//! Running a realignment end to end — stages A through D as one cancellable job.
//!
//! The stages exist as independent, separately-tested pieces
//! ([`revert`](navigator_analysis::revert), [`navigator_align`],
//! [`postprocess`](navigator_analysis::postprocess), and `crate::realign`). This is the part that
//! knows the order, what to do between them, and how to stop.
//!
//! ## Shape of the job
//!
//! ```text
//! preflight ──► revert ──► index ──► map ──► sort ──► markdup ──► finalize ──► register
//! ```
//!
//! It is hours of work on tens of GB of scratch, so three things matter more than they would in a
//! shorter job:
//!
//! - **Nothing is destroyed.** The source alignment is read and never written. If any stage fails
//!   or is cancelled, the workspace is exactly as it was — the new row is inserted last, after the
//!   final artifact exists, so there is no window where a half-built alignment is registered.
//! - **Scratch is cleaned up.** Intermediates are the size of the input, several times over.
//!   [`JobScratch`] removes them on the way out whether the job succeeded, failed, or was
//!   cancelled.
//! - **It can be stopped.** Every stage takes the cancel token, and a cancelled job reports itself
//!   as cancelled rather than as a failure — the user pressed the button.
//!
//! ## Preflight
//!
//! Checking disk and memory before starting is not politeness. Filling a disk partway through hour
//! three fails the job *and* leaves the machine unusable until someone finds the scratch directory;
//! and the index sizing needs the RAM figure anyway. See [`preflight`].

use std::path::{Path, PathBuf};

use navigator_align::{BatchSize, Preset};
use navigator_analysis::cancel::CancelToken;
use navigator_analysis::postprocess::{self, MarkDupParams, SortParams};
use navigator_analysis::revert::{self, RevertParams};
use navigator_domain::workspace::Alignment;

use crate::error::AppError;
use crate::App;

/// The stages, in order, for progress reporting.
///
/// Named rather than numbered in the UI: "sorting" tells a user waiting two hours something that
/// "step 5 of 8" does not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RealignStage {
    Preflight,
    Revert,
    Index,
    Map,
    Sort,
    MarkDuplicates,
    Finalize,
    Register,
}

impl RealignStage {
    pub const ALL: [RealignStage; 8] = [
        RealignStage::Preflight,
        RealignStage::Revert,
        RealignStage::Index,
        RealignStage::Map,
        RealignStage::Sort,
        RealignStage::MarkDuplicates,
        RealignStage::Finalize,
        RealignStage::Register,
    ];

    /// What to show the user while this stage runs.
    pub fn label(self) -> &'static str {
        match self {
            RealignStage::Preflight => "Checking space and memory",
            RealignStage::Revert => "Recovering the original reads",
            RealignStage::Index => "Preparing the reference index",
            RealignStage::Map => "Mapping reads to the new reference",
            RealignStage::Sort => "Sorting by position",
            RealignStage::MarkDuplicates => "Marking duplicates",
            RealignStage::Finalize => "Indexing the new alignment",
            RealignStage::Register => "Adding it to the workspace",
        }
    }

    pub fn step(self) -> usize {
        RealignStage::ALL.iter().position(|s| *s == self).unwrap_or(0) + 1
    }
}

/// Progress from a running job.
#[derive(Debug, Clone)]
pub struct RealignProgress {
    pub stage: RealignStage,
    pub total_stages: usize,
    /// Free-text detail — a record count, a part number. May be empty.
    pub detail: String,
}

/// What a finished job produced.
///
/// The counts are optional because a resumed job did not necessarily run the stage that produces
/// them. A run that picks up from a previous attempt's sorted BAM never reverted anything, so it
/// cannot report how many unmapped reads that revert saw; [`ScratchState`] carries the figure
/// across when the earlier attempt recorded it, and `None` says plainly that nobody measured it
/// rather than reporting a zero that reads like a finding.
#[derive(Debug, Clone)]
pub struct RealignOutcome {
    pub alignment: Alignment,
    /// Reads that were unmapped in the source and therefore had a chance to be placed by the new
    /// reference. This is the number the module exists to move, so it is surfaced rather than
    /// buried in a log.
    pub source_unmapped_reads: Option<u64>,
    pub reads_written: Option<u64>,
    pub duplicates_marked: Option<u64>,
}

/// Tuning a caller may override; the defaults are what the UI uses.
#[derive(Debug, Clone)]
pub struct RealignParams {
    /// Target build, e.g. `chm13v2.0`.
    pub target_build: String,
    /// The target build's FASTA, already resolved and indexed by `navigator-refgenome`.
    pub target_reference: PathBuf,
    /// `None` infers from the run's technology, which is what the UI does.
    pub preset: Option<Preset>,
    /// Where intermediates live. `None` puts them beside the output.
    pub scratch_root: Option<PathBuf>,
    /// Pick up from a previous attempt's intermediates instead of starting over.
    ///
    /// Off by default, because reusing files a *different* job left behind would be a correctness
    /// bug, and the scratch path alone cannot prove they came from this source and this target.
    /// The caller opting in is what supplies that knowledge. See [`Resumed`].
    pub resume: bool,
}

/// How far a previous attempt got, judged from what it left in the scratch directory.
///
/// Ordered by how much work it saves, and checked newest-artifact-first: a complete `marked.bam`
/// makes the sort behind it irrelevant, so there is no reason to ask about `sorted.bam` as well.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
enum Resumed {
    /// Start from the source alignment.
    Nothing,
    /// Reads were recovered and mapped; sorting is next.
    Mapped,
    /// The mapped BAM was sorted; duplicate marking is next.
    Sorted,
    /// Duplicates were marked; only finalising and registration remain.
    Marked,
}

impl Resumed {
    /// What to tell the user about the stages this skips.
    fn detail(self) -> &'static str {
        match self {
            Resumed::Nothing => "",
            Resumed::Mapped => "reusing the mapped reads from an earlier attempt",
            Resumed::Sorted => "reusing the sorted alignment from an earlier attempt",
            Resumed::Marked => "reusing the marked alignment from an earlier attempt",
        }
    }
}

/// Counts a resumed job cannot re-derive, left beside the intermediates they describe.
///
/// Each stage's contribution to [`RealignOutcome`] is measured while that stage runs and is gone
/// once it has. A resumed job skips stages by design, so the numbers are written down as they are
/// produced. Best-effort throughout: failing to write this file must not fail a job that has
/// otherwise done hours of correct work, and a missing or unreadable file simply means the counts
/// come back `None`.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
struct ScratchState {
    unmapped_reads: Option<u64>,
    reads_written: Option<u64>,
    duplicates_marked: Option<u64>,
    /// The furthest stage that *returned successfully*, as opposed to merely leaving a file behind.
    ///
    /// Written after the stage returns, never before, so it says something the file itself cannot:
    /// see [`discard_partial`] for why a finished-looking BAM is not proof on its own. A scratch
    /// directory predating this field has `None`, and then the marker is all there is to go on.
    completed_through: Option<Resumed>,
}

impl ScratchState {
    const FILE: &'static str = "state.json";

    fn load(scratch: &Path) -> Self {
        std::fs::read(scratch.join(Self::FILE))
            .ok()
            .and_then(|bytes| serde_json::from_slice(&bytes).ok())
            .unwrap_or_default()
    }

    fn store(&self, scratch: &Path) {
        if let Ok(bytes) = serde_json::to_vec_pretty(self) {
            let _ = std::fs::write(scratch.join(Self::FILE), bytes);
        }
    }
}

/// Remove a stage's output when the stage did not finish.
///
/// The BGZF end-of-file marker cannot carry the whole weight of "this file is complete", and
/// finding that out the hard way is what this function exists to prevent. noodles' multithreaded
/// writer finishes its stream from `Drop`, so a stage that *unwinds* — cancelled, failed, panicked
/// — leaves a partial file wearing a finished file's marker. Measured: a merge cancelled at 13.2 GB
/// of an expected ~30 GB, whose output `is_complete_bam` then accepted. Resuming from it would have
/// marked duplicates on a truncated alignment and registered the result, with nothing anywhere
/// saying the reads were missing.
///
/// A hard kill is the opposite case, and the one resume was built for: no `Drop` runs, no marker is
/// written, and the file is correctly refused. So the rule that makes the marker trustworthy again
/// is this one — whenever the job is still alive to notice a stage did not finish, its output does
/// not survive to be mistaken for one that did.
fn discard_partial(path: &Path) {
    let _ = std::fs::remove_file(path);
}

/// Run a stage, and make sure a failure leaves nothing that looks like a success.
async fn stage<T>(output: &Path, work: impl std::future::Future<Output = Result<T, AppError>>) -> Result<T, AppError> {
    match work.await {
        Ok(value) => Ok(value),
        Err(e) => {
            discard_partial(output);
            Err(e)
        }
    }
}

/// How far a previous attempt in `scratch` got, or [`Resumed::Nothing`] if it left nothing usable.
///
/// Two things have to agree. The file must carry the BGZF end-of-file marker, and — when the
/// previous attempt was new enough to have recorded one — its own account of the last stage it
/// completed must reach at least as far. The lower of the two wins, because each catches a case the
/// other cannot: the marker catches a scratch directory written by an older build, and the record
/// catches a file whose marker was written by an unwinding `Drop`. See [`discard_partial`].
fn resumable(scratch: &Path, state: &ScratchState) -> Resumed {
    let by_marker = resumable_by_marker(scratch);
    match state.completed_through {
        Some(claimed) => by_marker.min(claimed),
        None => by_marker,
    }
}

/// The marker-only half of [`resumable`], kept separate so the two rules can be tested apart.
fn resumable_by_marker(scratch: &Path) -> Resumed {
    if postprocess::is_complete_bam(&scratch.join("marked.bam")) {
        Resumed::Marked
    } else if postprocess::is_complete_bam(&scratch.join("sorted.bam")) {
        Resumed::Sorted
    } else if postprocess::is_complete_bam(&scratch.join("mapped.bam")) {
        Resumed::Mapped
    } else {
        Resumed::Nothing
    }
}

/// Clear a stage's working directory before that stage runs.
///
/// A resumed job re-runs the first stage it cannot skip, and that stage's leftovers from the
/// killed attempt are pure cost: the sort ignores run files it did not write itself, so stale ones
/// are not a correctness problem, but at WGS scale they are tens of GB held against a disk that
/// the same job is about to need. Best-effort — a directory that will not clear is not a reason to
/// refuse to start.
fn clear_stage_dir(dir: &Path) {
    let _ = std::fs::remove_dir_all(dir);
}

impl App {
    /// Realign `source_id` onto another reference, end to end.
    ///
    /// Long-running and cancellable. `progress` is called as each stage begins.
    pub async fn realign_alignment(
        &self,
        source_id: i64,
        params: RealignParams,
        cancel: CancelToken,
        progress: impl Fn(RealignProgress) + Send + 'static,
    ) -> Result<RealignOutcome, AppError> {
        let source = self.alignment_or_err(source_id).await?;
        let source_bam = source.bam_path.clone().ok_or(AppError::MissingPaths(source_id))?;
        let source_reference = source.reference_path.clone().map(PathBuf::from);

        let preset = match params.preset {
            Some(p) => p,
            None => self.infer_preset(&source).await?,
        };

        let output = self.realigned_output_path(source_id, &params.target_build);
        let scratch = params
            .scratch_root
            .clone()
            .unwrap_or_else(|| output.with_extension("scratch"));
        // Removes every intermediate on the way out, however the job ends — unless
        // [`keep_scratch_on_failure`] is set and the job dies.
        let mut scratch = JobScratch::new(scratch)?;

        let report = |stage: RealignStage, detail: &str| {
            progress(RealignProgress {
                stage,
                total_stages: RealignStage::ALL.len(),
                detail: detail.to_string(),
            });
        };

        let mapped = scratch.path().join("mapped.bam");
        let sorted = scratch.path().join("sorted.bam");
        let marked = scratch.path().join("marked.bam");

        // What a previous attempt left behind, and what it measured while it was running.
        let previous = ScratchState::load(scratch.path());
        let resumed = if params.resume {
            resumable(scratch.path(), &previous)
        } else {
            Resumed::Nothing
        };
        let mut state = if resumed == Resumed::Nothing {
            ScratchState::default()
        } else {
            previous
        };

        // Watch the machine for as long as the job runs. It reports; it never intervenes — see
        // `navigator_analysis::resource`. Started here so that it covers every stage, including the
        // ones a resumed job skips over quickly.
        let _watch = navigator_analysis::resource::ResourceWatch::start(
            navigator_analysis::resource::DEFAULT_INTERVAL,
            |sample| {
                // Anything short of trouble is a log line; the bands exist so that trouble is
                // greppable afterwards rather than buried in six hours of normal readings.
                if sample.pressure == navigator_analysis::resource::Pressure::Normal {
                    eprintln!("realign: {}", sample.summary());
                } else {
                    eprintln!("realign: WARNING {}", sample.summary());
                }
            },
        );

        // ---- preflight ----
        report(RealignStage::Preflight, resumed.detail());
        cancel.check()?;
        let source_size = std::fs::metadata(&source_bam).map(|m| m.len()).unwrap_or(0);
        // A resumed job is sized from the intermediate it is resuming from, not from the source.
        // The source's expansion — CRAM to FASTQ and back to BAM — has already happened and is
        // sitting on the disk; charging for it a second time would refuse a job that fits.
        let plan = if resumed == Resumed::Nothing {
            preflight(scratch.path(), Path::new(&source_bam), source_size)?
        } else {
            resume_preflight(scratch.path(), &mapped, &sorted, &marked)?
        };

        // ---- stage A: revert ----
        report(RealignStage::Revert, resumed.detail());
        let reverted = if resumed > Resumed::Nothing {
            None
        } else {
            let bam = PathBuf::from(&source_bam);
            let reference = source_reference.clone();
            let dir = scratch.path().join("revert");
            clear_stage_dir(&dir);
            let token = cancel.clone();
            let stats = tokio::task::spawn_blocking(move || {
                revert::revert_alignment(&bam, reference.as_deref(), &dir, &RevertParams::default(), &token)
            })
            .await
            .map_err(|e| AppError::Join(e.to_string()))??;

            state.unmapped_reads = Some(stats.stats.unmapped_reads);
            state.store(scratch.path());
            Some(stats)
        };

        // ---- stage B: index, then map ----
        report(RealignStage::Index, resumed.detail());
        if let Some(reverted) = &reverted {
            let index = {
                let build = params.target_build.clone();
                let reference = params.target_reference.clone();
                let batch = plan.batch;
                tokio::task::spawn_blocking(move || {
                    navigator_align::index::ensure_index(
                        &navigator_align::index::cache_root(),
                        &build,
                        &reference,
                        preset,
                        batch,
                        &mut |_, _| {},
                    )
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))??
            };

            report(RealignStage::Map, "");
            let (r1, r2, singles) = (
                reverted.read1.clone(),
                reverted.read2.clone(),
                reverted.singletons.clone(),
            );
            let out = mapped.clone();
            let dir = scratch.path().join("map");
            clear_stage_dir(&dir);
            let token = cancel.clone();
            let map_params = navigator_align::MapParams {
                preset,
                threads: 0,
                read_group: None,
                format: navigator_align::OutputFormat::Bam,
                reference: None,
            };
            stage(&mapped, async {
                tokio::task::spawn_blocking(move || -> Result<(), AppError> {
                    let cancelled = move || token.is_cancelled();
                    if preset.is_paired() {
                        navigator_align::map_pairs(
                            &index,
                            &r1,
                            &r2,
                            &out,
                            &dir,
                            &map_params,
                            &cancelled,
                            &mut |_, _, _| {},
                        )?;
                    } else {
                        // Long-read presets are single-ended; the revert stage puts everything that
                        // did not pair into the singletons file, which is the whole read set here.
                        navigator_align::map_reads(
                            &index,
                            &singles,
                            &out,
                            &dir,
                            &map_params,
                            &cancelled,
                            &mut |_, _, _| {},
                        )?;
                    }
                    Ok(())
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))?
            })
            .await?;

            state.completed_through = Some(Resumed::Mapped);
            state.store(scratch.path());
        } else {
            // Still reported, even though there is nothing to do: a progress display that skips
            // from stage 3 to stage 5 reads as a missing step rather than a saved one.
            report(RealignStage::Map, resumed.detail());
        }

        // Every stage's input is dead once the next stage has read it, and at WGS scale each is
        // tens of GB. Holding them all until the job ends — which is what this did first — roughly
        // doubles the peak and is the difference between fitting on a normal disk and not.
        let discard = discard_partial;
        if let Some(reverted) = &reverted {
            discard(&reverted.read1);
            discard(&reverted.read2);
            discard(&reverted.singletons);
        }

        // ---- stage C: sort, mark duplicates, compress ----
        report(RealignStage::Sort, resumed.detail());
        if resumed < Resumed::Sorted {
            let (input, out) = (mapped.clone(), sorted.clone());
            let dir = scratch.path().join("sort");
            // A resumed job re-sorts from the start, so the killed attempt's spilled runs are dead
            // weight — and at WGS scale they are tens of GB the sort is about to want back.
            clear_stage_dir(&dir);
            let token = cancel.clone();
            stage(&sorted, async {
                tokio::task::spawn_blocking(move || {
                    postprocess::sort_alignment(&input, &out, &dir, &SortParams::default(), &token, &mut |_| {})
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))?
                .map_err(AppError::from)
            })
            .await?;

            // The runs are consumed by the merge; the sorted BAM is what the next stage reads.
            clear_stage_dir(&scratch.path().join("sort"));
            state.completed_through = Some(Resumed::Sorted);
            state.store(scratch.path());
        }

        // Only what *this* run produced is disposable. A resumed run that skipped the sort has
        // re-derived nothing, and deleting the mapped BAM it resumed past would throw away the most
        // expensive artifact in the pipeline on the strength of a belief about a file it did not
        // write. That is not hypothetical: combined with the marker bug in `discard_partial`, it is
        // exactly how a 59 GB `mapped.bam` — 3 h 58 m of revert and mapping — was destroyed.
        if resumed < Resumed::Sorted {
            discard(&mapped);
        }

        report(RealignStage::MarkDuplicates, resumed.detail());
        if resumed < Resumed::Marked {
            let (input, out) = (sorted.clone(), marked.clone());
            let token = cancel.clone();
            // Long-read libraries are typically PCR-free and long reads rarely share endpoints, so
            // marking them would discard real coverage.
            let md_params = MarkDupParams {
                enabled: preset.is_paired(),
                ..Default::default()
            };
            let markdup = stage(&marked, async {
                tokio::task::spawn_blocking(move || {
                    postprocess::mark_duplicates(&input, &out, &md_params, &token, &mut |_| {})
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))?
                .map_err(AppError::from)
            })
            .await?;

            state.reads_written = Some(markdup.records);
            state.duplicates_marked = Some(markdup.duplicates);
            state.completed_through = Some(Resumed::Marked);
            state.store(scratch.path());
        }

        if resumed < Resumed::Marked {
            discard(&sorted);
        }

        report(RealignStage::Finalize, resumed.detail());
        let finalized = {
            let (input, out) = (marked.clone(), output.clone());
            tokio::task::spawn_blocking(move || postprocess::finalize_bam(&input, &out))
                .await
                .map_err(|e| AppError::Join(e.to_string()))??
        };

        // The output exists, so the intermediates behind it are disposable even if registration
        // then fails — re-registering is seconds of work, not hours. `marked` is not discarded
        // here because finalising *moved* it into place rather than copying it.
        scratch.completed();

        // ---- stage D: register ----
        //
        // Last, deliberately. The row is only inserted once the artifact it points at exists, so a
        // job that fails or is cancelled leaves no alignment referring to a file that was never
        // finished.
        report(RealignStage::Register, "");
        cancel.check()?;
        let alignment = self
            .register_realigned_alignment(
                source_id,
                &finalized.bam,
                &params.target_reference,
                &params.target_build,
                "minimap2",
                preset.as_str(),
            )
            .await?;

        Ok(RealignOutcome {
            alignment,
            source_unmapped_reads: state.unmapped_reads,
            reads_written: state.reads_written,
            duplicates_marked: state.duplicates_marked,
        })
    }

    /// The mapper preset for a source alignment, from its run's recorded technology.
    async fn infer_preset(&self, source: &Alignment) -> Result<Preset, AppError> {
        let run = navigator_store::sequence_run::get(self.store.pool(), source.sequence_run_id)
            .await?
            .ok_or(AppError::MissingPaths(source.id))?;
        Preset::infer(Some(&run.test_type), Some(&run.platform_name)).map_err(|e| {
            // A refusal here is the right outcome — mapping long reads under a short-read preset
            // produces plausible, wrong alignments rather than failing — so it is surfaced as-is.
            AppError::Import(format!(
                "cannot choose a mapper preset for alignment #{}: {e}",
                source.id
            ))
        })
    }
}

/// What preflight decided.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RealignPlan {
    /// Index batch size, chosen from the machine's memory.
    pub batch: BatchSize,
    /// Bytes of scratch the job is expected to need.
    pub scratch_needed: u64,
    /// Bytes free where the scratch will live.
    pub scratch_free: u64,
}

/// Scratch multiple of the source's **uncompressed** volume.
///
/// Calibrated against a measured run, not reasoned from stage sizes — the reasoned figure was 25%
/// low. Realigning WGS229 (17.3 GB CRAM) peaked at **276 GB of scratch**, i.e. 16x the source
/// file, which this reproduces as `4` here times the CRAM expansion factor below.
///
/// The peak is inside the revert, not where the stage list suggests: its spill runs are the
/// bespoke uncompressed encoding while its FASTQ output is gzipped, so the two coexist at very
/// different densities and the sum beats anything later in the pipeline.
const SCRATCH_MULTIPLE: u64 = 4;

/// How much larger the data is than the file holding it.
///
/// This is the correction that matters, and getting it wrong is not a rounding error. A CRAM is
/// reference-compressed: 17 GB of CRAM is ~70 GB of BAM-equivalent data and, reverted to FASTQ,
/// larger still — FASTQ spends an ASCII byte per base on quality where BAM packs. Estimating
/// scratch from the *file* size would have told a user that a job needing ~200 GB needed 69, and
/// the disk would have filled somewhere in the third hour.
fn expansion_factor(source: &Path) -> u64 {
    match source.extension().and_then(|e| e.to_str()) {
        Some(e) if e.eq_ignore_ascii_case("cram") => 4,
        // BAM is bgzf-compressed at roughly 4x, but the reverted FASTQ that comes out of it is
        // gzipped too, so the ratio between them is far closer to 1.
        _ => 2,
    }
}

/// Check that the job can finish before starting it.
///
/// Running out of disk in hour three fails the job *and* leaves the machine wedged until someone
/// finds the scratch directory. The RAM figure is needed regardless, to size the index.
pub fn preflight(scratch: &Path, source: &Path, source_size: u64) -> Result<RealignPlan, AppError> {
    let needed = source_size
        .saturating_mul(expansion_factor(source))
        .saturating_mul(SCRATCH_MULTIPLE);
    plan_for(scratch, needed, "realign")
}

/// Preflight for a job resuming from intermediates that already exist.
///
/// The full [`preflight`] estimate is anchored to the *source* file, and multiplies it by the
/// expansion the pipeline is about to perform. For a resumed job that expansion is already on the
/// disk, so the same sum charges twice for it — enough to refuse a job with room to spare.
///
/// The remaining stages are each roughly the size of the intermediate they read, and at most two
/// coexist: the sort holds its spilled runs while the mapped BAM is still there, and duplicate
/// marking writes its output while the sorted BAM is still there. Three times the largest
/// surviving intermediate covers that with a margin, and the intermediate's size is a measurement
/// rather than a guess.
fn resume_preflight(scratch: &Path, mapped: &Path, sorted: &Path, marked: &Path) -> Result<RealignPlan, AppError> {
    let size = |path: &Path| std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let largest = size(mapped).max(size(sorted)).max(size(marked));
    plan_for(scratch, largest.saturating_mul(3), "resume the realignment")
}

/// Measure the disk, refuse a job that cannot finish on it, and describe what was decided.
///
/// The two preflights differ only in how they size `needed`; everything after that — probing free
/// space, the refusal, the wording, the plan — was written out twice and had to be kept in step by
/// hand. `what` is the verb in the refusal, so the two messages stay exactly as they were.
fn plan_for(scratch: &Path, needed: u64, what: &str) -> Result<RealignPlan, AppError> {
    let free = free_space(scratch);

    if !has_room(needed, free) {
        return Err(AppError::Import(format!(
            "not enough room to {what}: about {} GB of working space is needed and {} GB is free \
             on {}",
            gb(needed),
            gb(free),
            scratch.display(),
        )));
    }

    Ok(RealignPlan {
        // A resumed job never reaches the index stage, but the plan is one shape and a caller may
        // still log the figure.
        batch: BatchSize::for_this_machine(),
        scratch_needed: needed,
        scratch_free: free,
    })
}

fn gb(bytes: u64) -> u64 {
    bytes / 1_000_000_000
}

/// Whether a job needing `needed` bytes may start given `free` bytes available.
///
/// Split from the syscall so the *decision* is testable on any machine; the probe itself is not.
/// `free == 0` means the platform would not say, and is treated as permission to proceed —
/// refusing a multi-hour job because a free-space call failed would be a worse outcome than
/// letting it run and fail honestly on a real write.
fn has_room(needed: u64, free: u64) -> bool {
    free == 0 || free >= needed
}

/// Free bytes on the filesystem holding `path`, or 0 when it cannot be determined.
///
/// Zero means "unknown", and preflight treats it as "do not block" — refusing a job because the
/// free-space call failed would be worse than letting it run and fail honestly on a real write.
fn free_space(path: &Path) -> u64 {
    // Walk up to the nearest existing ancestor: the scratch directory itself may not exist yet.
    let mut probe = path;
    loop {
        if probe.exists() {
            break;
        }
        match probe.parent() {
            Some(parent) => probe = parent,
            None => return 0,
        }
    }
    fs_free_space(probe)
}

#[cfg(unix)]
fn fs_free_space(path: &Path) -> u64 {
    use std::os::unix::ffi::OsStrExt;
    let Ok(c_path) = std::ffi::CString::new(path.as_os_str().as_bytes()) else {
        return 0;
    };
    // SAFETY: `statvfs` writes into a zeroed struct we own, and `c_path` is a valid NUL-terminated
    // string that outlives the call.
    unsafe {
        let mut stat: libc::statvfs = std::mem::zeroed();
        if libc::statvfs(c_path.as_ptr(), &mut stat) != 0 {
            return 0;
        }
        (stat.f_bavail as u64).saturating_mul(stat.f_frsize as u64)
    }
}

/// Free bytes on the volume holding `path`.
///
/// `lpFreeBytesAvailableToCaller` rather than the volume's total free space, which is the subtle
/// half of this call: on a volume with disk quotas the two differ, and the number preflight needs is
/// what *this user* may actually write. It is the same choice the Unix arm makes in taking
/// `f_bavail` (blocks available to an unprivileged process) over `f_bfree`.
#[cfg(windows)]
fn fs_free_space(path: &Path) -> u64 {
    use std::os::windows::ffi::OsStrExt;
    use windows_sys::Win32::Storage::FileSystem::GetDiskFreeSpaceExW;

    // `encode_wide` does not NUL-terminate, and the API requires it.
    let mut wide: Vec<u16> = path.as_os_str().encode_wide().collect();
    wide.push(0);

    let mut available: u64 = 0;
    // SAFETY: `wide` is a valid NUL-terminated UTF-16 buffer owned here and outliving the call, and
    // `available` is a `u64` we own. The two totals this does not need are passed as null, which the
    // API documents as permitted.
    let ok = unsafe {
        GetDiskFreeSpaceExW(
            wide.as_ptr(),
            &mut available,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        )
    };

    // A failure is reported as "unknown" rather than as zero free space, so that a probe that could
    // not answer does not masquerade as a full disk and refuse the job. See `has_room`.
    if ok == 0 {
        return 0;
    }
    available
}

#[cfg(not(any(unix, windows)))]
fn fs_free_space(_path: &Path) -> u64 {
    // No probe for this platform. Preflight reports "unknown" and declines to block — see
    // `free_space`.
    0
}

/// A scratch directory that removes itself.
///
/// Intermediates are several times the size of the input, so leaving them behind after a failed or
/// cancelled job would quietly consume a disk. Removal is best-effort: it must never mask the error
/// that is already unwinding.
struct JobScratch {
    path: PathBuf,
    /// Set once the job has produced its output, so [`Drop`] can tell "finished" from "died".
    completed: bool,
    /// Keep the intermediates behind when the job does *not* complete.
    keep_on_failure: bool,
}

/// Opt in to keeping a failed job's intermediates (`NAVIGATOR_REALIGN_KEEP_SCRATCH=1`).
///
/// Off by default, and that default is deliberate: at WGS scale this directory is hundreds of GB,
/// and a desktop application that leaves that behind after a failure has taken the user's disk
/// hostage over a bug they did not cause.
///
/// It exists because the opposite default has a matching failure of its own. A realignment that
/// dies in stage 7 of 8 discards the seven stages that worked — a measured 10.7 hours on a 30x WGS.
/// For anyone iterating on the pipeline, that trade is the wrong way round, so they can invert it.
///
/// This is now also the switch that makes [`RealignParams::resume`] worth anything: resume reads
/// the intermediates a previous attempt left, and by default a failed attempt does not leave any.
/// The two are meant to be used together when a run is expected to be interrupted — which, at six
/// hours a run on a machine someone is still using, is a reasonable thing to expect. A job killed
/// outright (a session teardown, a power loss) never runs `Drop` at all, so its scratch survives
/// regardless of this setting; that is the case resume was built for.
fn keep_scratch_on_failure() -> bool {
    std::env::var("NAVIGATOR_REALIGN_KEEP_SCRATCH")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false)
}

impl JobScratch {
    fn new(path: PathBuf) -> Result<Self, AppError> {
        std::fs::create_dir_all(&path)?;
        Ok(Self {
            path,
            completed: false,
            keep_on_failure: keep_scratch_on_failure(),
        })
    }

    fn path(&self) -> &Path {
        &self.path
    }

    /// The job got its output. Anything left here is now genuinely disposable.
    fn completed(&mut self) {
        self.completed = true;
    }
}

impl Drop for JobScratch {
    fn drop(&mut self) {
        if !self.completed && self.keep_on_failure {
            eprintln!(
                "realign: keeping intermediates at {} (NAVIGATOR_REALIGN_KEEP_SCRATCH)",
                self.path.display()
            );
            return;
        }
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stages_are_labelled_and_numbered_in_order() {
        assert_eq!(RealignStage::Preflight.step(), 1);
        assert_eq!(RealignStage::Register.step(), RealignStage::ALL.len());
        for stage in RealignStage::ALL {
            assert!(!stage.label().is_empty(), "{stage:?} needs a label");
        }
    }

    /// Preflight must refuse a job that cannot finish, and say how much is needed rather than
    /// leaving the user to guess.
    ///
    /// Runs wherever [`fs_free_space`] has a real implementation, which is now both desktop
    /// families. Platforms without one report "unknown" and cannot refuse anything — pinned
    /// separately by `preflight_cannot_refuse_where_free_space_is_unknown`.
    #[cfg(any(unix, windows))]
    #[test]
    fn preflight_refuses_when_the_disk_is_too_small() {
        let dir = std::env::temp_dir();
        // A source larger than any plausible free space.
        let err = preflight(&dir, Path::new("x.bam"), u64::MAX / 8).expect_err("must refuse");
        let message = format!("{err}");
        assert!(message.contains("working space"), "unhelpful: {message}");
        assert!(message.contains("free"), "unhelpful: {message}");
    }

    #[test]
    fn preflight_accepts_a_small_job_and_sizes_the_index() {
        let plan = preflight(&std::env::temp_dir(), Path::new("x.bam"), 1024).expect("a 1 KB source always fits");
        assert_eq!(plan.scratch_needed, 1024 * 2 * SCRATCH_MULTIPLE);
        assert!(plan.batch.bases() > 0);
    }

    /// The correction that matters: a CRAM holds several times its own size in read data, so
    /// estimating scratch from the file size understates a WGS job by well over 100 GB.
    #[test]
    fn a_cram_is_estimated_larger_than_a_bam_of_the_same_size() {
        let dir = std::env::temp_dir();
        let bam = preflight(&dir, Path::new("s.bam"), 1_000_000).unwrap();
        let cram = preflight(&dir, Path::new("s.cram"), 1_000_000).unwrap();
        assert!(
            cram.scratch_needed > bam.scratch_needed,
            "a CRAM expands more than a BAM: {} vs {}",
            cram.scratch_needed,
            bam.scratch_needed
        );
        assert_eq!(cram.scratch_needed, 1_000_000 * 4 * SCRATCH_MULTIPLE);
    }

    /// Free space is a check, not a gate on the check working: if the platform will not say, the
    /// job runs and fails honestly on a real write rather than being refused for an unrelated
    /// reason.
    #[test]
    fn an_unknown_free_space_does_not_block_the_job() {
        assert!(has_room(u64::MAX, 0), "0 means unknown, not full");
        assert!(has_room(100, 100), "exactly enough is enough");
        assert!(!has_room(101, 100), "one byte short is short");
    }

    /// The scratch directory does not exist when preflight runs — it is created by the job — so the
    /// probe has to answer for the filesystem that will hold it, which means walking up to the
    /// nearest existing ancestor rather than giving up.
    ///
    /// Both desktop families, for the same reason as above: the ancestor walk is platform-
    /// independent, but it is only meaningful where the call it ends at can answer.
    #[cfg(any(unix, windows))]
    #[test]
    fn free_space_resolves_through_a_directory_that_does_not_exist_yet() {
        let unborn = std::env::temp_dir().join("dun-not-created-yet").join("nor-this");
        assert!(!unborn.exists());
        assert!(
            free_space(&unborn) > 0,
            "must report the filesystem that will hold the scratch"
        );
    }

    /// What the platforms without a free-space probe actually do, stated as a test rather than left
    /// as the absence of one.
    ///
    /// `preflight` cannot refuse a job it has no measurement for, and refusing on an unknown would
    /// block every realignment on that platform. So a job that would obviously not fit is allowed
    /// to start and fail honestly on a real write. Windows used to be in this bucket; it now has
    /// `GetDiskFreeSpaceExW` and is tested by the two above.
    #[cfg(not(any(unix, windows)))]
    #[test]
    fn preflight_cannot_refuse_where_free_space_is_unknown() {
        let dir = std::env::temp_dir();
        assert_eq!(free_space(&dir), 0, "no probe on this platform yet");

        let plan = preflight(&dir, Path::new("x.bam"), u64::MAX / 8).expect("unknown must not block");
        assert_eq!(plan.scratch_free, 0);
        assert!(plan.scratch_needed > 0, "the estimate is still made and reported");
    }

    /// Scratch is several times the size of the input; a cancelled job must not leave it behind.
    #[test]
    fn scratch_is_removed_when_the_job_ends() {
        let dir = std::env::temp_dir().join(format!("dun-jobscratch-{}", std::process::id()));
        {
            let scratch = JobScratch::new(dir.clone()).unwrap();
            std::fs::write(scratch.path().join("intermediate"), b"big").unwrap();
            assert!(dir.exists());
        }
        assert!(!dir.exists(), "scratch outlived the job");
    }
}

#[cfg(test)]
mod resume_tests {
    use super::*;

    /// A scratch directory of this test's own, holding whatever files it names.
    fn scratch(tag: &str, files: &[(&str, bool)]) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("navigator-resume-{tag}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        for (name, complete) in files {
            let path = dir.join(name);
            // A "complete" BAM is one carrying the BGZF end-of-file marker; an incomplete one is
            // the same bytes with the marker cut short, which is what a killed writer leaves.
            let mut bytes = vec![0u8; 64];
            let eof: [u8; 28] = [
                0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b,
                0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            ];
            if *complete {
                bytes.extend_from_slice(&eof);
            } else {
                bytes.extend_from_slice(&eof[..20]);
            }
            std::fs::write(&path, &bytes).unwrap();
        }
        dir
    }

    #[test]
    fn an_empty_scratch_directory_resumes_nothing() {
        let dir = scratch("empty", &[]);
        assert_eq!(resumable(&dir, &ScratchState::default()), Resumed::Nothing);
    }

    /// The bug this pair of rules exists for, reproduced.
    ///
    /// A cancelled merge left 13.2 GB of an expected ~30 GB `sorted.bam` carrying a valid BGZF
    /// end-of-file marker, because noodles' multithreaded writer finishes its stream from `Drop`
    /// while the job unwinds. The marker check alone accepted it and the next run went straight to
    /// marking duplicates on a truncated alignment. The stage record is what refuses it: the sort
    /// never returned, so nothing ever claimed it completed.
    #[test]
    fn a_marker_written_by_an_unwinding_drop_is_not_enough() {
        let dir = scratch("unwound", &[("mapped.bam", true), ("sorted.bam", true)]);

        assert_eq!(
            resumable_by_marker(&dir),
            Resumed::Sorted,
            "the marker alone is fooled — that is the premise"
        );

        let state = ScratchState {
            completed_through: Some(Resumed::Mapped),
            ..Default::default()
        };
        assert_eq!(
            resumable(&dir, &state),
            Resumed::Mapped,
            "the stage record caps it at what actually finished"
        );
    }

    /// The record cannot promise more than the files deliver either — a scratch directory whose
    /// `marked.bam` was removed must not be resumed from just because a stale record mentions it.
    #[test]
    fn the_record_cannot_outrun_the_files() {
        let dir = scratch("stale-record", &[("mapped.bam", true)]);
        let state = ScratchState {
            completed_through: Some(Resumed::Marked),
            ..Default::default()
        };

        assert_eq!(resumable(&dir, &state), Resumed::Mapped);
    }

    /// A scratch directory written before the stage record existed still has to work: the Aug-13
    /// `mapped.bam` was left by a process killed outright, so no `Drop` ran and the marker means
    /// exactly what it says.
    #[test]
    fn a_scratch_without_a_record_falls_back_to_the_marker() {
        let dir = scratch("legacy", &[("mapped.bam", true)]);
        assert_eq!(resumable(&dir, &ScratchState::default()), Resumed::Mapped);
    }

    /// The 2026-08-13 case exactly: mapping finished, the sort was killed partway through writing
    /// its output. The mapped BAM is worth four hours and must be picked up; the half-written
    /// sorted BAM must not be.
    #[test]
    fn a_complete_map_and_a_truncated_sort_resumes_from_the_map() {
        let dir = scratch("mid-sort", &[("mapped.bam", true), ("sorted.bam", false)]);
        assert_eq!(resumable(&dir, &ScratchState::default()), Resumed::Mapped);
    }

    #[test]
    fn the_furthest_complete_stage_wins() {
        let dir = scratch(
            "furthest",
            &[("mapped.bam", true), ("sorted.bam", true), ("marked.bam", true)],
        );
        assert_eq!(resumable(&dir, &ScratchState::default()), Resumed::Marked);

        let dir = scratch("furthest-sorted", &[("mapped.bam", true), ("sorted.bam", true)]);
        assert_eq!(resumable(&dir, &ScratchState::default()), Resumed::Sorted);
    }

    /// The stage guards are written as `resumed < Resumed::Sorted`, so the ordering is load-bearing
    /// rather than cosmetic: getting it backwards would skip work that has not been done.
    #[test]
    fn the_stages_are_ordered_by_how_much_they_skip() {
        assert!(Resumed::Nothing < Resumed::Mapped);
        assert!(Resumed::Mapped < Resumed::Sorted);
        assert!(Resumed::Sorted < Resumed::Marked);
    }

    /// Resume must not inherit a *different* job's numbers, and it must not invent them either.
    #[test]
    fn counts_survive_a_round_trip_and_default_to_unmeasured() {
        let dir = scratch("state", &[]);
        assert_eq!(ScratchState::load(&dir).unmapped_reads, None);

        let state = ScratchState {
            unmapped_reads: Some(301_431),
            reads_written: Some(62_429_459),
            duplicates_marked: Some(2_122_650),
            completed_through: Some(Resumed::Marked),
        };
        state.store(&dir);

        let loaded = ScratchState::load(&dir);
        assert_eq!(loaded.unmapped_reads, Some(301_431));
        assert_eq!(loaded.reads_written, Some(62_429_459));
        assert_eq!(loaded.duplicates_marked, Some(2_122_650));
        assert_eq!(loaded.completed_through, Some(Resumed::Marked));
    }

    /// A resumed job is sized from what is already on disk. Sizing it from the source again would
    /// charge twice for an expansion that has already happened and refuse a job that fits.
    #[test]
    fn resume_preflight_sizes_from_the_surviving_intermediate() {
        let dir = scratch("preflight", &[("mapped.bam", true)]);
        let mapped = dir.join("mapped.bam");
        let plan = resume_preflight(&dir, &mapped, &dir.join("sorted.bam"), &dir.join("marked.bam"))
            .expect("a 92-byte intermediate fits anywhere");

        let size = std::fs::metadata(&mapped).unwrap().len();
        assert_eq!(plan.scratch_needed, size * 3);
    }
}
