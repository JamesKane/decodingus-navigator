//! One realignment from start to end. This module runs stage A to stage D as one job, and the user
//! can stop that job.
//!
//! Each stage is a separate piece of code with its own tests. They are
//! [`revert`](navigator_analysis::revert), [`navigator_align`],
//! [`postprocess`](navigator_analysis::postprocess), and `crate::realign`. This module knows their
//! order, the work between them, and the way to stop them.
//!
//! ## Shape of the job
//!
//! ```text
//! preflight ──► revert ──► index ──► map ──► sort ──► markdup ──► finalize ──► register
//! ```
//!
//! The job needs hours of work and tens of GB of scratch space. So three rules matter more here
//! than in a short job.
//!
//! - **The job destroys nothing.** It reads the source alignment and never writes to it. After a
//!   failed stage, and after a stop, the workspace holds what it held before. The code adds the new
//!   row last, after the final file exists. So the workspace never holds a row for an alignment
//!   that the job did not complete.
//! - **The job removes its scratch files.** Those files are some times the size of the input.
//!   [`JobScratch`] removes them at the end, after a success, after a failure, and after a stop.
//! - **The user can stop the job.** Each stage takes the cancel token. A job that the user stopped
//!   reports a stop and not a failure, because the user pressed the button.
//!
//! ## Preflight
//!
//! The job checks the disk and the memory before it starts, and that check is necessary.
//!
//! A disk that fills in the third hour fails the job. It also leaves the machine unusable until
//! somebody finds the scratch directory. The size of the index also needs the RAM value. See
//! [`preflight`].

use std::path::{Path, PathBuf};

use navigator_align::{BatchSize, Preset};
use navigator_analysis::cancel::CancelToken;
use navigator_analysis::postprocess::{self, MarkDupParams, SortParams};
use navigator_analysis::revert::{self, RevertParams};
use navigator_domain::workspace::Alignment;

use crate::error::AppError;
use crate::App;

/// The stages, in their order, for the progress report.
///
/// The UI shows a name and not a number. The word "sorting" tells a user who waits two hours more
/// than the text "step 5 of 8".
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

/// The progress of a job that is in operation.
#[derive(Debug, Clone)]
pub struct RealignProgress {
    pub stage: RealignStage,
    pub total_stages: usize,
    /// Text with more detail, such as a count of records or a part number. The value can be
    /// empty.
    pub detail: String,
}

/// The result of a job that is complete.
///
/// Each count is optional, because a job that continues an earlier job does not always run the stage
/// that makes that count.
///
/// A run that starts from the sorted BAM file of an earlier run does no revert step. So it can not
/// report the count of unmapped reads that the revert step saw.
///
/// [`ScratchState`] holds that value when the earlier run wrote it. A value of `None` states that no
/// code measured the count. A zero value looks like a measurement, and it would be wrong.
#[derive(Debug, Clone)]
pub struct RealignOutcome {
    pub alignment: Alignment,
    /// The count of reads with no place in the source alignment. The new reference can give each of
    /// them a place. This module exists to increase that count, so the report holds it. A log entry
    /// alone would hide it.
    pub source_unmapped_reads: Option<u64>,
    pub reads_written: Option<u64>,
    pub duplicates_marked: Option<u64>,
}

/// The values that a caller can change. The UI uses the default values.
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
    /// Use the intermediate files of an earlier run, and do not start again from the source.
    ///
    /// The default is off. A run that used the files of a *different* job gives a wrong result.
    /// The scratch path alone does not prove that those files came from this source and this
    /// target. The caller who sets this field supplies that knowledge. See [`Resumed`].
    pub resume: bool,
}

/// The last stage that an earlier run completed. The code decides from the files in the scratch
/// directory.
///
/// The order is the quantity of work that each value saves. The code tests the newest file first. A
/// complete `marked.bam` file makes the sort behind it unnecessary, so the code does not test
/// `sorted.bam` also.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
enum Resumed {
    /// Start from the source alignment.
    Nothing,
    /// The run recovered the reads and mapped them. The sort is next.
    Mapped,
    /// The run sorted the mapped BAM file. The duplicate mark step is next.
    Sorted,
    /// The run marked the duplicates. Only the finalize step and the register step remain.
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

/// The counts that a job which continues an earlier job can not calculate again. The file sits
/// beside the intermediate files that it describes.
///
/// Each stage measures its part of [`RealignOutcome`] while that stage runs, and that value is gone
/// after the stage ends. A job that continues an earlier job skips stages, by design. So the code
/// writes each number to this file as a stage produces it.
///
/// Each write is optional. A failed write of this file must not fail a job that did hours of correct
/// work. An absent file, and a file that the code can not read, give counts of `None`.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
struct ScratchState {
    unmapped_reads: Option<u64>,
    reads_written: Option<u64>,
    duplicates_marked: Option<u64>,
    /// The last stage that *returned a success*. A stage that only left a file behind does not
    /// count.
    ///
    /// The code writes this value after the stage returns, and never before. So the value states
    /// something that the file itself can not state. [`discard_partial`] gives the reason: a BAM
    /// file that looks complete is not proof on its own.
    ///
    /// A scratch directory from before this field holds `None`. The code then has the marker
    /// only.
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

/// Remove the output of a stage that did not complete.
///
/// The BGZF end-of-file marker alone does not prove that a file is complete. This function exists,
/// so that no user learns that fact from a broken result.
///
/// The writer of noodles uses many threads, and it closes its stream from `Drop`. So a stage that
/// *unwinds* leaves a partial file with the marker of a complete file. A stop, a failure, and a
/// panic each unwind.
///
/// One measurement showed this. A user stopped a merge at 13.2 GB of an expected 30 GB, and
/// `is_complete_bam` then accepted that output. A run from that file would mark the duplicates of a
/// short alignment and add the result to the workspace. No record anywhere would state that the
/// reads were absent.
///
/// A hard kill is the opposite case, and the resume feature exists for it. No `Drop` runs, the code
/// writes no marker, and it correctly refuses the file.
///
/// So this rule makes the marker reliable again. While the job runs and can see that a stage did not
/// complete, the output of that stage must not survive. If it survives, a later run reads it as the
/// output of a stage that did complete.
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

/// The last stage that an earlier run in `scratch` completed. The function returns
/// [`Resumed::Nothing`] when that run left nothing that this run can use.
///
/// Two values must agree. The file must hold the BGZF end-of-file marker. The earlier run must also
/// record a last stage that reaches at least as far. Only a run from a new enough build writes such
/// a record.
///
/// The lower of the two values wins, because each one finds a case that the other can not. The
/// marker finds a scratch directory from an older build. The record finds a file whose marker a
/// `Drop` call wrote while the code unwound. See [`discard_partial`].
fn resumable(scratch: &Path, state: &ScratchState) -> Resumed {
    let by_marker = resumable_by_marker(scratch);
    match state.completed_through {
        Some(claimed) => by_marker.min(claimed),
        None => by_marker,
    }
}

/// The half of [`resumable`] that reads the marker only. It is a separate function, so a test can
/// cover each of the two rules alone.
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

/// Empty the directory of a stage before that stage runs.
///
/// A job that continues an earlier job runs the first stage that it can not skip. The files of that
/// stage from the earlier run give nothing.
///
/// They are not a correctness problem, because the sort ignores a run file that it did not write.
/// But for a WGS sample they hold tens of GB on a disk that this same job needs.
///
/// The step is optional. A directory that the code can not empty is not a reason to refuse the
/// job.
fn clear_stage_dir(dir: &Path) {
    let _ = std::fs::remove_dir_all(dir);
}

impl App {
    /// Realign `source_id` onto another reference, from start to end.
    ///
    /// The job runs for a long time, and the user can stop it. The code calls `progress` at the
    /// start of each stage.
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
        // This value removes each intermediate file at the end of the job, after each result. The
        // one exception is a failed job with [`keep_scratch_on_failure`] set.
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

        // The files that an earlier run left, and the values that it measured.
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

        // Watch the machine while the job runs. This watch reports, and it never stops the job.
        // See `navigator_resource`. It starts here, so it covers each stage. That set holds the
        // stages that a job which continues an earlier job passes quickly.
        let _watch = navigator_resource::ResourceWatch::start(navigator_resource::DEFAULT_INTERVAL, |sample| {
            // A normal reading is one line in the log. The bands exist so that a user can find a
            // problem in that log later. Without them, six hours of normal readings hide it.
            if sample.pressure == navigator_resource::Pressure::Normal {
                eprintln!("realign: {}", sample.summary());
            } else {
                eprintln!("realign: WARNING {}", sample.summary());
            }
        });

        // ---- preflight ----
        report(RealignStage::Preflight, resumed.detail());
        cancel.check()?;
        let source_size = std::fs::metadata(&source_bam).map(|m| m.len()).unwrap_or(0);
        // A job that continues an earlier job takes its size from the intermediate file that it
        // starts from, and not from the source.
        //
        // The growth of the source, from CRAM to FASTQ and then back to BAM, already happened, and
        // those files are on the disk. A second count of that space refuses a job that fits.
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
            let params = RevertParams::default();
            log_buffer("collating", params.sort_buffer_bytes);
            let stats = tokio::task::spawn_blocking(move || {
                revert::revert_alignment(&bam, reference.as_deref(), &dir, &params, &token)
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
            // Report the stage, although it has no work. A progress display that goes from stage
            // 3 to stage 5 looks like a step that the code lost. It does not look like a step that
            // the code saved.
            report(RealignStage::Map, resumed.detail());
        }

        // The input of a stage has no more use after the next stage reads it. For a WGS sample,
        // each input is tens of GB.
        //
        // An earlier version kept each input until the end of the job. That version needed about
        // two times the peak space, and a normal disk was then too small.
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
            // A job that continues an earlier job sorts from the start. So the spilled runs of the
            // earlier run have no use. For a WGS sample they are tens of GB, and the sort needs
            // that space.
            clear_stage_dir(&dir);
            let token = cancel.clone();
            let params = SortParams::default();
            log_buffer("sorting", params.buffer_bytes);
            stage(&sorted, async {
                tokio::task::spawn_blocking(move || {
                    postprocess::sort_alignment(&input, &out, &dir, &params, &token, &mut |_| {})
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))?
                .map_err(AppError::from)
            })
            .await?;

            // The merge reads each run file. The next stage reads the sorted BAM file.
            clear_stage_dir(&scratch.path().join("sort"));
            state.completed_through = Some(Resumed::Sorted);
            state.store(scratch.path());
        }

        // The code can remove only the files that *this* run made.
        //
        // A run that continues an earlier run and skips the sort made nothing again. A delete of the
        // mapped BAM file that it started from removes the file that costs the most in this
        // pipeline. The code would remove a file that it did not write, on a belief about that file.
        //
        // That fault occurred. With the marker fault in `discard_partial`, it destroyed a 59 GB
        // `mapped.bam` file, which was 3 hours and 58 minutes of revert work and mapping work.
        if resumed < Resumed::Sorted {
            discard(&mapped);
        }

        report(RealignStage::MarkDuplicates, resumed.detail());
        if resumed < Resumed::Marked {
            let (input, out) = (sorted.clone(), marked.clone());
            let token = cancel.clone();
            // A long-read library usually needs no PCR step, and two long reads rarely have the
            // same end points. So a mark on those reads removes real coverage.
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

        // The output file exists, so the code can remove each intermediate file behind it. A failed
        // registration then costs seconds of work and not hours.
        //
        // The code does not remove `marked` here. The finalize stage *moved* that file into its
        // place, and it made no copy.
        scratch.completed();

        // ---- stage D: register ----
        //
        // This stage is last, by design. The code adds the row only after the file that it names
        // exists. So a failed job, and a job that the user stopped, leave no alignment that names a
        // file with no content.
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
            // A refusal here is the correct result. A map of long reads under a short-read preset
            // does not fail. It gives alignments that look correct and are wrong. So the code
            // reports the refusal to the user.
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
    /// The count of scratch bytes that the job needs.
    pub scratch_needed: u64,
    /// Bytes free where the scratch will live.
    pub scratch_free: u64,
}

/// The factor between the scratch space and the **uncompressed** volume of the source.
///
/// This value comes from a measured run. An estimate from the size of each stage gave a value that
/// was 25% too small.
///
/// A realignment of WGS229, from a CRAM file of 17.3 GB, reached a peak of **276 GB of scratch**.
/// That peak is 16 times the size of the source file. The value `4` here, times the CRAM expansion
/// factor below, gives that result.
///
/// The peak is inside the revert stage, and not at the place that the stage list suggests. The
/// spill runs of that stage use an uncompressed format of this project, and its FASTQ output uses
/// gzip. So the two files exist together at very different densities, and their sum is larger than
/// any later stage.
const SCRATCH_MULTIPLE: u64 = 4;

/// The factor between the size of the data and the size of the file that holds it.
///
/// This correction is important, and a wrong value here is not a small error.
///
/// A CRAM file uses the reference to compress its data. A CRAM file of 17 GB holds about 70 GB of
/// data in the BAM form. The revert step writes FASTQ, and that file is larger again. A FASTQ file
/// uses one ASCII byte for the quality of each base, and a BAM file packs those values.
///
/// An estimate from the size of the *file* tells a user that a job needs 69 GB, when it needs about
/// 200 GB. The disk then fills in the third hour.
fn expansion_factor(source: &Path) -> u64 {
    match source.extension().and_then(|e| e.to_str()) {
        Some(e) if e.eq_ignore_ascii_case("cram") => 4,
        // A BAM file uses bgzf compression, at a factor of about 4. But the revert step writes a
        // FASTQ file with gzip compression. So the ratio between the two files is near 1.
        _ => 2,
    }
}

/// Check that the job can complete, before it starts.
///
/// A disk that fills in the third hour fails the job. It also leaves the machine unusable until
/// somebody finds the scratch directory. The code also needs the RAM value in each case, to size the
/// index.
pub fn preflight(scratch: &Path, source: &Path, source_size: u64) -> Result<RealignPlan, AppError> {
    let needed = source_size
        .saturating_mul(expansion_factor(source))
        .saturating_mul(SCRATCH_MULTIPLE);
    plan_for(scratch, needed, "realign")
}

/// The preflight of a job that continues from intermediate files that exist.
///
/// The full [`preflight`] estimate starts from the *source* file. It then multiplies that size by
/// the growth that the pipeline gives.
///
/// For a job that continues an earlier job, that growth is already on the disk. So the same sum
/// counts it two times, and it refuses a job that has enough space.
///
/// Each remaining stage needs about the size of the intermediate file that it reads, and at most two
/// files exist together. The sort holds its spilled runs while the mapped BAM file is still there.
/// The duplicate mark step writes its output while the sorted BAM file is still there.
///
/// Three times the largest intermediate file covers that need, with a margin. The size of that file
/// is a measurement and not an estimate.
fn resume_preflight(scratch: &Path, mapped: &Path, sorted: &Path, marked: &Path) -> Result<RealignPlan, AppError> {
    let size = |path: &Path| std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let largest = size(mapped).max(size(sorted)).max(size(marked));
    plan_for(scratch, largest.saturating_mul(3), "resume the realignment")
}

/// Measure the disk, refuse a job that it can not hold, and report the decision.
///
/// The two preflight functions differ only in the way that they calculate `needed`. Each later step
/// was the same code two times, and a developer had to change both. Those steps are the read of the
/// free space, the refusal, the text, and the plan.
///
/// The `what` value is the verb of the refusal, so each message stays as it was.
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

/// Record the spill budget that a stage will use.
///
/// The machine now gives both budgets, and no constant holds them. So the version number no longer
/// gives the count of runs that a stage produces. Without this log line, the statement "it spilled
/// 400 runs" in a bug report has no context.
fn log_buffer(stage: &str, bytes: usize) {
    eprintln!("realign: {stage} with a {} MB buffer", bytes / (1024 * 1024));
}

/// Shows whether a job that needs `needed` bytes can start with `free` bytes available.
///
/// This function is separate from the system call, so a test can cover the *decision* on any
/// machine. No test can cover the system call itself.
///
/// A `free` value of 0 means that the platform gave no answer, and the function then permits the
/// job. A job of many hours must not stop because a call for the free space failed. It is better to
/// run that job and let it fail on a real write.
fn has_room(needed: u64, free: u64) -> bool {
    free == 0 || free >= needed
}

/// The count of free bytes on the file system of `path`. The function returns 0 when it can not
/// find that value.
///
/// A zero means "unknown", and the preflight then permits the job. A refusal, because a call for the
/// free space failed, is worse than a job that runs and then fails on a real write.
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

/// The count of free bytes on the volume of `path`.
///
/// The code reads `lpFreeBytesAvailableToCaller`, and not the total free space of the volume. That
/// choice is the important part of this call. On a volume with a disk quota, the two values differ,
/// and the preflight needs the quantity that *this user* can write.
///
/// The Unix code makes the same choice. It reads `f_bavail`, which is the count of blocks for a
/// process with no special rights, and not `f_bfree`.
#[cfg(windows)]
fn fs_free_space(path: &Path) -> u64 {
    use std::os::windows::ffi::OsStrExt;
    use windows_sys::Win32::Storage::FileSystem::GetDiskFreeSpaceExW;

    // `encode_wide` writes no NUL character at the end, and the API needs one.
    let mut wide: Vec<u16> = path.as_os_str().encode_wide().collect();
    wide.push(0);

    let mut available: u64 = 0;
    // SAFETY: `wide` is a correct UTF-16 buffer with a NUL character at the end. This function owns
    // it, and it lives longer than the call. This function also owns `available`, which is a `u64`.
    // The two other totals are null, and the API permits that value.
    let ok = unsafe {
        GetDiskFreeSpaceExW(
            wide.as_ptr(),
            &mut available,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        )
    };

    // A failure gives "unknown" and not a free space of zero. So a call that gave no answer does
    // not look like a full disk, and it does not refuse the job. See `has_room`.
    if ok == 0 {
        return 0;
    }
    available
}

#[cfg(not(any(unix, windows)))]
fn fs_free_space(_path: &Path) -> u64 {
    // This platform has no such call. The preflight reports "unknown" and permits the job. See
    // `free_space`.
    0
}

/// A scratch directory that removes itself.
///
/// The intermediate files are some times the size of the input. After a failed job, and after a job
/// that the user stopped, those files fill a disk and the user sees no message.
///
/// The removal is optional. It must never hide the error that already unwinds the code.
struct JobScratch {
    path: PathBuf,
    /// Set once the job has produced its output, so [`Drop`] can tell "finished" from "died".
    completed: bool,
    /// Keep the intermediates behind when the job does *not* complete.
    keep_on_failure: bool,
}

/// Keep the intermediate files of a failed job. Set `NAVIGATOR_REALIGN_KEEP_SCRATCH=1`.
///
/// The default is off, by design. For a WGS sample this directory holds hundreds of GB. A desktop
/// application must not leave that data after a failure. The user then loses that disk space for a
/// fault that they did not cause.
///
/// The setting exists because the other default has its own fault. A realignment that fails in
/// stage 7 of 8 removes the work of the seven stages that succeeded. On a 30x WGS sample, one
/// measurement gave 10.7 hours of such work. A developer who works on the pipeline needs the
/// opposite behaviour, so this setting gives it.
///
/// This setting also makes [`RealignParams::resume`] useful. A new run reads the intermediate files
/// of an earlier run, and by default a failed run leaves none.
///
/// Use the two together when something can stop a run. A run of six hours, on a machine that
/// somebody also uses, meets that condition often.
///
/// A job that stops at once runs no `Drop` call, so its scratch files stay in each case. The causes
/// are the end of a login session and a power loss. The resume feature exists for that case.
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

    /// The preflight must refuse a job that it can not complete. It must also state the quantity of
    /// space that the job needs, so the user makes no estimate.
    ///
    /// This test runs on each platform where [`fs_free_space`] has real code, and that set now
    /// holds both desktop families. A platform with no such code reports "unknown" and can refuse
    /// nothing. The test `preflight_cannot_refuse_where_free_space_is_unknown` covers that case.
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

    /// The important correction. A CRAM file holds some times its own size in read data. So an
    /// estimate of the scratch space from the size of that file is more than 100 GB too small for a
    /// WGS job.
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

    /// The free space is a check. The check itself must not depend on it. When the platform gives
    /// no answer, the job runs and fails on a real write. It must not stop for a reason that has no
    /// connection to the data.
    #[test]
    fn an_unknown_free_space_does_not_block_the_job() {
        assert!(has_room(u64::MAX, 0), "0 means unknown, not full");
        assert!(has_room(100, 100), "exactly enough is enough");
        assert!(!has_room(101, 100), "one byte short is short");
    }

    /// The scratch directory does not exist while the preflight runs, because the job makes it. So
    /// the call must answer for the file system that will hold that directory. The code moves up the
    /// path to the nearest directory that exists, and it does not stop.
    ///
    /// The test covers both desktop families, for the reason above. The move up the path does not
    /// depend on the platform, but it has a result only where the call at the end can answer.
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

    /// The behaviour of a platform with no call for the free space. This test states that
    /// behaviour, and the absence of a test does not.
    ///
    /// `preflight` can refuse no job that it did not measure. A refusal on an unknown value stops
    /// each realignment on that platform. So the code starts a job that clearly does not fit, and
    /// that job fails on a real write.
    ///
    /// Windows was in this group. It now has `GetDiskFreeSpaceExW`, and the two tests above cover
    /// it.
    #[cfg(not(any(unix, windows)))]
    #[test]
    fn preflight_cannot_refuse_where_free_space_is_unknown() {
        let dir = std::env::temp_dir();
        assert_eq!(free_space(&dir), 0, "no probe on this platform yet");

        let plan = preflight(&dir, Path::new("x.bam"), u64::MAX / 8).expect("unknown must not block");
        assert_eq!(plan.scratch_free, 0);
        assert!(plan.scratch_needed > 0, "the estimate is still made and reported");
    }

    /// The scratch files are some times the size of the input. A job that the user stopped must
    /// remove them.
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

    /// A scratch directory for this test alone. It holds the files that the test names.
    fn scratch(tag: &str, files: &[(&str, bool)]) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("navigator-resume-{tag}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        for (name, complete) in files {
            let path = dir.join(name);
            // A "complete" BAM file holds the BGZF end-of-file marker. An incomplete file holds
            // the same bytes with a short marker, and a writer that stopped leaves such a file.
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

    /// The fault that these two rules exist for. This test reproduces it.
    ///
    /// A user stopped a merge. It left a `sorted.bam` file of 13.2 GB, and the full file is about
    /// 30 GB.
    ///
    /// That partial file held a correct BGZF end-of-file marker. The writer of noodles uses many
    /// threads, and it closes its stream from `Drop` while the code unwinds.
    ///
    /// The marker test alone accepted that file. The next run then marked the duplicates of a short
    /// alignment.
    ///
    /// The stage record refuses the file. The sort never returned, so no code stated that it
    /// completed.
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

    /// The record can also promise no more than the files hold. A scratch directory with no
    /// `marked.bam` file must not start a new run, and an old record that names that file changes
    /// nothing.
    #[test]
    fn the_record_cannot_outrun_the_files() {
        let dir = scratch("stale-record", &[("mapped.bam", true)]);
        let state = ScratchState {
            completed_through: Some(Resumed::Marked),
            ..Default::default()
        };

        assert_eq!(resumable(&dir, &state), Resumed::Mapped);
    }

    /// A scratch directory from before the stage record must still work. A process that stopped at
    /// once left the `mapped.bam` file of 13 August. No `Drop` call ran, so the marker states the
    /// truth.
    #[test]
    fn a_scratch_without_a_record_falls_back_to_the_marker() {
        let dir = scratch("legacy", &[("mapped.bam", true)]);
        assert_eq!(resumable(&dir, &ScratchState::default()), Resumed::Mapped);
    }

    /// The exact case of 2026-08-13. The mapping stage completed, and the sort stopped while it
    /// wrote its output. The mapped BAM file holds four hours of work, and a new run must use it.
    /// That run must not use the partial sorted BAM file.
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

    /// Each stage guard is a test of the form `resumed < Resumed::Sorted`. So the order of these
    /// values controls the code, and it is not only for the reader. A wrong order skips a stage that
    /// did not run.
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

    /// A job that continues an earlier job takes its size from the files on the disk. A size from
    /// the source counts the growth of that source two times, and it refuses a job that fits.
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
