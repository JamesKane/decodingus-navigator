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
//! preflight ──► revert ──► index ──► map ──► sort ──► markdup ──► CRAM ──► register
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
    Compress,
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
        RealignStage::Compress,
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
            RealignStage::Compress => "Compressing to CRAM",
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
#[derive(Debug, Clone)]
pub struct RealignOutcome {
    pub alignment: Alignment,
    /// Reads that were unmapped in the source and therefore had a chance to be placed by the new
    /// reference. This is the number the module exists to move, so it is surfaced rather than
    /// buried in a log.
    pub source_unmapped_reads: u64,
    pub reads_written: u64,
    pub duplicates_marked: u64,
    /// Secondary/supplementary records the mapper emitted without `SEQ`, which CRAM cannot encode.
    /// No read is lost — the bases live on the primary — but the count is reported rather than
    /// hidden, like every other drop in this pipeline.
    pub sequenceless_dropped: u64,
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

        // ---- preflight ----
        report(RealignStage::Preflight, "");
        cancel.check()?;
        let source_size = std::fs::metadata(&source_bam).map(|m| m.len()).unwrap_or(0);
        let plan = preflight(scratch.path(), Path::new(&source_bam), source_size)?;

        // ---- stage A: revert ----
        report(RealignStage::Revert, "");
        let reverted = {
            let bam = PathBuf::from(&source_bam);
            let reference = source_reference.clone();
            let dir = scratch.path().join("revert");
            let token = cancel.clone();
            tokio::task::spawn_blocking(move || {
                revert::revert_alignment(&bam, reference.as_deref(), &dir, &RevertParams::default(), &token)
            })
            .await
            .map_err(|e| AppError::Join(e.to_string()))??
        };

        // ---- stage B: index, then map ----
        report(RealignStage::Index, "");
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
        let mapped = scratch.path().join("mapped.bam");
        {
            let (r1, r2, singles) = (
                reverted.read1.clone(),
                reverted.read2.clone(),
                reverted.singletons.clone(),
            );
            let out = mapped.clone();
            let dir = scratch.path().join("map");
            let token = cancel.clone();
            let map_params = navigator_align::MapParams {
                preset,
                threads: 0,
                read_group: None,
                format: navigator_align::OutputFormat::Bam,
                reference: None,
            };
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
            .map_err(|e| AppError::Join(e.to_string()))??;
        }

        // Every stage's input is dead once the next stage has read it, and at WGS scale each is
        // tens of GB. Holding them all until the job ends — which is what this did first — roughly
        // doubles the peak and is the difference between fitting on a normal disk and not.
        let discard = |path: &Path| {
            let _ = std::fs::remove_file(path);
        };
        discard(&reverted.read1);
        discard(&reverted.read2);
        discard(&reverted.singletons);

        // ---- stage C: sort, mark duplicates, compress ----
        report(RealignStage::Sort, "");
        let sorted = scratch.path().join("sorted.bam");
        {
            let (input, out) = (mapped.clone(), sorted.clone());
            let dir = scratch.path().join("sort");
            let token = cancel.clone();
            tokio::task::spawn_blocking(move || {
                postprocess::sort_alignment(&input, &out, &dir, &SortParams::default(), &token, &mut |_| {})
            })
            .await
            .map_err(|e| AppError::Join(e.to_string()))??
        };

        discard(&mapped);

        report(RealignStage::MarkDuplicates, "");
        let marked = scratch.path().join("marked.bam");
        let markdup = {
            let (input, out) = (sorted.clone(), marked.clone());
            let token = cancel.clone();
            // Long-read libraries are typically PCR-free and long reads rarely share endpoints, so
            // marking them would discard real coverage.
            let md_params = MarkDupParams {
                enabled: preset.is_paired(),
                ..Default::default()
            };
            tokio::task::spawn_blocking(move || {
                postprocess::mark_duplicates(&input, &out, &md_params, &token, &mut |_| {})
            })
            .await
            .map_err(|e| AppError::Join(e.to_string()))??
        };

        discard(&sorted);

        report(RealignStage::Compress, "");
        let cram = {
            let (input, out) = (marked.clone(), output.clone());
            let reference = params.target_reference.clone();
            let token = cancel.clone();
            tokio::task::spawn_blocking(move || postprocess::write_cram(&input, &out, &reference, &token, &mut |_| {}))
                .await
                .map_err(|e| AppError::Join(e.to_string()))??
        };

        // The CRAM exists, so the intermediates behind it are disposable even if registration then
        // fails — re-registering is seconds of work, not hours.
        scratch.completed();

        // ---- stage D: register ----
        //
        // Last, deliberately. The row is only inserted once the artifact it points at exists, so a
        // job that fails or is cancelled leaves no alignment referring to a file that was never
        // finished.
        discard(&marked);

        report(RealignStage::Register, "");
        cancel.check()?;
        let alignment = self
            .register_realigned_alignment(
                source_id,
                &cram.cram,
                &params.target_reference,
                &params.target_build,
                "minimap2",
                preset.as_str(),
            )
            .await?;

        Ok(RealignOutcome {
            alignment,
            source_unmapped_reads: reverted.stats.unmapped_reads,
            reads_written: cram.records,
            duplicates_marked: markdup.duplicates,
            sequenceless_dropped: cram.sequenceless_dropped,
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
    let free = free_space(scratch);

    if !has_room(needed, free) {
        return Err(AppError::Import(format!(
            "not enough room to realign: about {} GB of working space is needed and {} GB is free \
             on {}",
            gb(needed),
            gb(free),
            scratch.display(),
        )));
    }

    Ok(RealignPlan {
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

#[cfg(not(unix))]
fn fs_free_space(_path: &Path) -> u64 {
    // Windows free-space needs `GetDiskFreeSpaceExW`. Until that is wired, preflight reports
    // "unknown" and declines to block — see `free_space`.
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
/// dies in stage 7 of 8 discards the seven stages that worked — a measured 10.7 hours on a 30x WGS
/// — and there is no way to resume, because the input each stage needed has been deleted. For
/// anyone iterating on the pipeline, that trade is the wrong way round, so they can invert it.
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
    #[test]
    fn free_space_resolves_through_a_directory_that_does_not_exist_yet() {
        let unborn = std::env::temp_dir().join("dun-not-created-yet").join("nor-this");
        assert!(!unborn.exists());
        assert!(
            free_space(&unborn) > 0,
            "must report the filesystem that will hold the scratch"
        );
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
