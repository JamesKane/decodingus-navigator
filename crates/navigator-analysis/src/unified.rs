//! Unified quality-metrics walker — one coordinate-ordered pass over a BAM/CRAM that
//! collects coverage + callable loci, read-level QC metrics, and sex inference together.
//!
//! The rewrite already had three focused single-pass walkers ([`crate::coverage`],
//! [`crate::read_metrics`], [`crate::sex`]); run separately they read a BAM end-to-end
//! **twice** (coverage pileup + read-metrics scan) and a CRAM **three times** (plus the
//! sex scan, since `.crai` carries no per-reference counts). This walker fuses them into a
//! single record loop — 2→1 for BAM, 3→1 for CRAM (CRAM decode being the expensive case).
//!
//! There is **no metric change**: each record is dispatched to the same `*State` accumulators
//! the standalone walkers use (the single source of truth), so the numbers are byte-for-byte
//! identical to running the three separately. The only subtlety is filtering — coverage
//! pre-filters hard (mapped/primary/main-assembly), but read-metrics needs *every* record and
//! sex needs per-contig mapped tallies, so the loop hands every record to all three states and
//! each applies its own filtering internally.
//!
//! Sex here is tallied directly from the record stream (no BAI dependency), matching the
//! standalone CRAM path's math; the standalone [`crate::sex::infer_from_bam`] keeps its BAI
//! fast path for the cheap à-la-carte "Sex inference" command.

use std::collections::HashSet;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use noodles::core::Region;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::contig;
use crate::coverage::{
    merge_coverage_partials, CallableLociParams, ContigCoverageAccum, ContigCoveragePartial, CoverageResult,
    CoverageState,
};
use crate::error::AnalysisError;
use crate::read_metrics::{ReadMetrics, ReadMetricsState};
use crate::reader;
use crate::sex::{self, SexInferenceResult, SexState};

/// Algorithm version for the unified artifact cache key; bump on any change that alters output.
/// (The three sub-results are persisted under their own existing keys; this is for completeness.)
pub const UNIFIED_VERSION: &str = "unified-1";

/// The three quality-metric results collected in one pass. Sex is `None` when inference
/// can't be computed for the input (no autosomes/chrX, or no autosomal reads — e.g. a
/// targeted panel or chrY-only test); coverage + read-metrics are unaffected, mirroring the
/// pipeline where sex is an independent step whose failure doesn't kill the others.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnifiedMetricsResult {
    pub coverage: CoverageResult,
    pub read_metrics: ReadMetrics,
    pub sex: Option<SexInferenceResult>,
}

/// Single-pass coverage + read-metrics + sex over a coordinate-sorted BAM/CRAM. `reference`
/// is required (CRAM decode + reference-N detection). Equivalent to running
/// [`crate::coverage::collect_coverage_callable`], [`crate::read_metrics::collect_read_metrics`],
/// and tallying sex from the same records — but one read of the file.
pub fn collect_unified_metrics(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<UnifiedMetricsResult, AnalysisError> {
    collect_unified_metrics_with_progress(bam_path, reference_path, params, contig_allowlist, &mut |_, _| {})
}

/// Like [`collect_unified_metrics`], reporting `progress(contigs_done, contigs_total)` as the
/// coverage pass finalizes each tracked contig (the slow whole-genome step — so a progress bar
/// can advance instead of sitting frozen for minutes). Assumes a coordinate-sorted BAM/CRAM.
pub fn collect_unified_metrics_with_progress(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &mut dyn FnMut(usize, usize),
) -> Result<UnifiedMetricsResult, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, Some(reference_path))?;
    let mut cov = CoverageState::new(&header, reference_path, *params, contig_allowlist)?;
    let mut rm = ReadMetricsState::default();
    let mut sx = SexState::new(&header);
    progress(0, cov.total_tracked());

    for result in reader.records(&header) {
        let record = result?;
        // Every record to all three; each state filters internally (see module docs).
        rm.accept(&record);
        sx.accept(&record);
        cov.accept(&record, progress)?;
    }

    let coverage = cov.finish(progress)?;
    let read_metrics = rm.finish();
    // Sex inference is best-effort: `None` (not a hard error) when the input lacks the
    // autosomes/chrX it needs, so coverage + read-metrics still come back.
    let sex = sx.finish().ok();
    Ok(UnifiedMetricsResult { coverage, read_metrics, sex })
}

/// Per-contig parallel unified metrics — the same result as [`collect_unified_metrics`] but
/// computed concurrently across contigs. Coverage is embarrassingly parallel per contig; the
/// per-position pileup compute (not decompression) is the bottleneck a sequential pass hits.
///
/// Requires an **indexed BAM** (per-contig region queries + an unmapped-tail sweep). Anything
/// else — CRAM (no `.crai` unmapped query), or a BAM without a `.bai` — transparently falls
/// back to the sequential [`collect_unified_metrics`], so callers can always prefer this.
///
/// Output is byte-identical to the sequential walker: per contig it runs the same `*State`
/// accumulators, and the merge is over commutative sums / header-ordered per-contig outputs.
/// Read-metrics covers **every** contig (not just main-assembly) plus the unmapped tail — the
/// same record set the sequential pass sees — so totals match exactly.
pub fn collect_unified_metrics_parallel(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<UnifiedMetricsResult, AnalysisError> {
    collect_unified_metrics_parallel_with_progress(bam_path, reference_path, params, contig_allowlist, &|_, _| {})
}

/// Worker threads for the per-contig fan-out. Defaults to all available cores capped at 12 —
/// past that the wall time is floored by the largest contig + the unmapped sweep, so more
/// threads only add memory. Override with `NAVIGATOR_ANALYSIS_THREADS`. Shared with the
/// de-novo caller's region fan-out.
pub(crate) fn analysis_thread_count() -> usize {
    std::env::var("NAVIGATOR_ANALYSIS_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(12))
        .max(1)
}

/// A token from the reference-load semaphore; returns itself to the pool on drop (including on
/// the error path). Bounds how many contigs hold their full reference buffer at once — the peak
/// memory driver, since the per-contig N-mask is tiny once built.
struct LoadPermit<'a> {
    tx: &'a std::sync::mpsc::Sender<()>,
}

impl Drop for LoadPermit<'_> {
    fn drop(&mut self) {
        let _ = self.tx.send(());
    }
}

/// One contig's partial result from the parallel fan-out.
struct ContigPartial {
    rm: ReadMetricsState,
    cov: Option<ContigCoveragePartial>,
    autosome_reads: u64,
    x_reads: u64,
}

/// Like [`collect_unified_metrics_parallel`], reporting `progress(contigs_done, contigs_total)`
/// as each tracked (main-assembly) contig finishes. The progress callback is `Fn + Sync`
/// because it's invoked concurrently from worker threads.
pub fn collect_unified_metrics_parallel_with_progress(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &(dyn Fn(usize, usize) + Sync),
) -> Result<UnifiedMetricsResult, AnalysisError> {
    // The parallel path needs an indexed BAM; everything else takes the sequential walker.
    if !reader::has_bai_index(bam_path) {
        return collect_unified_metrics_with_progress(
            bam_path,
            reference_path,
            params,
            contig_allowlist,
            &mut |d, t| progress(d, t),
        );
    }

    let header = reader::read_header(bam_path, Some(reference_path))?;

    // Work items: one per reference sequence (read-metrics + sex span all contigs). Coverage
    // runs only for tracked = main-assembly ∩ allowlist contigs (matching the sequential walker).
    struct Work {
        ref_id: usize,
        name: String,
        length: usize,
        tracked: bool,
        class: u8, // 0 = other, 1 = autosome, 2 = chrX (for the sex tally)
    }
    let mut works: Vec<Work> = Vec::new();
    let (mut autosome_length, mut x_length) = (0u64, None);
    for (ref_id, (name_bytes, map)) in header.reference_sequences().iter().enumerate() {
        let name = String::from_utf8_lossy(name_bytes.as_ref()).into_owned();
        let length = map.length().get();
        let tracked = contig::is_main_assembly(&name)
            && contig_allowlist.map_or(true, |s| s.contains(&name));
        let class = if contig::is_autosome(&name) {
            autosome_length += length as u64;
            1
        } else if contig::is_chr_x(&name) {
            x_length = Some(length as u64);
            2
        } else {
            0
        };
        works.push(Work { ref_id, name, length, tracked, class });
    }

    let total_cov = works.iter().filter(|w| w.tracked).count();
    let done = AtomicUsize::new(0);
    progress(0, total_cov);

    let n_threads = analysis_thread_count();
    // Bound concurrent full-reference loads (the peak-memory driver) independently of compute
    // parallelism: at most a few contigs hold their raw reference at once while building the
    // compact N-mask. A token pool implements the counting semaphore.
    let load_permits = n_threads.min(4);
    let (perm_tx, perm_rx) = std::sync::mpsc::channel::<()>();
    for _ in 0..load_permits {
        let _ = perm_tx.send(());
    }
    let perm_rx = std::sync::Mutex::new(perm_rx);

    let process_contig = |w: &Work| -> Result<ContigPartial, AnalysisError> {
        let (h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
        let region = Region::new(w.name.as_bytes().to_vec(), ..); // whole contig

        let mut cov_accum = if w.tracked {
            // Hold a load permit only across the raw-reference load + mask build; release before
            // the long pileup (which keeps just the small mask).
            let _permit = {
                let _ = perm_rx.lock().unwrap().recv();
                LoadPermit { tx: &perm_tx }
            };
            let ref_bases = reader::read_contig_sequence(reference_path, &w.name)?;
            Some(ContigCoverageAccum::new(w.name.clone(), w.length, ref_bases, *params))
        } else {
            None
        };
        let mut rm = ReadMetricsState::default();
        let (mut autosome_reads, mut x_reads) = (0u64, 0u64);

        {
            let q = idx.query(&h, &region)?;
            for r in q {
                let record = r?;
                rm.accept(&record);
                if let Some(acc) = cov_accum.as_mut() {
                    acc.accept(&record);
                }
                if w.class != 0 && !record.flags().is_unmapped() {
                    if w.class == 1 {
                        autosome_reads += 1;
                    } else {
                        x_reads += 1;
                    }
                }
            }
        }

        let cov = cov_accum.map(|a| a.finish(w.ref_id));
        if w.tracked {
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            progress(d, total_cov);
        }
        Ok(ContigPartial { rm, cov, autosome_reads, x_reads })
    };

    // The unmapped tail (no reference position) is invisible to region queries but the
    // sequential read-metrics counts it (total/pf reads, read-length) — sweep it separately.
    let process_unmapped = || -> Result<ReadMetricsState, AnalysisError> {
        let (h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
        let mut rm = ReadMetricsState::default();
        {
            let q = idx.query_unmapped(&h)?;
            for r in q {
                rm.accept(&r?);
            }
        }
        Ok(rm)
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|e| AnalysisError::Message(format!("thread pool: {e}")))?;

    let (contig_results, unmapped_rm) = pool.install(|| {
        rayon::join(
            || works.par_iter().map(&process_contig).collect::<Result<Vec<_>, AnalysisError>>(),
            process_unmapped,
        )
    });
    let contig_results = contig_results?;
    let unmapped_rm = unmapped_rm?;

    // Merge: read-metrics is a commutative fold; coverage merges per-contig (header order);
    // sex sums per-contig class counts into one tally.
    let mut rm_total = ReadMetricsState::default();
    let mut cov_partials: Vec<ContigCoveragePartial> = Vec::new();
    let (mut autosome_reads, mut x_reads) = (0u64, 0u64);
    for p in contig_results {
        rm_total.merge(p.rm);
        if let Some(c) = p.cov {
            cov_partials.push(c);
        }
        autosome_reads += p.autosome_reads;
        x_reads += p.x_reads;
    }
    rm_total.merge(unmapped_rm);

    let coverage = merge_coverage_partials(cov_partials);
    let read_metrics = rm_total.finish();
    let sex = sex::result_from_tally((autosome_reads, autosome_length, x_reads, x_length)).ok();
    progress(total_cov, total_cov);
    Ok(UnifiedMetricsResult { coverage, read_metrics, sex })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coverage, read_metrics, sex};
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name)
    }

    /// The fused walker yields exactly what the three standalone walkers do, field for field.
    #[test]
    fn unified_matches_standalone_walkers() {
        let bam = fixture("coverage.bam");
        let reference = fixture("ref.fa");
        let params = CallableLociParams::default();

        let unified = collect_unified_metrics(&bam, &reference, &params, None).unwrap();

        let cov = coverage::collect_coverage_callable(&bam, &reference, &params, None).unwrap();
        let rm = read_metrics::collect_read_metrics(&bam, Some(&reference)).unwrap();
        // The fixture is chrM-only (no autosomes/chrX), so sex inference can't be computed —
        // the fused walker reports it as `None` (best-effort) while still returning coverage +
        // read-metrics, and the standalone walker errors. Both agree sex is unavailable here.
        assert!(sex::infer_from_bam(&bam, Some(&reference)).is_err());

        assert_eq!(unified.coverage, cov, "coverage diverged");
        assert_eq!(unified.read_metrics, rm, "read metrics diverged");
        assert_eq!(unified.sex, None, "expected chrM-only fixture to lack autosomes");
    }
}
