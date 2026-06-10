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

use serde::{Deserialize, Serialize};

use crate::coverage::{CallableLociParams, CoverageResult, CoverageState};
use crate::error::AnalysisError;
use crate::read_metrics::{ReadMetrics, ReadMetricsState};
use crate::reader;
use crate::sex::{SexInferenceResult, SexState};

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
