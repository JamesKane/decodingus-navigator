//! SV caller orchestration — port of the Scala `SvCaller`: evidence collection ->
//! depth segmentation -> clustering -> result. Artifact/VCF writing (`SvVcfWriter`) is
//! deferred, like the coverage walker's BED output.

use std::collections::BTreeMap;
use std::path::Path;

use super::types::{SvAnalysisResult, SvCallerConfig};
use super::{clusterer, segmenter, walker};
use crate::error::AnalysisError;

const MERGE_MAX_GAP: i64 = 50_000;

/// Run the full SV pipeline on a BAM. Requires >= 10x mean coverage (Scala threshold).
#[allow(clippy::too_many_arguments)]
pub fn call_structural_variants(
    bam_path: &Path,
    contig_lengths: &BTreeMap<String, i64>,
    reference_build: &str,
    mean_coverage: f64,
    mean_insert_size: f64,
    insert_size_sd: f64,
    mean_read_length: f64,
    config: &SvCallerConfig,
) -> Result<SvAnalysisResult, AnalysisError> {
    if mean_coverage < 10.0 {
        return Err(AnalysisError::Message(format!(
            "coverage too low for SV calling ({mean_coverage}x, minimum 10x required)"
        )));
    }

    let evidence = walker::collect_evidence(bam_path, contig_lengths, mean_insert_size, insert_size_sd, config)?;

    let raw_segments = segmenter::segment(
        &evidence.depth_bins,
        contig_lengths,
        mean_coverage,
        mean_read_length,
        config,
    );
    let merged_segments = segmenter::merge_nearby_segments(&raw_segments, MERGE_MAX_GAP);

    let sv_calls = clusterer::cluster(&evidence, &merged_segments, config);

    Ok(SvAnalysisResult {
        sv_calls: sv_calls.into_iter().filter(|c| c.filter == "PASS").collect(),
        total_discordant_pairs: evidence.total_discordant_pairs(),
        total_split_reads: evidence.total_split_reads(),
        cnv_segments: merged_segments.len(),
        reference_build: reference_build.to_string(),
        mean_coverage,
    })
}
