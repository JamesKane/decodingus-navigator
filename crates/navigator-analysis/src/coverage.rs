//! Coverage + callable-loci walker — the Rust port of the Scala
//! `CoverageCallableWalker` (which itself replaces GATK `CollectWgsMetrics` +
//! `CallableLoci` over htsjdk). Single pass over a coordinate-sorted BAM/CRAM builds,
//! per main-assembly contig: a depth histogram, per-position callable state, and
//! samtools-style coverage stats.
//!
//! Parity target is the Scala walker, not samtools — notably mean base/mapping quality
//! are averaged **per base observation** (Σ quality / Σ depth), where samtools averages
//! per read. See the crate fixture tests for hand-computed expected values.
//!
//! Memory: this slice accumulates dense per-position arrays for each tracked contig
//! (fine for mtDNA/Y or any `contig_allowlist`-bounded run; a streaming pileup is the
//! whole-genome optimization). BED-interval output and progress callbacks from the
//! Scala walker are deferred.

use std::collections::HashSet;
use std::path::Path;

use noodles::bam;
use noodles::core::Region;
use noodles::fasta;

use serde::{Deserialize, Serialize};

use crate::contig;
use crate::error::AnalysisError;

/// Algorithm version for the coverage artifact cache key; bump on any change that
/// alters output (plan §6 cache versioning).
pub const COVERAGE_VERSION: &str = "coverage-1";

/// Callable-loci parameters. Defaults match GATK `CallableLoci` (and the Scala walker).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallableLociParams {
    pub min_depth: u32,
    pub max_depth: Option<u32>,
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub max_low_mapq: u8,
    pub max_fraction_low_mapq: f64,
}

impl Default for CallableLociParams {
    fn default() -> Self {
        CallableLociParams {
            min_depth: 4,
            max_depth: None,
            min_mapping_quality: 10,
            min_base_quality: 20,
            max_low_mapq: 1,
            max_fraction_low_mapq: 0.1,
        }
    }
}

/// Per-position callable classification (GATK `CallableLoci` states). Hierarchical:
/// the first failing condition wins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CallableState {
    RefN,
    NoCoverage,
    PoorMappingQuality,
    LowCoverage,
    ExcessiveCoverage,
    Callable,
}

/// Per-contig callable-state base counts.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContigCallableMetrics {
    pub contig: String,
    pub ref_n: u64,
    pub callable: u64,
    pub no_coverage: u64,
    pub low_coverage: u64,
    pub excessive_coverage: u64,
    pub poor_mapping_quality: u64,
}

/// Per-contig samtools-`coverage`-style stats (averaged per base observation; see module docs).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContigCoverageStats {
    pub contig: String,
    pub start_pos: u64, // always 1
    pub end_pos: u64,   // contig length
    pub num_reads: u64,
    pub cov_bases: u64,
    pub coverage: f64, // percent of contig with depth > 0
    pub mean_depth: f64,
    pub mean_base_q: f64,
    pub mean_map_q: f64,
}

/// Combined coverage + callable result (replaces the Scala `CoverageCallableResult`'s
/// global-metrics + callable-summary + samtools-stats fields).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CoverageResult {
    pub genome_territory: u64,
    pub mean_coverage: f64,
    pub median_coverage: f64,
    pub sd_coverage: f64,
    /// Depth histogram, clamped at index 255.
    pub coverage_histogram: Vec<u64>,
    pub pct_1x: f64,
    pub pct_5x: f64,
    pub pct_10x: f64,
    pub pct_15x: f64,
    pub pct_20x: f64,
    pub pct_25x: f64,
    pub pct_30x: f64,
    pub pct_40x: f64,
    pub pct_50x: f64,
    pub callable_bases: u64,
    pub contig_callable: Vec<ContigCallableMetrics>,
    pub contig_coverage_stats: Vec<ContigCoverageStats>,
}

const HIST_LEN: usize = 256;

/// Dense per-position accumulators for one contig.
struct ContigAccum {
    name: String,
    length: usize,
    depth: Vec<u32>,
    base_q_sum: Vec<u64>,
    map_q_sum: Vec<u64>,
    qc_pass: Vec<u32>,
    low_mapq: Vec<u32>,
    read_count: u64,
}

impl ContigAccum {
    fn new(name: String, length: usize) -> Self {
        ContigAccum {
            name,
            length,
            depth: vec![0; length],
            base_q_sum: vec![0; length],
            map_q_sum: vec![0; length],
            qc_pass: vec![0; length],
            low_mapq: vec![0; length],
            read_count: 0,
        }
    }
}

/// Walk a coordinate-sorted BAM and collect coverage + callable metrics for the
/// main-assembly contigs (optionally restricted to `contig_allowlist`).
pub fn collect_coverage_callable(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<CoverageResult, AnalysisError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .map_err(|e| AnalysisError::io(bam_path, e))?;
    let header = reader
        .read_header()
        .map_err(|e| AnalysisError::io(bam_path, e))?;

    // Resolve which reference ids we track, and their names/lengths (header order).
    let ref_seqs = header.reference_sequences();
    let mut tracked: Vec<Option<ContigAccum>> = Vec::with_capacity(ref_seqs.len());
    for (name_bytes, map) in ref_seqs.iter() {
        let name = String::from_utf8_lossy(name_bytes.as_ref()).into_owned();
        let keep = contig::is_main_assembly(&name)
            && contig_allowlist.map_or(true, |set| set.contains(&name));
        tracked.push(keep.then(|| ContigAccum::new(name, map.length().get())));
    }

    // Single pass over records.
    for result in reader.records() {
        let record = result.map_err(|e| AnalysisError::io(bam_path, e))?;
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        let ref_id = match record.reference_sequence_id() {
            Some(r) => r.map_err(|e| AnalysisError::io(bam_path, e))?,
            None => continue,
        };
        let accum = match tracked.get_mut(ref_id).and_then(|o| o.as_mut()) {
            Some(a) => a,
            None => continue,
        };

        let start = match record.alignment_start() {
            Some(p) => p.map_err(|e| AnalysisError::io(bam_path, e))?.get(),
            None => continue,
        };
        let mapq = record.mapping_quality().map_or(255u8, |m| m.get());
        let quals = record.quality_scores();
        let quals = quals.as_ref();

        accum.read_count += 1;

        // Walk the CIGAR over reference positions; only M/=/X contribute a base.
        let mut ref_pos = start; // 1-based
        let mut query_off = 0usize;
        for op in record.cigar().iter() {
            let op = op.map_err(|e| AnalysisError::io(bam_path, e))?;
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let pos = ref_pos + i; // 1-based
                        if pos >= 1 && pos <= accum.length {
                            let idx = pos - 1;
                            let base_q = quals.get(query_off + i).copied().unwrap_or(0);
                            accum.depth[idx] += 1;
                            accum.base_q_sum[idx] += base_q as u64;
                            accum.map_q_sum[idx] += mapq as u64;
                            if mapq >= params.min_mapping_quality
                                && base_q >= params.min_base_quality
                            {
                                accum.qc_pass[idx] += 1;
                            }
                            if mapq <= params.max_low_mapq {
                                accum.low_mapq[idx] += 1;
                            }
                        }
                    }
                    ref_pos += len;
                    query_off += len;
                }
                (true, false) => ref_pos += len,  // D / N
                (false, true) => query_off += len, // I / S
                (false, false) => {}               // H / P
            }
        }
    }

    finalize(reference_path, params, tracked)
}

fn finalize(
    reference_path: &Path,
    params: &CallableLociParams,
    tracked: Vec<Option<ContigAccum>>,
) -> Result<CoverageResult, AnalysisError> {
    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)
        .map_err(|e| AnalysisError::io(reference_path, e))?;

    let mut global_hist = vec![0u64; HIST_LEN];
    let mut global_n: u64 = 0;
    let mut global_sum_depth: u128 = 0;
    let mut global_sum_sq: u128 = 0;

    let mut contig_callable = Vec::new();
    let mut contig_stats = Vec::new();

    for accum in tracked.into_iter().flatten() {
        // Load the contig reference bases (uppercased) for N detection.
        let region: Region = format!("{}", accum.name)
            .parse()
            .map_err(|_| AnalysisError::Message(format!("bad region for contig {}", accum.name)))?;
        let ref_record = fasta_reader
            .query(&region)
            .map_err(|e| AnalysisError::io(reference_path, e))?;
        let ref_bases = ref_record.sequence().as_ref();

        let mut cm = ContigCallableMetrics {
            contig: accum.name.clone(),
            ref_n: 0,
            callable: 0,
            no_coverage: 0,
            low_coverage: 0,
            excessive_coverage: 0,
            poor_mapping_quality: 0,
        };
        let mut covered: u64 = 0;
        let mut total_base_obs: u64 = 0;
        let mut base_q_total: u64 = 0;
        let mut map_q_total: u64 = 0;
        let mut contig_sum_depth: u128 = 0;

        for idx in 0..accum.length {
            let depth = accum.depth[idx];
            let clamped = depth.min(255) as usize;
            global_hist[clamped] += 1;
            global_n += 1;
            global_sum_depth += depth as u128;
            global_sum_sq += (depth as u128) * (depth as u128);
            contig_sum_depth += depth as u128;

            if depth > 0 {
                covered += 1;
                total_base_obs += depth as u64;
                base_q_total += accum.base_q_sum[idx];
                map_q_total += accum.map_q_sum[idx];
            }

            let ref_base = ref_bases.get(idx).copied().unwrap_or(b'N');
            match determine_state(ref_base, depth, accum.qc_pass[idx], accum.low_mapq[idx], params) {
                CallableState::RefN => cm.ref_n += 1,
                CallableState::NoCoverage => cm.no_coverage += 1,
                CallableState::PoorMappingQuality => cm.poor_mapping_quality += 1,
                CallableState::LowCoverage => cm.low_coverage += 1,
                CallableState::ExcessiveCoverage => cm.excessive_coverage += 1,
                CallableState::Callable => cm.callable += 1,
            }
        }

        let length = accum.length as f64;
        contig_stats.push(ContigCoverageStats {
            contig: accum.name.clone(),
            start_pos: 1,
            end_pos: accum.length as u64,
            num_reads: accum.read_count,
            cov_bases: covered,
            coverage: if accum.length == 0 { 0.0 } else { covered as f64 / length * 100.0 },
            mean_depth: if accum.length == 0 { 0.0 } else { contig_sum_depth as f64 / length },
            mean_base_q: if total_base_obs == 0 { 0.0 } else { base_q_total as f64 / total_base_obs as f64 },
            mean_map_q: if total_base_obs == 0 { 0.0 } else { map_q_total as f64 / total_base_obs as f64 },
        });
        contig_callable.push(cm);
    }

    let mean = if global_n == 0 { 0.0 } else { global_sum_depth as f64 / global_n as f64 };
    let sd = if global_n < 2 {
        0.0
    } else {
        let var = global_sum_sq as f64 / global_n as f64 - mean * mean;
        var.max(0.0).sqrt()
    };
    let callable_bases = contig_callable.iter().map(|c| c.callable).sum();

    Ok(CoverageResult {
        genome_territory: global_n,
        mean_coverage: mean,
        median_coverage: median_from_hist(&global_hist, global_n),
        sd_coverage: sd,
        pct_1x: pct_at_least(&global_hist, global_n, 1),
        pct_5x: pct_at_least(&global_hist, global_n, 5),
        pct_10x: pct_at_least(&global_hist, global_n, 10),
        pct_15x: pct_at_least(&global_hist, global_n, 15),
        pct_20x: pct_at_least(&global_hist, global_n, 20),
        pct_25x: pct_at_least(&global_hist, global_n, 25),
        pct_30x: pct_at_least(&global_hist, global_n, 30),
        pct_40x: pct_at_least(&global_hist, global_n, 40),
        pct_50x: pct_at_least(&global_hist, global_n, 50),
        coverage_histogram: global_hist,
        callable_bases,
        contig_callable,
        contig_coverage_stats: contig_stats,
    })
}

/// GATK `CallableLoci` hierarchy — first failing condition wins. Mirrors the Scala
/// `determineCallableState`.
fn determine_state(
    ref_base: u8,
    depth: u32,
    qc_pass: u32,
    low_mapq: u32,
    params: &CallableLociParams,
) -> CallableState {
    if ref_base == b'N' || ref_base == b'n' {
        return CallableState::RefN;
    }
    if depth == 0 {
        return CallableState::NoCoverage;
    }
    let low_frac = low_mapq as f64 / depth as f64;
    if low_frac > params.max_fraction_low_mapq {
        return CallableState::PoorMappingQuality;
    }
    if qc_pass < params.min_depth {
        return CallableState::LowCoverage;
    }
    if params.max_depth.is_some_and(|m| qc_pass > m) {
        return CallableState::ExcessiveCoverage;
    }
    CallableState::Callable
}

fn pct_at_least(hist: &[u64], total: u64, min_depth: usize) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let at_least: u64 = hist[min_depth..].iter().sum();
    at_least as f64 / total as f64
}

fn median_from_hist(hist: &[u64], total: u64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let half = total / 2;
    let mut cumulative = 0u64;
    for (depth, &count) in hist.iter().enumerate() {
        cumulative += count;
        if cumulative >= half {
            return depth as f64;
        }
    }
    255.0
}

