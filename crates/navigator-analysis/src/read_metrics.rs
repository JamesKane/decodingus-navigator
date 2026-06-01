//! Read-level metrics walker — Rust port of the Scala `UnifiedMetricsWalker`
//! (replaces GATK `CollectAlignmentSummaryMetrics` + `CollectInsertSizeMetrics`,
//! no R dependency). Single pass over the BAM collects alignment-summary counts,
//! read-length and insert-size distributions, pair orientation, and mean MAPQ.
//!
//! Parity target is the Scala walker. Primary metrics exclude secondary/supplementary
//! records; insert size is taken from first-of-pair proper pairs only (no double count).

use std::collections::BTreeMap;
use std::path::Path;

use noodles::bam;

use crate::error::AnalysisError;

const MAX_INSERT_SIZE: i32 = 10_000;

/// Dominant pair orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairOrientation {
    Fr,
    Rf,
    Tandem,
}

impl PairOrientation {
    pub fn as_str(self) -> &'static str {
        match self {
            PairOrientation::Fr => "FR",
            PairOrientation::Rf => "RF",
            PairOrientation::Tandem => "TANDEM",
        }
    }
}

/// Read-level metrics (mirrors the Scala `ReadMetrics`).
#[derive(Debug, Clone, PartialEq)]
pub struct ReadMetrics {
    pub total_reads: u64,
    pub pf_reads: u64,
    pub pf_reads_aligned: u64,
    pub reads_aligned_in_pairs: u64,
    pub proper_pairs: u64,

    pub pct_pf_reads_aligned: f64,
    pub pct_reads_aligned_in_pairs: f64,
    pub pct_proper_pairs: f64,

    pub median_read_length: f64,
    pub mean_read_length: f64,
    pub std_read_length: f64,
    pub min_read_length: u32,
    pub max_read_length: u32,
    pub read_length_histogram: BTreeMap<u32, u64>,

    pub median_insert_size: f64,
    pub mean_insert_size: f64,
    pub std_insert_size: f64,
    pub min_insert_size: u32,
    pub max_insert_size: u32,
    pub insert_size_histogram: BTreeMap<u32, u64>,

    pub pair_orientation: PairOrientation,

    pub pct_chimeras: f64,
    pub mean_mapping_quality: f64,
}

/// Online accumulators for a distribution summarized by a histogram.
#[derive(Default)]
struct DistAccum {
    hist: BTreeMap<u32, u64>,
    sum: u128,
    sum_sq: u128,
    count: u64,
    min: u32,
    max: u32,
}

impl DistAccum {
    fn add(&mut self, value: u32) {
        *self.hist.entry(value).or_insert(0) += 1;
        self.sum += value as u128;
        self.sum_sq += (value as u128) * (value as u128);
        if self.count == 0 || value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }
        self.count += 1;
    }

    /// (median, mean, std) using population variance, matching the Scala walker.
    fn stats(&self) -> (f64, f64, f64) {
        if self.count == 0 {
            return (0.0, 0.0, 0.0);
        }
        let n = self.count as f64;
        let mean = self.sum as f64 / n;
        let var = (self.sum_sq as f64 / n) - mean * mean;
        let std = var.max(0.0).sqrt();
        (median_from_hist(&self.hist, self.count), mean, std)
    }
}

fn median_from_hist(hist: &BTreeMap<u32, u64>, total: u64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let half = total / 2;
    let mut cumulative = 0u64;
    let mut last = 0.0;
    for (&value, &count) in hist {
        cumulative += count;
        last = value as f64;
        if cumulative >= half {
            return value as f64;
        }
    }
    last
}

/// Single pass over the BAM collecting read-level metrics.
pub fn collect_read_metrics(bam_path: &Path) -> Result<ReadMetrics, AnalysisError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .map_err(|e| AnalysisError::io(bam_path, e))?;
    let _header = reader
        .read_header()
        .map_err(|e| AnalysisError::io(bam_path, e))?;

    let mut total_reads = 0u64;
    let mut pf_reads = 0u64;
    let mut pf_reads_aligned = 0u64;
    let mut reads_aligned_in_pairs = 0u64;
    let mut proper_pairs = 0u64;
    let mut chimeric_reads = 0u64;

    let mut total_mapq = 0u64;
    let mut mapped_for_mq = 0u64;

    let mut read_len = DistAccum::default();
    let mut insert = DistAccum::default();
    let (mut fr, mut rf, mut tandem) = (0u64, 0u64, 0u64);

    for result in reader.records() {
        let record = result.map_err(|e| AnalysisError::io(bam_path, e))?;
        let flags = record.flags();

        // Primary metrics exclude secondary/supplementary.
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        total_reads += 1;
        if flags.is_qc_fail() {
            continue;
        }
        pf_reads += 1;

        read_len.add(record.sequence().len() as u32);

        if flags.is_unmapped() {
            continue;
        }
        pf_reads_aligned += 1;

        if let Some(mq) = record.mapping_quality() {
            total_mapq += mq.get() as u64; // None == 255 (unavailable), excluded
            mapped_for_mq += 1;
        }

        if !flags.is_segmented() {
            continue; // not paired
        }
        if !flags.is_mate_unmapped() {
            reads_aligned_in_pairs += 1;
        }

        if flags.is_properly_segmented() {
            proper_pairs += 1;
            if flags.is_first_segment() {
                let isize = record.template_length().abs();
                if isize > 0 && isize < MAX_INSERT_SIZE {
                    insert.add(isize as u32);
                    match detect_orientation(&record, bam_path)? {
                        PairOrientation::Fr => fr += 1,
                        PairOrientation::Rf => rf += 1,
                        PairOrientation::Tandem => tandem += 1,
                    }
                }
            }
        }

        if !flags.is_mate_unmapped() {
            let r = opt_usize(record.reference_sequence_id(), bam_path)?;
            let m = opt_usize(record.mate_reference_sequence_id(), bam_path)?;
            if r != m {
                chimeric_reads += 1;
            }
        }
    }

    let (median_rl, mean_rl, std_rl) = read_len.stats();
    let (median_is, mean_is, std_is) = insert.stats();

    let pair_orientation = if fr >= rf && fr >= tandem {
        PairOrientation::Fr
    } else if rf >= tandem {
        PairOrientation::Rf
    } else {
        PairOrientation::Tandem
    };

    Ok(ReadMetrics {
        total_reads,
        pf_reads,
        pf_reads_aligned,
        reads_aligned_in_pairs,
        proper_pairs,
        pct_pf_reads_aligned: ratio(pf_reads_aligned, pf_reads),
        pct_reads_aligned_in_pairs: ratio(reads_aligned_in_pairs, pf_reads_aligned),
        pct_proper_pairs: ratio(proper_pairs, pf_reads_aligned),
        median_read_length: median_rl,
        mean_read_length: mean_rl,
        std_read_length: std_rl,
        min_read_length: read_len.min,
        max_read_length: read_len.max,
        read_length_histogram: read_len.hist,
        median_insert_size: median_is,
        mean_insert_size: mean_is,
        std_insert_size: std_is,
        min_insert_size: insert.min,
        max_insert_size: insert.max,
        insert_size_histogram: insert.hist,
        pair_orientation,
        pct_chimeras: ratio(chimeric_reads, reads_aligned_in_pairs),
        mean_mapping_quality: ratio(total_mapq, mapped_for_mq),
    })
}

fn ratio(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

fn opt_usize(
    v: Option<std::io::Result<usize>>,
    path: &Path,
) -> Result<Option<usize>, AnalysisError> {
    match v {
        Some(r) => Ok(Some(r.map_err(|e| AnalysisError::io(path, e))?)),
        None => Ok(None),
    }
}

/// Pair orientation from a first-of-pair proper pair, mirroring the Scala logic.
fn detect_orientation(record: &bam::Record, path: &Path) -> Result<PairOrientation, AnalysisError> {
    let read_neg = record.flags().is_reverse_complemented();
    let mate_neg = record.flags().is_mate_reverse_complemented();
    if read_neg == mate_neg {
        return Ok(PairOrientation::Tandem);
    }
    let start = record
        .alignment_start()
        .transpose()
        .map_err(|e| AnalysisError::io(path, e))?
        .map_or(0, |p| p.get());
    let mate_start = record
        .mate_alignment_start()
        .transpose()
        .map_err(|e| AnalysisError::io(path, e))?
        .map_or(0, |p| p.get());

    let fr = if start < mate_start {
        !read_neg && mate_neg
    } else {
        read_neg && !mate_neg
    };
    Ok(if fr { PairOrientation::Fr } else { PairOrientation::Rf })
}
