//! Read-level metrics walker — Rust port of the Scala `UnifiedMetricsWalker`
//! (replaces GATK `CollectAlignmentSummaryMetrics` + `CollectInsertSizeMetrics`,
//! no R dependency). Single pass over the BAM collects alignment-summary counts,
//! read-length and insert-size distributions, pair orientation, and mean MAPQ.
//!
//! Parity target is the Scala walker. Primary metrics exclude secondary/supplementary
//! records; insert size is taken from first-of-pair proper pairs only (no double count).

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;
use crate::reader;
use crate::readview::AlnRead;

const MAX_INSERT_SIZE: i32 = 10_000;

/// Dominant pair orientation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum PairOrientation {
    #[default]
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
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
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

    /// Fold another accumulator in (for the parallel per-contig merge). All fields are
    /// commutative sums / set-min-max / histogram unions, so the merged distribution is
    /// independent of how records were partitioned across contigs.
    fn merge(&mut self, other: DistAccum) {
        for (value, count) in other.hist {
            *self.hist.entry(value).or_insert(0) += count;
        }
        self.sum += other.sum;
        self.sum_sq += other.sum_sq;
        if other.count > 0 {
            if self.count == 0 {
                self.min = other.min;
                self.max = other.max;
            } else {
                self.min = self.min.min(other.min);
                self.max = self.max.max(other.max);
            }
        }
        self.count += other.count;
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
pub fn collect_read_metrics(bam_path: &Path, reference: Option<&Path>) -> Result<ReadMetrics, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, reference)?;
    let mut state = ReadMetricsState::default();
    for result in reader.records_lazy(&header) {
        state.accept(&result?);
    }
    Ok(state.finish())
}

/// Read-level metrics accumulator shared by the standalone walker and the fused
/// [`crate::unified`] walker (one source of truth → byte-identical numbers). Feed every
/// record via [`ReadMetricsState::accept`], then call [`ReadMetricsState::finish`].
#[derive(Default)]
pub(crate) struct ReadMetricsState {
    total_reads: u64,
    pf_reads: u64,
    pf_reads_aligned: u64,
    reads_aligned_in_pairs: u64,
    proper_pairs: u64,
    chimeric_reads: u64,
    total_mapq: u64,
    mapped_for_mq: u64,
    read_len: DistAccum,
    insert: DistAccum,
    fr: u64,
    rf: u64,
    tandem: u64,
}

impl ReadMetricsState {
    /// Fold another state in (parallel per-contig + unmapped-sweep merge). Every field is a
    /// commutative count / histogram union, so the result equals a single sequential pass
    /// regardless of how records were split across contigs.
    pub(crate) fn merge(&mut self, other: ReadMetricsState) {
        self.total_reads += other.total_reads;
        self.pf_reads += other.pf_reads;
        self.pf_reads_aligned += other.pf_reads_aligned;
        self.reads_aligned_in_pairs += other.reads_aligned_in_pairs;
        self.proper_pairs += other.proper_pairs;
        self.chimeric_reads += other.chimeric_reads;
        self.total_mapq += other.total_mapq;
        self.mapped_for_mq += other.mapped_for_mq;
        self.fr += other.fr;
        self.rf += other.rf;
        self.tandem += other.tandem;
        self.read_len.merge(other.read_len);
        self.insert.merge(other.insert);
    }

    pub(crate) fn accept(&mut self, record: &impl AlnRead) {
        let flags = record.flags();

        // Primary metrics exclude secondary/supplementary.
        if flags.is_secondary() || flags.is_supplementary() {
            return;
        }
        self.total_reads += 1;
        if flags.is_qc_fail() {
            return;
        }
        self.pf_reads += 1;

        self.read_len.add(record.sequence_len() as u32);

        if flags.is_unmapped() {
            return;
        }
        self.pf_reads_aligned += 1;

        if let Some(mq) = record.mapping_quality() {
            self.total_mapq += mq as u64; // None == 255 (unavailable), excluded
            self.mapped_for_mq += 1;
        }

        if !flags.is_segmented() {
            return; // not paired
        }
        if !flags.is_mate_unmapped() {
            self.reads_aligned_in_pairs += 1;
        }

        if flags.is_properly_segmented() {
            self.proper_pairs += 1;
            if flags.is_first_segment() {
                let isize = record.template_length().abs();
                if isize > 0 && isize < MAX_INSERT_SIZE {
                    self.insert.add(isize as u32);
                    match detect_orientation(record) {
                        PairOrientation::Fr => self.fr += 1,
                        PairOrientation::Rf => self.rf += 1,
                        PairOrientation::Tandem => self.tandem += 1,
                    }
                }
            }
        }

        if !flags.is_mate_unmapped() && record.reference_sequence_id() != record.mate_reference_sequence_id() {
            self.chimeric_reads += 1;
        }
    }

    pub(crate) fn finish(self) -> ReadMetrics {
        let (median_rl, mean_rl, std_rl) = self.read_len.stats();
        let (median_is, mean_is, std_is) = self.insert.stats();

        let pair_orientation = if self.fr >= self.rf && self.fr >= self.tandem {
            PairOrientation::Fr
        } else if self.rf >= self.tandem {
            PairOrientation::Rf
        } else {
            PairOrientation::Tandem
        };

        ReadMetrics {
            total_reads: self.total_reads,
            pf_reads: self.pf_reads,
            pf_reads_aligned: self.pf_reads_aligned,
            reads_aligned_in_pairs: self.reads_aligned_in_pairs,
            proper_pairs: self.proper_pairs,
            pct_pf_reads_aligned: ratio(self.pf_reads_aligned, self.pf_reads),
            pct_reads_aligned_in_pairs: ratio(self.reads_aligned_in_pairs, self.pf_reads_aligned),
            pct_proper_pairs: ratio(self.proper_pairs, self.pf_reads_aligned),
            median_read_length: median_rl,
            mean_read_length: mean_rl,
            std_read_length: std_rl,
            min_read_length: self.read_len.min,
            max_read_length: self.read_len.max,
            read_length_histogram: self.read_len.hist,
            median_insert_size: median_is,
            mean_insert_size: mean_is,
            std_insert_size: std_is,
            min_insert_size: self.insert.min,
            max_insert_size: self.insert.max,
            insert_size_histogram: self.insert.hist,
            pair_orientation,
            pct_chimeras: ratio(self.chimeric_reads, self.reads_aligned_in_pairs),
            mean_mapping_quality: ratio(self.total_mapq, self.mapped_for_mq),
        }
    }
}

fn ratio(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

/// Pair orientation from a first-of-pair proper pair, mirroring the Scala logic.
fn detect_orientation(record: &impl AlnRead) -> PairOrientation {
    let read_neg = record.flags().is_reverse_complemented();
    let mate_neg = record.flags().is_mate_reverse_complemented();
    if read_neg == mate_neg {
        return PairOrientation::Tandem;
    }
    let start = record.alignment_start().unwrap_or(0);
    let mate_start = record.mate_alignment_start().unwrap_or(0);

    let fr = if start < mate_start {
        !read_neg && mate_neg
    } else {
        read_neg && !mate_neg
    };
    if fr {
        PairOrientation::Fr
    } else {
        PairOrientation::Rf
    }
}
