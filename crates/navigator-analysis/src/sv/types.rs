//! The core SV types. This is the port of the Scala `SvTypes`, which held SvType, SvCall, the
//! configuration and the confidence.

use serde::{Deserialize, Serialize};

/// Structural variant type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SvType {
    Del, // deletion
    Dup, // duplication
    Inv, // inversion
    Bnd, // breakend (translocation)
    Ins, // insertion
}

impl SvType {
    pub fn as_str(self) -> &'static str {
        match self {
            SvType::Del => "DEL",
            SvType::Dup => "DUP",
            SvType::Inv => "INV",
            SvType::Bnd => "BND",
            SvType::Ins => "INS",
        }
    }
}

/// A called structural variant.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SvCall {
    pub id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub sv_type: SvType,
    pub sv_len: i64,
    pub ci_pos: (i32, i32),
    pub ci_end: (i32, i32),
    pub quality: f64,
    pub paired_end_support: u32,
    pub split_read_support: u32,
    pub relative_depth: Option<f64>,
    pub mate_chrom: Option<String>,
    pub mate_pos: Option<i64>,
    pub filter: String,
    pub genotype: String,
}

/// The confidence, in [0,1]. It weighs the PE evidence, the SR evidence and the depth evidence. It
/// has the same shape as `calculateConfidence`.
pub fn calculate_confidence(call: &SvCall) -> f64 {
    let pe_weight = 0.3;
    let sr_weight = 0.4;
    let depth_weight = 0.3;

    let pe_score = (call.paired_end_support as f64 / 10.0).min(1.0);
    let sr_score = (call.split_read_support as f64 / 5.0).min(1.0);
    let depth_score = call.relative_depth.map_or(0.0, |rd| {
        let deviation = (1.0 - rd).abs();
        (deviation / 0.5).min(1.0)
    });

    pe_score * pe_weight + sr_score * sr_weight + depth_score * depth_weight
}

/// SV-calling configuration. Defaults match the Scala `SvCallerConfig`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SvCallerConfig {
    pub bin_size: i64,
    pub min_depth_z_score: f64,
    pub min_cnv_size: i64,
    pub insert_size_z_threshold: f64,
    pub min_mapq: u8,
    pub max_cluster_distance: i64,
    pub min_paired_end_support: u32,
    pub min_split_read_support: u32,
    pub min_total_support: u32,
    pub min_quality: f64,
    /// The upper limit on the discordant pairs that one walk keeps, and, apart from that, on the
    /// split reads that it keeps.
    ///
    /// It is a safety valve, and not a filter. One library of an unusual shape must not take the
    /// whole process down with an out-of-memory error. In a batch, that would take the other 147
    /// samples with it. See the default.
    pub max_evidence_records: u64,
}

impl Default for SvCallerConfig {
    fn default() -> Self {
        SvCallerConfig {
            bin_size: 1000,
            min_depth_z_score: 2.5,
            min_cnv_size: 10_000,
            insert_size_z_threshold: 4.0,
            min_mapq: 20,
            max_cluster_distance: 500,
            min_paired_end_support: 2,
            min_split_read_support: 1,
            min_total_support: 3,
            min_quality: 10.0,
            // This value comes from the real spread, and it is not a round number. Across 33
            // alignments that somebody analysed, the discordant-pair counts run from 0.05M to
            // 15.6M, and the split reads from 0 to 1.1M. So 32M is about 2x the worst case that
            // anybody has seen, and it fires on none of them.
            //
            // What it does give is an upper limit of about 2 GB on evidence that the size of the
            // file alone used to bound.
            max_evidence_records: 32_000_000,
        }
    }
}

/// Result of SV analysis (timestamp/VCF output handled by the orchestrator/caller).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SvAnalysisResult {
    pub sv_calls: Vec<SvCall>,
    pub total_discordant_pairs: u64,
    pub total_split_reads: u64,
    pub cnv_segments: usize,
    pub reference_build: String,
    pub mean_coverage: f64,
}
