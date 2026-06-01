//! SV core types — port of the Scala `SvTypes` (SvType, SvCall, config, confidence).

/// Structural variant type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
#[derive(Debug, Clone, PartialEq)]
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

/// Confidence in [0,1] weighting PE / SR / depth evidence (mirrors `calculateConfidence`).
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
        }
    }
}

/// Result of SV analysis (timestamp/VCF output handled by the orchestrator/caller).
#[derive(Debug, Clone, PartialEq)]
pub struct SvAnalysisResult {
    pub sv_calls: Vec<SvCall>,
    pub total_discordant_pairs: u64,
    pub total_split_reads: u64,
    pub cnv_segments: usize,
    pub reference_build: String,
    pub mean_coverage: f64,
}
