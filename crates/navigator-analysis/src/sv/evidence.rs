//! SV evidence models — port of the Scala `SvEvidence` (discordant pairs, split reads,
//! depth segments, the evidence collection, and breakpoint clusters).

use std::collections::BTreeMap;
use std::sync::Arc;

use super::types::SvType;

/// Why a read pair is considered discordant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscordantReason {
    InsertSizeOutlier,
    WrongOrientation,
    InterChromosomal,
}

/// A discordant read pair (potential SV breakpoint evidence).
///
/// Contig names are `Arc<str>` rather than `String`, and there is no read name. A 30x WGS retains
/// 3–16 M of these (measured across the workspace), and both choices are about that scale: the
/// walker interns one `Arc` per contig and clones a pointer instead of allocating a name per
/// record, and the read name — which nothing downstream ever read — cost an allocation and ~55
/// bytes each to carry evidence that clustering identifies purely by position.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscordantPair {
    pub chrom1: Arc<str>,
    pub pos1: i64,
    pub strand1: char,
    pub chrom2: Arc<str>,
    pub pos2: i64,
    pub strand2: char,
    pub insert_size: i32,
    pub mapq: u8,
    pub reason: DiscordantReason,
}

/// A split read (alignment split across two loci, from the SA tag). Interned contig names and no
/// read name, for the reasons in [`DiscordantPair`].
#[derive(Debug, Clone, PartialEq)]
pub struct SplitRead {
    pub primary_chrom: Arc<str>,
    pub primary_pos: i64,
    pub primary_strand: char,
    pub supp_chrom: Arc<str>,
    pub supp_pos: i64,
    pub supp_strand: char,
    pub clip_length: i32,
    pub mapq: u8,
}

/// A genome segment with abnormal copy number.
#[derive(Debug, Clone, PartialEq)]
pub struct DepthSegment {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub mean_depth: f64,
    pub log2_ratio: f64,
    pub z_score: f64,
    pub num_bins: u32,
    pub sv_type: SvType,
}

/// All SV evidence gathered from a BAM. `depth_bins` maps contig -> per-bin read counts.
#[derive(Debug, Clone, PartialEq)]
pub struct SvEvidenceCollection {
    pub discordant_pairs: Vec<DiscordantPair>,
    pub split_reads: Vec<SplitRead>,
    pub depth_bins: BTreeMap<String, Vec<u32>>,
    pub sample_name: String,
    pub expected_insert_size: f64,
    pub insert_size_sd: f64,
    /// Evidence seen but not retained, because `SvCallerConfig::max_evidence_records` was already
    /// met. Zero in every normal run — see that field. Kept so the `total_*` counts below stay the
    /// number of items *found*, which is what the walker's stats mean, capped or not.
    pub discordant_pairs_dropped: u64,
    pub split_reads_dropped: u64,
}

impl SvEvidenceCollection {
    pub fn total_discordant_pairs(&self) -> u64 {
        self.discordant_pairs.len() as u64 + self.discordant_pairs_dropped
    }

    pub fn total_split_reads(&self) -> u64 {
        self.split_reads.len() as u64 + self.split_reads_dropped
    }

    pub fn inter_chromosomal_pairs(&self) -> Vec<DiscordantPair> {
        self.discordant_pairs
            .iter()
            .filter(|p| p.reason == DiscordantReason::InterChromosomal)
            .cloned()
            .collect()
    }
}

/// Grouped evidence supporting a single breakpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct BreakpointCluster {
    pub chrom: String,
    pub position: i64,
    pub ci_low: i32,
    pub ci_high: i32,
    pub discordant_pairs: Vec<DiscordantPair>,
    pub split_reads: Vec<SplitRead>,
    pub mate_chrom: Option<String>,
    pub mate_position: Option<i64>,
}

impl BreakpointCluster {
    pub fn total_support(&self) -> u32 {
        (self.discordant_pairs.len() + self.split_reads.len()) as u32
    }

    pub fn pe_support(&self) -> u32 {
        self.discordant_pairs.len() as u32
    }

    pub fn sr_support(&self) -> u32 {
        self.split_reads.len() as u32
    }

    pub fn mean_mapq(&self) -> f64 {
        let n = self.discordant_pairs.len() + self.split_reads.len();
        if n == 0 {
            return 0.0;
        }
        let sum: u64 = self.discordant_pairs.iter().map(|p| p.mapq as u64).sum::<u64>()
            + self.split_reads.iter().map(|s| s.mapq as u64).sum::<u64>();
        sum as f64 / n as f64
    }
}
