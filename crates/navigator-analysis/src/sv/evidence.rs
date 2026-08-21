//! The models of the SV evidence. This is the port of the Scala `SvEvidence`. It holds the
//! discordant pairs, the split reads, the depth segments, the collection of the evidence, and the
//! clusters of breakpoints.

use std::collections::BTreeMap;
use std::sync::Arc;

use super::types::SvType;

/// The reason that a read pair counts as discordant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscordantReason {
    InsertSizeOutlier,
    WrongOrientation,
    InterChromosomal,
}

/// A discordant read pair (potential SV breakpoint evidence).
///
/// A contig name is an `Arc<str>`, and not a `String`. There is also no read name. A 30x WGS keeps
/// 3M to 16M of these, as a measurement over the workspace showed. Both choices come from that
/// scale.
///
/// The walker interns one `Arc` for each contig, and it then clones a pointer. It does not allocate
/// a name at each record.
///
/// The read name cost an allocation, and about 55 bytes, at each record. Nothing after the walker
/// ever read it. The step that groups the evidence finds a breakpoint from the positions alone.
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

/// All of the SV evidence that the walker collected from a BAM. `depth_bins` maps a contig to the
/// read count of each bin.
#[derive(Debug, Clone, PartialEq)]
pub struct SvEvidenceCollection {
    pub discordant_pairs: Vec<DiscordantPair>,
    pub split_reads: Vec<SplitRead>,
    pub depth_bins: BTreeMap<String, Vec<u32>>,
    pub sample_name: String,
    pub expected_insert_size: f64,
    pub insert_size_sd: f64,
    /// The count of evidence items that the walker saw and did not keep, because the count had
    /// already reached `SvCallerConfig::max_evidence_records`. It is zero in every usual run. See
    /// that field.
    ///
    /// This field exists so that the `total_*` counts below stay the count of items that the walker
    /// *found*. That is what the statistics of the walker mean, whether the cap fired or not.
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

/// The evidence behind one breakpoint, in a group.
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
