//! The call of a structural variant. It is the port of the Scala `analysis.sv` subsystem, which is
//! a caller of its own, in the style of BreakDancer, Pindel and CNV-seq. It is not a GATK tool.
//!
//! The pipeline is:
//!
//! - [`walker`] collects the evidence in one pass over the BAM. It takes the depth of each bin. It
//!   takes the discordant pairs, from the insert size, from the orientation, and from a pair that
//!   crosses two chromosomes. And it takes the split reads, from the SA tag.
//! - [`segmenter`] turns the depth bins into CNV segments, through a z-score.
//! - [`clusterer`] puts the PE and SR evidence into groups at each breakpoint, infers the SV type,
//!   and brings in the depth segments.
//! - [`caller::call_structural_variants`] drives all of the above.
//!
//! The write of a VCF or an artifact, which `SvVcfWriter` does, waits for later work. The parity
//! target is the Scala caller.

pub mod caller;
pub mod clusterer;
pub mod evidence;
pub mod segmenter;
pub mod types;
pub mod walker;

pub use caller::call_structural_variants;
pub use evidence::{
    BreakpointCluster, DepthSegment, DiscordantPair, DiscordantReason, SplitRead, SvEvidenceCollection,
};
pub use types::{calculate_confidence, SvAnalysisResult, SvCall, SvCallerConfig, SvType};
