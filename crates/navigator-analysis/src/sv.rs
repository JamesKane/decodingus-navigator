//! Structural variant calling — port of the Scala `analysis.sv` subsystem (a custom
//! BreakDancer/Pindel/CNV-seq-style caller, not a GATK tool). Pipeline:
//!
//! - [`walker`] gathers evidence in one BAM pass: per-bin depth, discordant pairs
//!   (insert-size / orientation / inter-chromosomal), and SA-tag split reads.
//! - [`segmenter`] turns depth bins into CNV segments via z-score analysis.
//! - [`clusterer`] groups PE/SR evidence into breakpoints, infers SV type, and
//!   integrates depth segments.
//! - [`caller::call_structural_variants`] orchestrates the above.
//!
//! VCF/artifact output (`SvVcfWriter`) is deferred. Parity target is the Scala caller.

pub mod caller;
pub mod clusterer;
pub mod evidence;
pub mod segmenter;
pub mod types;
pub mod walker;

pub use caller::call_structural_variants;
pub use evidence::{
    BreakpointCluster, DepthSegment, DiscordantPair, DiscordantReason, SplitRead,
    SvEvidenceCollection,
};
pub use types::{calculate_confidence, SvAnalysisResult, SvCall, SvCallerConfig, SvType};
