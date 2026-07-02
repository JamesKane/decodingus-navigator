//! Y-variant profile — the **Y-DNA adapter** over the generic [`crate::consensus`] engine.
//!
//! The reconciliation machinery (quality-weighted voting, status taxonomy, summary) is DNA-type
//! agnostic and lives in [`crate::consensus`]; this module is the Y-DNA view of it. The app gathers
//! each Y-bearing source's per-SNP calls (a WGS alignment's haplogroup placement, the chip/BISDNA
//! placement, the private-Y bucket), groups them **by SNP name** (build-independent — M269 is M269
//! whether the source aligned to GRCh37 or GRCh38) via [`reconcile_y`], and classifies each SNP as
//! confirmed / novel / conflict / single-source. The mtDNA (variants vs rCRS) and autosomal consumers
//! reuse the same engine through their own thin adapters.
//!
//! The Y-flavored aliases below keep call sites read as Y-specific while sharing one implementation.

pub use crate::consensus::{
    interpret, obs_weight, reconcile as reconcile_y, summarize, to_observed, CallableState as YCallableState,
    ConsensusObs as YObsInput, ConsensusState as YState, ConsensusStatus as YVariantStatus,
    ConsensusSummary as YProfileSummary, ConsensusVariant as YProfileVariant, ObservedProfile, ObservedSource,
    ObservedVariant, SourceObs as YSourceObs, SourceSummary,
};
