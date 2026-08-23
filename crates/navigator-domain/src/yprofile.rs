//! Y-variant profile: the **Y-DNA adapter** over the generic [`crate::consensus`] engine.
//!
//! The reconciliation code is not specific to one DNA type. It holds a vote weighted by quality, a
//! status taxonomy, and a summary, and it lives in [`crate::consensus`]. This module is the Y-DNA
//! view of it.
//!
//! The app collects the calls at each SNP from every source that carries Y. Those sources are the
//! haplogroup placement of a WGS alignment, the chip or BISDNA placement, and the private-Y
//! bucket. It then groups them **by SNP name** through [`reconcile_y`]. The name is independent of
//! the build, because M269 is M269 whether the source aligned to GRCh37 or to GRCh38. It
//! classifies each SNP as confirmed, novel, conflict, or single-source.
//!
//! The mtDNA consumer (variants vs rCRS) and the autosomal consumer reuse the same engine through
//! their own thin adapters. The Y-flavored aliases below keep call sites Y-specific, and there is
//! still one implementation.

pub use crate::consensus::{
    interpret, obs_weight, reconcile as reconcile_y, summarize, to_observed, CallableState as YCallableState,
    ConsensusObs as YObsInput, ConsensusState as YState, ConsensusStatus as YVariantStatus,
    ConsensusSummary as YProfileSummary, ConsensusVariant as YProfileVariant, ObservedProfile, ObservedSource,
    ObservedVariant, SourceObs as YSourceObs, SourceSummary,
};
