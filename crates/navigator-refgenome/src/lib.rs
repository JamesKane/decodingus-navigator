//! Reference-genome + liftover-chain asset retrieval and on-disk cache (plan §4f).
//!
//! It resolves a reference *build*, such as `chm13v2.0`, to a local file that the app can use.
//! That file is a decompressed FASTA with a `.fai` index. On a cache miss it fetches the file and
//! caches it. It also caches the UCSC liftover chains, for `du-bio` to parse. The index step is
//! pure Rust, through `noodles::fasta::fs::index`, so it needs no samtools and no GATK.
//!
//! This crate sits below `navigator-app`. It depends on `du-bio` alone, plus reqwest, noodles, and
//! flate2.

pub mod cache;
pub mod download;
pub mod error;
pub mod gateway;
pub mod index;
pub mod regions;
pub mod registry;
pub mod vcf_lift;

pub use error::RefgenomeError;
pub use gateway::{LiftedPos, RefStatus, ReferenceGateway, VerifyOutcome};
pub use regions::{ChromosomeRegions, Cytoband, GenomeRegions, RegionAnnotation};
pub use registry::{canonical_build, Build, BuildOverride, ReferencePolarity, UserConfig};
pub use vcf_lift::{VcfLiftOpts, VcfLiftStats};
