//! Reference-genome + liftover-chain asset retrieval and on-disk cache (plan §4f).
//!
//! Resolves a reference *build* (e.g. `chm13v2.0`) to a usable local file: a decompressed,
//! `.fai`-indexed FASTA, fetched + cached on a miss. Also caches UCSC liftover chains for
//! `du-bio` to parse. Indexing is in-Rust (`noodles::fasta::fs::index`) — no samtools/GATK.
//!
//! Layered below `navigator-app`; depends only on `du-bio` + reqwest/noodles/flate2.

pub mod cache;
pub mod download;
pub mod error;
pub mod gateway;
pub mod index;
pub mod regions;
pub mod registry;

pub use error::RefgenomeError;
pub use gateway::{LiftedPos, ReferenceGateway, RefStatus};
pub use regions::{ChromosomeRegions, Cytoband, GenomeRegions, RegionAnnotation};
pub use registry::{canonical_build, Build, BuildOverride, ReferencePolarity, UserConfig};
