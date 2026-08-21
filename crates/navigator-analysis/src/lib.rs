//! The analysis crate of Navigator. It replaces htsjdk and GATK, on the Navigator side.
//!
//! It owns the I/O layer over `noodles`, for BAM, CRAM, FASTA, BGZF and the index files. That layer
//! stays out of the shared `du-bio` crate. `du-bio` holds coordinate arithmetic, and it reads text,
//! and it does little I/O.
//!
//! It also owns the GATK walkers that this project ported: `coverage`, `read_metrics`, `sv` and
//! `sex`. And it owns the haploid variant caller that somebody built for this purpose.
//!
//! That caller does two things. It force-calls a genotype at a known site. And it discovers
//! de-novo on the Y and the mtDNA, for private-variant matching and for a new branch.
//!
//! It stands on `du-bio` for the primitives of liftover, callability and coordinates. A parity
//! harness, of GATK against Rust, over a golden truth, gates the cutover. Phases 2 and 3 of the
//! roadmap built it.

pub mod ancestry;
pub mod archaic;
pub mod archaic_match;
pub mod archaic_segments;
pub mod caller;
pub mod callset;
pub mod cancel;
pub mod contig;
pub mod coverage;
pub mod error;
pub mod genotype;
pub mod gvcf;
pub mod gzio;
pub mod haplo;
pub mod heteroplasmy;
pub mod ibd;
pub mod ibd_attest;
pub mod ibd_panel;
pub mod index;
pub mod lai;
pub mod library_stats;
pub mod manifest;
pub mod mask;
pub mod mastervar;
/// The derivation of the mtDNA variants, and the liftover between the CHM13 `chrM` and the rCRS.
/// This code moved to the shared `du-bio` crate, so that the AppView and Navigator share one
/// implementation. This module exports it again, under its original path, so that a call site that
/// says `navigator_analysis::mtvariants::…` does not change.
pub use du_bio::mt as mtvariants;
pub mod parity;
pub mod phasing;
pub mod postprocess;
pub mod preflight;
pub mod probe;
pub mod read_metrics;
pub mod reader;
pub mod readview;
pub mod realign;
pub mod reassembly;
pub mod revert;
pub mod roh;
pub mod scan;
pub mod sex;
pub mod sidecar;
pub mod strcaller;
pub mod strmarker;
pub mod strref;
pub mod sv;
pub mod testtype;
pub mod unified;
pub mod vcf;

pub use cancel::CancelToken;
pub use error::{guard_walk, AnalysisError};
