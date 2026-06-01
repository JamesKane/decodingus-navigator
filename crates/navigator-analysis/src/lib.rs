//! Navigator analysis — the htsjdk/GATK replacement, Navigator-side.
//!
//! Owns the `noodles` BAM/CRAM/FASTA/BGZF/index I/O layer (kept out of shared
//! `du-bio`, which stays IO-light coordinate math + text parsing), the ported GATK
//! walkers (`coverage`, `read_metrics`, `sv`, `sex`), and the purpose-built haploid
//! variant caller: force-call genotyping at known sites plus de-novo Y/mtDNA discovery
//! for private-variant matching and branch creation.
//!
//! Built on `du-bio` for liftover/callable/coordinate primitives. A GATK-vs-Rust
//! golden-truth parity harness gates cutover. Implemented in roadmap phases 2–3.

pub mod caller;
pub mod contig;
pub mod coverage;
pub mod error;
pub mod parity;
pub mod read_metrics;
pub mod realign;
pub mod sex;
pub mod sv;

pub use error::AnalysisError;
