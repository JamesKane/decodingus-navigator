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

pub mod ancestry;
pub mod caller;
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
pub mod library_stats;
pub mod manifest;
pub mod mask;
pub mod mastervar;
/// mtDNA variant derivation + CHM13 `chrM`↔rCRS liftover. Moved to the shared `du-bio` crate
/// so the AppView and Navigator share one implementation; re-exported here under the original
/// path so existing `navigator_analysis::mtvariants::…` call sites are unchanged.
pub use du_bio::mt as mtvariants;
pub mod parity;
pub mod preflight;
pub mod probe;
pub mod read_metrics;
pub mod reader;
pub mod readview;
pub mod realign;
pub mod reassembly;
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

pub use error::{guard_walk, AnalysisError};
