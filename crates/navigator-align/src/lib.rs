//! Read mapping for the realignment module — stage B of
//! `documents/design/realignment-module.md`.
//!
//! Takes the reads [`navigator-analysis`'s revert stage](../navigator_analysis/revert/index.html)
//! recovered and maps them to a new reference. Everything here is about doing that within a
//! desktop's memory, because that — not accuracy, and no longer platform support — is the module's
//! binding constraint.
//!
//! ## The backend
//!
//! The default mapper is [`minimap2-pure-rs`], a pure-Rust translation of minimap2 v2.31. It was
//! chosen over linking the C library through FFI because it needs no C toolchain, which means
//! Windows and every other Rust target build unchanged, and because a parity test measured it
//! 99.74% byte-identical to the C implementation with **zero disagreements at MAPQ > 0**.
//!
//! That crate describes itself as an "LLM-mediated faithful translation" and asks users to stay
//! alert to bugs. The mitigation is the `ffi` feature: the same work can be run through the C
//! minimap2 and the results diffed. It is off by default because enabling it requires a C
//! toolchain — and note the feature spelling in `Cargo.toml`, because the obvious one silently
//! links htslib.
//!
//! ## Memory
//!
//! A reference is indexed in *parts* of at most [`BatchSize`] bases, and exactly one part is
//! resident at a time. Building CHM13's index as a single part costs ~19 GiB; building it in
//! 1 Gbase parts costs 11.7 GiB, produces the same output, and takes the same wall time. See
//! [`batch`] for the measured table and for why bigger parts are preferred within the budget.
//!
//! [`BatchSize::for_this_machine`] reads the machine's RAM and picks for itself. That is the
//! intended entry point: this module's users click a button, and "bases per index part" is not a
//! question they can be asked — a wrong answer is an out-of-memory failure, not a preference.
//!
//! ## What is here so far
//!
//! [`preset`] (which mapper preset a run's reads need), [`batch`] (the memory control), and
//! [`index`] (the `.mmi` cache and the part-by-part build). The mapping pass itself — running
//! reads against each part and merging the per-part results, which is what makes a split index
//! produce the same alignments as a whole one — is the next piece.

pub mod batch;
pub mod error;
pub mod index;
pub mod map;
pub mod preset;

pub use batch::BatchSize;
pub use error::AlignError;
pub use map::{map_reads, MapParams, MapStats};
pub use preset::Preset;
