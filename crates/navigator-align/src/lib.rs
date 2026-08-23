//! Read mapping for the realignment module. This is stage B of
//! `documents/design/realignment-module.md`.
//!
//! This module takes the reads that
//! [`navigator-analysis`'s revert stage](../navigator_analysis/revert/index.html) recovered, and
//! maps them to a new reference. It does that work inside a desktop's memory. Memory is the
//! constraint that limits the module. Accuracy is not, and platform support is no longer.
//!
//! ## The backend
//!
//! The default mapper is [`minimap2-pure-rs`], a pure-Rust translation of minimap2 v2.31. The
//! module uses it, and does not link the C library through FFI. It needs no C toolchain, so
//! Windows and every other Rust target build unchanged. A parity test also measured it 99.74%
//! byte-identical to the C implementation, with **zero disagreements at MAPQ > 0**.
//!
//! That crate describes itself as an "LLM-mediated faithful translation", and it asks users to
//! look for bugs. So its output needs a check against the original. Do that check **outside the
//! application**: run upstream minimap2 over the same reads and compare the two outputs. Do not
//! keep a second backend here.
//!
//! A link to the C library would put a C toolchain, an `unsafe` surface, and a Windows-unproven
//! dependency into the artifact that users install. That is a poor trade for a development
//! activity that no user ever runs.
//!
//! ## Memory
//!
//! This module indexes a reference in *parts* of not more than [`BatchSize`] bases. Exactly one
//! part stays in memory at a time. An index of CHM13 that is one part costs about 19 GiB. An index
//! in parts of 1 Gbase costs 11.7 GiB, gives the same output, and takes the same wall time. See
//! [`batch`] for the measured table, and for why a larger part is better inside the budget.
//!
//! [`BatchSize::for_this_machine`] reads the machine's RAM and chooses for itself. That is the
//! intended entry point. Users of this module click a button. Nobody can ask them how many bases
//! an index part holds, because a wrong answer is an out-of-memory failure and not a preference.
//!
//! ## What is here so far
//!
//! - [`preset`] finds which mapper preset a run's reads need.
//! - [`batch`] is the memory control.
//! - [`index`] holds the `.mmi` cache and the part-by-part build.
//! - [`map`] does single-end mapping. It also does the cross-part merge, which makes a split index
//!   give the same alignments as a whole one.
//! - [`pe`] does paired-end mapping, which is what `sr` and most vendor WGS need.
//! - [`output`] writes SAM, BAM and CRAM through noodles.
//!
//! A job starts at [`index::ensure_cached_index`], then calls [`map::map_reads`] or
//! [`pe::map_pairs`]. The first resolves the cache location and sets the index size for the
//! machine. The other two write BAM by default.
//!
//! Stage B does not sort, and it does not mark duplicates. That is deliberate: those two steps are
//! stage C. CRAM belongs after the sort, and not here.

pub mod batch;
pub mod error;
pub mod index;
pub mod map;
pub mod output;
pub mod pe;
pub mod preset;

pub use batch::BatchSize;
pub use error::AlignError;
pub use map::{map_reads, MapParams, MapStats};
pub use output::OutputFormat;
pub use pe::map_pairs;
pub use preset::Preset;
