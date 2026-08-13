//! Stage C of the realignment pipeline — turning the mapper's output into an alignment the rest
//! of Navigator can use (`documents/design/realignment-module.md`).
//!
//! The mapper emits records in the order the reads arrived, which is useless to every consumer:
//! coverage walks, region queries, and the variant caller all assume coordinate order, and CRAM's
//! compression assumes it too. So stage C is:
//!
//! ```text
//! mapped BAM (read order)
//!     │
//!  sort ──────────► coordinate order          <- [`sort`]
//!     │
//!  mark duplicates (short read only)          <- next
//!     │
//!  CRAM + .crai                               <- next
//! ```
//!
//! Duplicate marking is deliberately short-read-only. It identifies reads that start and end at
//! the same place with the same orientation, on the reasoning that two such fragments are almost
//! certainly PCR copies of one original. For HiFi and ONT that reasoning does not hold — long
//! reads genuinely share endpoints far less often, and the libraries are usually PCR-free — so
//! marking them would discard real coverage. This matches standard practice and the design's
//! Stage C note.

pub mod bamio;
mod cram;
pub mod markdup;
pub mod sort;

pub use cram::{write_cram, CramOutput};
pub use markdup::{mark_duplicates, MarkDupParams, MarkDupStats};
pub use sort::{sort_alignment, SortParams, SortStats};

#[cfg(test)]
mod tests;
