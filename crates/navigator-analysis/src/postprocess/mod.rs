//! Stage C of the realignment pipeline. It turns the output of the mapper into an alignment that
//! the rest of Navigator can use. See `documents/design/realignment-module.md`.
//!
//! The mapper gives its records in the order that the reads arrived, and no consumer can use that
//! order. A coverage walk, a region query and the variant caller all need coordinate order, and
//! the compression of a CRAM needs it too. So stage C does this:
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
//! The mark on a duplicate goes onto a short read alone, and that is deliberate. The step finds
//! two reads that start and end at the same place, with the same orientation. Two such fragments
//! are almost surely PCR copies of one original.
//!
//! That reasoning does not hold for HiFi and ONT. Two long reads share their endpoints far less
//! often, and such a library usually has no PCR step. A mark on them would then throw away real
//! coverage. This agrees with standard practice, and with the Stage C note of the design.

pub mod bamio;
mod cram;
mod finalize;
pub mod markdup;
pub mod sort;

pub use bamio::is_complete_bam;
pub use cram::{write_cram, CramOutput};
pub use finalize::{bai_path, finalize_bam, index_bam, FinalizedAlignment};
pub use markdup::{mark_duplicates, MarkDupParams, MarkDupStats};
pub use sort::{sort_alignment, SortParams, SortStats};

#[cfg(test)]
mod tests;
