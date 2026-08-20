//! Turn an aligned BAM or CRAM back into the unaligned reads that made it. This is stage A of the
//! realignment pipeline. See `documents/design/realignment-module.md`.
//!
//! It is the job of GATK `RevertSam` plus `SamToFastq`, or of `samtools collate | fastq`, in Rust
//! and on noodles. It does not know about any backend, and that is deliberate. Nothing here knows
//! which mapper the reads go to, so it holds its value even if somebody changes the aligner under
//! it.
//!
//! ## Why this is the hard part
//!
//! To recover the original reads from an alignment in coordinate order is not a filter. It is a
//! regroup. Four things make it awkward:
//!
//! 1. **A read and its mate are far apart.** In coordinate order they can sit gigabases from each
//!    other, so the code can rebuild the pair only from a group by name. At WGS scale, which is
//!    about 10⁹ records, a hash map from a read name to a record is not possible. So [`collate`]
//!    does an external merge sort on disk. It fills a memory budget, sorts, spills a run, and
//!    merges the runs back k at a time. The memory then stays flat at any input size.
//! 2. **An aligner rewrites the read.** A reverse-strand alignment stores SEQ and QUAL as the
//!    reverse complement of what the sequencer gave, so the code must restore both.
//! 3. **Only a primary record carries the full read.** The code drops a secondary record and a
//!    supplementary one. A supplementary record usually carries a hard clip, which means that the
//!    aligner already threw sequence away.
//! 4. **An unmapped read must survive.** It is not an edge case to accept. A read that did not map
//!    on GRCh38 is exactly the read that may land in sequence that CHM13 resolves. Those reads are
//!    the gain of the whole realignment, and they must reach the FASTQ.
//!
//! ## What comes out
//!
//! Paired FASTQ, as `_1.fastq` and `_2.fastq`, and the two stay in step. A file of singletons sits
//! beside them, for anything with no pair. That covers a library with no pairs, a mate whose
//! partner the code dropped, and a read whose flags disagree with themselves.
//!
//! The read names go out bare, with no `/1` and `/2` at the end, as `samtools fastq` writes them.
//! The position in the file says which two reads make a pair, and those suffixes confuse some
//! later tools more than they help.
//!
//! The design records one other option: a uBAM output, which would keep the `@RG` at each read,
//! and not in the header alone. FASTQ is the default, because every aligner takes it. [`writer`]
//! holds the writer on its own, so a uBAM writer can go beside it, and the collation does not
//! change.

mod collate;
mod transform;
mod writer;

use std::path::{Path, PathBuf};

use noodles::sam::alignment::RecordBuf;

use crate::cancel::CancelToken;
use crate::error::AnalysisError;
use crate::reader;

// The reverted-read vocabulary is public because stage B (the mapper) consumes it. `Skipped` is
// not: it is the private reason-code that pairs with `revert_record`, which stays internal.
pub use transform::{Mate, RevertedRead};

/// How often the record loop asks whether somebody cancelled it. It asks often enough that a click
/// feels immediate, and rarely enough that the atomic load never shows in a profile. That is the
/// same reasoning as in the other walkers. See [`crate::cancel`].
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// What to do with a **primary** record whose CIGAR contains a hard clip.
///
/// A hard clip means that the aligner threw sequence away from the record, so nothing can recover
/// the whole read.
///
/// A mainstream aligner puts a hard clip on a supplementary record alone, and the code drops those
/// in any case. But some pipelines give a primary record with a hard clip. To emit such a record
/// as if it were whole would give the mapper a short read, and nobody would see it happen.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HardClipPolicy {
    /// Drop the read and count it. The default: a missing read is visible in the stats, whereas a
    /// truncated one is not visible anywhere.
    #[default]
    Skip,
    /// Emit whatever sequence survived. For inputs where hard-clipped primaries are the norm and a
    /// partial read beats no read.
    Emit,
}

/// The controls of [`revert_alignment`].
#[derive(Debug, Clone)]
pub struct RevertParams {
    /// How many bytes of reverted reads stay in memory before the code spills a sorted run to the
    /// scratch space.
    pub sort_buffer_bytes: usize,
    /// Treatment of hard-clipped primary records.
    pub hard_clipped: HardClipPolicy,
    /// Prefer the `OQ` tag (original, pre-BQSR qualities) over the record's stored qualities.
    pub prefer_original_qualities: bool,
}

impl Default for RevertParams {
    /// The size of the collator comes from the machine, by the same rule as the coordinate sort.
    /// See [`navigator_resource::spill_budget`], which also documents `NAVIGATOR_REVERT_SORT_MB`.
    ///
    /// The constant before it was 256 MB, and its comment gave the reason: it kept the run count
    /// low for a WGS. A reverted read is about 340 bytes, so 256 MB is one run in every million
    /// reads. A 30x WGS then spilled some hundreds of runs, and the merge opened every one.
    fn default() -> Self {
        Self {
            sort_buffer_bytes: navigator_resource::spill_budget("NAVIGATOR_REVERT_SORT_MB") as usize,
            hard_clipped: HardClipPolicy::default(),
            prefer_original_qualities: true,
        }
    }
}

/// What the revert did. It goes into the job log, and into the honest report that the design asks
/// for.
///
/// Every count here exists because the code *dropped or changed* something. A revert that loses
/// reads where nobody sees it is the failure that this whole struct makes impossible.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RevertStats {
    /// The count of records in the input, before any filter.
    pub records_read: u64,
    /// Secondary alignments dropped (`0x100`).
    pub secondary_dropped: u64,
    /// Supplementary alignments dropped (`0x800`).
    pub supplementary_dropped: u64,
    /// The count of primary records that the code dropped or cut short because of a hard clip.
    /// [`HardClipPolicy`] decides which.
    pub hard_clipped: u64,
    /// The count of records that the code dropped because they hold no sequence at all, where
    /// `SEQ` is `*`. There is nothing to revert.
    pub no_sequence_dropped: u64,
    /// The count of records that carried no qualities, where `QUAL` is `*`, and for which the code
    /// **made** a flat phred 40. It keeps those reads, because a mapper can still map a read with
    /// no qualities. But the code invents the qualities that come out, and this count is what
    /// keeps that visible.
    pub qualities_synthesized: u64,
    /// The count of records whose qualities came from the `OQ` tag, and not from `QUAL`.
    pub original_qualities_used: u64,
    /// The count of reads that had no mapping in the input, at flag `0x4`. Those are the gain of
    /// the realignment. See the module documentation.
    pub unmapped_reads: u64,
    /// Reads written across all three output files.
    pub reads_emitted: u64,
    /// Complete pairs written to the `_1`/`_2` files.
    pub pairs: u64,
    /// Reads written to the singletons file.
    pub singletons: u64,
    /// Sorted runs spilled to scratch. One means everything fit in memory.
    pub runs_spilled: usize,
}

/// Where the reverted reads landed, plus what happened on the way.
#[derive(Debug, Clone)]
pub struct RevertOutput {
    pub read1: PathBuf,
    pub read2: PathBuf,
    pub singletons: PathBuf,
    pub stats: RevertStats,
}

/// Revert `path` (BAM or CRAM) into paired FASTQ under `out_dir`.
///
/// A CRAM needs `reference`, and a BAM ignores it, as in [`reader::open_seq`]. The three FASTQ
/// files go into `out_dir`. The spill files of the sort go there too, and the code removes those
/// before it returns.
pub fn revert_alignment(
    path: &Path,
    reference: Option<&Path>,
    out_dir: &Path,
    params: &RevertParams,
    cancel: &CancelToken,
) -> Result<RevertOutput, AnalysisError> {
    let (header, mut sr) = reader::open_seq(path, reference)?;
    let records = sr.records(&header);
    revert_records(records, out_dir, params, cancel)
}

/// The container-independent core of [`revert_alignment`], over any source of records.
///
/// It is separate so that a test can run the pipeline on records that somebody built by hand, and
/// write no BAM fixture. The layer that reads a file format is the job of [`reader`], and the
/// tests there already cover it.
pub fn revert_records(
    records: impl Iterator<Item = Result<RecordBuf, AnalysisError>>,
    out_dir: &Path,
    params: &RevertParams,
    cancel: &CancelToken,
) -> Result<RevertOutput, AnalysisError> {
    std::fs::create_dir_all(out_dir).map_err(|e| AnalysisError::io(out_dir, e))?;

    let mut stats = RevertStats::default();
    let mut collator = collate::Collator::new(out_dir, params.sort_buffer_bytes);

    for (i, record) in records.enumerate() {
        if i as u64 % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
        }
        let record = record?;
        stats.records_read += 1;

        match transform::revert_record(&record, params, &mut stats) {
            Ok(read) => collator.push(read)?,
            Err(skipped) => transform::count_skip(skipped, &mut stats),
        }
    }

    let merged = collator.finish()?;
    stats.runs_spilled = merged.run_count();
    let out = writer::write_fastq(merged, out_dir, &mut stats, cancel)?;

    Ok(RevertOutput {
        read1: out.0,
        read2: out.1,
        singletons: out.2,
        stats,
    })
}

#[cfg(test)]
mod tests;
