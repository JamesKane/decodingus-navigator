//! Revert an aligned BAM/CRAM back to the unaligned reads it was built from — stage A of the
//! realignment pipeline (`documents/design/realignment-module.md`).
//!
//! This is the GATK `RevertSam` + `SamToFastq` / `samtools collate | fastq` job, in Rust on
//! noodles. It is deliberately **backend-agnostic**: nothing here knows which mapper the reads are
//! headed for, so it is worth having even if the aligner decision changes underneath it.
//!
//! ## Why this is the hard part
//!
//! Recovering the original reads from a coordinate-sorted alignment is not a filter, it is a
//! regrouping. Four things make it awkward:
//!
//! 1. **Mates are far apart.** In coordinate order a read and its mate can sit gigabases away, so
//!    pairing requires grouping by name. A read-name→record hash map is not an option at WGS scale
//!    (~10⁹ records), so [`collate`] does a disk-backed external merge sort instead: fill a memory
//!    budget, sort, spill a run, then k-way merge the runs back. Memory stays flat and bounded
//!    regardless of input size, which is the whole point.
//! 2. **Aligners rewrite the read.** A reverse-strand alignment stores SEQ and QUAL
//!    reverse-complemented relative to the sequencer's output, so both must be restored.
//! 3. **Only primaries carry the full read.** Secondary and supplementary records are dropped;
//!    supplementaries are typically hard-clipped, meaning sequence has already been discarded.
//! 4. **Unmapped reads must survive.** They are not an edge case to tolerate — reads that failed to
//!    map on GRCh38 are exactly the ones that may land in CHM13-resolved sequence, so they are the
//!    realignment payoff and have to reach the FASTQ.
//!
//! ## What comes out
//!
//! Paired FASTQ (`_1.fastq` / `_2.fastq`, kept in lockstep) plus a singletons file for anything
//! that did not pair — an unpaired library, a mate whose partner was dropped, or a read whose
//! flags disagree with themselves. Read names are written bare, without `/1` and `/2` suffixes,
//! matching `samtools fastq`: the pairing is carried by file position, and suffixes confuse some
//! downstream tools more than they help.
//!
//! uBAM output (which would preserve `@RG` per-read rather than only in the header) is the other
//! option the design records; FASTQ is the default because every aligner takes it. The writer is
//! isolated in [`writer`] so a uBAM sibling can be added without disturbing collation.

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

/// How often the record loop asks whether it has been cancelled. Frequent enough that a click
/// feels immediate, rare enough that the atomic load never shows up in a profile — the same
/// reasoning as the other walkers (see [`crate::cancel`]).
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// What to do with a **primary** record whose CIGAR contains a hard clip.
///
/// Hard clipping means the aligner discarded sequence from the record, so the read cannot be fully
/// recovered. Mainstream aligners hard-clip only supplementary records (which we drop anyway), but
/// some pipelines emit hard-clipped primaries, and emitting those as if whole would silently feed
/// a truncated read to the mapper.
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

/// Tuning for [`revert_alignment`].
#[derive(Debug, Clone)]
pub struct RevertParams {
    /// Bytes of reverted reads held in memory before a sorted run is spilled to scratch.
    pub sort_buffer_bytes: usize,
    /// Treatment of hard-clipped primary records.
    pub hard_clipped: HardClipPolicy,
    /// Prefer the `OQ` tag (original, pre-BQSR qualities) over the record's stored qualities.
    pub prefer_original_qualities: bool,
}

impl Default for RevertParams {
    /// The collator is sized from the machine by the same rule as the coordinate sort — see
    /// [`navigator_resource::spill_budget`], which also documents `NAVIGATOR_REVERT_SORT_MB`. The
    /// constant this replaced was 256 MB, described as keeping the run count low for a WGS; at
    /// ~340 bytes a reverted read that is a run every million reads, so a 30x WGS spilled several
    /// hundred of them and the merge opened every one.
    fn default() -> Self {
        Self {
            sort_buffer_bytes: navigator_resource::spill_budget("NAVIGATOR_REVERT_SORT_MB") as usize,
            hard_clipped: HardClipPolicy::default(),
            prefer_original_qualities: true,
        }
    }
}

/// What the revert did, for the job log and for the honest reporting the design asks for.
///
/// Every count here exists because something was *dropped or changed*; a revert that silently
/// loses reads is the failure mode this whole struct is here to make impossible.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RevertStats {
    /// Records seen in the input, before any filtering.
    pub records_read: u64,
    /// Secondary alignments dropped (`0x100`).
    pub secondary_dropped: u64,
    /// Supplementary alignments dropped (`0x800`).
    pub supplementary_dropped: u64,
    /// Primary records dropped or truncated because of hard clipping, per [`HardClipPolicy`].
    pub hard_clipped: u64,
    /// Records dropped for having no sequence at all (`SEQ` = `*`) — nothing to revert.
    pub no_sequence_dropped: u64,
    /// Records whose qualities were absent (`QUAL` = `*`) and were therefore **synthesized** at a
    /// flat phred 40. The reads are kept because a read with no qualities is still mappable, but
    /// the qualities that come out are invented and this count is how that stays visible.
    pub qualities_synthesized: u64,
    /// Records whose qualities came from the `OQ` tag rather than `QUAL`.
    pub original_qualities_used: u64,
    /// Reads that were unmapped in the input (`0x4`) — the realignment payoff; see the module docs.
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
/// `reference` is required for CRAM and ignored for BAM, matching [`reader::open_seq`]. `out_dir`
/// receives the three FASTQ files and is also used for the sort's spill files, which are removed
/// before returning.
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
/// Split out so the pipeline can be tested on hand-built records without writing BAM fixtures —
/// the file-format layer is [`reader`]'s job and is already covered there.
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
