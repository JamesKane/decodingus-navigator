//! Coordinate-sort a BAM, on disk.
//!
//! Same shape as [`crate::revert::collate`] and for the same reason: a WGS alignment does not fit
//! in memory, so this fills a budget, sorts it, spills a run, and k-way merges the runs at the end.
//! Peak memory is the budget plus one buffered block per run, independent of input size.
//!
//! The budget comes from the machine ([`navigator_resource::spill_budget`]), not from a constant.
//! It was 512 MB for everyone, which on a 30x WGS spilled **688 runs** that the merge then opened
//! at once — bounded memory by design, and a great deal of fan-in to buy on a machine with 128 GB
//! sitting idle.
//!
//! ## Runs are ordinary BAM files
//!
//! The revert stage spills a bespoke binary encoding because its records are four small fields.
//! An alignment record is not — it has a CIGAR, a tag dictionary, and per-base sequence and
//! quality — so inventing a serialization for it would mean re-deriving BAM's encoding badly.
//! Each run is written as a real BAM instead, using the same noodles encoder as the final output.
//! It costs a header per run (a few kB against runs of hundreds of MB) and buys a format that is
//! already correct, already tested, and inspectable with any tool when something looks wrong.
//!
//! ## Order
//!
//! SAM's coordinate order is by reference sequence, then by alignment start, with **unplaced reads
//! last**. Those unplaced reads are not an afterthought here: recovering reads the old reference
//! could not place is a large part of why realignment exists, so they have to survive the sort and
//! land somewhere a reader will find them, rather than being dropped for having no coordinate.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::path::{Path, PathBuf};

use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::RecordBuf;

use super::bamio;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

/// How often the record loop asks whether it has been cancelled — same cadence as the walkers.
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// Tuning for [`sort_alignment`].
#[derive(Debug, Clone)]
pub struct SortParams {
    /// Approximate bytes of records held before a sorted run is spilled to scratch.
    pub buffer_bytes: usize,
}

impl Default for SortParams {
    /// Sized from the machine, not from a constant — see [`navigator_resource::spill_budget`], which
    /// also documents `NAVIGATOR_SORT_MB`. The constant this replaced was 512 MB, which spilled 688
    /// runs on a 30x WGS regardless of whether the machine had 8 GB or 128 GB to work with.
    fn default() -> Self {
        Self {
            buffer_bytes: navigator_resource::spill_budget("NAVIGATOR_SORT_MB") as usize,
        }
    }
}

/// What the sort did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SortStats {
    /// Records read, and written — the sort never drops one.
    pub records: u64,
    /// Records with no coordinate, which sort to the end.
    pub unplaced: u64,
    /// Sorted runs spilled. One means everything fit in the budget.
    pub runs: usize,
}

/// Coordinate-sort `input` into `output`, using `scratch` for spilled runs.
///
/// The output header is the input's, with `@HD SO:coordinate` set — the claim readers rely on to
/// decide whether they may binary-search or must scan.
pub fn sort_alignment(
    input: &Path,
    output: &Path,
    scratch: &Path,
    params: &SortParams,
    cancel: &CancelToken,
    progress: &mut dyn FnMut(u64),
) -> Result<SortStats, AnalysisError> {
    std::fs::create_dir_all(scratch).map_err(|e| AnalysisError::io(scratch, e))?;
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AnalysisError::io(parent, e))?;
    }

    let mut reader = bamio::open(input)?;
    let header = reader.read_header().map_err(|e| AnalysisError::io(input, e))?;

    let mut stats = SortStats::default();
    let mut spiller = Spiller::new(scratch, &header, params.buffer_bytes);

    for (i, result) in reader.record_bufs(&header).enumerate() {
        if i as u64 % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
            progress(stats.records);
        }
        let record = result.map_err(|e| AnalysisError::io(input, e))?;
        if sort_key(&record).0 == UNPLACED {
            stats.unplaced += 1;
        }
        stats.records += 1;
        spiller.push(record)?;
    }

    let runs = spiller.finish()?;
    stats.runs = runs.paths.len();

    let sorted_header = with_coordinate_sort_order(header);
    merge(runs, &sorted_header, output, cancel, &mut stats, progress)?;
    progress(stats.records);
    Ok(stats)
}

// ---- spilling -------------------------------------------------------------

/// Accumulates records, spilling a sorted BAM run whenever the budget is reached.
struct Spiller {
    dir: PathBuf,
    header: sam::Header,
    budget: usize,
    buffered: Vec<RecordBuf>,
    buffered_bytes: usize,
    paths: Vec<PathBuf>,
}

impl Spiller {
    fn new(dir: &Path, header: &sam::Header, budget: usize) -> Self {
        Self {
            dir: dir.to_path_buf(),
            header: header.clone(),
            budget: budget.max(1),
            buffered: Vec::new(),
            buffered_bytes: 0,
            paths: Vec::new(),
        }
    }

    fn push(&mut self, record: RecordBuf) -> Result<(), AnalysisError> {
        self.buffered_bytes += heap_bytes(&record);
        self.buffered.push(record);
        if self.buffered_bytes >= self.budget {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> Result<(), AnalysisError> {
        if self.buffered.is_empty() {
            return Ok(());
        }
        // `sort_by_key` (stable) rather than `sort_unstable_by_key`: records that share a
        // coordinate keep their input order, so the whole sort is deterministic. A duplicate
        // marker downstream picks a representative from a group of identical positions, and it
        // should pick the same one on every run.
        self.buffered.sort_by_key(sort_key);

        let path = self.dir.join(format!("sort-run-{:05}.bam", self.paths.len()));
        let mut writer = bamio::create(&path)?;
        writer
            .write_header(&self.header)
            .map_err(|e| AnalysisError::io(&path, e))?;
        for record in &self.buffered {
            writer
                .write_alignment_record(&self.header, record)
                .map_err(|e| AnalysisError::io(&path, e))?;
        }
        bamio::finish(writer, &path)?;

        self.buffered.clear();
        self.buffered_bytes = 0;
        self.paths.push(path);
        Ok(())
    }

    fn finish(mut self) -> Result<Runs, AnalysisError> {
        self.spill()?;
        Ok(Runs { paths: self.paths })
    }
}

/// The spilled runs, deleted on drop so a failed or cancelled sort leaves no scratch behind — at
/// WGS scale these are the size of the alignment itself.
struct Runs {
    paths: Vec<PathBuf>,
}

impl Drop for Runs {
    fn drop(&mut self) {
        for path in &self.paths {
            let _ = std::fs::remove_file(path);
        }
    }
}

// ---- merging --------------------------------------------------------------

/// One record plus the run it came from, ordered by coordinate.
struct Entry {
    key: (u32, u64),
    record: RecordBuf,
    run: usize,
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}
impl Eq for Entry {}
impl Ord for Entry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Ties break on run index so the merge is deterministic; combined with the stable sort
        // within each run, equal-coordinate records come out in a repeatable order.
        self.key.cmp(&other.key).then(self.run.cmp(&other.run))
    }
}
impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

fn merge(
    runs: Runs,
    header: &sam::Header,
    output: &Path,
    cancel: &CancelToken,
    stats: &mut SortStats,
    progress: &mut dyn FnMut(u64),
) -> Result<(), AnalysisError> {
    let mut writer = bamio::create(output)?;
    writer.write_header(header).map_err(|e| AnalysisError::io(output, e))?;

    // Each run is read through its own reader; only one record per run is resident at a time.
    let mut readers = Vec::with_capacity(runs.paths.len());
    for path in &runs.paths {
        let mut reader = bamio::open_many(path)?;
        let run_header = reader.read_header().map_err(|e| AnalysisError::io(path, e))?;
        readers.push(RunReader {
            path: path.clone(),
            header: run_header,
            reader,
        });
    }

    let mut heap = BinaryHeap::with_capacity(readers.len());
    for (run, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next_record()? {
            heap.push(Reverse(Entry {
                key: sort_key(&record),
                record,
                run,
            }));
        }
    }

    let mut written = 0u64;
    while let Some(Reverse(entry)) = heap.pop() {
        if written % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
            progress(written);
        }
        writer
            .write_alignment_record(header, &entry.record)
            .map_err(|e| AnalysisError::io(output, e))?;
        written += 1;

        if let Some(record) = readers[entry.run].next_record()? {
            heap.push(Reverse(Entry {
                key: sort_key(&record),
                record,
                run: entry.run,
            }));
        }
    }

    bamio::finish(writer, output)?;

    // A sort that loses records is the failure this whole stage must not have, and it would be
    // invisible downstream — coverage would simply read low.
    if written != stats.records {
        return Err(AnalysisError::Message(format!(
            "sort lost records: read {}, wrote {written}",
            stats.records
        )));
    }
    Ok(())
}

struct RunReader {
    path: PathBuf,
    header: sam::Header,
    /// Single-threaded by design — see [`bamio::open_many`]. Every run in the merge is open at
    /// once, so a worker pool per run is thousands of threads for work that is already parallel.
    reader: bamio::PlainBamReader,
}

impl RunReader {
    fn next_record(&mut self) -> Result<Option<RecordBuf>, AnalysisError> {
        let mut record = RecordBuf::default();
        let n = self
            .reader
            .read_record_buf(&self.header, &mut record)
            .map_err(|e| AnalysisError::io(&self.path, e))?;
        Ok((n != 0).then_some(record))
    }
}

// ---- keys and helpers -----------------------------------------------------

/// Reference id for a record with no place on the reference. `u32::MAX` sorts after every real
/// reference sequence, which is exactly where SAM wants unplaced reads.
const UNPLACED: u32 = u32::MAX;

fn sort_key(record: &RecordBuf) -> (u32, u64) {
    match record.reference_sequence_id() {
        Some(id) => (id as u32, record.alignment_start().map(|p| p.get() as u64).unwrap_or(0)),
        None => (UNPLACED, 0),
    }
}

/// What one record costs the buffer.
///
/// This used to be the variable-length parts plus a flat 256, which read as a reasonable stand-in
/// for "fixed fields and allocator overhead" and was not one. It left out the tag dictionary
/// entirely, and a mapped record carries a dozen tags — `NM`, `MD`, `AS`, `ms`, `nn`, `tp`, `cm`,
/// `s1`, `s2`, `de`, `rl` from minimap2 alone — each an entry in a `Vec<(Tag, Value)>`. The buffer
/// therefore held well over its stated budget, which mattered little against a constant picked with
/// an unwritten margin and matters a great deal now that the budget is a fraction of the machine
/// (see [`navigator_resource::spill_budget`]).
///
/// Still an estimate: it does not chase a `Value`'s own heap (a string tag's bytes) or a `Vec`'s
/// spare capacity. It is close enough to size a buffer by, and it no longer omits a whole field.
pub(super) fn heap_bytes(record: &RecordBuf) -> usize {
    use noodles::sam::alignment::record::Cigar as _;
    // The record itself sits inline in the buffer's `Vec`, so its size is part of what a record
    // costs — not something to approximate around.
    std::mem::size_of::<RecordBuf>()
        + record.name().map(|n| n.len()).unwrap_or(0)
        + record.sequence().len()
        + record.quality_scores().len()
        + record.cigar().len() * 4
        + record.data().len() * TAG_ENTRY_BYTES
        // Name, sequence, qualities, CIGAR, tags: five vectors, five allocations.
        + 5 * navigator_resource::ALLOCATION_OVERHEAD
}

/// Bytes one tag occupies in a record's `Vec<(Tag, Value)>`.
///
/// Pinned by `a_tag_entry_is_not_larger_than_the_estimate_assumes`, so a noodles upgrade that grows
/// `Value` fails a test here rather than quietly halving the buffer's honesty.
pub(super) const TAG_ENTRY_BYTES: usize = 48;

/// Stamp `@HD SO:coordinate` on the header.
///
/// Not cosmetic: an index is only valid for a coordinate-sorted file, and readers decide whether
/// they may query a region by looking at this. A correctly sorted file that fails to say so is
/// treated as unsorted and re-scanned.
fn with_coordinate_sort_order(mut header: sam::Header) -> sam::Header {
    use noodles::sam::header::record::value::map::header::tag;
    use noodles::sam::header::record::value::{map, Map};

    match header.header_mut() {
        Some(hd) => {
            hd.other_fields_mut()
                .insert(tag::SORT_ORDER, b"coordinate".as_slice().into());
        }
        none => {
            let mut hd = Map::<map::Header>::default();
            hd.other_fields_mut()
                .insert(tag::SORT_ORDER, b"coordinate".as_slice().into());
            *none = Some(hd);
        }
    }
    header
}
