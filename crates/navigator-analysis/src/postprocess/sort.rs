//! Coordinate-sort a BAM, on disk.
//!
//! This has the same shape as [`crate::revert::collate`], and for the same reason. A WGS alignment
//! does not fit in memory. So this code fills a budget, sorts it, spills a run, and then merges
//! the runs k at a time at the end. The peak memory is the budget, plus one buffered block for
//! each run, and it does not grow with the input.
//!
//! The budget comes from the machine. See [`navigator_resource::spill_budget`]. It is not a
//! constant. It was 512 MB for everybody, and on a 30x WGS that spilled **688 runs**, which the
//! merge then opened at one time. The memory had a bound by design. But that is a great deal of
//! fan-in to buy on a machine with 128 GB that nothing uses.
//!
//! ## A run is an ordinary BAM file
//!
//! The revert stage spills a binary encoding of its own, because its records hold four small
//! fields. An alignment record does not. It holds a CIGAR, a tag dictionary, and a sequence and a
//! quality at each base. To invent a serialization for that would mean a second, worse version of
//! the BAM encoding.
//!
//! So each run goes out as a real BAM, through the same noodles encoder as the final output. That
//! costs one header for each run, which is a few kB against runs of hundreds of MB. It buys a
//! format that is already correct, already tested, and that any tool can open when something looks
//! wrong.
//!
//! ## The order
//!
//! The coordinate order of SAM goes by reference sequence, then by alignment start, and it puts
//! **the reads with no place last**. Those reads matter here. To recover a read that the old
//! reference could not place is a large part of why realignment exists. So they must survive the
//! sort, and land where a reader will find them. The sort must not drop a read because it has no
//! coordinate.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::path::{Path, PathBuf};

use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::RecordBuf;

use super::bamio;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

/// How often the record loop asks whether somebody cancelled it. This is the same rate as in the
/// walkers.
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// The controls of [`sort_alignment`].
#[derive(Debug, Clone)]
pub struct SortParams {
    /// About how many bytes of records the code holds before it sorts them and spills a run to
    /// the scratch space.
    pub buffer_bytes: usize,
}

impl Default for SortParams {
    /// The size comes from the machine, and not from a constant. See
    /// [`navigator_resource::spill_budget`], which also documents `NAVIGATOR_SORT_MB`. The constant
    /// before it was 512 MB. That spilled 688 runs on a 30x WGS, whether the machine had 8 GB or
    /// 128 GB to work with.
    fn default() -> Self {
        Self {
            buffer_bytes: navigator_resource::spill_budget("NAVIGATOR_SORT_MB") as usize,
        }
    }
}

/// What the sort did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SortStats {
    /// The count of records that the code read, and wrote. The sort never drops one.
    pub records: u64,
    /// Records with no coordinate, which sort to the end.
    pub unplaced: u64,
    /// Sorted runs spilled. One means everything fit in the budget.
    pub runs: usize,
}

/// Sort `input` into `output`, in coordinate order. It uses `scratch` for the runs that it spills.
///
/// The output header is the header of the input, with `@HD SO:coordinate` set. A reader depends on
/// that claim to decide whether it may do a binary search, or must scan.
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

/// It collects records. Whenever they reach the budget, it sorts them and spills a BAM run.
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
        // This uses `sort_by_key`, which is stable, and not `sort_unstable_by_key`. Records that
        // share a coordinate then keep their input order, so the whole sort is deterministic. A
        // later step marks the duplicates. It takes one record from a group at the same position,
        // and it must take the same one on every run.
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

/// The runs that the code spilled. They go away when this drops, so a sort that failed, or that
/// somebody cancelled, leaves no scratch space behind. At WGS scale they are as large as the
/// alignment itself.
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
        // A tie breaks on the run index, so the merge is deterministic. With the stable sort
        // inside each run, two records at the same coordinate always come out in the same
        // order.
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

    // Each run has its own reader. Only one record of each run sits in memory at a time.
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

    // A sort that loses records is the one failure that this stage must not have. No later step
    // would see it, because the coverage would only read low.
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
    /// It runs on one thread, and that is deliberate. See [`bamio::open_many`]. Every run of the
    /// merge is open at one time. A worker pool for each run is then thousands of threads, for
    /// work that already runs in parallel.
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
/// This was the parts of variable length, plus a flat 256. That 256 read as a reasonable stand-in
/// for the fixed fields and the allocator overhead, and it was not one. It left the tag dictionary
/// out completely.
///
/// A mapped record carries a dozen tags. minimap2 alone gives `NM`, `MD`, `AS`, `ms`, `nn`, `tp`,
/// `cm`, `s1`, `s2`, `de` and `rl`. Each one is an entry in a `Vec<(Tag, Value)>`. So the buffer
/// held well over the budget that it stated.
///
/// That mattered little against a constant that somebody chose with a margin nobody wrote down. It
/// matters a great deal now that the budget is a fraction of the machine. See
/// [`navigator_resource::spill_budget`].
///
/// This is still an estimate. It does not follow the own heap of a `Value`, which holds the bytes
/// of a string tag. It also does not count the spare capacity of a `Vec`. It is close enough to
/// size a buffer by, and it no longer leaves a whole field out.
pub(super) fn heap_bytes(record: &RecordBuf) -> usize {
    use noodles::sam::alignment::record::Cigar as _;
    // The record itself sits inside the `Vec` of the buffer. Its size is part of what a record
    // costs, and not something to work around.
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
/// `a_tag_entry_is_not_larger_than_the_estimate_assumes` holds this value. So an upgrade of
/// noodles that makes `Value` larger fails a test here. Without that test, the buffer would hold
/// twice what it says, and nobody would see it.
pub(super) const TAG_ENTRY_BYTES: usize = 48;

/// Stamp `@HD SO:coordinate` on the header.
///
/// This is not for appearance. An index is correct only for a file in coordinate order, and a
/// reader reads this to decide whether it may query a region. A file that the code sorted
/// correctly, and that does not say so, reads as unsorted, and every reader scans it again.
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
