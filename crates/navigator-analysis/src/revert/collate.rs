//! Collate reverted reads by name, on disk.
//!
//! A WGS BAM in coordinate order holds about 10⁹ records, and its mates lie far apart. To put two
//! mates together, the code needs the records in name order.
//!
//! Two things rule out the two direct methods. A hash map from a name to a record does not fit in
//! memory at that scale. And a second pass over the index costs a full extra decode of the file.
//!
//! So this is a standard **external merge sort**. Fill a fixed memory budget, sort it, spill a run
//! to the scratch space, and repeat. Then merge the runs with a heap, k at a time.
//!
//! The property that matters is the peak memory: the budget, plus one buffered block for each run,
//! and *nothing that grows with the input*. That is what lets one code path revert a 5 GB exome
//! and a 200 GB WGS on the same laptop. It is also why this code writes to disk, instead of a
//! method that looks smarter. The budget itself comes from the machine. See
//! [`navigator_resource::spill_budget`].
//!
//! A run holds a plain binary encoding, where a length comes before each field. It uses no
//! serialization framework. There are three reasons. This one file writes the format and reads it.
//! It is a hot path, over billions of records. And an encoder that you can read answers the
//! question "how many bytes did that read cost", where one that a macro derived does not.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use super::transform::{Mate, RevertedRead};
use crate::error::AnalysisError;

/// The buffer size for a spill of a run, and for the read back. It is large enough that the reads
/// of the merge, on each run, stay sequential. It is also small enough that the runs of a WGS do
/// not add up to real memory. The merge holds all of them open at one time.
const RUN_IO_BUFFER: usize = 256 * 1024;

/// It collects the reverted reads. When they reach the budget, it sorts them and spills a run to
/// the scratch space.
pub struct Collator {
    dir: PathBuf,
    buffer: Vec<RevertedRead>,
    buffered_bytes: usize,
    budget: usize,
    runs: Vec<PathBuf>,
}

impl Collator {
    pub fn new(dir: &Path, budget: usize) -> Self {
        Self {
            dir: dir.to_path_buf(),
            buffer: Vec::new(),
            buffered_bytes: 0,
            budget: budget.max(1),
            runs: Vec::new(),
        }
    }

    /// Add a read, spilling a sorted run first if this one would exceed the budget.
    pub fn push(&mut self, read: RevertedRead) -> Result<(), AnalysisError> {
        self.buffered_bytes += read.heap_bytes();
        self.buffer.push(read);
        if self.buffered_bytes >= self.budget {
            self.spill()?;
        }
        Ok(())
    }

    /// Sort the in-memory buffer and write it out as one run.
    fn spill(&mut self) -> Result<(), AnalysisError> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer.sort_unstable();

        let path = self.dir.join(format!("revert-run-{:05}.bin", self.runs.len()));
        let file = File::create(&path).map_err(|e| AnalysisError::io(&path, e))?;
        // The code paces this write. The spill runs are the largest thing that this pipeline
        // writes, and the peak of the scratch space is here, not in a later stage. So these are
        // exactly the writes that must not collect in the page cache with nothing behind them.
        // See `navigator_resource::PacedFile`.
        let mut w = BufWriter::with_capacity(RUN_IO_BUFFER, navigator_resource::PacedFile::new(file));
        for read in &self.buffer {
            write_read(&mut w, read).map_err(|e| AnalysisError::io(&path, e))?;
        }
        w.flush().map_err(|e| AnalysisError::io(&path, e))?;

        self.buffer.clear();
        self.buffered_bytes = 0;
        self.runs.push(path);
        Ok(())
    }

    /// Spill the remaining reads, and then open the merge over all of the runs.
    pub fn finish(mut self) -> Result<Merged, AnalysisError> {
        self.spill()?;
        Merged::open(self.runs)
    }
}

/// A merge over the spilled runs, k at a time. It gives the reads back in `(name, mate)` order.
///
/// It owns its run files, and it deletes them when it drops. So a revert that somebody cancelled,
/// or that failed, leaves no scratch space behind. That space is tens of gigabytes.
pub struct Merged {
    runs: Vec<PathBuf>,
    readers: Vec<BufReader<File>>,
    heap: BinaryHeap<Reverse<Entry>>,
}

/// A read plus the run it came from. `Ord` is by the read alone so the heap orders by
/// `(name, mate)`; the run index only says where to pull the replacement from.
#[derive(Debug, PartialEq, Eq)]
struct Entry {
    read: RevertedRead,
    run: usize,
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.read.cmp(&other.read)
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Merged {
    fn open(runs: Vec<PathBuf>) -> Result<Self, AnalysisError> {
        let mut readers = Vec::with_capacity(runs.len());
        for path in &runs {
            let file = File::open(path).map_err(|e| AnalysisError::io(path, e))?;
            readers.push(BufReader::with_capacity(RUN_IO_BUFFER, file));
        }

        // Prime the heap with the first read of every run.
        let mut heap = BinaryHeap::with_capacity(readers.len());
        for (run, reader) in readers.iter_mut().enumerate() {
            if let Some(read) = read_read(reader).map_err(|e| AnalysisError::io(&runs[run], e))? {
                heap.push(Reverse(Entry { read, run }));
            }
        }

        Ok(Self { runs, readers, heap })
    }

    /// The count of runs that the code spilled. A value of one means that the input fit in the
    /// memory budget.
    pub fn run_count(&self) -> usize {
        self.runs.len()
    }

    /// The next read in `(name, mate)` order.
    pub fn next_read(&mut self) -> Result<Option<RevertedRead>, AnalysisError> {
        let Some(Reverse(Entry { read, run })) = self.heap.pop() else {
            return Ok(None);
        };
        if let Some(next) = read_read(&mut self.readers[run]).map_err(|e| AnalysisError::io(&self.runs[run], e))? {
            self.heap.push(Reverse(Entry { read: next, run }));
        }
        Ok(Some(read))
    }

    /// The next group of reads that share a name. That group is one template.
    ///
    /// The reads arrive in name order, so a group is the run of equal names at the front, and
    /// nothing more. This is where two mates can come together at last. The records of a template
    /// finally sit beside each other, and coordinate order had taken that away.
    pub fn next_group(&mut self, group: &mut Vec<RevertedRead>) -> Result<bool, AnalysisError> {
        group.clear();
        let Some(first) = self.next_read()? else {
            return Ok(false);
        };
        group.push(first);

        // Look at the top of the heap, which is the next read of all. Compare the names, and do
        // not take it off.
        while let Some(Reverse(entry)) = self.heap.peek() {
            if entry.read.name != group[0].name {
                break;
            }
            let read = self.next_read()?.expect("peek said there is a read");
            group.push(read);
        }
        Ok(true)
    }
}

impl Drop for Merged {
    fn drop(&mut self) {
        // Close the handles before the code removes the files. This is best effort. A cleanup
        // that fails must not hide an error that is already on its way out.
        self.readers.clear();
        for path in &self.runs {
            let _ = std::fs::remove_file(path);
        }
    }
}

// ---- run encoding ---------------------------------------------------------

fn mate_byte(mate: Mate) -> u8 {
    match mate {
        Mate::Unpaired => 0,
        Mate::One => 1,
        Mate::Two => 2,
    }
}

fn mate_from_byte(b: u8) -> Mate {
    match b {
        1 => Mate::One,
        2 => Mate::Two,
        _ => Mate::Unpaired,
    }
}

/// `name_len:u32 | name | mate:u8 | seq_len:u32 | sequence | qualities`
///
/// The qualities carry no length before them. [`super::transform::revert_record`] makes sure that
/// their length matches the sequence. To write a number that the code already knows would cost 4
/// bytes at each read, over billions of reads, and it would buy nothing.
fn write_read<W: Write>(w: &mut W, read: &RevertedRead) -> std::io::Result<()> {
    w.write_all(&(read.name.len() as u32).to_le_bytes())?;
    w.write_all(&read.name)?;
    w.write_all(&[mate_byte(read.mate)])?;
    w.write_all(&(read.sequence.len() as u32).to_le_bytes())?;
    w.write_all(&read.sequence)?;
    w.write_all(&read.qualities)
}

fn read_read<R: Read>(r: &mut R) -> std::io::Result<Option<RevertedRead>> {
    let Some(name_len) = read_u32_or_eof(r)? else {
        return Ok(None);
    };
    let mut name = vec![0u8; name_len as usize];
    r.read_exact(&mut name)?;

    let mut mate = [0u8; 1];
    r.read_exact(&mut mate)?;

    let mut len = [0u8; 4];
    r.read_exact(&mut len)?;
    let seq_len = u32::from_le_bytes(len) as usize;

    let mut sequence = vec![0u8; seq_len];
    r.read_exact(&mut sequence)?;
    let mut qualities = vec![0u8; seq_len];
    r.read_exact(&mut qualities)?;

    Ok(Some(RevertedRead {
        name,
        mate: mate_from_byte(mate[0]),
        sequence,
        qualities,
    }))
}

/// A clean end of a run is the one place where a short read is correct. Anywhere else, the data
/// carries damage, and `read_exact` says so.
fn read_u32_or_eof<R: Read>(r: &mut R) -> std::io::Result<Option<u32>> {
    let mut buf = [0u8; 4];
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..])? {
            0 if filled == 0 => return Ok(None),
            0 => return Err(std::io::ErrorKind::UnexpectedEof.into()),
            n => filled += n,
        }
    }
    Ok(Some(u32::from_le_bytes(buf)))
}
