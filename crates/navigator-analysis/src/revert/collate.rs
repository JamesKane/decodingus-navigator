//! Collate reverted reads by name, on disk.
//!
//! A coordinate-sorted WGS BAM holds ~10⁹ records and its mates are scattered, so pairing needs the
//! records in name order. Two things rule out the obvious approaches: a name→record hash map does
//! not fit in memory at that scale, and a second indexed pass costs a full extra decode of the
//! file. So this is a textbook **external merge sort** — fill a fixed memory budget, sort it, spill
//! a run to scratch, repeat; then merge the runs with a k-way heap.
//!
//! The property that matters is that peak memory is the budget plus one buffered block per run,
//! *independent of input size*. That is what lets the same code path revert a 5 GB exome and a
//! 200 GB WGS on the same laptop, and it is why this is disk-backed rather than clever.
//!
//! Runs use a plain length-prefixed binary encoding rather than a serialization framework: the
//! format is written and read in this one file, it is a hot path measured in billions of records,
//! and an explicit encoder is easier to reason about than a derived one when the question is "how
//! many bytes did that read cost".

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use super::transform::{Mate, RevertedRead};
use crate::error::AnalysisError;

/// Buffer size for run spill/read-back. Large enough that the merge's per-run reads stay
/// sequential, small enough that a few dozen concurrent runs don't add up to real memory.
const RUN_IO_BUFFER: usize = 256 * 1024;

/// Accumulates reverted reads, spilling sorted runs to scratch when the budget is reached.
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
        // Paced: the spill runs are the biggest thing this pipeline writes — the scratch peak is
        // here, not in the post-processing stages — so they are exactly what must not be allowed to
        // pile up dirty in the page cache. See `crate::resource::PacedFile`.
        let mut w = BufWriter::with_capacity(RUN_IO_BUFFER, crate::resource::PacedFile::new(file));
        for read in &self.buffer {
            write_read(&mut w, read).map_err(|e| AnalysisError::io(&path, e))?;
        }
        w.flush().map_err(|e| AnalysisError::io(&path, e))?;

        self.buffer.clear();
        self.buffered_bytes = 0;
        self.runs.push(path);
        Ok(())
    }

    /// Spill whatever is left and open the merge over all runs.
    pub fn finish(mut self) -> Result<Merged, AnalysisError> {
        self.spill()?;
        Merged::open(self.runs)
    }
}

/// A k-way merge over the spilled runs, yielding reads in `(name, mate)` order.
///
/// Owns its run files and deletes them on drop, so a cancelled or failed revert does not leave
/// tens of gigabytes of scratch behind.
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

    /// How many runs were spilled. One means the input fit in the memory budget.
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

    /// The next group of reads sharing a name — i.e. one template.
    ///
    /// Reads arrive name-ordered, so a group is just the run of equal names at the front. This is
    /// where pairing becomes possible: a template's records are finally adjacent, which is the
    /// thing coordinate order took away.
    pub fn next_group(&mut self, group: &mut Vec<RevertedRead>) -> Result<bool, AnalysisError> {
        group.clear();
        let Some(first) = self.next_read()? else {
            return Ok(false);
        };
        group.push(first);

        // Peek: the heap's top is the next read overall, so compare names without consuming.
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
        // Close the handles before unlinking; best-effort, since a failed cleanup must not mask
        // whatever error is already unwinding.
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
/// Qualities are not length-prefixed: [`super::transform::revert_record`] guarantees they match the
/// sequence length, and re-encoding a number we already know would cost 4 bytes per read across
/// billions of reads for nothing.
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

/// A clean end-of-run is the only place a short read is expected; anywhere else it is corruption
/// and `read_exact` will say so.
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
