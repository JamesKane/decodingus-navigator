//! Multithreaded BAM I/O for the post-processing stages.
//!
//! Stage C moves the whole alignment through BGZF four times over — the sort reads the mapped BAM,
//! writes its spilled runs, reads them back to merge, and writes the sorted output — and duplicate
//! marking and CRAM emission each read it once more. On a 30x WGS that is hundreds of GB of
//! inflate and deflate, and every pass of it was running on one thread: the first WGS-scale run
//! measured the sort at 4 h 44 m, the most expensive stage in the pipeline, ahead of the mapping
//! it feeds.
//!
//! BGZF is a *block*-gzip stream, so compression and decompression parallelize across blocks while
//! the record stream stays sequential and byte-identical. That is the same reasoning
//! [`crate::reader::open_seq`] already applies to reading vendor BAMs; these stages simply never
//! got it, because they were written against the plain constructors.
//!
//! Compression *level* is deliberately unchanged. These are all intermediates consumed once, so a
//! faster level is tempting, but it would inflate the scratch footprint that
//! `navigator-app::realign_job`'s preflight is calibrated against — a separate change, with a
//! recalibration attached.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use noodles::{bam, bgzf};

use crate::error::AnalysisError;

/// A BAM reader whose block decompression runs on a worker pool.
pub(crate) type BamReader = bam::io::Reader<bgzf::io::MultithreadedReader<File>>;

/// A BAM writer whose block compression runs on a worker pool.
pub(crate) type BamWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<PacedFile>>>;

const WRITE_BUFFER: usize = 1 << 20;

/// Bytes written between forced flushes to disk. See [`PacedFile`].
const DEFAULT_SYNC_MB: u64 = 256;

/// How much may be left dirty before this stream pushes it to disk.
///
/// `NAVIGATOR_IO_SYNC_MB=0` turns the pacing off and restores the old behaviour.
fn sync_interval() -> u64 {
    std::env::var("NAVIGATOR_IO_SYNC_MB")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(DEFAULT_SYNC_MB)
        * 1024
        * 1024
}

/// A file that will not let an unbounded amount of its output sit dirty in the page cache.
///
/// Without this, a stage writes as fast as it can into the page cache and leaves write-back to the
/// operating system, which sounds like the right division of labour and is, until the volume gets
/// far enough out of scale. The 2026-08-13 WGS sort dirtied 549 GB of file-backed memory over one
/// run — enough that macOS filed a disk-writes resource notice against the process for exceeding
/// its sustained write-back limit by 1.4x, and enough that WindowServer's main thread missed a
/// 40-second watchdog check-in and was killed, taking the login session and this six-hour job with
/// it.
///
/// Flushing on a byte cadence caps how much can be outstanding at once. The write path pays for
/// its own I/O as it goes instead of handing the machine a debt to settle later, which is also why
/// this is not obviously a slowdown: the same bytes reach the same disk, in steadier instalments
/// rather than in storms.
///
/// `sync_data` rather than `sync_all` — the file's contents must be durable, its metadata need not
/// be, and on a stream this size the difference is many thousands of inode updates. It is
/// std's portable spelling: `fdatasync` where there is one, `FlushFileBuffers` on Windows.
pub(crate) struct PacedFile {
    file: File,
    since_sync: u64,
    interval: u64,
}

impl PacedFile {
    fn new(file: File) -> Self {
        Self {
            file,
            since_sync: 0,
            interval: sync_interval(),
        }
    }
}

impl Write for PacedFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.file.write(buf)?;
        crate::resource::record_bytes_written(written as u64);

        if self.interval > 0 {
            self.since_sync += written as u64;
            if self.since_sync >= self.interval {
                self.file.sync_data()?;
                self.since_sync = 0;
            }
        }

        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

/// The BGZF end-of-file marker: an empty deflate block, written last by every conforming writer.
///
/// Not exported by noodles, so it is spelled out here. It is fixed by the BAM specification, which
/// is why a 28-byte literal is the right way to hold it rather than something derived.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Whether `path` is a BAM whose writer finished.
///
/// The only difference between a complete BGZF stream and one whose process was killed partway is
/// this marker, so it is the whole test. It answers a question that matters at a scale where
/// re-deriving the file costs hours: `navigator-app`'s realignment resume uses it to decide whether
/// a previous attempt's intermediate can be picked up or has to be thrown away.
///
/// Cheap by construction — one open and a 28-byte read from the end, regardless of a file that may
/// be 60 GB. It does not validate the records inside; a writer that finished wrote them, and a
/// deeper check would mean reading the whole file, which is the cost this exists to avoid.
pub fn is_complete_bam(path: &Path) -> bool {
    use std::io::{Read, Seek, SeekFrom};

    let Ok(mut file) = File::open(path) else {
        return false;
    };
    let Ok(len) = file.metadata().map(|m| m.len()) else {
        return false;
    };
    if len < BGZF_EOF.len() as u64 {
        return false;
    }
    if file.seek(SeekFrom::End(-(BGZF_EOF.len() as i64))).is_err() {
        return false;
    }

    let mut tail = [0u8; BGZF_EOF.len()];
    file.read_exact(&mut tail).is_ok() && tail == BGZF_EOF
}

/// Open `path` for a sequential read with threaded block inflation.
pub(crate) fn open(path: &Path) -> Result<BamReader, AnalysisError> {
    let file = File::open(path).map_err(|e| AnalysisError::io(path, e))?;
    let inner = bgzf::io::MultithreadedReader::with_worker_count(crate::reader::bgzf_worker_count(), file);
    Ok(bam::io::Reader::from(inner))
}

/// Create `path` for writing with threaded block deflation.
pub(crate) fn create(path: &Path) -> Result<BamWriter, AnalysisError> {
    let file = File::create(path).map_err(|e| AnalysisError::io(path, e))?;
    let inner = bgzf::io::MultithreadedWriter::with_worker_count(
        crate::reader::bgzf_worker_count(),
        BufWriter::with_capacity(WRITE_BUFFER, PacedFile::new(file)),
    );
    Ok(bam::io::Writer::from(inner))
}

/// Finish the stream, including the BGZF end-of-file block.
///
/// A BGZF file without its EOF block is indistinguishable from a truncated one, and the plain
/// writer's `try_finish` does not exist on this type — the equivalent is draining the workers.
///
/// The finished file is then synced, which matters more here than it looks: that EOF block is
/// exactly what [`navigator_app::realign_job`]'s resume uses to tell a stage output it can trust
/// from one a killed job left half-written. A marker still sitting in the page cache would be a
/// promise the disk has not made.
pub(crate) fn finish(mut writer: BamWriter, path: &Path) -> Result<(), AnalysisError> {
    let mut buffered = writer.get_mut().finish().map_err(|e| AnalysisError::io(path, e))?;
    buffered.flush().map_err(|e| AnalysisError::io(path, e))?;
    buffered
        .get_ref()
        .file
        .sync_data()
        .map_err(|e| AnalysisError::io(path, e))
}
