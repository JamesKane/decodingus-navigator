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
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use noodles::{bam, bgzf};

use crate::error::AnalysisError;
use crate::resource::PacedFile;

/// A BAM reader whose block decompression runs on a worker pool.
pub(crate) type BamReader = bam::io::Reader<bgzf::io::MultithreadedReader<File>>;

/// A BAM reader that inflates on the calling thread.
pub(crate) type PlainBamReader = bam::io::Reader<bgzf::io::Reader<BufReader<File>>>;

/// Read buffer for a stream that is one of many open at once.
const READ_BUFFER: usize = 1 << 18;

/// A BAM writer whose block compression runs on a worker pool.
pub(crate) type BamWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<PacedFile>>>;

const WRITE_BUFFER: usize = 1 << 20;

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

/// Open `path` as one of many streams that will be read concurrently.
///
/// [`open`] is right for a stream that is alone: threading its inflation is the difference between
/// one core and six on a file the stage reads end to end. It is badly wrong for a stream that is
/// one of hundreds. The sort's merge opens every spilled run at once — 688 of them on a 30x WGS at
/// the default budget — and giving each its own worker pool spawned **4,843 threads** in a measured
/// run, each independently reading ahead. The machine did 15,000 IOPS and 6.6 GB/s of disk reads to
/// produce 5 MB/s of merged output, WindowServer could not get scheduled, and the run was killed by
/// the watchdog that noticed.
///
/// There is nothing for a worker pool to do here anyway: the merge is already parallel across runs,
/// and it consumes one record at a time from each. So this inflates on the calling thread, behind a
/// modest buffer — 688 of these cost 688 buffers and no threads at all.
pub(crate) fn open_many(path: &Path) -> Result<PlainBamReader, AnalysisError> {
    let file = File::open(path).map_err(|e| AnalysisError::io(path, e))?;
    let inner = bgzf::io::Reader::new(BufReader::with_capacity(READ_BUFFER, file));
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
    buffered.get_ref().sync().map_err(|e| AnalysisError::io(path, e))
}
