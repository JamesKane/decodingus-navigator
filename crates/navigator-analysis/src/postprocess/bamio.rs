//! Multithreaded BAM I/O for the post-processing stages.
//!
//! Stage C moves the whole alignment through BGZF four times. The sort reads the mapped BAM,
//! writes its spilled runs, reads them back to merge, and writes the sorted output. The duplicate
//! mark and the CRAM output then each read it once more.
//!
//! On a 30x WGS that is hundreds of GB of inflate and deflate, and every pass of it ran on one
//! thread. The first run at WGS scale measured the sort at 4 h 44 m. That made it the stage in the
//! pipeline that cost the most, ahead of the mapping that it feeds.
//!
//! BGZF is a stream of gzip *blocks*. So the compression and the decompression run in parallel
//! across those blocks, while the stream of records stays sequential, and its bytes do not
//! change.
//! That is the same reasoning that [`crate::reader::open_seq`] already applies to a read of a
//! vendor BAM. These stages never got it, because somebody wrote them against the plain
//! constructors.
//!
//! The compression *level* does not change, and that is deliberate. These are all intermediates
//! that one step reads once, so a faster level looks attractive. But a faster level makes the
//! scratch space larger, and the preflight of `navigator-app::realign_job` holds a calibration
//! against that space. That is a separate change, and it carries a new calibration with it.

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use noodles::{bam, bgzf};

use crate::error::AnalysisError;
use navigator_resource::PacedFile;

/// A BAM reader whose block decompression runs on a worker pool.
pub(crate) type BamReader = bam::io::Reader<bgzf::io::MultithreadedReader<File>>;

/// A BAM reader that inflates on the calling thread.
pub(crate) type PlainBamReader = bam::io::Reader<bgzf::io::Reader<BufReader<File>>>;

/// Read buffer for a stream that is one of many open at once.
const READ_BUFFER: usize = 1 << 18;

/// A BAM writer whose block compression runs on a worker pool.
pub(crate) type BamWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<PacedFile>>>;

const WRITE_BUFFER: usize = 1 << 20;

/// The end-of-file marker of BGZF. It is an empty deflate block, and every writer that follows the
/// specification writes it last.
///
/// The code writes it out here, because noodles does not export it. The BAM specification fixes
/// its value, and that is why a literal of 28 bytes is the correct way to hold it. Nothing needs
/// to derive it.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Whether `path` is a BAM whose writer finished.
///
/// This marker is the only difference between a complete BGZF stream and one whose process
/// somebody killed in the middle. So it is the whole test.
///
/// It answers a question that matters at a scale where a new derivation of the file costs hours.
/// The realignment resume in `navigator-app` uses it. That resume decides whether it can take up
/// the intermediate of an earlier try, or must throw that file away.
///
/// It costs almost nothing by construction: one open, and a read of 28 bytes from the end, even
/// on a file of 60 GB. It does not check the records inside. A writer that reached its end wrote
/// them, and a deeper check would mean a read of the whole file. That cost is what this function
/// exists to avoid.
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

/// Open `path` as one of many streams that the code reads at the same time.
///
/// [`open`] is correct for a stream that stands alone. There, threads on its inflation are the
/// difference between one core and six, on a file that the stage reads from end to end. It is very
/// wrong for a stream that is one of hundreds.
///
/// The merge of the sort opens every spilled run at one time. That is 688 of them on a 30x WGS, at
/// the default budget. A worker pool for each one started **4,843 threads** in a measured run, and
/// each of those read ahead on its own. The machine did 15,000 IOPS, and 6.6 GB/s of reads from
/// the disk, to make 5 MB/s of merged output. WindowServer could not get onto a core, and its
/// watchdog saw that and killed the run.
///
/// A worker pool has nothing to do here in any case. The merge already runs in parallel across the
/// runs, and it takes one record at a time from each. So this function inflates on the thread that
/// calls it, behind a buffer of modest size. There, 688 of these cost 688 buffers, and no
/// threads.
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
/// Nothing can separate a BGZF file with no EOF block from one that stops early. This type has no
/// `try_finish`, which the plain writer has. Here the equivalent is to run the workers dry.
///
/// The code then syncs the finished file, and that matters more than it looks. The resume in
/// [`navigator_app::realign_job`] uses that EOF block. It separates a stage output that it can
/// trust from one that a killed job left half written. A marker that still sits in the page cache
/// is a promise that the disk has not made.
pub(crate) fn finish(mut writer: BamWriter, path: &Path) -> Result<(), AnalysisError> {
    let mut buffered = writer.get_mut().finish().map_err(|e| AnalysisError::io(path, e))?;
    buffered.flush().map_err(|e| AnalysisError::io(path, e))?;
    buffered.get_ref().sync().map_err(|e| AnalysisError::io(path, e))
}
