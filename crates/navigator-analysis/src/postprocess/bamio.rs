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
use std::io::BufWriter;
use std::path::Path;

use noodles::{bam, bgzf};

use crate::error::AnalysisError;

/// A BAM reader whose block decompression runs on a worker pool.
pub(crate) type BamReader = bam::io::Reader<bgzf::io::MultithreadedReader<File>>;

/// A BAM writer whose block compression runs on a worker pool.
pub(crate) type BamWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<File>>>;

const WRITE_BUFFER: usize = 1 << 20;

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
        BufWriter::with_capacity(WRITE_BUFFER, file),
    );
    Ok(bam::io::Writer::from(inner))
}

/// Finish the stream, including the BGZF end-of-file block.
///
/// A BGZF file without its EOF block is indistinguishable from a truncated one, and the plain
/// writer's `try_finish` does not exist on this type — the equivalent is draining the workers.
pub(crate) fn finish(mut writer: BamWriter, path: &Path) -> Result<(), AnalysisError> {
    writer
        .get_mut()
        .finish()
        .map(|_| ())
        .map_err(|e| AnalysisError::io(path, e))
}
