//! A read of an alignment that does not depend on the format. The walkers, which are coverage,
//! the caller and read-metrics, must read a record the same way from a BAM and from a CRAM.
//!
//! But `noodles` gives two different families of reader. A BAM gives a borrowed `bam::Record`. A
//! CRAM gives an owned `sam::alignment::RecordBuf`, and it needs the reference FASTA to decode.
//! This module brings both to a `RecordBuf`. That is one owned allocation at each record, which is
//! the same order of cost that a CRAM pays in any case. So the hot loops over the bases do not
//! know the format, and they allocate nothing.
//!
//! noodles stays inside this crate, and that is deliberate. See lib.rs. This is the one place that
//! knows about the repository of reference sequences of a CRAM.

use std::fs::File;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use noodles::core::region::Interval;
use noodles::core::{Position, Region};
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, sam};

/// The count of worker threads that decompress bgzf, for a sequential read of a BAM.
///
/// bgzf is a stream of gzip blocks. So the inflation of those blocks runs in parallel, while the
/// parse of the records stays sequential. The output does not change, to the last byte, because
/// the threads only decompress.
///
/// The default is the available parallelism, less one for the consumer that parses the records,
/// and up to 6. Above a few inflate workers, the one consumer thread is the limit.
/// `NAVIGATOR_BGZF_THREADS` overrides it, clamped to 1 or more. Set it to 1 to use one thread.
pub(crate) fn bgzf_worker_count() -> NonZeroUsize {
    if let Some(n) = std::env::var("NAVIGATOR_BGZF_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
    {
        return NonZeroUsize::new(n.max(1)).unwrap();
    }
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    NonZeroUsize::new(cores.saturating_sub(1).clamp(1, 6)).unwrap()
}

use crate::cancel::CancelToken;
use crate::error::AnalysisError;
use crate::readview::{AlnRead, SeqRecord};

/// The stack size of each thread, in bytes, for any thread that decodes a BAM or **CRAM** record.
///
/// The CRAM decoder of noodles recurses in proportion to the data. The CRAM **3.1** codecs are the
/// ones that matter: the range and arithmetic coder, fqzcomp, and the name tokenizer. An older 3.0
/// file never uses those. The decoder can recurse deep enough to overflow a default thread stack
/// of 2 MiB, and even the pools of rayon.
///
/// A stack overflow **aborts the process**. It is not a panic that the code can catch. One file
/// with a deep encoding would else take down the whole app or batch. So give a decode thread a
/// large stack. `NAVIGATOR_DECODE_STACK_MB` overrides the size, in whole MiB, clamped to 8 or
/// more.
pub fn decode_stack_size() -> usize {
    let mb = std::env::var("NAVIGATOR_DECODE_STACK_MB")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(64)
        .max(8);
    mb * 1024 * 1024
}

/// Build a rayon pool whose worker threads have a stack that is safe for a decode. See
/// [`decode_stack_size`]. Use it for any parallel work that decodes a CRAM or BAM record. The
/// rayon default of 2 MiB is not enough for a CRAM 3.1 file with a deep encoding, and neither is a
/// small fixed increase.
pub fn decode_pool(threads: usize) -> Result<rayon::ThreadPool, AnalysisError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(decode_stack_size())
        .build()
        .map_err(|e| AnalysisError::Message(format!("thread pool: {e}")))
}

/// On-disk alignment container, by extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Bam,
    Cram,
}

/// Detect the alignment format from the path extension (`.cram` → CRAM, else BAM).
pub fn detect_format(path: &Path) -> Format {
    match path.extension().and_then(|e| e.to_str()) {
        Some(e) if e.eq_ignore_ascii_case("cram") => Format::Cram,
        _ => Format::Bam,
    }
}

/// Build a caching FASTA sequence repository from an indexed reference (needs a `.fai`).
/// Required to decode CRAM; reused for any reference-backed reading.
pub fn build_repository(reference: &Path) -> Result<fasta::Repository, AnalysisError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .map_err(|e| AnalysisError::io(reference, e))?;
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(
        reader,
    )))
}

/// CRAM needs a reference; surface a clear error if one was not supplied.
fn require_reference<'a>(path: &Path, reference: Option<&'a Path>) -> Result<&'a Path, AnalysisError> {
    reference.ok_or_else(|| AnalysisError::Message(format!("CRAM {} requires a reference FASTA", path.display())))
}

// ---- sequential (whole-file) reading --------------------------------------

/// A reader over a whole BAM or CRAM file. Hold it, and call [`SeqReader::records`]. The BAM path
/// uses a bgzf reader with more than one thread. So the block decompression runs on a worker pool,
/// while the parse of the records stays sequential. See [`bgzf_worker_count`].
pub enum SeqReader {
    Bam {
        inner: bam::io::Reader<bgzf::io::MultithreadedReader<File>>,
        path: PathBuf,
    },
    Cram {
        inner: cram::io::Reader<File>,
        path: PathBuf,
    },
}

/// Open `path` for a sequential pass. Returns the header and the reader. A CRAM needs
/// `reference`, and a BAM ignores it.
pub fn open_seq(path: &Path, reference: Option<&Path>) -> Result<(sam::Header, SeqReader), AnalysisError> {
    match detect_format(path) {
        Format::Bam => {
            let file = File::open(path).map_err(|e| AnalysisError::io(path, e))?;
            let mt = bgzf::io::MultithreadedReader::with_worker_count(bgzf_worker_count(), file);
            let mut inner = bam::io::Reader::from(mt);
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((
                header,
                SeqReader::Bam {
                    inner,
                    path: path.to_path_buf(),
                },
            ))
        }
        Format::Cram => {
            let repo = build_repository(require_reference(path, reference)?)?;
            let mut inner = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(path)
                .map_err(|e| AnalysisError::io(path, e))?;
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((
                header,
                SeqReader::Cram {
                    inner,
                    path: path.to_path_buf(),
                },
            ))
        }
    }
}

impl SeqReader {
    /// Iterate every record as a `RecordBuf`. `header` must be the one returned by
    /// [`open_seq`] (CRAM decodes against it; BAM converts its records through it).
    pub fn records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = Result<RecordBuf, AnalysisError>> + 'a> {
        match self {
            SeqReader::Bam { inner, path } => {
                let path = path.clone();
                Box::new(inner.records().map(move |r| {
                    let rec = r.map_err(|e| AnalysisError::io(&path, e))?;
                    RecordBuf::try_from_alignment_record(header, &rec).map_err(|e| AnalysisError::io(&path, e))
                }))
            }
            SeqReader::Cram { inner, path } => {
                let path = path.clone();
                Box::new(
                    inner
                        .records(header)
                        .map(move |r| r.map_err(|e| AnalysisError::io(&path, e))),
                )
            }
        }
    }

    /// Walk every record as a [`SeqRecord`]. This is the **lazy** counterpart to
    /// [`SeqReader::records`].
    ///
    /// The BAM path gives the zero-copy `bam::Record`. It does no decode into an owned
    /// `RecordBuf`, and it parses no tag, and that is the gain on the hot path. The CRAM path gives
    /// the decoded `RecordBuf`, because there is no cheaper form. The walkers take an
    /// `&impl AlnRead`, so a `SeqRecord` drives them with no allocation on the BAM path.
    pub fn records_lazy<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = Result<SeqRecord, AnalysisError>> + 'a> {
        match self {
            SeqReader::Bam { inner, path } => {
                let path = path.clone();
                Box::new(
                    inner
                        .records()
                        .map(move |r| r.map(SeqRecord::Bam).map_err(|e| AnalysisError::io(&path, e))),
                )
            }
            SeqReader::Cram { inner, path } => {
                let path = path.clone();
                Box::new(
                    inner
                        .records(header)
                        .map(move |r| r.map(SeqRecord::Cram).map_err(|e| AnalysisError::io(&path, e))),
                )
            }
        }
    }
}

// ---- indexed (region) reading ---------------------------------------------

/// An indexed reader over BAM or CRAM. Hold it and call [`IdxReader::query`].
pub enum IdxReader {
    Bam {
        inner: bam::io::IndexedReader<bgzf::io::Reader<File>>,
        path: PathBuf,
    },
    Cram {
        inner: cram::io::IndexedReader<File>,
        repo: fasta::Repository,
        path: PathBuf,
    },
}

/// The file offsets of the `.crai` containers that can hold a record inside `interval` on
/// `ref_id`.
///
/// **This is the whole reason that a CRAM region query is usable.** A CRAM container is the unit
/// of a decode, and you can not decode part of one. To limit *which containers the code decodes*
/// is the only place where a region query can save work.
///
/// The `Query` of noodles, and our own `for_each` before this code, chose the containers by
/// reference sequence alone. They then threw away the records outside the region, *after* they
/// decoded those records. So every query cost a whole chromosome, at any size of region.
/// A measurement gave 20.9 s for a 1 bp query on chr21, and 116 s on chr1. The same query on a BAM
/// takes 4 to 6 ms. chr21 of a 30x WGS holds 1,140 containers, and a point query needs exactly one
/// of them.
///
/// The code **keeps** a container whose `alignment_start` is absent. This index can not place such
/// a container, and to drop it would lose records where nobody sees it happen. The code skips a
/// container only on positive evidence that the container lies outside the interval.
fn cram_container_offsets(index: &cram::crai::Index, ref_id: usize, interval: Interval) -> Vec<u64> {
    index
        .iter()
        .filter(|r| r.reference_sequence_id() == Some(ref_id))
        .filter(|r| match r.alignment_start() {
            Some(start) => {
                // Span 0 would make an empty range, which intersects nothing; treat it as one base.
                let span = r.alignment_span().max(1);
                let end = Position::new(usize::from(start).saturating_add(span - 1)).unwrap_or(start);
                interval.intersects((start..=end).into())
            }
            None => true,
        })
        .map(|r| r.offset())
        .collect()
}

/// Open `path` for a region query over its index. It loads the `.bai` or `.crai` itself. A CRAM
/// needs `reference`.
pub fn open_indexed(path: &Path, reference: Option<&Path>) -> Result<(sam::Header, IdxReader), AnalysisError> {
    match detect_format(path) {
        Format::Bam => {
            let mut inner = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .map_err(|e| AnalysisError::io(path, e))?;
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((
                header,
                IdxReader::Bam {
                    inner,
                    path: path.to_path_buf(),
                },
            ))
        }
        Format::Cram => {
            let repo = build_repository(require_reference(path, reference)?)?;
            let mut inner = cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repo.clone())
                .build_from_path(path)
                .map_err(|e| AnalysisError::io(path, e))?;
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((
                header,
                IdxReader::Cram {
                    inner,
                    repo,
                    path: path.to_path_buf(),
                },
            ))
        }
    }
}

impl IdxReader {
    /// Walk the records inside `region` as `RecordBuf` values.
    pub fn query<'a>(
        &'a mut self,
        header: &'a sam::Header,
        region: &Region,
    ) -> Result<Box<dyn Iterator<Item = Result<RecordBuf, AnalysisError>> + 'a>, AnalysisError> {
        match self {
            IdxReader::Bam { inner, path } => {
                let path = path.clone();
                let q = inner.query(header, region).map_err(|e| AnalysisError::io(&path, e))?;
                Ok(Box::new(q.records().map(move |r| {
                    let rec = r.map_err(|e| AnalysisError::io(&path, e))?;
                    RecordBuf::try_from_alignment_record(header, &rec).map_err(|e| AnalysisError::io(&path, e))
                })))
            }
            // This code does the work itself, and it does not call `inner.query(...)`. The
            // `Query` of noodles decodes every container of the contig, and it filters the records
            // after that. A 1 bp query then costs a whole chromosome. See
            // [`cram_container_offsets`].
            //
            // This code decodes only the containers that can hold a record in the region, and it
            // does that lazily, one container at a time. So a caller that stops early, such as a
            // `.take(n)` probe or a walk that somebody cancelled, does not pay for the rest.
            IdxReader::Cram { inner, repo, path } => {
                use std::io::{Seek, SeekFrom};

                use noodles::sam::alignment::Record as _; // alignment_start/_end on cram::Record

                let path = path.clone();
                let repo = repo.clone();
                let ref_id = header
                    .reference_sequences()
                    .get_index_of(region.name())
                    .ok_or_else(|| {
                        AnalysisError::Message(format!(
                            "contig {} not in {} header",
                            String::from_utf8_lossy(region.name()),
                            path.display()
                        ))
                    })?;
                let interval = region.interval();
                let mut offsets = cram_container_offsets(inner.index(), ref_id, interval).into_iter();

                let mut pending: std::vec::IntoIter<RecordBuf> = Vec::new().into_iter();
                let mut container = cram::io::reader::Container::default();
                Ok(Box::new(std::iter::from_fn(move || {
                    loop {
                        if let Some(rec) = pending.next() {
                            return Some(Ok(rec));
                        }
                        // Next container that can overlap; `None` ends the iterator.
                        let offset = offsets.next()?;
                        let io_err = |e| AnalysisError::io(&path, e);
                        if let Err(e) = inner.get_mut().seek(SeekFrom::Start(offset)).map_err(io_err) {
                            return Some(Err(e));
                        }
                        match inner.read_container(&mut container).map_err(io_err) {
                            Ok(0) => continue,
                            Ok(_) => {}
                            Err(e) => return Some(Err(e)),
                        }
                        let compression_header = match container.compression_header().map_err(io_err) {
                            Ok(h) => h,
                            Err(e) => return Some(Err(e)),
                        };
                        let mut buf = Vec::new();
                        for slice in container.slices() {
                            let slice = match slice.map_err(io_err) {
                                Ok(s) => s,
                                Err(e) => return Some(Err(e)),
                            };
                            let (core, external) = match slice.decode_blocks().map_err(io_err) {
                                Ok(b) => b,
                                Err(e) => return Some(Err(e)),
                            };
                            let records = match slice
                                .records(repo.clone(), header, &compression_header, &core, &external)
                                .map_err(io_err)
                            {
                                Ok(r) => r,
                                Err(e) => return Some(Err(e)),
                            };
                            for rec in &records {
                                // The same test at each record that noodles applies after its
                                // decode. The container filter is a coarse first pass, and it does
                                // not replace this test.
                                if let (Some(Ok(start)), Some(Ok(end))) = (rec.alignment_start(), rec.alignment_end()) {
                                    if !interval.intersects((start..=end).into()) {
                                        continue;
                                    }
                                } else {
                                    continue;
                                }
                                match RecordBuf::try_from_alignment_record(header, rec).map_err(io_err) {
                                    Ok(r) => buf.push(r),
                                    Err(e) => return Some(Err(e)),
                                }
                            }
                        }
                        pending = buf.into_iter();
                    }
                })))
            }
        }
    }

    /// Walk the unmapped records that have no place, which are the tail of a BAM, as `RecordBuf`
    /// values. This works for a BAM alone. The `.crai` of a CRAM gives no query for the unmapped
    /// records, so this returns an error there. A caller that needs the unmapped tail of a CRAM
    /// must take a sequential pass.
    pub fn query_unmapped<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Result<Box<dyn Iterator<Item = Result<RecordBuf, AnalysisError>> + 'a>, AnalysisError> {
        match self {
            IdxReader::Bam { inner, path } => {
                let path = path.clone();
                let q = inner.query_unmapped().map_err(|e| AnalysisError::io(&path, e))?;
                Ok(Box::new(q.map(move |r| {
                    let rec = r.map_err(|e| AnalysisError::io(&path, e))?;
                    RecordBuf::try_from_alignment_record(header, &rec).map_err(|e| AnalysisError::io(&path, e))
                })))
            }
            IdxReader::Cram { path, .. } => Err(AnalysisError::Message(format!(
                "unmapped-record query unsupported for CRAM {}",
                path.display()
            ))),
        }
    }
}

/// A consumer of one record, which the indexed reader drives over a region.
///
/// The `accept` method is generic over [`AlnRead`], so the compiler makes one copy for each record
/// type. The BAM path gives it the **lazy, zero-copy** `bam::Record`. There is no owned
/// `RecordBuf` allocation at each read, and that is the gain on the hot path. The CRAM path gives
/// it the decoded `RecordBuf`. One sink serves both.
pub trait RecordSink {
    fn accept(&mut self, record: &impl AlnRead);
}

/// How often a record loop polls the cancel token. The value keeps the check out of a profile,
/// and it keeps the worst delay between a click and a stop well below one frame. At about 1M
/// records/s that is a check every few milliseconds. The check itself is one relaxed atomic
/// load.
const CANCEL_CHECK_RECORDS: u32 = 4096;

impl IdxReader {
    /// Drive `sink` over every record inside `region`. A BAM gives a lazy record, and a CRAM gives
    /// a `RecordBuf`. A record that the code can not read stops the walk with an error. This is the
    /// counterpart of [`IdxReader::query`] that allocates nothing, where that one copies each
    /// record into an owned `RecordBuf`.
    ///
    /// The loop polls `cancel` every [`CANCEL_CHECK_RECORDS`] records. So a walk that somebody
    /// cancelled stops in the middle of a contig, and not at the next contig boundary. On
    /// chr1 that is the difference between a stop in milliseconds and a stop in minutes. Pass
    /// [`CancelToken::none`] when there is nothing to cancel.
    pub fn for_each<S: RecordSink>(
        &mut self,
        header: &sam::Header,
        region: &Region,
        sink: &mut S,
        cancel: &CancelToken,
    ) -> Result<(), AnalysisError> {
        match self {
            IdxReader::Bam { inner, path } => {
                let path = path.clone();
                let q = inner.query(header, region).map_err(|e| AnalysisError::io(&path, e))?;
                let mut seen = 0u32;
                for r in q.records() {
                    sink.accept(&r.map_err(|e| AnalysisError::io(&path, e))?);
                    seen += 1;
                    if seen % CANCEL_CHECK_RECORDS == 0 {
                        cancel.check()?;
                    }
                }
                Ok(())
            }
            IdxReader::Cram { inner, repo, path } => {
                // Decode the CRAM containers of the region down to borrowed `cram::Record`
                // values, and drive the sink from those directly. That leaves out the `RecordBuf`
                // copy at each read, which the high-level `query` iterator pays. On a 30x WGS CRAM
                // that copy costs about 1.74 times the decode of one read.
                //
                // This has the same shape as the `Query` of noodles. Seek each `.crai` container
                // whose reference matches, decode the slices of that container, and keep the
                // records inside the query interval.
                use std::io::{Seek, SeekFrom};

                use noodles::sam::alignment::Record as _; // alignment_start/_end on cram::Record

                use crate::readview::CramRead;

                let path = path.clone();
                let repo = repo.clone();
                let io_err = |e| AnalysisError::io(&path, e);

                // Resolve the query contig to its @SQ index, and capture the query interval.
                let ref_id = header
                    .reference_sequences()
                    .get_index_of(region.name())
                    .ok_or_else(|| {
                        AnalysisError::Message(format!(
                            "contig {} not in {} header",
                            String::from_utf8_lossy(region.name()),
                            path.display()
                        ))
                    })?;
                let interval = region.interval();

                // Collect the file offsets of the containers that can hold a record in the query.
                // Do that before the code takes a mutable borrow of `inner` to seek and read. The
                // borrow of the `.crai` index can not overlap the borrow for the read.
                //
                // The choice goes on the interval, and not on the contig alone. That is what keeps
                // the cost proportional to the region, and not to the chromosome. See
                // [`cram_container_offsets`].
                let offsets = cram_container_offsets(inner.index(), ref_id, interval);

                let mut container = cram::io::reader::Container::default();
                for offset in offsets {
                    // The check goes at each container, and not at each record. The code decodes a
                    // CRAM container as a unit, so this is the smallest step at which a stop saves
                    // work.
                    cancel.check()?;
                    inner.get_mut().seek(SeekFrom::Start(offset)).map_err(io_err)?;
                    if inner.read_container(&mut container).map_err(io_err)? == 0 {
                        continue;
                    }
                    let compression_header = container.compression_header().map_err(io_err)?;
                    for slice in container.slices() {
                        let slice = slice.map_err(io_err)?;
                        let (core, external) = slice.decode_blocks().map_err(io_err)?;
                        let records = slice
                            .records(repo.clone(), header, &compression_header, &core, &external)
                            .map_err(io_err)?;
                        for rec in &records {
                            // Same overlap test noodles' `Query` applies post-decode.
                            if let (Some(Ok(start)), Some(Ok(end))) = (rec.alignment_start(), rec.alignment_end()) {
                                if interval.intersects((start..=end).into()) {
                                    sink.accept(&CramRead { rec, header });
                                }
                            }
                        }
                    }
                }
                Ok(())
            }
        }
    }

    /// Drive `sink` over the unplaced unmapped records (BAM only; CRAM errors, as in
    /// [`IdxReader::query_unmapped`]).
    pub fn for_each_unmapped<S: RecordSink>(
        &mut self,
        sink: &mut S,
        cancel: &CancelToken,
    ) -> Result<(), AnalysisError> {
        match self {
            IdxReader::Bam { inner, path } => {
                let path = path.clone();
                let q = inner.query_unmapped().map_err(|e| AnalysisError::io(&path, e))?;
                let mut seen = 0u32;
                for r in q {
                    sink.accept(&r.map_err(|e| AnalysisError::io(&path, e))?);
                    seen += 1;
                    if seen % CANCEL_CHECK_RECORDS == 0 {
                        cancel.check()?;
                    }
                }
                Ok(())
            }
            IdxReader::Cram { path, .. } => Err(AnalysisError::Message(format!(
                "unmapped-record query unsupported for CRAM {}",
                path.display()
            ))),
        }
    }
}

/// True when a BAM index sits beside `path`, as `foo.bam.bai` or as `foo.bai`. The parallel walker
/// over the contigs needs one for its region queries. A caller falls back to a sequential pass
/// when this is false. A CRAM is out, because its `.crai` gives no query for the unmapped
/// records.
pub fn has_bai_index(path: &Path) -> bool {
    if detect_format(path) != Format::Bam {
        return false;
    }
    let dotted = path.with_extension("bam.bai"); // foo.bam -> foo.bam.bai
    let replaced = path.with_extension("bai"); // foo.bam -> foo.bai
    dotted.exists() || replaced.exists()
}

/// Whether a CRAM `.crai` coordinate index is present (`foo.cram.crai` or `foo.crai`).
pub fn has_crai_index(path: &Path) -> bool {
    if detect_format(path) != Format::Cram {
        return false;
    }
    path.with_extension("cram.crai").exists() || path.with_extension("crai").exists()
}

/// True when the file has a coordinate index that supports a **region query on one contig**. That
/// is a `.bai` for a BAM, or a `.crai` for a CRAM. The parallel walker over the contigs needs this.
/// A CRAM also has no region query for its unmapped tail, and a caller handles that on its own.
pub fn has_region_index(path: &Path) -> bool {
    has_bai_index(path) || has_crai_index(path)
}

// ---- header-only ----------------------------------------------------------

/// Read the SAM header alone, for example to find the length of a contig. A CRAM needs
/// `reference`.
pub fn read_header(path: &Path, reference: Option<&Path>) -> Result<sam::Header, AnalysisError> {
    open_seq(path, reference).map(|(header, _)| header)
}

/// The names of the reference sequences of the alignment, which are the contigs, in header order.
/// A CRAM needs `reference`.
///
/// Use it to reconcile the contig of a panel or a site against the names that the file uses. A
/// GRCh37 alignment can use a bare `1` where a panel locus holds `chr1`, and the two can also be
/// the other way round.
pub fn contig_names(path: &Path, reference: Option<&Path>) -> Result<Vec<String>, AnalysisError> {
    let header = read_header(path, reference)?;
    Ok(header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect())
}

/// Read one contig's full sequence from an indexed FASTA (needs a `.fai`). Used to pull
/// `chrM` out of a reference for the rCRS↔chrM mtDNA coordinate map.
pub fn read_contig_sequence(reference: &Path, contig: &str) -> Result<Vec<u8>, AnalysisError> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .map_err(|e| AnalysisError::io(reference, e))?;
    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    let record = reader.query(&region).map_err(|e| AnalysisError::io(reference, e))?;
    Ok(record.sequence().as_ref().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_cram_by_extension() {
        assert_eq!(detect_format(Path::new("x/HG00096.chm13.cram")), Format::Cram);
        assert_eq!(detect_format(Path::new("x/HG00096.CRAM")), Format::Cram);
        assert_eq!(detect_format(Path::new("x/sample.bam")), Format::Bam);
        assert_eq!(detect_format(Path::new("x/sample")), Format::Bam);
    }

    /// The fields a walker reads off a record, captured comparably from any [`AlnRead`].
    #[derive(Debug, PartialEq)]
    struct Captured {
        flags: u16,
        start: Option<usize>,
        mate_start: Option<usize>,
        ref_id: Option<usize>,
        mate_ref_id: Option<usize>,
        mapq: Option<u8>,
        tlen: i32,
        seq_len: usize,
        quals: Vec<u8>,
        cigar: Vec<(u8, usize)>,
    }

    fn capture(r: &impl AlnRead) -> Captured {
        let (quals, cigar) = r.pileup_with(|q, ops| (q.to_vec(), ops.map(|(k, l)| (k as u8, l)).collect::<Vec<_>>()));
        Captured {
            flags: r.flags().bits(),
            start: r.alignment_start(),
            mate_start: r.mate_alignment_start(),
            ref_id: r.reference_sequence_id(),
            mate_ref_id: r.mate_reference_sequence_id(),
            mapq: r.mapping_quality(),
            tlen: r.template_length(),
            seq_len: r.sequence_len(),
            quals,
            cigar,
        }
    }

    /// The new CRAM `for_each` path works at the level of a slice, with a borrowed `cram::Record`.
    /// It must give records whose fields all match those from the high-level `query` path, which
    /// gives an owned `RecordBuf`. This test guards our copy of the internal noodles API, which is
    /// the crai seek and the slice decode, against a change of version.
    #[test]
    fn cram_for_each_matches_query_recordbuf() {
        struct CollectSink(Vec<Captured>);
        impl RecordSink for CollectSink {
            fn accept(&mut self, record: &impl AlnRead) {
                self.0.push(capture(record));
            }
        }

        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
        let cram = dir.join("coverage.cram");
        let reference = dir.join("ref.fa");
        let region = Region::new(b"chrM".to_vec(), ..);

        // New path: for_each over borrowed cram::Record.
        let (header, mut idx) = open_indexed(&cram, Some(&reference)).expect("open");
        let mut sink = CollectSink(Vec::new());
        idx.for_each(&header, &region, &mut sink, &CancelToken::none())
            .expect("for_each");

        // Old path: query yields RecordBuf.
        let (header2, mut idx2) = open_indexed(&cram, Some(&reference)).expect("open2");
        let via_query: Vec<Captured> = idx2
            .query(&header2, &region)
            .expect("query")
            .map(|r| capture(&r.expect("rec")))
            .collect();

        assert!(!sink.0.is_empty(), "fixture should have chrM records");
        assert_eq!(sink.0, via_query, "cram::Record path must match RecordBuf path");
    }

    /// [`cram_container_offsets`] decides which containers the code decodes at all. Its behaviour
    /// at a boundary *is* the correctness of every CRAM region query. A container that it skips
    /// wrongly means reads that a variant call never sees, and no later test would trace that back
    /// to the reader. This test looks at the edges, where an off-by-one error lives.
    #[test]
    fn container_offsets_select_only_overlapping_containers() {
        let p = |n: usize| Position::new(n).unwrap();
        // Three containers on ref 0, which cover [1000,1099], [2000,2099] and [3000,3099]. One
        // container on ref 1. And one that the index can not place.
        let idx: cram::crai::Index = vec![
            cram::crai::Record::new(Some(0), Some(p(1000)), 100, 10, 0, 0),
            cram::crai::Record::new(Some(0), Some(p(2000)), 100, 20, 0, 0),
            cram::crai::Record::new(Some(0), Some(p(3000)), 100, 30, 0, 0),
            cram::crai::Record::new(Some(1), Some(p(2000)), 100, 40, 0, 0),
            cram::crai::Record::new(Some(0), None, 0, 50, 0, 0),
        ];
        let sel = |a: usize, b: usize| cram_container_offsets(&idx, 0, (p(a)..=p(b)).into());

        // A point inside one container decodes that container, and not the contig. This one
        // assertion is the difference between 8 ms and 21 s on a real chr21.
        assert_eq!(sel(2050, 2050), vec![20, 50]);
        // The boundaries. A query that reaches the first base or the last base of a container
        // counts as an overlap.
        assert_eq!(sel(2099, 2099), vec![20, 50], "last base of a container overlaps");
        assert_eq!(sel(2000, 2000), vec![20, 50], "first base of a container overlaps");
        assert_eq!(sel(2100, 2100), vec![50], "one past the end does not");
        assert_eq!(sel(1999, 1999), vec![50], "one before the start does not");
        // A span that crosses more than one container takes exactly the containers that it
        // crosses.
        assert_eq!(sel(1050, 2050), vec![10, 20, 50]);
        // The other reference is never selected, even at identical coordinates.
        assert_eq!(cram_container_offsets(&idx, 1, (p(2050)..=p(2050)).into()), vec![40]);
        // An unbounded interval keeps every container on the reference.
        assert_eq!(
            cram_container_offsets(&idx, 0, Region::new(b"x".to_vec(), ..).interval()),
            vec![10, 20, 30, 50]
        );
    }

    /// The `query` in this module must return exactly what the `Query` of noodles returns. We
    /// wrote our own for speed. A version of it that drops records where nobody looks would seem
    /// to be a faster caller, and not a broken one. That failure needs a test.
    #[test]
    fn cram_query_matches_noodles_query() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
        let cram = dir.join("coverage.cram");
        let reference = dir.join("ref.fa");

        for region in [
            Region::new(b"chrM".to_vec(), ..),
            Region::new(
                b"chrM".to_vec(),
                Position::new(1).unwrap()..=Position::new(200).unwrap(),
            ),
            Region::new(
                b"chrM".to_vec(),
                Position::new(50).unwrap()..=Position::new(60).unwrap(),
            ),
        ] {
            let (header, mut ours) = open_indexed(&cram, Some(&reference)).expect("open");
            let mine: Vec<Captured> = ours
                .query(&header, &region)
                .expect("query")
                .map(|r| capture(&r.expect("rec")))
                .collect();

            // noodles' own indexed query, unmodified, as the oracle.
            let repo = build_repository(&reference).expect("repo");
            let mut theirs = cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram)
                .expect("noodles open");
            let nheader = theirs.read_header().expect("noodles header");
            let reference_impl: Vec<Captured> = theirs
                .query(&nheader, &region)
                .expect("noodles query")
                .map(|r| capture(&r.expect("rec")))
                .collect();

            assert_eq!(
                mine, reference_impl,
                "region {region:?}: must match noodles' Query exactly"
            );
        }
    }
}
