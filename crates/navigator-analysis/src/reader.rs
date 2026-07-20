//! Format-agnostic alignment reading. The walkers (coverage, caller, read-metrics) need
//! to read records the same way whether the file is BAM or CRAM, but `noodles` exposes
//! two different reader families: BAM yields borrowed `bam::Record`s, CRAM yields owned
//! `sam::alignment::RecordBuf`s and needs the reference FASTA to decode. This module
//! normalizes both to `RecordBuf` (one owned allocation per record — the same order CRAM
//! pays anyway) so the hot per-base loops stay format-blind and allocation-free.
//!
//! noodles is intentionally confined to this crate (see lib.rs); this is the single place
//! that knows about CRAM's reference-sequence repository.

use std::fs::File;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use noodles::core::Region;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, sam};

/// Worker threads for multithreaded bgzf decompression of BAM sequential reads. bgzf is a
/// block-gzip stream, so block inflation parallelizes while record parsing stays sequential
/// (output is byte-identical — only decompression is threaded). Defaults to the available
/// parallelism minus one (the record-parsing consumer), capped at 6 — beyond a handful of
/// inflate workers the single consumer thread is the limit. Override with
/// `NAVIGATOR_BGZF_THREADS` (clamped to >= 1; set to 1 to disable threading).
fn bgzf_worker_count() -> NonZeroUsize {
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

/// Per-thread stack size (bytes) for any thread that decodes a BAM/**CRAM** record. noodles' CRAM
/// decoder recurses proportionally to the data — notably the CRAM **3.1** codecs (range/arithmetic
/// coder, fqzcomp, name tokenizer), which older 3.0 files never exercise — and can recurse deep
/// enough to blow a default thread stack (2 MiB) or even rayon's pools. A stack overflow **aborts
/// the process** (it is not a catchable panic), so a single deeply-encoded file would otherwise take
/// down the whole app/batch. Give decode threads a generous stack. Override with
/// `NAVIGATOR_DECODE_STACK_MB` (whole MiB; clamped to >= 8).
pub fn decode_stack_size() -> usize {
    let mb = std::env::var("NAVIGATOR_DECODE_STACK_MB")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(64)
        .max(8);
    mb * 1024 * 1024
}

/// Build a rayon pool whose worker threads have a decode-safe stack ([`decode_stack_size`]).
/// Use this for any parallel work that decodes CRAM/BAM records — the rayon default (2 MiB) and
/// even a modest fixed bump are not enough for deeply-encoded CRAM 3.1 files.
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

/// CRAM needs a reference; surface a clear error if one wasn't supplied.
fn require_reference<'a>(path: &Path, reference: Option<&'a Path>) -> Result<&'a Path, AnalysisError> {
    reference.ok_or_else(|| AnalysisError::Message(format!("CRAM {} requires a reference FASTA", path.display())))
}

// ---- sequential (whole-file) reading --------------------------------------

/// A whole-file reader over BAM or CRAM. Hold it and call [`SeqReader::records`]. The BAM
/// path uses a multithreaded bgzf reader so block decompression runs on a worker pool while
/// records are parsed sequentially (see [`bgzf_worker_count`]).
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

/// Open `path` for a sequential pass, returning the header and reader. `reference` is
/// required for CRAM (ignored for BAM).
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

    /// Iterate every record as a [`SeqRecord`] — the **lazy** counterpart to [`SeqReader::records`].
    /// The BAM path yields the zero-copy `bam::Record` (no owned `RecordBuf` decode/tag-parse, the
    /// hot-path win); the CRAM path yields the decoded `RecordBuf` (no cheaper form). The walkers
    /// consume `&impl AlnRead`, so `SeqRecord` drives them with no allocation on the BAM path.
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

/// Open `path` for indexed region queries (autoloads the `.bai`/`.crai`). `reference` is
/// required for CRAM.
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
    /// Iterate the records overlapping `region` as `RecordBuf`s.
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
            IdxReader::Cram { inner, path, .. } => {
                let path = path.clone();
                let q = inner.query(header, region).map_err(|e| AnalysisError::io(&path, e))?;
                Ok(Box::new(q.map(move |r| r.map_err(|e| AnalysisError::io(&path, e)))))
            }
        }
    }

    /// Iterate the unplaced unmapped records (the BAM tail) as `RecordBuf`s. BAM only —
    /// CRAM's `.crai` exposes no unmapped query, so it returns an error (callers needing the
    /// unmapped tail for CRAM should take a sequential pass instead).
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

/// A per-record consumer the indexed reader drives over a region. The `accept` method is generic
/// over [`AlnRead`], so it monomorphizes for each record type: the BAM path hands it the **lazy,
/// zero-copy** `bam::Record` (no per-read owned `RecordBuf` allocation — the hot-path win) and the
/// CRAM path hands it the decoded `RecordBuf`. A single sink serves both.
pub trait RecordSink {
    fn accept(&mut self, record: &impl AlnRead);
}

/// How often the record loops poll the cancel token. Tuned to be invisible in a profile while
/// keeping the worst-case delay between a click and a stop well under a frame: at ~1M records/s
/// this is a check every few milliseconds. The check itself is one relaxed atomic load.
const CANCEL_CHECK_RECORDS: u32 = 4096;

impl IdxReader {
    /// Drive `sink` over every record overlapping `region` (BAM: lazy record; CRAM: `RecordBuf`).
    /// A record that fails to read aborts with an error. The allocation-free counterpart to
    /// [`IdxReader::query`] (which copies each record into an owned `RecordBuf`).
    ///
    /// `cancel` is polled every [`CANCEL_CHECK_RECORDS`] records, so a cancelled walk stops
    /// mid-contig instead of at the next contig boundary — on chr1 that is the difference between
    /// stopping in milliseconds and stopping in minutes. Pass [`CancelToken::none`] when there is
    /// nothing to cancel.
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
                // Decode the region's CRAM containers down to borrowed `cram::Record`s and drive the
                // sink off them directly — skipping the per-read `RecordBuf` copy the high-level
                // `query` iterator pays (~1.74× the per-read decode on a 30× WGS CRAM). This mirrors
                // noodles' own `Query`: seek each `.crai` container whose reference matches, decode
                // its slices, and keep the records overlapping the query interval.
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

                // Collect the file offsets of this contig's containers before borrowing `inner`
                // mutably to seek/read (the `.crai` index borrow can't overlap the read borrow).
                let offsets: Vec<u64> = inner
                    .index()
                    .iter()
                    .filter(|r| r.reference_sequence_id() == Some(ref_id))
                    .map(|r| r.offset())
                    .collect();

                let mut container = cram::io::reader::Container::default();
                for offset in offsets {
                    // Per container rather than per record: a CRAM container is decoded as a unit,
                    // so this is the finest granularity at which stopping actually saves work.
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

/// Whether a sibling BAM index (`.bai`, as `foo.bam.bai` or `foo.bai`) exists for `path`.
/// The per-contig parallel walker needs one for region queries; callers fall back to a
/// sequential pass when this is false. CRAM is excluded (its `.crai` has no unmapped query).
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

/// Whether the file has a coordinate index supporting **per-contig region queries** — a BAM `.bai`
/// or a CRAM `.crai`. The prerequisite for the parallel per-contig walker (CRAM additionally can't
/// region-query the unmapped tail; callers handle that separately).
pub fn has_region_index(path: &Path) -> bool {
    has_bai_index(path) || has_crai_index(path)
}

// ---- header-only ----------------------------------------------------------

/// Read just the SAM header (e.g. to resolve a contig length). `reference` is required
/// for CRAM.
pub fn read_header(path: &Path, reference: Option<&Path>) -> Result<sam::Header, AnalysisError> {
    open_seq(path, reference).map(|(header, _)| header)
}

/// The alignment's reference-sequence (contig) names, in header order. `reference` is required for
/// a CRAM. Used to reconcile a panel/site contig against the file's naming convention — a GRCh37
/// alignment may use bare `1` where a panel locus stores `chr1` (or vice versa).
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

    /// The new slice-level CRAM `for_each` path (borrowed `cram::Record`) must yield records
    /// field-identical to the high-level `query` path (owned `RecordBuf`) — guards the noodles
    /// internal-API replication (crai seek + slice decode) against version drift.
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
}
