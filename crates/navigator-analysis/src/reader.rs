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
    if let Some(n) = std::env::var("NAVIGATOR_BGZF_THREADS").ok().and_then(|s| s.parse::<usize>().ok()) {
        return NonZeroUsize::new(n.max(1)).unwrap();
    }
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    NonZeroUsize::new(cores.saturating_sub(1).clamp(1, 6)).unwrap()
}

use crate::error::AnalysisError;
use crate::readview::AlnRead;

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
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader)))
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
    Bam { inner: bam::io::Reader<bgzf::MultithreadedReader<File>>, path: PathBuf },
    Cram { inner: cram::io::Reader<File>, path: PathBuf },
}

/// Open `path` for a sequential pass, returning the header and reader. `reference` is
/// required for CRAM (ignored for BAM).
pub fn open_seq(path: &Path, reference: Option<&Path>) -> Result<(sam::Header, SeqReader), AnalysisError> {
    match detect_format(path) {
        Format::Bam => {
            let file = File::open(path).map_err(|e| AnalysisError::io(path, e))?;
            let mt = bgzf::MultithreadedReader::with_worker_count(bgzf_worker_count(), file);
            let mut inner = bam::io::Reader::from(mt);
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((header, SeqReader::Bam { inner, path: path.to_path_buf() }))
        }
        Format::Cram => {
            let repo = build_repository(require_reference(path, reference)?)?;
            let mut inner = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(path)
                .map_err(|e| AnalysisError::io(path, e))?;
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((header, SeqReader::Cram { inner, path: path.to_path_buf() }))
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
                Box::new(inner.records(header).map(move |r| r.map_err(|e| AnalysisError::io(&path, e))))
            }
        }
    }
}

// ---- indexed (region) reading ---------------------------------------------

/// An indexed reader over BAM or CRAM. Hold it and call [`IdxReader::query`].
pub enum IdxReader {
    Bam { inner: bam::io::IndexedReader<bgzf::Reader<File>>, path: PathBuf },
    Cram { inner: cram::io::IndexedReader<File>, path: PathBuf },
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
            Ok((header, IdxReader::Bam { inner, path: path.to_path_buf() }))
        }
        Format::Cram => {
            let repo = build_repository(require_reference(path, reference)?)?;
            let mut inner = cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(path)
                .map_err(|e| AnalysisError::io(path, e))?;
            let header = inner.read_header().map_err(|e| AnalysisError::io(path, e))?;
            Ok((header, IdxReader::Cram { inner, path: path.to_path_buf() }))
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
                Ok(Box::new(q.map(move |r| {
                    let rec = r.map_err(|e| AnalysisError::io(&path, e))?;
                    RecordBuf::try_from_alignment_record(header, &rec).map_err(|e| AnalysisError::io(&path, e))
                })))
            }
            IdxReader::Cram { inner, path } => {
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

impl IdxReader {
    /// Drive `sink` over every record overlapping `region` (BAM: lazy record; CRAM: `RecordBuf`).
    /// A record that fails to read aborts with an error. The allocation-free counterpart to
    /// [`IdxReader::query`] (which copies each record into an owned `RecordBuf`).
    pub fn for_each<S: RecordSink>(
        &mut self,
        header: &sam::Header,
        region: &Region,
        sink: &mut S,
    ) -> Result<(), AnalysisError> {
        match self {
            IdxReader::Bam { inner, path } => {
                let path = path.clone();
                let q = inner.query(header, region).map_err(|e| AnalysisError::io(&path, e))?;
                for r in q {
                    sink.accept(&r.map_err(|e| AnalysisError::io(&path, e))?);
                }
                Ok(())
            }
            IdxReader::Cram { inner, path } => {
                let path = path.clone();
                let q = inner.query(header, region).map_err(|e| AnalysisError::io(&path, e))?;
                for r in q {
                    sink.accept(&r.map_err(|e| AnalysisError::io(&path, e))?);
                }
                Ok(())
            }
        }
    }

    /// Drive `sink` over the unplaced unmapped records (BAM only; CRAM errors, as in
    /// [`IdxReader::query_unmapped`]).
    pub fn for_each_unmapped<S: RecordSink>(&mut self, sink: &mut S) -> Result<(), AnalysisError> {
        match self {
            IdxReader::Bam { inner, path } => {
                let path = path.clone();
                let q = inner.query_unmapped().map_err(|e| AnalysisError::io(&path, e))?;
                for r in q {
                    sink.accept(&r.map_err(|e| AnalysisError::io(&path, e))?);
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

// ---- header-only ----------------------------------------------------------

/// Read just the SAM header (e.g. to resolve a contig length). `reference` is required
/// for CRAM.
pub fn read_header(path: &Path, reference: Option<&Path>) -> Result<sam::Header, AnalysisError> {
    open_seq(path, reference).map(|(header, _)| header)
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
}
