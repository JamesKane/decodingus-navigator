//! Alignment output — SAM, BAM, or CRAM, through noodles.
//!
//! ## Why this sits between the mapper and the file
//!
//! `minimap2-pure-rs` formats records as SAM *text* and nothing else. Two problems follow. Every
//! downstream stage in Navigator reads BAM/CRAM through noodles, so SAM text would have to be
//! converted somewhere anyway; and the paired fields the mapper does not fill in were being
//! patched into that text by column position, which is fragile in a way that fails silently if the
//! formatter's layout ever shifts.
//!
//! So each record makes one hop through a type: the mapper's line is parsed into a
//! [`RecordBuf`], the paired fields are set on the *typed* record, and noodles writes it. The
//! parse costs something, but the mapper's own formatting already allocates a string per record,
//! and both are noise next to alignment itself — mapping a WGS is hours, serializing it is
//! minutes. What it buys is that `RNEXT` can no longer be written into the `TLEN` column.
//!
//! The alternative — building records from `AlignReg` directly and skipping SAM text — was
//! rejected on purpose. It would mean reimplementing CIGAR emission, clipping, `SEQ` orientation,
//! and every tag, which is the delicate part of the mapper's output and exactly the code most
//! worth *not* rewriting.
//!
//! ## Which format
//!
//! **BAM is the right choice for this stage.** The mapper emits reads in input order, and CRAM's
//! compression assumes coordinate-sorted, reference-adjacent records — writing it here would be
//! both slow and large. CRAM belongs after the sort in stage C. It is supported anyway because
//! the cost is one enum arm and callers past the sort will want it.

use std::io::{BufWriter, Write};
use std::path::Path;

use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, cram, fasta, sam};

use crate::error::AlignError;

/// On-disk container for the mapper's output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    /// Uncompressed SAM text. Useful for tests and eyeballing; wasteful for a real run.
    Sam,
    /// The default, and what stage C expects to sort.
    #[default]
    Bam,
    /// Reference-compressed. Needs `reference`, and really wants coordinate-sorted input, so it
    /// belongs after the sort rather than here.
    Cram,
}

impl OutputFormat {
    /// Guess from the file extension, defaulting to BAM.
    pub fn from_path(path: &Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some(e) if e.eq_ignore_ascii_case("sam") => OutputFormat::Sam,
            Some(e) if e.eq_ignore_ascii_case("cram") => OutputFormat::Cram,
            _ => OutputFormat::Bam,
        }
    }
}

/// Writes mapped records in the chosen container.
pub struct AlignmentWriter {
    header: sam::Header,
    inner: Inner,
}

enum Inner {
    Sam(sam::io::Writer<BufWriter<std::fs::File>>),
    Bam(bam::io::Writer<noodles::bgzf::io::Writer<BufWriter<std::fs::File>>>),
    Cram(Box<cram::io::Writer<std::fs::File>>),
}

impl AlignmentWriter {
    /// Open `path` and write the header.
    ///
    /// `header_text` is the mapper's `@HD`/`@SQ`/`@RG`/`@PG` block, parsed here so the typed
    /// header can be handed to the encoders — BAM stores reference names as indices into it, so
    /// this is not merely decorative.
    pub fn create(
        path: &Path,
        format: OutputFormat,
        header_text: &str,
        reference: Option<&Path>,
    ) -> Result<Self, AlignError> {
        let header = parse_header(header_text)?;

        let inner = match format {
            OutputFormat::Sam => {
                let mut w = sam::io::Writer::new(BufWriter::new(create_file(path)?));
                w.write_header(&header).map_err(|e| AlignError::io(path, e))?;
                Inner::Sam(w)
            }
            OutputFormat::Bam => {
                let mut w = bam::io::Writer::new(BufWriter::new(create_file(path)?));
                w.write_header(&header).map_err(|e| AlignError::io(path, e))?;
                Inner::Bam(w)
            }
            OutputFormat::Cram => {
                let reference = reference.ok_or_else(|| {
                    AlignError::Message("CRAM output needs the reference FASTA it will be compressed against".into())
                })?;
                let repository = fasta_repository(reference)?;
                let mut w = cram::io::writer::Builder::default()
                    .set_reference_sequence_repository(repository)
                    .build_from_path(path)
                    .map_err(|e| AlignError::io(path, e))?;
                w.write_header(&header).map_err(|e| AlignError::io(path, e))?;
                Inner::Cram(Box::new(w))
            }
        };

        Ok(Self { header, inner })
    }

    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Parse one SAM line from the mapper and hand it to `edit` before writing.
    ///
    /// `edit` is where paired fields get set. It sees a typed record, so it cannot write a value
    /// into the wrong column — which was the entire failure mode this module removes.
    pub fn write_line_with(
        &mut self,
        line: &str,
        path: &Path,
        edit: impl FnOnce(&mut RecordBuf, &sam::Header),
    ) -> Result<(), AlignError> {
        let mut record = parse_record(line, &self.header, path)?;
        edit(&mut record, &self.header);
        self.write_record(&record, path)
    }

    pub fn write_record(&mut self, record: &RecordBuf, path: &Path) -> Result<(), AlignError> {
        let header = &self.header;
        match &mut self.inner {
            Inner::Sam(w) => w.write_alignment_record(header, record),
            Inner::Bam(w) => w.write_alignment_record(header, record),
            Inner::Cram(w) => w.write_alignment_record(header, record),
        }
        .map_err(|e| AlignError::io(path, e))
    }

    /// Flush and close. CRAM in particular must be finished explicitly — its final container is
    /// only written on shutdown, so a dropped writer yields a truncated file.
    pub fn finish(self, path: &Path) -> Result<(), AlignError> {
        match self.inner {
            Inner::Sam(mut w) => w.get_mut().flush().map_err(|e| AlignError::io(path, e)),
            // BAM is BGZF, which ends with a specific empty block. Flushing alone leaves the file
            // without it, and readers treat that as truncated.
            Inner::Bam(mut w) => w.try_finish().map_err(|e| AlignError::io(path, e)),
            Inner::Cram(mut w) => w.try_finish(&self.header).map_err(|e| AlignError::io(path, e)),
        }
    }
}

fn create_file(path: &Path) -> Result<std::fs::File, AlignError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AlignError::io(parent, e))?;
    }
    std::fs::File::create(path).map_err(|e| AlignError::io(path, e))
}

fn fasta_repository(reference: &Path) -> Result<fasta::Repository, AlignError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .map_err(|e| AlignError::io(reference, e))?;
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(
        reader,
    )))
}

fn parse_header(text: &str) -> Result<sam::Header, AlignError> {
    // The mapper hands back the block without a trailing newline; the parser wants line-terminated
    // input, and an unterminated last line is silently dropped.
    let mut owned = text.to_string();
    if !owned.ends_with('\n') {
        owned.push('\n');
    }
    let mut reader = sam::io::Reader::new(std::io::Cursor::new(owned.into_bytes()));
    reader
        .read_header()
        .map_err(|e| AlignError::Message(format!("could not parse the SAM header: {e}")))
}

fn parse_record(line: &str, header: &sam::Header, path: &Path) -> Result<RecordBuf, AlignError> {
    let raw = sam::Record::try_from(line.as_bytes())
        .map_err(|e| AlignError::Message(format!("could not parse a SAM record: {e}")))?;
    RecordBuf::try_from_alignment_record(header, &raw).map_err(|e| AlignError::io(path, e))
}

/// Read every record from a SAM file. Test-facing: it lets a test assert on typed fields rather
/// than on column positions, which is the same reason the writer exists.
pub fn read_all(path: &Path) -> Result<(sam::Header, Vec<RecordBuf>), AlignError> {
    let file = std::fs::File::open(path).map_err(|e| AlignError::io(path, e))?;
    let mut reader = sam::io::Reader::new(std::io::BufReader::new(file));
    let header = reader.read_header().map_err(|e| AlignError::io(path, e))?;
    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        records.push(result.map_err(|e| AlignError::io(path, e))?);
    }
    Ok((header, records))
}

/// Read every record from a BAM file, for the same reason as [`read_all`].
pub fn read_all_bam(path: &Path) -> Result<(sam::Header, Vec<RecordBuf>), AlignError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|e| AlignError::io(path, e))?;
    let header = reader.read_header().map_err(|e| AlignError::io(path, e))?;
    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        records.push(result.map_err(|e| AlignError::io(path, e))?);
    }
    Ok((header, records))
}
