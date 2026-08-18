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

use navigator_resource::PacedFile;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, sam};

use crate::error::AlignError;

/// Write buffer under the container encoders.
///
/// BGZF hands down ~64 KB blocks, so `BufWriter`'s 8 KB default coalesced nothing at all: this
/// stage's output is the largest file the pipeline produces and it was reaching the disk in
/// block-sized dribs. Matches the post-processing writers.
const WRITE_BUFFER: usize = 1 << 20;

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

/// Every arm writes through a [`PacedFile`], and that is not incidental.
///
/// This stage produces the pipeline's largest file — ~60 GB of `mapped.bam` for a 30x WGS — as fast
/// as sixteen cores can compress it, and left to itself that goes into the page cache and becomes
/// the operating system's problem to write back. On macOS it became everyone's problem: a
/// realignment dirtied 549 GB of file-backed memory, exceeded the sustained write-back limit by
/// 1.4x, and WindowServer's watchdog took the login session down with the job. Pacing caps what can
/// be outstanding; the accounting is what makes the stage visible to
/// [`navigator_resource::ResourceWatch`] at all, which until now reported `0 MB/s` through the
/// longest stage in the job because the only writer it has was unwrapped.
enum Inner {
    Sam(sam::io::Writer<BufWriter<PacedFile>>),
    Bam(bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<PacedFile>>>),
    Cram(Box<cram::io::Writer<BufWriter<PacedFile>>>),
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
                let mut w = sam::io::Writer::new(paced(path)?);
                w.write_header(&header).map_err(|e| AlignError::io(path, e))?;
                Inner::Sam(w)
            }
            OutputFormat::Bam => {
                // Threaded BGZF, because this is where the mapping stage's wall clock went. A
                // profile of the stage attributes ~60% of the serial phase to zlib deflate —
                // `longest_match` alone is a third of it — while sixteen cores wait for the next
                // batch. Block compression parallelizes; the byte stream is unchanged.
                let inner = bgzf::io::MultithreadedWriter::with_worker_count(bgzf_worker_count(), paced(path)?);
                let mut w = bam::io::Writer::from(inner);
                w.write_header(&header).map_err(|e| AlignError::io(path, e))?;
                Inner::Bam(w)
            }
            OutputFormat::Cram => {
                let reference = reference.ok_or_else(|| {
                    AlignError::Message("CRAM output needs the reference FASTA it will be compressed against".into())
                })?;
                let repository = fasta_repository(reference)?;
                // `build_from_writer`, not `build_from_path`: the latter opens the file itself, and
                // an encoder holding its own raw `File` is exactly the writer that goes uncounted.
                let mut w = cram::io::writer::Builder::default()
                    .set_reference_sequence_repository(repository)
                    .build_from_writer(paced(path)?);
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
    /// `edit` is where paired fields get set. It sees a typed record, so it can not write a value
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
    ///
    /// Each arm then syncs, which matters more here than it looks. A resumed realignment decides
    /// whether it can pick this file up by checking for the BGZF end-of-file block on the end of it
    /// (`navigator_analysis::postprocess::bamio::is_complete_bam`), and a marker still sitting in
    /// the page cache is a promise the disk has not made. Getting that wrong once already cost a
    /// 59 GB intermediate: a truncated file that looked complete was resumed past and the real one
    /// deleted.
    pub fn finish(self, path: &Path) -> Result<(), AlignError> {
        match self.inner {
            Inner::Sam(mut w) => sync(w.get_mut(), path),
            // BAM is BGZF, which ends with a specific empty block. Flushing alone leaves the file
            // without it, and readers treat that as truncated. On the threaded writer that means
            // draining the workers, which is what `finish` does.
            Inner::Bam(mut w) => {
                let mut buffered = w.get_mut().finish().map_err(|e| AlignError::io(path, e))?;
                sync(&mut buffered, path)
            }
            Inner::Cram(mut w) => {
                w.try_finish(&self.header).map_err(|e| AlignError::io(path, e))?;
                sync(w.get_mut(), path)
            }
        }
    }
}

/// Flush the buffer and push the file itself to disk.
fn sync(buffered: &mut BufWriter<PacedFile>, path: &Path) -> Result<(), AlignError> {
    buffered.flush().map_err(|e| AlignError::io(path, e))?;
    buffered.get_ref().sync().map_err(|e| AlignError::io(path, e))
}

/// Worker threads for BGZF block compression.
///
/// Compression is the mapping stage's serial bottleneck, so this wants more workers than the
/// read-side default: the consumer here is the mapper's pool, not a single parsing thread. Shares
/// `NAVIGATOR_ALIGN_THREADS` with the mapper so one knob still controls the stage.
fn bgzf_worker_count() -> std::num::NonZeroUsize {
    let n = std::env::var("NAVIGATOR_ALIGN_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|n| *n > 0)
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4));
    std::num::NonZeroUsize::new(n.clamp(1, 8)).expect("clamped above zero")
}

/// Create `path` — parents included — behind a buffer and the write pacer.
fn paced(path: &Path) -> Result<BufWriter<PacedFile>, AlignError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AlignError::io(parent, e))?;
    }
    let file = std::fs::File::create(path).map_err(|e| AlignError::io(path, e))?;
    Ok(BufWriter::with_capacity(WRITE_BUFFER, PacedFile::new(file)))
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

/// Parse one of the mapper's SAM lines into a typed record.
///
/// Public to the crate so records can be built *off* the writing thread. Formatting an alignment
/// to SAM text and parsing it back is per-record independent work, and at WGS scale it dominates
/// the mapping stage — see the pipeline note on `pe::map_pairs_single_part`.
pub(crate) fn parse_record(line: &str, header: &sam::Header, path: &Path) -> Result<RecordBuf, AlignError> {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("dun-output-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n";
    const RECORD: &str = "read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII";

    /// The regression this crate's dependency on `navigator-resource` exists for.
    ///
    /// The mapping stage writes the largest file in the pipeline, and for the whole of its first
    /// WGS run it wrote that file through a bare `File` — so the resource watch, which reports what
    /// the pipeline is doing to the machine, logged `0 MB/s` for hours while ~60 GB went to disk.
    /// The counter is process-global precisely so that a writer in *this* crate lands in the same
    /// total as the sort's, and the only way to keep that true is to assert it from here.
    #[test]
    fn the_mappers_output_reaches_the_shared_byte_counter() {
        let dir = scratch("counted");
        let path = dir.join("out.bam");

        let before = navigator_resource::bytes_written();
        let mut writer = AlignmentWriter::create(&path, OutputFormat::Bam, HEADER, None).unwrap();
        writer.write_line_with(RECORD, &path, |_, _| {}).unwrap();
        writer.finish(&path).unwrap();

        // Strictly greater, not an exact figure: the counter is shared with anything else running
        // in this binary, so the claim under test is that these bytes were counted at all.
        assert!(
            navigator_resource::bytes_written() > before,
            "the mapper's BAM output was not accounted for"
        );

        let (_, records) = read_all_bam(&path).unwrap();
        assert_eq!(records.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
