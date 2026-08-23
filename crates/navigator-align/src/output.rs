//! Alignment output: SAM, BAM, or CRAM, through noodles.
//!
//! ## Why this sits between the mapper and the file
//!
//! `minimap2-pure-rs` writes records as SAM *text* and nothing else. Two problems come from that.
//! First, every downstream stage in Navigator reads BAM or CRAM through noodles, so something must
//! convert the SAM text anyway. Second, the code patched the paired fields that the mapper does
//! not fill into that text by column position. That is fragile, and it fails with no warning if
//! the layout of the formatter ever moves.
//!
//! So each record goes through a type one time. The code parses the mapper's line into a
//! [`RecordBuf`], sets the paired fields on the *typed* record, and lets noodles write it. The
//! parse has a cost. But the mapper already allocates a string for each record when it makes the
//! SAM text. Both costs are small next to alignment itself: a WGS takes hours to map, and minutes
//! to write. The gain is that `RNEXT` can no longer go into the `TLEN` column.
//!
//! There is an alternative: build records from `AlignReg` directly, and do not make SAM text at
//! all. This module refuses that on purpose. It would need a new implementation of CIGAR emission,
//! clipping, `SEQ` orientation, and every tag. That is the delicate part of the mapper's output,
//! and it is the code this module most wants to leave alone.
//!
//! ## Which format
//!
//! **BAM is the right choice for this stage.** The mapper emits reads in input order. CRAM
//! compression expects records in coordinate order, and near to the reference. CRAM here would be
//! both slow and large. CRAM belongs after the sort, in stage C. This module still has a CRAM arm,
//! because the cost is one enum arm, and a caller after the sort will want it.

use std::io::{BufWriter, Write};
use std::path::Path;

use navigator_resource::PacedFile;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, sam};

use crate::error::AlignError;

/// Write buffer under the container encoders.
///
/// BGZF hands down ~64 KB blocks, so the 8 KB default of `BufWriter` combined nothing at all.
/// This stage writes the largest file in the pipeline, and it went to the disk in pieces the size
/// of one block. This matches the writers in post-processing.
const WRITE_BUFFER: usize = 1 << 20;

/// On-disk container for the mapper's output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    /// Uncompressed SAM text. Useful for tests and for a look by a person. Wasteful for a real run.
    Sam,
    /// The default, and what stage C expects to sort.
    #[default]
    Bam,
    /// Reference-compressed. Needs `reference`, and needs input in coordinate order, so it
    /// belongs after the sort and not here.
    Cram,
}

impl OutputFormat {
    /// Guess from the file extension. The default is BAM.
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

/// Every arm writes through a [`PacedFile`], and that is not an accident.
///
/// This stage makes the largest file in the pipeline: about 60 GB of `mapped.bam` for a 30x WGS.
/// It writes as fast as sixteen cores can compress it. With no control, all of that goes into the
/// page cache, and the write-back becomes the problem of the operating system. On macOS it became
/// the problem of everybody. One realignment made 549 GB of file-backed memory dirty. It went 1.4x
/// above the sustained write-back limit, and the watchdog of WindowServer took the login session
/// down with the job.
///
/// The pacing puts a limit on how many bytes can wait. The count is what makes the stage visible
/// to [`navigator_resource::ResourceWatch`] at all. Until now that watch reported `0 MB/s` through
/// the longest stage in the job, because nothing wrapped the one writer it has.
enum Inner {
    Sam(sam::io::Writer<BufWriter<PacedFile>>),
    Bam(bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<PacedFile>>>),
    Cram(Box<cram::io::Writer<BufWriter<PacedFile>>>),
}

impl AlignmentWriter {
    /// Open `path` and write the header.
    ///
    /// `header_text` is the mapper's `@HD`/`@SQ`/`@RG`/`@PG` block. This function parses it, so
    /// that the encoders get a typed header. BAM stores reference names as indices into that
    /// header, so the parse is not only decoration.
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
                // Threaded BGZF, because this is where the wall clock of the mapping stage went.
                // A profile of the stage gives ~60% of the serial phase to zlib deflate, and
                // `longest_match` alone is a third of that. Sixteen cores wait for the next batch
                // while this happens. Block compression runs in parallel, and the byte stream does
                // not change.
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
                // `build_from_writer`, not `build_from_path`. The second one opens the file
                // itself, and an encoder that holds its own raw `File` is exactly the writer that
                // no counter sees.
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
    /// `edit` is where the code sets the paired fields. It sees a typed record, so it can not
    /// write a value into the wrong column. That was the whole failure mode this module removes.
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

    /// Flush and close. CRAM above all needs an explicit finish. It writes its final container
    /// only on shutdown, so a writer that drops leaves a truncated file.
    ///
    /// Each arm then syncs, which matters more here than it looks. A realignment that resumes
    /// looks for the BGZF end-of-file block at the end of this file
    /// (`navigator_analysis::postprocess::bamio::is_complete_bam`). That is how it decides whether
    /// it can use the file. A marker that is still in the page cache is a promise the disk has not
    /// made. This was wrong one time, and it cost a 59 GB intermediate. A truncated file looked
    /// complete, the resume went past it, and the code deleted the real one.
    pub fn finish(self, path: &Path) -> Result<(), AlignError> {
        match self.inner {
            Inner::Sam(mut w) => sync(w.get_mut(), path),
            // BAM is BGZF, which ends with a specific empty block. A flush alone leaves the file
            // without it, and readers treat that as truncated. On the threaded writer the code
            // must also empty the workers, which is what `finish` does.
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
/// Compression is what makes the mapping stage serial, so this wants more workers than the
/// read-side default. The consumer here is the pool of the mapper, and not one thread that parses.
/// This shares `NAVIGATOR_ALIGN_THREADS` with the mapper, so one control still sets the stage.
fn bgzf_worker_count() -> std::num::NonZeroUsize {
    let n = std::env::var("NAVIGATOR_ALIGN_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|n| *n > 0)
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4));
    std::num::NonZeroUsize::new(n.clamp(1, 8)).expect("clamped above zero")
}

/// Create `path`, and its parent directories, behind a buffer and the write pacer.
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
    // The mapper hands back the block with no newline at the end. The parser needs each line to
    // end with a newline, and it drops a last line that does not, with no warning.
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
/// This is public to the crate, so that the code can make records *off* the write thread. To make
/// SAM text from an alignment, and to parse it back, is independent work for each record. At WGS
/// scale it is the largest part of the mapping stage. See the pipeline note on
/// `pe::map_pairs_single_part`.
pub(crate) fn parse_record(line: &str, header: &sam::Header, path: &Path) -> Result<RecordBuf, AlignError> {
    let raw = sam::Record::try_from(line.as_bytes())
        .map_err(|e| AlignError::Message(format!("could not parse a SAM record: {e}")))?;
    RecordBuf::try_from_alignment_record(header, &raw).map_err(|e| AlignError::io(path, e))
}

/// Read every record from a SAM file. This is for tests: it lets a test assert on typed fields,
/// and not on column positions. That is the same reason the writer exists.
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
    /// The mapping stage writes the largest file in the pipeline. For the whole of its first WGS
    /// run it wrote that file through a bare `File`. So the resource watch, which reports what the
    /// pipeline does to the machine, logged `0 MB/s` for hours while ~60 GB went to disk.
    ///
    /// The counter is process-global, so that a writer in *this* crate lands in the same total as
    /// the writer of the sort. The only way to keep that true is to assert it from here.
    #[test]
    fn the_mappers_output_reaches_the_shared_byte_counter() {
        let dir = scratch("counted");
        let path = dir.join("out.bam");

        let before = navigator_resource::bytes_written();
        let mut writer = AlignmentWriter::create(&path, OutputFormat::Bam, HEADER, None).unwrap();
        writer.write_line_with(RECORD, &path, |_, _| {}).unwrap();
        writer.finish(&path).unwrap();

        // Strictly greater, and not an exact number. Anything else in this binary shares the
        // counter, so the claim under test is only that a counter saw these bytes.
        assert!(
            navigator_resource::bytes_written() > before,
            "the mapper's BAM output was not accounted for"
        );

        let (_, records) = read_all_bam(&path).unwrap();
        assert_eq!(records.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
