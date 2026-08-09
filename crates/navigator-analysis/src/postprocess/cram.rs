//! Compress a sorted, marked alignment to CRAM and index it — the last step of stage C.
//!
//! CRAM stores each read as its *difference* from the reference rather than as bases, which is why
//! it is dramatically smaller than BAM and why it needs the reference to read back. Navigator
//! already reads CRAM this way (`reader::open_seq` takes a reference for exactly this), so a
//! realigned alignment lands in the same shape as a vendor one.
//!
//! ## The reference is part of the file
//!
//! A CRAM is unreadable without the reference it was written against — not "degraded", unreadable.
//! That makes the reference argument here a correctness matter rather than a tuning knob, and it is
//! why the realigned alignment's row records `reference_path` alongside `bam_path` in stage D.
//! Compressing against the wrong reference produces a file that decodes to wrong bases rather than
//! failing, so the caller must pass the reference the reads were actually mapped to.
//!
//! ## Order matters, and it is checked
//!
//! CRAM's compression assumes coordinate-sorted, reference-adjacent records; feeding it read-order
//! input produces a file that is both slow to write and larger than the BAM it replaced. Rather
//! than trust the caller to have sorted first, [`write_cram`] reads the `@HD SO` the sort stamped
//! and refuses input that does not claim coordinate order.
//!
//! ## The index
//!
//! A `.crai` is what turns a region query from a full-file scan into a seek. It is built by reading
//! the finished CRAM back and recording where each container starts, so it necessarily happens
//! after the file is complete rather than during the write.

use std::path::{Path, PathBuf};

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, cram, fasta, sam};

use crate::cancel::CancelToken;
use crate::error::AnalysisError;

const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// What the CRAM step produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CramOutput {
    pub cram: PathBuf,
    pub index: PathBuf,
    /// Records written — the same count that went in.
    pub records: u64,
}

/// Compress `input` (a coordinate-sorted BAM) to CRAM against `reference`, then write its `.crai`.
pub fn write_cram(
    input: &Path,
    output: &Path,
    reference: &Path,
    cancel: &CancelToken,
    progress: &mut dyn FnMut(u64),
) -> Result<CramOutput, AnalysisError> {
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AnalysisError::io(parent, e))?;
    }

    let file = std::fs::File::open(input).map_err(|e| AnalysisError::io(input, e))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(|e| AnalysisError::io(input, e))?;
    require_coordinate_sorted(&header, input)?;

    let repository = fasta_repository(reference)?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(output)
        .map_err(|e| AnalysisError::io(output, e))?;
    writer.write_header(&header).map_err(|e| AnalysisError::io(output, e))?;

    let mut records = 0u64;
    for (i, result) in reader.record_bufs(&header).enumerate() {
        if i as u64 % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
            progress(records);
        }
        let record = result.map_err(|e| AnalysisError::io(input, e))?;
        writer
            .write_alignment_record(&header, &record)
            .map_err(|e| AnalysisError::io(output, e))?;
        records += 1;
    }

    // CRAM buffers records into containers and only writes the last one — and the end-of-file
    // marker — on shutdown. A dropped writer leaves a file that looks complete and is not.
    writer.try_finish(&header).map_err(|e| AnalysisError::io(output, e))?;
    progress(records);

    let index = index_cram(output)?;
    Ok(CramOutput {
        cram: output.to_path_buf(),
        index,
        records,
    })
}

/// Build the `.crai` for a finished CRAM, returning its path.
///
/// Separate from [`write_cram`] because indexing reads the completed file back, so it is also what
/// a "reindex this alignment" repair would call.
pub fn index_cram(cram_path: &Path) -> Result<PathBuf, AnalysisError> {
    let index = cram::fs::index(cram_path).map_err(|e| AnalysisError::io(cram_path, e))?;
    let index_path = crai_path(cram_path);
    cram::crai::fs::write(&index_path, &index).map_err(|e| AnalysisError::io(&index_path, e))?;
    Ok(index_path)
}

/// `foo.cram` → `foo.cram.crai`, the name every reader looks for.
pub fn crai_path(cram_path: &Path) -> PathBuf {
    let mut name = cram_path.as_os_str().to_os_string();
    name.push(".crai");
    PathBuf::from(name)
}

/// Refuse input that does not declare coordinate order.
///
/// The sort stamps `@HD SO:coordinate`; anything else here means the caller skipped the sort, and
/// the resulting CRAM would be slow to write, larger than the BAM, and unusable for region queries.
/// Failing now beats discovering that after hours of compression.
fn require_coordinate_sorted(header: &sam::Header, input: &Path) -> Result<(), AnalysisError> {
    use noodles::sam::header::record::value::map::header::tag;

    let sort_order = header
        .header()
        .and_then(|hd| hd.other_fields().get(&tag::SORT_ORDER))
        .map(|v| v.to_vec());

    match sort_order.as_deref() {
        Some(b"coordinate") => Ok(()),
        other => Err(AnalysisError::Message(format!(
            "{} is not coordinate-sorted (@HD SO is {}); sort it before writing CRAM",
            input.display(),
            other
                .map(|v| String::from_utf8_lossy(v).to_string())
                .unwrap_or_else(|| "absent".to_string()),
        ))),
    }
}

fn fasta_repository(reference: &Path) -> Result<fasta::Repository, AnalysisError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .map_err(|e| AnalysisError::io(reference, e))?;
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(
        reader,
    )))
}
