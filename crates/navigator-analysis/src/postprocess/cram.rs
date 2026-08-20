//! Compress a sorted alignment that carries its duplicate marks to a CRAM, and make its index.
//! This is the last step of stage C.
//!
//! A CRAM stores each read as its *difference* from the reference, and not as bases. That is why a
//! CRAM is much smaller than a BAM, and why a reader needs the reference to read it back.
//! Navigator already reads a CRAM this way, and `reader::open_seq` takes a reference for exactly
//! that. So a realigned alignment comes out in the same shape as a vendor one.
//!
//! ## The reference is part of the file
//!
//! Nothing can read a CRAM without the reference that it went out against. The file is not
//! "worse": no reader can read it at all. So the reference argument here is a matter of
//! correctness, and not a control that somebody tunes. It is also why the row of a realigned
//! alignment records `reference_path` beside `bam_path`, in stage D.
//!
//! A compression against the wrong reference makes a file that decodes to the wrong bases. It does
//! not fail. So the caller must give the reference that the mapper mapped the reads to.
//!
//! ## The order matters, and the code checks it
//!
//! The compression of a CRAM needs records in coordinate order, and near to the reference. Read
//! order instead makes a file that is slow to write, and larger than the BAM that it replaced.
//! [`write_cram`] does not trust the caller to have sorted first. It reads the `@HD SO` that the
//! sort wrote, and it refuses input that does not say coordinate order.
//!
//! ## The index
//!
//! A `.crai` turns a region query from a scan of the whole file into a seek. The code builds it
//! from a read of the finished CRAM, and it notes where each container starts. That must happen
//! after the file is complete, and not during the write.

use std::path::{Path, PathBuf};

use noodles::sam::alignment::io::Write as _;
use noodles::{cram, fasta, sam};

use super::bamio;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// Write buffer under the CRAM encoder. Matches [`bamio`]'s, for the same reason: containers arrive
/// far larger than `BufWriter`'s 8 KB default, which coalesces nothing.
const CRAM_WRITE_BUFFER: usize = 1 << 20;

/// What the CRAM step produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CramOutput {
    pub cram: PathBuf,
    pub index: PathBuf,
    /// Records written.
    pub records: u64,
    /// The count of records that are not primary, and that the code dropped because they hold no
    /// `SEQ`. See [`is_unencodable`].
    pub sequenceless_dropped: u64,
}

/// Can the code not write this record as a difference from the reference?
///
/// A CRAM stores a read as its difference against the reference. To encode one means a walk over
/// its CIGAR, and a comparison of the bases. A record with an aligned CIGAR, and a `SEQ` of `*`,
/// has no base to compare. There, noodles indexes the empty sequence, and it does not check first.
/// The result is a panic from inside the writer, ten hours into a WGS job:
/// `range end index N out of range for slice of length 0`.
///
/// This input is not malformed. SAM lets a secondary alignment carry a `SEQ` of `*`, and minimap2
/// uses that. The primary alignment alone holds the bases. A secondary one points at another locus
/// that the same read could have come from.
///
/// So to drop these loses no read, because the primary holds the sequence. But a CRAM can not
/// represent these records, and they must not reach the writer.
fn is_unencodable(record: &noodles::sam::alignment::RecordBuf) -> bool {
    !record.cigar().as_ref().is_empty() && record.sequence().as_ref().is_empty()
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

    let mut reader = bamio::open(input)?;
    let header = reader.read_header().map_err(|e| AnalysisError::io(input, e))?;
    require_coordinate_sorted(&header, input)?;

    let repository = fasta_repository(reference)?;
    // The final CRAM is tens of GB. It was the last writer in the pipeline that still gave its
    // output straight to the page cache. `build_from_path` opens the file itself, and that is how
    // an encoder comes to hold a raw `File` that nothing paces and nothing counts.
    let file = std::fs::File::create(output).map_err(|e| AnalysisError::io(output, e))?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_writer(std::io::BufWriter::with_capacity(
            CRAM_WRITE_BUFFER,
            navigator_resource::PacedFile::new(file),
        ));
    writer.write_header(&header).map_err(|e| AnalysisError::io(output, e))?;

    let mut records = 0u64;
    let mut sequenceless_dropped = 0u64;
    for (i, result) in reader.record_bufs(&header).enumerate() {
        if i as u64 % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
            progress(records);
        }
        let record = result.map_err(|e| AnalysisError::io(input, e))?;

        if is_unencodable(&record) {
            let flags = record.flags();
            // A secondary record, and a supplementary one, each give a second opinion about where
            // a read could go. The read itself lives on its primary record, so to drop one of
            // those costs no sequence.
            //
            // A *primary* record with no SEQ is a real read that goes missing. Nothing may take
            // that in without a word, after the hours that this run took to reach here.
            if !flags.is_secondary() && !flags.is_supplementary() {
                return Err(AnalysisError::Message(format!(
                    "{}: primary record {} has an aligned CIGAR but no SEQ, so it cannot be encoded \
                     as CRAM and dropping it would lose the read",
                    input.display(),
                    record
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()))
                        .unwrap_or_default(),
                )));
            }
            sequenceless_dropped += 1;
            continue;
        }

        writer
            .write_alignment_record(&header, &record)
            .map_err(|e| AnalysisError::io(output, e))?;
        records += 1;
    }

    // A CRAM collects records into containers. It writes the last container, and the end-of-file
    // marker, at shutdown alone. A writer that drops leaves a file that looks complete and is
    // not.
    writer.try_finish(&header).map_err(|e| AnalysisError::io(output, e))?;
    {
        use std::io::Write as _;
        let buffered = writer.get_mut();
        buffered.flush().map_err(|e| AnalysisError::io(output, e))?;
        // Sync it before the code makes its index and gives it to the workspace. `index_cram`
        // reads the file back at once, and every later step takes this path as the finished
        // alignment.
        buffered.get_ref().sync().map_err(|e| AnalysisError::io(output, e))?;
    }
    progress(records);

    let index = index_cram(output)?;
    Ok(CramOutput {
        cram: output.to_path_buf(),
        index,
        records,
        sequenceless_dropped,
    })
}

/// Build the `.crai` of a finished CRAM. Returns its path.
///
/// It sits apart from [`write_cram`] because it reads the completed file back. So it is also the
/// function that a repair which makes the index again would call.
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
/// The sort writes `@HD SO:coordinate`. Anything else here means that the caller left the sort
/// out. The CRAM would then be slow to write, larger than the BAM, and of no use for a region
/// query. To fail now is better than to find that out after hours of compression.
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
