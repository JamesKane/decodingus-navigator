//! Build the coordinate index that a BAM or CRAM needs for a **region query**.
//!
//! Three paths seek to a `contig:start-end`: the walker over the contigs ([`crate::unified`]), the
//! scan for callable intervals ([`crate::coverage::callable_intervals`]), and the de-novo and STR
//! callers. They all go through the indexed readers of `noodles`, which load a `.bai` beside a
//! BAM, or a `.crai` beside a CRAM.
//!
//! An aligned file usually comes with one. But a file that somebody exported again, or imported on
//! its own, often does not. Without an index, those paths either fail outright, or fall back to a
//! linear scan of the whole file. This module builds the missing index once, at the start, so that
//! every query path is fast.
//!
//! To build an index is one sequential pass over the file, which costs the same as one analysis
//! read. So it is worth a progress report to the user. The BAM path gives a true fraction of the
//! bytes, from the compressed offset of the bgzf stream. The CRAM path hands the work to the
//! container walk of `noodles`, which offers no hook for an offset. So that path reports progress
//! with no end, at `total = None`, and the UI shows a spinner for it.

use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::bam::{self, bai};
use noodles::core::Position;
use noodles::cram::{self, crai};
use noodles::csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles::sam::{
    self,
    alignment::Record as _,
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};

use crate::error::AnalysisError;
use crate::reader::{detect_format, has_region_index, Format};

/// The sink for the progress of an index build, as `(done_bytes, total_bytes)`. `total` is `Some`
/// for a BAM, where it holds the compressed length of the file. It is `None` for a CRAM, which
/// offers no hook for an offset, so the progress there has no end.
pub type ProgressFn<'a> = &'a mut dyn FnMut(u64, Option<u64>);

/// The index path that this module writes beside `path`: `foo.bam` gives `foo.bam.bai`, and
/// `foo.cram` gives `foo.cram.crai`. That is the `.bam.bai` and `.cram.crai` form. The query
/// readers also accept the `.bai` and `.crai` form, but this module always writes the dotted form
/// that `samtools` writes.
pub fn index_path_for(path: &Path) -> PathBuf {
    let ext = match detect_format(path) {
        Format::Bam => "bai",
        Format::Cram => "crai",
    };
    let mut file_name = path.file_name().map(|n| n.to_os_string()).unwrap_or_default();
    file_name.push(".");
    file_name.push(ext);
    path.with_file_name(file_name)
}

/// Build the coordinate index of `path`, **if one is not there already**. Returns the path of the
/// index that it wrote. It returns `Ok(None)` when a `.bai` or `.crai` already existed, and there
/// was nothing to do.
///
/// A BAM input must be in coordinate order, with `SO:coordinate` in its header. A file that is not
/// sorted gives a clear error, and not a corrupt index.
///
/// Nothing uses `reference`. The BAM path never needs it. The CRAM path *would* need it for a
/// slice that holds more than one reference, and it has no way to get it there. See
/// [`build_crai`].
///
/// The argument stays in the signature for two reasons. It keeps this function interchangeable
/// with the reader and decode helpers, and a CRAM indexer that somebody fixes would want it.
pub fn ensure_index(
    path: &Path,
    _reference: Option<&Path>,
    progress: ProgressFn,
) -> Result<Option<PathBuf>, AnalysisError> {
    if has_region_index(path) {
        return Ok(None);
    }
    let dst = index_path_for(path);
    match detect_format(path) {
        Format::Bam => build_bai(path, &dst, progress)?,
        Format::Cram => build_crai(path, &dst, progress)?,
    }
    Ok(Some(dst))
}

/// Index a BAM that is in coordinate order. It reports a fraction of the bytes, from the
/// compressed offset of the bgzf stream.
///
/// It has the same shape as `bam::fs::index` in `noodles`. But it drives the record loop itself,
/// so that it can emit progress against the compressed length of the file on disk.
fn build_bai(path: &Path, dst: &Path, progress: ProgressFn) -> Result<(), AnalysisError> {
    let total = File::open(path)
        .and_then(|f| f.metadata())
        .map(|m| m.len())
        .map_err(|e| AnalysisError::io(path, e))?;

    // `bam::io::Reader::new` takes the raw `File` and wraps it in a bgzf reader internally, so
    // `get_ref()` below yields the bgzf reader whose `virtual_position` drives progress.
    let file = File::open(path).map_err(|e| AnalysisError::io(path, e))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(|e| AnalysisError::io(path, e))?;

    if !is_coordinate_sorted(&header) {
        return Err(AnalysisError::Message(format!(
            "cannot index {}: the BAM is not coordinate-sorted (need SO:coordinate). Sort it first, \
             or import a coordinate-sorted alignment.",
            path.display()
        )));
    }

    let mut record = bam::Record::default();
    let mut builder = Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();
    let mut last_reported = 0u64;

    loop {
        let n = reader
            .read_record(&mut record)
            .map_err(|e| AnalysisError::io(path, e))?;
        if n == 0 {
            break;
        }
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match alignment_context(&record).map_err(|e| AnalysisError::io(path, e))? {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };
        builder
            .add_record(alignment_context, chunk)
            .map_err(|e| AnalysisError::io(path, e))?;

        // Report on ~32 MB of compressed progress so a multi-GB BAM does not flood the channel.
        let done = end_position.compressed();
        if done.saturating_sub(last_reported) >= 32_000_000 {
            last_reported = done;
            progress(done, Some(total));
        }
        start_position = end_position;
    }

    let index: bai::Index = builder.build(header.reference_sequences().len());
    bai::fs::write(dst, &index).map_err(|e| AnalysisError::io(dst, e))?;
    progress(total, Some(total));
    Ok(())
}

/// Index a CRAM. It hands the work to the container walk of `noodles`. That walk offers no offset
/// as it goes, so the progress has no end: one `(0, None)` beat at the start, and then the
/// finish.
///
/// **This does not work on every CRAM.** A slice with *one* reference costs little to index,
/// because its span comes straight from the slice header. A slice with *more than one* reference
/// has no single span, so noodles decodes its records to get one.
///
/// To build the sequence of a mapped record again needs the reference bases. `cram::fs::index`
/// gives the decoder an empty `fasta::Repository`, which carries its own `// TODO`, still open at
/// noodles-cram 0.95. So it panics there, and it does not return an error.
///
/// An aligner writes its unmapped and decoy tail as slices with more than one reference. So most
/// real whole-genome CRAMs reach this. They reach it at the very end of the file, after the walk
/// has done almost all of the work. [`multi_reference_panic`] turns that panic into an error that
/// a user can act on, and not a crash.
fn build_crai(path: &Path, dst: &Path, progress: ProgressFn) -> Result<(), AnalysisError> {
    progress(0, None);
    let index = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| cram::fs::index(path))) {
        Ok(result) => result.map_err(|e| AnalysisError::io(path, e))?,
        Err(payload) => return Err(index_panic_error(path, &*payload)),
    };
    crai::fs::write(dst, &index).map_err(|e| AnalysisError::io(dst, e))?;
    progress(1, None);
    Ok(())
}

/// The panic text that noodles gives when the record decode of a slice with more than one
/// reference asks the empty repository for reference bases.
///
/// It is different from "invalid **slice** reference sequence name", which comes from a slice with
/// one reference. That message means something else: a reference that truly does not hold a
/// contig. So match on the text at the record level alone.
const MULTI_REFERENCE_PANIC: &str = "invalid reference sequence name";

/// Explain a panic that comes out of the CRAM index walk.
///
/// Somebody has already found the case of a slice with more than one reference, and there is a
/// concrete way around it. So this function names that case and gives the command. Any other
/// panic reports its own text, and this code does not guess.
fn index_panic_error(path: &Path, payload: &(dyn std::any::Any + Send)) -> AnalysisError {
    let text = crate::error::panic_text(payload).unwrap_or("no further detail");
    if multi_reference_panic(text) {
        AnalysisError::Message(format!(
            "cannot index {p}: this CRAM has multi-reference slices, which the built-in indexer \
             cannot span. Build the index with `samtools index {p}` and re-import — an existing \
             .crai is used as-is.",
            p = path.display()
        ))
    } else {
        AnalysisError::Message(format!(
            "cannot index {}: the CRAM reader hit a case it does not handle ({text})",
            path.display()
        ))
    }
}

fn multi_reference_panic(text: &str) -> bool {
    text.contains(MULTI_REFERENCE_PANIC) && !text.contains("slice reference sequence name")
}

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .map(|sort_order| sort_order == COORDINATE)
        .unwrap_or_default()
}

#[allow(clippy::type_complexity)]
fn alignment_context(record: &bam::Record) -> std::io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two panics of noodles differ by one word, and they mean opposite things. The one at
    /// the record level is a limit of the indexer itself. The one at the slice level means that
    /// the reference truly does not hold a contig. To tell a user to run `samtools index` for that
    /// second one would be wrong.
    #[test]
    fn distinguishes_the_indexer_limitation_from_a_genuinely_missing_contig() {
        assert!(multi_reference_panic("invalid reference sequence name"));
        assert!(!multi_reference_panic("invalid slice reference sequence name"));
        assert!(!multi_reference_panic("not yet implemented"));
    }

    #[test]
    fn index_panic_error_names_the_cause_and_the_workaround() {
        let path = Path::new("/data/sample.cram");

        let known = index_panic_error(path, &"invalid reference sequence name");
        let known = known.to_string();
        assert!(known.contains("multi-reference slices"), "names the cause: {known}");
        assert!(
            known.contains("samtools index /data/sample.cram"),
            "gives the command: {known}"
        );

        // A panic that this code does not know reports its own text. It does not take the
        // diagnosis of the known case.
        let other = index_panic_error(path, &String::from("not yet implemented")).to_string();
        assert!(other.contains("not yet implemented"), "quotes the panic: {other}");
        assert!(!other.contains("samtools"), "no bogus workaround: {other}");
    }
}
