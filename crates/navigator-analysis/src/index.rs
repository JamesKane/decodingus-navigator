//! Building the coordinate index a BAM/CRAM needs for **region queries**.
//!
//! The per-contig walker ([`crate::unified`]), the callable-interval scan
//! ([`crate::coverage::callable_intervals`]), and the de-novo / STR callers all seek to a
//! `contig:start-end` via `noodles`' indexed readers, which autoload a sibling `.bai` (BAM) or
//! `.crai` (CRAM). Aligned files usually ship with one, but a re-exported or single-imported file
//! often does not — without it those paths either error outright or fall back to a whole-file linear
//! scan. This module builds the missing index once, up front, so every query path is fast.
//!
//! Index construction is a single sequential pass over the file (the same cost as one analysis
//! read), so it is worth surfacing to the user with progress. The BAM path reports a true byte
//! fraction (the compressed offset of the bgzf stream); the CRAM path delegates to `noodles`'
//! container walk, which exposes no offset hook, so it reports *indeterminate* progress
//! (`total = None`) — the UI shows a spinner for it.

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

/// Progress sink for index construction: `(done_bytes, total_bytes)`. `total` is `Some` for BAM
/// (compressed file length) and `None` for CRAM (no offset hook — indeterminate).
pub type ProgressFn<'a> = &'a mut dyn FnMut(u64, Option<u64>);

/// The sibling index path this module writes for `path`: `foo.bam` → `foo.bam.bai`,
/// `foo.cram` → `foo.cram.crai`. (This is the `.bam.bai` / `.cram.crai` spelling; the query readers
/// also accept the `.bai` / `.crai` spelling, but we always write the dotted form `samtools` does.)
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

/// Build the coordinate index for `path` **if one is not already present**, returning the path of
/// the index that was written (`Ok(None)` when a `.bai`/`.crai` already existed — nothing to do).
///
/// The BAM input must be coordinate-sorted (`SO:coordinate` in the header); an unsorted file yields
/// a clear error rather than a corrupt index. `reference` is unused: the BAM path never needs it,
/// and the CRAM path *would* need it for multi-reference slices but has no way to supply it — see
/// [`build_crai`]. It is kept in the signature to stay drop-in with the reader/decode helpers, and
/// because threading it through is what a fixed CRAM indexer would want.
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

/// Index a coordinate-sorted BAM, reporting a byte fraction from the bgzf compressed offset. This
/// mirrors `noodles`' `bam::fs::index`, but drives the record loop ourselves so we can emit
/// progress against the on-disk (compressed) file length.
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

        // Report on ~32 MB of compressed progress so a multi-GB BAM doesn't flood the channel.
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

/// Index a CRAM by delegating to `noodles`' container walk. It exposes no incremental offset, so
/// progress is reported as indeterminate: one `(0, None)` heartbeat at the start, then completion.
///
/// **This does not work on every CRAM.** A *single*-reference slice is cheap to index — its span
/// comes straight from the slice header — but a *multi*-reference slice has no one span, so noodles
/// decodes its records to derive one. Reconstructing a mapped record's sequence needs the reference
/// bases, and `cram::fs::index` hands the decoder an empty `fasta::Repository` (its own `// TODO`,
/// still open as of noodles-cram 0.95), so it panics there instead of erroring. Aligners write
/// their unmapped/decoy tail as multi-reference slices, so most real whole-genome CRAMs hit this —
/// and only at the very end of the file, after the walk has already done nearly all the work.
/// [`multi_reference_panic`] turns that panic into an actionable error rather than a crash.
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

/// The panic text noodles emits when a multi-reference slice's record decode asks the (empty)
/// repository for reference bases. Distinct from the single-reference slice's "invalid **slice**
/// reference sequence name", which would mean something genuinely different — a reference that
/// really is missing a contig — so match on the record-level wording only.
const MULTI_REFERENCE_PANIC: &str = "invalid reference sequence name";

/// Explain a panic escaping the CRAM index walk. The multi-reference case is known and has a
/// concrete workaround, so name it and give the command; anything else reports its own text
/// instead of a guess.
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
fn alignment_context(
    record: &bam::Record,
) -> std::io::Result<(Option<usize>, Option<Position>, Option<Position>)> {
    Ok((
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two noodles panics differ by one word, and they mean opposite things: the record-level
    /// one is the indexer's own limitation, the slice-level one means the reference really is
    /// missing a contig. Telling a user to run `samtools index` for the latter would be wrong.
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
        assert!(known.contains("samtools index /data/sample.cram"), "gives the command: {known}");

        // An unclassified panic reports its own text rather than borrowing the known diagnosis.
        let other = index_panic_error(path, &String::from("not yet implemented")).to_string();
        assert!(other.contains("not yet implemented"), "quotes the panic: {other}");
        assert!(!other.contains("samtools"), "no bogus workaround: {other}");
    }
}
