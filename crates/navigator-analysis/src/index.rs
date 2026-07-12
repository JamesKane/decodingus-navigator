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
/// a clear error rather than a corrupt index. CRAM needs no reference to index (only alignment
/// spans, which the container headers carry), so `reference` is unused today but kept in the
/// signature to stay drop-in with the reader/decode helpers.
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

/// Index a CRAM by delegating to `noodles`' container walk (which needs no reference — only the
/// alignment spans the container headers carry). It exposes no incremental offset, so progress is
/// reported as indeterminate: one `(0, None)` heartbeat at the start, then completion.
fn build_crai(path: &Path, dst: &Path, progress: ProgressFn) -> Result<(), AnalysisError> {
    progress(0, None);
    let index = cram::fs::index(path).map_err(|e| AnalysisError::io(path, e))?;
    crai::fs::write(dst, &index).map_err(|e| AnalysisError::io(dst, e))?;
    progress(1, None);
    Ok(())
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
