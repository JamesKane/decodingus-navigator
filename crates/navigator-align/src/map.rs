//! Map reads against a cached index — the pass that turns reverted reads back into an alignment.
//!
//! ## The part-by-part problem
//!
//! [`crate::index`] deliberately builds the index in parts so no single one has to fit in memory.
//! That buys the memory bound but creates an obligation here: a read must be mapped against
//! *every* part, and the per-part results merged, or it will be placed against whichever fraction
//! of the genome happened to be resident. Worse, MAPQ is a statement about how much better the
//! best hit is than the second best — a claim that is only meaningful genome-wide. Merging is
//! therefore not an optimization, it is what makes a split index produce the same answer as a
//! whole one.
//!
//! So:
//!
//! ```text
//! part 0 ──map all reads──> part-0 hits ─┐
//! part 1 ──map all reads──> part-1 hits ─┼──> merge per read ──> re-rank, recompute MAPQ ──> SAM
//! part 2 ──map all reads──> part-2 hits ─┘
//! ```
//!
//! Each pass holds one part; the per-part hits go to scratch rather than memory. The reads are
//! streamed once per part, which is the same trade minimap2's own `--split-prefix` makes.
//!
//! The merge itself — re-ranking across parts and recomputing MAPQ — is
//! `minimap2::index::split::merge_split_query_records`, reused rather than reimplemented. That is
//! the subtlest arithmetic in the pipeline and the place an independent implementation would most
//! likely be quietly wrong.
//!
//! ## Why not the upstream file-level entry points
//!
//! `minimap2-pure-rs` ships `map_file_sam_split` and friends, which look like exactly this. They
//! can not be used: they write to **stdout** (unusable from a desktop app) and they take
//! `parts: &[MmIdx]`, holding every part resident — giving up the entire memory bound this design
//! exists to buy. What is reused is the per-part record format and the merge; the loop is ours.
//!
//! ## Scope
//!
//! Single-end, which is what the long-read presets need. Paired-end — `sr`, and so most vendor
//! WGS — lives in [`crate::pe`], which reuses the part-by-part machinery here and adds fragment
//! mapping, pairing, and the mate-facing SAM fields.

use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use minimap2::bseq::BseqFile;
use minimap2::flags::{IdxFlags, MapFlags};
use minimap2::index::reader::IdxReader;
use minimap2::index::split;
use minimap2::index::MmIdx;
use minimap2::map::map_query;
use minimap2::options::{mapopt_update, MapOpt};

use crate::error::AlignError;
use crate::index::path_str;
use crate::output::{AlignmentWriter, OutputFormat};
use crate::preset::Preset;

/// How often the read loop asks whether it has been cancelled — same reasoning as the analysis
/// walkers: often enough that a click feels immediate, rarely enough to stay off the profile.
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// Tuning for [`map_reads`].
#[derive(Debug, Clone)]
pub struct MapParams {
    pub preset: Preset,
    /// Mapping threads. Zero means "ask the machine".
    pub threads: usize,
    /// An `@RG` line to stamp into the header and onto each record, if the source had one.
    pub read_group: Option<String>,
    /// Output container. BAM by default — see [`crate::output`] for why not CRAM here.
    pub format: OutputFormat,
    /// The reference FASTA, required only for CRAM output.
    pub reference: Option<std::path::PathBuf>,
}

impl Default for MapParams {
    fn default() -> Self {
        Self {
            preset: Preset::ShortRead,
            threads: 0,
            read_group: None,
            format: OutputFormat::default(),
            reference: None,
        }
    }
}

impl MapParams {
    fn thread_count(&self) -> usize {
        if self.threads > 0 {
            return self.threads;
        }
        std::env::var("NAVIGATOR_ALIGN_THREADS")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|n| *n > 0)
            .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4))
    }
}

/// What the mapping pass did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MapStats {
    /// Reads read from the input.
    pub queries: u64,
    /// Reads with at least one alignment.
    pub mapped: u64,
    /// Reads with none — written as unmapped SAM records, never dropped.
    pub unmapped: u64,
    /// Index parts the reference was split into. More than one means the merge path ran.
    pub parts: usize,
}

/// Cancellation, as a callback rather than a shared token type.
///
/// This crate is a leaf — it deliberately does not depend on `navigator-analysis`, so it can not
/// take that crate's `CancelToken` without inverting the layering. A closure lets the caller wire
/// whatever cancellation it already has, and costs this crate no dependency.
pub type CancelFn<'a> = &'a dyn Fn() -> bool;

/// Progress: `(reads_done, parts_done, parts_total)`. `parts_total` is only known once the index
/// has been walked, so it is zero during the first pass.
pub type ProgressFn<'a> = &'a mut dyn FnMut(u64, usize, usize);

/// Map `reads` against the index at `index_path`, writing SAM to `out`.
///
/// `scratch` holds the per-part intermediates when the index is split; it is cleaned up before
/// returning, on success or failure.
pub fn map_reads(
    index_path: &Path,
    reads: &Path,
    out: &Path,
    scratch: &Path,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    std::fs::create_dir_all(scratch).map_err(|e| AlignError::io(scratch, e))?;
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AlignError::io(parent, e))?;
    }

    let (_idx_opt, mut map_opt) = minimap2::prelude::preset(params.preset.as_str())
        .map_err(|e| AlignError::Message(format!("preset {}: {e}", params.preset.as_str())))?;

    // What `-ax <preset>` sets, and both flags are load-bearing rather than cosmetic.
    //
    // `CIGAR` is what runs base-level alignment. Without it a mapping stops at chaining, so records
    // carry coordinates but no CIGAR — and, less obviously, `map_query` skips the block that
    // assigns primary/secondary status, leaving every region with `sam_pri` unset so that *every*
    // record is emitted flagged supplementary (0x800). The split path hid this, because the merge
    // re-runs that ranking unconditionally; only the whole-index fast path was affected.
    map_opt.flag |= MapFlags::OUT_SAM | MapFlags::CIGAR;

    let mut reader = open_index(index_path, params.preset)?;

    // Read the first part, then ask whether that was all of it. Knowing this up front matters: on
    // a machine large enough to hold a whole index the split machinery is pure overhead, and the
    // reads would be streamed twice for nothing.
    let Some(first) = reader.read_next().map_err(|e| AlignError::io(index_path, e))? else {
        return Err(AlignError::Message(format!(
            "{} contains no index parts",
            index_path.display()
        )));
    };
    let single_part = reader.is_eof().map_err(|e| AlignError::io(index_path, e))?;

    if single_part {
        return map_single_part(&first, reads, out, &map_opt, params, cancel, progress);
    }
    map_split(
        first,
        &mut reader,
        index_path,
        reads,
        out,
        scratch,
        &map_opt,
        params,
        cancel,
        progress,
    )
}

/// Open the cached `.mmi`. `is_idx = true` — this is a prebuilt index, not a FASTA to sketch, so
/// the sketching parameters are read back from the file rather than supplied.
pub(crate) fn open_index(index_path: &Path, preset: Preset) -> Result<IdxReader, AlignError> {
    let (idx_opt, _) = minimap2::prelude::preset(preset.as_str())
        .map_err(|e| AlignError::Message(format!("preset {}: {e}", preset.as_str())))?;
    IdxReader::open(
        &path_str(index_path)?,
        true,
        idx_opt.w as i32,
        idx_opt.k as i32,
        idx_opt.bucket_bits,
        IdxFlags::empty(),
        idx_opt.mini_batch_size,
        idx_opt.batch_size,
    )
    .map_err(|e| AlignError::io(index_path, e))
}

// ---- the whole-index fast path --------------------------------------------

fn map_single_part(
    index: &MmIdx,
    reads: &Path,
    out: &Path,
    map_opt: &MapOpt,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    let mut opt = map_opt.clone();
    mapopt_update(&mut opt, index);

    let mut writer = open_output(out, index, params)?;
    let mut queries = BseqFile::open(&path_str(reads)?).map_err(|e| AlignError::io(reads, e))?;
    let pool = thread_pool(params)?;
    let mut stats = MapStats {
        parts: 1,
        ..Default::default()
    };

    loop {
        let batch = queries
            .read_batch(opt.mini_batch_size, true)
            .map_err(|e| AlignError::io(reads, e))?;
        if batch.is_empty() {
            break;
        }
        if cancel() {
            return Err(AlignError::Cancelled);
        }

        for (record, result) in batch.iter().zip(map_batch(&pool, index, &opt, &batch)) {
            emit(
                &mut writer,
                index,
                record,
                &result.regs,
                result.rep_len,
                &opt,
                out,
                &mut stats,
            )?;
            stats.queries += 1;
        }
        progress(stats.queries, 1, 1);
    }

    writer.finish(out)?;
    progress(stats.queries, 1, 1);
    Ok(stats)
}

/// Map a batch across the pool, **preserving input order**.
///
/// Order is not cosmetic. The split path joins each part's hits to a read by position in the file,
/// so a reordered batch would silently attach one read's hits to another — the kind of corruption
/// that produces plausible alignments at wrong loci. `par_iter().collect()` preserves order, which
/// is why results are collected rather than written as they finish.
fn map_batch(
    pool: &rayon::ThreadPool,
    index: &MmIdx,
    opt: &MapOpt,
    batch: &[minimap2::bseq::BseqRecord],
) -> Vec<minimap2::map::MapResult> {
    use rayon::prelude::*;
    pool.install(|| {
        batch
            .par_iter()
            .map(|record| map_query(index, opt, &record.name, &record.seq))
            .collect()
    })
}

/// The mapping options a preset implies, with the flags SAM output requires.
///
/// `CIGAR` is load-bearing: without it `map_query` stops at chaining, emitting records with no
/// CIGAR *and* skipping the step that assigns primary/secondary status.
pub(crate) fn prepared_map_opt(preset: Preset) -> Result<MapOpt, AlignError> {
    let (_idx, mut opt) = minimap2::prelude::preset(preset.as_str())
        .map_err(|e| AlignError::Message(format!("preset {}: {e}", preset.as_str())))?;
    opt.flag |= MapFlags::OUT_SAM | MapFlags::CIGAR;
    Ok(opt)
}

/// A path as the mapper's `&str` API wants it.
pub(crate) fn path_str_of(path: &Path) -> Result<String, AlignError> {
    crate::index::path_str(path)
}

/// Mapping is the pipeline's dominant cost and is per-read independent, so it gets a pool sized to
/// the machine (or to `NAVIGATOR_ALIGN_THREADS`).
pub(crate) fn thread_pool(params: &MapParams) -> Result<rayon::ThreadPool, AlignError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(params.thread_count())
        .build()
        .map_err(|e| AlignError::Message(format!("mapping thread pool: {e}")))
}

// ---- the split path -------------------------------------------------------

/// Per-part pass, then a merge pass. `first` has already been read from `reader`.
#[allow(clippy::too_many_arguments)]
fn map_split(
    first: MmIdx,
    reader: &mut IdxReader,
    index_path: &Path,
    reads: &Path,
    out: &Path,
    scratch: &Path,
    map_opt: &MapOpt,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    let prefix = path_str(&scratch.join("part"))?;
    // Every intermediate is removed before returning, whatever happens — these are per-read hit
    // blocks for a whole WGS and would otherwise be left behind at genome scale.
    let cleanup = ScratchGuard {
        prefix: prefix.clone(),
        parts: 0,
    };
    let mut cleanup = cleanup;

    // Header-only accumulation of every part's sequences. This is the one thing that must span all
    // parts, and it is safe to: names and lengths for a few hundred contigs, not index data.
    let mut merged_header = header_only(&first);
    let mut part_opts = vec![part_opt(map_opt, &first)];
    let mut rid_shifts = vec![0u32];

    let pool = thread_pool(params)?;
    let mut part = first;
    let mut parts = 0usize;
    loop {
        // Every part sees the same reads, so this is the same number each pass — reported for
        // progress, not accumulated.
        let queries_seen = write_part_hits(&prefix, parts, &part, &part_opts[parts], reads, &pool, cancel)?;
        parts += 1;
        cleanup.parts = parts;
        progress(queries_seen, parts, 0);

        // Drop this part before pulling the next: this is the line that keeps peak memory at one
        // part rather than the whole index.
        drop(part);

        let Some(next) = reader.read_next().map_err(|e| AlignError::io(index_path, e))? else {
            break;
        };
        rid_shifts.push(merged_header.seqs.len() as u32);
        append_header(&mut merged_header, &next);
        part_opts.push(part_opt(map_opt, &next));
        part = next;
    }

    let stats = merge_parts(
        &prefix,
        parts,
        &merged_header,
        &rid_shifts,
        reads,
        out,
        map_opt,
        params,
        cancel,
        progress,
    )?;
    Ok(stats)
}

/// One pass over the reads against one part, appending each read's hits to that part's scratch
/// file. Read order is the join key for the merge, so the file is positional: record *n* here is
/// read *n* of the input, in every part.
fn write_part_hits(
    prefix: &str,
    part_index: usize,
    part: &MmIdx,
    opt: &MapOpt,
    reads: &Path,
    pool: &rayon::ThreadPool,
    cancel: CancelFn<'_>,
) -> Result<u64, AlignError> {
    let path = split::split_tmp_path(prefix, part_index);
    let file = split::create_split_tmp(prefix, part_index, part).map_err(|e| AlignError::io(&path, e))?;
    let mut writer = BufWriter::with_capacity(1 << 20, file);

    let mut queries = BseqFile::open(&path_str(reads)?).map_err(|e| AlignError::io(reads, e))?;
    let with_cigar = opt.flag.contains(MapFlags::CIGAR);
    let mut seen = 0u64;
    loop {
        let batch = queries
            .read_batch(opt.mini_batch_size, true)
            .map_err(|e| AlignError::io(reads, e))?;
        if batch.is_empty() {
            break;
        }
        if cancel() {
            return Err(AlignError::Cancelled);
        }

        // Order-preserving, and it must be: this file is joined to the reads positionally.
        for result in map_batch(pool, part, opt, &batch) {
            let block = split::SplitQueryRecord {
                n_reg: result.regs.len() as i32,
                rep_len: result.rep_len,
                frag_gap: result.frag_gap,
                regs: result.regs,
            };
            split::write_split_query_record(&mut writer, &block, with_cigar).map_err(|e| AlignError::io(&path, e))?;
            seen += 1;
        }
    }
    writer.flush().map_err(|e| AlignError::io(&path, e))?;
    Ok(seen)
}

/// Read one hit block per part per read, merge them, and emit SAM.
#[allow(clippy::too_many_arguments)]
fn merge_parts(
    prefix: &str,
    parts: usize,
    merged_header: &MmIdx,
    rid_shifts: &[u32],
    reads: &Path,
    out: &Path,
    map_opt: &MapOpt,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    let mut opt = map_opt.clone();
    mapopt_update(&mut opt, merged_header);
    let with_cigar = opt.flag.contains(MapFlags::CIGAR);

    let mut part_readers = Vec::with_capacity(parts);
    for part in 0..parts {
        let path = split::split_tmp_path(prefix, part);
        let file = std::fs::File::open(&path).map_err(|e| AlignError::io(&path, e))?;
        let mut r = BufReader::with_capacity(1 << 20, file);
        // Step past the header this part's writer stamped, leaving the reader on record 0.
        split::read_split_header(&mut r).map_err(|e| AlignError::io(&path, e))?;
        part_readers.push((path, r));
    }

    let mut writer = open_output(out, merged_header, params)?;
    let mut queries = BseqFile::open(&path_str(reads)?).map_err(|e| AlignError::io(reads, e))?;
    let mut stats = MapStats {
        parts,
        ..Default::default()
    };

    while let Some(record) = queries
        .read_record_with_qual(true)
        .map_err(|e| AlignError::io(reads, e))?
    {
        check_cancel(stats.queries, cancel)?;

        let mut blocks = Vec::with_capacity(parts);
        for (path, r) in part_readers.iter_mut() {
            blocks.push(split::read_split_query_record(r, with_cigar).map_err(|e| AlignError::io(&*path, e))?);
        }

        let merged = split::merge_split_query_records(&blocks, rid_shifts, &opt, merged_header.k, record.l_seq as i32);
        emit(
            &mut writer,
            merged_header,
            &record,
            &merged.regs,
            merged.rep_len,
            &opt,
            out,
            &mut stats,
        )?;
        stats.queries += 1;
        if stats.queries % CANCEL_CHECK_INTERVAL == 0 {
            progress(stats.queries, parts, parts);
        }
    }

    writer.finish(out)?;
    progress(stats.queries, parts, parts);
    Ok(stats)
}

// ---- shared helpers -------------------------------------------------------

/// Write one read's SAM records: the primary (or an unmapped record) plus any supplementaries.
///
/// An unmapped read gets a record rather than silence. Realignment exists partly to find reads the
/// old reference could not place, so which reads failed *here* is information, not noise.
#[allow(clippy::too_many_arguments)]
fn emit(
    writer: &mut AlignmentWriter,
    index: &MmIdx,
    record: &minimap2::bseq::BseqRecord,
    regs: &[minimap2::types::AlignReg],
    rep_len: i32,
    opt: &MapOpt,
    out: &Path,
    stats: &mut MapStats,
) -> Result<(), AlignError> {
    if regs.is_empty() {
        let line = minimap2::format::sam::write_sam_record(
            index,
            &record.name,
            &record.seq,
            &record.qual,
            None,
            0,
            regs,
            opt.flag,
            rep_len,
        );
        writer.write_line_with(&line, out, |_, _| {})?;
        stats.unmapped += 1;
        return Ok(());
    }

    for reg in regs.iter().filter(|reg| crate::pe::emits_record(opt, reg)) {
        let line = minimap2::format::sam::write_sam_record(
            index,
            &record.name,
            &record.seq,
            &record.qual,
            Some(reg),
            regs.len(),
            regs,
            opt.flag,
            rep_len,
        );
        writer.write_line_with(&line, out, |_, _| {})?;
    }
    stats.mapped += 1;
    Ok(())
}

/// Open the output container and write the header the mapper describes.
pub(crate) fn open_output(out: &Path, index: &MmIdx, params: &MapParams) -> Result<AlignmentWriter, AlignError> {
    // `@PG` args: the design asks the realigned header to carry a program record for this step.
    let args = vec!["navigator-align".to_string(), format!("-x{}", params.preset.as_str())];
    let header = minimap2::format::sam::write_sam_hdr(index, params.read_group.as_deref(), &args);
    AlignmentWriter::create(out, params.format, &header, params.reference.as_deref())
}

pub(crate) fn part_opt(base: &MapOpt, part: &MmIdx) -> MapOpt {
    // Per-part thresholds: `mapopt_update` derives occurrence cutoffs from the index's own
    // statistics, so a part must be scored against its own, not the whole reference's.
    let mut opt = base.clone();
    mapopt_update(&mut opt, part);
    opt
}

/// A metadata-only copy of an index part: sequence names and lengths, no minimizers.
///
/// SAM needs `RNAME` and `@SQ` for every contig across every part, and after the merge a region's
/// `rid` indexes that concatenation. Carrying names and lengths for a few hundred contigs costs
/// nothing; carrying the parts themselves would undo the whole design.
pub(crate) fn header_only(part: &MmIdx) -> MmIdx {
    let mut header = MmIdx::new(part.w, part.k, part.bucket_bits, IdxFlags::empty());
    append_header(&mut header, part);
    header
}

pub(crate) fn append_header(header: &mut MmIdx, part: &MmIdx) {
    let mut offset = header.seqs.last().map(|s| s.offset + s.len as u64).unwrap_or(0);
    for seq in &part.seqs {
        let mut seq = seq.clone();
        seq.offset = offset;
        offset += seq.len as u64;
        if seq.is_alt {
            header.n_alt += 1;
        }
        header.seqs.push(seq);
    }
}

fn check_cancel(seen: u64, cancel: CancelFn<'_>) -> Result<(), AlignError> {
    if seen % CANCEL_CHECK_INTERVAL == 0 && cancel() {
        return Err(AlignError::Cancelled);
    }
    Ok(())
}

/// Removes the per-part scratch on the way out, including on an error or a cancel.
pub(crate) struct ScratchGuard {
    prefix: String,
    pub(crate) parts: usize,
}

impl ScratchGuard {
    pub(crate) fn new(prefix: String) -> Self {
        Self { prefix, parts: 0 }
    }
}

impl Drop for ScratchGuard {
    fn drop(&mut self) {
        if self.parts > 0 {
            let _ = split::remove_split_tmps(&self.prefix, self.parts);
        }
    }
}

/// The scratch path a caller should hand [`map_reads`], under a job directory.
pub fn scratch_dir(job_dir: &Path) -> PathBuf {
    job_dir.join("align-scratch")
}

#[cfg(test)]
mod tests;
