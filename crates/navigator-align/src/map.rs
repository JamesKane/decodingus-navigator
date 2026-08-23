//! Map reads against a cached index. This is the pass that makes an alignment from reverted reads.
//!
//! ## The part-by-part problem
//!
//! [`crate::index`] builds the index in parts on purpose, so that no one part has to fit in
//! memory. That gives the memory limit, but it also makes an obligation here. The mapper must map
//! a read against *every* part, and the merge must join the results from each part. If it does
//! not, the read lands against whatever fraction of the genome was in memory.
//!
//! MAPQ is worse. It says how much better the best hit is than the second best. That claim
//! is true only over the whole genome. So the merge is not an optimization. It is what makes a
//! split index give the same answer as a whole one.
//!
//! So:
//!
//! ```text
//! part 0 ──map all reads──> part-0 hits ─┐
//! part 1 ──map all reads──> part-1 hits ─┼──> merge per read ──> re-rank, recompute MAPQ ──> SAM
//! part 2 ──map all reads──> part-2 hits ─┘
//! ```
//!
//! Each pass holds one part, and the hits for that part go to scratch and not to memory. The code
//! reads the reads one time for each part. That is the same trade that minimap2's own
//! `--split-prefix` makes.
//!
//! The merge itself is `minimap2::index::split::merge_split_query_records`. It ranks the hits
//! again over all parts, and calculates MAPQ again. This module uses that function, and does not
//! write its own. That code is the most delicate arithmetic in the pipeline. A separate
//! implementation would most probably be wrong there, and give no sign of it.
//!
//! ## Why not the upstream file-level entry points
//!
//! `minimap2-pure-rs` ships `map_file_sam_split` and the functions like it, which look exactly
//! like this. This module can not use them. They write to **stdout**, which a desktop app can not
//! use, and they take `parts: &[MmIdx]`, which keeps every part in memory. That gives up the whole
//! memory limit this design exists for. This module uses the record format for each part, and the
//! merge. The loop is its own.
//!
//! ## Scope
//!
//! Single-end, which is what the long-read presets need. Paired-end lives in [`crate::pe`], which
//! is `sr` and so most vendor WGS. That module uses the part-by-part code here, and adds fragment
//! mapping, pairing, and the SAM fields for the mate.

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

/// How often the read loop asks whether the user cancelled. The reason is the same as for the
/// analysis walkers: often enough that a click feels immediate, and rarely enough to stay off the
/// profile.
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// Options for [`map_reads`].
#[derive(Debug, Clone)]
pub struct MapParams {
    pub preset: Preset,
    /// Mapping threads. Zero means "ask the machine".
    pub threads: usize,
    /// An `@RG` line to stamp into the header and onto each record, if the source had one.
    pub read_group: Option<String>,
    /// Output container. BAM by default. See [`crate::output`] for why CRAM is not the default.
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
    /// Reads with none. The writer gives them an unmapped SAM record, and never drops them.
    pub unmapped: u64,
    /// Index parts the reference was split into. More than one means the merge path ran.
    pub parts: usize,
}

/// Cancellation, as a callback and not as a shared token type.
///
/// This crate is a leaf. It does not depend on `navigator-analysis`, and that is deliberate, so it
/// can not take that crate's `CancelToken`. That would invert the layers. A closure lets the
/// caller connect whatever cancellation it already has, and it costs this crate no dependency.
pub type CancelFn<'a> = &'a dyn Fn() -> bool;

/// Progress: `(reads_done, parts_done, parts_total)`. The code knows `parts_total` only after it
/// walks the index, so the value is zero during the first pass.
pub type ProgressFn<'a> = &'a mut dyn FnMut(u64, usize, usize);

/// Map `reads` against the index at `index_path`, writing SAM to `out`.
///
/// `scratch` holds the intermediates for each part when the index is split. This function removes
/// it before it returns, on success and on failure.
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

    // What `-ax <preset>` sets. Both flags do real work, and neither is cosmetic.
    //
    // `CIGAR` is what runs base-level alignment. Without it a mapping stops at chaining, so a
    // record carries coordinates but no CIGAR. Less obviously, `map_query` then skips the block
    // that assigns primary status and secondary status. Every region keeps `sam_pri` unset, so
    // *every* record comes out with the supplementary flag (0x800). The split path hid this,
    // because the merge always ranks the hits again. Only the whole-index fast path had the
    // fault.
    map_opt.flag |= MapFlags::OUT_SAM | MapFlags::CIGAR;

    let mut reader = open_index(index_path, params.preset)?;

    // Read the first part, then ask whether that was all of it. The answer matters at this
    // point. On a machine large enough to hold a whole index, the split code is only overhead.
    // The code would then read the reads two times for nothing.
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

/// Open the cached `.mmi`. `is_idx = true`, because this is an index that already exists, and not
/// a FASTA to sketch. So the sketch parameters come back from the file, and no caller gives them.
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

/// Map a batch across the pool, and **keep the input order**.
///
/// Order is not cosmetic. The split path joins the hits of each part to a read by position in the
/// file. A batch in a different order would attach the hits of one read to another read, with no
/// warning. That corruption gives plausible alignments at wrong loci. `par_iter().collect()` keeps
/// the order, and that is why this collects the results instead of a write as each one ends.
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
/// `CIGAR` does real work here. Without it `map_query` stops at chaining. It then emits records
/// with no CIGAR, *and* it skips the step that assigns primary status and secondary status.
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

/// Mapping is the largest cost in the pipeline, and the work for each read is independent. So it
/// gets a pool that fits the machine, or that `NAVIGATOR_ALIGN_THREADS` sets.
pub(crate) fn thread_pool(params: &MapParams) -> Result<rayon::ThreadPool, AlignError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(params.thread_count())
        .build()
        .map_err(|e| AlignError::Message(format!("mapping thread pool: {e}")))
}

// ---- the split path -------------------------------------------------------

/// One pass for each part, then a merge pass. The caller already read `first` from `reader`.
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
    // This removes every intermediate before it returns, whatever happens. These are hit blocks
    // for each read of a whole WGS, and they would otherwise stay on disk at genome scale.
    let cleanup = ScratchGuard {
        prefix: prefix.clone(),
        parts: 0,
    };
    let mut cleanup = cleanup;

    // A header-only collection of the sequences of every part. This is the one thing that must
    // cover all parts, and that is safe. It holds names and lengths for a few hundred contigs,
    // and no index data.
    let mut merged_header = header_only(&first);
    let mut part_opts = vec![part_opt(map_opt, &first)];
    let mut rid_shifts = vec![0u32];

    let pool = thread_pool(params)?;
    let mut part = first;
    let mut parts = 0usize;
    loop {
        // Every part sees the same reads, so this is the same number in each pass. It goes to
        // progress, and nothing adds it up.
        let queries_seen = write_part_hits(&prefix, parts, &part, &part_opts[parts], reads, &pool, cancel)?;
        parts += 1;
        cleanup.parts = parts;
        progress(queries_seen, parts, 0);

        // Drop this part before the next one arrives. This is the line that keeps peak memory at
        // one part, and not at the whole index.
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

/// One pass over the reads against one part. It adds the hits of each read to the scratch file of
/// that part. Read order is the join key for the merge, so the file is positional. Record *n*
/// here holds the hits of read *n* of the input, in every part.
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

        // This keeps the order, and it must: the merge joins this file to the reads by position.
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

/// Read one hit block for each part and each read, merge them, and emit SAM.
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
        // Step past the header that the writer of this part wrote. The reader is then on record 0.
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
/// An unmapped read gets a record, and not silence. Realignment exists in part to find reads that
/// the old reference could not place. So which reads failed *here* is information, and not noise.
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
    // Thresholds for each part: `mapopt_update` derives occurrence cutoffs from the statistics of
    // the index itself. So a part must score against its own statistics, and not against those of
    // the whole reference.
    let mut opt = base.clone();
    mapopt_update(&mut opt, part);
    opt
}

/// A metadata-only copy of an index part: sequence names and lengths, no minimizers.
///
/// SAM needs `RNAME` and `@SQ` for every contig in every part, and after the merge the `rid` of a
/// region indexes that concatenation. Names and lengths for a few hundred contigs cost nothing to
/// keep. To keep the parts themselves would undo the whole design.
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

/// Removes the scratch of each part on the way out, and also on an error or a cancel.
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
