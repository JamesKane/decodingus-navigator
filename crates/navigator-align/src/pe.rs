//! Paired-end mapping: the `sr` path, and so most vendor WGS.
//!
//! A pair is not two independent reads. Mapping them together lets a mate with a confident
//! position rescue a mate that is ambiguous. The expected span of the fragment is also evidence
//! about where the second end belongs. Both of those feed MAPQ.
//!
//! So the mapper maps the two ends as one *fragment* ([`minimap2::map::map_frag_queries`]), and
//! then pairs them ([`minimap2::pe::pair`]). That step sets `proper_frag`, adjusts MAPQ, and
//! decides which region of each end is the primary.
//!
//! ## What this module has to build itself
//!
//! The pieces above are public in `minimap2-pure-rs`. **The code that makes its PE SAM text is
//! not.** That code lives in private `pipeline.rs` helpers. So this module builds the mate half of
//! each record itself: the paired flags, `RNEXT`/`PNEXT`, and `TLEN`.
//!
//! The method is the same as upstream. Make the single-end line with the public writer, which
//! already knows CIGAR, clipping, and tags. Then fill in the paired fields. That is string surgery
//! on a SAM line, and this module says so instead of a quiet name for it. But the alternative is a
//! new implementation of CIGAR emission and tag emission, which is more of the delicate work, and
//! not less. [`set_pair_fields`] is small on purpose, and it has many tests for that reason.
//!
//! ## Split indexes
//!
//! Same shape as the single-end path. Map the fragment against each part, and spill the hit
//! blocks for each segment. Then merge each end over all parts, and pair them again. The pair step
//! after the merge is necessary. Pairing that the code decides for each part would use only a
//! fraction of the genome. That is exactly the error the merge exists to stop.

use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use minimap2::bseq::{BseqFile, BseqRecord};
use minimap2::flags::MapFlags;
use minimap2::index::split;
use minimap2::index::MmIdx;
use minimap2::map::MapResult;
use minimap2::options::{mapopt_update, MapOpt};
use minimap2::types::AlignReg;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::RecordBuf;

use crate::error::AlignError;
use crate::map::{
    append_header, header_only, open_index, open_output, part_opt, path_str_of, prepared_map_opt, thread_pool,
    CancelFn, MapParams, MapStats, ProgressFn, ScratchGuard,
};
use crate::output::AlignmentWriter;

/// SAM flag bits this module sets. They have names because a bare `0x20` in flag arithmetic is
/// hard to read. The difference between `0x20` and `0x10` is a wrong strand with no warning.
mod flag {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const MATE_REVERSE: u16 = 0x20;
    pub const FIRST: u16 = 0x40;
    pub const LAST: u16 = 0x80;
}

/// Map `reads1`/`reads2` as pairs against the index, writing SAM to `out`.
///
/// The two files must be in lockstep: record *n* of each is one template. The revert stage of
/// `navigator-analysis` guarantees that for the FASTQ it makes.
#[allow(clippy::too_many_arguments)]
pub fn map_pairs(
    index_path: &Path,
    reads1: &Path,
    reads2: &Path,
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

    let map_opt = prepared_map_opt(params.preset)?;
    let mut reader = open_index(index_path, params.preset)?;

    let Some(first) = reader.read_next().map_err(|e| AlignError::io(index_path, e))? else {
        return Err(AlignError::Message(format!(
            "{} contains no index parts",
            index_path.display()
        )));
    };
    if reader.is_eof().map_err(|e| AlignError::io(index_path, e))? {
        return map_pairs_single_part(&first, reads1, reads2, out, &map_opt, params, cancel, progress);
    }

    let mut parts = vec![first];
    while let Some(next) = reader.read_next().map_err(|e| AlignError::io(index_path, e))? {
        parts.push(next);
    }
    map_pairs_split(&parts, reads1, reads2, out, scratch, &map_opt, params, cancel, progress)
}

// ---- whole-index path -----------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn map_pairs_single_part(
    index: &MmIdx,
    reads1: &Path,
    reads2: &Path,
    out: &Path,
    map_opt: &MapOpt,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    let mut opt = map_opt.clone();
    mapopt_update(&mut opt, index);

    let writer = open_output(out, index, params)?;
    let pool = thread_pool(params)?;
    let mut stats = MapStats {
        parts: 1,
        ..Default::default()
    };
    let header = writer.header().clone();
    let chunk = opt.mini_batch_size;

    // Read, map and write run as three stages that overlap, and not one after the other.
    //
    // One stage after the other left the pool idle through two phases of every cycle. A profile of
    // the stage put ~60% of the serial phase in BGZF compression, and ~15% in BAM record
    // encoding. Another ~9% was gzip inflate. All of it ran while sixteen cores waited.
    //
    // Threads on the compression removed the largest piece. Stages that overlap remove the *idle
    // time*. The reader can inflate batch n+1, and the writer can encode batch n-1, while the pool
    // maps batch n.
    //
    // The design keeps the order: one reader, one writer, and batches that cross each channel in
    // sequence. The depth is small on purpose. A batch is tens of MB. A deep queue would only
    // let a fast mapper build a backlog in memory ahead of a slow disk.
    let (read_tx, read_rx) = std::sync::mpsc::sync_channel::<Vec<(BseqRecord, BseqRecord)>>(PIPELINE_DEPTH);
    let (write_tx, write_rx) = std::sync::mpsc::sync_channel::<Vec<Vec<RecordBuf>>>(PIPELINE_DEPTH);

    let outcome = std::thread::scope(|scope| -> Result<MapStats, AlignError> {
        let reader = scope.spawn(move || -> Result<(), AlignError> {
            let mut pairs = PairReader::open(reads1, reads2)?;
            while let Some(batch) = pairs.next_batch(chunk)? {
                // A closed channel means the mapper stopped: the user cancelled it, or a stage
                // failed. The error of the mapper is the one to report, so this exits with no
                // message.
                if read_tx.send(batch).is_err() {
                    break;
                }
            }
            Ok(())
        });

        let scribe = scope.spawn(move || -> Result<(), AlignError> {
            let mut writer = writer;
            for batch in write_rx {
                for records in batch {
                    for record in &records {
                        writer.write_record(record, out)?;
                    }
                }
            }
            // The finish here, on the thread that owns the writer, is what writes the BGZF
            // end-of-file block. A writer that drops in mid-stream leaves a file that readers call
            // truncated.
            writer.finish(out)
        });

        let mapping = || -> Result<(), AlignError> {
            for batch in read_rx {
                if cancel() {
                    return Err(AlignError::Cancelled);
                }
                let mapped = map_frag_batch(&pool, index, &opt, &batch);
                let built = build_batch_records(&pool, index, &opt, &header, &batch, mapped, out)?;

                let mut records = Vec::with_capacity(built.len());
                for (batch_records, batch_stats) in built {
                    stats.queries += batch_stats.queries;
                    stats.mapped += batch_stats.mapped;
                    stats.unmapped += batch_stats.unmapped;
                    records.push(batch_records);
                }
                if write_tx.send(records).is_err() {
                    // The writer died; join below surfaces why.
                    break;
                }
                progress(stats.queries, 1, 1);
            }
            Ok(())
        };
        let mapped_result = mapping();

        // Dropped so the writer's loop can end; it is what triggers `finish`.
        drop(write_tx);

        let scribe_result = scribe
            .join()
            .map_err(|_| AlignError::Message("writer thread panicked".into()))?;
        let reader_result = reader
            .join()
            .map_err(|_| AlignError::Message("reader thread panicked".into()))?;

        // A downstream failure usually shows first as a send error upstream, so this reports the
        // real cause before the symptom.
        scribe_result?;
        mapped_result?;
        reader_result?;
        Ok(stats)
    })?;

    progress(outcome.queries, 1, 1);
    Ok(outcome)
}

/// Batches in progress in each pipeline stage.
///
/// One in hand and one in the queue is enough to keep a stage off its neighbour. More than that
/// only costs memory. A batch is `mini_batch_size` bases of reads, plus the records the code makes
/// from them.
const PIPELINE_DEPTH: usize = 2;

// ---- split path -----------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn map_pairs_split(
    parts: &[MmIdx],
    reads1: &Path,
    reads2: &Path,
    out: &Path,
    scratch: &Path,
    map_opt: &MapOpt,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    let prefix = path_str_of(&scratch.join("pe-part"))?;
    let mut cleanup = ScratchGuard::new(prefix.clone());

    // Header-only view across every part, for RNAME/@SQ and for the merge's rid arithmetic.
    let mut merged_header = header_only(&parts[0]);
    let mut rid_shifts = vec![0u32];
    for part in &parts[1..] {
        rid_shifts.push(merged_header.seqs.len() as u32);
        append_header(&mut merged_header, part);
    }

    let pool = thread_pool(params)?;
    // One scratch file for each part. Each holds two blocks for each template, R1 then R2, so
    // the file is positional in exactly the way the single-end file is.
    for (index, part) in parts.iter().enumerate() {
        let opt = part_opt(map_opt, part);
        let path = split::split_tmp_path(&prefix, index);
        let file = split::create_split_tmp(&prefix, index, part).map_err(|e| AlignError::io(&path, e))?;
        let mut writer = BufWriter::with_capacity(1 << 20, file);
        let with_cigar = opt.flag.contains(MapFlags::CIGAR);

        let mut pairs = PairReader::open(reads1, reads2)?;
        while let Some(batch) = pairs.next_batch(opt.mini_batch_size)? {
            if cancel() {
                return Err(AlignError::Cancelled);
            }
            for results in map_frag_batch(&pool, part, &opt, &batch) {
                for result in results {
                    let block = split::SplitQueryRecord {
                        n_reg: result.regs.len() as i32,
                        rep_len: result.rep_len,
                        frag_gap: result.frag_gap,
                        regs: result.regs,
                    };
                    split::write_split_query_record(&mut writer, &block, with_cigar)
                        .map_err(|e| AlignError::io(&path, e))?;
                }
            }
        }
        writer.flush().map_err(|e| AlignError::io(&path, e))?;
        cleanup.parts = index + 1;
        progress(0, index + 1, parts.len());
    }

    merge_pairs(
        &prefix,
        parts.len(),
        &merged_header,
        &rid_shifts,
        reads1,
        reads2,
        out,
        map_opt,
        params,
        cancel,
        progress,
    )
}

/// Merge each end across parts, re-pair, then emit.
#[allow(clippy::too_many_arguments)]
fn merge_pairs(
    prefix: &str,
    parts: usize,
    merged_header: &MmIdx,
    rid_shifts: &[u32],
    reads1: &Path,
    reads2: &Path,
    out: &Path,
    map_opt: &MapOpt,
    params: &MapParams,
    cancel: CancelFn<'_>,
    progress: ProgressFn<'_>,
) -> Result<MapStats, AlignError> {
    let mut opt = map_opt.clone();
    mapopt_update(&mut opt, merged_header);
    let with_cigar = opt.flag.contains(MapFlags::CIGAR);

    let mut readers = Vec::with_capacity(parts);
    for part in 0..parts {
        let path = split::split_tmp_path(prefix, part);
        let file = std::fs::File::open(&path).map_err(|e| AlignError::io(&path, e))?;
        let mut r = BufReader::with_capacity(1 << 20, file);
        split::read_split_header(&mut r).map_err(|e| AlignError::io(&path, e))?;
        readers.push((path, r));
    }

    let mut writer = open_output(out, merged_header, params)?;
    let mut pairs = PairReader::open(reads1, reads2)?;
    let mut stats = MapStats {
        parts,
        ..Default::default()
    };

    while let Some(batch) = pairs.next_batch(opt.mini_batch_size)? {
        if cancel() {
            return Err(AlignError::Cancelled);
        }
        for (r1, r2) in &batch {
            // Two blocks for each template and each part, in the order the writer made them.
            let mut blocks1 = Vec::with_capacity(parts);
            let mut blocks2 = Vec::with_capacity(parts);
            for (path, reader) in readers.iter_mut() {
                blocks1
                    .push(split::read_split_query_record(reader, with_cigar).map_err(|e| AlignError::io(&*path, e))?);
                blocks2
                    .push(split::read_split_query_record(reader, with_cigar).map_err(|e| AlignError::io(&*path, e))?);
            }

            let m1 = split::merge_split_query_records(&blocks1, rid_shifts, &opt, merged_header.k, r1.l_seq as i32);
            let m2 = split::merge_split_query_records(&blocks2, rid_shifts, &opt, merged_header.k, r2.l_seq as i32);
            let mut res1 = MapResult {
                regs: m1.regs,
                rep_len: m1.rep_len,
                frag_gap: m1.frag_gap,
            };
            let mut res2 = MapResult {
                regs: m2.regs,
                rep_len: m2.rep_len,
                frag_gap: m2.frag_gap,
            };

            // Restore orientation only now. The blocks for each part are in the flipped space
            // that the mapper worked in, and the merge operates on those coordinates.
            let (rev1, rev2) = orient_flags(&opt);
            restore_orientation(&mut res1, r1.l_seq as i32, rev1);
            restore_orientation(&mut res2, r2.l_seq as i32, rev2);

            // Pair again *after* the merge. The merge built the regions of each end again from
            // the start, so any pairing that the passes for each part made is gone. Pairing that
            // the code decides for each part would use only a fraction of the genome anyway. That
            // is the error the merge exists to undo.
            repair(&opt, &mut res1, &mut res2, r1, r2);
            emit_pair(&mut writer, merged_header, &opt, r1, r2, &res1, &res2, out, &mut stats)?;
        }
        progress(stats.queries, parts, parts);
    }

    writer.finish(out)?;
    progress(stats.queries, parts, parts);
    Ok(stats)
}

/// Map a batch of templates as fragments, and keep the input order.
///
/// Each element is the result for each segment of that template, R1 then R2. Both ends go in
/// together, because that is what lets a mate with a confident position inform a mate that is
/// ambiguous. To map them apart and then join the results throws that away.
///
/// Results come back in **flipped** coordinate space when the library orientation needs it. A
/// caller restores them at the right moment, which is different in the whole-index path and the
/// split path.
fn map_frag_batch(
    pool: &rayon::ThreadPool,
    index: &MmIdx,
    opt: &MapOpt,
    batch: &[(BseqRecord, BseqRecord)],
) -> Vec<Vec<MapResult>> {
    use rayon::prelude::*;
    pool.install(|| {
        batch
            .par_iter()
            .map(|(r1, r2)| {
                let (s1, s2, _, _) = orient(opt, &r1.seq, &r2.seq);
                minimap2::map::map_frag_queries(index, opt, &r1.name, &[&s1, &s2])
            })
            .collect()
    })
}

/// Turn a mapped batch into finished records, in the pool, and keep the batch order.
///
/// The orientation restore moves in here with the record build. It is work for each template that
/// also ran on the write thread. `collect` into a `Result` keeps the first error, and keeps the
/// output in order. So a failure reads the same as it did when this was a serial loop.
#[allow(clippy::too_many_arguments)]
fn build_batch_records(
    pool: &rayon::ThreadPool,
    index: &MmIdx,
    opt: &MapOpt,
    header: &sam::Header,
    batch: &[(BseqRecord, BseqRecord)],
    mapped: Vec<Vec<MapResult>>,
    out: &Path,
) -> Result<Vec<(Vec<RecordBuf>, MapStats)>, AlignError> {
    use rayon::prelude::*;
    pool.install(|| {
        batch
            .par_iter()
            .zip(mapped.into_par_iter())
            .map(|((r1, r2), mut results)| {
                // Fragment mapping returns one result for each segment, R1 then R2. A pop in
                // reverse keeps that association. A segment that is absent (which must not happen)
                // becomes "this end mapped nowhere", and does not move the pairing.
                let mut res2 = results.pop().unwrap_or_else(empty_result);
                let mut res1 = results.pop().unwrap_or_else(empty_result);
                let (rev1, rev2) = orient_flags(opt);
                restore_orientation(&mut res1, r1.l_seq as i32, rev1);
                restore_orientation(&mut res2, r2.l_seq as i32, rev2);

                // No `repair` here, and that is deliberate. `map_frag_queries` already ran
                // `pe::pair` over the fragment, so the ends arrive in a pair, with `proper_frag`,
                // MAPQ and `sam_pri` all set. A second pair step *clears* `proper_frag`. It scores
                // against a fragment gap that matters only in the split path, where the merge
                // discarded the original pairing. A call here looks like the correct change, and
                // it drops 0x2 from every record with no warning.
                build_pair_records(index, opt, header, r1, r2, &res1, &res2, out)
            })
            .collect()
    })
}

/// Which ends the code flips for mapping, from the library orientation of the preset.
fn orient_flags(opt: &MapOpt) -> (bool, bool) {
    ((opt.pe_ori >> 1) & 1 != 0, opt.pe_ori & 1 != 0)
}

/// Put both ends into the orientation the fragment mapper expects, from the library orientation
/// of the preset (`pe_ori`).
///
/// This is easy to miss, and it fails with no message. `sr` sets `pe_ori = 1`, which is FR. R2
/// arrives as the reverse complement of R1. The code must flip it, so that both ends read the same
/// way before the chaining step.
///
/// Skip that flip and the ends still *map*: coordinates, strands and mate fields all come out
/// right. But no pair is ever concordant, so nothing ever sets `proper_frag`, and every record
/// loses its 0x2 flag.
///
/// Returns the sequences, which it may have flipped, and which ends it flipped. That second value
/// is for [`restore_orientation`].
fn orient(opt: &MapOpt, seq1: &[u8], seq2: &[u8]) -> (Vec<u8>, Vec<u8>, bool, bool) {
    let rev1 = (opt.pe_ori >> 1) & 1 != 0;
    let rev2 = opt.pe_ori & 1 != 0;
    (flip(seq1, rev1), flip(seq2, rev2), rev1, rev2)
}

fn flip(seq: &[u8], revcomp: bool) -> Vec<u8> {
    let mut out = seq.to_vec();
    if revcomp {
        minimap2::seq::revcomp_ascii(&mut out);
    }
    out
}

/// Undo [`orient`] on the results. Coordinates and strands then describe the read as the caller
/// gave it, and not the flipped copy the mapper saw.
fn restore_orientation(result: &mut MapResult, qlen: i32, was_flipped: bool) {
    if !was_flipped {
        return;
    }
    for r in &mut result.regs {
        let old_qs = r.qs;
        r.qs = qlen - r.qe;
        r.qe = qlen - old_qs;
        r.rev = !r.rev;
        if let Some(extra) = &mut r.extra {
            extra.trans_strand = match extra.trans_strand {
                1 => 2,
                2 => 1,
                other => other,
            };
        }
    }
}

/// A result with no alignments. `MapResult` has no `Default`, and the fields it needs are not
/// clearly zero, so this file writes the value out one time.
fn empty_result() -> MapResult {
    MapResult {
        regs: Vec::new(),
        rep_len: 0,
        frag_gap: 0,
    }
}

/// Run the pairing step over an already-mapped pair.
fn repair(opt: &MapOpt, res1: &mut MapResult, res2: &mut MapResult, r1: &BseqRecord, r2: &BseqRecord) {
    // `pe::pair` reads the alignment extras of each end. Without base-level alignment on both
    // ends there is nothing to score a pairing against, and upstream skips it on the same
    // condition.
    if res1.regs.is_empty() || res2.regs.is_empty() || res1.regs[0].extra.is_none() || res2.regs[0].extra.is_none() {
        return;
    }
    let qlens = [r1.l_seq as i32, r2.l_seq as i32];
    let mut n_regs = [res1.regs.len(), res2.regs.len()];
    let mut regs = [std::mem::take(&mut res1.regs), std::mem::take(&mut res2.regs)];
    minimap2::pe::pair(
        res2.frag_gap,
        opt.pe_bonus,
        opt.a * 2 + opt.b,
        opt.a,
        &qlens,
        &mut n_regs,
        &mut regs,
    );
    let [regs1, regs2] = regs;
    res1.regs = regs1;
    res2.regs = regs2;
}

// ---- SAM emission ---------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn emit_pair(
    writer: &mut AlignmentWriter,
    index: &MmIdx,
    opt: &MapOpt,
    r1: &BseqRecord,
    r2: &BseqRecord,
    res1: &MapResult,
    res2: &MapResult,
    out: &Path,
    stats: &mut MapStats,
) -> Result<(), AlignError> {
    emit_end(writer, index, opt, r1, res1, res2, true, out, stats)?;
    emit_end(writer, index, opt, r2, res2, res1, false, out, stats)?;
    Ok(())
}

/// Every SAM record for one template, and the writer does not take part.
///
/// This is the work that used to happen on the write thread. minimap2 makes a SAM line from each
/// alignment. The code parses that line back into a typed record, and then fills in the paired
/// fields. The work for each record is independent, so it belongs in the mapping pool.
///
/// On the write thread it made mapping a strict read → map → write cycle, and the pool was idle
/// for two of the three phases. A WGS run measured ~26% pool use, and the mapping itself was only
/// a quarter of a stage that took 3 h 40 m.
#[allow(clippy::too_many_arguments)]
fn build_pair_records(
    index: &MmIdx,
    opt: &MapOpt,
    header: &sam::Header,
    r1: &BseqRecord,
    r2: &BseqRecord,
    res1: &MapResult,
    res2: &MapResult,
    out: &Path,
) -> Result<(Vec<RecordBuf>, MapStats), AlignError> {
    let mut records = Vec::new();
    let mut stats = MapStats::default();
    build_end_records(index, opt, header, r1, res1, res2, true, out, &mut records, &mut stats)?;
    build_end_records(index, opt, header, r2, res2, res1, false, out, &mut records, &mut stats)?;
    Ok((records, stats))
}

/// [`build_pair_records`] for one end. Mirrors [`emit_end`], which the split path still uses.
#[allow(clippy::too_many_arguments)]
fn build_end_records(
    index: &MmIdx,
    opt: &MapOpt,
    header: &sam::Header,
    record: &BseqRecord,
    own: &MapResult,
    mate: &MapResult,
    is_first: bool,
    out: &Path,
    records: &mut Vec<RecordBuf>,
    stats: &mut MapStats,
) -> Result<(), AlignError> {
    let qname = strip_mate_suffix(&record.name);
    let mate_primary = primary(mate);
    stats.queries += 1;

    if own.regs.is_empty() {
        let line = minimap2::format::sam::write_sam_record(
            index,
            qname,
            &record.seq,
            &record.qual,
            None,
            0,
            &[],
            opt.flag,
            own.rep_len,
        );
        let mut rec = crate::output::parse_record(&line, header, out)?;
        set_pair_fields(&mut rec, None, mate_primary, is_first);
        records.push(rec);
        stats.unmapped += 1;
        return Ok(());
    }

    for reg in own.regs.iter().filter(|reg| emits_record(opt, reg)) {
        let line = minimap2::format::sam::write_sam_record(
            index,
            qname,
            &record.seq,
            &record.qual,
            Some(reg),
            own.regs.len(),
            &own.regs,
            opt.flag,
            own.rep_len,
        );
        let mut rec = crate::output::parse_record(&line, header, out)?;
        set_pair_fields(&mut rec, Some(reg), mate_primary, is_first);
        records.push(rec);
    }
    stats.mapped += 1;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_end(
    writer: &mut AlignmentWriter,
    index: &MmIdx,
    opt: &MapOpt,
    record: &BseqRecord,
    own: &MapResult,
    mate: &MapResult,
    is_first: bool,
    out: &Path,
    stats: &mut MapStats,
) -> Result<(), AlignError> {
    let qname = strip_mate_suffix(&record.name);
    let mate_primary = primary(mate);
    stats.queries += 1;

    if own.regs.is_empty() {
        let line = minimap2::format::sam::write_sam_record(
            index,
            qname,
            &record.seq,
            &record.qual,
            None,
            0,
            &[],
            opt.flag,
            own.rep_len,
        );
        writer.write_line_with(&line, out, |rec, _| set_pair_fields(rec, None, mate_primary, is_first))?;
        stats.unmapped += 1;
        return Ok(());
    }

    for reg in own.regs.iter().filter(|reg| emits_record(opt, reg)) {
        let line = minimap2::format::sam::write_sam_record(
            index,
            qname,
            &record.seq,
            &record.qual,
            Some(reg),
            own.regs.len(),
            &own.regs,
            opt.flag,
            own.rep_len,
        );
        writer.write_line_with(&line, out, |rec, _| {
            set_pair_fields(rec, Some(reg), mate_primary, is_first)
        })?;
    }
    stats.mapped += 1;
    Ok(())
}

/// The mapper's own rule for which regions reach the output.
///
/// A region whose `parent` is not itself is a *secondary* alignment. It is another place the read
/// could have gone. `NO_PRINT_2ND` says not to emit those, and the `sr` preset sets it, because
/// for short reads MAPQ already carries the ambiguity.
///
/// This module must apply the rule here, because Navigator builds records itself.
/// `minimap2-pure-rs` keeps its PE SAM assembly private. So this code does not go through the
/// pipeline that would apply the rule: `minimap2::map`, on the `r.id != r.parent` test.
///
/// Without the rule, every alternative placement went into the file. On a targeted-Y sample that
/// was 404 million secondary records against 62 million primaries, or 86.6% of the file, and none
/// carried SEQ. It grew the alignment to 17 GB, and every later stage had to read all of them.
///
/// This keeps supplementary alignments. Their `parent` *is* themselves, they carry sequence, and
/// they are how the format shows a split read.
pub(crate) fn emits_record(opt: &MapOpt, reg: &AlignReg) -> bool {
    !(opt.flag.contains(MapFlags::NO_PRINT_2ND) && reg.id != reg.parent)
}

/// The region a mate is "at" for `RNEXT` and `PNEXT`: its primary alignment.
fn primary(result: &MapResult) -> Option<&AlignReg> {
    result.regs.iter().find(|r| r.sam_pri).or_else(|| result.regs.first())
}

/// Fill in the paired half of a record: flags, `RNEXT`, `PNEXT`, `TLEN`.
///
/// The single-end writer made everything else. This code used to patch the SAM text by column
/// position. It now changes a typed [`RecordBuf`], so a mate position can not land in the
/// template-length field, whatever the layout of the formatter does.
///
/// `own` is this record's region (`None` for an unmapped read) and `mate` is the mate's primary.
fn set_pair_fields(record: &mut RecordBuf, own: Option<&AlignReg>, mate: Option<&AlignReg>, is_first: bool) {
    let mut flags = record.flags().bits();
    flags |= flag::PAIRED | if is_first { flag::FIRST } else { flag::LAST };
    match mate {
        Some(m) => {
            if m.rev {
                flags |= flag::MATE_REVERSE;
            }
        }
        None => flags |= flag::MATE_UNMAPPED,
    }
    if let (Some(o), Some(m)) = (own, mate) {
        // `proper_frag` is `pe::pair`'s verdict that these two ends form a concordant fragment; it
        // is the only thing entitled to set 0x2.
        if o.proper_frag && m.proper_frag {
            flags |= flag::PROPER_PAIR;
        }
    }
    *record.flags_mut() = Flags::from(flags);

    match (own, mate) {
        (Some(_), Some(m)) => {
            *record.mate_reference_sequence_id_mut() = Some(m.rid as usize);
            *record.mate_alignment_start_mut() = Position::new(m.rs as usize + 1);
        }
        (Some(o), None) => {
            // By convention, an unmapped mate goes at the locus of this record, so the pair
            // stays together after a coordinate sort of the file.
            *record.mate_reference_sequence_id_mut() = Some(o.rid as usize);
            *record.mate_alignment_start_mut() = Position::new(o.rs as usize + 1);
        }
        (None, Some(m)) => {
            // Likewise in reverse: place the unmapped read at its mapped mate's locus.
            *record.reference_sequence_id_mut() = Some(m.rid as usize);
            *record.alignment_start_mut() = Position::new(m.rs as usize + 1);
            *record.mate_reference_sequence_id_mut() = Some(m.rid as usize);
            *record.mate_alignment_start_mut() = Position::new(m.rs as usize + 1);
        }
        (None, None) => {
            *record.mate_reference_sequence_id_mut() = None;
            *record.mate_alignment_start_mut() = None;
        }
    }

    *record.template_length_mut() = match (own, mate) {
        // Only meaningful when both ends sit on the same reference sequence.
        (Some(o), Some(m)) if o.rid == m.rid => tlen(o, m) as i32,
        _ => 0,
    };
}

/// Signed observed template length: the span from the leftmost start to the rightmost end,
/// negative for whichever end is rightmost. Zero when both ends start at the same base, since
/// neither is leftmost and SAM has no way to break the tie consistently.
fn tlen(own: &AlignReg, mate: &AlignReg) -> i64 {
    let start = own.rs.min(mate.rs) as i64;
    let end = own.re.max(mate.re) as i64;
    let span = end - start;
    match own.rs.cmp(&mate.rs) {
        std::cmp::Ordering::Less => span,
        std::cmp::Ordering::Greater => -span,
        std::cmp::Ordering::Equal => 0,
    }
}

/// Both ends of a template must share a QNAME, so a `/1` or `/2` at the end has to go.
///
/// Our own revert stage writes bare names. But vendor FASTQ often carries the suffix, and a QNAME
/// that does not match breaks pairing for every downstream tool, with no warning.
fn strip_mate_suffix(name: &str) -> &str {
    let bytes = name.as_bytes();
    if bytes.len() > 2 && bytes[bytes.len() - 2] == b'/' {
        let last = bytes[bytes.len() - 1];
        if last == b'1' || last == b'2' {
            return &name[..name.len() - 2];
        }
    }
    name
}

// ---- paired input ---------------------------------------------------------

/// Reads two FASTQ files in lockstep.
struct PairReader {
    left: BseqFile,
    right: BseqFile,
    left_path: std::path::PathBuf,
    right_path: std::path::PathBuf,
}

impl PairReader {
    fn open(reads1: &Path, reads2: &Path) -> Result<Self, AlignError> {
        Ok(Self {
            left: BseqFile::open(&path_str_of(reads1)?).map_err(|e| AlignError::io(reads1, e))?,
            right: BseqFile::open(&path_str_of(reads2)?).map_err(|e| AlignError::io(reads2, e))?,
            left_path: reads1.to_path_buf(),
            right_path: reads2.to_path_buf(),
        })
    }

    /// The next batch of templates, or `None` at end of input.
    ///
    /// This reads the two files **one record at a time in step**, and collects records until it
    /// reaches the base budget.
    ///
    /// There is a simpler implementation: ask each file for a batch, and zip the results. It is
    /// wrong, and it looks right on tidy data. The reader below batches by *bases*, so two files
    /// whose reads differ in length give different record counts from the same budget. Real data
    /// has a tail of shorter reads, which adapter trimming and quality trimming make. So a
    /// mismatch of 332,653 against 332,722 appears on a real WGS, and never on a fixture where
    /// every read is the same length.
    ///
    /// A file that truly stops before the other is still an error, and not a truncation. R1 and R2
    /// that have lost their order would pair every later read with the wrong mate. That is much
    /// worse than a refusal to run.
    fn next_batch(&mut self, chunk: i64) -> Result<Option<Vec<(BseqRecord, BseqRecord)>>, AlignError> {
        let mut batch = Vec::new();
        let mut bases: i64 = 0;

        loop {
            let left = self
                .left
                .read_record_with_qual(true)
                .map_err(|e| AlignError::io(&self.left_path, e))?;
            let right = self
                .right
                .read_record_with_qual(true)
                .map_err(|e| AlignError::io(&self.right_path, e))?;

            match (left, right) {
                (Some(a), Some(b)) => {
                    bases += a.l_seq as i64 + b.l_seq as i64;
                    batch.push((a, b));
                }
                (None, None) => break,
                (a, _) => {
                    return Err(AlignError::Message(format!(
                        "{} and {} are not in lockstep — {} ran out first, so the remaining mates \
                         would be paired wrongly",
                        self.left_path.display(),
                        self.right_path.display(),
                        if a.is_none() {
                            "the first file"
                        } else {
                            "the second file"
                        },
                    )));
                }
            }

            if bases >= chunk {
                break;
            }
        }

        Ok((!batch.is_empty()).then_some(batch))
    }
}

#[cfg(test)]
mod tests;
