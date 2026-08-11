//! Paired-end mapping — the `sr` path, and so most vendor WGS.
//!
//! A pair is not two independent reads. Mapping them together lets a confidently-placed mate
//! rescue an ambiguous one, and the fragment's expected span is evidence about where the second
//! end belongs; both feed MAPQ. So the two ends are mapped as one *fragment*
//! ([`minimap2::map::map_frag_queries`]) and then paired ([`minimap2::pe::pair`]), which is what
//! sets `proper_frag`, adjusts MAPQ, and decides which region of each end is the primary.
//!
//! ## What this module has to build itself
//!
//! The pieces above are public in `minimap2-pure-rs`. Its **PE SAM formatting is not** — that
//! lives in private `pipeline.rs` helpers, so the mate-facing half of each record is assembled
//! here: the paired flags, `RNEXT`/`PNEXT`, and `TLEN`.
//!
//! The approach mirrors upstream's: format the single-end line with the public writer (which
//! already knows CIGAR, clipping, and tags), then fill in the paired fields. That is string
//! surgery on a formatted SAM line, which is worth naming rather than hiding — but the
//! alternative is reimplementing CIGAR and tag emission, which is far more of the delicate work,
//! not less. [`set_pair_fields`] is deliberately small and heavily tested for that reason.
//!
//! ## Split indexes
//!
//! Same shape as the single-end path: map the fragment against each part, spill per-segment hit
//! blocks, then merge each end across parts and re-pair. The re-pair after the merge is essential
//! — pairing decided per part would be based on a fraction of the genome, exactly the error the
//! merge exists to prevent.

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
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::RecordBuf;

use crate::error::AlignError;
use crate::map::{
    append_header, header_only, open_index, open_output, part_opt, path_str_of, prepared_map_opt, thread_pool,
    CancelFn, MapParams, MapStats, ProgressFn, ScratchGuard,
};
use crate::output::AlignmentWriter;

/// SAM flag bits this module sets. Named because a bare `0x20` in flag arithmetic is unreadable
/// and the difference between `0x20` and `0x10` is a silently wrong strand.
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
/// The two files must be in lockstep — record *n* of each is one template — which is what
/// `navigator-analysis`'s revert stage guarantees for the FASTQ it produces.
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

    let mut writer = open_output(out, index, params)?;
    let mut pairs = PairReader::open(reads1, reads2)?;
    let pool = thread_pool(params)?;
    let mut stats = MapStats {
        parts: 1,
        ..Default::default()
    };

    while let Some(batch) = pairs.next_batch(opt.mini_batch_size)? {
        if cancel() {
            return Err(AlignError::Cancelled);
        }
        let mapped = map_frag_batch(&pool, index, &opt, &batch);
        for ((r1, r2), mut results) in batch.iter().zip(mapped) {
            // Fragment mapping returns one result per segment, R1 then R2. Popping in reverse
            // keeps that association; an absent segment (which should not happen) degrades to
            // "this end mapped nowhere" rather than shifting the pairing.
            let mut res2 = results.pop().unwrap_or_else(empty_result);
            let mut res1 = results.pop().unwrap_or_else(empty_result);
            let (rev1, rev2) = orient_flags(&opt);
            restore_orientation(&mut res1, r1.l_seq as i32, rev1);
            restore_orientation(&mut res2, r2.l_seq as i32, rev2);

            // No `repair` here, deliberately. `map_frag_queries` already ran `pe::pair` over the
            // fragment, so the ends arrive paired: `proper_frag`, MAPQ, and `sam_pri` are set.
            // Pairing them a second time *clears* `proper_frag` — the re-pair is scored against a
            // fragment gap that only means something in the split path, where merging discarded
            // the original pairing. Adding a call here is the obvious-looking change that silently
            // drops 0x2 from every record.
            emit_pair(&mut writer, index, &opt, r1, r2, &res1, &res2, out, &mut stats)?;
        }
        progress(stats.queries, 1, 1);
    }

    writer.finish(out)?;
    progress(stats.queries, 1, 1);
    Ok(stats)
}

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
    // One scratch file per part; each holds two blocks per template, R1 then R2, so the file is
    // positional in exactly the way the single-end path's is.
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
            // Two blocks per template per part, in the order they were written.
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

            // Restore orientation only now: the per-part blocks were written in the flipped space
            // the mapper worked in, and the merge operates on those coordinates.
            let (rev1, rev2) = orient_flags(&opt);
            restore_orientation(&mut res1, r1.l_seq as i32, rev1);
            restore_orientation(&mut res2, r2.l_seq as i32, rev2);

            // Re-pair *after* merging: the merge rebuilt each end's regions from scratch, so
            // whatever pairing the per-part passes established is gone. Pairing decided per part
            // would rest on a fraction of the genome anyway — the error the merge exists to undo.
            repair(&opt, &mut res1, &mut res2, r1, r2);
            emit_pair(&mut writer, merged_header, &opt, r1, r2, &res1, &res2, out, &mut stats)?;
        }
        progress(stats.queries, parts, parts);
    }

    writer.finish(out)?;
    progress(stats.queries, parts, parts);
    Ok(stats)
}

/// Map a batch of templates as fragments, preserving input order.
///
/// Each element is that template's per-segment results, R1 then R2. Both ends go in together
/// because that is what lets a confidently-placed mate inform an ambiguous one — mapping them
/// separately and reconciling afterwards throws that away.
///
/// Results come back in **flipped** coordinate space when the library orientation calls for it;
/// callers restore at the right moment, which differs between the whole-index and split paths.
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

/// Whether each end is flipped for mapping, from the preset's library orientation.
fn orient_flags(opt: &MapOpt) -> (bool, bool) {
    ((opt.pe_ori >> 1) & 1 != 0, opt.pe_ori & 1 != 0)
}

/// Put both ends into the orientation the fragment mapper expects, per the preset's library
/// orientation (`pe_ori`).
///
/// This is easy to miss and fails quietly. `sr` sets `pe_ori = 1`, meaning FR: R2 arrives
/// reverse-complemented relative to R1, and must be flipped so both ends read the same way before
/// the fragment is chained. Skip it and the ends still *map* — coordinates, strands, and mate
/// fields all come out right — but no pair is ever judged concordant, so `proper_frag` is never
/// set and every record loses its 0x2 flag.
///
/// Returns the possibly-flipped sequences and whether each was flipped, for [`restore_orientation`].
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

/// Undo [`orient`] on the results, so coordinates and strands describe the read as it was given
/// to us rather than the flipped copy the mapper saw.
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

/// A result with no alignments. `MapResult` has no `Default`, and the fields it would need are
/// not obviously zero, so this is spelled out once.
fn empty_result() -> MapResult {
    MapResult {
        regs: Vec::new(),
        rep_len: 0,
        frag_gap: 0,
    }
}

/// Run the pairing step over an already-mapped pair.
fn repair(opt: &MapOpt, res1: &mut MapResult, res2: &mut MapResult, r1: &BseqRecord, r2: &BseqRecord) {
    // `pe::pair` reads each end's alignment extras; without base-level alignment on both ends
    // there is nothing to score a pairing against, and upstream skips it on the same condition.
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

    for reg in &own.regs {
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

/// The region a mate is "at" for the purposes of `RNEXT`/`PNEXT` — its primary alignment.
fn primary(result: &MapResult) -> Option<&AlignReg> {
    result.regs.iter().find(|r| r.sam_pri).or_else(|| result.regs.first())
}

/// Fill in the paired half of a record: flags, `RNEXT`, `PNEXT`, `TLEN`.
///
/// The single-end writer produced everything else. This used to patch the formatted SAM text by
/// column position; it now mutates a typed [`RecordBuf`], so a mate position cannot end up in the
/// template-length field however the formatter's layout changes.
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
            // An unmapped mate is conventionally reported at this record's own locus, so the pair
            // stays together once the file is coordinate-sorted.
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

/// Both ends of a template must share a QNAME, so a trailing `/1` or `/2` has to go.
///
/// Our own revert stage writes bare names, but vendor FASTQ frequently carries the suffix, and a
/// mismatched QNAME silently breaks pairing for every downstream tool.
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
    /// Reads the two files **one record at a time in step**, accumulating until the base budget is
    /// reached. The obvious implementation — ask each file for a batch and zip the results — is
    /// wrong, and wrong in a way that looks right on tidy data: the underlying reader batches by
    /// *bases*, so two files whose reads differ in length yield different record counts from the
    /// same budget. Real data has a tail of shorter reads from adapter and quality trimming, so a
    /// 332,653-against-332,722 mismatch appears on a genuine WGS and never on a fixture where every
    /// read is the same length.
    ///
    /// A file genuinely ending before the other is still an error rather than a truncation: R1/R2
    /// that have drifted out of step would pair every later read with the wrong mate, which is far
    /// worse than refusing to run.
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
