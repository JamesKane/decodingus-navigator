//! Coverage + callable-loci walker — the Rust port of the Scala
//! `CoverageCallableWalker` (which itself replaces GATK `CollectWgsMetrics` +
//! `CallableLoci` over htsjdk). Single pass over a coordinate-sorted BAM/CRAM builds,
//! per main-assembly contig: a depth histogram, per-position callable state, and
//! samtools-style coverage stats.
//!
//! Parity target is the Scala walker, not samtools — notably mean base/mapping quality
//! are averaged **per base observation** (Σ quality / Σ depth), where samtools averages
//! per read. See the crate fixture tests for hand-computed expected values.
//!
//! Memory: a **sliding-window pileup** finalizes each position once the read frontier
//! passes it, so peak memory is the span of currently-open reads — not the contig
//! length. The only contig-sized allocation is the reference bases of the contig being
//! walked (one at a time, for N detection); whole-genome HG002 peaks ~2 GB vs the
//! ~84 GB a dense per-contig-arrays approach would need. Requires a coordinate-sorted
//! BAM. (Streaming the reference in windows too is a further optimization.) BED-interval
//! output and progress callbacks from the Scala walker are deferred.

use std::collections::{HashMap, HashSet, VecDeque};
use std::path::Path;

use noodles::core::Region;
use noodles::fasta;

use serde::{Deserialize, Serialize};

use crate::contig;
use crate::error::AnalysisError;
use crate::reader;

/// Algorithm version for the coverage artifact cache key; bump on any change that
/// alters output (plan §6 cache versioning).
pub const COVERAGE_VERSION: &str = "coverage-1";

/// Callable-loci parameters. Defaults match GATK `CallableLoci` (and the Scala walker).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallableLociParams {
    pub min_depth: u32,
    pub max_depth: Option<u32>,
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub max_low_mapq: u8,
    pub max_fraction_low_mapq: f64,
}

impl Default for CallableLociParams {
    fn default() -> Self {
        CallableLociParams {
            min_depth: 4,
            max_depth: None,
            min_mapping_quality: 10,
            min_base_quality: 20,
            max_low_mapq: 1,
            max_fraction_low_mapq: 0.1,
        }
    }
}

/// Per-position callable classification (GATK `CallableLoci` states). Hierarchical:
/// the first failing condition wins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CallableState {
    RefN,
    NoCoverage,
    PoorMappingQuality,
    LowCoverage,
    ExcessiveCoverage,
    Callable,
}

/// Per-contig callable-state base counts.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContigCallableMetrics {
    pub contig: String,
    pub ref_n: u64,
    pub callable: u64,
    pub no_coverage: u64,
    pub low_coverage: u64,
    pub excessive_coverage: u64,
    pub poor_mapping_quality: u64,
}

/// Per-contig samtools-`coverage`-style stats (averaged per base observation; see module docs).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContigCoverageStats {
    pub contig: String,
    pub start_pos: u64, // always 1
    pub end_pos: u64,   // contig length
    pub num_reads: u64,
    pub cov_bases: u64,
    pub coverage: f64, // percent of contig with depth > 0
    pub mean_depth: f64,
    pub mean_base_q: f64,
    pub mean_map_q: f64,
}

/// Combined coverage + callable result (replaces the Scala `CoverageCallableResult`'s
/// global-metrics + callable-summary + samtools-stats fields).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CoverageResult {
    pub genome_territory: u64,
    pub mean_coverage: f64,
    pub median_coverage: f64,
    pub sd_coverage: f64,
    /// Depth histogram, clamped at index 255.
    pub coverage_histogram: Vec<u64>,
    pub pct_1x: f64,
    pub pct_5x: f64,
    pub pct_10x: f64,
    pub pct_15x: f64,
    pub pct_20x: f64,
    pub pct_25x: f64,
    pub pct_30x: f64,
    pub pct_40x: f64,
    pub pct_50x: f64,
    pub callable_bases: u64,
    pub contig_callable: Vec<ContigCallableMetrics>,
    pub contig_coverage_stats: Vec<ContigCoverageStats>,
}

const HIST_LEN: usize = 256;

/// Per-position accumulator within the sliding pileup window.
#[derive(Clone, Default)]
struct Col {
    depth: u32,
    base_q_sum: u64,
    map_q_sum: u64,
    qc_pass: u32,
    low_mapq: u32,
}

/// Global accumulators, folded as positions finalize.
struct Globals {
    hist: Vec<u64>,
    n: u64,
    sum_depth: u128,
    sum_sq: u128,
}

impl Globals {
    fn new() -> Self {
        Globals { hist: vec![0; HIST_LEN], n: 0, sum_depth: 0, sum_sq: 0 }
    }
}

/// Finished per-contig output (kept tiny; assembled into the result in header order).
struct ContigOut {
    callable: ContigCallableMetrics,
    stats: ContigCoverageStats,
}

/// Streaming state for the contig currently being walked. Memory is bounded by the
/// sliding window (the span of currently-open reads), not the contig length.
struct CurContig {
    name: String,
    length: usize,
    ref_bases: Vec<u8>,
    window: VecDeque<Col>,
    emit_cursor: usize, // 1-based position next to finalize; window front aligns here
    read_count: u64,
    cm: ContigCallableMetrics,
    covered: u64,
    total_base_obs: u64,
    base_q_total: u64,
    map_q_total: u64,
    sum_depth: u128,
}

impl CurContig {
    fn new(name: String, length: usize, ref_bases: Vec<u8>) -> Self {
        let cm = ContigCallableMetrics {
            contig: name.clone(),
            ref_n: 0,
            callable: 0,
            no_coverage: 0,
            low_coverage: 0,
            excessive_coverage: 0,
            poor_mapping_quality: 0,
        };
        CurContig {
            name,
            length,
            ref_bases,
            window: VecDeque::new(),
            emit_cursor: 1,
            read_count: 0,
            cm,
            covered: 0,
            total_base_obs: 0,
            base_q_total: 0,
            map_q_total: 0,
            sum_depth: 0,
        }
    }

    /// Fold one finalized column (covered or empty) into global + contig accumulators.
    fn finalize_col(&mut self, pos: usize, col: Col, params: &CallableLociParams, g: &mut Globals) {
        let depth = col.depth;
        let clamped = depth.min(255) as usize;
        g.hist[clamped] += 1;
        g.n += 1;
        g.sum_depth += depth as u128;
        g.sum_sq += (depth as u128) * (depth as u128);
        self.sum_depth += depth as u128;
        if depth > 0 {
            self.covered += 1;
            self.total_base_obs += depth as u64;
            self.base_q_total += col.base_q_sum;
            self.map_q_total += col.map_q_sum;
        }
        let ref_base = self.ref_bases.get(pos - 1).copied().unwrap_or(b'N');
        match determine_state(ref_base, depth, col.qc_pass, col.low_mapq, params) {
            CallableState::RefN => self.cm.ref_n += 1,
            CallableState::NoCoverage => self.cm.no_coverage += 1,
            CallableState::PoorMappingQuality => self.cm.poor_mapping_quality += 1,
            CallableState::LowCoverage => self.cm.low_coverage += 1,
            CallableState::ExcessiveCoverage => self.cm.excessive_coverage += 1,
            CallableState::Callable => self.cm.callable += 1,
        }
    }

    /// Finalize all positions strictly before `target` (clamped to the contig end).
    fn advance_to(&mut self, target: usize, params: &CallableLociParams, g: &mut Globals) {
        while self.emit_cursor < target && self.emit_cursor <= self.length {
            let col = self.window.pop_front().unwrap_or_default();
            let pos = self.emit_cursor;
            self.finalize_col(pos, col, params, g);
            self.emit_cursor += 1;
        }
    }

    /// Add one covered base at 1-based `pos` (>= `emit_cursor`) to the window.
    fn add(&mut self, pos: usize, base_q: u8, mapq: u8, params: &CallableLociParams) {
        let idx = pos - self.emit_cursor;
        while self.window.len() <= idx {
            self.window.push_back(Col::default());
        }
        let col = &mut self.window[idx];
        col.depth += 1;
        col.base_q_sum += base_q as u64;
        col.map_q_sum += mapq as u64;
        if mapq >= params.min_mapping_quality && base_q >= params.min_base_quality {
            col.qc_pass += 1;
        }
        if mapq <= params.max_low_mapq {
            col.low_mapq += 1;
        }
    }

    /// Flush the remaining window + uncovered tail, then produce the contig output.
    fn finish(mut self, params: &CallableLociParams, g: &mut Globals) -> ContigOut {
        self.advance_to(self.length + 1, params, g);
        let length = self.length as f64;
        let stats = ContigCoverageStats {
            contig: self.name.clone(),
            start_pos: 1,
            end_pos: self.length as u64,
            num_reads: self.read_count,
            cov_bases: self.covered,
            coverage: if self.length == 0 { 0.0 } else { self.covered as f64 / length * 100.0 },
            mean_depth: if self.length == 0 { 0.0 } else { self.sum_depth as f64 / length },
            mean_base_q: if self.total_base_obs == 0 { 0.0 } else { self.base_q_total as f64 / self.total_base_obs as f64 },
            mean_map_q: if self.total_base_obs == 0 { 0.0 } else { self.map_q_total as f64 / self.total_base_obs as f64 },
        };
        ContigOut { callable: self.cm, stats }
    }
}

/// Single coordinate-ordered pass with a sliding pileup window: positions finalize once
/// the read frontier passes them, so peak memory is the open-read span, not the contig
/// length. Requires a coordinate-sorted BAM/CRAM (the standard genomics layout). The
/// reference is needed both to detect reference-N positions and to decode CRAM.
pub fn collect_coverage_callable(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<CoverageResult, AnalysisError> {
    collect_coverage_callable_with_progress(bam_path, reference_path, params, contig_allowlist, &mut |_, _| {})
}

/// Like [`collect_coverage_callable`], reporting `progress(contigs_done, contigs_total)` as each
/// tracked contig is finalized — so a whole-genome pass (minutes on a real WGS BAM) can drive a
/// progress bar instead of looking stalled. Assumes a coordinate-sorted BAM (contigs in order).
pub fn collect_coverage_callable_with_progress(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &mut dyn FnMut(usize, usize),
) -> Result<CoverageResult, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, Some(reference_path))?;

    // ref_id -> (name, length) for tracked (main-assembly, allowlisted) contigs.
    let tracked: Vec<Option<(String, usize)>> = header
        .reference_sequences()
        .iter()
        .map(|(name_bytes, map)| {
            let name = String::from_utf8_lossy(name_bytes.as_ref()).into_owned();
            let keep = contig::is_main_assembly(&name)
                && contig_allowlist.map_or(true, |set| set.contains(&name));
            keep.then(|| (name, map.length().get()))
        })
        .collect();

    // Total contigs we'll walk, for the progress denominator.
    let total_tracked = tracked.iter().filter(|o| o.is_some()).count();
    let mut contigs_done = 0usize;
    progress(0, total_tracked);

    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)
        .map_err(|e| AnalysisError::io(reference_path, e))?;

    let mut g = Globals::new();
    let mut finished: HashMap<usize, ContigOut> = HashMap::new();
    let mut cur: Option<(usize, CurContig)> = None;

    for result in reader.records(&header) {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }
        let ref_id = match record.reference_sequence_id() {
            Some(r) => r,
            None => continue,
        };
        let Some((name, length)) = tracked.get(ref_id).and_then(|o| o.as_ref()) else {
            continue;
        };

        // Contig transition: finalize the previous contig, load the new contig's ref.
        if cur.as_ref().map(|(id, _)| *id) != Some(ref_id) {
            if let Some((id, c)) = cur.take() {
                finished.insert(id, c.finish(params, &mut g));
                contigs_done += 1;
                progress(contigs_done, total_tracked);
            }
            let region: Region = name
                .parse()
                .map_err(|_| AnalysisError::Message(format!("bad region for contig {name}")))?;
            let rec = fasta_reader
                .query(&region)
                .map_err(|e| AnalysisError::io(reference_path, e))?;
            let ref_bases = rec.sequence().as_ref().to_vec();
            cur = Some((ref_id, CurContig::new(name.clone(), *length, ref_bases)));
        }

        let start = match record.alignment_start() {
            Some(p) => p.get(),
            None => continue,
        };
        let (_, c) = cur.as_mut().unwrap();
        c.advance_to(start, params, &mut g);
        c.read_count += 1;

        let mapq = record.mapping_quality().map_or(255u8, |m| m.get());
        let quals = record.quality_scores();
        let quals = quals.as_ref();

        let mut ref_pos = start; // 1-based
        let mut query_off = 0usize;
        for op in record.cigar().as_ref() {
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        if pos >= 1 && pos <= c.length {
                            let base_q = quals.get(query_off + i).copied().unwrap_or(0);
                            c.add(pos, base_q, mapq, params);
                        }
                    }
                    ref_pos += len;
                    query_off += len;
                }
                (true, false) => ref_pos += len,
                (false, true) => query_off += len,
                (false, false) => {}
            }
        }
    }
    if let Some((id, c)) = cur.take() {
        finished.insert(id, c.finish(params, &mut g));
        contigs_done += 1;
        progress(contigs_done, total_tracked);
    }

    // Assemble in header (ref_id) order; tracked contigs with no reads are zero-coverage.
    let mut contig_callable = Vec::new();
    let mut contig_stats = Vec::new();
    for (ref_id, opt) in tracked.iter().enumerate() {
        let Some((name, length)) = opt else { continue };
        if let Some(out) = finished.remove(&ref_id) {
            contig_callable.push(out.callable);
            contig_stats.push(out.stats);
        } else {
            // No reads: every position is depth 0 (RefN where the reference is N).
            let region: Region = name
                .parse()
                .map_err(|_| AnalysisError::Message(format!("bad region for contig {name}")))?;
            let rec = fasta_reader
                .query(&region)
                .map_err(|e| AnalysisError::io(reference_path, e))?;
            let ref_bases = rec.sequence();
            let ref_bases = ref_bases.as_ref();
            let mut ref_n: u64 = 0;
            for idx in 0..*length {
                let b = ref_bases.get(idx).copied().unwrap_or(b'N');
                if b == b'N' || b == b'n' {
                    ref_n += 1;
                }
            }
            g.hist[0] += *length as u64;
            g.n += *length as u64;
            contig_callable.push(ContigCallableMetrics {
                contig: name.clone(),
                ref_n,
                callable: 0,
                no_coverage: *length as u64 - ref_n,
                low_coverage: 0,
                excessive_coverage: 0,
                poor_mapping_quality: 0,
            });
            contig_stats.push(ContigCoverageStats {
                contig: name.clone(),
                start_pos: 1,
                end_pos: *length as u64,
                num_reads: 0,
                cov_bases: 0,
                coverage: 0.0,
                mean_depth: 0.0,
                mean_base_q: 0.0,
                mean_map_q: 0.0,
            });
        }
    }

    let mean = if g.n == 0 { 0.0 } else { g.sum_depth as f64 / g.n as f64 };
    let sd = if g.n < 2 {
        0.0
    } else {
        let var = g.sum_sq as f64 / g.n as f64 - mean * mean;
        var.max(0.0).sqrt()
    };
    let callable_bases = contig_callable.iter().map(|c| c.callable).sum();

    Ok(CoverageResult {
        genome_territory: g.n,
        mean_coverage: mean,
        median_coverage: median_from_hist(&g.hist, g.n),
        sd_coverage: sd,
        pct_1x: pct_at_least(&g.hist, g.n, 1),
        pct_5x: pct_at_least(&g.hist, g.n, 5),
        pct_10x: pct_at_least(&g.hist, g.n, 10),
        pct_15x: pct_at_least(&g.hist, g.n, 15),
        pct_20x: pct_at_least(&g.hist, g.n, 20),
        pct_25x: pct_at_least(&g.hist, g.n, 25),
        pct_30x: pct_at_least(&g.hist, g.n, 30),
        pct_40x: pct_at_least(&g.hist, g.n, 40),
        pct_50x: pct_at_least(&g.hist, g.n, 50),
        coverage_histogram: g.hist,
        callable_bases,
        contig_callable,
        contig_coverage_stats: contig_stats,
    })
}

/// Mean read length and mean fragment (template) length, sampled from the first ~50k
/// primary mapped reads. The molecule-length proxy for the self-referential callable
/// run-length gate (long reads → long molecules → long callable runs). Fragment falls
/// back to read length when templates are unpaired (e.g. long-read single-end).
pub fn estimate_molecule_lengths(bam_path: &Path, reference: Option<&Path>) -> Result<(f64, f64), AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, reference)?;

    let (mut n, mut read_sum, mut frag_n, mut frag_sum) = (0u64, 0u64, 0u64, 0u64);
    for result in reader.records(&header) {
        let record = result?;
        let f = record.flags();
        if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
            continue;
        }
        let len = record.sequence().len() as u64;
        if len == 0 {
            continue;
        }
        read_sum += len;
        n += 1;
        // Fragment length only from properly-paired reads, with a sanity cap — chimeric or
        // improper pairs carry enormous |TLEN| that would otherwise blow up the mean (and
        // the run-length gate). Single-end / long reads have no proper pairs -> read-length.
        if f.is_properly_segmented() {
            let tlen = record.template_length().unsigned_abs() as u64;
            if tlen > 0 && tlen < 100_000 {
                frag_sum += tlen;
                frag_n += 1;
            }
        }
        if n >= 50_000 {
            break;
        }
    }
    if n == 0 {
        return Ok((0.0, 0.0));
    }
    let read_len = read_sum as f64 / n as f64;
    let frag_len = if frag_n > 0 { frag_sum as f64 / frag_n as f64 } else { read_len };
    Ok((read_len, frag_len))
}

/// CALLABLE intervals (BED 0-based half-open) on one `contig`, coalesced and kept only
/// when the run is at least `min_run_len` bases. Reference-free: positions are classified
/// by depth / QC / MAPQ via the GATK hierarchy (reference-N regions carry no reads and
/// fall out as no-coverage). Memory is bounded by the open-read window; needs a BAM index.
pub fn callable_intervals(
    bam_path: &Path,
    contig: &str,
    params: &CallableLociParams,
    min_run_len: u32,
    reference: Option<&Path>,
) -> Result<Vec<(i64, i64)>, AnalysisError> {
    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    let query = reader.query(&header, &region)?;

    let mut window: VecDeque<Col> = VecDeque::new();
    let mut emit_cursor: usize = 1;
    let mut intervals: Vec<(i64, i64)> = Vec::new();
    let mut run_start: Option<usize> = None;
    let mut run_end: usize = 0;

    let mut step = |pos: usize, col: &Col| {
        let callable = matches!(
            determine_state(b'A', col.depth, col.qc_pass, col.low_mapq, params),
            CallableState::Callable
        );
        if callable {
            if run_start.is_none() {
                run_start = Some(pos);
            }
            run_end = pos;
        } else if let Some(s) = run_start.take() {
            if (run_end - s + 1) as u32 >= min_run_len {
                intervals.push(((s - 1) as i64, run_end as i64));
            }
        }
    };

    for result in query {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() || flags.is_duplicate() || flags.is_qc_fail() {
            continue;
        }
        let start = match record.alignment_start() {
            Some(p) => p.get(),
            None => continue,
        };
        while emit_cursor < start {
            let col = window.pop_front().unwrap_or_default();
            step(emit_cursor, &col);
            emit_cursor += 1;
        }
        let mapq = record.mapping_quality().map_or(255u8, |m| m.get());
        let quals = record.quality_scores();
        let quals = quals.as_ref();
        let mut ref_pos = start;
        let mut query_off = 0usize;
        for op in record.cigar().as_ref() {
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        if pos >= emit_cursor {
                            let idx = pos - emit_cursor;
                            while window.len() <= idx {
                                window.push_back(Col::default());
                            }
                            let col = &mut window[idx];
                            let base_q = quals.get(query_off + i).copied().unwrap_or(0);
                            col.depth += 1;
                            if mapq >= params.min_mapping_quality && base_q >= params.min_base_quality {
                                col.qc_pass += 1;
                            }
                            if mapq <= params.max_low_mapq {
                                col.low_mapq += 1;
                            }
                        }
                    }
                    ref_pos += len;
                    query_off += len;
                }
                (true, false) => ref_pos += len,
                (false, true) => query_off += len,
                (false, false) => {}
            }
        }
    }
    while let Some(col) = window.pop_front() {
        step(emit_cursor, &col);
        emit_cursor += 1;
    }
    if let Some(s) = run_start {
        if (run_end - s + 1) as u32 >= min_run_len {
            intervals.push(((s - 1) as i64, run_end as i64));
        }
    }
    Ok(intervals)
}

/// GATK `CallableLoci` hierarchy — first failing condition wins. Mirrors the Scala
/// `determineCallableState`.
fn determine_state(
    ref_base: u8,
    depth: u32,
    qc_pass: u32,
    low_mapq: u32,
    params: &CallableLociParams,
) -> CallableState {
    if ref_base == b'N' || ref_base == b'n' {
        return CallableState::RefN;
    }
    if depth == 0 {
        return CallableState::NoCoverage;
    }
    let low_frac = low_mapq as f64 / depth as f64;
    if low_frac > params.max_fraction_low_mapq {
        return CallableState::PoorMappingQuality;
    }
    if qc_pass < params.min_depth {
        return CallableState::LowCoverage;
    }
    if params.max_depth.is_some_and(|m| qc_pass > m) {
        return CallableState::ExcessiveCoverage;
    }
    CallableState::Callable
}

fn pct_at_least(hist: &[u64], total: u64, min_depth: usize) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let at_least: u64 = hist[min_depth..].iter().sum();
    at_least as f64 / total as f64
}

fn median_from_hist(hist: &[u64], total: u64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let half = total / 2;
    let mut cumulative = 0u64;
    for (depth, &count) in hist.iter().enumerate() {
        cumulative += count;
        if cumulative >= half {
            return depth as f64;
        }
    }
    255.0
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name)
    }

    #[test]
    fn callable_intervals_cover_the_fixture_and_honor_run_length() {
        let bam = fixture("coverage.bam"); // chrM, 50 bp, well covered
        let params = CallableLociParams::default();

        // No run-length gate: some callable bases, all within the 50 bp contig, intervals
        // sorted and non-overlapping (BED 0-based half-open).
        let ivs = callable_intervals(&bam, "chrM", &params, 1, None).unwrap();
        assert!(!ivs.is_empty(), "expected callable intervals on the fixture");
        let callable_bases: i64 = ivs.iter().map(|(s, e)| e - s).sum();
        assert!((1..=50).contains(&callable_bases), "callable bases in range: {callable_bases}");
        for w in ivs.windows(2) {
            assert!(w[0].1 <= w[1].0, "intervals sorted & disjoint");
        }
        assert!(ivs.iter().all(|(s, e)| *s >= 0 && *e <= 50));

        // An impossibly long run-length gate drops everything (fixture is only 50 bp).
        let none = callable_intervals(&bam, "chrM", &params, 10_000, None).unwrap();
        assert!(none.is_empty(), "no run clears a 10 kb gate on a 50 bp contig");
    }

    #[test]
    fn estimate_molecule_lengths_on_fixture() {
        let (read_len, frag_len) = estimate_molecule_lengths(&fixture("coverage.bam"), None).unwrap();
        assert!(read_len > 0.0);
        assert!(frag_len >= read_len || frag_len == read_len); // fragment >= read, or == when unpaired
    }
}
