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
use crate::readview::AlnRead;

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
    /// This contig's depth histogram (bin `d` = bases at depth `d`, clamped at index 255 — same
    /// convention as the genome-wide [`CoverageResult::coverage_histogram`]). Empty for fast-path
    /// (pipeline-sidecar) imports, which have no per-depth histogram. `#[serde(default)]` keeps
    /// coverage blobs cached before this field was added loading (the histogram repopulates on the
    /// next analysis) — so `COVERAGE_VERSION` does not need a bump.
    #[serde(default)]
    pub histogram: Vec<u64>,
}

/// Combined coverage + callable result (replaces the Scala `CoverageCallableResult`'s
/// global-metrics + callable-summary + samtools-stats fields).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct CoverageResult {
    pub genome_territory: u64,
    pub mean_coverage: f64,
    pub median_coverage: f64,
    pub sd_coverage: f64,
    /// Median absolute deviation of depth (Picard MAD_COVERAGE). `#[serde(default)]` → coverage
    /// blobs cached before this field load with 0.0 and repopulate on the next analysis.
    #[serde(default)]
    pub mad_coverage: f64,
    /// Fraction of observed bases excluded for low mapping quality (Picard PCT_EXC_MAPQ).
    #[serde(default)]
    pub pct_exc_mapq: f64,
    /// Fraction of observed bases excluded for low base quality, MAPQ having passed (PCT_EXC_BASEQ).
    #[serde(default)]
    pub pct_exc_baseq: f64,
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
    /// Bases excluded for low mapping quality (mutually exclusive with `exc_baseq`; MAPQ checked first).
    exc_mapq: u32,
    /// Bases that passed MAPQ but were excluded for low base quality.
    exc_baseq: u32,
}

/// Global accumulators, folded as positions finalize.
struct Globals {
    hist: Vec<u64>,
    n: u64,
    sum_depth: u128,
    sum_sq: u128,
    /// Total observed bases excluded for low MAPQ / low base-Q (Picard PCT_EXC_{MAPQ,BASEQ}).
    sum_exc_mapq: u128,
    sum_exc_baseq: u128,
}

impl Globals {
    fn new() -> Self {
        Globals {
            hist: vec![0; HIST_LEN],
            n: 0,
            sum_depth: 0,
            sum_sq: 0,
            sum_exc_mapq: 0,
            sum_exc_baseq: 0,
        }
    }
}

/// Finished per-contig output (kept tiny; assembled into the result in header order).
struct ContigOut {
    callable: ContigCallableMetrics,
    stats: ContigCoverageStats,
}

/// Reference-N mask for one contig: one bit per base (set = the reference base is N), so
/// callable classification needs ~1/8 the memory of holding the raw reference bytes. That
/// matters for the parallel walker, where each concurrent contig task would otherwise pin its
/// full reference (chr1 ≈ 248 MB) for the whole pileup.
struct NMask {
    bits: Vec<u64>,
}

impl NMask {
    fn from_bases(bases: &[u8]) -> Self {
        let mut bits = vec![0u64; bases.len().div_ceil(64)];
        for (i, &b) in bases.iter().enumerate() {
            if b == b'N' || b == b'n' {
                bits[i >> 6] |= 1u64 << (i & 63);
            }
        }
        NMask { bits }
    }

    /// Whether the reference base at 0-based `idx` is N. Out-of-range reads as N, matching the
    /// old `ref_bases.get(..).unwrap_or(b'N')` defensiveness.
    fn is_n(&self, idx: usize) -> bool {
        self.bits.get(idx >> 6).map_or(true, |w| (w >> (idx & 63)) & 1 == 1)
    }
}

/// Streaming state for the contig currently being walked. Memory is bounded by the
/// sliding window (the span of currently-open reads), not the contig length.
struct CurContig {
    name: String,
    length: usize,
    ref_n_mask: NMask,
    window: VecDeque<Col>,
    emit_cursor: usize, // 1-based position next to finalize; window front aligns here
    read_count: u64,
    cm: ContigCallableMetrics,
    /// This contig's own depth histogram (kept alongside the global `Globals::hist` so the
    /// per-contig histogram can be surfaced; both walker paths finalize through here).
    hist: Vec<u64>,
    covered: u64,
    total_base_obs: u64,
    base_q_total: u64,
    map_q_total: u64,
    sum_depth: u128,
}

impl CurContig {
    /// Builds from the contig's raw reference bytes, retaining only the compact N-mask (the
    /// `ref_bases` buffer can be dropped by the caller right after).
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
            ref_n_mask: NMask::from_bases(&ref_bases),
            window: VecDeque::new(),
            emit_cursor: 1,
            read_count: 0,
            cm,
            hist: vec![0; HIST_LEN],
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
        self.hist[clamped] += 1;
        g.n += 1;
        g.sum_depth += depth as u128;
        g.sum_sq += (depth as u128) * (depth as u128);
        g.sum_exc_mapq += col.exc_mapq as u128;
        g.sum_exc_baseq += col.exc_baseq as u128;
        self.sum_depth += depth as u128;
        if depth > 0 {
            self.covered += 1;
            self.total_base_obs += depth as u64;
            self.base_q_total += col.base_q_sum;
            self.map_q_total += col.map_q_sum;
        }
        let ref_is_n = self.ref_n_mask.is_n(pos - 1);
        match determine_state(ref_is_n, depth, col.qc_pass, col.low_mapq, params) {
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
        // Mutually-exclusive exclusion attribution (MAPQ first, then base-Q), mirroring Picard.
        if mapq < params.min_mapping_quality {
            col.exc_mapq += 1;
        } else if base_q < params.min_base_quality {
            col.exc_baseq += 1;
        } else {
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
            coverage: if self.length == 0 {
                0.0
            } else {
                self.covered as f64 / length * 100.0
            },
            mean_depth: if self.length == 0 {
                0.0
            } else {
                self.sum_depth as f64 / length
            },
            mean_base_q: if self.total_base_obs == 0 {
                0.0
            } else {
                self.base_q_total as f64 / self.total_base_obs as f64
            },
            mean_map_q: if self.total_base_obs == 0 {
                0.0
            } else {
                self.map_q_total as f64 / self.total_base_obs as f64
            },
            histogram: std::mem::take(&mut self.hist),
        };
        ContigOut {
            callable: self.cm,
            stats,
        }
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
    let mut state = CoverageState::new(&header, reference_path, *params, contig_allowlist)?;
    progress(0, state.total_tracked());
    for result in reader.records_lazy(&header) {
        let record = result?;
        state.accept(&record, progress)?;
    }
    state.finish(progress)
}

/// Streaming coverage + callable accumulator shared by the standalone walker and the fused
/// [`crate::unified`] walker, so both produce byte-identical numbers from one source of
/// truth. Feed it every record via [`CoverageState::accept`] (it applies coverage's own
/// flag/contig filtering internally), then call [`CoverageState::finish`].
pub(crate) struct CoverageState {
    /// ref_id -> (name, length) for tracked (main-assembly, allowlisted) contigs; `None` elsewhere.
    tracked: Vec<Option<(String, usize)>>,
    fasta_reader: fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,
    reference_path: std::path::PathBuf,
    params: CallableLociParams,
    g: Globals,
    finished: HashMap<usize, ContigOut>,
    cur: Option<(usize, CurContig)>,
    total_tracked: usize,
    contigs_done: usize,
}

impl CoverageState {
    pub(crate) fn new(
        header: &noodles::sam::Header,
        reference_path: &Path,
        params: CallableLociParams,
        contig_allowlist: Option<&HashSet<String>>,
    ) -> Result<Self, AnalysisError> {
        let tracked: Vec<Option<(String, usize)>> = header
            .reference_sequences()
            .iter()
            .map(|(name_bytes, map)| {
                let name = String::from_utf8_lossy(name_bytes.as_ref()).into_owned();
                let keep = contig::is_main_assembly(&name) && contig_allowlist.map_or(true, |set| set.contains(&name));
                keep.then(|| (name, map.length().get()))
            })
            .collect();
        let total_tracked = tracked.iter().filter(|o| o.is_some()).count();
        let fasta_reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(reference_path)
            .map_err(|e| AnalysisError::io(reference_path, e))?;
        Ok(CoverageState {
            tracked,
            fasta_reader,
            reference_path: reference_path.to_path_buf(),
            params,
            g: Globals::new(),
            finished: HashMap::new(),
            cur: None,
            total_tracked,
            contigs_done: 0,
        })
    }

    /// Contigs this state will walk — the progress denominator.
    pub(crate) fn total_tracked(&self) -> usize {
        self.total_tracked
    }

    /// Feed one record. Records that coverage doesn't care about (unmapped/secondary/
    /// supplementary/duplicate/qc-fail, or off a tracked contig) are ignored, so the fused
    /// walker can hand every record here unfiltered. Fires `progress` on contig finalization.
    pub(crate) fn accept(
        &mut self,
        record: &impl AlnRead,
        progress: &mut dyn FnMut(usize, usize),
    ) -> Result<(), AnalysisError> {
        if !coverage_passes_filter(record) {
            return Ok(());
        }
        let ref_id = match record.reference_sequence_id() {
            Some(r) => r,
            None => return Ok(()),
        };
        if !matches!(self.tracked.get(ref_id), Some(Some(_))) {
            return Ok(());
        }

        // Contig transition: finalize the previous contig, load the new contig's ref.
        if self.cur.as_ref().map(|(id, _)| *id) != Some(ref_id) {
            if let Some((id, c)) = self.cur.take() {
                self.finished.insert(id, c.finish(&self.params, &mut self.g));
                self.contigs_done += 1;
                progress(self.contigs_done, self.total_tracked);
            }
            let (name, length) = {
                let t = self.tracked[ref_id].as_ref().unwrap();
                (t.0.clone(), t.1)
            };
            let region: Region = name
                .parse()
                .map_err(|_| AnalysisError::Message(format!("bad region for contig {name}")))?;
            let rec = self
                .fasta_reader
                .query(&region)
                .map_err(|e| AnalysisError::io(&self.reference_path, e))?;
            let ref_bases = rec.sequence().as_ref().to_vec();
            self.cur = Some((ref_id, CurContig::new(name, length, ref_bases)));
        }

        let (_, c) = self.cur.as_mut().unwrap();
        feed_into_contig(c, record, &self.params, &mut self.g);
        Ok(())
    }

    /// Finalize the last contig, then assemble the result in header order (tracked contigs
    /// with no reads are zero-coverage, ref-N counted from the reference).
    pub(crate) fn finish(mut self, progress: &mut dyn FnMut(usize, usize)) -> Result<CoverageResult, AnalysisError> {
        if let Some((id, c)) = self.cur.take() {
            self.finished.insert(id, c.finish(&self.params, &mut self.g));
            self.contigs_done += 1;
            progress(self.contigs_done, self.total_tracked);
        }

        let mut contig_callable = Vec::new();
        let mut contig_stats = Vec::new();
        for ref_id in 0..self.tracked.len() {
            let Some((name, length)) = self.tracked[ref_id].clone() else {
                continue;
            };
            if let Some(out) = self.finished.remove(&ref_id) {
                contig_callable.push(out.callable);
                contig_stats.push(out.stats);
            } else {
                // No reads: every position is depth 0 (RefN where the reference is N).
                let region: Region = name
                    .parse()
                    .map_err(|_| AnalysisError::Message(format!("bad region for contig {name}")))?;
                let rec = self
                    .fasta_reader
                    .query(&region)
                    .map_err(|e| AnalysisError::io(&self.reference_path, e))?;
                let ref_bases = rec.sequence();
                let ref_bases = ref_bases.as_ref();
                let mut ref_n: u64 = 0;
                for idx in 0..length {
                    let b = ref_bases.get(idx).copied().unwrap_or(b'N');
                    if b == b'N' || b == b'n' {
                        ref_n += 1;
                    }
                }
                self.g.hist[0] += length as u64;
                self.g.n += length as u64;
                // Per-contig histogram: every position at depth 0 (matches the parallel path,
                // where an unseen contig finalizes all positions at depth 0 via CurContig::finish).
                let mut hist = vec![0u64; HIST_LEN];
                hist[0] = length as u64;
                contig_callable.push(ContigCallableMetrics {
                    contig: name.clone(),
                    ref_n,
                    callable: 0,
                    no_coverage: length as u64 - ref_n,
                    low_coverage: 0,
                    excessive_coverage: 0,
                    poor_mapping_quality: 0,
                });
                contig_stats.push(ContigCoverageStats {
                    contig: name.clone(),
                    start_pos: 1,
                    end_pos: length as u64,
                    num_reads: 0,
                    cov_bases: 0,
                    coverage: 0.0,
                    mean_depth: 0.0,
                    mean_base_q: 0.0,
                    mean_map_q: 0.0,
                    histogram: hist,
                });
            }
        }

        Ok(assemble_coverage_result(
            self.g.hist,
            self.g.n,
            self.g.sum_depth,
            self.g.sum_sq,
            self.g.sum_exc_mapq,
            self.g.sum_exc_baseq,
            contig_callable,
            contig_stats,
        ))
    }
}

/// Coverage's read filter — skip unmapped / secondary / supplementary / duplicate / qc-fail.
/// Shared by the sequential [`CoverageState`] and the per-contig [`ContigCoverageAccum`] so
/// both pileups see the identical read set.
fn coverage_passes_filter(record: &impl AlnRead) -> bool {
    let f = record.flags();
    !(f.is_unmapped() || f.is_secondary() || f.is_supplementary() || f.is_duplicate() || f.is_qc_fail())
}

/// Feed one (already filter-passing) record into a contig's sliding-window pileup: advance
/// the finalize frontier to the read's start, then add each reference-consuming base. Shared
/// by the sequential and per-contig coverage paths so the per-base accounting is identical.
fn feed_into_contig(c: &mut CurContig, record: &impl AlnRead, params: &CallableLociParams, g: &mut Globals) {
    let start = match record.alignment_start() {
        Some(p) => p,
        None => return,
    };
    c.advance_to(start, params, g);
    c.read_count += 1;

    let mapq = record.mapping_quality().unwrap_or(255u8);
    record.pileup_with(|quals, ops| {
        let mut ref_pos = start; // 1-based
        let mut query_off = 0usize;
        for (kind, len) in ops {
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
    });
}

/// Assemble the genome-wide [`CoverageResult`] from merged histogram/territory/depth sums and
/// per-contig outputs (already in header order). Single source of truth for the result tail,
/// used by both the sequential `finish` and the parallel `merge_coverage_partials`.
#[allow(clippy::too_many_arguments)] // one accumulator per metric — a coverage roll-up, not a refactor target
fn assemble_coverage_result(
    hist: Vec<u64>,
    n: u64,
    sum_depth: u128,
    sum_sq: u128,
    sum_exc_mapq: u128,
    sum_exc_baseq: u128,
    contig_callable: Vec<ContigCallableMetrics>,
    contig_coverage_stats: Vec<ContigCoverageStats>,
) -> CoverageResult {
    let mean = if n == 0 { 0.0 } else { sum_depth as f64 / n as f64 };
    let sd = if n < 2 {
        0.0
    } else {
        (sum_sq as f64 / n as f64 - mean * mean).max(0.0).sqrt()
    };
    let median = median_from_hist(&hist, n);
    // Exclusion fractions over total observed bases (Picard PCT_EXC_{MAPQ,BASEQ}). `sum_depth`
    // already counts every observed base (excluded ones included), so it is the denominator. Other
    // exclusion reasons (dup/unpaired/overlap/capped) aren't tallied, so these don't sum to a total.
    let (pct_exc_mapq, pct_exc_baseq) = if sum_depth == 0 {
        (0.0, 0.0)
    } else {
        (
            sum_exc_mapq as f64 / sum_depth as f64,
            sum_exc_baseq as f64 / sum_depth as f64,
        )
    };
    let callable_bases = contig_callable.iter().map(|c| c.callable).sum();
    CoverageResult {
        genome_territory: n,
        mean_coverage: mean,
        median_coverage: median,
        sd_coverage: sd,
        mad_coverage: mad_from_hist(&hist, n, median),
        pct_exc_mapq,
        pct_exc_baseq,
        pct_1x: pct_at_least(&hist, n, 1),
        pct_5x: pct_at_least(&hist, n, 5),
        pct_10x: pct_at_least(&hist, n, 10),
        pct_15x: pct_at_least(&hist, n, 15),
        pct_20x: pct_at_least(&hist, n, 20),
        pct_25x: pct_at_least(&hist, n, 25),
        pct_30x: pct_at_least(&hist, n, 30),
        pct_40x: pct_at_least(&hist, n, 40),
        pct_50x: pct_at_least(&hist, n, 50),
        coverage_histogram: hist,
        callable_bases,
        contig_callable,
        contig_coverage_stats,
    }
}

/// Per-contig coverage accumulator for the parallel walker — one contig's sliding-window
/// pileup with a local copy of the genome-wide accumulators (summed across contigs at merge
/// time). Built per contig in the rayon fan-out; feed it that contig's region-query records.
pub(crate) struct ContigCoverageAccum {
    c: CurContig,
    g: Globals,
    params: CallableLociParams,
}

/// One contig's finished coverage contribution: its per-contig output plus this contig's
/// share of the genome-wide histogram / territory / depth sums.
pub(crate) struct ContigCoveragePartial {
    ref_id: usize,
    callable: ContigCallableMetrics,
    stats: ContigCoverageStats,
    hist: Vec<u64>,
    n: u64,
    sum_depth: u128,
    sum_sq: u128,
    sum_exc_mapq: u128,
    sum_exc_baseq: u128,
}

impl ContigCoverageAccum {
    pub(crate) fn new(name: String, length: usize, ref_bases: Vec<u8>, params: CallableLociParams) -> Self {
        ContigCoverageAccum {
            c: CurContig::new(name, length, ref_bases),
            g: Globals::new(),
            params,
        }
    }

    /// Feed one record; the coverage read filter is applied internally, so off-filter records
    /// are ignored and the caller can pass every record from the contig's region query.
    pub(crate) fn accept(&mut self, record: &impl AlnRead) {
        if coverage_passes_filter(record) {
            feed_into_contig(&mut self.c, record, &self.params, &mut self.g);
        }
    }

    /// Finalize the contig (flushing its window + uncovered tail) into a partial tagged with
    /// `ref_id` for header-order reassembly. A contig that saw no reads still finalizes every
    /// position at depth 0 — counting ref-N and no-coverage exactly as the sequential walker's
    /// zero-coverage branch does.
    pub(crate) fn finish(mut self, ref_id: usize) -> ContigCoveragePartial {
        let out = self.c.finish(&self.params, &mut self.g);
        ContigCoveragePartial {
            ref_id,
            callable: out.callable,
            stats: out.stats,
            hist: self.g.hist,
            n: self.g.n,
            sum_depth: self.g.sum_depth,
            sum_sq: self.g.sum_sq,
            sum_exc_mapq: self.g.sum_exc_mapq,
            sum_exc_baseq: self.g.sum_exc_baseq,
        }
    }
}

/// Merge per-contig coverage partials into the genome-wide [`CoverageResult`]: sum the
/// histogram/territory/depth accumulators and order the per-contig outputs by `ref_id` (header
/// order), so the result is byte-identical to the sequential walker's.
pub(crate) fn merge_coverage_partials(mut partials: Vec<ContigCoveragePartial>) -> CoverageResult {
    partials.sort_by_key(|p| p.ref_id);
    let mut hist = vec![0u64; HIST_LEN];
    let (mut n, mut sum_depth, mut sum_sq) = (0u64, 0u128, 0u128);
    let (mut sum_exc_mapq, mut sum_exc_baseq) = (0u128, 0u128);
    let mut contig_callable = Vec::with_capacity(partials.len());
    let mut contig_stats = Vec::with_capacity(partials.len());
    for p in partials {
        for (i, v) in p.hist.iter().enumerate() {
            hist[i] += v;
        }
        n += p.n;
        sum_depth += p.sum_depth;
        sum_sq += p.sum_sq;
        sum_exc_mapq += p.sum_exc_mapq;
        sum_exc_baseq += p.sum_exc_baseq;
        contig_callable.push(p.callable);
        contig_stats.push(p.stats);
    }
    assemble_coverage_result(
        hist,
        n,
        sum_depth,
        sum_sq,
        sum_exc_mapq,
        sum_exc_baseq,
        contig_callable,
        contig_stats,
    )
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
    let frag_len = if frag_n > 0 {
        frag_sum as f64 / frag_n as f64
    } else {
        read_len
    };
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
            determine_state(false, col.depth, col.qc_pass, col.low_mapq, params),
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
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
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
/// `determineCallableState`. `ref_is_n` is whether the reference base is N (non-callable).
fn determine_state(
    ref_is_n: bool,
    depth: u32,
    qc_pass: u32,
    low_mapq: u32,
    params: &CallableLociParams,
) -> CallableState {
    if ref_is_n {
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

/// Median absolute deviation of depth: the median of `|depth − median|` over the depth histogram.
/// Depths are clamped at index 255 in the histogram, so deviations in that tail are lower bounds
/// (negligible for typical WGS coverage).
fn mad_from_hist(hist: &[u64], total: u64, median: f64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    // Histogram of absolute deviations (same 0..=255 domain bound as depth).
    let mut dev: Vec<u64> = vec![0; hist.len()];
    for (depth, &count) in hist.iter().enumerate() {
        let d = (depth as f64 - median).abs().round() as usize;
        dev[d.min(hist.len() - 1)] += count;
    }
    median_from_hist(&dev, total)
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
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name)
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
        assert!(
            (1..=50).contains(&callable_bases),
            "callable bases in range: {callable_bases}"
        );
        for w in ivs.windows(2) {
            assert!(w[0].1 <= w[1].0, "intervals sorted & disjoint");
        }
        assert!(ivs.iter().all(|(s, e)| *s >= 0 && *e <= 50));

        // An impossibly long run-length gate drops everything (fixture is only 50 bp).
        let none = callable_intervals(&bam, "chrM", &params, 10_000, None).unwrap();
        assert!(none.is_empty(), "no run clears a 10 kb gate on a 50 bp contig");
    }

    #[test]
    fn mad_from_histogram() {
        // Depths {0,2,10,12}, one position each. median_from_hist uses the lower-median (cumulative
        // ≥ total/2) convention → median 2. |0-2|,|2-2|,|10-2|,|12-2| = {2,0,8,10}; sorted {0,2,8,10},
        // lower-median → 2.
        let mut hist = vec![0u64; HIST_LEN];
        for d in [0usize, 2, 10, 12] {
            hist[d] += 1;
        }
        let median = median_from_hist(&hist, 4);
        assert_eq!(median, 2.0);
        assert_eq!(mad_from_hist(&hist, 4, median), 2.0);
        // Constant depth → MAD 0.
        let mut flat = vec![0u64; HIST_LEN];
        flat[30] = 100;
        assert_eq!(mad_from_hist(&flat, 100, median_from_hist(&flat, 100)), 0.0);
    }

    #[test]
    fn per_contig_histograms_sum_to_genome_wide() {
        let params = CallableLociParams::default();

        // chrM-only and a multi-contig (autosomes + chrX) fixture, so the sum invariant is
        // exercised across more than one contig.
        for (bam_name, ref_name) in [("coverage.bam", "ref.fa"), ("sex.bam", "sexref.fa")] {
            let cov = collect_coverage_callable(&fixture(bam_name), &fixture(ref_name), &params, None).unwrap();
            assert!(
                !cov.contig_coverage_stats.is_empty(),
                "{bam_name}: expected tracked contigs"
            );

            let width = cov.coverage_histogram.len();
            let mut summed = vec![0u64; width];
            for s in &cov.contig_coverage_stats {
                // Every contig carries a full-width histogram...
                assert_eq!(s.histogram.len(), width, "{bam_name}/{}: histogram width", s.contig);
                // ...and every finalized position lands in exactly one depth bin, so the bins
                // total the contig length.
                let contig_total: u64 = s.histogram.iter().sum();
                assert_eq!(
                    contig_total, s.end_pos,
                    "{bam_name}/{}: histogram bins should total the contig length",
                    s.contig
                );
                for (acc, v) in summed.iter_mut().zip(&s.histogram) {
                    *acc += v;
                }
            }

            // The strong invariant: per-contig histograms reconstruct the genome-wide histogram
            // exactly (the genome-wide one is just their bin-wise sum).
            assert_eq!(
                summed, cov.coverage_histogram,
                "{bam_name}: per-contig histograms must sum bin-for-bin to the genome-wide histogram"
            );
        }
    }

    #[test]
    fn estimate_molecule_lengths_on_fixture() {
        let (read_len, frag_len) = estimate_molecule_lengths(&fixture("coverage.bam"), None).unwrap();
        assert!(read_len > 0.0);
        assert!(frag_len >= read_len || frag_len == read_len); // fragment >= read, or == when unpaired
    }
}
