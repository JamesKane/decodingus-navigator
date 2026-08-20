//! The walker over the coverage and the callable loci. It is the Rust port of the Scala
//! `CoverageCallableWalker`, which itself replaces GATK `CollectWgsMetrics` and `CallableLoci`
//! over htsjdk. It makes one pass over a BAM or CRAM in coordinate order. It builds three things
//! for each contig on the main assembly. Those are a depth histogram, the callable state at each
//! position, and coverage statistics in the style of samtools.
//!
//! The parity target is the Scala walker, and not samtools. The difference that matters is the
//! mean base quality and the mean mapping quality. This walker takes the mean **over the base
//! observations**, as Σ quality / Σ depth. samtools takes the mean over the reads. The fixture
//! tests of this crate carry the expected values, computed by hand.
//!
//! On memory: a **pileup in a window that slides** finalizes each position once the read frontier
//! passes it. The peak memory is then the span of the reads that are open, and not the length of
//! the contig.
//!
//! One allocation has the size of a contig. It holds the reference bases of the contig that the
//! walker is on, one contig at a time, to find the N bases. Over the whole genome, HG002 peaks at
//! about 2 GB. A method with dense arrays for each contig would need about 84 GB.
//!
//! This walker needs a BAM in coordinate order. To stream the reference in windows as well is a
//! further improvement, and nobody has done it. The BED-interval output and the progress
//! callbacks of the Scala walker wait for later work.

use std::collections::{HashMap, HashSet, VecDeque};
use std::path::Path;

use noodles::core::Region;
use noodles::fasta;

use serde::{Deserialize, Serialize};

use crate::cancel::CancelToken;
use crate::contig;
use crate::error::AnalysisError;
use crate::reader;
use crate::readview::AlnRead;

/// The algorithm version, for the cache key of the coverage artifact. Raise it after any change
/// that alters the output. See plan §6, on the version of a cache.
///
/// At `coverage-2`, the parallel walker attributes a record by its reference id. Before that, the
/// coverage count was too low on a CRAM whose slices hold more than one reference, such as an
/// FTDNA Big Y.
///
/// The new version makes every cached `coverage-1` result stale, and that is correct. Nobody knows
/// which sort order and slice layout each of those files came from. Each alignment computes the
/// correct value again at its next analysis. That run also overwrites the stale
/// read-metrics and sex from the same fused walk.
pub const COVERAGE_VERSION: &str = "coverage-2";

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

/// The callable class at one position. These are the states of GATK `CallableLoci`. They form a
/// hierarchy, and the first condition that fails wins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CallableState {
    RefN,
    NoCoverage,
    PoorMappingQuality,
    LowCoverage,
    ExcessiveCoverage,
    Callable,
}

/// The count of bases in each callable state, for each contig.
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

/// The coverage statistics of each contig, in the style of samtools `coverage`. The mean is over
/// the base observations. See the module documentation.
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
    /// The depth histogram of this contig. Bin `d` holds the count of bases at depth `d`, and the
    /// index stops at 255. That is the same convention as the genome-wide
    /// [`CoverageResult::coverage_histogram`].
    ///
    /// It is empty for an import on the fast path, from a pipeline sidecar, because such an import
    /// has no histogram over the depths. `#[serde(default)]` lets a coverage blob that the cache
    /// holds from before this field still load. The histogram fills again at the next analysis, so
    /// `COVERAGE_VERSION` does not need a new value.
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
    /// The fraction of the observed bases that the walker removed for a low base quality, after
    /// those bases passed the MAPQ filter. This is PCT_EXC_BASEQ.
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

/// The accumulator at one position, inside the pileup window that slides.
#[derive(Clone, Default)]
struct Col {
    depth: u32,
    base_q_sum: u64,
    map_q_sum: u64,
    qc_pass: u32,
    low_mapq: u32,
    /// Bases excluded for low mapping quality (mutually exclusive with `exc_baseq`; MAPQ checked first).
    exc_mapq: u32,
    /// The count of bases that passed MAPQ, and that a low base quality then removed.
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

/// The finished output of one contig. It stays very small, and the code puts it into the result
/// in header order.
struct ContigOut {
    callable: ContigCallableMetrics,
    stats: ContigCoverageStats,
}

/// The mask of the N bases in the reference of one contig. It holds one bit for each base, and a
/// set bit means that the reference base is N. The callable class then needs about one eighth of
/// the memory that the raw reference bytes need.
///
/// That matters for the parallel walker. Each contig task that runs at the same time would else
/// hold its full reference for the whole pileup, and chr1 is about 248 MB.
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

/// The streaming state of the contig that the walker is on. The window that slides bounds the
/// memory, and that window is the span of the reads that are open. The length of the contig does
/// not bound it.
struct CurContig {
    name: String,
    length: usize,
    ref_n_mask: NMask,
    window: VecDeque<Col>,
    emit_cursor: usize, // 1-based position next to finalize; window front aligns here
    read_count: u64,
    cm: ContigCallableMetrics,
    /// The own depth histogram of this contig. It sits beside the global `Globals::hist`, so that
    /// the UI can show the histogram of each contig. Both walker paths finish through here.
    hist: Vec<u64>,
    covered: u64,
    total_base_obs: u64,
    base_q_total: u64,
    map_q_total: u64,
    sum_depth: u128,
    /// The count of bases that the walker dropped, because they mapped *before* the finalize
    /// frontier. That can happen only when the input is not strictly in coordinate order. See
    /// [`CurContig::add`]. [`CurContig::finish`] gives a warning about it. It stays 0 for the
    /// standard sorted layout.
    dropped_unsorted: u64,
}

impl CurContig {
    /// Build from the raw reference bytes of the contig. This keeps the compact N-mask alone, so
    /// the caller can drop the `ref_bases` buffer immediately after.
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
            dropped_unsorted: 0,
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

    /// Add one covered base at the 1-based `pos` to the window.
    ///
    /// The window holds only the positions at the finalize frontier or after it. A base *before*
    /// the frontier, where `pos < emit_cursor`, belongs to a column that the walker already
    /// emitted. The code drops that base, and it does not let `pos - emit_cursor` go below zero.
    ///
    /// This can happen only when the input is not strictly in coordinate order, which occurs on
    /// some vendor CRAMs, such as an FTDNA Big Y. The streaming pileup assumes sorted input at its
    /// root. So the few bases that are out of order count as dropped, which `finish` reports, and
    /// the walk does not crash. On the standard sorted layout, `pos >= emit_cursor` always holds,
    /// and this guard never fires.
    fn add(&mut self, pos: usize, base_q: u8, mapq: u8, params: &CallableLociParams) {
        if pos < self.emit_cursor {
            self.dropped_unsorted += 1;
            return;
        }
        let idx = pos - self.emit_cursor;
        while self.window.len() <= idx {
            self.window.push_back(Col::default());
        }
        let col = &mut self.window[idx];
        col.depth += 1;
        col.base_q_sum += base_q as u64;
        col.map_q_sum += mapq as u64;
        // A base goes to one exclusion reason and no more than one: MAPQ first, then base-Q.
        // Picard does the same.
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
        if self.dropped_unsorted > 0 {
            eprintln!(
                "warning: coverage on {}: dropped {} base(s) that mapped before the finalize \
                 frontier — the input is not strictly coordinate-sorted, so this contig's depth \
                 may be slightly under-counted",
                self.name, self.dropped_unsorted
            );
        }
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

/// One pass in coordinate order, with a pileup window that slides. A position finalizes once the
/// read frontier passes it. The peak memory is then the span of the open reads, and not the length
/// of the contig.
///
/// This needs a BAM or CRAM in coordinate order, which is the standard layout in genomics. It
/// needs the reference for two things: to find the N positions, and to decode a CRAM.
pub fn collect_coverage_callable(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<CoverageResult, AnalysisError> {
    collect_coverage_callable_with_progress(
        bam_path,
        reference_path,
        params,
        contig_allowlist,
        &mut |_, _| {},
        &CancelToken::none(),
    )
}

/// The same as [`collect_coverage_callable`], and it also reports
/// `progress(contigs_done, contigs_total)` as it finishes each contig that it tracks.
///
/// A pass over the whole genome takes minutes on a real WGS BAM. It can then drive a progress bar,
/// and it does not look stopped. This function needs a BAM in coordinate order, with the contigs
/// in order.
pub fn collect_coverage_callable_with_progress(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &mut dyn FnMut(usize, usize),
    cancel: &CancelToken,
) -> Result<CoverageResult, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, Some(reference_path))?;
    let mut state = CoverageState::new(&header, reference_path, *params, contig_allowlist)?;
    progress(0, state.total_tracked());
    let mut seen = 0u32;
    for result in reader.records_lazy(&header) {
        let record = result?;
        state.accept(&record, progress)?;
        seen += 1;
        if seen % 4096 == 0 {
            cancel.check()?;
        }
    }
    state.finish(progress)
}

/// The streaming accumulator for the coverage and the callable state. The separate walker and the
/// fused [`crate::unified`] walker share it, so both give the same numbers, to the last digit,
/// from one source of truth. Give it every record through [`CoverageState::accept`], which applies
/// the flag filter and the contig filter of the coverage pass inside. Then call
/// [`CoverageState::finish`].
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

    /// The contigs that this state walks. It is the denominator of the progress.
    pub(crate) fn total_tracked(&self) -> usize {
        self.total_tracked
    }

    /// Give it one record. It ignores a record that the coverage pass does not want. There are
    /// five of those: a record with no mapping, a secondary record, a supplementary record, a
    /// duplicate, and a qc-fail. It also ignores a record on a contig that the walker does not
    /// track. The fused walker can then give every record here, with no filter of its own. It
    /// calls `progress` when a contig finishes.
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
                // The histogram of this contig: every position sits at depth 0. That matches the
                // parallel path, where a contig that the walker never saw finalizes every position
                // at depth 0, through CurContig::finish.
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

/// The read filter of the coverage pass. It skips a read with no mapping, a secondary or
/// supplementary read, a duplicate, and a qc-fail. The sequential [`CoverageState`] and the
/// [`ContigCoverageAccum`] of one contig share it, so both pileups see the same set of reads.
fn coverage_passes_filter(record: &impl AlnRead) -> bool {
    let f = record.flags();
    !(f.is_unmapped() || f.is_secondary() || f.is_supplementary() || f.is_duplicate() || f.is_qc_fail())
}

/// Give one record, which already passed the filter, to the pileup window of a contig. It moves
/// the finalize frontier to the start of the read, and then adds each base that the read takes
/// from the reference. The sequential path and the path over one contig share it, so the tally at
/// each base is the same in both.
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

/// Build the genome-wide [`CoverageResult`]. It takes the merged sums of the histogram, the
/// territory and the depth, and the output of each contig, which is already in header order.
///
/// It is the single source of truth for the tail of the result. The sequential `finish` and the
/// parallel `merge_coverage_partials` both use it.
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
    // exclusion reasons (dup/unpaired/overlap/capped) are not tallied, so these do not sum to a total.
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

/// The coverage accumulator of one contig, for the parallel walker. It holds the pileup window of
/// that contig, and a local copy of the genome-wide accumulators. The merge adds those copies
/// across the contigs. The rayon fan-out builds one for each contig. Give it the records from the
/// region query of that contig.
pub(crate) struct ContigCoverageAccum {
    c: CurContig,
    g: Globals,
    params: CallableLociParams,
}

/// The finished coverage of one contig. It holds the output of that contig. It also holds the
/// share that the contig adds to the genome-wide sums of the histogram, the territory and the
/// depth.
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

    /// Give it one record. It applies the read filter of the coverage pass inside, so it ignores a
    /// record that the filter rejects. The caller can then give it every record from the region
    /// query of the contig.
    pub(crate) fn accept(&mut self, record: &impl AlnRead) {
        if coverage_passes_filter(record) {
            feed_into_contig(&mut self.c, record, &self.params, &mut self.g);
        }
    }

    /// Finalize the contig into a partial result that carries `ref_id`, so that a later step can
    /// put the partials back into header order. This sends out the window of the contig, and its
    /// tail that no read covered.
    ///
    /// A contig that saw no read still finalizes every position at depth 0. It counts the ref-N
    /// bases and the bases with no coverage exactly as the zero-coverage branch of the sequential
    /// walker does.
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

/// Merge the coverage partials of the contigs into the genome-wide [`CoverageResult`]. It adds up
/// the accumulators of the histogram, the territory and the depth. It then puts the output of each
/// contig into `ref_id` order, which is the header order. The result then matches that of the
/// sequential walker, to the last digit.
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

/// The mean read length and the mean fragment length, which is the template length. The code
/// samples the first 50k primary mapped reads, or about that many.
///
/// This is the proxy for the molecule length, which the callable run-length gate needs, and that
/// gate refers to itself. Long reads mean long molecules, and long molecules mean long callable
/// runs. The fragment length falls back to the read length when the templates have no pair, as in
/// single-end long-read data.
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
        // Take the fragment length from a read with a correct pair alone, and put an upper limit
        // on it. A chimeric pair, or a pair that is not correct, carries a very large |TLEN|.
        // That would move the mean far away, and the run-length gate with it. A single-end read
        // and a long read have no correct pair, so those fall back to the read length.
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

/// The CALLABLE intervals on one `contig`, in BED form, which is 0-based and half-open. The code
/// joins the intervals that touch, and it keeps a run only when that run holds `min_run_len` bases
/// or more.
///
/// It needs no reference. It classifies a position by the depth, the QC flags and the MAPQ,
/// through the GATK hierarchy. A region of N bases in the reference carries no read, and it comes
/// out as no-coverage. The window of open reads bounds the memory. This needs a BAM index.
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

/// The hierarchy of GATK `CallableLoci`. The first condition that fails wins. It has the same
/// shape as the Scala `determineCallableState`. `ref_is_n` is true when the reference base is N,
/// and such a base is not callable.
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

/// The median absolute deviation of the depth. It is the median of `|depth − median|` over the
/// depth histogram. The histogram stops the depth at index 255, so a deviation in that tail is a
/// lower bound. At a usual WGS coverage that makes no difference.
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
    fn add_tolerates_a_base_before_the_finalize_frontier() {
        // An input that is not in coordinate order can give a base whose position is *behind* the
        // finalize frontier of the window. The code must drop that base, and count it. It must not
        // let `pos - emit_cursor` go below zero. That was the bug: coverage.rs panicked with
        // `attempt to subtract with overflow`.
        let params = CallableLociParams::default();
        let mut g = Globals::new();
        let mut c = CurContig::new("chrT".into(), 100, vec![b'A'; 100]);

        // Advance the frontier to 50 (positions 1..50 finalize at depth 0).
        c.advance_to(50, &params, &mut g);
        // A base at pos 30 is behind the frontier → dropped, no panic.
        c.add(30, 30, 60, &params);
        assert_eq!(c.dropped_unsorted, 1, "the stale base is counted as dropped");
        // A base at/after the frontier is still recorded normally.
        c.add(60, 30, 60, &params);
        assert_eq!(c.dropped_unsorted, 1, "an in-order base is not dropped");

        // The finish runs to its end, it does not panic, and it gives the statistics of the
        // contig.
        let out = c.finish(&params, &mut g);
        assert_eq!(out.stats.end_pos, 100);
    }

    #[test]
    fn callable_intervals_cover_the_fixture_and_honor_run_length() {
        let bam = fixture("coverage.bam"); // chrM, 50 bp, well covered
        let params = CallableLociParams::default();

        // With no run-length gate, there are some callable bases, and all of them lie inside the
        // 50 bp contig. The intervals come in sorted order, and none of them overlaps another. The
        // form is BED, which is 0-based and half-open.
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

        // One fixture holds chrM alone. The other holds more than one contig, which are the
        // autosomes and chrX. The test thereby covers the sum invariant across more than one
        // contig.
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

            // The strong invariant. The histograms of the contigs build the genome-wide histogram
            // exactly. The genome-wide one is their sum, bin by bin.
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
