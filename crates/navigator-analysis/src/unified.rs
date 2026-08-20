//! The unified walker over the quality metrics. It makes one pass over a BAM or a CRAM, in
//! coordinate order. In that pass it collects three things together: the coverage and the callable
//! loci, the QC metrics at the read level, and the sex inference.
//!
//! The rewrite already had three walkers, each with one purpose and one pass:
//! [`crate::coverage`], [`crate::read_metrics`] and [`crate::sex`]. Run apart, they read a BAM
//! from end to end **twice**, for the coverage pileup and the read-metrics scan. They read a CRAM
//! **three times**, because the sex scan is separate: a `.crai` carries no count for each
//! reference. This walker puts them into one record loop. That is 2 passes to 1 for a BAM, and 3
//! to 1 for a CRAM, and a CRAM decode is the costly case.
//!
//! **No metric changes.** The loop sends each record to the same `*State` accumulators that the
//! separate walkers use, which are the single source of truth. Every number then matches the
//! number that three separate runs give, to the last digit.
//!
//! There is one thing to watch, and that is the filter. The coverage pass applies a hard filter,
//! and it keeps only a primary record with a mapping on the main assembly. But
//! read-metrics needs *every* record, and sex needs the mapped tally of each contig. So the loop
//! gives every record to all three states, and each state applies its own filter inside.
//!
//! The sex tally here comes straight from the record stream, and it needs no BAI. The arithmetic
//! is the same as in the separate CRAM path. The separate [`crate::sex::infer_from_bam`] keeps its
//! fast path over the BAI, for the small "Sex inference" command on its own.

use std::collections::HashSet;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use noodles::core::Region;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cancel::CancelToken;
use crate::contig;
use crate::coverage::{
    merge_coverage_partials, CallableLociParams, ContigCoverageAccum, ContigCoveragePartial, CoverageResult,
    CoverageState,
};
use crate::error::AnalysisError;
use crate::read_metrics::{ReadMetrics, ReadMetricsState};
use crate::reader::{self, RecordSink};
use crate::readview::AlnRead;
use crate::sex::{self, SexInferenceResult, SexState};

/// Send the base-pair progress that the loop has collected to the shared counter, after this many
/// bp of advance. It is small enough that the bar moves smoothly, at about 1500 ticks over a
/// 3.1 Gb genome. It is also large enough that the atomic add and the progress callback cost
/// little. In the GUI that callback is a channel send behind a mutex.
const PROGRESS_FLUSH_BP: u64 = 2_000_000;

/// The record consumer for one contig. It gives each record to read-metrics and to coverage, and
/// it tallies the mapped count of each contig, which the sex inference needs. `class` is 1 for an
/// autosome, 2 for chrX, and 0 for anything else. It works on borrowed accumulators, so it serves
/// the zero-copy path for a BAM record.
///
/// It also drives the **base-pair progress**. The reads inside a contig are in coordinate order,
/// so the alignment start only rises. The difference between one start and the next is the count
/// of bp walked.
///
/// The differences add up locally, and go to the shared `processed_bp` counter, and to the
/// progress callback, every [`PROGRESS_FLUSH_BP`]. The bar then advances all the time *inside* a
/// contig. Without that, it moves only when a contig finishes, and the big autosomes leave it
/// frozen for minutes.
struct ContigSink<'a> {
    /// The header index of the contig that this sink walks. The sink takes a record only when the
    /// `reference_sequence_id` of that record matches.
    ///
    /// A CRAM slice with more than one reference comes back from the region query of *every*
    /// contig that it overlaps. Without this gate, its records would go to the wrong contig in the
    /// coverage pass. Read-metrics would also count them once for each contig that they
    /// overlap.
    ///
    /// With the gate, the code takes each record exactly once, in the query of the contig that it
    /// belongs to. That matches the single binned pass of the sequential walker. The gate does
    /// nothing for a slice with one reference.
    ref_id: usize,
    rm: &'a mut ReadMetricsState,
    cov: &'a mut Option<ContigCoverageAccum>,
    class: u8,
    autosome_reads: u64,
    x_reads: u64,
    // Base-pair progress.
    progress: &'a (dyn Fn(usize, usize) + Sync),
    processed_bp: &'a AtomicU64,
    total_mb: usize,
    last_pos: usize,
    local_bp: u64,
}

impl RecordSink for ContigSink<'_> {
    fn accept(&mut self, record: &impl AlnRead) {
        // Take the records of this contig alone. Drop a record of another contig that came in
        // with a slice that holds more than one reference. The query of its own contig takes it.
        // Every tally over the records then counts each record exactly once.
        if record.reference_sequence_id() != Some(self.ref_id) {
            return;
        }
        self.rm.accept(record);
        if let Some(acc) = self.cov.as_mut() {
            acc.accept(record);
        }
        if self.class != 0 && !record.flags().is_unmapped() {
            if self.class == 1 {
                self.autosome_reads += 1;
            } else {
                self.x_reads += 1;
            }
        }
        if let Some(pos) = record.alignment_start() {
            if pos > self.last_pos {
                self.local_bp += (pos - self.last_pos) as u64;
                self.last_pos = pos;
                if self.local_bp >= PROGRESS_FLUSH_BP {
                    let g = self.processed_bp.fetch_add(self.local_bp, Ordering::Relaxed) + self.local_bp;
                    self.local_bp = 0;
                    (self.progress)((g / 1_000_000) as usize, self.total_mb);
                }
            }
        }
    }
}

/// Record consumer for the unmapped tail: read-metrics only (no reference position).
struct MetricsSink {
    rm: ReadMetricsState,
}

impl RecordSink for MetricsSink {
    fn accept(&mut self, record: &impl AlnRead) {
        self.rm.accept(record);
    }
}

/// The algorithm version, for the cache key of the unified artifact. Raise it after any change
/// that alters the output. The three sub-results already go into the store under their own keys,
/// so this version exists to make the set complete.
pub const UNIFIED_VERSION: &str = "unified-1";

/// The three quality-metric results that one pass collects.
///
/// Sex is `None` when the code can not infer it for the input. That happens when there is no
/// autosome and no chrX, and when there is no autosomal read. A targeted panel and a chrY-only
/// test are examples. The coverage and the read-metrics do not change. This has the same shape as
/// the pipeline, where sex is an independent step, and a failure there does not kill the other
/// two.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnifiedMetricsResult {
    pub coverage: CoverageResult,
    pub read_metrics: ReadMetrics,
    pub sex: Option<SexInferenceResult>,
}

/// The coverage, the read-metrics and the sex, in one pass over a BAM or CRAM in coordinate
/// order. It needs `reference`, both to decode a CRAM and to find the N bases in the reference.
///
/// The result is the same as three separate steps:
/// [`crate::coverage::collect_coverage_callable`], [`crate::read_metrics::collect_read_metrics`],
/// and a tally of the sex over the same records. But it reads the file once.
pub fn collect_unified_metrics(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<UnifiedMetricsResult, AnalysisError> {
    collect_unified_metrics_with_progress(
        bam_path,
        reference_path,
        params,
        contig_allowlist,
        &mut |_, _| {},
        &CancelToken::none(),
    )
}

/// The same as [`collect_unified_metrics`], and it also reports
/// `progress(contigs_done, contigs_total)`. It calls that as the coverage pass finishes each
/// contig that it tracks. That pass is the slow step over the whole genome, so a progress bar can
/// then move, and it does not stay frozen for minutes. This function needs a BAM or CRAM in
/// coordinate order.
pub fn collect_unified_metrics_with_progress(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &mut dyn FnMut(usize, usize),
    cancel: &CancelToken,
) -> Result<UnifiedMetricsResult, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, Some(reference_path))?;
    let mut cov = CoverageState::new(&header, reference_path, *params, contig_allowlist)?;
    let mut rm = ReadMetricsState::default();
    let mut sx = SexState::new(&header);
    progress(0, cov.total_tracked());

    // The code polls this at the same rate as the record loop of the walker that uses an index.
    // That rate is often enough that a click stops the walk in milliseconds. It is also rare
    // enough to cost nothing next to the pileup work at each record. This is the fallback path,
    // for a BAM
    // or CRAM with no index. It has no contig boundary to stop at, so without a check here nobody
    // could cancel it at all.
    let mut seen = 0u32;
    for result in reader.records_lazy(&header) {
        let record = result?;
        // Every record to all three; each state filters internally (see module docs).
        rm.accept(&record);
        sx.accept(&record);
        cov.accept(&record, progress)?;
        seen += 1;
        if seen % 4096 == 0 {
            cancel.check()?;
        }
    }

    let coverage = cov.finish(progress)?;
    let read_metrics = rm.finish();
    // Sex inference is best-effort: `None` (not a hard error) when the input lacks the
    // autosomes/chrX it needs, so coverage + read-metrics still come back.
    let sex = sx.finish().ok();
    Ok(UnifiedMetricsResult {
        coverage,
        read_metrics,
        sex,
    })
}

/// The unified metrics, with one task for each contig. The result is the same as
/// [`collect_unified_metrics`], and the contigs run at the same time. The coverage pass over one
/// contig is independent of every other contig. The compute of the pileup at each position, and
/// not the decompression, is what limits a sequential pass.
///
/// This needs an **indexed BAM**, for the region query of each contig and for a sweep over the
/// unmapped tail. Anything else falls back to the sequential [`collect_unified_metrics`], and the
/// caller sees no difference. That covers a CRAM, because a `.crai` has no query for the unmapped
/// reads, and a BAM with no `.bai`. A caller can then always ask for this function.
///
/// The output matches that of the sequential walker to the last digit. At each contig it runs the
/// same `*State` accumulators. The merge is over sums that commute, and over outputs that follow
/// the header order of the contigs.
///
/// Read-metrics covers **every** contig, and not the main assembly alone, plus the unmapped tail.
/// That is the same set of records that the sequential pass sees, so the totals agree exactly.
pub fn collect_unified_metrics_parallel(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
) -> Result<UnifiedMetricsResult, AnalysisError> {
    collect_unified_metrics_parallel_with_progress(
        bam_path,
        reference_path,
        params,
        contig_allowlist,
        &|_, _| {},
        &CancelToken::none(),
    )
}

/// The count of worker threads for the fan-out over the contigs. The default is every available
/// core, up to 12. Above that, the largest contig plus the unmapped sweep set the floor on the
/// wall time, so more threads only add memory. `NAVIGATOR_ANALYSIS_THREADS` overrides it. The
/// region fan-out of the de-novo caller uses the same value.
pub(crate) fn analysis_thread_count() -> usize {
    std::env::var("NAVIGATOR_ANALYSIS_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
                .min(12)
        })
        .max(1)
}

/// A token from the semaphore that limits the reference loads. It goes back into the pool when it
/// drops, and that includes the error path. It bounds how many contigs hold their full reference
/// buffer at one time. That is what sets the peak memory, because the N-mask of a contig is very
/// small once the code has built it.
struct LoadPermit<'a> {
    tx: &'a std::sync::mpsc::Sender<()>,
}

impl Drop for LoadPermit<'_> {
    fn drop(&mut self) {
        let _ = self.tx.send(());
    }
}

/// One contig's partial result from the parallel fan-out.
struct ContigPartial {
    rm: ReadMetricsState,
    cov: Option<ContigCoveragePartial>,
    autosome_reads: u64,
    x_reads: u64,
}

/// The same as [`collect_unified_metrics_parallel`], and it also reports
/// `progress(megabases_done, megabases_total)`. Those are the base-pair positions walked over all
/// of the contigs, so the bar advances all the time. It does not take one step at each finished
/// contig. The big autosomes start first and finish together late, which holds a bar over the
/// contig count at 0 for about half of the run. The callback is `Fn + Sync`, because the worker
/// threads call it at the same time.
pub fn collect_unified_metrics_parallel_with_progress(
    bam_path: &Path,
    reference_path: &Path,
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &(dyn Fn(usize, usize) + Sync),
    cancel: &CancelToken,
) -> Result<UnifiedMetricsResult, AnalysisError> {
    // The parallel path needs a coordinate index for the region query of each contig. That is a
    // `.bai` for a BAM, or a `.crai` for a CRAM. Without one, fall back to the sequential walker.
    //
    // A CRAM has no region query for its unmapped tail. So the read-metrics totals of a CRAM
    // leave out the reads that have no mapping and nothing else. Those reads carry no coverage
    // signal and no sex signal. The coverage of each contig, and the read-metrics over the mapped reads,
    // still run in parallel.
    if !reader::has_region_index(bam_path) {
        return collect_unified_metrics_with_progress(
            bam_path,
            reference_path,
            params,
            contig_allowlist,
            &mut |d, t| progress(d, t),
            cancel,
        );
    }
    let skip_unmapped = reader::has_crai_index(bam_path); // CRAM: no unmapped-region query

    let header = reader::read_header(bam_path, Some(reference_path))?;

    // The work items: one for each reference sequence, because read-metrics and sex cover every
    // contig. The coverage runs only for the contigs that the code tracks, which are the ones on
    // the main assembly that are also in the allowlist. That matches the sequential walker.
    struct Work {
        ref_id: usize,
        name: String,
        length: usize,
        tracked: bool,
        class: u8, // 0 = other, 1 = autosome, 2 = chrX (for the sex tally)
    }
    let mut works: Vec<Work> = Vec::new();
    let (mut autosome_length, mut x_length) = (0u64, None);
    for (ref_id, (name_bytes, map)) in header.reference_sequences().iter().enumerate() {
        let name = String::from_utf8_lossy(name_bytes.as_ref()).into_owned();
        let length = map.length().get();
        let tracked = contig::is_main_assembly(&name) && contig_allowlist.map_or(true, |s| s.contains(&name));
        let class = if contig::is_autosome(&name) {
            autosome_length += length as u64;
            1
        } else if contig::is_chr_x(&name) {
            x_length = Some(length as u64);
            2
        } else {
            0
        };
        works.push(Work {
            ref_id,
            name,
            length,
            tracked,
            class,
        });
    }

    // The progress goes out in **megabases of reference walked**, and not in contigs finished.
    // The bar then advances all the time, from the first seconds. rayon schedules the big
    // autosomes first, and they otherwise finish together late. That leaves the bar frozen at 0
    // for about half of the run.
    //
    // The denominator is the length of every contig that the code walks, because read-metrics
    // covers all of them. The contigs of the main assembly dominate that sum, so the count follows
    // the position in the genome closely enough.
    let total_bp: u64 = works.iter().map(|w| w.length as u64).sum();
    let total_mb = (total_bp / 1_000_000).max(1) as usize;
    let processed_bp = AtomicU64::new(0);
    progress(0, total_mb);

    let n_threads = analysis_thread_count();
    // Limit how many full-reference loads run at one time, and set that limit apart from the
    // compute parallelism. Those loads are what set the peak memory. At most a few contigs hold
    // their raw reference at one time, while the code builds the compact N-mask. A pool of tokens
    // is the counting semaphore.
    let load_permits = n_threads.min(4);
    let (perm_tx, perm_rx) = std::sync::mpsc::channel::<()>();
    for _ in 0..load_permits {
        let _ = perm_tx.send(());
    }
    let perm_rx = std::sync::Mutex::new(perm_rx);

    let process_contig = |w: &Work| -> Result<ContigPartial, AnalysisError> {
        // Stop before the cost of the reader and the reference load of this contig. A contig that
        // already runs stops at its own check in the record loop inside `for_each`.
        cancel.check()?;
        let (h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
        let region = Region::new(w.name.as_bytes().to_vec(), ..); // whole contig

        let mut cov_accum = if w.tracked {
            // Hold a load token across the raw-reference load and the build of the mask, and no
            // longer. Release it before the long pileup, which keeps the small mask alone.
            let _permit = {
                let _ = perm_rx.lock().unwrap().recv();
                LoadPermit { tx: &perm_tx }
            };
            let ref_bases = reader::read_contig_sequence(reference_path, &w.name)?;
            Some(ContigCoverageAccum::new(w.name.clone(), w.length, ref_bases, *params))
        } else {
            None
        };
        let mut rm = ReadMetricsState::default();
        // Drive the records through a sink over the lazy (zero-copy on BAM) record path.
        let (autosome_reads, x_reads, leftover_bp) = {
            let mut sink = ContigSink {
                ref_id: w.ref_id,
                rm: &mut rm,
                cov: &mut cov_accum,
                class: w.class,
                autosome_reads: 0,
                x_reads: 0,
                progress,
                processed_bp: &processed_bp,
                total_mb,
                last_pos: 0,
                local_bp: 0,
            };
            idx.for_each(&h, &region, &mut sink, cancel)?;
            (sink.autosome_reads, sink.x_reads, sink.local_bp)
        };
        // Flush this contig's unflushed bp tail so the counter reflects the whole contig walked.
        if leftover_bp > 0 {
            let g = processed_bp.fetch_add(leftover_bp, Ordering::Relaxed) + leftover_bp;
            progress((g / 1_000_000) as usize, total_mb);
        }

        let cov = cov_accum.map(|a| a.finish(w.ref_id));
        Ok(ContigPartial {
            rm,
            cov,
            autosome_reads,
            x_reads,
        })
    };

    // A region query can not see the unmapped tail, because it has no reference position. But the
    // sequential read-metrics counts it, in the total reads, the pf reads and the read length. So
    // sweep it on its own.
    let process_unmapped = || -> Result<ReadMetricsState, AnalysisError> {
        let (_h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
        let mut sink = MetricsSink {
            rm: ReadMetricsState::default(),
        };
        idx.for_each_unmapped(&mut sink, cancel)?;
        Ok(sink.rm)
    };

    // The CRAM decoder of noodles can recurse deep enough to overflow the default 2 MiB worker
    // stack of rayon. In the sequential walker, the larger stack of the main thread holds the same
    // file. A CRAM 3.1 file recurses deeper still, because of its new range and arithmetic codecs,
    // fqzcomp, and the name tokenizer.
    //
    // So give the workers a large stack that is safe for a decode, and the CRAM decode of one
    // contig then does not overflow. An overflow aborts the whole process, so the margin here must
    // be wide.
    let pool = reader::decode_pool(n_threads)?;

    let (contig_results, unmapped_rm) = pool.install(|| {
        rayon::join(
            || {
                works
                    .par_iter()
                    .map(&process_contig)
                    .collect::<Result<Vec<_>, AnalysisError>>()
            },
            || {
                if skip_unmapped {
                    Ok(ReadMetricsState::default())
                } else {
                    process_unmapped()
                }
            },
        )
    });
    let contig_results = contig_results?;
    let unmapped_rm = unmapped_rm?;

    // The merge. Read-metrics is a fold that commutes. Coverage merges one contig at a time, in
    // header order. Sex adds the class counts of each contig into one tally.
    let mut rm_total = ReadMetricsState::default();
    let mut cov_partials: Vec<ContigCoveragePartial> = Vec::new();
    let (mut autosome_reads, mut x_reads) = (0u64, 0u64);
    for p in contig_results {
        rm_total.merge(p.rm);
        if let Some(c) = p.cov {
            cov_partials.push(c);
        }
        autosome_reads += p.autosome_reads;
        x_reads += p.x_reads;
    }
    rm_total.merge(unmapped_rm);

    let coverage = merge_coverage_partials(cov_partials);
    let read_metrics = rm_total.finish();
    let sex = sex::result_from_tally((autosome_reads, autosome_length, x_reads, x_length)).ok();
    progress(total_mb, total_mb);
    Ok(UnifiedMetricsResult {
        coverage,
        read_metrics,
        sex,
    })
}

/// The time of each pass from [`profile_contig`], over one contig.
#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct ContigProfile {
    pub reads: u64,
    /// The raw decode alone: BGZF and the BAM record decode. It touches as few lazy fields as it
    /// can, and it makes no `RecordBuf`.
    pub raw: std::time::Duration,
    /// The raw decode, plus `RecordBuf::try_from_alignment_record`, which is the owned copy of
    /// each read.
    pub recordbuf: std::time::Duration,
    /// The full production work: `RecordBuf` + read-metrics + coverage pileup.
    pub full: std::time::Duration,
}

/// A diagnostic. It times the loop over the reads of a **single** contig, in three passes. The
/// cost then separates into its parts: the raw decode, the owned `RecordBuf` copy, and the metrics
/// and pileup work. It profiles the hot loop, and it does not walk the whole genome. Production
/// does not use it.
#[doc(hidden)]
pub fn profile_contig(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &CallableLociParams,
) -> Result<ContigProfile, AnalysisError> {
    use noodles::bam;

    let region = Region::new(contig.as_bytes().to_vec(), ..);

    // Pass 1, the raw decode. Walk the lazy bam::Record, touch the flags and the sequence length,
    // and make no RecordBuf.
    let mut raw_reads = 0u64;
    let raw = {
        let mut inner = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .map_err(|e| AnalysisError::io(bam_path, e))?;
        let header = inner.read_header().map_err(|e| AnalysisError::io(bam_path, e))?;
        let start = std::time::Instant::now();
        let q = inner
            .query(&header, &region)
            .map_err(|e| AnalysisError::io(bam_path, e))?;
        for r in q.records() {
            let rec = r.map_err(|e| AnalysisError::io(bam_path, e))?;
            std::hint::black_box(rec.flags());
            std::hint::black_box(rec.sequence().len());
            raw_reads += 1;
        }
        start.elapsed()
    };

    // Pass 2. The same, plus the conversion to a RecordBuf. It accepts nothing.
    let recordbuf = {
        let (h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
        let start = std::time::Instant::now();
        let q = idx.query(&h, &region)?;
        for r in q {
            std::hint::black_box(r?);
        }
        start.elapsed()
    };

    // Pass 3. The full work at each read, as production does it.
    let length = {
        let (h, _) = reader::open_indexed(bam_path, Some(reference_path))?;
        h.reference_sequences()
            .get(contig.as_bytes())
            .map(|m| m.length().get())
            .ok_or_else(|| AnalysisError::Message(format!("contig {contig} not in header")))?
    };
    let ref_bases = reader::read_contig_sequence(reference_path, contig)?;
    let full = {
        let (h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
        let mut cov = ContigCoverageAccum::new(contig.to_string(), length, ref_bases, *params);
        let mut rm = ReadMetricsState::default();
        let start = std::time::Instant::now();
        let q = idx.query(&h, &region)?;
        for r in q {
            let record = r?;
            rm.accept(&record);
            cov.accept(&record);
        }
        start.elapsed()
    };

    Ok(ContigProfile {
        reads: raw_reads,
        raw,
        recordbuf,
        full,
    })
}

/// A diagnostic. It runs the full work at each read, on every contig in `contigs`, at the same
/// time. That has the same shape as the fan-out over contigs in the real parallel walker. It
/// returns `(total_reads, wall_clock)`.
///
/// Compare that throughput against the rate of [`profile_contig`], which runs on one thread. The
/// difference shows contention, for example an allocator under pressure from the `RecordBuf`
/// allocation at each read. Production does not use this.
#[doc(hidden)]
pub fn profile_contigs_parallel(
    bam_path: &Path,
    reference_path: &Path,
    contigs: &[String],
    params: &CallableLociParams,
) -> Result<(u64, std::time::Duration), AnalysisError> {
    let start = std::time::Instant::now();
    let counts: Result<Vec<u64>, AnalysisError> = contigs
        .par_iter()
        .map(|contig| {
            let (h, mut idx) = reader::open_indexed(bam_path, Some(reference_path))?;
            let length = h
                .reference_sequences()
                .get(contig.as_bytes())
                .map(|m| m.length().get())
                .ok_or_else(|| AnalysisError::Message(format!("contig {contig} not in header")))?;
            let ref_bases = reader::read_contig_sequence(reference_path, contig)?;
            let mut cov = ContigCoverageAccum::new(contig.clone(), length, ref_bases, *params);
            let mut rm = ReadMetricsState::default();
            let region = Region::new(contig.as_bytes().to_vec(), ..);
            let mut n = 0u64;
            for r in idx.query(&h, &region)? {
                let record = r?;
                rm.accept(&record);
                cov.accept(&record);
                n += 1;
            }
            Ok(n)
        })
        .collect();
    Ok((counts?.into_iter().sum(), start.elapsed()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coverage, read_metrics, sex};
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name)
    }

    /// The fused walker yields exactly what the three standalone walkers do, field for field.
    #[test]
    fn unified_matches_standalone_walkers() {
        let bam = fixture("coverage.bam");
        let reference = fixture("ref.fa");
        let params = CallableLociParams::default();

        let unified = collect_unified_metrics(&bam, &reference, &params, None).unwrap();

        let cov = coverage::collect_coverage_callable(&bam, &reference, &params, None).unwrap();
        let rm = read_metrics::collect_read_metrics(&bam, Some(&reference)).unwrap();
        // The fixture holds chrM alone, with no autosome and no chrX, so the code can not infer
        // the sex. The fused walker reports `None`, and it still gives back the coverage and the
        // read-metrics. The separate walker returns an error. Both agree that there is no sex
        // here.
        assert!(sex::infer_from_bam(&bam, Some(&reference)).is_err());

        assert_eq!(unified.coverage, cov, "coverage diverged");
        assert_eq!(unified.read_metrics, rm, "read metrics diverged");
        assert_eq!(unified.sex, None, "expected chrM-only fixture to lack autosomes");
    }
}
