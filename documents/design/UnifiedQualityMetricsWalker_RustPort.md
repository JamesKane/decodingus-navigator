# Unified Quality Metrics Walker — Rust port plan

Last updated: 2026-06-10. Branch: `rust-rewrite`. Companion to the Scala-era
`UnifiedQualityMetricsWalker.md` (which describes the *old* htsjdk implementation —
Phases 1–3 "complete" refer to the Scala codebase, not this rewrite).

## Why (the actual gap in the rewrite)

The rewrite already split the work into focused single-pass walkers:

- `coverage::collect_coverage_callable` — one coordinate-ordered pileup pass over a
  coordinate-sorted BAM/CRAM. Sliding-window pileup → depth histogram + GATK callable
  states + samtools-style per-contig stats. This is **the slow whole-genome pass**
  (single-threaded full-genome pileup, minutes on a real WGS BAM).
- `read_metrics::collect_read_metrics` — one record-level pass: alignment-summary counts,
  read-length + insert-size distributions, pair orientation, mean MAPQ.
- `sex::infer_from_bam` — BAM reads per-reference counts from **BAI metadata** (O(contigs),
  no scan); CRAM has no per-reference counts in `.crai` so it does a **full record scan**.

In `run_full_analysis_streaming` (worker.rs) these run as separate steps. Net file reads:

| Input | coverage | sex | read_metrics | total end-to-end passes |
|-------|----------|-----|--------------|-------------------------|
| BAM   | 1 (pileup) | 0 (BAI) | 1 (scan) | **2** |
| CRAM  | 1 (pileup) | 1 (scan) | 1 (scan) | **3** |

CRAM decode is the expensive case (per-record reference reconstruction). Fusing the three
into one record loop is **3→1 for CRAM, 2→1 for BAM**, and removes the only place CRAM is
scanned solely for sex.

This is the same I/O-bound win the Scala doc targeted (2.5–3×), realized for the rewrite's
already-factored walkers rather than re-collapsing three GATK tools.

## What unification requires (the crux)

`collect_coverage_callable` pre-filters hard: it skips unmapped / secondary / supplementary
/ duplicate / qc-fail records and **only** processes main-assembly contigs (`is_main_assembly`).

`collect_read_metrics` needs the opposite — it must see **every** record (unmapped included,
for `total_reads`/`pf_reads`/`pct_pf_reads_aligned`; off-main-assembly included, for accurate
alignment + pairing rates). Sex needs per-contig **mapped** read tallies across autosomes +
chrX (header lengths give the denominators).

So the unified walker must **not** pre-filter at the loop head. One record loop:

1. **Every record** → feed read-metrics accumulators (the existing `read_metrics` classify
   logic verbatim: secondary/supp excluded from primary metrics, qc-fail → not pf, etc.).
2. **Mapped records** → tally per-contig read counts by class (autosome / chrX) for sex
   (works for BAM and CRAM alike — no BAI dependency in the fused path; see note below).
3. **Primary + mapped + main-assembly + passes coverage filters** → feed the sliding-window
   pileup (the existing `coverage` per-position machinery verbatim).

The two accumulator sets are independent; the only shared work is decode + the cigar walk
(which only the pileup branch needs). No metric changes — byte-for-byte the same numbers,
just collected together.

## Design

New module `crates/navigator-analysis/src/unified.rs`:

```rust
pub struct UnifiedMetricsResult {
    pub coverage: CoverageResult,        // reuse coverage.rs type as-is
    pub read_metrics: ReadMetrics,       // reuse read_metrics.rs type as-is
    pub sex: SexInferenceResult,         // reuse sex.rs type as-is
}

pub fn collect_unified_metrics(
    bam_path: &Path,
    reference_path: &Path,                // required: CRAM decode + ref-N detection
    params: &CallableLociParams,
    contig_allowlist: Option<&HashSet<String>>,
    progress: &mut dyn FnMut(usize, usize),
) -> Result<UnifiedMetricsResult, AnalysisError>;
```

Reuse, don't duplicate. Refactor the per-record bodies of the three existing walkers into
small `pub(crate)` helpers so both the standalone and fused entry points call the same code:

- `coverage.rs`: lift the per-record pileup feed (the cigar match + `CurContig::add`/
  `advance_to`/contig-transition + `Globals`) behind a `CoverageState` struct with
  `fn accept(&mut self, &RecordBuf)` + `fn finish(self) -> CoverageResult`. The existing
  `collect_coverage_callable*` become thin wrappers that filter, then drive `CoverageState`.
- `read_metrics.rs`: lift the per-record classify into `ReadMetricsState::accept` +
  `finish() -> ReadMetrics`. `collect_read_metrics` becomes a wrapper.
- `sex.rs`: add `SexState` accumulating the `Tally` from records (the CRAM `tally_via_scan`
  inner loop), plus `finish(header) -> SexInferenceResult` reusing `determine_sex`.

The fused loop holds all three states, dispatches each record, then assembles the result.
This keeps the three standalone functions (still used à la carte by `run_sv`, the
per-contig de-novo path, single-step UI commands) and adds the fused path with **zero**
metric divergence — the helpers are the single source of truth.

### Sex in the fused path

The fused walker tallies per-contig mapped reads directly (it's already touching every
record), so it does **not** use BAI. That's strictly more robust (no BAI dependency) and
matches the CRAM path's existing math. Standalone `infer_from_bam` keeps its BAI fast path
for the cheap single-step "Sex inference" command — only the full-analysis pipeline uses the
fused tally. Validate the two agree on the GFX BAM (they should: same reads, same classes).

### Memory + ordering

Unchanged from `coverage.rs`: requires a coordinate-sorted BAM/CRAM; peak memory is the
open-read span (sliding window) + one contig's reference bases. Read-metrics + sex
accumulators are O(1)/O(contigs). Progress callback fires per finalized contig as today.

## Wiring

- `lib.rs`: add `run_unified_metrics(alignment_id)` that runs the fused walker on a blocking
  thread and persists **all three** artifacts (`coverage`/`COVERAGE_VERSION`, `read_metrics`/
  `"1"`, `sex`/`"1"`) in one shot — same cache keys, so `cached_coverage`/`cached_read_metrics`/
  `cached_sex` and `run_sv`'s reuse logic all keep working untouched. Keep the inferred-sex
  write-back to the biosample (currently in `run_sex`).
- `worker.rs` `run_full_analysis_streaming`: replace steps 1–3 (Coverage / Sex / Read metrics)
  with a single "Quality metrics" step driving `run_unified_metrics` with the per-contig
  progress callback. Total drops 8→6. Steps 4–8 (SV, chrM de-novo, Y/mt haplogroup, ancestry)
  are unchanged and still read the cached artifacts. **Reuse the cached coverage** short-circuit
  exactly as today (don't re-scan if a cached result exists).
- Feature toggle parity with the Scala `experimental.unified-metrics-enabled` is unnecessary
  — there's no behavioral change to gate (identical numbers), so wire it directly.

## Validation

- Unit: fused result == running the three standalone walkers separately, on the existing
  fixtures (`coverage.bam`, `paired.bam`, `sex.bam` + their CRAM twins). Add to
  `tests/` a `unified.rs` asserting field-by-field equality.
- Live (`#[ignore]`, real data): extend `parity_real.rs` — assert `collect_unified_metrics`
  on the GFX0457637 CHM13 BAM matches the three standalone calls and still yields sex=Male,
  European-consistent coverage. Run via the HANDOFF env-var recipe.
- `cargo test --workspace` green; `cargo build` warning-free (workspace policy).

## Out of scope (deferred, as in the current walkers)

BED-interval output from the fused path (use `callable_intervals` à la carte), multi-threaded
per-contig parallelism (Scala Phase 4), direct SVG binning, CRAM reference-cache tuning.

## Estimated shape

~1 new module + 3 small refactors (extract per-record helpers) + 1 lib method + 1 worker
edit + 2 test files. No new deps, no schema/artifact changes, no metric changes. The risk is
entirely in the refactor preserving exact numbers — the equality tests are the guardrail.

---

## Update — threading (implemented 2026-06-10)

Two layers, both byte-identical to the sequential walker (guarded by live parity tests in
`crates/navigator-analysis/tests/parity_real.rs`):

1. **MT bgzf decompression** (`reader.rs`, commit 900381d). BAM sequential reads wrap the file
   in `bgzf::MultithreadedReader` — parallel block inflation, sequential record parsing. Only
   ~8% on the GFX BAM: **decompression is not the bottleneck, the per-position pileup compute
   is.** `NAVIGATOR_BGZF_THREADS` (default cores−1, cap 6; 1 disables).

2. **Per-contig parallel walker** (`unified::collect_unified_metrics_parallel`, commit 7a9d7f2).
   Coverage is independent per contig, so fan out over contigs with rayon, each running the same
   `CoverageState`/`ReadMetricsState`/sex accumulators, then merge (commutative folds +
   header-ordered per-contig outputs). Read-metrics covers **every** contig (region queries) plus
   an **unmapped-tail sweep** (`query_unmapped`) so totals equal the sequential pass exactly.
   **BAM + `.bai` only** (region + unmapped queries); CRAM (no `.crai` unmapped query) and
   unindexed BAM transparently fall back to the sequential walker.

   Memory (the parallel path's real constraint, since each contig task would otherwise pin its
   full reference): a **1-bit-per-base reference N-mask** replaces the raw bytes (coverage only
   needs N-detection) — ~31 MB for chr1 vs ~248 MB — and a **load semaphore** caps concurrent
   full-reference loads (≤4) independently of compute threads (default `min(cores, 12)` — the
   knee; past it wall time is floored by chr1 + the serial unmapped sweep).
   `NAVIGATOR_ANALYSIS_THREADS` tunes it.

   **Measured on the real 9 GB GFX0457637 pbmm2 CHM13 BAM: unified 64.7s → 12.6s (5.15×), peak
   RSS 2.8 GB** (≈ the sequential walker's footprint).

Further headroom (not done): the serial unmapped sweep and the single largest contig bound the
tail; splitting big contigs into sub-regions (with per-region coverage merge) or parallelizing
the unmapped sweep would push past ~5×, at more complexity.
