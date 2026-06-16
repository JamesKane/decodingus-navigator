# SIMD Optimization Targets — `navigator-analysis`

**Status:** Analysis only. No code changes. Report scoped to the `rust-rewrite` branch.
**Date:** 2026-06-15
**Scope:** Hot numeric / per-base loops in `crates/navigator-analysis/`, plus the build
configuration that gates autovectorization.

---

## TL;DR

1. **Do the cheap thing first.** There is **no `.cargo/config.toml`** and no `target-cpu`
   setting, so release builds target baseline `x86-64` (SSE2-era). Some wins are available
   *for free* from the autovectorizer once it targets the actual baseline CPU
   (`target-cpu=ivybridge` — see below; **not** `x86-64-v3`, which the baseline can't run).
   Measure that before hand-writing a single intrinsic.
2. **The genuinely hot loops are per-base genome-scale walkers** (coverage pileup, caller
   tally). They are good SIMD shapes (u8 threshold → accumulate) but are **partly bottlenecked
   by noodles record decode and the sliding-window data structure**, so arithmetic SIMD has a
   bounded ceiling. Profile to confirm the arithmetic is actually on the critical path.
3. **One loop has a bigger *algorithmic* win than a SIMD win:** the IBD sliding-window IBS
   count (`ibd.rs`) recomputes the whole window at every position — O(n·w). A rolling update
   makes it O(n) and dwarfs any vectorization gain.
4. **The ancestry math is hand-rolled scalar** (nalgebra is only linked into the offline
   `navigator-panelbuild`, not the runtime crate), but those loops run **once per sample**, so
   they are correctness-clean low-priority targets, not hot paths.

---

## 0. Build configuration (prerequisite, highest leverage / lowest risk)

Current state:

- No `.cargo/config.toml` / no `RUSTFLAGS` → default codegen for the host triple, which on
  x86-64 means **SSE2 only**. AVX is *not* emitted by the autovectorizer.
- `[profile.release]` already sets `lto = "thin"`, `codegen-units = 1` — good, keep it.
- Toolchain is nightly (`rustc 1.96.0-nightly`), so `std::simd` (portable SIMD) **is**
  available if we want it. Workspace `edition = 2021`, `rust-version = "1.80"`.

### Baseline hardware: Mac Pro 2013 (Xeon E5 **v2 = Ivy Bridge-EP**)

The chosen reference machine pins the ISA ceiling, so this is no longer an open product call —
it is a hard constraint. Ivy Bridge-EP supports:

- ✅ **AVX** — 256-bit *floating-point* vectors (`f64x4`, `f32x8`), SSE4.2, POPCNT, F16C
- ❌ **AVX2** — no 256-bit *integer* vectors
- ❌ **FMA** — no fused multiply-add
- ❌ BMI2

In `target-cpu` terms that is **`x86-64-v2` + AVX**, i.e. `-C target-cpu=ivybridge`. It does
**not** satisfy `x86-64-v3` (which requires AVX2 + FMA + BMI2). **Do not** use
`target-cpu=x86-64-v3` — it would emit AVX2/FMA instructions that fault (`SIGILL`) on the
baseline machine.

**The asymmetry this creates is the single most important planning fact in this report:**

| Loop family | Element type | Width on Ivy Bridge | Consequence |
|-------------|--------------|---------------------|-------------|
| Per-base coverage / caller (T1.1–T1.3) — the *hot* loops | `u8` / `i8` **integer** | **128-bit SSE only** (16× u8) — AVX integer needs AVX2 | Capped at SSE width; AVX buys nothing here |
| IBD IBS count (T2.1) | `i8` **integer** | **128-bit SSE only** | Same cap — algorithmic fix matters more |
| Ancestry / genotype likelihood (T3) | `f64` **float** | **256-bit AVX** (4× f64), but no FMA fusion | Full AVX width, run once per sample |

So on this baseline the **hot integer loops are the ones AVX can't widen**, while the loops that
*do* get 256-bit AVX (the float ancestry math) are the cold per-sample ones. That inverts the
usual intuition and pushes the priority order toward **(a) loop hygiene + SSE autovectorization**
and **(b) rayon parallelism** rather than wide-vector integer SIMD.

### Rayon is the bigger lever on this box

The Mac Pro 2013 range spans **4 to 12 physical cores** — E5-1620 v2 (4C/8T) at the floor,
through E5-1650 v2 (6C), E5-1680 v2 (8C), up to E5-2697 v2 (12C/24T) — all on quad-channel
DDR3-1866. Because the baseline must assume the **4-core floor**, but the top unit gives 12C/24T,
the parallel speedup curve is what varies most across the supported fleet — far more than any
SIMD width does. The existing per-contig parallel walker already exploits this, and adding
cores scales the hot *integer* loops further than a 128-bit-capped SSE pass would. Spend effort
keeping the rayon fan-out saturated and memory-bandwidth-friendly (quad-channel DDR3 is the
shared ceiling — a 12-core part can starve on a bandwidth-bound per-base scan) before
hand-writing integer SIMD. Validate scaling on the 4-core floor, not just a 12-core dev box.

### Idle GPU compute (out of scope, but noted)

Every Mac Pro 2013 also ships **dual AMD FirePro GPUs** — D300/D500/D700 depending on config
(the D700 is Tahiti-class, ~3.5 TFLOPS FP32 and notably strong FP64 at ~1/4 rate, 6 GB GDDR5
each). On the premium units that is substantial data-parallel compute sitting completely idle
during analysis. This is **explicitly out of scope** for a SIMD report and is a research
direction, not a near-term lever, but worth recording because it changes the long-run calculus:
the embarrassingly-parallel-over-sites work (genotype-likelihood tallies, panel genotyping,
all-pairs IBS) is a more natural GPU fit than the 128-bit-capped integer SSE paths the CPU is
stuck with here. **Caveats that make it a hard sell:** these are old GCN1 cards on a platform
Apple is sunsetting; Apple deprecated **OpenCL** (last good support ~Monterey) and never gave
these cards strong **Metal** compute support; noodles BAM/CRAM decode stays CPU-bound either
way; and the sliding-window data dependencies in the coverage/IBD loops don't map cleanly to a
GPU. Treat as a "someday, if a GPU compute path is ever justified" footnote, not a TODO.

Recommended experiments, in order:

| Step | Action | Risk | Expected |
|------|--------|------|----------|
| 0a | Build with `RUSTFLAGS="-C target-cpu=ivybridge"` and benchmark a WGS run | Low — matches the baseline machine exactly | Free SSE autovec + AVX on the float loops; 1.1–1.3× where loops are clean |
| 0b | Try `target-cpu=native` on the dev machine to see the *upper* bound (then ignore any gain that came from AVX2/FMA — the baseline can't use it) | None (local only) | Diagnostic: shows how much is AVX2-gated and therefore *unavailable* in production |
| 0c | Restructure the hottest loop bodies so LLVM can vectorize them at SSE width (hoist slice bounds, drop per-element `unwrap_or`, expose contiguous slices) | Low | Unlocks 0a on the integer loops — the main realistic win |

> **Distribution note:** pinning to `ivybridge` is safe for the stated baseline and for every
> newer x86 Mac. If you later want newer x86 Macs (Haswell+) to use AVX2/FMA, do **runtime
> feature detection** (`is_x86_feature_detected!("avx2")`) with an Ivy-Bridge SSE/AVX fallback —
> a single binary that lights up the wider path only where present. **Apple Silicon caveat:** a
> universal binary's ARM slice uses NEON natively, but the x86 slice **under Rosetta 2 does not
> support AVX at all** — so AVX paths must always have an SSE-or-scalar fallback regardless.

**Recommendation:** treat §0 as the actual first deliverable. Establish a repeatable benchmark
(a fixed BAM + `examples/profile_metrics` already exists per the metrics-walker work) and
record the autovectorized baseline before committing to intrinsics.

---

## Tier 1 — Hot, genome-scale, good SIMD shape

### T1.1 Coverage per-base pileup — `coverage.rs:574` (`feed_into_contig`)

```rust
for i in 0..len {                         // len = CIGAR M/=/X run
    let pos = ref_pos + i;
    if pos >= 1 && pos <= c.length {
        let base_q = quals.get(query_off + i).copied().unwrap_or(0);
        c.add(pos, base_q, mapq, params); // threshold + per-position accumulate
    }
}
```

- **Hotness:** every aligned base of every passing read, whole genome. This is *the* hottest
  arithmetic loop in the crate.
- **Shape:** `u8` base-quality scan, threshold (`base_q >= min_base_quality`,
  `mapq >= min_mapping_quality`), conditional accumulate into a sliding window column.
- **SIMD potential:** the threshold/compare half vectorizes cleanly (compare 16/32 `u8` quals
  against a constant → mask). The *accumulate* half scatters into `c.add` per position, which
  is the harder part.
- **Caveats / honest ceiling:**
  - `quals.get(...).copied().unwrap_or(0)` does a bounds-check + branch **per base** — this
    alone likely blocks autovectorization. Hoisting the slice bound out of the loop (§0c) is a
    prerequisite and may capture much of the win with zero intrinsics.
  - `c.add` mutates a `VecDeque`/window structure; inlining it and exposing a contiguous
    `&mut [Col]` slice for the run would let the compiler vectorize the column updates.
  - Real bottleneck is shared with noodles decode + `pileup_with` closure overhead — **profile
    first** (the prior `metrics-walker-perf` work found `RecordBuf::try_from` was 50% of CPU,
    not the arithmetic).

### T1.2 Callable-interval per-base dual-threshold — `coverage.rs:834` (`callable_intervals`)

- Same shape as T1.1 but **two** threshold comparisons per base (mapq + baseq) and three
  counter updates (`depth`, `qc_pass`, `low_mapq`) into a `VecDeque<Col>` window.
- **Hotness:** per-base genome-wide whenever callable BED is requested (query-time critical).
- **SIMD potential:** higher arithmetic density than T1.1 (two compares) → better ratio of
  vectorizable work to scatter. Same `VecDeque`-contiguity caveat.

### T1.3 Caller ACGT base-count histogram — `caller.rs:403` (`tally_region`) & `:311` (`tally_targets`)

```rust
for i in 0..len {
    let pos = ref_pos + i;
    if pos >= lo && pos <= hi {
        let base_q = quals.get(query_off + i).copied().unwrap_or(0);
        if base_q >= params.min_base_quality {
            if let Some(bi) = seq.get(query_off + i).and_then(base_index) {
                counts[pos - lo][bi] += 1;     // [u32; 4] per position
            }
        }
    }
}
```

- **Hotness:** dense de-novo calling walks every base of every read in 8 MB chunks (rayon-parallel
  across chunks). `tally_targets` is the same body but sparse (haplogroup sites).
- **Shape:** `u8` quality threshold + `base_index` table-lookup (A/C/G/T → 0..3) + scatter into a
  4-wide histogram keyed by position.
- **SIMD potential:**
  - The `base_index` lookup vectorizes via a small shuffle/`pshufb`-style table.
  - The threshold mask vectorizes.
  - The scatter into `counts[pos-lo][bi]` is the obstacle — it is a per-position, per-base
    histogram. A SIMD-friendly restructure: separate the 4 bases into 4 SoA count planes and
    accumulate with masked adds, or transpose so the inner loop is over a fixed small width.
- **Note:** already parallel by chunk, so SIMD stacks multiplicatively with the existing rayon
  fan-out.

---

## Tier 2 — Algorithmic win > SIMD win

### T2.1 IBD sliding-window IBS fraction — `ibd.rs:332` (`find_candidate_segments`)

```rust
for i in 0..n {                                   // n = SNPs on chromosome
    let look_back = i.saturating_sub(half);
    let look_forward = (i + half).min(n - 1);
    let (mut local_ibs2, mut local_ibs0, mut local_total) = (0, 0, 0);
    for &s in &ibs[look_back..=look_forward] {    // window of ~window_size SNPs
        match s { 2 => local_ibs2 += 1, 0 => local_ibs0 += 1, _ => {} }
        local_total += 1;
    }
    ...
}
```

- **Hotness:** per chromosome per pair; the inner window re-scans `window_size` elements at every
  position → **O(n · window_size)**. For database-scale all-pairs IBD this is the dominant cost.
- **The real fix is algorithmic, not SIMD:** this is a sliding window with a *fixed* stride of 1.
  Maintain rolling `ibs0`/`ibs2`/`total` counters, subtracting the element leaving the window and
  adding the one entering → **O(n)**. That is a constant-factor *and* complexity win that no SIMD
  can match.
- **If SIMD is still wanted** after the rolling rewrite: the IBS state array is `i8 ∈ {0,1,2}`;
  counting matches against a constant across the whole array is a textbook `pcmpeqb` + popcount.
  But do the rolling-counter rewrite first; it likely removes this from the hot list entirely.

---

## Tier 3 — Hand-rolled scalar, but per-sample (low priority)

All of these are in `ancestry.rs`, run **once per sample**, and use no nalgebra (nalgebra is
linked only into the offline `navigator-panelbuild::pca`, not this runtime crate). They are
clean, correct, and small — vectorize only if profiling a single-sample ancestry report shows
them, which is unlikely.

| Loop | Location | Shape | Why low priority |
|------|----------|-------|------------------|
| AF likelihood | `ancestry.rs:492` (`estimate_by_allele_frequency`) | `f64` reduction over `sites × pops` | ~10–20K ops/sample, one-shot |
| Admixture EM | `ancestry.rs:572` (`estimate_admixture`) | `f64` dot products, 500 EM iters × sites × k | Heaviest of the group (~10M FLOP/sample) — the *one* here worth a look if ancestry latency matters |
| PCA projection | `ancestry.rs:141` (`project_pca`) | `f64` matrix-vector, `sites × components` | Iterator overhead dominates; a nalgebra `DVector` would help more than intrinsics |
| Mahalanobis | `ancestry.rs:162` (`mahalanobis_sq`) | `f64` vector reduction over components | ~10 calls/sample |
| nMonte Frank–Wolfe | `ancestry.rs:675` (`nmonte_fit`) | repeated `f64` dot products, ≤1000 iters | Only on the G25_NMONTE method path |

> If ancestry latency ever does matter, the highest-value move is **not** intrinsics but linking
> `navigator-analysis` against nalgebra (already a workspace dep) and expressing the admixture
> E-step / PCA projection as matrix ops — nalgebra’s `matrixmultiply` backend autovectorizes and
> is far less error-prone than hand SIMD.

---

## Tier 4 — Not worth SIMD (documented to forestall re-investigation)

- **Genotype likelihood** `genotype.rs:57` / `:119`: inner dimension is ploidy (≤2) or
  `n*(n+1)/2` for a handful of alleles, and the body has a `10f64.powf()` per observation. The
  win here is a **phred→error lookup table** (256-entry `[f64; 256]`), not SIMD — and that table
  amortizes far more than vectorizing a ploidy-2 loop.
- **Banded DP alignment** `realign.rs:58`, `mtvariants.rs:213`: anti-diagonal data dependencies
  + control flow on the argmin move. SIMD-able only with a wavefront/striped rewrite (Farrar-style),
  which is a large effort for code that runs once per mtDNA / per indel window. Not justified.
- **N-mask construction** `coverage.rs:189`: once per contig; the per-position `is_n` *query* is
  the hot side and is already a single shift+mask. Fine as is.
- **Histogram/MAD finalization** `coverage.rs:918`/`:931`, **heteroplasmy top-two**
  `heteroplasmy.rs:73`, **right-align/rotation** `mtvariants.rs`: all once-per-result or
  fixed-tiny-width. Leave alone.

---

## Recommended sequence

1. **§0 build flags + benchmark harness.** Establish an autovectorized baseline
   (`target-cpu=x86-64-v3`) on a fixed BAM. This is the control group for everything else.
2. **§0c loop hygiene** on T1.1/T1.2/T1.3 — hoist slice bounds, expose contiguous slices,
   drop per-element `unwrap_or` — and re-measure. Much of the "SIMD" win is really "let LLVM
   vectorize what it already wants to."
3. **T2.1 rolling-counter rewrite** — pure algorithmic, no portability cost, likely the single
   biggest IBD speedup.
4. **Only then**, if profiling still shows arithmetic on the critical path, hand-write SIMD for
   the chosen 2–3 kernels (T1.3 base histogram is the best-shaped) behind runtime feature
   detection with a scalar fallback.
5. **Phred→error LUT** in `genotype.rs` as an independent, low-risk cleanup.

### Portability decision — RESOLVED

Minimum target CPU is fixed by the recommended baseline machine: **Mac Pro 2013, Xeon E5 v2
(Ivy Bridge-EP)** → `target-cpu=ivybridge` (`x86-64-v2` + AVX, **no AVX2/FMA**). Consequences,
already folded into §0:

- The **hot integer loops (T1, T2) are capped at 128-bit SSE** on this hardware — AVX does not
  widen integer ops. Treat wide integer SIMD as out of scope for the baseline; lean on loop
  hygiene + SSE autovec + rayon instead.
- The **float loops (T3) get 256-bit AVX** but no FMA, and they're per-sample/cold — low ROI.
- Any AVX2/FMA path for newer x86 Macs must be **runtime-detected with an Ivy-Bridge fallback**,
  and any AVX path needs an SSE/scalar fallback because **Rosetta 2 (x86-on-Apple-Silicon) has
  no AVX**.

---

## Appendix — files referenced

- `crates/navigator-analysis/src/coverage.rs` — pileup, callable intervals, N-mask, finalization
- `crates/navigator-analysis/src/caller.rs` — `tally_region`, `tally_targets`, `collect_bases`
- `crates/navigator-analysis/src/genotype.rs` — `call_genotype`, `call_genotype_multi`
- `crates/navigator-analysis/src/ibd.rs` — `find_candidate_segments`, `intersect_positions`
- `crates/navigator-analysis/src/ancestry.rs` — AF likelihood, admixture, PCA, Mahalanobis, nMonte
- `crates/navigator-analysis/src/realign.rs`, `mtvariants.rs`, `heteroplasmy.rs` — DP / fixed-width
- `crates/navigator-panelbuild/src/pca.rs` — the *only* nalgebra user (offline asset build)
- Build: root `Cargo.toml` `[profile.release]`; **no `.cargo/config.toml`** exists today
