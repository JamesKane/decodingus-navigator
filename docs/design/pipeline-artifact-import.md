# Pipeline-artifact batch import (fast-path ingest)

Status: design / proposed
Scope: `navigator-analysis` (scan + new `gvcf`/`sidecar` readers), `navigator-app`
(import + analyze orchestration), `navigator-ui` (worker job model), `navigator-store`
(provenance on cached analyses). Branch `rust-rewrite`.

## Motivation

The `ytree` workflow (`/Users/jkane/Genomics/ytree`) joint-genotypes chrY + mtDNA
across ~3,350 CHM13v2 samples. As a side effect it archives, **next to each source
CRAM on the NAS**, the per-sample intermediates that Navigator's own analysis would
otherwise recompute by walking the (often 10–12 GB) CRAM:

```
PRJEB31736/HG00096/
  HG00096.chm13.chrY.g.vcf.gz (+.tbi)   GATK HaplotypeCaller ploidy-1, non-PAR chrY:2,458,321-62,122,809 (males)
  HG00096.chm13.chrM.g.vcf.gz (+.tbi)   GATK HaplotypeCaller ploidy-1, full chrM
  HG00096.chm13.chrYM.callable.bed      CallableLoci --min-depth 4 --min-mapping-quality 20
  HG00096.chm13.chrYM.callable.summary.txt   per-state nBases (REF_N/CALLABLE/NO_COVERAGE/…)
  HG00096.chm13.sex                     "male" | "female"
  coverage.txt                          `samtools coverage` per-contig (meandepth, coverage%, covbases…)
  stats.txt                             `samtools stats` (SN summary lines)
  HG00096.chm13.cram (+.crai)
```

The single most expensive thing Navigator does per sample — `App::base_calls`
genotyping the CRAM at every Y-tree / mt-tree position (`caller::call_bases_at`) — is
**exactly** what the chrY/chrM GVCFs already contain. Reading a 3 MB GVCF instead of
walking a 12 GB CRAM turns per-sample haplogroup placement from minutes into
milliseconds, and makes a 3,350-sample import tractable.

### Goals

1. **Fast path for haplogroups.** When the sidecar GVCFs are present, place Y and mt
   from them; never touch the CRAM. This is the headline win.
2. **Fast path for the cheap metrics** that the pipeline already computed: sex
   (`.sex`), read metrics (`stats.txt`), and a *lite* coverage roll-up
   (`coverage.txt` + `callable.summary.txt`).
3. **Deferred deep analysis is additive.** Autosomal ancestry, the full coverage
   histogram (pct_10x/20x/median), SV, and IBD panel genotyping still require the
   CRAM. Running them later must **not** wipe the fast-path results. A re-analyze only
   recomputes a result when its *input fingerprint* changed.
4. **Non-blocking.** Import + fast-path ingest is cheap and returns quickly; the deep
   pass runs as a cancellable background job that yields, so the rest of the UI stays
   responsive.

### Decisions (locked)

- **Fast-path ingest is default-on** in `import_project_dir` when sidecars are present
  and the build matches. The old CRAM-only behavior stays reachable via a flag
  (`fast_path: bool` / a `--no-sidecars` CLI option).
- **The deep pass is manually triggered**, never auto-queued. Import does only the fast
  path and returns; the expensive CRAM walk (coverage histogram, ancestry, SV, IBD)
  runs only when the user invokes "Run deep analysis" (per-project or per-subject). This
  keeps a 3,350-sample import from silently launching a multi-day walk.

### Non-goals

- Re-deriving the joint tree / mask (that's the pipeline's job, off-box).
- Heteroplasmy, private-Y, and IBD from sidecars — those still want the CRAM.
- Trusting sidecars blindly across builds: the GVCFs are `chm13`/`hs1`; we record the
  build and only take the liftover-free path when tree build == GVCF build.

## Current state (what we're changing)

- `scan::scan` classifies `coverage.txt`/`stats.txt` as recognized-but-ignored, and
  lumps all `*.vcf.gz` into `variant_files`. It has no notion of "the chrY GVCF for
  this sample".
- `import_project_dir` creates Project→Biosample→SequenceRun→Alignment rows from
  `alignment_files` only. Variant/coverage/stats sidecars are dropped on the floor.
- `analyze_project` loops samples and, per sample, calls `run_coverage_for_alignment`,
  `assign_y_haplogroup`, `run_sex`, `run_read_metrics`, `run_sv` — **every one walks
  the CRAM**. mtDNA is intentionally skipped here.
- `assign_y_haplogroup` / `assign_mtdna_haplogroup_from_alignment` already have an
  input-fingerprint cache (`y_score_fingerprint` = alignment content hash ⊕ tree
  hash). We extend, not replace, this.
- Cached analyses are keyed `(alignment_id, kind, version)` via `save_analysis` /
  `load_analysis`. No provenance column today.

## Design

### 1. Scanner: typed sidecars (`navigator-analysis/src/scan.rs`)

Add an optional `sidecars: SampleSidecars` to `DiscoveredSample`, matched by suffix
against the sample's files (build segment is wildcard — `*.chrY.g.vcf.gz` etc.):

```rust
pub struct SampleSidecars {
    pub chr_y_gvcf:        Option<PathBuf>,  // *.chrY.g.vcf.gz
    pub chr_m_gvcf:        Option<PathBuf>,  // *.chrM.g.vcf.gz
    pub callable_bed:      Option<PathBuf>,  // *.callable.bed
    pub callable_summary:  Option<PathBuf>,  // *.callable.summary.txt
    pub sex:               Option<PathBuf>,  // *.sex
    pub coverage:          Option<PathBuf>,  // coverage.txt
    pub stats:             Option<PathBuf>,  // stats.txt
    pub build_hint:        Option<String>,   // parsed from the GVCF name segment (chm13 → hs1)
}
```

The general `variant_files` list stays (other importers use it). Sidecar detection is
purely additive and name-based; no file is opened during scan.

### 2. GVCF reader (`navigator-analysis/src/gvcf.rs`, new)

The core new primitive. Given a tabixed ploidy-1 GVCF, a contig, and a set of target
positions, return the *observed haploid base* at each target — the same shape
`base_calls` produces from the pileup, so it drops straight into `haplo::score`.

```rust
pub struct CalledBases {
    pub variant_bases: HashMap<i64, char>, // SNP ALT at sites with GT=1 (uppercase)
    pub callable:      HashSet<i64>,        // covered by a passing ref-block OR a passing variant
}
pub fn read_called_bases(
    gvcf: &Path, contig: &str, targets: &HashSet<i64>, p: &GvcfReadParams,
) -> Result<CalledBases, AnalysisError>;
```

GVCF semantics (verified against the real files):
- **Variant record** (`ALT != <NON_REF>`, `GT=1`): SNP → `variant_bases[pos] = alt[0]`;
  indel (`len(REF)!=len(ALT)`) → skip the ALT but still mark `callable`. Honor
  `DP`/`GQ` thresholds.
- **Ref block** (`ALT=<NON_REF>`, `GT=0`, `END=` in INFO): every position in
  `[POS, END]` that passes `MIN_DP`/`GQ` is `callable` (hom-ref). We do **not** synthesize
  a base here — see assembly below.

Reading: tabix `.tbi` query over the (clamped) target range. Use the existing VCF
plumbing — `du_bio::vcf` for line parsing, noodles tabix for the index — to avoid a
full-file scan. Targets outside the GVCF's region (e.g. PAR, off-mito) simply never
match → no-call, which is correct.

**Assembly into `calls: HashMap<i64,char>` (app side).** For each tree position `t`:
- `variant_bases[t]` present → `calls[t] = variant_bases[t]` (derived observed)
- else `t ∈ callable` → `calls[t] = ref_base[t]` — the **reference genome base** at `t`
- else omit (`no-call`)

> **Correction (validated against the real HG00096 GVCF).** The first design assumed a
> callable hom-ref site could take the tree's *ancestral* allele "because the reference
> base == ancestral on the native build, so no FASTA is needed." **This is wrong.** A
> reference genome is a *real* human Y deep in the tree — CHM13's Y is HG002 (haplogroup
> J1). At every backbone SNP that J1 and the sample both carry as derived, the sample
> matches the reference → the GVCF emits a ref block (hom-ref), and the true base there
> is the *derived* allele, not ancestral. Assuming ancestral mis-set the whole shared
> backbone and collapsed placement to the root. So the fast path **does** read the
> reference FASTA at the callable tree positions (`App::reference_bases`, one off-thread
> contig read) — exactly the base `caller::call_bases_at` reads off the reads. The
> reference is therefore required (recorded path, else gateway-resolved/cached).

This reconstructs exactly what `caller::call_bases_at` yields **for native-build
placement** (DecodingUs tree on `hs1`, mt FTDNA tree rCRS/direct — liftover-free for
these chm13 GVCFs). Cross-build (GVCF build ≠ tree build) falls back to the existing
CRAM path; we don't lift GVCF coordinates in v1. The CHM13 mt case *does* lift (tree
rCRS → `chrM`), reusing the existing rCRS↔chrM map and reading the reference base at the
lifted positions.

**Terminal selection — robust, not strict.** A joint-genotyped GVCF gives confident
calls that include a few stray ancestral contradictions on the deep backbone (recurrent
sites, the J1 reference, joint hard-filters). The strict `path_admissible` guard then
vetoes the genuine deep lineage and drops to a shallow node (HG00096 → A1b instead of
its true R1b1a1b1a1a, which `score` ranks top at 344/364). Same regime as BISDNA chip
data, so the GVCF path uses `assemble_assignment_robust`. **Validated: HG00096 →
R1b1a1b1a1a (344/364) in ~5 s vs a ~22 min CRAM walk.**

### 3. Sidecar metric parsers (`navigator-analysis/src/sidecar.rs`, new)

Pure text parsers → the existing result structs, so caching/UI/report are unchanged:

- **`.sex`** → `SexInferenceResult { inferred_sex, confidence: High, .. }` (ratios
  unknown → 0/sentinel; the label is what matters).
- **`stats.txt`** (samtools stats `SN` lines) → `ReadMetrics`: `total_reads` (raw total
  sequences), `pf_reads_aligned` (reads mapped), `pct_pf_reads_aligned`,
  `mean_read_length` (average length), `mean_insert_size` (insert size average),
  `std_insert_size` (insert size standard deviation), proper-pair %. Histograms
  (`read_length_histogram`, insert histogram) left empty — distributions need the BAM;
  the scalar QC the report shows is fully covered.
- **`coverage.txt` + `callable.summary.txt`** → a **lite** `CoverageResult`:
  `mean_coverage` = depth-weighted mean of per-contig `meandepth` (chr1..22,X,Y),
  per-contig stats from the table, `callable_bases` from the summary's CALLABLE rows
  (chrY+chrM only in this pipeline). `median_coverage`, `sd_coverage`, the depth
  histogram and `pct_Nx` are **not** derivable from samtools coverage → left as
  `f64::NAN`/empty and flagged partial (see provenance). The full histogram is a
  deep-pass product.

### 4. Provenance on cached analyses (`navigator-store`)

**Built (P5).** Two nullable columns on the analysis cache (migration `0017`):
`source TEXT` (`"pipeline-sidecar"` | `"navigator-walk"`) and
`completeness TEXT` (`"full"` | `"partial"`). `save_analysis_with_provenance` writes
them; `save_analysis` defaults to `navigator-walk`/`full` (existing call-sites
unchanged); `analysis_provenance(aln, kind, version)` reads them, defaulting `None`
columns (pre-provenance rows) to `navigator-walk`/`full`. This lets:
- the UI badge a "lite" coverage ("from pipeline; run deep coverage for histograms"),
- the deep pass know a `partial` coverage is upgradeable while a `full` one is not.

**Additive haplogroups — via the consensus gate, not cross-source fingerprints.** The
fast path records the Y/mt call (P4) with a `gv:`-prefixed fingerprint. The deep pass
(`analyze_project`) already skips Y when `haplogroup_consensus(Y).is_some()`, so a
sidecar-placed Y/mt is never re-walked. (We deliberately do *not* teach the CRAM-path
`assign_y_haplogroup` to recognize the `gv:` fingerprint — it has no GVCF path to
re-hash; a direct, explicit CRAM Y assignment is *meant* to supersede.) So the
fingerprint distinguishes provenance for the UI/audit; the skip is structural.

### 5. App orchestration (`navigator-app/src/lib.rs`)

New entry points, composing the pieces:

- `assign_y_from_gvcf(alignment_id, gvcf, build)` / `assign_mt_from_gvcf(...)`:
  parse tree → `read_called_bases` → assemble `calls` → `assemble_assignment[_robust]`
  → `record_call_fp` with the GVCF-based fingerprint. Mirrors `assign_y_decodingus`
  with the GVCF swapped in for `base_calls`.
- `ingest_sidecars(alignment_id, &SampleSidecars)`: best-effort, each independent —
  Y gvcf → Y call; M gvcf → mt call; `.sex` → sex; `stats.txt` → read metrics;
  `coverage.txt`(+summary) → lite coverage. Records what it could, logs what it
  couldn't. **This is the fast path.**
- `import_project_dir` is extended: after creating each Alignment, if the sample has
  sidecars and the alignment build matches the GVCF build, call `ingest_sidecars`
  inline (it's cheap — small text/GVCF reads, no CRAM, no hashing). Gated by a param
  `fast_path: bool` (default on) so the old behavior is still reachable.
- `analyze_project` (the deep pass) is additive (P5): coverage is re-run only when the
  cached result is `partial` (the lite sidecar coverage) — `analysis_provenance` gates
  it; a `full` walk is skipped. Y is skipped when a consensus call already exists (so a
  sidecar-placed Y is not re-walked); sex/read-metrics from sidecars are `full` and
  skipped; SV/ancestry/IBD are never sidecar-sourced so they always run. The lite
  coverage is overwritten in place by the full walk (same key, provenance → `full`).

### 6. Worker / non-blocking (`navigator-ui/src/worker.rs`)

- `ImportProjectDir` stays a single command but now returns quickly *with* a fast-path
  summary (`{samples, alignments, y_placed, mt_placed, sex, metrics, lite_coverage,
  deep_pending}`), because ingest is cheap.
- Add a **manually-triggered background deep-analyze job** (a "Run deep analysis"
  button, per-project and per-subject): a new streaming command (modeled on the
  existing `RunFullAnalysis`/`EstimateAncestry` spawn-loop handlers) that processes the
  target samples one at a time on a `tokio::spawn`'d task, emits per-sample progress
  events, honors a cancel flag, and `await`s between samples so other UI commands
  interleave. The fast path having already filled the report, this job is purely the
  opt-in "now do the expensive autosomal/coverage/SV work" pass — never auto-started,
  cancellable, and it leaves the haplogroups alone.

## Data-coverage matrix

| Result            | Sidecar source                         | Fast path | Completeness | Deep pass still needed for |
|-------------------|----------------------------------------|-----------|--------------|----------------------------|
| Y haplogroup      | `*.chrY.g.vcf.gz`                       | ✅        | full         | —                          |
| mt haplogroup     | `*.chrM.g.vcf.gz`                       | ✅        | full         | —                          |
| sex               | `*.sex`                                 | ✅        | full         | —                          |
| read metrics      | `stats.txt`                            | ✅        | full*        | length/insert histograms   |
| coverage roll-up  | `coverage.txt` + `callable.summary`    | ✅        | partial      | median, pct_Nx, histogram, genome-wide callable |
| autosomal ancestry| —                                      | ❌        | —            | CRAM genotyping @ AIM panel |
| SV                | —                                      | ❌        | —            | CRAM (≥10×)                |
| IBD               | —                                      | ❌        | —            | CRAM panel genotyping       |

\* the scalars the report/QC view shows; only the distributions are deferred.

## Validation

- Unit: GVCF reader against committed fixtures (ref-block span, SNP ALT, indel skip,
  out-of-region target, DP/GQ filtering). `stats.txt`/`coverage.txt`/`.sex` parsers
  against trimmed real fixtures.
- **Parity:** place HG00096 (and the known GFX donor) Y+mt from the sidecar GVCF and
  from the CRAM walk; assert the same terminal. This is the correctness gate for the
  whole fast path.
- Idempotence: fast-path import, then deep analyze, then re-import — assert the Y/mt
  calls are byte-identical across runs (fingerprint cache holds; no wipe).

## Rollout / phasing

1. **GVCF reader + parity test** (no wiring) — proves the fast path is correct in
   isolation.
2. Sidecar metric parsers + their unit tests.
3. Scanner sidecar detection.
4. App `assign_*_from_gvcf` + `ingest_sidecars`; wire into `import_project_dir`
   behind `fast_path`.
5. Provenance columns + fingerprint extension; make `analyze_project` additive.
6. Worker background deep-analyze job + UI summary/badges.

Each phase is independently testable and commits cleanly; phases 1–4 already deliver
the headline win (haplogroups + sex + metrics with no CRAM walk) before any of the
deep-pass refactor lands.
```

`/Volumes/nas/.../HG00096` is the running parity fixture.
