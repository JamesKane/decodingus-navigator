# DUNavigator → Rust Rewrite Plan

**Status:** Planning
**Date:** 2026-06-01
**Decisions locked:** egui/eframe GUI · purpose-built pure-Rust haploid caller (no GATK/JVM) · shared crates extracted, Navigator in its own Cargo workspace

---

## 1. Why this is more tractable than it looks

The sister repo `/Users/jkane/Development/decodingus` (web rewrite) has already built
the foundation we need. It is a clean Cargo workspace whose lower layers are
**desktop-reusable as-is**:

| decodingus crate | What it gives Navigator |
|:---|:---|
| `du-domain` | Pure domain types, strongly-typed IDs, serde enums — **zero IO**. Directly fixes our type-triplication debt. |
| `du-atproto` | DID/handle parsing, did:key verification, PDS resolution, **OAuth: PKCE, ES256, DPoP, PAR** — the exact work from the AT Proto OAuth plan. |
| `du-bio` | **Coordinate math + text-format parsing only** (pure Rust, IO-light): `callable` (BED interval merge), `liftover` (UCSC chain), `vcf` (text reader for catalog ingest), `ybrowse`. **BAM/CRAM I/O and variant calling are explicitly out of scope** — those are Navigator-side (`navigator-analysis`). Depends only on `du-domain` + `thiserror`; no noodles. |
| `du-db` (pattern) | sqlx query-module-per-aggregate + JSONB→domain mapping pattern to copy for SQLite. |

> **Scope note (resolved at extraction, 2026-06-01):** an earlier draft put the noodles
> BAM/CRAM I/O and the walker ports inside shared `du-bio`. That was narrowed: shared `du-bio`
> stays coordinate-math + text-parsing (so the decodingus server never pulls in noodles), and
> **all noodles raw-read I/O + the walkers live in Navigator's `navigator-analysis` crate.** The
> haploid caller was already Navigator-side. See §2 and §4a for the corrected ownership.

The web UI (`du-web`: Axum + Askama + **HTMX**, server-rendered) is **not** reusable as a
desktop UI — but its domain/business logic and the crates above are. Hence: **egui for
the UI, reuse the lower crates.**

Master architecture reference for the sister repo: `~/.claude/plans/robust-knitting-lampson.md`.

---

## 2. Target topology — shared crates + Navigator workspace

The genuinely shared crates are extracted from the decodingus repo into
**`decodingus-shared`**, a sibling repo at the Development root. All three projects sit
side-by-side so both apps depend on the shared crates and fixes flow both ways:

```
/Users/jkane/Development/
├─ decodingus/                  # web app (Axum + Askama + HTMX)
├─ decodingus-shared/           # extracted shared crates, versioned independently
│    du-domain/                 # types, IDs, enums  (no IO)
│    du-atproto/                # DID + OAuth/DPoP    (identity)
│    du-bio/                    # coordinate math + text parsing: liftover + callable + vcf/ybrowse (NO BAM/CRAM, NO caller)
└─ DUNavigator/                 # desktop app — moved up from the legacy scala/ parent
     crates/
       navigator-domain/        # desktop-only types not in du-domain
                                #   (SequenceRun, Alignment, AnalysisArtifact,
                                #    YProfile, IBD, Workspace/Project aggregate)
       navigator-store/         # local persistence: SQLite via sqlx
       navigator-analysis/      # noodles BAM/CRAM I/O + the walkers + the haploid caller
                                #   (uses du-bio for liftover/callable/coordinate math)
       navigator-sync/          # PDS push/pull + AsyncSync (completed, not stubbed)
       navigator-app/           # application/command layer (no UI, no ScalaFX analog)
       navigator-ui/            # egui/eframe front end (thin: render + dispatch)
     Cargo.toml                 # workspace; depends on ../decodingus-shared/* (path or git)
```

**Dependency rule (the antidote to today's god object):** dependencies point *down only*
— `ui → app → {analysis, store, sync} → domain`. The UI never touches the DB or a
processor directly; it dispatches commands and renders state. No crate below `app`
imports anything UI-related.

---

## 3. Persistence — SQLite via sqlx (replaces H2/Slick)

- **Engine:** SQLite (embedded, file-based, perfect for a desktop app), accessed through
  `sqlx` (already in the workspace; supports the `sqlite` feature). Cross-platform, no
  server, single file under `~/.decodingus/`.
- **Kill the Slick-22-tuple JSONB workaround.** sqlx has no tuple-arity limit, so model
  complex children as **proper rows** (e.g. `file_info`, `contig_metrics`, `ibd_segment`
  tables) instead of `asJson.noSpaces` blobs. Reserve JSON columns for genuinely
  schema-less payloads (AT Proto record snapshots), using `sqlx::types::Json<T>` —
  exactly the `du-db` pattern.
- **Migrations:** `sqlx migrate` (versioned, checked-in SQL), replacing the 12 forward-only
  Flyway-style files. Author down-migrations.
- **One source of truth.** Persisted state is authoritative; the UI reads from a
  projection. No more `_workspace.value = …` imperative mutation racing async H2 writes
  (the current race-condition hotspot in `WorkbenchViewModel`).

---

## 4. GATK / HTSJDK replacement strategy

`noodles` (pure Rust) replaces htsjdk for BAM/CRAM/VCF/FASTA/BGZF/index I/O — and crucially
gives us a **clean Windows build** that `rust-htslib` (C bindings) would not. **This noodles
I/O lives in Navigator's `navigator-analysis` crate, not in shared `du-bio`** (which stays
IO-light coordinate math + text parsing, per §1). Navigator depends on `du-bio` for the
liftover/callable/coordinate primitives and adds the raw-read layer on top.

### 4a. The walkers → mechanical port onto noodles pileups

These already reimplement GATK logic in Scala-over-htsjdk; they are pileup loops with
statistics and port almost directly. **All targets are in `navigator-analysis`**; the
`du-bio` references are the shared coordinate/text primitives it builds on:

| Scala component | Replaces (GATK) | Port target |
|:---|:---|:---|
| `CoverageCallableWalker` | `CollectWgsMetrics` + `CallableLoci` | `navigator-analysis::coverage` over `noodles` pileup, using `du-bio::callable` for the BED merge/summary |
| `UnifiedMetricsWalker` | `CollectAlignmentSummaryMetrics` + `CollectInsertSizeMetrics` | `navigator-analysis::read_metrics` |
| `SvEvidenceWalker` | (custom BreakDancer/Pindel-style) | `navigator-analysis::sv` |
| `SexInference` | (index-based ratio) | `navigator-analysis::sex` (noodles BAI metadata) |
| `LiftoverProcessor` | `LiftoverVcf` | `du-bio::liftover` (UCSC chain parse — pure coordinate math, shared) |
| `ReferenceQuerier` | htsjdk FASTA | `navigator-analysis` FASTA via noodles, contig caching |

### 4b. The genuine gap: variant calling — a purpose-built **haploid caller**

There is no pure-Rust GATK, and your flows need three things from it. **Confirmed
requirement:** de-novo Y/mtDNA calls are needed for **private-variant matching and new
branch creation** — so the caller is not just a force-call genotyper. It has two modes:

1. **Force-call genotyping at known sites** (haplogroup tree sites, ancestry-informative
   SNPs). Pileup at each target position → call ref/alt by depth + base-quality majority.
   Straightforward.
2. **De-novo discovery across the contig** (Y and mtDNA) → **private variants** for branch
   creation. Because Y and mtDNA are **haploid (ploidy = 1)**, this is *pileup-consensus
   calling*, not the diploid local-reassembly that makes `HaplotypeCaller`/`Mutect2`
   expensive to reproduce. Algorithm: walk callable positions, compute the consensus
   non-reference allele where depth ≥ min and base-quality/MAPQ filters pass and allele
   fraction ≥ threshold; subtract known tree positions → private set. mtDNA runs at high
   depth over 16.5 kb (cheap); Y is large but limited to reliably callable regions.

**The risk to validate, stated honestly:** SNPs are tractable; **indels and homopolymer
runs are where a naive pileup caller diverges from GATK** — local reassembly is exactly
what improves indel accuracy, and mtDNA homopolymer/length-heteroplasmy sites
(e.g. ~309, ~16193, the 3107 "N") are notorious. Mitigations, in order of preference:
   - Add a **light local realignment** around candidate indels before calling.
   - Restrict branch-defining/private calls to **SNPs** initially; treat indels as
     advisory until validated.
   - Keep a **pinned external caller fallback** (see §4d) usable for indel-sensitive
     mtDNA work during validation only.

`CollectWgsMetrics`, `BuildBamIndex`, `IndexFeatureFile`: not separate algorithms —
metrics fold into the coverage walker; indexing is a noodles call.

### 4c. Validation — parity is a first-class deliverable

Before the Rust caller is trusted, build a **golden-truth harness**: run the existing
GATK pipeline and the new Rust caller on a panel of samples (varied depth, platforms,
builds; include known-hard mtDNA homopolymer cases) and assert agreement on called
genotypes, private-variant sets, and haplogroup assignment. This harness is also the
regression guard for the rewrite. No flow flips to the Rust caller until it passes.

### 4d. Transition bridge (optional, time-boxed)

If the caller's validation lags the rest of the rewrite, the JVM GATK can run as a
**subprocess sidecar** behind the same `navigator-analysis` interface so UI/domain work
proceeds unblocked — explicitly temporary, removed once §4c passes. (This is *not* the
shipped architecture; the end state is JVM-free.)

### 4e. Self-referential callable loci (design note — not yet built)

**Idea.** Derive the callable-region mask from *the sample's own alignment* (a per-sample
callable-loci BED), rather than gating against a fixed external mask. The CallableLoci BED
is already an artifact of the Scala/GATK pipeline; the rewrite currently only *summarizes*
callability (per-position `CallableState` counts) — it should also **emit the CALLABLE
runs as a stored BED artifact** per alignment, and consume *that* downstream.

**Criterion.** Replace the fixed per-base aggregate (depth ≥ N) gate with a **run-length
requirement proportional to molecule length**: a CALLABLE run counts only if its length
`≥ f · fragment_length` (insert size for paired short reads; read length for single-end /
long reads — both already estimated by the `read_metrics` walker; `f` tunable, start ≈1.0).
A callable island shorter than the molecules covering it sits inside a repeat those
molecules straddle ambiguously, so it's dropped — a mappability proxy with **no external
mappability track**.

**Why it matters.**
- *Self-referential* — callability tracks the sample's real evidence and **sequencing tech**
  automatically; supersedes external masks (e.g. `b38_sites.bed`, which is a ~143 k-site
  allowlist — not discovery-capable). The per-sample callable *region* BED both removes Y
  repeat noise from the private bucket **and** lets genuinely novel positions inside
  reliable regions surface.
- *Rewards long reads* — short reads (~150 bp / ~500 bp fragment) clear the proportional
  bar over only modest stretches of chrY; **PacBio HiFi** (~15 kb molecules map uniquely
  across the repeats that defeat short reads → long CALLABLE runs) clears it over far more,
  so HiFi yields materially more callable Y and more callable SNPs. The advantage is
  emergent from the criterion, not special-cased.
- *Perf* — restricting the de-novo scan to the callable BED also collapses the whole-chrY
  pass (~13 min observed) to the reliable subset.

**Ingredients in place:** per-position `CallableState` classification (`coverage` walker);
read/fragment-length estimates (`read_metrics` walker). **Missing:** coalesce CALLABLE
runs → BED, the run-length gate, store as a per-alignment artifact, and route it into
`private_y_variants` / de-novo region restriction / force-call site selection in place of
external masks. `ExcessiveCoverage` (collapsed-repeat pileups) stays non-callable; the
run-length gate is additive. Validation: compare callable-chrY size for a short-read vs a
HiFi sample of the same individual.

> Note: `combBED.b38.bed` (Poznik ∩ Big Y original design) is **only** for replicating
> YFULL's age algorithms — not a general calling mask.

**Built + validated (status update).** `coverage::callable_intervals` (run-length-gated
CALLABLE BED, reference-free) + `estimate_molecule_lengths`; `app::callable_chr_intervals`
derives tech-adaptive params; `private_y_variants_self_masked` consumes the sample's own
BED; UI exposes self-referential / external-BED / none. The haploid **caller** depth is
now tech-adaptive too (long reads → `min_depth 2`), mirroring the mask — the two no longer
disagree.

Real-data validation (same individual, GRCh38):
- *Callable chrY:* WGS229 short-read ~13× → 13.9 Mb (3336 runs); GFX0457637 HiFi ~1.6× →
  1.1 Mb (70 runs ≈ one ~11 kb read each).
- *Private-Y bucket (WGS229), by mask:* none 4583 / external `b38_sites.bed` 177 /
  self-referential 1151 novel — self sits between noise and the over-restrictive allowlist
  and stays discovery-capable.
- *Low-coverage CCS for haploid calling:* with the caller still at `min_depth 4`, HiFi Y
  matched 283/1919 (score 0.306); after the adaptive `min_depth 2`, **1097/1919 (0.537)** —
  ≈4× more backbone SNPs, supporting "≈1× CCS suffices on single-ploidy regions." mtDNA
  resolved `U5a1b1g` exactly on both technologies.

Open knobs: `f` (run-length proportion) is hardcoded 1.0 — at low-coverage HiFi the gate is
coupled to coverage continuity, so `f < 1.0` may help; needs more HiFi data points (N is
small, and the HiFi here is the same person as WGS229). A depth-1 floor for pure CCS is
plausible but unvalidated.

### 4f. Reference & liftover management (BUILT — Y + mtDNA haplogroups validated on CHM13)

**Prerequisite for haplogroups on CHM13-aligned data.** The FTDNA Y haplotree is in
**GRCh38** coordinates and mtDNA is **rCRS**, while the sample CRAMs here are
**T2T-CHM13v2.0** (and CHM13 `chrM` is not necessarily rCRS). Closing that needs either
liftover of tree positions GRCh38→CHM13, or a CHM13 reference that carries rCRS mito for the
mtDNA path.

**BUILT (2026-06-03): the `navigator-refgenome` crate** — retrieval + on-disk cache of
reference FASTAs and liftover chains. Resolves a build name → a cached, decompressed,
**in-Rust-indexed** `.fa` (`noodles::fasta::fs::index` — no samtools/GATK), fetching on a
miss via streaming reqwest + `flate2` decompress. `ReferenceGateway` exposes
`reference_status`/`cached_reference`/`resolve_reference` (+ `resolve_chain`/`load_liftover`
over `du-bio`). `App` holds it; `import_project_dir`'s reference is now optional — on a cache
miss it returns `AppError::ReferenceNeeded`, and the UI prompts → downloads with a progress
bar → auto-retries. Registry defaults from `docs/chm13-reference-resources.md` (+ a
`reference_sources.json` override); cache under `~/.decodingus` (`$NAVIGATOR_REFGENOME_DIR`).

**Applying liftover — BUILT + validated for Y *and* mtDNA (2026-06-04):**
`assign_haplogroup_from_alignment` lifts each tree's positions onto the alignment's build,
queries the lifted coords, maps observed bases back to tree positions (scoring unchanged):
- **chrY**: GRCh38→build via the cached, auto-downloaded nuclear chain (`gateway.lift_positions`).
  Minus-strand lifts are reverse-complemented (large CHM13-Y tracts are inverted vs GRCh38).
- **chrM**: no chrM chain exists (and CHM13's `chrM` is a *circular permutation* of rCRS —
  origin at rCRS ~577), so we **self-generate** the map: `mtvariants::mt_position_map` detects
  the rotation, rotates into the rCRS frame, banded-aligns (bundled rCRS vs this reference's
  `chrM`), and composes the offset back. No chain, no download, no cache (~16.5 kb align).
Validated live on GFX0457637 (CHM13 HiFi): **Y = R-FGC29071** (1092/1919) and
**mtDNA = U5a1b1g** (53/55), both matching the GRCh38 truth. Related caches if/when ancestry
lands: `AncestryReferenceGateway`/`Cache`, `TreeCache`, `AnalysisCache`.
- **Resource catalog:** `docs/chm13-reference-resources.md` lists the concrete CHM13v2.0 URLs —
  reference FASTAs (incl. **`chm13v2.0_maskedY_rCRS.fa.gz`**, the most relevant for this app),
  GRCh38↔CHM13 and hg19↔CHM13 **1:1 liftover chains** + `unique_to_*` BEDs (unliftable regions),
  lifted variant catalogs (dbSNP155, ClinVar, gnomAD, 1000G/SGDP on CHM13), and accessibility/
  repeat/censat masks. All on the public `human-pangenomics` bucket (CC0).

---

## 5. GUI — egui / eframe

- **Single static binary** per OS, GPU-rendered, no system webview or JRE — the simplest
  robust cross-platform story for Win/Lin/macOS, which is the stated priority.
- **Architecture:** `navigator-ui` is thin. It renders immutable view-state and dispatches
  **commands** to `navigator-app`; long-running analysis runs on a worker thread (or
  `tokio` runtime) and streams progress back via a channel, which the egui repaint loop
  reads. No business logic, no DB calls, no domain decisions in the UI — directly fixing
  the "dialogs make domain decisions" debt (`FingerprintMatchDialog`, etc., become app-layer
  policy with a UI prompt only when truly needed).
- **What maps over:** the 37 ScalaFX dialogs/panels become egui panels/windows; data
  tables and coverage/haplogroup charts are an egui strength. `SubjectDetailView`
  (3,104 LOC) decomposes into per-tab widgets backed by app-layer queries.
- Immediate-mode is awkward for very large forms — accept that; the workbench is
  table/chart/report-heavy, which suits egui.

---

## 6. Tech-debt remediation — explicit mapping

| Today (Scala) | Rewrite fix |
|:---|:---|
| `WorkbenchViewModel` god object (4,021 LOC) | Split across `navigator-app` (commands), `navigator-store`, `navigator-analysis`; `navigator-ui` holds only view-state + dispatch. |
| `HaplogroupResult` / `ScoredHaplogroup` / `HaplogroupAssignments` triplication | One domain type in `du-domain`/`navigator-domain`; serialization variants via serde, not parallel types. |
| Slick 22-tuple → JSONB blobs | Proper SQLite tables; `Json<T>` only for AT Proto snapshots (§3). |
| `EntityConversions` (549 LOC) fragile mapping | sqlx `FromRow` + small `into_domain` per aggregate (the `du-db` pattern). |
| Mixed error handling (`Either`/`Try`/`Option`/exceptions; swallowed failures) | `thiserror` enums per layer (`DomainError`/`StoreError`/`AppError`), propagate with `?`; no proceed-after-failed-persist. |
| `AsyncSyncService` stubbed | `navigator-sync` completed: explicit retry/backoff, conflict policy, offline indicator; reuse `du-atproto` OAuth session. |
| Config fragmentation (3 sources, deprecated overrides) | One validated config (env > file > default). |
| Cache with no invalidation/versioning | Cache key includes **algorithm version**; explicit invalidation + hit metrics. |
| 7-table Y-profile schema + manual concordance | Simplify to a `YDnaProfile { snps, strs, sources, reconciliation }` aggregate. |

"Port as-is" (logic is sound, just translate): haplogroup **scoring** (tree DP traversal),
external clients (ENA, facility, tree providers), IBD detection/relationship math.
"Redesign": workspace state model, persistence, sync, UI, analysis-processor interface
(define one `Processor` trait with `init/run/cleanup` + standard progress signature).

---

## 7. AT Protocol / OAuth

Reuse `du-atproto` — but **partially**. The primitives reuse directly: DID/handle
resolution, PDS discovery, PKCE (S256), DPoP proofs. The confidential-client pieces it also
implements (`private_key_jwt` ES256 client assertion, served client-metadata/JWKS, cookie
session) are for the **decodingus web** client and do **not** apply here: Navigator is a
desktop app, hence a **public/native client** — PKCE only, its own native
`client-metadata.json`, a **loopback redirect** (`http://127.0.0.1:<port>/callback`), and
tokens in the **OS keychain**. Confirm `du-atproto`'s token-exchange builder runs
PKCE-without-client-assertion (small add if not). See
`documents/atmosphere/11-Auth-and-Permissions.md` §6–7 and the server-side companion
`decodingus/rust/docs/atproto-oauth-findings.md`.

This is the OAuth migration the design review called for: app-password `createSession` is
gone; Navigator authenticates via OAuth and requests the `navigatorCore` **write** scope.

**AppView scope reduction (2026-06) shifts work to Navigator.** The AppView no longer
mirrors the network (see `documents/atmosphere/08-AppView-Lifecycle.md`); per-sample data is
authoritative in the Navigator workspace. This adds two Navigator responsibilities:
- **Publish per-sample coverage summaries** as public PDS records (under `navigatorCore`)
  so the AppView can aggregate them on demand. Coverage reporting for the researcher's own
  cohort is also a Navigator-local feature over the local workspace.
- **A variant-proposal submission client** — Navigator posts variant/branch proposals
  directly to the AppView curation API (decoupled from PDS records), replacing the old
  firehose-harvest of private variants.

These live in `navigator-sync` (publishing) and a small `navigator-app` command (proposal
submission). The desktop app remains a **writer** (direct-to-PDS) plus a reader of its own
AppView records; no broad PDS read scope is needed at this stage.

---

## 8. Data migration (H2 → SQLite)

Follow the `du-migrate` precedent: a one-time CLI (`navigator-migrate`) reads the existing
`~/.decodingus/data/workspace.mv.db` (H2) and writes the new SQLite schema. Map
biosamples/projects/sequence-runs/alignments/haplogroups; re-derive anything cheap rather
than migrating it (analysis caches can be recomputed). Keep the H2 file untouched as a
rollback.

---

## 9. Phased roadmap

1. **Extract shared crates.** Promote `du-domain`/`du-atproto`/`du-bio` to the shared
   location; stand up the `navigator` workspace skeleton with the §2 crates compiling empty.
2. **Raw-read layer + finish du-bio.** Finish shared `du-bio` `liftover`/`callable`/`vcf`
   (coordinate math + text parsing). In **`navigator-analysis`**, implement the `noodles`
   BAM/CRAM/FASTA/BGZF/index I/O layer and port the **walkers** (§4a) on top of it, with
   unit tests.
3. **Caller + parity harness.** Build the haploid caller (§4b) and the golden-truth
   validation harness (§4c). Gate on parity. (Sidecar bridge §4d available if needed.)
4. **Store + app layer.** SQLite schema + migrations; `navigator-app` command/query layer;
   `navigator-migrate` from H2.
5. **UI.** egui shell → dashboard, subject detail, haplogroup/coverage/ancestry views,
   reports. Replace dialogs with app-layer policy + thin prompts.
6. **Sync + OAuth.** Wire `du-atproto` OAuth; complete `navigator-sync` (the part that is
   stubbed today).
7. **IBD + Y-profile.** Port detection/relationship math; simplified Y-profile aggregate.
   IBD detection/relationship is done. The simplified Y/mtDNA **multi-source reconciliation**
   aggregate (combine runs across platforms + Sanger confirmation, with identity
   verification) is specified in `MultiSource_Reconciliation.md` — built (all 6 phases).
8. **Reference & liftover management (§4f).** Retrieval + on-disk cache **built** as
   `navigator-refgenome` (2026-06-03): resolve a build → cached, in-Rust-indexed `.fa`
   (fetch/decompress/index on miss); chains cached for `du-bio`; import resolves references
   from the cache with a download prompt. Applying liftover is **built + validated for Y and
   mtDNA** on CHM13: Y via the nuclear chain (reverse-complementing inverted tracts), mtDNA via
   a self-generated rotation-aware rCRS↔`chrM` map. GFX0457637 → R-FGC29071 + U5a1b1g, matching
   GRCh38.
9. **Ancestry / PCA.** Phase 1 **built** (2026-06-04): allele-frequency likelihood at super-
   population granularity (AFR/AMR/EAS/EUR/SAS). No GATK — `caller::genotype_sites` (diploid GL
   model) genotypes the sample at an AIMs panel, `analysis::ancestry::estimate_by_allele_frequency`
   scores each population by the binomial likelihood of the diploid genotypes. The panel
   (per-population alt-allele frequencies at high-Fst sites) is built offline by the new
   `navigator-panelbuild` binary from the 1000G-on-CHM13 VCFs (per-super-pop counts come straight
   from INFO `AC_<POP>_unrel`/`AN_<POP>_unrel`); selected by Nei Fst (`--max-sites`/`--min-fst`).
   `App::estimate_ancestry` loads the build-matched panel (`$NAVIGATOR_ANCESTRY_PANEL` or
   `<cache>/ancestry/`), genotypes, estimates, persists (`ancestry_result` table); UI shows a
   super-population proportion bar. **Phase 2 — PCA — built + validated (2026-06-04):**
   `navigator-panelbuild pca` builds `PcaLoadings` (per-SNP loadings+means, per-pop
   centroids+variances) from the 1KGP genotype matrix via a sample-space Gram eigendecomposition
   (nalgebra, pure Rust). The genotype matrix is extracted with `bcftools query -R <panel sites>`
   against the remote `unrelated_samples_2504` VCFs (run on an in-AWS EC2 instance — the ~1 TB
   files are tabix-indexed, so only the ~20k panel sites are pulled). `analysis::project_pca`
   (+`classify_pca`) projects a sample; `App::estimate_ancestry` fills `pca_coordinates`, and the
   UI draws a PC1×PC2 scatter (reference centroids + sample). Validated on GFX0457637: 18,634
   sites × 10 PCs, classic 1000G structure (PC1=AFR axis, PC2=EAS, PC3=SAS), sample projects
   onto the EUR centroid. **Fine-grained (26-population) — built + validated (2026-06-04):**
   `navigator-panelbuild fine-panel` computes per-fine-population AF from the genotype matrix, and
   `pca` keys centroids on the 26 fine populations; the domain carries the fine catalog +
   fine→super mapping, and the estimator rolls fine components up into the super-population
   summary. The UI leads with the (robust) super-population rollup and lists the top fine
   populations; the scatter shows all 26 centroids (colored by continent). Validated on
   GFX0457637: European 100% (rolled up), nearest fine centroid IBS — within-noise, as expected
   for a continental-AIMs panel.
10. **Cutover.** Feature-parity check against the golden harness; ship.

---

## 10. Risks & open questions

**Risks**
- **Indel/homopolymer calling parity** (§4b) — the single biggest technical risk; owned by
  the §4c harness. May constrain v1 branch-creation to SNP evidence.
- **Numerical drift** in metrics/coverage vs GATK — covered by the same harness.
- **Shared-crate coordination** — extraction adds a release axis between two teams/repos.
- **egui polish ceiling** — acceptable for a workbench; revisit only if a report view needs
  rich native widgets.

**Open questions**
1. ~~Where do the shared crates live?~~ **Resolved:** `decodingus-shared`, a sibling repo at
   the Development root (DUNavigator/decodingus/decodingus-shared as siblings). Remaining
   sub-question: path deps vs git deps vs published registry for CI.
2. Does Navigator share the **same SQLite schema shape** as the decodingus Postgres schema
   where domains overlap (biosample/variant/haplogroup), or only the domain *types*?
3. Long-read specifics (PacBio/ONT) the walkers must preserve (e.g. `COUNT_UNPAIRED`,
   read-length handling) — enumerate so parity tests cover them.
4. mtDNA indel fallback: is a pinned external caller acceptable *temporarily* for
   validation, or must everything be pure-Rust from day one?
