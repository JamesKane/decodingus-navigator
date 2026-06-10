# Rust rewrite — handoff / resume notes

Last updated: 2026-06-07. Branch: `rust-rewrite`. Pick up here next session.

## Three-repo topology

- **DUNavigator** (this repo, `rust-rewrite`) — the desktop edge app (egui). Pulls the shared
  crates from the sibling repo by **path**, so working-tree changes are picked up immediately.
- **decodingus-shared** (`/Development/decodingus-shared`, branch `feat/fed-report-records`,
  pushed) — `du-domain` / `du-atproto` / `du-bio`. The federated record contracts live in
  `du-domain::fed`.
- **decodingus** (`/Development/decodingus`, the AppView, branch `rust-rewrite-foundation`) — the
  PostgreSQL hub + web (`rust/` is the Rust rewrite of the Play app). Pulls the shared crates by
  pinned **git rev** (`du-domain` currently pinned at `f975a08`). **Bump the rev (or add a local
  `[patch]`) to pick up newer shared changes** — e.g. the `fitDistance` field added after f975a08
  is NOT yet visible to the AppView.

## State

Workspace builds clean (no warnings); `cargo test` green. Major work since 2026-06-04, by theme
(newest first; all on `rust-rewrite` unless noted):

- **i18n** — `0bb19e9` scaffolding (Lang En/Es, `key=value` catalogs in `crates/navigator-ui/
  locales/`, `tr()` with active→En→key fallback, app-bar language switcher); `0cf412b` migrated
  card titles, dashboard, empty states, primary buttons. Deeper forms/headers/status still English
  (fall back fine). See `memory/navigator-i18n.md`.
- **Subject-centric model** (`documents/design/SubjectCentricModel.md`, P1–P3 done) — `ba6ffc2`
  auto-select default alignment; `9259f7a` donor STR consensus + ancestry provenance; `d928af8`
  donor ancestry (best-of) + private-Y union. Tabs now present the **donor**, per-run detail in
  Data Sources. Remaining: true genotype-level ancestry pooling.
- **Analysis QoL / fixes** — `463e938` per-contig de-novo (Y-DNA=chrY / mtDNA=chrM), mtDNA
  haplogroup in full analysis, inferred-sex write-back, table Y/mt/sex; `4d54a8d` standalone
  "Assign mtDNA haplogroup"; `046b654` dropped the UNMASKED raw chrY de-novo (use "Find private Y
  variants" instead) + Y card shows persisted consensus; `f258636` full analysis reuses cached
  coverage+ancestry, private-Y persisted.
- **Header auto-detect** — `9a2062b` `navigator-analysis::probe` reads BAM/CRAM headers → build /
  aligner / platform / test-type; `add_data` auto-imports a BAM (creates run+alignment, no
  questions); reference resolved from the build via the gateway (never asked).
- **Full Analysis** — `e94f8ec` modal + `RunFullAnalysis` pipeline (8 steps, cancellable);
  `dc48f98` per-contig coverage progress + alive modal (it was slow, not hung).
- **UI Workbench redesign** — `c64cd50` dark theme + nav tabs + subjects table + detail sub-tabs;
  `133dfdf` Data Sources cards; later all detail tabs carded. See `memory/ui-workbench-redesign.md`.
- **Ancestry methods** — `7a67c22` `estimate_pca_gmm` + `estimate_nmonte` (Frank-Wolfe distance
  fit, `fit_distance`); `method` field captured not inferred; estimate_ancestry computes+persists+
  publishes 3 methods (ADMIXTURE / PCA_PROJECTION_GMM / G25_NMONTE). `e5b502f` projection-mode PCA
  (`--basis-pops`). See `memory/ancestry-pca-gmm.md` + `documents/design/AncestryAnalysis.md`.
- **Federated report records** — `du-domain::fed` (shared) defines the atproto wire contracts
  (WireF64 strings + numeric storage projections); Navigator publishes alignment/biosample/
  sequencerun/populationBreakdown (`7a67c22`/`dca5b5c`); AppView `du-jobs` ingests (`83bfc0f`).
  Shared: `f975a08` + `0aa960c` (fitDistance). See `memory/fed-report-records.md`.
- **Global ancestry panel pipeline** — `scripts/ancestry-panel/` (`67a5185`): fetch AADR/HGDP/
  SGDP/1000G-CHM13 → liftover → 1240k-restricted AIMs → matrices → PCA+freq assets → CDN. Decoupled
  from raw alignments (needs genotypes, not reads). URLs/scale verified (see
  `memory/ancestry-panel-pipeline.md`): a naive pull is ~5 TB — slice panel sites, don't bulk-DL.

Uncommitted of note: none. Untracked: `CLAUDE.md`, `GEMINI.md` (leave). A stray
`crates/.claude/settings.local.json` is recreated by the environment and would break the cargo
glob — handled by `exclude = ["crates/.claude"]` in the root `Cargo.toml`.

## Build / validate

```bash
cargo build && cargo test --workspace
cargo run -p navigator-ui            # desktop app
```

### Live (`#[ignore]`) tests — real data, run explicitly
Test sample: `/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam` (CHM13 HiFi, ~4×,
male, Y=R-FGC29071, mtDNA=U5a1b1g, European). Reference: `/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa`.

```bash
GFX_CHM13_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam \
GFX_CHM13_REF=/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa \
NAVIGATOR_ANCESTRY_PANEL=/Users/jkane/.decodingus/ancestry/ancestry_panel_chm13v2.0.bin \
NAVIGATOR_ANCESTRY_PCA=/Users/jkane/.decodingus/ancestry/ancestry_pca_chm13v2.0.bin \
  cargo test -p navigator-app --release \
  validate_gfx_chm13_ancestry local_ancestry_paints_gfx gfx_sex_is_male -- --ignored --nocapture
```
Expected: European ~98% (admixture), DNA painting EUR-dominant, sex=Male. Other ignored live tests:
`validate_gfx_chm13_haplogroups` (Y/mt), parity_real.rs (HG002 env), PDS publish (PDS_TEST_URL).
`NAVIGATOR_ANCESTRY_PCA_ANCIENT` points the PCA-GMM at an ancient-component asset when present.

## Ancestry assets (regenerable; not committed)

Installed at `~/.decodingus/ancestry/`:
- `ancestry_panel_chm13v2.0.bin` — AF panel (genotyping + admixture; the default the app loads)
- `ancestry_pca_chm13v2.0.bin` — PCA loadings + per-pop centroids (drives PCA-GMM + nMonte)

Today's assets come from the archived genotype matrix `~/Genomics/archive/1kgp_chm13_pca_build/`
(`gt_all.tsv.gz` 1000G + `sgdp_gt.tsv.gz` + `combined_pops.txt`):
```bash
A=~/Genomics/archive/1kgp_chm13_pca_build; O=~/.decodingus/ancestry
navigator-panelbuild fine-panel --matrix $A/gt_all.tsv.gz,$A/sgdp_gt.tsv.gz \
  --samples $A/samples.txt,$A/sgdp_subset_samples.txt --pops $A/combined_pops.txt --out $O/ancestry_panel_chm13v2.0.bin
navigator-panelbuild pca        --matrix $A/gt_all.tsv.gz,$A/sgdp_gt.tsv.gz \
  --samples $A/samples.txt,$A/sgdp_subset_samples.txt --pops $A/combined_pops.txt [--basis-pops modern.txt] --out $O/ancestry_pca_chm13v2.0.bin
```
The **next-gen** asset path is the global pipeline in `scripts/ancestry-panel/` (modern + ancient
deep components over a 1240k-restricted panel, projection-mode PCA, CDN publish) — needs the
datasets fetched (verify `# VERIFY` URLs; slice panel sites to avoid the multi-TB pull).

## EC2 (genotype extraction)

`admin@ec2-3-132-31-28.us-east-2.compute.amazonaws.com`, key `~/Decoding-Us.pem` (chmod 600),
Debian, bcftools. Used to tabix-fetch panel-site genotypes from the ~1 TB 1000G/SGDP VCFs (in-AWS
S3 is fast). Recipe + region files archived in `1kgp_chm13_pca_build/ancestry_build.tar.gz`. The
matrices are already pulled; re-extraction only needed to add reference samples.

## Key gotchas (also in agent memory)

- CHM13 `chrM` is a circular permutation of rCRS → rotation-aware self-generated map; CHM13 Y has
  inverted tracts → liftover reverse-complements minus-strand lifts.
- **Raw chrY de-novo is unmasked** — calls ~13k mostly-artifact variants. The validated Y discovery
  is **"Find private Y variants"** (callable-mask + backbone-classified). chrM de-novo is fine (small,
  fully callable). Don't re-add a raw chrY de-novo button.
- Full Analysis **reuses cached** coverage + ancestry (the slow whole-genome steps); coverage is a
  single-threaded full-genome pileup (minutes on WGS — slow, not hung; per-contig progress shows it).
- PCA projection of a low-coverage sample is rescaled by `total/used` sites (else it shrinks toward
  origin). SV needs ≥10× — won't run on 4× GFX. AIMs were super-pop-Fst selected → fine resolution noisy.
- **i18n borrow gotcha**: `TextEdit::singleline(&mut self.x).hint_text(self.tr(k))` fails — bind
  `let hint = self.tr(k);` first. `tr()` returns `&'static`.
- **AppView pins the shared crate by git rev** — bump it (or `[patch]`) for new `du-domain::fed`
  fields; additive optional fields keep the old rev compiling.

## Architecture / design pointers

- `documents/design/SubjectCentricModel.md` — donor-centric tab model (P1–P3 implemented).
- `documents/design/AncestryAnalysis.md` — the 3 estimators + ancient-asset build + nMonte/G25.
- `documents/atmosphere/` — the lexicon spec; `du-domain::fed` is the implemented write subset.
- `documents/BACKLOG.md` — **Scala-era** feature inventory (March 2026, pre-rewrite); use as the
  master feature list, not current status.
- Agent memory (`~/.claude/projects/.../memory/`) is the most current running record.

## Remaining gaps (from the 2026-06-07 audit)

**Unported from Scala**: i18n (scaffolded, partial migration), feature toggles/config, the
DecodingUs haplogroup tree provider (FTDNA-only today), haplogroup/comparison report export (CSV
only), full Y-STR concordance subsystem.
**Designed, not built**: IBD **network** matching (local detector + UI tab exist; consent/match-
discovery/chromosome-browser/relationship records do not — the biggest remaining feature), Genome
Regions API, sequence-run fingerprinting, academic-ENA / FTDNA-project / pangenome-GAM imports,
imputation, instrument-observation + ancestral-STR records, interactive tree viz.

**Done since this audit**: **Unified Quality Metrics Walker** (2026-06-10, committed 842b2bb) —
fused coverage + read-metrics + sex into one record-loop pass (`navigator-analysis::unified`); BAM
2→1 / CRAM 3→1 file reads; full-analysis steps 8→6. Shared `pub(crate)` `*State` accept/finish
helpers in coverage/read_metrics/sex → byte-identical numbers (standalone fns are now wrappers).
**Threaded** (committed 900381d + 7a9d7f2): MT bgzf decompression (~8% — decode isn't the
bottleneck), then a **per-contig parallel walker** (rayon) — **5.15× on the GFX BAM (64.7s→12.6s),
peak 2.8 GB, byte-identical** to sequential. Reference N-mask (1 bit/base) + a load semaphore bound
memory; `NAVIGATOR_ANALYSIS_THREADS` / `NAVIGATOR_BGZF_THREADS` tune it. BAM+indexed uses the
parallel path; CRAM / unindexed BAM fall back to sequential. Live parity + perf-smoke tests in
`tests/parity_real.rs`. See `documents/design/UnifiedQualityMetricsWalker_RustPort.md` +
`memory/uqmw-rust-port.md`.
**AppView side**: pick which ancestry method to surface; ingest `fitDistance` (needs rev bump);
IBD-matching AppView backlog; AppView backfeed.
**Small**: Edit/Delete + Add-to-Project UI stubs; Compare needs multi-select; table Y/mt only fills
the selected subject; ancestry genotype-pooling deferred; global panel asset needs data.
**Perf backlog** (unified walker, 2026-06-10): the per-contig parallel walker plateaus at ~5×
(knee ~12 threads) because the serial **unmapped-tail sweep** and the **single largest contig**
(chr1) floor the wall time. To push further: split big contigs into sub-regions with a per-region
coverage merge, and/or parallelize the unmapped sweep. Lower priority than functional gaps. See the
"Further headroom" note in `documents/design/UnifiedQualityMetricsWalker_RustPort.md`.

## Recommended next steps (pick one)

1. ~~**Unified Quality Metrics Walker**~~ — DONE 2026-06-10 (committed, incl. per-contig
   parallelism — 5.15× on the GFX BAM); see "Done since this audit" above.
2. **i18n** — keep migrating (forms, grid headers, status messages) and/or persist the chosen lang.
3. **DecodingUs tree provider** — you control that tree; FTDNA-only is a real limitation for Y work.
4. **IBD network matching** — biggest user-facing feature; detection + identity math built, the
   consent/discovery/browser UI + records are not (needs AppView work too).
5. **Edit/Delete + Add-to-Project** — small but visible UI stubs that need backend commands.
