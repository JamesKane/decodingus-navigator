# Ancestry fine-population wiring + IBD reference artifacts — design

Status: **design / specification only** (no code). Branch: `rust-rewrite`.
Scope if built: `navigator-panelbuild` (new asset products), `navigator-analysis`
(fine-admixture + genetic-map loader), `navigator-app` (asset loading + IBD wiring),
`navigator-ui` (fine-pop display + asset transparency).

## Problem

The `navigator-panelbuild` pipeline now produces a rich asset set in `~/.decodingus/ancestry/`
(1000G fine pops + SGDP + AADR ancient deep components on CHM13), but the app consumes only a
fraction of it, and the IBD modules run on stand-in inputs.

What the app actually loads today (verified):
- `ancestry_panel_<build>.bin` — **required**; drives the 5-way super-pop admixture donut
  (`navigator-app/src/lib.rs:545,3655`; `estimate_admixture`).
- `ancestry_pca_<build>.bin` (+ optional `ancestry_pca_ancient_<build>.bin`) — **optional**; drives
  PCA-GMM / nMonte and the PC1×PC2 scatter (`lib.rs:557,3676`).
- **Not referenced anywhere:** `ancestry_freq_global_<build>.bin` (the fine 173-population AF),
  the `ancestry_panel_{fine,super}` / `ancestry_pca_{fine,super}` variants, and
  `ancestry_manifest_<build>.json`.

Consequences:
1. The admixture donut is only 5 super-populations even though we built fine (26 1000G) + SGDP +
   ancient frequencies. The fine/SGDP/ancient pops appear *only* as centroids in the PCA scatter,
   not as admixture proportions.
2. IBD (`navigator-analysis/src/ibd.rs`, `App::compare_ibd` at `lib.rs:5217`) reuses the **sparse,
   Fst-selected AIM panel** as its SNP set and a **hardcoded uniform genetic map**
   (`GeneticMap::uniform(1.0, …)` at `lib.rs:5249`) — both suboptimal for IBD.
3. Published asset integrity isn't verified (the manifest sha256s are written but never checked on
   load).

This doc specifies (A) wiring fine-population ancestry, and (B) the IBD reference artifacts.

## Goals

- Surface fine-population + ancient ancestry as **proportions**, not just scatter centroids.
- Give IBD a real genetic map and a fit-for-purpose SNP panel, replacing the stand-ins.
- Make the IBD panel **consumer-chip compatible** (23andMe / AncestryDNA / MyHeritage / FTDNA /
  LivingDNA) through an **allele-aware, multi-build liftover** — WGS is the target, but chip kits
  outnumber WGS by orders of magnitude, so chip↔chip and chip↔WGS matching is the volume case and
  the binding design constraint.
- Keep the offline build (`navigator-panelbuild`) the single source of these assets; keep the app a
  consumer with graceful degradation when an asset is absent.

## Non-goals

- Re-deriving reference panels from raw alignments (the build already does this off public call sets).
- Phased / IBD2 detection, chromosome-painting accuracy work, the federated match-discovery UI
  (separate tracks).
- Changing the CHM13-only build-match gate (`lib.rs:3659`) — see cross-cutting note.

---

## Part A — Fine-population ancestry

### A1. The asset is already built; the gap is consumption

`ancestry_freq_global_<build>.bin` holds per-fine-population allele frequencies **at the same panel
sites** the app already genotypes for admixture (it's produced by `panelbuild fine-panel` from the
same matrices). So **no extra genotyping is needed** — the existing `genotype_panel()` dosages over
`ancestry_panel`'s sites already cover it (the two share the site set).

### A2. Loader + path resolution

Mirror the existing `ancestry_panel_path` / `ancestry_pca_path` helpers:
`ancestry_freq_global_path(build)` → `$NAVIGATOR_ANCESTRY_FREQ` override, else
`<base>/ancestry/ancestry_freq_global_<build>.bin`. Optional, like the PCA asset — absent ⇒ fine
admixture silently skipped.

### A3. Estimator — and the 173-population problem

`estimate_admixture` is a supervised EM over reference allele frequencies. A naïve 173-way EM is
ill-posed (many collinear modern pops + sparse pseudo-haploid ancient refs), so **don't** run one
flat admixture over every population in `freq_global`. Recommended structure:

| Layer | Reference set | Method | Asset |
|-------|---------------|--------|-------|
| Super-pop (today) | 5 super-pops | supervised EM | `ancestry_panel` |
| **Fine modern (new)** | 26 1000G + curated SGDP regions | supervised EM, restricted | `ancestry_freq_global` (modern subset) |
| **Ancient deep (new framing)** | Steppe/ANF/WHG/EHG/Iran_N | distance / nMonte | `ancestry_pca` (projection coords) |

- **Fine modern admixture:** run the same EM as `estimate_admixture`, but over a *curated modern
  subset* of `freq_global` (the 26 1000G fine pops, optionally a small set of well-sampled SGDP
  regional pops). The asset carries all 173; the estimator selects the modern reference subset (a
  config/embedded list) to keep the EM well-conditioned.
- **Ancient deep components:** the ancient pops were built in *projection mode* and belong with the
  distance estimators (nMonte / PCA-GMM), which already run off `ancestry_pca`. Surface ancient as a
  separate "deep ancestry" estimate, not mixed into the modern EM. (`freq_global`'s ancient columns
  are mainly useful for an optional supervised "which deep source" pass later.)

### A4. Data model + persistence

`AncestryResult` already carries a `method` field and is persisted/keyed per method
(`estimate_ancestry` stores ADMIXTURE / PCA_PROJECTION_GMM / G25_NMONTE separately). Add a
`FINE_ADMIXTURE` method record; the UI hierarchy grid already renders fine sub-populations under
each super-pop when present (`ui.rs:4917-4950`), so the display path largely exists — it just needs
the fine result populated and attached.

### A5. UI

- Populate the existing super→fine hierarchy from the new fine result (donut stays super-pop;
  expandable rows show fine modern pops > 0.5%).
- Add a small **asset-transparency** affordance (tooltip / "data sources" line): which assets are
  loaded and at what resolution — e.g. "super-pop panel + fine frequencies (173 ref pops) + PCA
  (modern basis, ancient projected)". Today nothing tells the user whether PCA/fine assets are
  present; absence is silent.

### A6. Cross-cutting: the CHM13-only gate

Ancestry requires `alignment.reference_build == panel.build` (`chm13v2.0`) — a vendor GRCh37/38 BAM
can't run *any* of this. That's the direct tie-in to **`realignment-module.md`**: realignment to
hs1 is the prerequisite that lets off-build samples reach these panels. Worth a one-line pointer in
both docs.

---

## Part B — IBD reference artifacts

IBD detection is pure pairwise IBS over diploid dosages; it needs two reference inputs the build
doesn't yet provide well: a SNP set and a genetic map.

### B1. Genetic map (the bigger correctness win)

Today: `GeneticMap::uniform(1.0, …)` — a flat 1 cM/Mb. Segment lengths in cM (and therefore the
relationship thresholds in `RelationshipEstimate`) are only as good as the map, so a real
recombination map matters more than panel density.

- **Source:** a sex-averaged genetic map (e.g. deCODE 2019, or HapMap II) on GRCh38, **lifted to
  CHM13** — coordinates only, no alleles, so the allele-aware concern from the AADR work doesn't
  apply; the stage-2 BED `CrossMap` lift is sufficient (interpolate cM at lifted positions).
- **Asset:** `genetic_map_<build>.bin` — per-chromosome `(bp, cM)` arrays. New `panelbuild`
  subcommand (`genetic-map`) or a step that fetches + lifts + serializes it.
- **Loader:** `genetic_map_path(build)` + `GeneticMap::from_bytes`; `compare_ibd` loads it and
  passes it to `detect_segments` instead of `GeneticMap::uniform`. Absent ⇒ fall back to uniform
  with a logged warning (don't hard-fail).

### B2. The IBD panel must be consumer-chip compatible (the volume case)

The 20k Fst-selected AIM panel is the wrong shape for IBD twice over: it's sparse and
ancestry-*biased* (bad for neutral IBS), **and** it isn't aligned to what consumer arrays assay.
WGS is the target, but chip kits dominate by orders of magnitude, so the panel must be the set
where **chip and WGS overlap** — otherwise chip-only samples (the majority) can't match anything.
WGS covers everything; the binding constraint is the arrays' content.

**Composition** — `ibd_panel_<build>.bin`, built by a new `panelbuild ibd-panel` step:
- Intersect the **common backbone of the major consumer arrays** (the OmniExpress / Illumina-GSA
  cores the vendors share) with common 1000G-CHM13 SNPs (biallelic, MAF ≥ ~5%). Denser than the
  AIMs (hundreds of k) and, crucially, chip-assayed.
- **Exclude strand-ambiguous palindromes (A/T, C/G)** — see B2b.

### B2a. Multi-build, allele-aware coordinate map (the key new structure)

Chip raw data arrives on **GRCh37** (23andMe v5, AncestryDNA) or **GRCh38** (newer kits), on the
array's design strand. So the CHM13 IBD panel must be **multi-build**: each site carries its rsID
and its `(contig, pos, REF, ALT)` on **CHM13, GRCh37, and GRCh38**, so a chip genotype on any build
resolves to the right CHM13 site and orientation. This generalizes the **BISDNA build-keyed
coordinate dictionary already in the codebase** (SNP→locus per build, CHM13 via liftover).

Building the cross-build coordinates uses the **allele-aware liftover the population-panel work
established — GATK `LiftoverVcf`, not CrossMap.** This is the same reverse-complement lesson from
the AADR build: CrossMap blanked the ALT wherever the reference base differed, silently corrupting
~3/4 of sites; GATK reverse-complements on inverted/minus-strand chain blocks and swaps REF/ALT.
A missed strand flip here corrupts IBD dosages exactly the same way it corrupted the ancient PCA.

### B2b. Chip strand reconciliation + the palindrome trap

A chip reports two observed alleles per SNP (e.g. `AG`) on the array's strand, which may be + or −
relative to the panel's REF/ALT. Reconcile to the CHM13 panel orientation by comparing the chip's
allele pair against `{REF, ALT}` and against `{rc(REF), rc(ALT)}`:
- chip alleles ⊆ `{REF, ALT}` → same strand; dosage = count of ALT.
- chip alleles ⊆ `{rc(REF), rc(ALT)}` → opposite strand → **reverse-complement**, then count ALT.
- neither → drop (genotyping error / multiallelic / wrong site).

This is the `strand_reconcile` logic already used in chip haplogroup import, generalized to the IBD
panel. **Palindromic SNPs (A/T, C/G) cannot be resolved by allele comparison** — `rc(A)=T` is also a
valid allele, so strand is ambiguous. These are the classic strand-flip trap, so **exclude A/T and
C/G sites from the chip-compatible panel** (standard for cross-array IBD/imputation). MAF-based
resolution is possible but error-prone and not worth it for IBD.

### B2c. Genotype sourcing per data type (uniform output)

- **WGS / CRAM:** genotype the IBD panel sites from the BAM (caller), cached per alignment.
- **Chip:** the imported `ChipProfile` genotypes → mapped through the multi-build map → strand
  reconcile (B2b) → dosage at the CHM13 panel sites. **No alignment needed** — works for chip-only
  samples, which is the whole point.

Both paths emit the *same* dosage vector over the CHM13 IBD panel, so the comparison stays
data-type-agnostic. (For chip↔chip, IBD is limited to the panel sites both kits assayed; for
chip↔WGS, to the chip's sites. Record the effective overlapping-site count per comparison — sparse
overlap weakens short-segment calls and must be surfaced, not hidden.)

### B3. App wiring

`compare_ibd` (`lib.rs:5217`) changes to: load the `ibd_panel` + each sample's IBD-panel dosages
(from CRAM **or** chip) + the real `genetic_map`, then `detect_segments`. The dosage *source*
differs by sample type; the comparison is uniform. Identity verification keeps its
autosomal-concordance + optional Y-STR path (no new asset). Surface the per-comparison overlapping
site count and (for chip pairs) which array generations were involved.

### B4. Status note

The IBD math (`PairwiseIbdDetector`, `GeneticMap`, `RelationshipEstimate`) is built and
`compare_ibd` is wired, but the UI uses it to *validate federated suggestions*
(`network_suggestions_section`), not as a standalone local compare. This asset work is independent
of — and a prerequisite for — making local IBD trustworthy whether or not the match-discovery UI
lands.

---

## Cross-cutting: manifest verification

`ancestry_manifest_<build>.json` records per-asset sha256 and is published for clients to verify,
but the app never reads it. When the asset loaders are touched, have them (optionally) check the
loaded `.bin` against the manifest and refuse a checksum mismatch — cheap integrity guard for CDN-
delivered assets.

## Validation plan

- **Fine ancestry:** on the validated GFX donor (EUR), the fine modern admixture should resolve to
  plausible European fine pops (e.g. NW-European 1000G pops) and sum to ~100%; the super-pop donut
  must be unchanged. Ancient deep estimate should not contaminate the modern EM.
- **Genetic map:** segment cM totals for a known relative pair (or the same sample vs itself →
  ~full genome IBD) should match expectation far better than the 1 cM/Mb stand-in; spot-check
  against a published map's chromosome cM lengths.
- **IBD panel:** parent/child and sibling fixtures land in the correct `RelationshipEstimate`
  bands; self-vs-self ≈ fully IBD.
- **Chip compatibility (the critical cross-check):** the *same individual* genotyped by both a chip
  (23andMe/AncestryDNA) and WGS must come out as self/identical IBD — this is the end-to-end proof
  the multi-build map + strand reconciliation are correct. A strand-flip bug would instead make the
  chip and WGS of one person look *unrelated* (dosages anti-correlated), the same silent-corruption
  signature as the AADR PCA. Also verify a known chip↔chip relative pair, and confirm A/T,C/G sites
  were excluded (no palindromes in the panel).
- **Graceful degradation:** with each new asset absent, the app degrades (super-pop only / uniform
  map) without error.

## Phasing / rollout

1. **Fine-admixture wiring** (asset already exists): loader + modern-subset EM + persist
   `FINE_ADMIXTURE` + populate the existing hierarchy UI. Smallest, highest-visibility win.
2. **Asset transparency + manifest check** in the loaders.
3. **Genetic map**: `panelbuild genetic-map` (fetch + coordinate-only lift + serialize) + loader;
   swap `compare_ibd` off `GeneticMap::uniform`.
4. **Chip-compatible IBD panel**: `panelbuild ibd-panel` building the **multi-build, palindrome-free
   panel** (CHM13 + GRCh37 + GRCh38 coords via GATK allele-aware liftover); WGS genotyping at its
   sites cached in deep analysis; point `compare_ibd` at it.
5. **Chip → IBD dosage path**: map `ChipProfile` genotypes through the multi-build map + strand
   reconcile (reuse/generalize the chip-import `strand_reconcile`), so chip-only and chip↔WGS
   comparisons work. This is the volume unlock.
6. (Later) ancient supervised "deep source" pass from `freq_global`; AF-aware / IBD2 / phased IBD.

Phase 1 surfaces the richness already built; phases 3–5 make IBD quantitatively trustworthy and —
critically — usable for the chip-majority population, not just WGS.

## Open questions

- Which modern fine reference subset for the EM — 26 1000G only, or 1000G + a curated SGDP regional
  set? (Trade-off: resolution vs EM conditioning / collinearity.)
- Genetic-map source: deCODE 2019 vs HapMap II vs a CHM13-native map if one is available — and
  sex-averaged vs sex-specific.
- **IBD panel array backbone:** which vendor generations to intersect (23andMe v4/v5, AncestryDNA
  v1/v2, MyHeritage, FTDNA, LivingDNA — they diverged toward GSA over time, so the common core
  shrinks with newer kits). Target density vs cross-vendor overlap is the core trade-off.
- **Palindrome handling:** exclude A/T,C/G outright (recommended) vs MAF-resolve a subset to recover
  density — and the minimum overlapping-site count below which a chip↔chip IBD call is suppressed.
- Whether to fold IBD-panel genotyping into the standard deep-analysis pass (cost per WGS sample)
  or keep it on-demand.
- Should the multi-build IBD panel and the **BISDNA build-keyed coordinate dictionary** be unified
  into one shared cross-build SNP-coordinate facility (both need rsID→per-build locus + strand)?
- Should `ancestry_freq_global` ever feed a flat fine admixture, or is the super-pop EM +
  distance-based ancient the right division of labor permanently?
