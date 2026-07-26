# Archaic Ancestry Report (Neanderthal / Denisovan) — Design

**Status:** Design + **implementation plan (§10)**. No code yet. Drafted 2026-07-23; plan added
2026-07-26 on branch `feat/archaic-ancestry`. Two of the three §9 open questions are now settled;
the third is a gating checkpoint inside M1.
**Goal:** Reconstruct a 23andMe-style Neanderthal report — and go beyond it with a Denisovan
estimate and a true whole-genome introgression map — from public archaic reference genomes and
recent methods, using the app's existing ancestry/panel/HMM machinery.

---

## 1. What we're building, and why it splits in two

23andMe's report (see the reference UI: a headline **variant count**, a **percentile vs other
customers**, trait associations, an educational section) is deliberately *not* a "percent
Neanderthal." It is a **marker count**: 23andMe assays 3,731 (v5) / 1,436 (pre-v5) of the 135,171
Sankararaman-2014 Neanderthal-derived SNPs and reports how many archaic-allele copies you carry,
because a sparse chip "captures only a fraction of the total Neanderthal ancestry an individual
carries" (their words). That limitation is fundamental to chip data and dictates the architecture.

We ingest **both** consumer chip data *and* whole-genome BAM/CRAM/VCF. So the feature is two tiers:

- **Tier A — Marker count (all inputs).** A precomputed archaic-informative SNP panel; count the
  subject's archaic-derived allele copies; report a count + population percentile. This is the
  *only* scientifically honest option for ~600–700k-SNP chip data, and it is the 23andMe-equivalent.
  It also runs on WGS as a fast summary.
- **Tier B — Segment HMM (WGS/VCF only).** An hmmix-style archaic-segment caller producing an actual
  **introgressed Mb / % of genome**, a per-segment **Neanderthal-vs-Denisovan** attribution, and a
  chromosome browser. Requires dense genome-wide variants a chip cannot provide, so it is gated to
  WGS/VCF inputs.

This mirrors how the codebase already gates deep ancestry: `estimate_qpadm_ancestry` returns `None`
when the model doesn't apply rather than fabricating numbers. We adopt the same discipline for
Denisovan-in-Europeans (see §7).

---

## 2. Data availability (verified) — what we can actually download

All four archaic genomes are **openly downloadable, no access application**, from Max Planck EVA:

| Genome | Path | Coverage | Build |
|---|---|---|---|
| Altai Neanderthal (Prüfer 2014) | `cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/` | ~52× | hg19/GRCh37 |
| Vindija 33.19 (Prüfer 2017) | `cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/` | ~30× | hg19/GRCh37 |
| Chagyrskaya 8 (Mafessoni 2020) | `ftp.eva.mpg.de/neandertal/Chagyrskaya/` (VCF/, FilterBed/) | ~28× | hg19/GRCh37 |
| Denisova 3 (Meyer 2012) | `cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/` | ~30× | hg19/GRCh37 |

Per-chromosome genotype VCFs + BAMs + BED quality masks (`FilterBed/`). **All hg19/GRCh37** — the
Vindija README says so verbatim; Altai filenames carry the `hg19_1000g` token.

**Licensing:** governed only by the **Ft. Lauderdale principles** (a norm reserving *first*
genome-wide analysis to the producers), *not* a formal open-source license. Because Prüfer 2014/2017,
Meyer 2012, and Mafessoni 2020 are all long-published, the first-analysis reservation is moot.
Practically redistributable, but there is **no explicit CC/open grant** — same caution we logged for
AADR ([[germanic-panel-sources]]). **Design decision:** we do *not* redistribute the raw archaic
VCFs. We fetch them at panel-build time and redistribute only our **derived** artifacts (computed
marker sites + archaic allele calls at those sites), which are our own small computed asset — cleaner
footing and a few MB instead of tens of GB.

**Build alignment is a solved problem here.** The archaic VCFs are GRCh37; the panel pipeline already
lifts GRCh37→CHM13 with GATK allele-aware liftover offline (`scripts/ancestry-panel/02_liftover_panel_sites.sh`),
and subject genotyping is build-agnostic via `resolve_chip` / `resolve_alignment` (the IBD panel
pre-computes GRCh37/38/CHM13 coordinates; no runtime liftover). Archaic sites lift to CHM13 at build
time exactly like every other panel locus.

---

## 3. Method choice (verified)

**Tier A** = the 23andMe marker-count model: a fixed list of archaic-derived SNPs, 0/1/2 archaic
copies each, summed. Simple, exact, chip-compatible.

**Tier B** = the **Skov et al. 2018 HMM (hmmix)** — "Detecting archaic introgression using an
unadmixed outgroup." Chosen because it is the best fit for a solo-dev per-individual tool:

- **No archaic reference genome required** to *find* segments — it removes variants shared with an
  African outgroup and flags regions dense in the remaining "private" derived variants.
- **Works on unphased diploid genomes** ("heterozygous archaic segments still stand out"). No phasing
  step needed — critical, since we don't have a phasing pipeline.
- **Per-individual** — no population-level dataset needed at runtime.
- It is exactly what **deCODE/Skov ran on 27,566 Icelandic genomes** (Skov 2020, *Nature*), assigning
  **84.5% Neanderthal / 3.3% Denisovan / 12.2% unknown**.

Crucial nuance: the HMM **cannot itself distinguish Neanderthal from Denisovan** (they share a common
ancestor before either coalesces with humans). Archaic reference genomes are used **downstream only**,
to classify each called segment by matching its derived alleles against Altai/Vindija/Chagyrskaya
(Neanderthal) vs Denisova. So Tier B is: *HMM finds segments → archaic VCFs label them.*

Rejected alternatives: **S\*** and **Sprime** need an introgression-free modern reference panel and
are population-scale; **IBDmix** is reference-based and heavier; **DAIseg** (jointly separates
Nea/Den, no post-processing) is attractive but an unreviewed 2025 preprint with self-benchmarks — a
possible future upgrade, not the v1 foundation.

---

## 3a. Marker-list sourcing — RESOLVED: compute our own

The question was whether to ingest an existing public archaic-informative SNP list or compute one.
Verified answer:

- **The Sankararaman 2014 list 23andMe used is NOT publicly available** — not on the Reich Lab
  datasets page, no mirror surfaced; obtainable only by request to the authors, as 23andMe did. Dead
  end for an open-source project.
- **Three openly-licensed alternatives exist**, but each has a catch:
  - **hmmix/Skov callsets** (Zenodo DOI 10.5281/zenodo.14136628, **CC BY 4.0**) — per-SNP files with
    ancestral/derived bases and Neanderthal/Denisovan sharing flags. The best *ready-made* option, but
    **hg38** and *per-individual probabilistic* calls (would need aggregating across individuals into a
    panel).
  - **Sprime/Browning 2018** (Mendeley DOI 10.17632/y7hyt83vxr.1, **CC BY 4.0**) — archaic-match SNPs
    over 1kGP non-Africans + SGDP Papuans; schema only in the in-archive README.
  - **Vernot & Akey** (Akey Lab Google Drive) — segment/BED level, **no stated license**. Weakest.
- **Compute-our-own is fully supported and scientifically standard**, and the build math is ideal:
  the EVA archaic VCFs and the polarity resource are **both GRCh37**, so ancestral/derived assignment
  needs *no* liftover — only the final GRCh37→CHM13 lift the pipeline already runs.

**Decision: compute our own panel** (clean provenance, full control, no private-list or
unclear-license dependency), and use the **hmmix Zenodo set as a citable cross-validation reference**
(CC BY 4.0). Concrete polarity resource: **Ensembl release-75 human ancestral sequence**
(`ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2`,
766 MB, 6-primate EPO / Ortheus, GRCh37 — matches the EVA VCF build exactly). This directly answers
the previous §9 open question.

---

## 4. Offline assets (built by `navigator-panelbuild`)

Following the established asset pattern (bincode `.bin`, `(contig, pos)` CHM13 keys, SHA-256 manifest,
download-on-first-use via `ensure_ancestry_asset`). New `Cmd::ArchaicPanel` + builder module, mirroring
`build_ancient_panel`. New offline stage `scripts/ancestry-panel/08_build_archaic.sh`.

**Asset 1 — `archaic_markers_<build>.bin` (Tier A).** The archaic-informative marker panel, computed
by us (§3a). All inputs GRCh37; polarity assigned pre-lift. Construction (offline):
1. Intersect the four archaic VCFs at biallelic SNP sites passing their BED quality masks (`FilterBed/`).
2. Assign ancestral/derived from the **Ensembl release-75 EPO ancestral sequence** (GRCh37); keep sites
   where an archaic genome is **homozygous-derived** (the introgression donor state).
3. Require the derived allele **rare/absent in a Sub-Saharan African outgroup** (1kGP AFR: YRI/LWK/GWD/MSL/ESN,
   e.g. AFR freq < 1%) and **present in non-Africans** — the introgression signature (Sankararaman/Vernot logic).
   Exact thresholds (AFR cutoff, min non-African freq) to be tuned against the hmmix cross-validation set.
4. Lift GRCh37→CHM13; drop palindromic A/T, C/G sites (already excluded by `is_palindromic`).
5. Store per site: `contig, pos, ref, alt, archaic_derived_allele`, per-archaic genotype
   (Altai/Vindija/Chagyrskaya/Denisova), and a `diagnostic_class` (Neanderthal-diagnostic /
   Denisovan-diagnostic / shared-archaic).

**Asset 2 — `archaic_outgroup_af_<build>.bin` (Tier B baseline).** Per-site African-outgroup allele
frequencies (or a "variable-in-Africans" site set) used to strip shared variants from the test genome
before the HMM. Derived from the 1kGP-on-CHM13 VCFs the panel pipeline already fetches.

**Asset 3 — `archaic_classify_<build>.bin` (Tier B downstream).** Genome-wide archaic diagnostic
alleles (Neanderthal set vs Denisovan set) for labeling called segments. Can be the genome-wide
superset of Asset 1's per-archaic calls.

**Asset 4 — percentile reference (`archaic_marker_dist_<build>.bin`).** Distribution of Tier-A counts
across 1kGP super-populations, so we can render "more than X% of \<EUR\> samples" honestly against a
real cohort instead of an unfalsifiable "other customers" pool.

All four added to the JSON manifest (`AssetManifest` / `.verify()`), path helpers in
`navigator-app/src/lib.rs` (following `ancestry_qpadm_path`), and `ensure_ancestry_asset` /
`ancestry_asset_status`.

---

## 5. Runtime analysis (`navigator-analysis`)

New module `crates/navigator-analysis/src/archaic.rs`.

**Tier A — `count_archaic_markers(&[SiteGenotype], &ArchaicMarkerPanel) -> ArchaicMarkerResult`.**
Intersect the subject's consensus dosages (`consensus_genotypes`, already build-normalized) with the
marker sites; for each covered site add archaic-derived-allele copies (0/1/2). Emit total count,
call rate, count split by `diagnostic_class`, and a percentile from Asset 4. Pure dosage arithmetic —
reuses the exact machinery `estimate_by_allele_frequency` already uses. Works for chip *and* WGS.

**Tier B — `call_archaic_segments(variant_calls, &OutgroupAf, &ArchaicClassify, &GeneticMap, &ArchaicConfig)
-> ArchaicSegmentResult`.** Gated to WGS/VCF (needs genome-wide variants from the existing diploid
SNV caller, [[diploid-snv-caller]]). Pipeline:
1. **Strip shared variants:** drop subject variants present in the African outgroup (Asset 2) and
   sites fixed in humans → keep "private" derived variants.
2. **Windowed Poisson HMM** (e.g. 1 kb windows) over private-variant density. Two states
   (non-archaic / archaic) with different Poisson rates; distance-scaled transitions via `GeneticMap`
   cM. **Viterbi** MAP path + **forward/backward** posteriors — the exact log-space idiom of `roh.rs`
   and `paint_local_ancestry`. Emission differs (Poisson point-process over private SNVs, not het/hom),
   everything else is the established pattern.
3. **Classify each segment:** count derived-allele matches to the Neanderthal set vs Denisovan set
   (Asset 3) → assign `Neanderthal` / `Denisovan` / `Unknown` with a confidence, per Skov 2020.
4. **Aggregate:** total archaic Mb, % of callable genome, Nea/Den/unknown split, per-segment records.

Config (`ArchaicConfig`) carries window size, Poisson rates, min-segment length, classification
thresholds — same shape as `RohConfig` / `PaintParams`.

---

## 6. Domain, persistence, orchestration, UI

- **Domain** (`navigator-domain/src/ancestry.rs`): `ArchaicMarkerResult { total_copies, possible_copies,
  call_rate, neanderthal_copies, denisovan_copies, percentile, cohort }` and `ArchaicSegmentResult {
  segments: Vec<ArchaicSegment>, total_mb, pct_genome, neanderthal_pct, denisovan_pct, unknown_pct }`
  with `ArchaicSegment { contig, start, end, source: ArchaicSource, posterior, n_private_snps }`.
- **Store**: migration `00XX_consensus_archaic` — a `consensus_archaic` table keyed by
  `biosample_guid` + `consensus_sig` (staleness = consensus `last_reconciled_at`), results as JSON,
  exactly like `consensus_roh` / `consensus_painting`.
- **App** (`navigator-app/src/haplogroup.rs`): `estimate_archaic_from_consensus(guid)` runs Tier A
  always; runs Tier B when the subject has a WGS/VCF source with genome-wide calls; persists under
  `CONSENSUS_SOURCE_ID`. Command `EstimateArchaicFromConsensus` wired in `worker.rs`.
- **UI** (`navigator-ui/src/ui/central.rs`): a new **"Archaic Ancestry"** card in
  `DetailTab::Ancestry`, following the fine/qpAdm/painting/ROH cards. Presentation, matching the
  23andMe reference but honest:
  - **Headline:** Neanderthal variant **count** + **percentile** ("more than X% of European samples")
    — a count, never a fake "% Neanderthal", exactly per 23andMe's own guidance.
  - **WGS bonus:** when Tier B ran, a second line with **% of genome** + a **chromosome painting** of
    archaic segments (reuse `draw_chromosome_painting` / `draw_roh` in `charts.rs`), colored by
    Nea/Den/unknown.
  - **Denisovan:** shown only when it clears a confidence floor (§7); otherwise an explicit
    "No reliable Denisovan signal — expected for European ancestry" note.
  - An educational "decoded" blurb + optional trait-association list (deferred, §8).

---

## 7. Setting honest expectations (the Denisovan trap)

Verified population baselines: a typical **European carries ~1.5–2% Neanderthal** (single pulse
~45–60 kya; Europeans ~1.8%, East Asians ~2.3–2.6%). **Denisovan is effectively zero in Europeans** —
in the Icelandic data Denisovan was 3.3% *of* the ~2% archaic total, i.e. a European's Denisovan
signal sits **at the noise floor**. Denisovan ancestry concentrates in Oceanians / East & South Asians.

Design consequence: **never emit a fabricated small Denisovan number for a European.** This is the
same failure mode that got ancient ancestry disabled ([[ancient-ancestry-broken]] — fabricated numbers
from centroids that sat on top of each other). Follow the `estimate_qpadm_ancestry` precedent: report
Denisovan only when Tier B's classified Denisovan fraction clears a confidence threshold; otherwise say
"none reliably detected" and explain why. For the ground-truth European sample (GFX0457637 = James),
the expectation is ~1.5–2% Neanderthal and no Denisovan — a natural validation target, cross-checkable
against his real 23andMe Neanderthal count if available.

---

## 8. Phasing

- **Phase 1 (MVP):** Asset 1 + Asset 4, Tier A `count_archaic_markers`, domain/store/UI card. Delivers
  the full 23andMe-equivalent for **both chip and WGS**, reusing the ancestry-panel machinery almost
  verbatim. Low risk. Validate the marker count against James's actual 23andMe report.
- **Phase 2:** Assets 2 + 3, Tier B segment HMM for WGS, % genome + chromosome browser + Nea/Den
  attribution. Higher effort (Poisson-HMM emission, outgroup strip, classification) but reuses the
  `roh.rs` / painting HMM idiom and the existing diploid caller.
- **Phase 3:** Trait associations (23andMe-style introgressed-variant → phenotype list; each needs a
  curated GWAS-backed table — treat as optional/curatorial), export/PDS records
  (`export.rs`/`publish.rs`), and possibly a DAIseg-style joint Nea/Den upgrade if it matures.

## 9. Open questions

### Resolved

**Marker-list sourcing — VERDICT: compute our own (see §3a).** The exact Sankararaman 2014 list
23andMe used is *not* publicly downloadable (not on the Reich Lab page, no mirror; author-request
only), so it is unusable for an open-source project. Three openly-licensed alternatives exist but each
is disqualified for direct use: the **hmmix/Skov** callsets (Zenodo, CC BY 4.0) carry clean per-SNP
polarity + Nea/Den sharing but are hg38 and per-individual probabilistic; **Sprime/Browning 2018**
(Mendeley, CC BY 4.0) is usable but schema-opaque; **Vernot & Akey** (Google Drive) is segment-level
with no stated license. We therefore **compute our own panel** from the EVA archaic VCFs + **Ensembl
release-75 EPO ancestral alleles** (both GRCh37 — no liftover for polarity) + a **1kGP AFR outgroup**,
and redistribute only our derived sites (avoids every private-list and unclear-license dependency).
The **hmmix Zenodo set (CC BY 4.0)** is retained as a citable cross-validation reference. One sub-task
carries into Phase 1: finalize the exact site-selection thresholds (AFR-frequency cutoff, min
non-African frequency) by calibrating against that cross-validation set.

### Still open before Phase 1

*(Reviewed 2026-07-26 when the implementation plan below was written; two of the three are now settled.)*

1. **Chip overlap — MEASUREMENT SCHEDULED (M1 checkpoint B).** How many informative sites actually
   overlap 23andMe/AncestryDNA v5 chip content after lift to CHM13? Unanswerable until candidate
   sites exist, so it becomes a gating checkpoint *inside* M1 rather than a blocker before it. It
   sets the realistic chip call rate and the honest ceiling on the count.
2. **Tier B inputs — RESOLVED.** The AFR outgroup needs no new data source: stage 1 of the panel
   pipeline (`01_fetch.sh`) already fetches the 1000G-on-CHM13 VCFs carrying **per-super-pop INFO
   AC/AN** (this is what `03_select_panel.sh` Fst-ranks AIMs from), which is exactly what Asset 1
   step 3 (AFR freq < 1%) and Asset 2 require. The callability mask reuses the existing
   1kGP-on-CHM13 assets, and the private-variant source is the diploid caller's genome-wide VCF
   ([[diploid-snv-caller]]).
3. **Percentile cohort — DECIDED for v1: 1kGP super-population.** Not the user's own inferred
   fine-ancestry group. Rationale: keying the percentile to the inferred ancestry couples this report
   to the ancestry estimate, so an ancestry error would silently shift the archaic percentile — and
   an unexplained shift in a headline number is the failure mode §7 exists to prevent. A fixed
   super-pop cohort is falsifiable and independently checkable. A fine-pop percentile is a
   worthwhile **Phase 3** refinement once Tier A is trusted; Asset 4 should therefore store the
   distribution per population, not a pre-reduced summary, so the cohort can be re-keyed later
   without an asset rebuild.

---

## 10. Implementation plan

Written 2026-07-26 on branch `feat/archaic-ancestry`. Milestones map onto the §8 phases: **M1–M2 =
Phase 1 (MVP)**, **M3 = Phase 2**, **M4 = Phase 3**. Each milestone's touchpoints were verified
against the tree before this plan was written.

### M0 — Decisions (no code) — DONE

The §9 review above: Q2 resolved, Q3 decided, Q1 converted into an M1 checkpoint.

### M1 — Offline assets (the bulk of the work)

New offline stage `scripts/ancestry-panel/08_build_archaic.sh` plus a `navigator-panelbuild` module
and `Cmd::ArchaicPanel`, mirroring `pca::build_ancient_panel` (`main.rs:130`).

**Asset 1 — `archaic_markers_<build>.bin`**, per §4:
1. Intersect the four EVA archaic VCFs at biallelic SNP sites passing their `FilterBed/` masks.
2. Assign ancestral/derived from the Ensembl release-75 EPO ancestral sequence — both inputs are
   GRCh37, so **polarity is assigned pre-lift with no liftover**. Keep sites where an archaic genome
   is homozygous-derived.
3. Require the derived allele rare in the AFR outgroup (freq < 1 %) and present in non-Africans,
   using the per-super-pop INFO AC/AN established in §9 Q2.
4. Lift GRCh37→CHM13 through the existing `02_liftover_panel_sites.sh`; drop palindromic sites
   (`is_palindromic`).
5. Emit per site: `contig, pos, ref, alt, archaic_derived_allele`, the four per-archaic genotypes,
   and a `diagnostic_class` (Neanderthal-diagnostic / Denisovan-diagnostic / shared-archaic).

**Asset 4 — `archaic_marker_dist_<build>.bin`**: Tier-A counts across 1kGP samples, stored **per
population** (see §9 Q3) and reduced to a super-pop percentile at read time.

Both assets follow the established pattern: bincode `.bin`, `(contig, pos)` CHM13 keys, an entry in
the SHA-256 manifest (`AssetManifest::verify`), a path helper beside `ancestry_qpadm_path`
(`navigator-app/src/lib.rs:1098`), and download-on-first-use via `ensure_ancestry_asset`
(`import_unified.rs:916`). Asset 1 is small enough to bundle; check its size against the
`ON_DEMAND_PREFIXES` policy in `packaging/stage-assets.sh` before deciding.

**Two gating checkpoints inside M1:**

- **Checkpoint A — threshold calibration.** Fix the AFR cutoff and the min-non-African frequency by
  cross-validating against the hmmix Zenodo set (CC BY 4.0, §3a). This is the scientific core of the
  milestone: the site list *is* the product, and every downstream number inherits its errors. Do not
  treat the §4 example thresholds as final.
- **Checkpoint B — chip overlap (§9 Q1).** Intersect the final site list with a real 23andMe v5 raw
  file to measure the actual chip call rate. This number goes into the UI copy as the honest ceiling.

**Licensing guard rail (§2):** the raw EVA archaic VCFs and the 766 MB Ensembl ancestral tarball are
**fetch-at-build-time only, never redistributed**. Only our derived sites ship. `08_build_archaic.sh`
must fetch to the pipeline's raw/ area, like the AADR handling in [[germanic-panel-sources]].

### M2 — Tier A runtime + UI (the shippable MVP)

- **Analysis:** new `crates/navigator-analysis/src/archaic.rs` with
  `count_archaic_markers(&[SiteGenotype], &ArchaicMarkerPanel) -> ArchaicMarkerResult` — pure dosage
  arithmetic over `consensus_genotypes` (`lib.rs:2415`), so chip and WGS both work with no decode.
- **Domain:** `ArchaicMarkerResult` in `navigator-domain/src/ancestry.rs` per §6.
- **Store:** migration `0039_consensus_archaic` (0038 is the current head), keyed by
  `biosample_guid` + `consensus_sig`, modelled on `consensus_roh` / `consensus_painting`.
- **App/worker:** `estimate_archaic_from_consensus(guid)` persisted under `CONSENSUS_SOURCE_ID`, with
  an `EstimateArchaicFromConsensus` command.
- **UI:** an "Archaic Ancestry" card in `DetailTab::Ancestry`. Headline is a **count + percentile**,
  never a "% Neanderthal" (§1, §7).

**Validation gate:** the marker count for GFX0457637 must be checkable against that person's real
23andMe Neanderthal count, and the chip-derived and WGS-derived counts must agree within the call-rate
limit measured at checkpoint B. *This requires the actual 23andMe report figure as an external
oracle — without it M2 verifies internal consistency only, not correctness.*

### M3 — Tier B segment HMM (WGS/VCF)

Assets 2 + 3, then `call_archaic_segments` per §5: strip AFR-shared variants → windowed Poisson HMM
over private-variant density (two states, cM-scaled transitions via `GeneticMap`) → Viterbi MAP path
+ forward/backward posteriors in log space, following the `roh.rs` / `paint_local_ancestry` idiom →
classify each segment against the Nea vs Den diagnostic sets → aggregate to Mb, % of callable genome,
and a Nea/Den/unknown split. UI adds the % line and a chromosome painting via `draw_roh` /
`draw_chromosome_painting` (`navigator-ui/src/charts.rs:117`, `:227`).

**Ship Tier B behind a feature gate until validated on real data.** The precedent is
[[ancient-ancestry-broken]]: a plausible-looking but fabricated breakdown shipped and had to be
disabled. Tier B has more moving parts than Tier A and the same blast radius, so it gets the same
discipline — a constant gate that also covers the read *and* publish paths, flipped only once the
validation targets below pass.

**Validation targets:** ~1.5–2 % Neanderthal and **no** Denisovan for the European ground-truth
sample (§7); the overall Nea/Den/unknown split should be in the neighbourhood of Skov 2020's
84.5 / 3.3 / 12.2 on Icelanders. The Denisovan confidence floor must be exercised by a test that
asserts "none reliably detected" rather than a small number.

### M4 — Phase 3 (optional)

Trait associations (each needs a curated GWAS-backed table — curatorial work, not engineering),
export + PDS records, fine-pop percentile re-keying, and a possible DAIseg-style joint Nea/Den
upgrade if that method matures past preprint.

### Sequencing note

M1 is the long pole and everything else depends on Asset 1, but M2 is where the user-visible feature
lands. M3 is separable — if effort runs short, **M1 + M2 is a complete, honest, shippable report**
(it is exactly what 23andMe ships), and Tier B can follow later without rework.
