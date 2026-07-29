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

**Validation gate — the external oracle (obtained 2026-07-26).** For GFX0457637 (= James), the real
23andMe v5 report reads **191 Neanderthal variants out of 7,462**.

Two things follow directly:

- **The denominator is *copies*, not sites.** 7,462 = 2 × 3,731, and 3,731 is exactly the v5 assayed
  site count in §1 — so 23andMe reports archaic-allele *copies* carried out of copies assayed. This
  confirms the §6 domain shape: `total_copies` = 191, `possible_copies` = 7,462. The observed rate is
  **2.56 % of assayed copies** (0.051 archaic copies per assayed site).
- **Do NOT treat 191 as a number our panel must reproduce.** We compute our own marker panel (§3a)
  from different inputs, so it will have a different size and different membership; the raw count is
  panel-relative and a direct equality check is meaningless. Tuning M1's thresholds until our count
  hits 191 would be fitting the panel to one sample — precisely the kind of circular validation that
  §3.2 of the ancient-ancestry investigation already burned us on (the A′ "pass" that turned out to
  be a circular ∩chip comparison).

**The comparison that is actually valid** is the per-site archaic rate **on the intersection of our
panel with the v5 chip content** — i.e. restrict both to shared sites and compare copies-per-site
(expect ≈ 0.051 on the 23andMe side). Checkpoint B already computes that intersection, so this costs
nothing extra. Secondary, weaker checks: our whole-panel rate should land in the same neighbourhood,
and Tier B's independent % -of-genome estimate should be consistent with a ~1.5–2 % European
Neanderthal fraction (§7) — two different methods agreeing is worth more than either matching a
vendor's count.

Also required at M2: the chip-derived and WGS-derived counts for the same person must agree within
the call-rate limit measured at checkpoint B. That is the cross-source stability gate that ancient
ancestry failed for a long time ([[ancient-ancestry-broken]]), and it is a genuinely independent
check because it needs no vendor number at all.

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

### M1 — as built (2026-07-28)

The pipeline ran end to end. Inputs: 202.8 GB of EVA all-sites VCFs (four genomes × 22 autosomes)
plus the EPO ancestral sequence. Funnel at the shipped thresholds:

| stage | count |
|---|---|
| candidates (GRCh37, polarized, archaic hom-derived) | 2,032,698 |
| lifted to CHM13 | 2,031,406 (99.94 %) |
| dropped — no outgroup AF | 627,063 |
| dropped — too common in AFR (>1 %) | 968,493 |
| dropped — too rare outside AFR (<1 %) | 247,165 |
| dropped — palindromic | 29,295 |
| dropped — failed CHM13 orientation | 6 |
| **kept** | **159,384** (118,522 Nea / 6,406 Den / 34,456 shared) |

**The orientation step earns its keep:** 8,844 sites were ref/alt-swapped relative to CHM13. Without
it they would have been silently inverted — the §7.16 defect from the ancient-ancestry work. Only 6
sites failed orientation outright, so the lift itself is clean.

The 627 k "no outgroup AF" losses are benign: those sites are absent from the 1000G AF VCFs because
1000G is monomorphic there, so their derived allele is at 0 % outside Africa and fails the floor
regardless.

### Checkpoint B — chip overlap, and what it revealed

Measured against a real 23andMe v5 raw file (1,407,553 markers). **Chip overlap is not the
constraint**: ~10,000 panel sites are assayed, comfortably more than the 3,731 Sankararaman sites
23andMe uses, so Tier A is viable on chip data.

The *rate* is the interesting part. Sweeping the non-African floor, measuring archaic copies per
called site on the chip intersection:

| `min_non_afr_freq` | panel | on chip | copies/site | vs 23andMe (0.0512) |
|---|---|---|---|---|
| 0.05 (initial) | 48,687 | 4,807 | 0.1545 | 3.02× |
| 0.02 | 113,194 | 8,236 | 0.1140 | 2.23× |
| **0.01 (adopted)** | **159,384** | **10,006** | **0.0991** | **1.94×** |
| 0.005 | 199,072 | 10,721 | 0.0950 | 1.86× |
| 0.001 | 267,852 | 10,907 | 0.0936 | 1.83× |

The floor was set far too high. At 0.05 the surviving sites piled up against it (p10 = 0.054,
mean 0.089, median 0.077), keeping the common tail and discarding the rare variants that make up
most of a real introgression panel. A *ceiling* cannot help: while the floor sets the distribution
it cannot push the rate below ~0.139.

Two inferences drawn here were **wrong**, and checkpoint A overturned both — recorded because the
reasoning is a trap worth not repeating:

- *"Below ~0.005 the Denisovan-diagnostic count inflates implausibly, so that is noise."* It is not.
  hmmix independently contains 40,408 Denisova-only diagnostic SNPs, so a large Denisovan count is
  signal. That was proxy reasoning with no oracle behind it.
- *"The residual ~1.83× rate gap is structural — a frequency filter cannot separate introgressed
  from out-of-Africa-specific variants."* Too pessimistic: measured precision against hmmix is
  ~78 %, so the panel is mostly genuinely introgressed sites. The rate gap reflects *which*
  introgressed sites each panel selects, not junk in ours.

### Checkpoint A — hmmix calibration (2026-07-28)

Source: the hmmix Zenodo callset (DOI 10.5281/zenodo.14136628, **CC BY 4.0**),
`hg38_1000g_SNPS.txt`. Its 370,960 `DAV` (directly diagnostic) SNPs were lifted hg38 → CHM13 with
the existing chain (553,597 of 553,763 lifted) and used as the positive set. `linkedDAV` rows are
excluded — they are LD-linked, not independent evidence.

**Three independent validations, all passing:**

| check | result |
|---|---|
| polarity (derived base) agreement on the overlap | **99.99 %** (125,304 vs 12) |
| diagnostic class on the diagonal | **98.5 %** |
| Neanderthal ↔ Denisovan confusions | **zero** |

Every off-diagonal classification is "we said lineage-specific, hmmix said shared" — the
conservative direction `classify_diagnostic` predicts, since a masked-out Denisova reads as
no-evidence rather than as absence. The polarity result is the important one: it independently
confirms the EPO-based ancestral/derived assignment, which is the one thing this panel cannot
afford to get wrong.

**Threshold sweep, scored by F1 against the hmmix positives:**

| `min_non_afr_freq` : `max_afr_freq` | panel | precision | recall | F1 |
|---|---|---|---|---|
| 0.05 : 0.01 (initial) | 48,687 | 74.7 % | 9.8 % | 0.173 |
| 0.01 : 0.01 | 159,384 | 78.6 % | 33.8 % | 0.473 |
| 0.001 : 0.01 | 267,852 | 79.8 % | 57.6 % | 0.669 |
| **0.0005 : 0.01 (ADOPTED)** | **299,958** | **78.4 %** | **63.4 %** | **0.701** |
| 0.0 : 0.01 | 371,839 | 64.7 % | 64.9 % | 0.648 |
| 0.001 : 0.02 | 291,283 | 73.8 % | 57.9 % | 0.649 |
| 0.001 : 0.05 | 336,338 | 63.9 % | 58.0 % | 0.608 |

Both bounds are load-bearing. Precision is flat (~75–80 %) while a floor exists at all, so lowering
it is nearly free recall — but removing it entirely (0.0) collapses precision to 64.7 % for almost
no recall gain. Relaxing the **AFR ceiling** costs precision without buying recall, which identifies
it as the criterion actually doing the archaic-specificity work; the non-African floor is only
suppressing the very rarest noise.

**Shipped panel: 299,958 sites** (202,097 Neanderthal / 36,268 Denisovan / 61,593 shared).

### Asset 4 — percentile reference (built 2026-07-28)

`archaic_marker_dist_chm13v2.0.bin`, built by `panelbuild archaic-dist`: every 1kGP sample scored
through the same marker arithmetic the app runs, from the 3202-sample CHM13 BCF (native build, no
liftover), grouped per fine population. 2,504 labelled samples over 26 populations; the 698
unlabelled (related) samples are dropped.

**This is an independent biological validation of the whole panel**, on data nothing was tuned to:

| super-pop | n | mean archaic copies |
|---|---|---|
| AFR | 661 | 1,223 |
| EUR | 503 | 13,234 |
| AMR | 347 | 13,748 |
| SAS | 489 | 15,639 |
| EAS | 504 | 18,201 |

Africans carry ~10× fewer archaic alleles, and East Asians exceed Europeans by ~38 % — both are the
textbook result (EAS ~2.3–2.6 % vs EUR ~1.8–2.0 % Neanderthal). The non-zero AFR floor is expected:
the AFR ceiling is 1 %, not 0.

### The percentile — but NOT the count — is incomparable across data types

**Scope first, because this is easy to over-read.** The headline is a *count over what was actually
assayed* — copies carried of copies possible, exactly the shape of 23andMe's "191 of 7,462"
(= 2 × 3,731 sites), and exactly what `ArchaicMarkerResult` already models via `total_copies` /
`possible_copies` / `call_rate`. That form is self-normalising: whatever subset of the panel a test
covers becomes the denominator, so **chip and WGS each get an honest headline with no cross-data-type
comparison and nothing below is a blocker for it.** For the ground-truth sample the chip figure is
1,005 of 21,478 (4.68 %), against 23andMe's 191 of 7,462 (2.56 %) — different because the panels
ascertain different sites, which §1 and §10 already forbid presenting as equivalent.

What follows applies **only to the percentile**, which ranks that number against a cohort and so
does require the cohort to have been scored on comparable sites.

Measured, not anticipated. A 23andMe v5 chip covers **3.6 %** of the calibrated panel (10,739 of
299,958 sites), so the ground-truth sample's raw chip count is 1,005 against a EUR WGS mean of
13,234. Rendering that against the cohort distribution puts **every chip user at the 0th
percentile**, purely as a call-rate artifact.

Naive rate-scaling does not fix it either, and this is the part that would be easy to miss: the
subject's chip *rate* is 0.0936 copies/site against a EUR WGS mean rate of 0.0441 — **2.1× too
high**, because array content is deliberately biased toward common variants while the calibrated
panel is mostly rare ones. The chip-overlapping sites are the panel's common tail, so scaling a rate
measured there across the whole panel over-estimates by roughly 2×.

So neither the raw count nor a scaled rate can be compared to this cohort. Options for M2, in
rough order of preference:

1. **Per-site, per-population derived-allele frequencies** in Asset 4 instead of (or alongside)
   per-sample totals. Compact, and it lets the expected count and variance be computed analytically
   for *any* subset of called sites — so chip, WGS and partial data all get an honest percentile.
2. **Two cohorts** — one scored over all panel sites (WGS) and one over the chip-overlap subset,
   picking whichever matches the subject's data. Pragmatic, mirrors how a vendor's chip-vs-chip
   percentile works, but needs a cohort per chip build.
3. Report the percentile **only** for WGS and suppress it for chip data, showing the bare count.
   Honest but weak, given chip is the common case.

Until one is implemented, the Tier A card renders the count and rate as normal and simply omits the
percentile for chip input — the count does not depend on any of this.

### Manifest (2026-07-28)

Both assets are now in `ancestry_manifest_chm13v2.0.json` (12 → 14 entries), so
`read_verified_asset` checks them instead of loading through the unverified passthrough:

| asset | bytes |
|---|---|
| `archaic_markers_chm13v2.0.bin` | 14,220,626 |
| `archaic_marker_dist_chm13v2.0.bin` | 10,829 |

`panelbuild manifest` globs `*_<build>.bin`, so no builder change was needed. Regenerating rewrites
every entry, so the pre-existing 12 were diffed against a backup
(`.bak-prearchaic`): **no hash changed and nothing was removed**, which matters because a changed
hash would invalidate an asset already published to the release. Both new hashes were then verified
independently of the builder.

**M1 is complete.** What remains before the feature is usable is M2, and two items belong with it
rather than here, because both would be dead code without a consumer:

- **App-side plumbing** (§4): `archaic_markers_path` / `archaic_marker_dist_path` beside
  `ancestry_qpadm_path`, plus `ensure_ancestry_asset` and `ancestry_asset_status` entries so the
  assets download on first use and appear in the Settings asset list.
- **The percentile comparability fix** — see the section above. It gates the *percentile* only; the
  count-and-rate headline can ship without it.

**Publishing is deliberately not done.** The raw archaic inputs must never be published (§2), and
the derived assets go to the GitHub asset release via `packaging/publish-assets.sh` as a separate,
explicit step. Note `archaic_markers_*.bin` is 14 MB — small enough to bundle, unlike the 133 MB
haplotype panel, so it does not need an `ON_DEMAND_PREFIXES` entry in `packaging/stage-assets.sh`.

### M2 — as built, and the coverage gap the real-data run exposed (2026-07-28)

Built end to end: `count_archaic_markers` → `consensus_archaic` (migration 0039) →
`estimate_archaic_from_consensus` → worker command/event → an Ancestry-tab card, plus a headless
`navigator archaic --subject <id>` mirroring `deep-ancestry`.

Running it on the ground-truth subject produced a working result and a **design gap that unit tests
could not have caught**:

```
Archaic markers (Tier A): 914 of 15354 copies
  7677 of 299958 panel sites called (2.6%)
```

**The autosomal consensus is 1240k-restricted.** §5 assumed the consensus dosages span the genome;
they do not — `DiploidProfile` carries genotypes only at the ~1.23 M 1240k/IBD panel loci. Our
archaic panel intersects that set at **7,679 sites**, and the app called 7,677 of them. So Tier A
sees 2.6 % of the panel *regardless of whether the subject has WGS*, and 97 % of the calibrated
panel is unreachable through this path.

**Resolved — option 1 taken.** `genotype_archaic_for_alignment` genotypes the panel directly from
the subject's best-callable alignment via `caller::genotype_sites_all_contigs`, cached per alignment
under a kind salted with the panel's manifest hash (recalibrating the panel invalidates stale
genotypes rather than mixing site sets). The consensus is still read first, so chips and
non-alignment sources keep contributing where the alignment has no depth; the direct call wins at any
site it covers. Measured on the ground-truth subject: **2.6 % → 100.0 % coverage** (299,932 of
299,958 sites) in a single 2 m 40 s decode pass.

The rejected alternatives, for the record: folding the archaic loci into the consensus target set
would be more efficient but changes the consensus contract and needs everyone re-genotyped;
restricting the panel to the 1240k intersection would discard the calibration and cap WGS at chip
resolution.

### Per-build loci — all three builds genotype natively (2026-07-28)

`ArchaicSite` now carries `grch37` / `grch38` loci beside its canonical CHM13 one, mirroring
`IbdPanelSite`, so the CHM13-only restriction is gone.

- **GRCh37 is exact, not lifted.** Those are the archaic VCFs' own hg19 coordinates and alleles
  carried straight through, so that build has no liftover and no strand risk whatsoever.
- **GRCh38 is lifted and oriented** against an hg38 reference exactly as the CHM13 pass is — the
  hg19→hg38 lift is not allele-aware either. 2,031,439 of 2,032,189 candidates carry an hg38 locus;
  750 were dropped where the reference base matched neither allele.
- A dosage measured on a non-CHM13 build is **re-keyed** to the CHM13 alleles
  (`ArchaicSite::rekey_dosage`) before counting, because that build's ref/alt may be swapped *or*
  strand-flipped. Feeding the raw dosage through would invert exactly those sites, silently — the
  same failure class the CHM13 orientation pass exists to prevent.

Site selection is unchanged (299,958); the asset grows 13.6 → 26.1 MB, still small enough to bundle.

**Cross-build validation** — the same person's GRCh38 alignment (`WGS229.b38.bam`) against the
CHM13 result, both at ~100 % coverage:

| | CHM13 | GRCh38 | Δ |
|---|---|---|---|
| copies | 12,126 | 12,112 | **0.12 %** |
| sites called | 299,932 | 299,832 | |
| Neanderthal | 9,460 | 9,449 | 0.12 % |
| Denisovan | 254 | 254 | 0 |

That agreement is the real test of the re-keying: had it been wrong, the swapped and strand-flipped
sites would have inverted and the two builds would have diverged grossly rather than by 0.12 %. It is
also the cross-source stability gate from the M2 plan, met on the cross-build axis.
`navigator archaic --alignment <id>` exists to run this check.

### Validation — the percentile problem is solvable, and the numbers are sane

Scoring the 1kGP cohort on **exactly the 7,679 sites the subject was scored on**:

| super-pop | n | mean copies |
|---|---|---|
| AFR | 661 | 80.2 |
| SAS | 489 | 940.9 |
| EUR | 503 | 1003.4 |
| AMR | 347 | 1036.2 |
| EAS | 504 | 1322.7 |

The population ordering survives on this subset (AFR ~12× lower; EAS above EUR), and the subject
lands at **914 copies — the 10th percentile within EUR** (mean 1003, range 819–1234): low-normal,
well inside the distribution.

That is the **percentile fix demonstrated**: comparing subject and cohort on the *same* called-site
basis is valid.

**Now wired.** With direct genotyping the subject reaches ~100 % of the panel — the same basis the
cohort was scored on — so the percentile is filled when
`call_rate >= ARCHAIC_PERCENTILE_MIN_CALL_RATE` (0.90) and left `None` below it. The cohort is keyed
to the subject's inferred **super**-population, not their fine population (§9 Q3): a wrong percentile
is worse than a coarse one, and no ancestry estimate means no percentile rather than a default.
The ground-truth subject reports **12,126 of 599,864 copies at 100 % coverage, "more than 7 % of
EUR"** — consistent with the 5th–10th percentile the two independent hand calculations gave.

### Chip percentiles — solved (2026-07-28)

Asset 4 now stores **per-site per-population derived frequencies** in panel order, so the cohort's
expected count and variance can be computed over exactly the sites a subject called. The call-rate
gate is gone: `percentile_for_called` is valid at any coverage, and a chip is compared against what
the cohort would have scored on those same ~3 % of sites.

The count is modelled as a sum of per-site binomials under Hardy-Weinberg, normal-approximated, with
the variance scaled by a **measured LD inflation** — archaic alleles travel in linked blocks, so the
independent-sites variance badly understates the spread.

Two things had to be got right, both found by measuring rather than assuming:

1. **The inflation is not a constant.** It is 52.4× on the full panel but 5.3× on a 2.6 % subset of
   the same panel, because a sparse subset samples fewer sites per linked block. Applying one
   full-panel factor to a chip would over-widen the deviation ~3× and squash every percentile toward
   50. No simple block model fits both measurements (solving for a common block size gives a
   negative size), so the inflation is measured at a ladder of densities and log-interpolated.
2. **Random subsets are not chips.** A random 3 % subset gives EUR 2.0× while the *real* 1240k
   intersection at 2.6 % gives 5.3× — array and capture content is spatially clustered and retains
   far more linkage. Mixing random and real rungs also made the ladder non-monotonic, which would
   have *lowered* the inflation as density rose. Real site sets (`--subset-sites`) therefore
   **displace** the synthetic rungs; only the full-panel rung, real by construction, is kept
   alongside them.

Validated against empirical rank on the same subsets:

| input | coverage | model | empirical rank |
|---|---|---|---|
| chip (23andMe v5) | 3.58 % | 7 % | 1 % |
| WGS, direct genotyping | 100 % | 12 % | ~5–7 % |

The model runs **conservative** — it reports closer to the middle than the truth — which is the
right direction to err for a consumer-facing number, but it should not be presented as a precise
rank. Exactness would need the cohort's per-sample genotypes at runtime (~188 MB bit-packed, the
shape `HaplotypeReference` already uses); that is the upgrade path if the approximation ever proves
insufficient.

### M4 — Phase 3 (optional)

Trait associations (each needs a curated GWAS-backed table — curatorial work, not engineering),
export + PDS records, fine-pop percentile re-keying, and a possible DAIseg-style joint Nea/Den
upgrade if that method matures past preprint.

### Sequencing note

M1 is the long pole and everything else depends on Asset 1, but M2 is where the user-visible feature
lands. M3 is separable — if effort runs short, **M1 + M2 is a complete, honest, shippable report**
(it is exactly what 23andMe ships), and Tier B can follow later without rework.
