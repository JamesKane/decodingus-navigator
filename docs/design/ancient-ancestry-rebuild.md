# Ancient ancestry — why it's disabled, and what a correct rebuild takes

**Status:** feature gated off (`navigator_app::ANCIENT_ANCESTRY_ENABLED = false`). Original PCA
implementation disabled 2026-07-13 (§1–§2). Option-A rebuild (§3) built and tested 2026-07;
**it fails the §3.4 stability gate — see §4** — so the feature stays off pending the §4 fix.
**Scope of this doc:** the evidence that the original implementation is unsound (§1–2), the rebuild
design (§3), and what the rebuild attempt found — including a reconsideration of the §3.1 method
choice (§4).

---

## 1. What was shipped, and why it's wrong

The "Ancient ancestry" card reported components like *Western Hunter-Gatherer 72.6%, Steppe
pastoralist 23.2%, Anatolian Farmer 2.5%* for a British subject. Those numbers are not merely
imprecise — they are **fabricated**, in the sense that no property of the input data produced them.
Three compounding defects, in increasing order of severity.

### (a) The reference set is contaminated

`ancestry_pca_ancient_<build>.bin` contains **173 populations: the 168 modern reference populations,
plus exactly 5 ancient ones** — `WHG`, `EHG`, `ANF`, `Iran_N`, `Steppe`. It is the modern panel with
five ancient populations appended, not an ancient reference.

Consequences:
- A **mixture** model trivially reconstructs a Briton as `English 55% + GBR 31%` — "British" is *in
  the reference set*, so the ancient sources are never needed. Correct math, wrong reference.
- Modern populations (`Colombian`, `Puerto Rican`) surface as *"ancient sources"* in the UI, because
  they are simply other entries in the same list.

### (b) The estimator is a classifier, not an admixture model

`ancestry::classify_pca` computes `exp(-0.5 · mahalanobis²)` from the sample to each population
centroid and normalises. That is a **GMM membership posterior** — *"which population **is** this
sample?"*. A modern European is not a **member** of WHG; they are a **mixture** of ancient sources.
Presenting a membership posterior as an admixture percentage is a category error.

It is also numerically brittle: `exp(-d²)` is exponentially sensitive to distance, so a small change
in the projected coordinates flips the answer. Observed on one subject, same person, same method:

| Genotype source | PCA sites used | Result |
|---|---|---|
| Consensus | 9,662 | WHG 72.6%, Steppe 23.2%, ANF 2.5% |
| 28× CHM13 BAM | 13,142 | ANF 49.5%, Steppe 40.8%, WHG 7.4% |

`project_pca` compounds this: it rescales the projection by `total_sites / used_sites` to "un-shrink"
for missing data, which amplifies noise precisely when sites are sparse.

### (c) **Fatal:** the ancient centroids carry no ancient signal

The five ancient centroids sit **on top of the modern European centroids** in PC space:

| Population | PC1 | PC2 |
|---|---|---|
| **WHG** *(ancient)* | +29.90 | −23.91 |
| English *(modern)* | +28.47 | −26.44 |
| GBR *(modern)* | +28.72 | −24.68 |
| **ANF** *(ancient)* | +25.76 | −24.72 |
| Sardinian *(modern)* | +25.45 | −25.56 |

WHG ≈ English. ANF ≈ Sardinian. In any valid population-genetics PCA, WHG is an **extreme outlier**
far from every modern European, and modern Europeans fall *between* the ancient sources — that
betweenness is exactly what makes admixture decomposition possible.

The PCA axes themselves are fine (YRI at PC1 −100, Han at +43/+31 — correctly separated), so this is
not a PCA bug. It is the classic **shrinkage artifact of projecting ancient samples onto a PCA built
from modern variation**: ancient samples are low-coverage, pseudo-haploid and highly divergent, so a
naive projection collapses them toward the centre of the modern cloud.

Because the ancient sources are nearly collinear and nearly co-located, **any** model over them is
ill-conditioned. Verified: filtering to ancient-only sources *and* switching to a proper simplex
mixture still yields `Steppe 80.5%, WHG 16.4%, ANF 3.1%` — Anatolian Farmer at 3% where ~30% is
expected for any European. **No code change rescues this.** The asset is the problem.

### What is *not* broken

- **Input signals are clean.** Chip-vs-WGS dosage concordance at the shared panel sites: 23andMe
  97.7%, AncestryDNA 97.7%, **0.0% orientation flips**. Nothing is strand-crossed.
- **Modern/fine admixture is sound** and stays enabled: `CEU 36.8 / GBR 23.6 / TSI 14.2 / IBS 13.1 /
  FIN 12.0` is a textbook British profile.

---

## 2. What "disabled" means today

`ANCIENT_ANCESTRY_ENABLED = false` gates four sites. The read paths are gated as well as the compute
path, so that **stale rows persisted by an earlier build cannot resurface**:

| Site | Effect |
|---|---|
| `haplogroup::estimate_ancestry_from_consensus` | stops computing/persisting `PCA_PROJECTION_GMM` + `G25_NMONTE` |
| `brief::subject_brief` | ancient components absent → Simple card, DNA-story HTML export, and LLM facts all go quiet |
| `worker` `LoadConsensusAncestryDetail` | Advanced detail report omits the ancient breakdowns |
| `publish::consensus_ancestry_results` | ancient methods are **not federated to the PDS** |

The publish filter matters most: fabricated breakdowns must not reach the network.

---

## 3. Rebuild design

### 3.1 Choose the model first — drop PCA-centroid distance

Distance-to-centroid in a modern PCA is the wrong tool. Two defensible replacements:

**Option A — Supervised admixture over allele frequencies (recommended).**
Model the sample's genotype at each AIM as drawn from a mixture of *K* ancient source populations
with known allele frequencies `f_k(site)`. Fit mixture weights `w` (non-negative, sum to 1) by
maximising the likelihood
`Σ_sites log Σ_k w_k · P(genotype | f_k(site))`
via EM or projected gradient. This is exactly the shape of the existing
`estimate_by_allele_frequency` / `estimate_admixture` code that already works for the super-population
estimate — so the machinery and its validation exist. **It needs ancient *allele frequencies*, not PCA
centroids.**

**Option B — qpAdm-style f4/f-statistics.** Rigorous and the field standard, but it needs outgroups,
a much heavier reference set, and a substantially larger implementation. Not warranted for a
consumer-facing three-way breakdown.

**Recommendation: Option A.** It reuses proven code, needs only a frequency table, and yields genuine
mixture proportions with a residual/fit statistic we can surface as a confidence.

### 3.2 Source populations

Target the classic three-way model, which maps 1:1 onto how the feature is already labelled (and onto
FTDNA's Hunter-Gatherer / Farmer / Metal Age Invader):

- **WHG** — Western Hunter-Gatherer
- **EEF / ANF** — Anatolian Neolithic Farmer
- **Steppe** — Yamnaya / Bronze Age steppe

Optionally a 4th (**CHG** or **Iran_N**) for non-European subjects. Avoid including both `Steppe` and
`EHG`+`CHG`: Steppe *is* approximately EHG+CHG, so they are collinear and the fit becomes
ill-conditioned — one of the failure modes we just diagnosed.

### 3.3 Building the frequency table (the real work)

1. **Source genomes.** Published ancient genomes (Allen Ancient DNA Resource / Reich Lab dataset is
   the standard source). Select the sample sets defining each source population.
2. **Genotypes.** Ancient data is pseudo-haploid and low coverage. Use the published pseudo-haploid
   calls; do **not** attempt diploid calling.
3. **Restrict to our AIM panel sites**, lifted to CHM13 (the panel already carries per-build
   coordinates, so this is a lookup, not a liftover — see `ibd_panel.rs`).
4. **Compute per-site alt-allele frequency per source population**, in the **canonical CHM13
   orientation** (reuse `dosage_from_alleles` so ref/alt swaps and strand flips are handled — the same
   orientation logic already validated for chips and for the b37/b38 alignment path).
5. **Emit** `ancestry_freq_ancient_<build>.bin` in the existing `AncestryPanel` shape
   (`PanelSite { contig, position, reference_allele, alternate_allele, freqs: Vec<f32> }`), with
   `populations = [WHG, ANF, Steppe]`. This is the *same* container the working super/fine panels use.
6. **Drop the ancient PCA asset entirely** — nothing should consume it again.

Build step lives in `navigator-panelbuild` alongside the existing panel builders; publish as a versioned
asset like the others.

### 3.4 Validation gates (do not ship without these)

Empirical, on real data — the diagnosis above only surfaced because we checked real numbers:

1. **Sanity band per region.** A NW-European sample must land roughly in
   Steppe 40–55%, ANF 25–40%, WHG 10–25%. A Sardinian must be ANF-dominant (>60%) with near-zero
   Steppe. A Yoruba must not fit the model at all (large residual → we report "not applicable"
   rather than a number).
2. **Stability.** The same subject genotyped from two different sources (e.g. 28× WGS vs. a chip)
   must agree within a few percent. The current implementation fails this badly (WHG 72.6% vs 7.4%);
   it is the single most diagnostic test.
3. **Density.** Results must be stable as sites are randomly downsampled (say 50%). If they are not,
   the model is over-fit or ill-conditioned.
4. **Fit residual.** Surface it. A poor fit must present as "we can't model this ancestry" rather than
   as confident percentages.

### 3.5 Effort

The code is the small part (a frequency-mixture estimator ≈ the existing `estimate_admixture`, plus a
panelbuild step). **The work is assembling and validating the ancient reference frequencies** —
sourcing the Reich dataset, selecting sample sets, and running the validation gates above.

Until that asset exists and passes §3.4, the feature stays off.

---

## 4. Rebuild attempt (2026-07): Option A built, and where it fails

§3 was implemented in full: the PCA-centroid classifier and the nMonte estimator were deleted;
`estimate_ancient_admixture` (the supervised frequency EM) and a dedicated `ancient-panel` builder
were added; the ancient AF asset `ancestry_freq_ancient_chm13v2.0.bin` was built from the AADR over
WHG (Villabruna-cluster Mesolithic, n≈75) / ANF / Steppe, with a per-source ≥8-called floor. Three of
the four §3.4 gates pass. **Gate 2 (stability) fails**, and the failure is diagnostic of a
method-level problem, not a data-sourcing one.

### 4.1 What passes, what fails

- **Recovery, sanity band, density** — pass *in simulation*. `validate-ancient` draws individuals from
  known reference frequencies and recovers mixtures exactly (20/30/50 → 19.5/30.2/50.3); a simulated
  GBR lands Steppe 50 / ANF 34 / WHG 16; Yoruba/East-Asian are rejected by the dispersion + European
  scope guards.
- **Stability (real data) — FAILS.** Subject `huF98AFD`, genotyped by both means, comes out
  **Steppe 58 / ANF 31 / WHG 10 from his two consumer chips** but **Steppe 80–90 / ANF 15 / WHG 4 from
  his WGS alignments** — a ~25-point disagreement on one person. The chip figures sit in the §3.4
  band; the WGS figures do not. Note the simulated gates cannot catch this: they draw genotypes from
  the same frequency space they fit against, so they are blind to a *mismatch between the sample's
  ascertainment and the panel's*.

### 4.2 The cause is ascertainment, established by elimination

The split tracks exactly one property: whether a site is on the subject's consumer chip. Restricting
either WGS alignment to chip-covered sites reproduces the chip answer (aln #9 full 80.6% Steppe →
∩chip 57.9%, dispersion 2.25 → 1.18); the WGS-only complement carries all the bias (∁chip ≈ 91%).
Everything else was ruled out with direct evidence:

- **Not genotyping.** The clean CHM13 alignments genotype at **99.9%** concordance with the chips at
  shared sites (and keep 99.9% of heterozygotes); they still read ~90% Steppe on the non-chip sites.
  (An *un*clean case turned up here — a GRCh38 alignment served corrupt never-het genotypes — but that
  was a stale pre-`3cf4956` genotype cache, fixed separately by bumping `IBD_PANEL_KIND`; re-genotyped
  it returns to 99.9% and still shows the chip/non-chip split.)
- **Not strand-ambiguity.** Dropping A/T + C/G sites moves Steppe < 3 points.
- **Not ref/alt polarity.** 0 of 19,727 ancient sites are swapped vs the (correctly-oriented) super
  panel.
- **Not aDNA transition damage.** Transversion-only sites still read ~86% Steppe.
- **Not low-MAF noise / the panel's differential zero-rate.** This is the decisive one: a modern-MAF
  floor makes the bias **worse**, not better (MAF ≥ 0.10 → full-WGS 89% Steppe). The bias lives at
  moderate-and-high-MAF non-chip sites, so it is not the sparse-source zero-frequency artifact.

What remains is that the estimate depends on **which ascertainment the sites came from** — the
consumer-chip ascertainment gives the correct, stable, in-band answer (two different chips agree to
2 points, dispersion ~1.2), and the AADR/1240k non-chip ascertainment does not. That dependence is the
textbook failure mode of allele-frequency ADMIXTURE, and it is precisely the problem qpAdm/f-statistics
(the §3.1 "Option B" we dismissed) were designed to be robust against. **§3.1's recommendation was made
without anticipating this; the diagnosis overturns it.**

### 4.3 Two fixes

**Option A′ — restrict the ancient panel to consumer-chip-ascertained sites.** Allele-frequency
admixture is only valid when the sample and the reference share ascertainment; forcing that match is
the *correct* way to run Option A, not a workaround. It already clears the stability gate on real data
(WGS-on-chip-sites = chips = ~58%, in band). Small change: a build-time site filter in `ancient-panel`
intersecting with a canonical consumer-chip manifest. Limits: it drops ~half the panel's sites, needs
that manifest, and — because `huF98AFD` is our only real dual-source subject — the non-European sanity
band can still only be checked in (ascertainment-blind) simulation until more real samples exist.

**Option B — implement qpAdm-style f4/f-statistics.** The principled, ascertainment-robust method,
now justified by the diagnosis despite §3.1's dismissal. Works on the full panel. Larger lift: an
outgroup reference set and an f4 solver.

**Decision (2026-07): pursue A′** as the near-term ship — it is a legitimate fix, it is small, and it
already passes the gate §3.4 calls the single most diagnostic — keeping Option B as the documented
fallback if A′ cannot hold the sanity band as real samples accumulate.

### 4.4 A′ implemented and validated

- **Builder.** `panelbuild ancient-panel` gained `--ascertain-sites <contig\tpos>`: it restricts the
  panel to a consumer-array manifest (mapped to CHM13). Off by default (builds the full panel);
  supplied by the pipeline when a manifest is configured.
- **Pipeline.** `05_build_assets.sh` maps `$CHIP_MANIFEST` (an array's rsID export) to CHM13 via the
  stage-02 1240k liftover and passes the result to the builder; the §3.4 validation step is a hard
  gate. Building without a manifest logs that the deep panel will fail stability.
- **Validation.** A chip-ascertained panel — the AADR ancient panel intersected with the 23andMe v5 +
  AncestryDNA v2 sites (9,971 of 19,727 sites) — clears every gate:
  - *Stability:* met by construction, and already measured — `huF98AFD`'s WGS restricted to array
    sites is 57.9% Steppe vs his chips at 56–58%.
  - *Recovery:* 20/30/50 → 20.2/29.2/50.6.
  - *Sanity band:* GBR 50.0 / CEU 49.6 / FIN 55.1 Steppe, TSI/IBS ANF-shifted, all in band; PJL by
    EUR-scope and CHB/JPT/YRI/LWK by dispersion all **rejected**. Dispersion is if anything *tighter*
    than the full panel (IBS 2.60 vs 3.02).
  - *Density:* GBR 50.0 → 51.8 at half the sites.

**Manifest (done).** `scripts/ancestry-panel/manifests/consumer_array_1240k_rsids.txt.gz` — the
23andMe v5 ∪ AncestryDNA v2 probe sets ∩ AADR 1240k (649,478 rsIDs; array *design*, not genotypes).
`config.sh` defaults `$CHIP_MANIFEST` to it, so a pipeline build is ascertained out of the box; the
whole path is verified (`.gz` → 05 rsID/1240k join → 9,971 ancient sites, the validated panel).

**Remaining before re-enabling** (`ANCIENT_ANCESTRY_ENABLED` stays `false`): (1) rebuild and
re-publish `ancestry_freq_ancient_<build>.bin` from the AADR through the ascertainment floor (the
current shipped `.bin` is still the full, unascertained panel); (2) re-run §3.4 on that rebuilt asset;
(3) flip the flag. The non-European sanity band remains simulation-only until more real dual-source
samples exist; a broader/vendor-neutral ascertainment (e.g. the AADR Human Origins array) is the
natural follow-up if wanted — see `manifests/README.md`.
