# Ancient ancestry — why it's disabled, and what a correct rebuild takes

**Status:** feature gated off (`navigator_app::ANCIENT_ANCESTRY_ENABLED = false`), 2026-07-13.
**Scope of this doc:** the evidence that the current implementation is unsound, and the design for a
replacement.

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
