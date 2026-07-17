# Ancient ancestry — why it's disabled, and what a correct rebuild takes

**Status:** feature gated off (`navigator_app::ANCIENT_ANCESTRY_ENABLED = false`). The original PCA
implementation was disabled 2026-07-13 for fabricating numbers (§1–§2). Three approaches have since
**failed the real-data stability gate**: raw allele-frequency admixture and a consumer-array
ascertainment floor (§3), and — prototyped and refuted — target pseudo-haploidization (§5.1). The
literature explains why and names the method that should work (§4); the design and the remaining path
(an f4/qpAdm estimator) are §5. The feature stays off until that is built and passes §5.4.

**Scope:** the evidence that the original is unsound (§1–2), what the rebuild attempts established and
why they fail (§3–4), and the design for a correct rebuild (§5–6).

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

## 3. What the rebuild attempts established (2026-07)

Both attempts were built end to end and tested; both fail. The failures are informative, so they are
the evidence base for the §5 design.

### 3.1 Attempt 1 — supervised allele-frequency admixture

The design originally chosen (the "Option A" of earlier drafts): model the sample's genotype at each
AIM as a mixture of source allele frequencies `f_k(site)` and fit simplex weights by EM — the same
machinery as the working super-population estimate, over a dedicated ancient frequency panel built
from the AADR (WHG / ANF / Steppe, per-source ≥8-called floor, no no-data-as-zero). Implemented as
`estimate_ancient_admixture`; asset built by `panelbuild ancient-panel`; gates by `validate-ancient`.

It **passes in simulation** — draw individuals from known frequencies and it recovers mixtures
exactly (20/30/50 → 19.5/30.2/50.3), lands a simulated GBR at Steppe ≈ 50, rejects non-Europeans. But
simulation draws genotypes from the same frequency space it fits against, so it is **blind to the
failure that matters**:

**It fails the stability gate on real data.** Subject `huF98AFD`, genotyped by both means, comes out
**Steppe 58 / ANF 31 / WHG 10 from his chips** but **Steppe 80–90 / ANF 15 / WHG 4 from his WGS**.
Same person, ~25 points apart. The chip figures sit in the expected band; the WGS figures do not.

The split tracks exactly one thing — the *data type* of the genotypes — and everything intrinsic was
ruled out with direct evidence:

- **Not genotyping.** Clean CHM13 alignments genotype at **99.9%** concordance with the chips at
  shared sites (99.9% of hets preserved) yet still read ~90% Steppe on non-chip sites. (A separate
  never-het corruption turned up on one GRCh38 alignment; it was a stale pre-`3cf4956` genotype cache,
  fixed by bumping `IBD_PANEL_KIND` — re-genotyped it returns to 99.9% and the split remains.)
- **Not strand-ambiguity.** Dropping A/T + C/G sites moves Steppe < 3 points.
- **Not ref/alt polarity.** 0 of 19,727 sites swapped vs the correctly-oriented super panel.
- **Not aDNA transition damage.** Transversion-only sites still read ~86% Steppe.
- **Not low-MAF / the panel's differential zero-rate.** A modern-MAF floor makes it *worse* (MAF ≥
  0.10 → 89% Steppe), so it is not the sparse-source zero-frequency artifact.

### 3.2 Attempt 2 — restrict the panel to consumer-array sites (A′)

Reasoning: allele-frequency admixture is only valid where sample and reference share ascertainment, so
restrict the panel to sites consumer arrays actually assay. Built as `ancient-panel --ascertain-sites`
+ a committed manifest (`manifests/consumer_array_1240k_rsids.txt.gz`, 23andMe v5 ∪ AncestryDNA v2 ∩
AADR 1240k). This looked like it worked, and **the "pass" was a mistake worth recording**: it rested
on the `∩chip` figure — WGS restricted to the subject's *own resolved chip calls* — which trivially
equals the chip because they are the same sites genotyped two ways.

The **end-to-end smoke test refutes it**: rebuild the ascertained asset (9,971 sites), publish it,
re-genotype the subject, fit the *actual app consensus* → **consensus 75.2% Steppe vs chips ~58%**.
The gate still fails. Restricting further to *called* (non-no-call) array sites barely moves it
(75 → 71). The residual bias is the array-manifest sites the chip path drops but a WGS caller
genotypes anyway (`∁chip`: 99% Steppe at dispersion 8.6 — the model fits them terribly). A coarse
site filter cannot isolate the well-behaved sites, and the only set that gives 58% for WGS is the
subject's own chip set, so validating it on our one dual-source subject is circular.

**Takeaway:** this is not a data-sourcing problem and not a site-filter problem. It is the method.

---

## 4. Why it fails — the literature

The instability is a known, named phenomenon, and the field long ago moved off the class of method we
were using. Four load-bearing findings:

1. **Frequency-mixture ADMIXTURE is ascertainment-sensitive by construction; f4/qpAdm is the
   standard because it is not.** The field estimates ancient-source proportions with **qpAdm**, which
   works on **f4-statistics** — allele-sharing measured *against a panel of outgroups* — not on the
   sample's raw frequencies. Differences-of-differences against outgroups cancel drift and are robust
   to SNP ascertainment (Harney et al. 2021; f4/admixture-graph reviews). This is exactly the property
   Attempts 1–2 lack.

2. **Our WGS-vs-chip split is a documented capture-vs-shotgun batch effect.** "Capture samples had
   higher attraction to the Anatolian Neolithic than shotgun samples," and that attraction shifts with
   hunter-gatherer admixture (*Testing times*, Genetics 2024). A subject's ANF/Steppe split moving
   with the data type — rather than with their ancestry — is precisely this.

3. **Do not co-analyze present-day and ancient data.** Harney et al. state it directly, and warn
   against mixing capture and shotgun in one model. Our feature does exactly this (a present-day
   individual against ancient sources), so even the gold-standard method flags the setup as hazardous.

4. **Reference bias + pseudo-haploid processing create spurious allele-sharing.** Ancient pseudo-
   haploid calls carry reference bias that "systematically creates spurious signals of allele sharing"
   in D/f-statistics (Günther & Nettelblad 2019). Our ancient references are pseudo-haploid-with-
   damage; our modern target is diploid-called — that mismatch is a documented source of "attraction."

And the original PCA/`nMonte` approach (§1) is a community distance method (Global25 + nMonte) known
to overfit and "not ideal for recent ancestry" — abandoning it was correct.

Sources: Harney et al., *Assessing the performance of qpAdm* (Genetics 2021,
`academic.oup.com/genetics/article/217/4/iyaa045`); *Testing times…* (Genetics 2024,
`.../228/1/iyae110`); Soraggi et al., f-statistics biased under all ascertainment schemes (PLOS Gen
2023, `journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010931`); Günther & Nettelblad,
reference bias (PLOS Gen 2019, `.../journal.pgen.1008302`); MetaGLIMPSE imputation
(`pmc.ncbi.nlm.nih.gov/articles/PMC12262289/`).

---

## 5. The rebuild design

Two levers, applied **in order**. The first is cheap, method-agnostic, and attacks the root cause the
literature identifies; the second is the correct estimator if the first is not sufficient.

### 5.1 Lever 1 — harmonize the target to the references (do this first)

Every finding in §4 reduces to one rule: **put the target and the references through one identical
pipeline**, so no data-type batch effect can exist. Concretely, reduce *every* subject — WGS or chip —
to the same representation as the AADR ancient references:

- **Pseudo-haploidize the target.** Sample a single allele (one read for WGS; one of the two reported
  alleles for a chip) per site, instead of feeding diploid dosages against pseudo-haploid references.
  Removes the diploid-vs-pseudo-haploid attraction.
- **Transversions only.** Drop A↔G, C↔T so post-mortem C→T/G→A damage in the references can't create
  attraction.
- **One fixed, orientation-checked site set**, identical for target and references, on CHM13.

**Prototype result (2026-07): Lever 1 is refuted.** `panelbuild --example phaploid_fit` genotyped a
clean CHM13 alignment at the ancient sites with allele depths and fit it read-level pseudo-haploid
(one allele per site, P(alt)=alt_depth/depth). It does **not** move the WGS estimate:
diploid 80.5% Steppe → pseudo-haploid **79.9%** (the transform changed the genotypes — dispersion
2.15 → 3.31 — but not the answer), still nowhere near the chip's ~58%. Transversions make it *worse*
(86–87%). So the split is **not** a target-ploidy batch effect: with the target now processed
identically to the pseudo-haploid references, it persists. What remains is the capture-vs-shotgun /
reference-side effect (§4.2/§4.4), which target-side transforms cannot touch — it points at Lever 2.

### 5.2 Lever 2 — estimator: qpAdm-style f4 (now the primary path, Lever 1 having failed)

Replace the frequency-mixture EM with an **f4-statistics** estimator (the qpAdm approach):

- **Left (sources):** the WHG / ANF / Steppe references. **Right (outgroups):** a diverse panel that
  is differentially related to the sources (e.g. an African outgroup such as Mbuti plus a spread of
  Eurasian/other populations) — the reference already carries enough populations to seed this.
- Estimate weights from the f4 covariance; accept by the standard criteria (model not rejected at
  p ≈ 5%, weights in [0,1]); use **rotating outgroups** rather than p-value ranking to choose among
  models.
- Robustness the literature documents and we need: works down to ~10k SNPs, tolerates high missing
  data and pseudo-haploid genotypes, and is ascertainment-robust — the property §3 lacked.

This is the larger lift (an f4/covariance solver + a curated outgroup set) and a genuine departure
from the "reuse `estimate_admixture`" plan, but it is the method the field actually trusts.

### 5.3 Source populations

The classic three-way model, which maps 1:1 onto the existing labels (and FTDNA's
Hunter-Gatherer / Farmer / Metal-Age-Invader):

- **WHG** — Western Hunter-Gatherer (Villabruna cluster; strengthened to n ≈ 75)
- **ANF / EEF** — Anatolian Neolithic Farmer
- **Steppe** — Yamnaya / Bronze-Age steppe

Optionally a 4th (**CHG** or **Iran_N**) for non-European subjects. Do **not** list `Steppe` alongside
`EHG`+`CHG`: Steppe ≈ EHG+CHG, so they are collinear and the fit becomes ill-conditioned.

### 5.4 Validation gates (do not ship without these)

Empirical, on real data — §3 only surfaced because we checked real numbers, and simulation is blind
to the failure that matters.

1. **Stability — the single most diagnostic test.** The same subject from two sources (WGS vs chip)
   must agree within a few percent. This is what both attempts fail. Measure it end-to-end through the
   app consensus (`debug-ancient`), **not** via a `∩chip` restriction — that restriction is circular
   and produced the false "pass" in §3.2.
2. **Pseudo-haploid consistency.** With Lever 1 in place, the diploid and pseudo-haploidized fits of
   the same WGS must agree — if they don't, the harmonization is incomplete.
3. **Sanity band per region.** NW-European ≈ Steppe 40–55 / ANF 25–40 / WHG 10–25; Sardinian
   ANF-dominant, near-zero Steppe; Yoruba does not fit (large residual → "not applicable").
4. **Density.** Stable under 50% site downsampling.
5. **Fit residual surfaced.** A poor fit presents as "we can't model this ancestry," never as
   confident percentages.

Non-European sanity remains simulation-only until more real dual-source subjects exist; `huF98AFD` is
currently the only one, which is also why the stability gate cannot be self-validated by ascertaining
on his own chips (§3.2).

### 5.5 What's already built and reusable

- `panelbuild ancient-panel` (+ `--ascertain-sites`) and the committed array manifest — the panel
  builder and site-restriction machinery, reusable for the harmonized site set.
- `panelbuild validate-ancient` — the simulated gates (recovery / band / density); keep, but treat as
  necessary-not-sufficient (it passed on both failed attempts).
- `debug-ancient` + `App::ancient_ancestry_stability` — the real-data stability diagnostic with the
  ∩chip/∁chip/density probes; this is the harness for gate 1 and the Lever-1 experiment.
- The AADR component map (`scripts/ancestry-panel/pops/aadr_component_map.tsv`) and WHG strengthening.

### 5.6 Effort and sequencing

1. ~~Pseudo-haploidize the target (Lever 1)~~ — **done, refuted** (§5.1): it doesn't move the WGS
   estimate, so target harmonization is not the fix.
2. **Next:** build the f4/qpAdm estimator (Lever 2) with a curated outgroup set — the ascertainment-
   robust method the diagnosis and the prototype both point to. Larger lift (an f4/covariance solver +
   outgroups). If that is judged too heavy for a consumer three-way, the fallback is **chip-only deep
   ancestry** — compute it from a chip-resolved site set for *every* source, WGS included, accepting
   the smaller site count, and validate non-circularly once a second dual-source subject exists.
3. Re-run §5.4 gates on real data end-to-end; only then flip `ANCIENT_ANCESTRY_ENABLED`.

---

## 6. Current status

`ANCIENT_ANCESTRY_ENABLED = false`. The shipped `ancestry_freq_ancient_<build>.bin` is the full,
unascertained panel (the A′ publish was rolled back). Modern/fine admixture stays enabled and is
unaffected. Three approaches have now failed the real-data stability gate: raw frequency-admixture
(§3.1), consumer-array ascertainment (§3.2), and target pseudo-haploidization (§5.1). The gate to
re-enable is §5.4 gate 1 passing end-to-end. Next concrete step: **§5.6 step 2 — the f4/qpAdm
estimator (Lever 2)**, the only remaining approach the literature supports, or the chip-only fallback.
