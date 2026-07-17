# Ancient ancestry — why it's disabled, and what a correct rebuild takes

**Status:** feature gated off (`navigator_app::ANCIENT_ANCESTRY_ENABLED = false`). The original PCA
implementation was disabled 2026-07-13 for fabricating numbers (§1–§2). Three approaches have since
**failed the real-data stability gate**: raw allele-frequency admixture and a consumer-array
ascertainment floor (§3), and — prototyped and refuted — target pseudo-haploidization (§5.1). The
literature explains why and names the method that should work (§4); the design and the remaining path
(an f4/qpAdm estimator) are §5, and its concrete implementation scope is §7. The feature stays off
until that is built and passes §5.4.

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
   outgroups); **full scope in §7** (the panel already stores the frequencies f4 needs and the
   outgroups are already extracted, so the real cost is the f4/covariance core, not the plumbing). If
   that is judged too heavy for a consumer three-way, the fallback is **chip-only deep ancestry** —
   compute it from a chip-resolved site set for *every* source, WGS included, accepting the smaller
   site count, and validate non-circularly once a second dual-source subject exists.
3. Re-run §5.4 gates on real data end-to-end; only then flip `ANCIENT_ANCESTRY_ENABLED`.

---

## 6. Current status

`ANCIENT_ANCESTRY_ENABLED = false`. The shipped `ancestry_freq_ancient_<build>.bin` is the full,
unascertained panel (the A′ publish was rolled back). Modern/fine admixture stays enabled and is
unaffected. **Four** approaches have now failed on real data: raw frequency-admixture (§3.1),
consumer-array ascertainment (§3.2), target pseudo-haploidization (§5.1), and the **qpAdm f4 estimator
(Lever 2)** — built, tested, and sound in isolation, but the real WGS returns *negative WHG*, and the
caller (GATK4 vs ours) and read technology (short-read vs PacBio) are both ruled out (**§7.9**). The
residual is a structural **capture (ancient sources) vs shotgun (modern target/outgroups) batch
effect** that runs between sources and target, which no target-side method or outgroup projection can
cancel. Remaining levers (§7.9): ancient-*shotgun* sources, chip-only (unproven), or accept deep
ancestry as not deliverable. The Lever-2 scope and results are **§7**.

---

## 7. Lever 2 implementation scope (f4 / qpAdm)

The good news the scoping surfaced: **most of the machinery already exists**, because f4 consumes
exactly what the panel already stores.

### 7.1 What we already have (and why f4 fits the existing asset)

- **The panel format is already the right shape.** `AncestryPanel` carries, per site, a `freqs:
  Vec<f32>` of per-population alt-allele frequencies. f4 is defined on population *frequencies*:
  `f4(A,B;C,D) = mean_site (a−b)(c−d)`. Crucially f4 is **unbiased from pooled frequencies** — the
  estimation noise in each population's `â` is independent across the four slots, so the cross-terms
  vanish in expectation. (This is unlike f2/f3, which need a sample-size hzcorr.) So **no new
  per-genotype asset format is required** — the stored frequencies are sufficient for the f4 point
  estimates. The target contributes its own frequency `g/2 ∈ {0, 0.5, 1}` in one slot.
- **The outgroups are already extracted.** `pops/aadr_component_map.tsv` already maps AADR groups to
  `African` (Mbuti, Yoruba), `EastAsian` (Han, DevilsCave_N), `Oceanian` (Papuan), `AASI` (Onge),
  `NativeAmerican` (Karitiana, Anzick), `ANE` (MA1), plus `CHG`, `Iran_N`, `EHG`. These are the
  classic Lazaridis/Harney "right" set. Today the ancient-panel builder keeps only `WHG,ANF,Steppe`
  and drops the rest; Lever 2 just stops dropping them.
- **The builder already parameterizes the population set.** `ancient-panel --components` takes an
  arbitrary ordered list with a per-source `--min-called` floor. Building the f4 asset is, at the
  matrix level, `--components WHG,ANF,Steppe,Mbuti,Yoruba,Han,Papuan,Onge,Karitiana,MA1,...`.
- **The gate harness exists.** `validate-ancient` (simulated recovery/band/density) and
  `App::ancient_ancestry_stability` / `debug-ancient` (real-data WGS-vs-chip) carry straight over as
  the acceptance tests — they just wrap the new estimator.

### 7.2 The estimand and the solve

qpAdm expresses the target as a weighted mix of the sources, using only allele-sharing measured
*against outgroups* (which is what cancels ascertainment/drift). Sources `S_1..S_n` (n=3:
WHG/ANF/Steppe), outgroups `R_1..R_m`, target `T`. Fix a base source `S_1` and a base outgroup `R_1`.
The admixture identity `T = Σ w_i S_i` implies, for every outgroup `R_j`:

```
  y[j]    = f4(T,   S_1; R_1, R_j)              j = 2..m         (target vs base source)
  A[i][j] = f4(S_i, S_1; R_1, R_j)              i = 2..n         (each source vs base source)
  model:  y = Σ_{i=2}^{n} w_i · A[i][·],   w_1 = 1 − Σ_{i≥2} w_i
```

Solve `(w_2..w_n)` by **GLS** with the block-jackknife covariance of `y` (§7.3), recover `w_1`, and
read the model-fit **χ²/p-value** from the weighted residual with `(m−1)−(n−1)` dof. Accept only when
the model is *not* rejected (p ≳ 0.05) **and** all `w_i ∈ [0,1]`. This is qpWave-rank-1-of-residual /
qpAdm; the exact ADMIXTOOLS index convention will be pinned against Harney et al. 2021 during
implementation, but the substance above is what gets built.

### 7.3 Covariance — the frequency block-jackknife (and its honest caveat)

The GLS weights, the standard errors, and the p-value all need `Cov(y)`. Compute it by **block
jackknife**: partition the genome into ~5 cM (or ~5 Mb) blocks by each site's `(contig,pos)` — data we
already have — recompute the full f4 vector leaving out one block at a time, and form the jackknife
covariance. This is pure arithmetic over the panel frequencies + target dosages; no new asset field.

**Caveat to document, not hide:** textbook ADMIXTOOLS jackknifes over *per-individual* genotypes; we
only retain *pooled* per-population frequencies, so ours is a frequency-level block jackknife — a
legitimate, lighter approximation that gets the point estimates exactly right and the covariance
approximately right. For a consumer three-way estimate this is the right tradeoff, but the SEs it
produces are approximate and the p-value gate must be treated as a screen, backstopped by the
empirical stability gate (§5.4 gate 1), which is the one that actually caught §3.

### 7.4 The outgroup set — the one expertise-driven artifact

As with `aadr_component_map.tsv`, outgroup choice is the judgment call that makes or breaks the
method. Requirements: outgroups must be **differentially related** to the three sources (so `A` has
full rank) and must **not** be downstream of the mixture (no gene flow *from* the sources after they
diverged). Safe distal starting set: **Mbuti, Yoruba, Han, Papuan, Onge, Karitiana, MA1** (+ ancient
Kostenki/Ust'-Ishim if we can source them). `EHG`/`CHG`/`Iran_N` are *ancestral inputs to Steppe* and
are the tempting-but-dangerous case — they can sharpen rank but risk violating the no-downstream-flow
assumption; **hold them out of the initial set** and add only if the stability gate improves with
them. Ship the chosen set as a committed, commented `rightpops`/`leftpops` file, parallel to the
component map. Verify each outgroup clears the `--min-called` floor with enough samples for a usable
frequency (several of the `.DG` groups may be thin — a build-time count check).

### 7.5 Build-pipeline changes

- Extend the ancient asset (or add `f4-panel`, likely just a thin wrapper over `ancient-panel`) so the
  built `AncestryPanel` carries **sources + outgroups** as its `populations`, each with the call
  floor. Source vs outgroup roles live in the **estimator inputs** (two code lists, like
  `AncestryPanel::subset`), not baked into the asset, so the asset stays a plain frequency panel.
- New committed artifact: `pops/qpadm_rightpops.txt` (+ `leftpops`) with provenance comments.
- Asset name: reuse `ancestry_freq_ancient_<build>.bin` (now carrying the extra populations) or mint
  `ancestry_f4_ancient_<build>.bin`; either way add it to the manifest/sha256 set.

### 7.6 Estimator API and integration

- New `navigator_analysis::ancestry::estimate_qpadm(genotypes, panel, sources, outgroups,
  reference_version) -> Option<AncestryResult>`, plus a raw `qpadm_fit` that also returns the
  p-value/SEs (mirroring the `estimate_ancient_admixture` / `ancient_admixture_fit` split so
  `validate-ancient` can report rejected fits).
- **Gate swap:** the acceptance test changes from `fit_distance ≤ ANCIENT_MAX_DISPERSION` to
  `p ≥ P_MIN` **and** `w_i ∈ [0,1]`; carry the p-value in `fit_distance` (or add a field) so the UI's
  "we can't model this ancestry" path still fires on a rejected model. Keep the `west_eurasian_share ≥
  50%` scope guard unchanged.
- **5 call sites, unchanged shape.** `haplogroup.rs` (×2), `brief.rs`, `worker.rs`,
  `validate_ancient.rs` all call `estimate_ancient_admixture(...) -> Option<AncestryResult>`. Point
  them at `estimate_qpadm` (extra `sources`/`outgroups` args from the committed files); the return
  type is identical, so the UI/brief/publish paths need no change. `ANCIENT_ANCESTRY_ENABLED` stays
  the kill switch and stays `false` until §5.4 passes.

### 7.7 Validation (unchanged bar; add the qpAdm-native screens)

Reuse every §5.4 gate. Add: **(a)** model-rejection behaves — Yoruba/non-European rejected by
*p-value* (not just dispersion); **(b)** weights land in `[0,1]` without truncation for a real NW
European; **(c)** stability gate 1 measured **end-to-end through the app consensus**, not `∩chip`
(the §3.2 circularity). Gate 1 remains the single most diagnostic and the only one that failed all
three prior attempts — a green simulation is necessary, never sufficient.

### 7.8 Task breakdown, effort, risk

1. **f4 core + block jackknife** (`ancestry.rs`, pure fn over frequencies + dosages; unit-tested on a
   synthetic admixture graph with known f4). ~The load-bearing piece.
2. **qpAdm GLS solve + rank/p-value** (nalgebra; sources/outgroups as code lists).
3. **Outgroup file + build wiring** (`qpadm_rightpops.txt`, extend `ancient-panel` components, rebuild
   asset, count-check outgroup coverage).
4. **Estimator API + swap the 5 call sites + gate swap** (p-value/weights replace dispersion).
5. **Validation**: `validate-ancient` over the new estimator + the real-data stability gate on
   `huF98AFD` end-to-end. **Only then** flip the flag.

**Risk / kill criteria.** The honest failure mode: even qpAdm can reject or destabilize on a single
diploid present-day target against pseudo-haploid ancient sources (§4.3 warns against exactly this
setup). If, after step 5, `huF98AFD` still splits WGS-vs-chip beyond a few points, the method is
exhausted for our data and the decision is **chip-only deep ancestry** (§5.6 fallback) — which needs a
second dual-source subject to validate non-circularly. Steps 1–2 are the real cost (a correct,
tested f4/covariance core); 3–4 are largely wiring over machinery that exists.

### 7.9 Empirical result (2026-07) — built, tested, and the WGS still won't fit

Steps 1–3 are done and committed. The estimator is **sound in isolation** and the caller is **ruled
out**, but the real WGS target does **not** fit the model. The evidence, in order:

- **The estimator recovers known mixtures.** `qpadm_selftest` draws a target from the panel's own
  `Σ wᵢ·sourceᵢ` frequencies and fits it back: NW-European 20/30/50 → WHG 23 / ANF 22 / Steppe 54
  (feasible), pure-WHG → 97% (feasible). Adding an **ANE outgroup** (MA1 + AfontovaGora + Yana, n≈5 —
  recovered by fixing the stale `Russia_MA1_HG.SG` group-ID to v66.p1's `Russia_Malta_UP` /
  `Russia_AfontovaGora_UP`) tightens the WHG/Steppe axis exactly as intended (WHG SE 21→14). The f4
  core, the GLS solve, the p-value, and the outgroup set all work.
- **The real WGS gives an infeasible fit.** `huF98AFD`'s WGS (aln 3713) → **WHG −34% to −68%**, ANF
  ~54–66, Steppe **80–102%** — a *negative* hunter-gatherer weight. The model is accepted by the
  p-value (p≈0.16) but the weights are outside [0,1]: the sources cannot express this genome.
- **Not the read technology.** Short-read (bwa/markdup) and **PacBio HiFi** (pbmm2, minimal reference
  bias) both give negative WHG (−68% / −46%). Long reads don't fix it.
- **Not our caller.** Re-genotyped the same WGS at the same 17 k sites with **GATK4 HaplotypeCaller**
  (force-called, local reassembly) → **WHG −34.2 / ANF 53.7 / Steppe 80.5**, concordant with our
  native caller (−34.1 / 53.7 / 80.4) to 0.1%. Two independent callers agree; the genotypes are right.

**Diagnosis.** The genotypes are correct and the estimator is correct, yet a real NW-European WGS
needs *negative WHG*. The remaining structural difference is the one the literature names (§4.2–4.3):
our **sources are 1240k-capture ancient pseudo-haploid** data while the **target and outgroups are
modern shotgun**. That capture-vs-shotgun / ancient-vs-modern batch line runs *between the sources and
the target*, so f4-against-outgroups cannot cancel it — the outgroups sit on the target's side of the
line. The self-test doesn't see it only because it draws the target from the capture-derived source
frequencies (same side of the line). This is consistent with every prior failure: frequency-EM (§3.1),
target pseudo-haploidization (§5.1), and now qpAdm all fail on the same WGS, because none removes a
source-vs-target batch effect.

**What this rules out / leaves open.**
- Ruled out: the estimator, the outgroup set, the read technology, and the variant caller.
- Open, and the only remaining levers: (a) **match data types** — rebuild the sources from ancient
  *shotgun* (`.SG`/`.DG`) genomes so sources and target are the same assay class (literature-consistent,
  but few shotgun WHG/ANF/Steppe genomes exist — data-limited); (b) **chip-only deep ancestry** (§5.6),
  though the chip is still modern-vs-ancient and its ~58% was only ever measured through the retired
  frequency-EM, never qpAdm, so it is unproven and possibly circular; (c) accept that deep ancestry is
  **not deliverable** on this data and keep it disabled, shipping only the sound modern/fine admixture.

`ANCIENT_ANCESTRY_ENABLED` stays `false`. Tasks 4–5 (wire the estimator into the app + flip the flag)
are **not** started: there is no point wiring an estimator that returns infeasible weights on real
data. The decision above is a genuine fork and is left to the maintainer.

### 7.10 Correction — the batch line was our *outgroup sourcing*, and it is largely fixable

§7.9 called the batch effect "structural." That was **too pessimistic**: comparing our build to the
standard qpAdm workflow (Reich-lab `ADMIXTOOLS`) exposed a real construction error, and fixing it
moved the real WGS most of the way back.

**The error.** In qpAdm the *right* (outgroup) set must live in the **same genotype callset as the
left (source) set**, so `f4(Left, Left; Right, Right)` is measured in one consistent space and only the
single target is "foreign." We instead sourced outgroups from **modern 1000G/SGDP shotgun** (a separate
download→CHM13-liftover→calling pipeline) while the sources were **AADR 1240k capture** — putting a
capture-vs-shotgun batch line straight through the middle of the f4 matrix. The "strengthen the thin
AADR outgroups from the workspace DB" move was backwards: a thin AADR outgroup (Mbuti n≈15) is
*correct* because it shares the callset; strengthening it from a foreign pipeline breaks qpAdm's core
assumption. We had also never used the canonical deep anchors **Ust'-Ishim** and **Kostenki14** (both
present in our AADR matrix) that resolve the West-Eurasian structure.

**The fix and its effect on `huF98AFD`'s WGS (GATK genotypes):**

| Outgroup set | WHG | ANF | Steppe | model p |
|---|---|---|---|---|
| Wrong-pipeline (1000G/SGDP + ANE) | −34% | 54 | 80 | 0.16 |
| **All-AADR** (Han/Papuan/Yoruba/Karitiana/Mbuti + Ust'-Ishim/Kostenki, deep anchors **pooled** for density) | **−8% (±42)** | 42 | 65 | **0.87** |

From *confidently wrong* (−34 to −68% WHG) to *statistically consistent with the truth* (−8% ± 42
includes the true ≈15%), with a near-perfect fit (p=0.87). Notes from the tuning: single-genome
anchors (Ust'-Ishim n=1) **bias** the frequency-based fit — their idiosyncratic allele-sharing doesn't
average out — so the deep anchors must be *pooled* (UP-Eurasian, ANE) for density; dropping them
entirely leaves the sources unresolved (SE 90%). Pseudo-haploidizing the target did **not** help
(again), so target representation is not the lever.

**What now blocks it: precision, not bias.** WHG comes back with **SE ≈ 40%**, where a real qpAdm run
resolves a Briton to **±2–3%**. ~10× too imprecise, so the point estimate wanders slightly negative
though it is consistent with the true value. The prime suspect is our **pooled-frequency + block-
jackknife-over-sites** approximation with **sparse ancient sources** (WHG ≈ 28 calls/site), versus
ADMIXTOOLS' **per-individual genotypes** with a per-sample jackknife. We could not settle
implementation-vs-data locally: no `ADMIXTOOLS` / R / conda is installed.

**Revised remaining levers.**
- (a) **Per-individual f4 machinery** — persist per-sample ancient genotypes (not just pooled `freqs`)
  and do the ADMIXTOOLS-style per-individual block jackknife. Most likely to recover the precision;
  the bigger lift.
- (b) **Validate against real `qpAdm`/`admixtools2`** (offline) on the same EIGENSTRAT to learn whether
  ±2–3% is even achievable with our sparse sources, or whether it is a data limit. Decisive, needs a
  C/R install.
- (c) **Denser sources** — the WHG source at ~28 calls/site is the sparse bottleneck; higher-coverage
  WHG genomes would tighten the SE.
- (d) Accept "directionally right, imprecise," keep disabled.

The committed `pops/qpadm_rightpops.txt` is updated to the AADR-native set (the methodologically
correct choice regardless of the precision question). `ANCIENT_ANCESTRY_ENABLED` remains `false`.
