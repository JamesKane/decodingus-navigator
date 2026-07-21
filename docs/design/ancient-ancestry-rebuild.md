# Ancient ancestry — why it's disabled, and what a correct rebuild takes

**Status:** feature still gated off (`navigator_app::ANCIENT_ANCESTRY_ENABLED = false`), but the method
is now **validated end to end** (§7.14) — the long investigation in §3–7 is largely superseded by that
result, kept for the record. The original PCA implementation was disabled 2026-07-13 for fabricating
numbers (§1–§2); several rebuild attempts then appeared to fail (§3–§7.12), but those "walls" were an
outgroup-sourcing mistake (§7.10), 16 k-SNP imprecision (§7.11), and finally a genotype-labelling bug
in a diagnostic example (§7.13). With the bug fixed, the full 1240k SNP set, and the **Patterson 2022
sister-outgroup qpAdm config**, `huF98AFD`'s real WGS resolves to **WHG 15 / EEF 45 / Steppe 41,
±1–2 %, model accepted** — a literature-grade British breakdown, reproduced identically by both
`admixtools2` and our own `qpadm_fit` (§7.14). What remains before enabling: the WGS-vs-chip stability
gate on this config, and wiring the full-1240k genotyping + `qpadm_fit` into the app path (tasks 4–5).
The estimator design is §5; the implementation + the working config are §7 (start at §7.14).

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
unaffected. The **qpAdm f4 estimator (Lever 2) now works end to end in the offline harness (§7.13)**:
at the full 1240k, `huF98AFD`'s WGS fits **WHG 0 / ANF 57 / Steppe 43, ±6 %, model accepted (p=0.70)**,
matching a 1000G-British AADR control — and a same-person concordance check (our genotyping vs the
AADR's own `HG00160`) is **99.84 %**. Getting here cleared a chain of false walls: not the caller
(GATK4 ≈ ours, §7.9), not a structural batch effect (that was our outgroup *sourcing* — fixed with
AADR-native outgroups, §7.10), not our estimator (`admixtools2` reproduces it, §7.11), and finally not
"target genotyping" either — §7.12 blamed the target, but that was a **genotype-labelling bug** in the
`genotype_bed` example (§7.13). The two real levers turned out to be **SNP count** (16 k AIM → full
1.15 M 1240k drops WHG SE ~8×) and **fixing that bug**. Remaining before enabling: the WGS-vs-chip
stability gate, a WHG≈0 source-config refinement, and wiring the full-1240k asset into the app
(tasks 4-5). Lever-2 scope + results are **§7**.

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

### 7.11 Decisive: real `admixtools2` gives the SAME huge SE — it's the SNP count, not our code

To settle implementation-vs-data (§7.10), installed **`admixtools2` 2.0.10** (R) and ran genuine
`qpadm` on the *identical* data: built an EIGENSTRAT of the 446 AADR source+outgroup individuals plus
the `huF98AFD` target (CHM13, our ~20 k panel sites, GATK genotypes; 15 961 SNPs after admixtools'
own filtering), same left (WHG/ANF/Steppe) and right (Han/Papuan/Yoruba/Karitiana/Mbuti/UPEuro/ANE).

Result — real qpAdm:

| left | weight | SE | z |
|---|---|---|---|
| WHG | 0.4% | **±43%** | 0.01 |
| ANF | 60.1% | ±28% | 2.17 |
| Steppe | 39.5% | ±40% | 1.00 |

Model **accepted** (χ²=2.66, p=0.62). **This is our own estimator's answer** (WHG ≈ 0 with SE ≈ 40%,
model accepted; ours gave WHG −17 ± 33 on the pooled panel). ADMIXTOOLS is marginally better-behaved
(feasible weights vs our slightly-negative point estimate) but its **SEs are just as large**. So:

- **Our f4/qpAdm implementation is validated** — it reproduces ADMIXTOOLS on the same input. Tasks 1–2
  were built correctly. The "per-individual machinery" lever (§7.10a) is **not** the fix — the pooled-
  frequency approximation is not what limits precision.
- **The precision limit is the SNP count.** We ran on **~16 k SNPs** — our high-Fst *AIM panel* subset.
  A real 1240k qpAdm uses **~1.15 M SNPs**, ~70× more. SE scales ≈ 1/√(n_SNP), so full 1240k is a
  ~8× SE reduction: WHG ≈ ±5% instead of ±43% — enough to resolve a Briton's ~15% WHG (the published
  studies get ±2–3% precisely this way). Our whole asset pipeline was built on a 20 k AIM panel
  (sized for the modern super-pop AF work); deep ancestry needs the full 1240k.

**Revised path — bounded and hopeful.** Deep ancestry is achievable with our validated estimator; it
needs the **full-1240k SNP set**, not the 20 k AIM subset:
1. Re-extract the AADR sources + AADR-native outgroups at the **full 1240k** (the AADR `.geno` already
   is 1240k; the pipeline currently down-selects to 20 k AIMs — build a *separate* dense ancient asset
   instead of subsetting).
2. Genotype the target at those ~1.15 M sites (CHM13), not just the AIM panel.
3. Re-run the (validated) `qpadm_fit`; expect WHG SE ≈ ±5%. Then §5.4 gates + flip.

This is a data-volume task, not a research dead-end. `ANCIENT_ANCESTRY_ENABLED` stays `false` until the
full-1240k rebuild passes the stability gate. The `admixtools2` cross-check (`examples`-adjacent R
script in the build scratch) is the reusable oracle for validating our estimator on any future asset.

### 7.12 Full-1240k rebuild — SNP count confirmed; the blocker is now our target genotyping

Did the rebuild. Genotyped the target at all **1,231,970** CHM13 1240k sites (our caller, 1.23 M called),
subset the AADR PLINK to the 446 source+outgroup individuals (`plink2`), reconciled the target's
alleles onto the AADR SNPs, and ran real `admixtools2` qpAdm on the combined ~1.15 M-SNP set. Added a
**positive control**: 12 **English (1000G British, `HG00xxx.DG`)** individuals already in the AADR —
the same population as `huF98AFD`, but genotyped through the AADR pipeline instead of ours.

| Target (full 1240k) | WHG | ANF | Steppe | SE (WHG) | model |
|---|---|---|---|---|---|
| **English (AADR pipeline)** | 0.0% | 50.1% | 49.8% | **±2.6%** | feasible weights, p=3e-4 |
| **huF98AFD (our CHM13 genotyping)** | −346% | +451% | −5% | ±20% | infeasible, p=2e-13 |

Two firm conclusions:

- **The SNP-count precision hypothesis is CONFIRMED.** The English control's SE collapsed to **±2.6%**
  (from ±43% at 16 k) — exactly the ~8× the arithmetic predicted, and the ~2–3% real studies report.
  Full 1240k delivers publication-grade precision, and the whole stack (our panel, the AADR-native
  outgroups, ADMIXTOOLS) resolves a real Briton feasibly and tightly.
- **The remaining blocker is entirely our target genotyping.** Identical sources/outgroups/SNPs/
  estimator: a real Briton *from the AADR callset* fits; the *same* person genotyped by us from a CHM13
  BAM gives −346% WHG. Flipping the target globally makes it worse (−5.4%), so it is **not** a uniform
  orientation bug — it is a *per-SNP* inconsistency between our CHM13-lifted genotyping and the AADR
  1240k allele space (the §3.1 "non-AIM sites read anomalously" finding, now proven at full scale: the
  16 k AIM panel is clean, the other ~1.13 M 1240k sites are not, for our target).

**Next step (target-side, bounded).** Genotype the target in a *standard 1240k build that matches the
AADR* rather than via the hg19→CHM13 liftover: the target has a local **GRCh38** BAM
(`WGS229.b38.bam`) and `1240k_sites.hg38.bed` exists, so genotype at hg38 1240k sites and key to the
AADR by rsID (hg38↔hg19 is the clean mapping 1000G itself uses — no exotic T2T liftover). If a
GRCh38-genotyped `huF98AFD` then fits like the English control, the feature is shippable; if it still
fails, the target-side genotyping/reference-bias problem is the true wall.

> **⚠️ §7.12 below reached a WRONG conclusion — see §7.13.** The "target genotyping is the wall"
> finding was a bug in the `genotype_bed` example (genotypes mis-labelled to rsIDs), not a real batch
> effect. With the bug fixed, the full-1240k qpAdm **works**: `huF98AFD` fits at WHG 0 / ANF 57 /
> Steppe 43, ±6 %, model accepted. §7.12 is kept for the record; §7.13 is the correction.

**Done — and GRCh38 fails identically. It is not the reference build.** Genotyped the target from the
**GRCh38** BAM at the hg38 1240k sites (1.15 M autosomal, 1.15 M called) and re-ran the same qpAdm:
**WHG −3.44, ANF +4.48, Steppe −0.04, model rejected (p=3e-13)** — essentially the CHM13 result
(−3.46 / +4.48 / −0.05) to two decimals. So the failure is **build-independent**: neither the T2T
liftover nor a strand/orientation bug (a global flip made it worse). The reconciliation is sound (the
same pipeline on the 16 k AIM subset gives the feasible-but-imprecise fit of §7.11, and the English
control on the *same* full-1240k SNPs fits at ±2.6 %). What is left is intrinsic: **our externally
genotyped modern WGS is systematically inconsistent with the AADR 1240 k data at the ~1.13 M non-AIM
sites, however we call it** — the high-Fst AIM sites we curated are clean, the rest are not. This is
the core "do not co-analyze data from different pipelines" wall (§4.3–4.4): the English control fits
*only because it lives inside the AADR callset*, processed identically to the sources. Any target we
genotype ourselves sits on the other side of that line.

**Where deep ancestry stands.** Every component is now proven: estimator (matches ADMIXTOOLS), sources,
AADR-native outgroups, and SNP-count precision (±2.6 % on a real Briton). The one unsolved piece is
target harmonization — making an externally sequenced genome statistically indistinguishable from an
AADR-processed one at all 1240k sites, not just the curated AIMs. That is a genuine research problem
(reference-bias correction / imputation onto the 1240k in the AADR's own representation), not a wiring
task. Absent it, the honest options are: **(a)** ship deep ancestry only for samples we can place on
the AADR's side of the line (none, currently), **(b)** accept the AIM-only fit (precise to ±40 %,
i.e. not shippable), or **(c)** keep it disabled. The decisive same-sample test — genotype a 1000G
individual that *is* in the AADR (`PRJEB36890`/`PRJEB31736` CRAMs) through our pipeline and diff
against its AADR genotypes — needs the NAS mounted; it would quantify exactly how far our calling
drifts from theirs. `ANCIENT_ANCESTRY_ENABLED` stays `false`.

**Secondary tuning note.** Even the English control comes back **WHG ≈ 0** and the model is technically
rejected (p=3e-4). A textbook British qpAdm usually retains ~10–20% WHG, so our specific source set
(the Villabruna WHG, the Yamnaya Steppe which already carries EHG/WHG-like ancestry) and outgroups are
not yet the canonical configuration — a real but *secondary* refinement, dwarfed by the target-
genotyping blocker above. Everything remains behind `ANCIENT_ANCESTRY_ENABLED = false`.

### 7.13 CORRECTION — it was a genotype-labelling bug, and full-1240k qpAdm WORKS

§7.12's "target genotyping is the wall" was **wrong**, and the way it was caught is the point of a
positive control. Copied a **1000G British genome that is itself in the AADR** (`HG00160`,
`PRJEB31736` CHM13 CRAM) locally and genotyped it through *our* pipeline, then compared to the AADR's
own `HG00160.DG` genotypes — the same person, two pipelines.

- First pass: **49.4 % concordant** — i.e. *random* (≈ the by-chance rate for this genotype
  distribution). No real genotyping is that bad; it meant the genotypes were **mis-labelled**.
- Root cause: `genotype_sites_all_contigs` returns genotypes **reordered** (per-contig rayon), so the
  `genotype_bed` example's `zip(gts, rsids)` paired each genotype with the *wrong* rsID — only ~14 %
  landed on the right one. The `(contig,pos,dosage)` on each row were correct; only the rsID column
  was scrambled. **All the full-1240k target files (§7.11–7.12) were keyed by that scrambled rsID.**
- Re-keyed by `(contig,pos)→rsID` from the BED: **99.84 % same-person concordance.** Our CHM13
  genotyping is correct. (Fixed in `examples/genotype_bed.rs`: key each returned genotype to its rsID
  by position, never by input order.)

Re-ran the full-1240k qpAdm for `huF98AFD` with the corrected keying:

| Target (full 1240k, corrected) | WHG | ANF | Steppe | SE | model |
|---|---|---|---|---|---|
| **huF98AFD (our WGS)** | −0.3 % | 57.3 % | 43.0 % | **±6 %** | **accepted, p=0.70** |
| English (AADR control) | 0.0 % | 50.1 % | 49.8 % | ±2.6 % | feasible |

**Feasible weights, ~6 % SE, model accepted — and it matches the English control.** So the entire
"batch effect / can't co-analyze" wall of §7.9–7.12 was an artifact of the labelling bug on top of the
16 k-SNP imprecision. With the bug fixed and the full 1240k, **the whole stack works end to end**:
validated estimator, AADR-native outgroups, ~6 % precision, and an externally sequenced WGS that fits
like a reference genome.

**Status flip.** Deep ancestry is now *technically working* in the offline harness. What remains before
`ANCIENT_ANCESTRY_ENABLED = true`:
1. **Stability gate (§5.4-1)** — run `huF98AFD`'s *chip* through the same full-1240k qpAdm and confirm
   WGS-vs-chip agree within a few percent (the WGS now gives WHG 0 / ANF 57 / Steppe 43).
2. **WHG ≈ 0 refinement** — both the target and the English control return WHG ≈ 0, where a textbook
   Briton keeps ~10–20 %. Likely the source config (Villabruna WHG; Yamnaya Steppe already carries
   EHG/WHG). A real refinement, but now the *only* substantive modelling question left, and the model
   is a valid accepted fit regardless.
3. **Wire the full-1240k asset into the app path** (tasks 4-5): the app currently genotypes the ~20 k
   AIM panel; deep ancestry needs the ~1.15 M-SNP genotyping + the validated `qpadm_fit`.

The `admixtools2` harness (`scripts/ancestry-panel/qpadm_validate.R`) + the local `HG00160` control are
the reusable oracles. `ANCIENT_ANCESTRY_ENABLED` stays `false` until gate 1 passes end-to-end.

### 7.14 The shipping config — Patterson 2022 sister-outgroups resolve WHG

§7.13's fix left one artifact: both the target *and* the English control returned **WHG ≈ 0** where a
textbook Briton keeps ~15–20 %, because our generic continental outgroups couldn't separate WHG from
the **EHG inside Steppe** (Yamnaya ≈ 50 % EHG, and EHG ≈ WHG + ANE). The fix is the qpAdm design from
**Patterson et al. 2022** (*Nature*, "Large-scale migration into Britain") — each outgroup is a
**sister of one source**:

| | source (Left) | sister outgroup (Right) |
|---|---|---|
| Hunter-gatherer | WHG (France_Mesolithic + Loschbour) | **Iron Gates HG** (Serbia/Romania) |
| Farmer | EEF (**Balkan** Neolithic, min-HG) | **Anatolia_N** (Turkey_N) |
| Steppe | Yamnaya + Poltavka | **Afanasievo** |
| base | — | **ancient** sub-Saharan (Malawi LSA + Mota) |

The sister structure gives qpAdm the differential relatedness to tease the three apart. Two of our
earlier choices invert: Anatolia_N was the *farmer source*; here it is the farmer *outgroup* (EEF =
Balkan_N is the source), and Iron Gates HG — excluded from the WHG *source* — is the WHG *outgroup*.
The African base is deliberately *ancient*, capture-processed (§4.4), not present-day Mbuti.

**Result on `huF98AFD`'s full-1240k WGS — WHG resolves, and our estimator matches ADMIXTOOLS:**

| | admixtools2 | **our `qpadm_fit`** | academic British | FTDNA |
|---|---|---|---|---|
| WHG | 14.6 % | **14.6 %** | ~15–20 | 43 |
| EEF/Farmer | 44.6 % | **44.8 %** | ~30–35 | 45 |
| Steppe | 40.8 % | **40.6 %** | ~45–50 | 12 |
| SE | 1–2 % | 1–2 % | — | — |
| model | accepted p=0.15 | accepted p=0.21 | — | — |

Identical to 0.2 % between ADMIXTOOLS and our own pooled-frequency estimator, tight SEs, model
accepted, WHG in the academic range. (FTDNA agrees on Farmer but is the off-consensus outlier on the
Steppe/HG axis — its 12 % Steppe contradicts the qpAdm literature, which puts British at ~45–50 %.)

**This is now the committed shipping config**: `pops/qpadm_leftpops.txt` (WHG/EEF/Steppe),
`pops/qpadm_rightpops.txt` (AnatoliaOG/Afanasievo/IronGates/African), and the group→component map in
`aadr_component_map.tsv`. Both the `admixtools2` oracle (`qpadm_validate.R`) and our estimator
(`examples/verify_qpadm_fit.rs`) reproduce it. The deep-ancestry method is **validated end to end**.
Remaining before `ANCIENT_ANCESTRY_ENABLED = true`: the WGS-vs-chip stability gate on this config, and
wiring the full-1240k genotyping + `qpadm_fit` into the app path (tasks 4-5). New component codes
(`EEF`, `AnatoliaOG`, etc.) will need catalog entries in `navigator-domain::ancestry` for display.

### 7.15 Stability gate PASSES — WGS vs chip agree to 0.3%

The one gate every earlier attempt failed (§5.4-1): the same person, genotyped two independent ways,
must agree. Ran `huF98AFD`'s 23andMe chip (606 k sites ∩ 1240k) through the *identical* Patterson
config + our `qpadm_fit`:

| component | WGS (full 1240k) | CHIP (23andMe) | Δ |
|---|---|---|---|
| WHG | 14.6 % | 14.3 % | 0.3 |
| EEF | 44.8 % | 44.9 % | 0.1 |
| Steppe | 40.6 % | 40.8 % | 0.2 |
| model | accepted p=0.21 | accepted p=0.16 | |

**WGS and chip agree to ≤0.3 % on every component.** The original frequency-EM (§3.1) had this exact
subject at 80 % Steppe from WGS vs 58 % from chip — a ~22-point split that defined the whole problem.
It is now ~0.3 points. The method is stable across data sources.

**Deep ancestry is fully validated on real data:** validated estimator (= ADMIXTOOLS), 99.84 %
same-person genotyping concordance, Patterson sister-outgroups resolving WHG to the academic ~15 %,
±1–2 % precision, our `qpadm_fit` = the oracle, and now the WGS-vs-chip stability gate passing to 0.3 %.
The remaining work is purely engineering — wiring the full-1240k genotyping + `qpadm_fit` + a persisted
1240k ancient asset into the app path, and catalog entries for the new component codes (tasks 4-5).
`ANCIENT_ANCESTRY_ENABLED` can flip once that path is in place and re-checks this gate end-to-end.

### 7.16 Panel orientation + the frontend/backend architecture (multi-ref + chip)

Wiring deep ancestry to GRCh37/38 alignments and chips surfaced how the ancestry assets actually fit
together — and a build bug in the qpAdm asset.

**The architecture is a frontend/backend split.** The `IbdPanel` asset is the full-1240k
(1,231,935 sites) **multi-build genotyping frontend**: it carries per-build loci (100% GRCh37, 100%
GRCh38, + CHM13) and rsIDs, with `resolve_alignment(build, …)` / `resolve_chip(build, …)` that re-key
any-build WGS *and* chips to canonical CHM13. It builds the **autosomal consensus** (`DiploidProfile`)
pooling every source. The frequency panels (`AncestryPanel`: super-pop 20k AIM, fine, qpAdm 1240k) are
**CHM13-only scoring backends**. Modern super-pop + fine ancestry already work across all three
references and chips *because they score the IBD-resolved consensus*. Deep ancestry was CHM13-WGS-only
only because `estimate_deep_ancestry` **bypassed the frontend** and genotyped one alignment directly at
the qpAdm sites.

**The orientation bug.** REF at a site is fixed by the reference base, so on one reference two panels
cannot legitimately disagree — yet the qpAdm asset was ref/alt-swapped vs the IBD panel at ~30% of
shared sites (339,297 of 1,149,621). Cause: `06_build_qpadm_panel.sh` read allele labels from the
1240k **CHM13 BED**, but that BED came from lifting the sites **hg19→CHM13** — the liftover moved the
*position* to CHM13 but kept the *hg19 allele labels*. At the ~30% of sites where CHM13 carries the
hg19-ALT base, REF/ALT are labelled backwards. Proven three ways: at every swapped site the actual
CHM13 FASTA base equals the **IBD** REF (not qpAdm's); the super/fine/IBD assets are mutually
0-swapped (all built from native-CHM13 VCFs, REF = reference base); only the lifted-BED qpAdm asset is
offset. It did **not** break the fit (14.6/44.8/40.6) because it was self-consistent — source freqs
and target dosages both counted the same (mis-labelled) ALT, and f4 is invariant under a consistent
flip — but the wrong labels break any cross-panel join.

**Fix.** `panelbuild ancient-panel --reference <chm13.fa>` orients every site so `reference_allele` is
the actual CHM13 base (swap ref↔alt and each freq→1−freq where reversed; drop sites matching neither).
Stage 6 passes the CHM13 FASTA. The rebuilt asset is CHM13-canonical (0-swapped vs the other assets)
and the fit is unchanged. There is **no** systemic "no canonical orientation" gap — it was one build
bug in the one asset built off a lifted sites file.

**Consequence.** With the qpAdm asset canonical, `estimate_deep_ancestry` consumes the autosomal
**consensus** (like modern ancestry) instead of genotyping one alignment — an orientation-free join
that delivers GRCh37/38 alignments *and* chip-only subjects for free, reusing the cached, proven
frontend. (The earlier single-best-callable-alignment path was a CHM13-WGS-only stopgap.)

### 7.17 Consensus contract + a note on progressive consensus (design direction)

Deep ancestry now **requires** the autosomal consensus (errors "build the autosomal consensus first"
when absent), the same contract as modern/fine ancestry and painting — rather than building it
silently on demand. The heavy full-1240k genotyping happens once, in the shared Autosomal-tab build
(which streams progress); every downstream estimate (modern, fine, deep, painting) is then a fast read
over the cached `DiploidProfile`.

**Design direction (not yet built): make the consensus *progressive*.** Today the autosomal consensus
is a single on-demand build that genotypes every source at 1240k in one pass. It would be better for
each **batch import** to incrementally contribute its genotypes as it processes — Y-DNA, mtDNA, and the
autosomal **panel** (1240k) sites — so the consensus accumulates during ingest instead of a separate
heavy step afterwards. The building blocks already exist: `ibd_panel_dosages` resolves any source
(any build, or a chip) to canonical CHM13 panel dosages, and `reconcile_diploid` pools per-source
observations by site — it is already an incremental reducer (add a source's obs, re-vote). The work is
to (a) genotype+persist each imported source's 1240k panel dosages during batch processing (as the Y/mt
placements already are), and (b) have the consensus be an *update* over the accumulated per-source
dosages rather than a from-scratch genotyping pass. Then deep/modern ancestry need no build step at all
after import — the consensus is always current. Scope: the import/batch pipeline + a per-source panel
dosage store; deferred pending a decision to take it on.

### 7.18 Progressive consensus — first increment (built)

Made the autosomal consensus **accumulate** instead of being a lazy from-scratch build:

- **Reduce vs genotype are now separate.** `refresh_autosomal_consensus` reduces over the per-source
  dosages that are *already available* — every chip / WGS-VCF (which resolve with no decode) plus any
  alignment whose IBD-panel dosages are *cached* (`cached_alignment_panel_dosages`) — and **never**
  decodes an uncached alignment. It is cheap and safe to call after every import.
- **Wired into batch import** (`add_sample_dir`): each imported sample refreshes the consensus with
  whatever is available. Y-only imports (D2C) contribute nothing autosomal and pay ~nothing; a chip
  folds in immediately.
- **Panel batch-process mode** (`genotype_panel_for_subject`, CLI `genotype-panel`): genotype the
  subject's single **best-callable** alignment (`genome_territory × pct_10x × (1−pct_exc_mapq)`, any
  build — the IBD panel re-keys GRCh37/38) at the 1240k panel + cache, then reconcile **once**. Chips
  and VCFs fold in via the refresh, so this pays at most **one** whole-genome decode per subject —
  not one per redundant same-person alignment. `genotype_panel_for_alignment` is the primitive
  (genotype + cache only; the caller reconciles at the batch boundary — reconciling millions of
  observations per source is wasted work if repeated per alignment).

Validated: `genotype-panel huF98AFD` → best alignment #3713, 1.23M sites, consensus refreshed (2m14s,
one decode); `deep-ancestry` then reads it → WHG 14.6 / EEF 44.7 / Steppe 40.6, p=0.30.

**Still to do** (deferred): fold the panel batch mode into the project-wide analyze/deep-analyze
streaming flow with progress; a GUI trigger; and — the broader §7.17 vision — Y/mt already accumulate
at import via the haplogroup-GVCF sidecar, so the autosomal panel is the piece this closes.
