# Why Tier B failed — the diagnostic scripts

These are the instruments behind
[`ArchaicAncestry_Design.md`](../../../documents/design/ArchaicAncestry_Design.md) § *Why it failed*.
They establish that the Tier B segment caller was built on the **wrong observable**, not mis-tuned —
and, just as importantly, they rule out the explanations that looked right first.

Kept because every one of the rejected hypotheses is plausible enough to be proposed again.

## Run order

All expect a scratch directory containing, for one individual (HG00096 in the recorded run):

| file | produced by |
|---|---|
| `HG00096.chr{21,22}.calls.json` | `sqlite3` on `analysis_artifact`, kind `diploid_denovo:chrN` |
| `private.chr{21,22}.tsv` | `cargo run --example archaic_private_dump` |
| `classify.tsv` | `cargo run --example archaic_classify_dump` |
| `callable.bed` | `cargo run --example archaic_callable_dump` (threshold 0.5) |
| `og_density.tsv` | `cargo run --example archaic_outgroup_density` |
| `truth_<sample>.chm13.bed` | hmmix segments lifted hg38→CHM13 — see [../README.md](../README.md) |

## What each one answers

| script | question | answer it gave |
|---|---|---|
| `background_variation.py` | is the background the flat Poisson the emission assumes? | **No** — 5.3× p10–p90 spread, **14.6× overdispersed**, larger than the 2.89× signal |
| `mutrate_check.py` | would a mutation-rate map fix that? | **Not enough** — the best available proxy explains 38% of the variance, leaving 7.4× |
| `quality_effect.py` | is the excess variance our caller's artifacts? | **No** — filtering lowers overdispersion only by discarding variants; enrichment falls with it |
| `callset_compare.py` | is our variant calling diluting the signal? | **No** — 1000G's own calls for the same person give the *same* contrast (1.98× vs 2.08×), though ours are 6× noisier |
| `truth_selfcheck.py` | is the truth set (or my lift) wrong? | **No** — hmmix's tracts are enriched 1.84× for their own archaic SNPs in native hg38, null 1.04× |
| `observable_compare.py` | is a more specific observable available? | private ∩ diagnostic is unusable (7 sites); the African strip removes almost all diagnostic sites |
| `haplotype_match.py` | does archaic-allele matching separate tracts? | **Yes** — 39.5% carrying inside tracts vs 13.0% elsewhere, over ~30 sites per tract |
| `discriminability.py` | how detectable is ONE tract under each observable? | density **14.3%** vs matching **95.1%** sensitivity at 5% false positives |

`discriminability.py` is the one to read first. It needs no inputs and is the whole argument:
both observables carry ~3× contrast, so contrast was never the problem — **evidence per tract** was.

## The trap these scripts exist to prevent

Two harness bugs in here produced confident wrong answers before being caught, both worth knowing:

- **`haplotype_match.py` originally conditioned on the subject having a variant call** at a
  diagnostic site, which samples only sites where he already has a variant and reported an
  impossible ~80% carrying rate against a known 4.3% background. Every diagnostic site in callable
  territory belongs in the denominator; no call means hom-reference.
- **Binning at 100 kb to measure in-tract contrast** dilutes a 36 kb tract with 64 kb of background
  and reported 1.14× where the correct measurement is 1.98×. Measure contrast at tract boundaries,
  not in fixed bins.

More generally: sensitivity alone is gameable by calling more sequence, because the
random-placement null rises with it. Always report the null at the extent actually called.
