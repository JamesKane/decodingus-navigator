# Tier B archaic-segment validation

The harness that produced the numbers in
[`documents/design/ArchaicAncestry_Design.md`](../../documents/design/ArchaicAncestry_Design.md)
§ *Tier B validation*, and the gate for ever setting `ARCHAIC_SEGMENTS_ENABLED` back to `true`.

## Why this exists

Tier B was shipped on one number: its total archaic extent landed at 1.01× the hmmix European mean.
That number was produced by three parameters fitted until it did. **A caller that emits a calibrated
constant passes a cohort-mean test**, so the mean says nothing about whether the caller measures the
person in front of it.

These scripts ask the two questions that do:

1. **Locations** — do we call archaic sequence *where* an independent callset does, better than
   random placement would? (`compare_locations.py`)
2. **Amounts** — across individuals, does our extent rise and fall with theirs? (`correlate_extent.py`)

Tier B scored *below* the random null on (1) and r = −0.018 on (2), which is why it is gated off.

## What you need

- **hmmix's published 1000G callset** (Zenodo, CC BY 4.0) — the same record checkpoint A used.
  `hmmix_segments_chr21_22.tsv` covers 2,310 individuals across four regions on chr21+22;
  `hmmix_eur_all.tsv` is genome-wide but Europeans only. Both are hg38.
- **CHM13 alignments in the workspace for those same individuals.** The comparison is only worth
  anything against *the same people* — a distribution comparison is what got us here.
- `CrossMap` and `~/.decodingus/liftover/hg38ToHs1.over.chain` for the lift.

## Procedure

**1. Call the contigs the truth covers.** chr21+22 is the right axis: the truth covers every
individual there, and `pct_callable` restricts its denominator to the contigs actually called, so a
partial run is directly comparable rather than silently wrong.

```sh
ls runs || mkdir runs
for s in HG00096 HG00133 ...; do ./run_one.sh "$s" "$PWD" /path/to/navigator; done   # ~90 s each
```

Run 2 at a time (`xargs -P 2`): one sample saturates ~6 cores of a 16-core machine, so two fill it.

**2. Build the truth, lifted to CHM13.** Emit one BED line per hmmix segment **with a unique id**,
CrossMap it, then reassemble.

Two traps, both of which produced wrong answers before they were caught:

- **Reassemble from the lifted fragments, not from min..max of them.** CrossMap splits a segment
  into many pieces, and a couple land across a rearranged region — taking the span inflated
  HG00096's 2.3 Mb to 26.6 Mb. The fragment-length sum should equal the hg38 input exactly; check it.
- **Merge fragments with a ~1 kb tolerance.** The lift leaves median 2 bp gaps, which a strict union
  will not close, turning 48 real tracts into 423 shards and making per-tract recovery meaningless.

**3. Union across haplotypes — never sum.** hmmix reports per haplotype; our caller is unphased.
Summing doubles their figure and would make a caller look correctly calibrated when it is not.
Unioning reproduces their published EUR mean of 2.09 Mb on chr21+22, which is the check that the
truth-side arithmetic is right.

**4. Compare.**

```sh
python3 compare_locations.py runs/HG00096.json truth_HG00096.chm13.bed chr21 chr22
python3 correlate_extent.py          # reads runs/*.json, prints r, rho, permutation p, spread
```

## Reading the output

- **Sensitivity must beat the null**, which `compare_locations.py` does not compute for you — draw
  segments of your own lengths at random within the truth's span and measure the same overlap. For
  Tier B that null was 5.0 % (p95 9.4 %) against an observed 2.1 %.
- **`r` near zero with a wide truth spread means a calibrated constant**, not a measurement. Check
  the spread ratio too: an estimator of a varying quantity should vary about as much as the quantity
  does. Tier B's was 0.63×.
- **Before concluding the caller is wrong, rule out the harness.** Coordinate offset (cross-correlate
  overlap against ±shift; a real bug peaks off zero), callability (what fraction of the truth is even
  reachable — use `cargo run --example archaic_callable_dump`), and the two lift traps above.
