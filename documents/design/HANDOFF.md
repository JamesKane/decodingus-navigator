# Rust rewrite — handoff / resume notes

Last updated: 2026-06-04. Branch: `rust-rewrite`. Pick up here next session.

## State

Whole workspace builds; `cargo test --workspace` and `cargo clippy --all-targets -- -D warnings`
are green. Recent work (newest first):

- `d3d4b6e` Wire SV calling into app (single + bulk + report)
- `cf85edd` Wire sex inference + read metrics (single + bulk + report)
- `8e21c93` Ancestry Phase C — local-ancestry "DNA painting" (AF-based HMM)
- `1ce7c95` Ancestry Phase D — donut + geographic map visuals
- `5136efb` Ancestry Phase B — SGDP panel (Middle East, Central Asia/Siberia, Oceania)
- `3595e05` Ancestry — supervised admixture composition + hierarchy
- `4ed5d43` Ancestry — fine-grained 26-population resolution
- `33c0353` Ancestry Phase 2 — PCA; `45e6506` Ancestry Phase 1 — AF-likelihood + panel tool

Nothing uncommitted of note. Untracked: `CLAUDE.md`, `GEMINI.md`, `.claude/` (agent config — leave).

## Build / validate

```bash
cargo build && cargo test --workspace && cargo clippy --all-targets -- -D warnings
cargo run -p navigator-ui            # desktop app
```

### Live (`#[ignore]`) tests — real data, run explicitly
Test sample: `/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam` (CHM13 HiFi, ~4×,
male, Y=R-FGC29071, mtDNA=U5a1b1g, European). Reference: `/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa`.

```bash
GFX_CHM13_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam \
GFX_CHM13_REF=/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa \
NAVIGATOR_ANCESTRY_PANEL=/Users/jkane/.decodingus/ancestry/ancestry_panel_chm13v2.0.bin \
NAVIGATOR_ANCESTRY_PCA=/Users/jkane/.decodingus/ancestry/ancestry_pca_chm13v2.0.bin \
  cargo test -p navigator-app --release \
  validate_gfx_chm13_ancestry local_ancestry_paints_gfx gfx_sex_is_male -- --ignored --nocapture
```
Expected: European ~98% (admixture), DNA painting EUR-dominant, sex=Male. Other ignored live
tests: `validate_gfx_chm13_haplogroups` (Y/mt), parity_real.rs (HG002 env), PDS publish (PDS_TEST_URL).

## Ancestry assets (regenerable; not committed)

Installed at `~/.decodingus/ancestry/`:
- `ancestry_panel_chm13v2.0.bin` — fine (29-pop) + SGDP AF panel (the default the app loads)
- `ancestry_pca_chm13v2.0.bin` — PCA loadings (29 pops)
- `*_super_*` / `*_fine_*` — earlier super-pop-only / 1000G-only backups

Rebuild from the archived genotype matrix `~/Genomics/archive/1kgp_chm13_pca_build/`
(`gt_all.tsv.gz` 1000G + `sgdp_gt.tsv.gz` + `combined_pops.txt` + README with exact commands):
```bash
A=~/Genomics/archive/1kgp_chm13_pca_build; O=~/.decodingus/ancestry
navigator-panelbuild fine-panel --matrix $A/gt_all.tsv.gz,$A/sgdp_gt.tsv.gz \
  --samples $A/samples.txt,$A/sgdp_subset_samples.txt --pops $A/combined_pops.txt --out $O/ancestry_panel_chm13v2.0.bin
navigator-panelbuild pca        --matrix $A/gt_all.tsv.gz,$A/sgdp_gt.tsv.gz \
  --samples $A/samples.txt,$A/sgdp_subset_samples.txt --pops $A/combined_pops.txt --out $O/ancestry_pca_chm13v2.0.bin
```
The super-pop AIMs panel itself was built from the sites-only 1KGP-CHM13 VCFs in
`~/Genomics/1kgp_chm13_af/` (9 GB, re-downloadable from the human-pangenomics S3 bucket).

## EC2 (genotype extraction)

`admin@ec2-3-132-31-28.us-east-2.compute.amazonaws.com`, key `~/Decoding-Us.pem` (chmod 600),
Debian, bcftools installed. Used to tabix-fetch panel-site genotypes from the ~1 TB 1000G/SGDP
genotype VCFs (in-AWS S3 is fast; a laptop fetch is .tbi/latency-bound). The recipe + region
files are archived in `1kgp_chm13_pca_build/ancestry_build.tar.gz`. The matrices are already
pulled, so re-extraction is only needed to add reference samples (e.g. a fine-Fst panel).

## Key gotchas (also in agent memory)

- CHM13 `chrM` is a circular permutation of rCRS → rotation-aware self-generated map.
- CHM13 Y has inverted tracts → liftover reverse-complements minus-strand lifts.
- `allele_freq/*withafinfo*` VCFs are **sites-only** (no genotypes despite the header); real
  genotypes are in `unrelated_samples_2504/` + `SGDP/` (anonymized IDs → `SGDP_sample_info.txt`).
- PCA projection of a low-coverage sample is rescaled by `total/used` sites (else it shrinks toward
  origin and mis-clusters).
- SV needs ≥10× — won't run on the 4× GFX sample (gate returns "coverage too low").
- AIMs were selected for 1000G **super-pop** Fst → within-continent fine resolution is noisy.

## Recommended next steps (pick one)

1. **IBD Matching + chromosome browser** — biggest remaining user-facing feature; detection +
   identity math already built, the consent/match-discovery/browser UI is not.
2. **Parity-harness automation** — make `parity_real.rs` a real CI gate (currently all `#[ignore]`),
   the stated cutover criterion (§4c).
3. Ancestry refinements — fine-Fst panel (sharpen sub-continental + SGDP continents), real CHM13
   genetic map for the painting HMM, pair-state (per-haplotype) painting, country-polygon map.
4. Validate SV on a ≥10× sample (e.g. HG002 in `~/Genomics/HG002`).
