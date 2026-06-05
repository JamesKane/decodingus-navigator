# Global ancestry-panel build pipeline (Option B)

Builds the global ancestry reference assets — a 1240k-restricted AIMs panel plus a global
PCA basis (modern + ancient deep components) — that feed all three estimators
(`estimate_admixture`, `estimate_pca_gmm`, `estimate_nmonte`). Output is published to a CDN
for Navigator clients to download.

**Why scripts, not raw alignments:** reference building needs *genotypes / allele
frequencies / PCA coordinates*, never raw reads. So this is decoupled from the multi-year
CHM13 alignment effort — it runs entirely off published call sets.

## Pipeline

| Stage | Script | Does |
|-------|--------|------|
| 0 | `config.sh` | all paths, versions, URLs, panel params, CDN target (override via env) |
| 1 | `01_fetch.sh` | download CHM13 FASTA + liftover chains, AADR, HGDP+1KG, SGDP, 1000G-CHM13 |
| 2 | `02_liftover_panel_sites.sh` | lift the AADR 1240k site universe hg19 → CHM13 |
| 3 | `03_select_panel.sh` | restrict 1000G-CHM13 to 1240k∩CHM13, Fst-select AIMs → `ancestry_panel_<build>.bin` |
| 4 | `04_build_matrices.sh` | per source: convert→VCF, liftover, align to CHM13 ref, cut to panel sites → matrices + pop map |
| 5 | `05_build_assets.sh` | `panelbuild pca` + `fine-panel` → global PCA + freq assets + provenance manifest |
| 6 | `06_publish_cdn.sh` | sha256 + upload assets + manifest to the CDN (`--apply` to actually upload) |

Run in order: `for s in 01 02 03 04 05; do ./"$s"_*.sh; done` then `./06_publish_cdn.sh --apply`.
Every stage is idempotent; intermediates live under `$WORK` (default `~/.decodingus/ancestry-build`).

## Prerequisites

External tools (not bundled): `curl`, `bcftools` + `tabix`, `CrossMap` (`pip install CrossMap`),
`convertf` + `plink2` (AADR EIGENSTRAT → VCF), `awk`, and `aws` **or** `rclone` (publish).
Plus this repo's `navigator-panelbuild` (run via `cargo run -p navigator-panelbuild`).

## You must provide / curate

- **`pops/aadr_component_map.tsv`** — AADR Group ID → deep-component map (the one
  expertise-driven artifact; a starter set is included — review it).
- **`$RAW/<src>.pops.tsv`** — `sample<TAB>population` for each modern source (1000G/HGDP/SGDP).
- **`# VERIFY` URLs in `config.sh`** — dataset versions/filenames roll forward; confirm the
  current AADR release, gnomAD HGDP+1KG subset names, SGDP location, and the CHM13 chains
  before a real run.

## Outputs (in `$ASSETS`, default `~/.decodingus/ancestry`)

- `ancestry_panel_<build>.bin` — AIMs AF panel (genotyping + super-pop admixture)
- `ancestry_pca_<build>.bin` — global PCA loadings + per-population centroids (GMM + nMonte)
- `ancestry_freq_global_<build>.bin` — global per-population AF (fine admixture)
- `ancestry_manifest_<build>.json` — provenance + sha256 (published; clients verify against it)

## Known refinements (tracked in `documents/design/AncestryAnalysis.md`)

- **Allele harmonization** to the CHM13 reference (`bcftools norm -c s`) is the highest-risk
  step — a missed strand flip silently corrupts dosages. Spot-check a known sample after a build.
- **Pseudo-haploid ancient genotypes** (AADR 1240k) inflate centroid variance — fine for a
  first asset, worth modelling later.
- **Projection-mode PCA** (basis = modern, project ancient) is a `navigator-panelbuild`
  enhancement; until it lands, keep very-low-coverage ancient samples out of the pop map.
- **App asset wiring**: the app currently loads `ancestry_pca_<build>.bin` (+ optional
  `ancestry_pca_ancient_<build>.bin`). The combined global asset *is* `ancestry_pca_<build>.bin`,
  so all three methods use it; consuming `ancestry_freq_global` for fine admixture is a
  follow-up in `navigator-app`.
