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
| 1 | `01_fetch.sh` | download CHM13 FASTA + chains, AADR (Dataverse), 1000G-CHM13 **AF** files; optionally SGDP. Big genotype sets are NOT bulk-pulled (see Download footprint) |
| 2 | `02_liftover_panel_sites.sh` | lift the AADR 1240k site universe hg19 → CHM13 (and → hg38 if GRCh38 sources are enabled) |
| 3 | `03_select_panel.sh` | restrict 1000G-CHM13 AF to 1240k∩CHM13, Fst-select AIMs → `ancestry_panel_<build>.bin` |
| 4 | `04_build_matrices.sh` | per source: get genotypes (slice 1000G BCF / slice gnomAD remotely / convert AADR+SGDP), liftover, align to CHM13 ref, cut to panel sites → matrices + pop map |
| 5 | `05_build_assets.sh` | `panelbuild pca` + `fine-panel` → global PCA + freq assets + provenance manifest |
| 6 | `06_publish_cdn.sh` | sha256 + upload assets + manifest to the CDN (`--apply` to actually upload) |

Run in order: `for s in 01 02 03 04 05; do ./"$s"_*.sh; done` then `./06_publish_cdn.sh --apply`.
Run **`./derive_modern_pops.sh`** once before stage 04 (it needs no stage-1 output — it pulls its
own public metadata) so the modern sources get labelled. Every stage is idempotent; intermediates
live under `$WORK` (default `~/.decodingus/ancestry-build`).

## Download footprint

The panel needs only ~20k AIM sites, so the pipeline **slices** rather than mirroring whole
genomes. The naive "download everything" pull would be ~5 TB; the default build pulls **~20 GB**:

| Source | What's pulled | Size |
|--------|---------------|------|
| 1000G-CHM13 AF (`unrelated_samples_2504/allele_freq`) | per-chrom `withafinfo` VCFs (carry `AC_<POP>_unrel` for the panel) | ~9.9 GB |
| 1000G-CHM13 genotypes (phased biallelic 3202 BCF) | **downloaded whole (~13 GB)** then sliced locally — remote `-R` over S3 is unreliable ("Illegal seek") | ~13 GB |
| AADR v66.p1 1240K (Harvard Dataverse) | geno/snp/ind/anno | ~7.3 GB |
| CHM13 FASTA + chains | reference + liftover | ~1 GB |
| **gnomAD HGDP+1KG** (`HGDP_1KG_ENABLE=1`) | **IMPRACTICAL** — chr22 alone ~60 GB (~2 TB total, no bulk DL); remote `-R` ~5 s/site (~28 h, fragile). Off; prefer SGDP/AADR | — |
| **SGDP** (`SGDP_ENABLE=1`) | GRCh37 PLINK fetched whole, panel-restricted + lifted hg19→CHM13 | ~3 GB |

The phased biallelic BCF (`KGP_GT_BCF_URL`) is one ~13 GB whole-genome file: point it at a local
mirror, or the pipeline downloads it once (htslib remote `-R` slicing of thousands of scattered
sites over the network is too slow/fragile). gnomAD/SGDP are **off by
default** (modern global resolution is an enhancement; 1000G-CHM13 + AADR give a working build).

## Prerequisites

External tools (not bundled): `curl`, `bcftools` + `tabix`, `samtools`, `CrossMap`
(`pip install CrossMap`), `plink2`, **AdmixTools' `convertf`**, **GATK** (`gatk` ≥4), `awk`, and
`aws` **or** `rclone` (publish). Plus this repo's `navigator-panelbuild`
(`cargo run -p navigator-panelbuild`). Notes:
- **convertf** must be AdmixTools', not EIGENSOFT's — recent AADR `.geno` ships as TGENO
  (transposed packed), which AdmixTools convertf reads but EIGENSOFT's does not.
- **VCF liftover uses GATK `LiftoverVcf`, not CrossMap** — it is allele-aware (reverse-complements
  on inverted/minus-strand chain blocks and swaps REF↔ALT + flips genotypes via
  `RECOVER_SWAPPED_REF_ALT`). CrossMap's `vcf` mode silently blanks the ALT wherever the CHM13 base
  differs from the source REF (most AIMs), corrupting ~3/4 of the lifted genotypes. CrossMap is
  still used for the allele-free **BED** site-lift in stage 02. Tune heap/spill with
  `GATK_JAVA_OPTS` / `GATK_MAX_RECORDS_IN_RAM` (defaults `-Xmx16g` / `2000`, sized for wide
  many-sample VCFs).
- bioconda has `admixtools` and `eigensoft`; install via conda/micromamba.

## You must provide / curate

- **`pops/aadr_component_map.tsv`** — AADR Group ID → deep-component map (the one
  expertise-driven artifact; a starter set is included — review it).
- **`$RAW/<src>.pops.tsv`** — `sample<TAB>population` for each modern source (1000G/HGDP/SGDP).
  Generate them from public metadata with **`./derive_modern_pops.sh`** (idempotent; `--force` to
  rebuild). It writes `1kgp.pops.tsv` (required — the 2504 unrelated → 26 fine pops), plus
  `sgdp.pops.tsv` and `hgdp1kg.pops.tsv` (only consumed when those optional sources are enabled).
  Columns are resolved by header name, so upstream reordering won't silently mis-map.
- **`# VERIFY` URLs in `config.sh`** — dataset versions/filenames roll forward. AADR is on
  Harvard Dataverse (release **v66.p1** = Dataverse version 14, fetched by numeric file id — re-pin
  `AADR_ID_*` from the Dataverse API when bumping `AADR_VERSION`). SGDP moved to
  `sharehost.hms.harvard.edu/.../variant_set` (its `.bim` ships as `.bim.zip`). (URLs web-verified
  2026-06-13.)

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
- **Projection-mode PCA** (basis = modern, project ancient) is wired: stage 5 passes
  `--basis-pops` so ancient deep components are projected onto the modern basis, not baked
  into the axes. Still keep *extremely* low-coverage ancient samples out of the pop map.
- **App asset wiring**: the app currently loads `ancestry_pca_<build>.bin` (+ optional
  `ancestry_pca_ancient_<build>.bin`). The combined global asset *is* `ancestry_pca_<build>.bin`,
  so all three methods use it; consuming `ancestry_freq_global` for fine admixture is a
  follow-up in `navigator-app`.
