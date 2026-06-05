# shellcheck shell=bash
# Central config for the global ancestry-panel build pipeline (Option B).
# Sourced by every stage script. Override any value via the environment, e.g.
#   AADR_VERSION=v62.0 WORK=/data/ancestry ./01_fetch.sh
#
# URLs marked `# VERIFY` are best-known landing/base locations as of writing — confirm
# the exact current filenames before a production run (dataset versions roll forward).

set -euo pipefail

# ── target build ────────────────────────────────────────────────────────────────
# Coordinate base for every asset. The app genotypes the sample on its own build and
# the panel must match; CHM13v2 is the project default.
BUILD="${BUILD:-chm13v2.0}"

# ── working tree ────────────────────────────────────────────────────────────────
# All downloads + intermediates live here; final assets land in $ASSETS.
WORK="${WORK:-$HOME/.decodingus/ancestry-build}"
RAW="$WORK/raw"            # untouched downloads
TMP="$WORK/tmp"            # intermediates (matrices, lifted VCFs, …)
ASSETS="${ASSETS:-$HOME/.decodingus/ancestry}"   # where the app + CDN read assets from
LOG="$WORK/log"

# ── panel parameters ────────────────────────────────────────────────────────────
MAX_SITES="${MAX_SITES:-20000}"   # AIMs kept (highest Fst first)
MIN_FST="${MIN_FST:-0.10}"        # Nei Fst floor across super-pops
PCA_COMPONENTS="${PCA_COMPONENTS:-25}"  # 25 → G25-comparable coordinate space
MIN_CALL_RATE="${MIN_CALL_RATE:-0.5}"   # ancient data is sparse; keep the floor low

# ── reference CHM13v2 FASTA + liftover chains ───────────────────────────────────
# CHM13v2.0 analysis-set FASTA (CrossMap needs the *target* reference).
CHM13_FASTA_URL="${CHM13_FASTA_URL:-https://s3.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz}"
# Liftover chains TO CHM13v2 (hs1). # VERIFY exact filenames at the base below.
#   T2T chains: https://s3.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/
#   UCSC hub:   https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/
CHAIN_HG19_TO_CHM13="${CHAIN_HG19_TO_CHM13:-https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg19-chm13v2.over.chain.gz}"   # VERIFY
CHAIN_GRCH38_TO_CHM13="${CHAIN_GRCH38_TO_CHM13:-https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz}" # VERIFY

# ── source datasets ─────────────────────────────────────────────────────────────
# 1000 Genomes recalibrated on CHM13v2 (native build; the project already mirrors these).
# Point at a LOCAL mirror if you have one; else the public bucket.
# VERIFY: https://s3.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/
KGP_CHM13_DIR="${KGP_CHM13_DIR:-$RAW/1kgp-chm13}"
KGP_CHM13_BASE_URL="${KGP_CHM13_BASE_URL:-https://s3.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0}" # VERIFY

# Allen Ancient DNA Resource (AADR) — ancient deep-component sources. EIGENSTRAT/packed.
# Landing page (pick the current version + the 1240k set):
#   https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
AADR_VERSION="${AADR_VERSION:-v62.0}"   # VERIFY current release
AADR_BASE_URL="${AADR_BASE_URL:-https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases}" # VERIFY
AADR_DATASET="${AADR_DATASET:-1240k}"   # 1240k (~1.23M, the ancient-capture set) | HO (~600k)

# HGDP + 1KG dense callset (modern global, non-European resolution). gnomAD v3, GRCh38.
# Landing: https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg
HGDP_1KG_BASE_URL="${HGDP_1KG_BASE_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes}" # VERIFY subset filenames

# Simons Genome Diversity Project (modern deep diversity). GRCh38.
# Reich lab mirror: https://reichdata.hms.harvard.edu/pub/datasets/sgdp/
SGDP_BASE_URL="${SGDP_BASE_URL:-https://reichdata.hms.harvard.edu/pub/datasets/sgdp}" # VERIFY

# Curated AADR population-label → deep-component map (edit this; ships in the repo).
AADR_COMPONENT_MAP="${AADR_COMPONENT_MAP:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pops/aadr_component_map.tsv}"

# ── CDN publish target ──────────────────────────────────────────────────────────
# Where 06_publish_cdn.sh uploads the released assets + manifest. s3:// (aws cli) or
# any rclone remote. The CDN (e.g. CloudFront) fronts this bucket/prefix.
CDN_REMOTE="${CDN_REMOTE:-s3://decodingus-assets}"
CDN_PREFIX="${CDN_PREFIX:-ancestry/$BUILD}"
ASSET_VERSION="${ASSET_VERSION:-1}"   # bump per published asset revision

# ── outputs (asset filenames the app/CDN expect) ────────────────────────────────
PANEL_OUT="$ASSETS/ancestry_panel_${BUILD}.bin"        # AF panel (genotyping + admixture)
PCA_OUT="$ASSETS/ancestry_pca_${BUILD}.bin"            # global PCA loadings + centroids
FINE_OUT="$ASSETS/ancestry_freq_global_${BUILD}.bin"   # global per-pop AF (fine admixture)
MANIFEST="$ASSETS/ancestry_manifest_${BUILD}.json"     # provenance + checksums
