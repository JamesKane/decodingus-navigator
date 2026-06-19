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
# Virtual-hosted-style S3 (portable from any region; path-style 301-redirects outside us-east-1).
CHM13_FASTA_URL="${CHM13_FASTA_URL:-https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz}"
# Liftover chains TO CHM13v2 (hs1). # VERIFY exact filenames at the base below.
#   T2T chains: https://s3.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/
#   UCSC hub:   https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/
CHAIN_HG19_TO_CHM13="${CHAIN_HG19_TO_CHM13:-https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg19-chm13v2.over.chain.gz}"   # VERIFY
CHAIN_GRCH38_TO_CHM13="${CHAIN_GRCH38_TO_CHM13:-https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz}" # VERIFY
# hg19 -> hg38 (standard UCSC chain) — used only to project the 1240k AIM universe into hg38 so
# GRCh38 sources (gnomAD/SGDP) can be sliced to just those sites before download.
CHAIN_HG19_TO_HG38="${CHAIN_HG19_TO_HG38:-https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz}"

# ── source datasets ─────────────────────────────────────────────────────────────
# 1000 Genomes recalibrated on CHM13v2 (native build — no liftover). Two distinct roles:
#   (a) PANEL (stage 03): per-super-pop AF in INFO (AC_<POP>_unrel/AN_<POP>_unrel). Those live
#       in the sites-only `withafinfo` VCFs under unrelated_samples_2504/allele_freq/ — ~9.9 GB
#       total, NOT the 1.6 TB per-genotype all_samples_3202 VCFs. panelbuild reads INFO only.
#   (b) GENOTYPES (stage 04 PCA basis): the phased biallelic 3202 BCF (~13 GB, one whole-genome
#       file). Stage 04 remote-slices it at the panel sites, so the actual pull is a fraction.
KGP_S3_BASE="${KGP_S3_BASE:-https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0}"
# (a) AF files — downloaded whole into KGP_CHM13_DIR; stage 03 restricts them to the 1240k sites.
KGP_CHM13_DIR="${KGP_CHM13_DIR:-$RAW/1kgp-chm13}"
KGP_AF_BASE_URL="${KGP_AF_BASE_URL:-$KGP_S3_BASE/unrelated_samples_2504/allele_freq}"
KGP_AF_PATTERN="${KGP_AF_PATTERN:-1KGP.CHM13v2.0.chr%s.recalibrated.snp_indel.pass.withafinfo.vcf.gz}"
# (b) Whole-genome phased biallelic genotype BCF — remote-sliced at panel sites (set to a local
#     path to use a mirror instead of streaming). Its .csi index sits next to it on S3.
KGP_GT_BCF_URL="${KGP_GT_BCF_URL:-$KGP_S3_BASE/Phased_SHAPEIT5_v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf.gz}"

# Allen Ancient DNA Resource (AADR) — ancient deep-component sources. EIGENSTRAT.
# Distributed via Harvard Dataverse (the old reichdata.hms.harvard.edu/.../curated_releases
# flat path is gone / access-restricted). Each file is fetched by its numeric Dataverse file
# id, pinned per release below. To move to a newer release, refresh AADR_VERSION + the id map:
#   curl 'https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId='"$AADR_DATAVERSE_DOI"
#   landing page: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
AADR_VERSION="${AADR_VERSION:-v66.p1}"  # current Dataverse release (v66.p1 = Dataverse dataset version 14; older notes said v62/v66 — now superseded)
AADR_DATASET="${AADR_DATASET:-1240K}"   # 1240K (~1.23M ancient-capture set, ~7.3 GB) | HO (~600k, ~4 GB)
AADR_DATAVERSE_DOI="${AADR_DATAVERSE_DOI:-doi:10.7910/DVN/FFIDCW}"
AADR_DOWNLOAD_BASE="${AADR_DOWNLOAD_BASE:-https://dataverse.harvard.edu/api/access/datafile}"
# Dataverse numeric file ids for the 1240K quartet (geno/snp/ind/anno). Re-pin per release —
# stale ids 403 via the access API. Verified against dataset version 14 (v66.p1) on 2026-06-13.
AADR_ID_GENO="${AADR_ID_GENO:-13994829}"
AADR_ID_SNP="${AADR_ID_SNP:-13994514}"
AADR_ID_IND="${AADR_ID_IND:-13994513}"
AADR_ID_ANNO="${AADR_ID_ANNO:-13994515}"
# Local filename stem the downloaded quartet shares (matches AADR's own naming).
AADR_FILE_PREFIX="${AADR_FILE_PREFIX:-${AADR_VERSION}.${AADR_DATASET}.aadr.PUB}"

# HGDP + 1KG dense callset (modern global, non-European resolution). gnomAD v3.1.2, GRCh38.
# OPTIONAL — but IMPRACTICAL on a workstation, off by default (see below). Anonymous gs:// reads of
# the public bucket DO work (no GCP project/auth needed; HGDP_1KG_GCP_PROJECT only if it ever flips
# to requester-pays). BUT: the per-chromosome callsets are dense whole-genome (chr22 alone is
# ~60 GB → ~2 TB total, so bulk download is out), and remote `-R` slicing is latency-bound at
# ~5 s/site over gs:// — ~28 h for the 20k-site panel and fragile (measured 2026-06-14). It also
# re-includes the 1KG samples (overlap with our 1000G source). RECOMMENDATION: leave disabled;
# 1000G + AADR + SGDP already give modern + ancient global coverage. If you truly need HGDP, source
# it from a smaller distribution (e.g. AADR present-day HGDP samples, already local) rather than the
# gnomAD gs:// callset. The pop map (hgdp1kg.pops.tsv) is the gnomAD meta `hgdp_tgp_meta.Population`.
HGDP_1KG_ENABLE="${HGDP_1KG_ENABLE:-0}"   # 1 to include (impractical via gnomAD gs://; see above)
HGDP_1KG_GCP_PROJECT="${HGDP_1KG_GCP_PROJECT:-}"
HGDP_1KG_BASE_URL="${HGDP_1KG_BASE_URL:-gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes}"
HGDP_1KG_PATTERN="${HGDP_1KG_PATTERN:-gnomad.genomes.v3.1.2.hgdp_tgp.chr%s.vcf.bgz}"

# Simons Genome Diversity Project (modern deep diversity). cteam_extended is **GRCh37/hg19**
# (verified by 1240k position overlap — NOT GRCh38), distributed as PLINK by the Reich lab.
# OPTIONAL and tractable (works). Served from sharehost.hms.harvard.edu (valid TLS — the old
# reichdata.hms.harvard.edu path 404s; the -k hack no longer applies). NOTE: at this host the .bim
# ships zipped as `${SGDP_PLINK_PREFIX}.bim.zip` (stage 01 unzips it). The cteam .fam IID is the
# metadata's `Sample_ID(Aliases)` (HGDP*/SGDP ids), NOT `SGDP_ID` — sgdp.pops.tsv is keyed on the
# alias accordingly. Stage 04 restricts to the panel SNPs (hg19 `chrom_pos` .bim ids) and lifts
# hg19->CHM13 before genotyping.
SGDP_ENABLE="${SGDP_ENABLE:-0}"           # 1 to include
SGDP_BASE_URL="${SGDP_BASE_URL:-https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set}"

# ── IBD asset inputs (genetic map + GRCh38 FASTA) ────────────────────────────────
# Recombination map (deCODE-derived, GRCh38 PLINK format) for the IBD genetic-map asset. Stage 01
# downloads + unzips it; stage 05 lifts GRCh38->CHM13 (coordinate-only) and serializes.
GMAP_URL="${GMAP_URL:-https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip}"
# GRCh38 FASTA (UCSC hg38 — matches the hg19->hg38 chain) for the IBD panel's GRCh38 allele column,
# so GRCh38 consumer chips resolve. OPTIONAL: set IBD_GRCH38=0 to skip the ~1 GB download; the IBD
# panel then carries CHM13 + GRCh37 only (covers GRCh37 chips like 23andMe v5 / AncestryDNA + WGS).
IBD_GRCH38="${IBD_GRCH38:-1}"
GRCH38_FASTA_URL="${GRCH38_FASTA_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz}"
SGDP_PLINK_PREFIX="${SGDP_PLINK_PREFIX:-cteam_extended.v4.maf0.1perc}" # verified at sharehost 2026-06-13 (.bim is .bim.zip)

# Curated AADR population-label → deep-component map (edit this; ships in the repo).
AADR_COMPONENT_MAP="${AADR_COMPONENT_MAP:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pops/aadr_component_map.tsv}"

# ── CDN publish target ──────────────────────────────────────────────────────────
# Where 06_publish_cdn.sh uploads the released assets + manifest. s3:// (aws cli) or
# any rclone remote. The CDN (e.g. CloudFront) fronts this bucket/prefix.
CDN_REMOTE="${CDN_REMOTE:-s3://decodingus-assets}"
CDN_PREFIX="${CDN_PREFIX:-ancestry/$BUILD}"
ASSET_VERSION="${ASSET_VERSION:-1}"   # bump per published asset revision

# ── outputs (asset filenames the app/CDN expect) ────────────────────────────────
PANEL_OUT="$ASSETS/ancestry_panel_${BUILD}.bin"            # AF panel (genotyping + admixture)
PCA_OUT="$ASSETS/ancestry_pca_${BUILD}.bin"                # modern PCA loadings + centroids (scatter)
PCA_ANCIENT_OUT="$ASSETS/ancestry_pca_ancient_${BUILD}.bin" # PCA w/ ancient deep components (GMM/nMonte)
FINE_OUT="$ASSETS/ancestry_freq_global_${BUILD}.bin"       # global per-pop AF (fine admixture)
GMAP_OUT="$ASSETS/genetic_map_${BUILD}.bin"               # IBD recombination map (bp->cM)
IBD_PANEL_OUT="$ASSETS/ibd_panel_${BUILD}.bin"            # chip-compatible multi-build IBD SNP panel
MANIFEST="$ASSETS/ancestry_manifest_${BUILD}.json"         # provenance + checksums
