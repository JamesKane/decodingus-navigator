# shellcheck shell=bash
# Config for the Y-SNP name→locus dictionary build. Sourced by every stage.
# Override any value via the environment, e.g. `WORK=/data/ysnp ./02_build.sh`.
#
# URLs marked `# VERIFY` are best-known locations as of writing — confirm before a run
# (YBrowse refreshes its master CSV ~weekly; UCSC chain filenames can roll).

set -euo pipefail

# ── working tree ────────────────────────────────────────────────────────────────
WORK="${WORK:-$HOME/.decodingus/ysnp-build}"
RAW="$WORK/raw"            # untouched downloads (YBrowse CSVs, chains, FASTA)
TMP="$WORK/tmp"            # intermediates (BEDs, lifted output)
ASSETS="${ASSETS:-$HOME/.decodingus/ysnp}"  # where the app reads dictionary.tsv / aliases.tsv
LOG="$WORK/log"

# ── builds emitted into the dictionary ──────────────────────────────────────────
# Native YBrowse extracts give GRCh38 + GRCh37 directly; hs1 (CHM13v2) is added by lifting
# the GRCh38 coordinate. Build keys MUST match the labels the app resolves against
# (navigator-domain ysnp_dict / haplo): "GRCh38", "GRCh37", "hs1".
BUILD_GRCH38="GRCh38"
BUILD_GRCH37="GRCh37"
BUILD_HS1="hs1"

# ── YBrowse master Y-SNP catalog (Thomas Krahn) ─────────────────────────────────
# Quoted CSV, GFF-shaped. Columns (1-based):
#   1 seqid  2 source  3 type  4 start  5 end  6 score  7 strand  8 phase
#   9 Name  10 ID  11 allele_anc  12 allele_der  13 YCC_haplogroup  14 ISOGG_haplogroup
#   15 mutation  16 count_tested  17 count_derived  18 ref  19 comment
YBROWSE_HG38_URL="${YBROWSE_HG38_URL:-https://ybrowse.org/gbrowse2/gff/snps_hg38.csv}"  # VERIFY
YBROWSE_HG19_URL="${YBROWSE_HG19_URL:-https://ybrowse.org/gbrowse2/gff/snps_hg19.csv}"  # VERIFY

# ── liftover GRCh38 → CHM13v2 (hs1) ─────────────────────────────────────────────
# CrossMap needs the chain + the TARGET FASTA. Same sources the ancestry pipeline uses.
CHAIN_GRCH38_TO_CHM13="${CHAIN_GRCH38_TO_CHM13:-https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz}"  # VERIFY
CHM13_FASTA_URL="${CHM13_FASTA_URL:-https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz}"

# The Y contig name in YBrowse output is `chrY`; CHM13 analysis-set also uses `chrY`.
Y_CONTIG="${Y_CONTIG:-chrY}"
