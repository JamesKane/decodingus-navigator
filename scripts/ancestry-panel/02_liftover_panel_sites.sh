#!/usr/bin/env bash
# Stage 2 — lift the AADR 1240k site universe (hg19) onto CHM13v2.
#
# Ancient samples only carry data at the 1240k capture sites, so 1240k is the universe
# the panel is selected from. The AADR `.snp` file lists those sites in hg19; we lift
# them to CHM13 as a BED, dropping sites that don't map. Output: $TMP/1240k_sites.<build>.bed
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool CrossMap "pip install CrossMap"
require_tool awk

SNP="$(ls "$RAW/${AADR_DATASET}"*.snp 2>/dev/null | head -1 || true)"
[[ -n "$SNP" ]] || die "AADR .snp not found in $RAW (run 01_fetch.sh / unpack AADR). Looked for ${AADR_DATASET}*.snp"

CHAIN="$(chain_for hg19)"
OUT="$TMP/1240k_sites.${BUILD}.bed"
mkdir -p "$TMP"

# EIGENSTRAT .snp columns: id  chrom  genpos  physpos  ref  alt  (whitespace-separated).
# Emit a 0-based BED (chrom, pos-1, pos, id\tref\talt) in hg19, then CrossMap to CHM13.
log "building hg19 BED from $(basename "$SNP")"
awk 'NF>=6 { chr=$2; if (chr=="23") chr="X"; if (chr=="24") chr="Y";
            printf "chr%s\t%d\t%d\t%s|%s|%s\n", chr, $4-1, $4, $1, $5, $6 }' "$SNP" \
  > "$TMP/1240k_sites.hg19.bed"

log "CrossMap bed -> $BUILD"
CrossMap bed "$CHAIN" "$TMP/1240k_sites.hg19.bed" "$OUT" || die "CrossMap bed failed"

# panelbuild/bcftools want a tab-separated CHROM<TAB>POS regions file too (1-based).
awk '{ split($4,a,"|"); printf "%s\t%d\n", $1, $3 }' "$OUT" | sort -k1,1 -k2,2n -u \
  > "$TMP/1240k_sites.${BUILD}.tsv"

log "stage 2 complete: $(wc -l < "$OUT") sites lifted -> $OUT (+ .tsv regions)."
