#!/usr/bin/env bash
# Stage 2 — build the dictionary asset from the YBrowse CSVs:
#   • GRCh38 + GRCh37 coordinates straight from the native extracts.
#   • hs1 (CHM13v2) coordinates by lifting the GRCh38 positions with CrossMap; alleles are
#     complemented when the lift inverts strand, so stored alleles are always the + strand
#     of the emitted build (matching how placement reads the sample's reference + strand base).
# Output: $ASSETS/dictionary.tsv (+ a header-only aliases.tsv — YBrowse lists synonymous SNP
# names as their own rows at the same locus, so each name is already a direct entry).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool awk
require_tool sort
require_tool CrossMap "pip install CrossMap"
ensure_dirs

HG38_CSV="$RAW/snps_hg38.csv"
HG19_CSV="$RAW/snps_hg19.csv"
[[ -s "$HG38_CSV" ]] || die "missing $HG38_CSV (run 01_fetch.sh)"

ROWS="$TMP/rows.tsv"; : > "$ROWS"

# ── native GRCh38 ────────────────────────────────────────────────────────────────
log "parsing GRCh38 rows"
ybrowse_to_rows "$HG38_CSV" "$BUILD_GRCH38" >> "$ROWS"

# ── native GRCh37 (optional) ──────────────────────────────────────────────────────
if [[ -s "$HG19_CSV" ]]; then
  log "parsing GRCh37 rows"
  ybrowse_to_rows "$HG19_CSV" "$BUILD_GRCH37" >> "$ROWS"
else
  log "no hg19 CSV — skipping GRCh37 coords"
fi

# ── lift GRCh38 → hs1 ──────────────────────────────────────────────────────────────
# CrossMap reads gzipped chains directly; use whichever 01_fetch left in place.
CHAIN="$RAW/$(basename "$CHAIN_GRCH38_TO_CHM13")"
[[ -s "$CHAIN" ]] || CHAIN="$RAW/$(basename "$CHAIN_GRCH38_TO_CHM13" .gz)"
[[ -s "$CHAIN" ]] || die "missing GRCh38→CHM13 chain in $RAW (run 01_fetch.sh)"

# BED6 from the GRCh38 rows: name field encodes NAME|ANC|DER; input strand is '+'.
BED_IN="$TMP/grch38_sites.bed"
awk -F'\t' -v b="$BUILD_GRCH38" '$2==b {
  printf "%s\t%d\t%d\t%s|%s|%s\t0\t+\n", $3, $4-1, $4, $1, $6, $7
}' "$ROWS" > "$BED_IN"
log "lifting $(wc -l < "$BED_IN") GRCh38 sites -> hs1"

BED_OUT="$TMP/hs1_sites.bed"
CrossMap bed "$CHAIN" "$BED_IN" "$BED_OUT" || die "CrossMap bed failed"

# Parse lifted BED: pos = end (1-based); complement alleles when the lift flipped to '-'.
awk -F'\t' -v build="$BUILD_HS1" '
  function comp(b){ return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b }
  {
    chrom=$1; pos=$3; strand=$6;
    n=split($4, a, "|"); name=a[1]; anc=a[2]; der=a[3];
    if (strand=="-") { anc=comp(anc); der=comp(der) }
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", name, build, chrom, pos, strand, anc, der;
  }' "$BED_OUT" >> "$ROWS"

# ── assemble the asset ─────────────────────────────────────────────────────────────
DICT="$ASSETS/dictionary.tsv"
{
  printf "name\tbuild\tchrom\tposition\tstrand\tancestral\tderived\n"
  # Sort by (name, build); dedup exact rows. First coord per (name,build) wins in the loader,
  # so a stable sort keeps the build deterministic.
  sort -t$'\t' -k1,1 -k2,2 -k4,4n -u "$ROWS"
} > "$DICT"

ALIASES="$ASSETS/aliases.tsv"
{
  printf "alias\tcanonical\n"
  printf "# YBrowse lists synonymous SNP names as separate rows at the same locus, so each\n"
  printf "# name is already a direct dictionary entry. Add manual alias overrides below if needed.\n"
} > "$ALIASES"

log "stage 2 complete: $(( $(wc -l < "$DICT") - 1 )) rows -> $DICT"
log "  entries by build:"; awk -F'\t' 'NR>1 {c[$2]++} END {for (b in c) printf "    %s: %d\n", b, c[b]}' "$DICT" >&2
