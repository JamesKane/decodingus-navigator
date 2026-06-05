#!/usr/bin/env bash
# Stage 3 — select the AIMs panel within the 1240k-restricted, CHM13-lifted site set.
#
# Restrict the 1000G-on-CHM13 VCFs to 1240k∩CHM13 (stage 2), then run `navigator-panelbuild
# panel` (Fst-ranked AIMs from per-super-pop INFO AC/AN). The result is ancient-compatible by
# construction: every panel site exists in 1240k, so AADR + modern + the user's sample all
# overlap. Output: $PANEL_OUT (+ the panel sites regions file consumed by later stages).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool bcftools
require_tool cargo

SITES="$TMP/1240k_sites.${BUILD}.tsv"
[[ -s "$SITES" ]] || die "missing $SITES (run 02_liftover_panel_sites.sh)"
FILTERED="$TMP/1kgp-chm13-1240k"; mkdir -p "$FILTERED"

# Restrict each 1000G-CHM13 VCF to the lifted 1240k sites (keeps the INFO AF the panel needs).
log "restricting 1000G-CHM13 to 1240k sites"
for vcf in "$KGP_CHM13_DIR"/*.vcf.gz; do
  [[ -e "$vcf" ]] || die "no 1000G-CHM13 VCFs in $KGP_CHM13_DIR (run 01_fetch.sh)"
  out="$FILTERED/$(basename "$vcf")"
  [[ -s "$out" ]] && { log "  have $(basename "$out")"; continue; }
  bcftools view -R "$SITES" -Oz -o "$out" "$vcf" && tabix -f -p vcf "$out"
done

# Build the AIMs panel from the restricted VCF dir.
log "navigator-panelbuild panel -> $PANEL_OUT"
cargo run --release -q -p navigator-panelbuild -- panel \
  --vcf-dir "$FILTERED" \
  --out "$PANEL_OUT" \
  --max-sites "$MAX_SITES" \
  --min-fst "$MIN_FST" \
  --tsv "$TMP/panel_sites.tsv"

# The chosen panel's (CHROM,POS) regions — every downstream matrix is cut to exactly these.
awk 'NR>1 { printf "%s\t%s\n", $1, $2 }' "$TMP/panel_sites.tsv" | sort -k1,1 -k2,2n -u \
  > "$TMP/panel_regions.${BUILD}.tsv"

log "stage 3 complete: $(wc -l < "$TMP/panel_regions.${BUILD}.tsv") AIMs -> $PANEL_OUT"
