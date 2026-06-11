#!/usr/bin/env bash
# Stage 3 (optional) — restrict the full dictionary to the SNP set a specific BISDNA chip
# probes, producing a small manifest that can be checked into the repo so import works
# offline with no 2M-row asset. Also reports which chip names the dictionary couldn't resolve.
#
#   ./03_restrict_panel.sh <bisdna-results.txt> [out.tsv]
# Default out: $ASSETS/chromo2-panel.tsv
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool awk

BISDNA="${1:-}"; OUT="${2:-$ASSETS/chromo2-panel.tsv}"
[[ -s "$BISDNA" ]] || die "usage: $0 <bisdna-results.txt> [out.tsv]"
DICT="$ASSETS/dictionary.tsv"
[[ -s "$DICT" ]] || die "missing $DICT (run 02_build.sh)"

# Chip SNP names: column 1 of the rows after the `SNPID\tgenotype\tresult` header.
NAMES="$TMP/chip_names.txt"
awk -F'\t' '
  seen { if (NF>=1 && $1!="") print $1; next }
  tolower($1)=="snpid" && tolower($2)=="genotype" { seen=1 }
' "$BISDNA" | sort -u > "$NAMES"
log "chip probes $(wc -l < "$NAMES") distinct SNP names"

# Keep dictionary rows whose name (lowercased) is in the chip set. Case-insensitive join.
{
  head -1 "$DICT"
  awk -F'\t' '
    NR==FNR { want[tolower($1)]=1; next }
    FNR==1  { next }                          # dictionary header
    (tolower($1) in want)
  ' "$NAMES" "$DICT"
} > "$OUT"

# Report chip names with NO dictionary entry on any build (unresolved — can never be placed).
RESOLVED="$TMP/resolved_names.txt"
awk -F'\t' 'NR>1 {print tolower($1)}' "$OUT" | sort -u > "$RESOLVED"
UNRESOLVED="$TMP/unresolved_names.txt"
comm -23 <(awk '{print tolower($0)}' "$NAMES" | sort -u) "$RESOLVED" > "$UNRESOLVED"

log "stage 3 complete: $(( $(wc -l < "$OUT") - 1 )) rows -> $OUT"
log "  resolved chip names: $(wc -l < "$RESOLVED") / $(wc -l < "$NAMES")"
if [[ -s "$UNRESOLVED" ]]; then
  log "  UNRESOLVED chip names ($(wc -l < "$UNRESOLVED")) listed in $UNRESOLVED (first 10):"
  head -10 "$UNRESOLVED" >&2
fi
