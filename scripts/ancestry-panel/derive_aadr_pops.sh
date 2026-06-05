#!/usr/bin/env bash
# Helper: emit `sample<TAB>component` for AADR samples, by mapping each individual's
# AADR "Group ID" through the curated label->component map. Only samples that (a) are in
# the extracted matrix sample list and (b) map to a component are emitted; everything else
# is dropped (most AADR individuals aren't useful source-population references).
#
#   derive_aadr_pops.sh <raw-dir> <component-map.tsv> <matrix-samples.txt>
set -euo pipefail
RAW="$1"; MAP="$2"; SAMPLES="$3"
ANNO="$(ls "$RAW"/*.anno 2>/dev/null | head -1 || true)"
[[ -n "$ANNO" ]] || { echo "WARN: no AADR .anno in $RAW" >&2; exit 0; }

# .anno is a large TSV; column 1 is the individual ID, and a "Group ID" / "Population"
# column holds the label. Resolve that column by header, then join via the map.
awk -F'\t' -v mapf="$MAP" -v sampf="$SAMPLES" '
  BEGIN {
    while ((getline l < mapf) > 0) { if (l ~ /^#/ || l=="") continue; split(l,m,"\t"); comp[m[1]]=m[2] }
    while ((getline s < sampf) > 0) { keep[s]=1 }
  }
  NR==1 {
    for (i=1;i<=NF;i++) { h=tolower($i); if (h ~ /group id/ || h=="population") gcol=i }
    if (!gcol) gcol=2; next
  }
  { id=$1; grp=$gcol; if (id in keep && grp in comp) printf "%s\t%s\n", id, comp[grp] }
' "$ANNO"
