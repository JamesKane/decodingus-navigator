#!/usr/bin/env bash
# Build the FINAL fine-admixture + PCA assets at the chosen depth cap from the combined modern
# reference set: 1000G (26 pops) + the AADR present-day CONTINENTAL-European pops (add_continental_pops.sh).
# Default cap = 200000 (the wide set); set CAP to sub-select the top-N Fst sites.
#   CAP=200000 ./build_final_fine_pca.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool cargo

CAP="${CAP:-200000}"
SWEEP="$WORK/sweep"
WIDE_SITES_TSV="$SWEEP/wide_sites.tsv"
KGP_MATRIX="$SWEEP/1kgp.wide.matrix.tsv.gz";        KGP_SAMPLES="$SWEEP/1kgp.wide.samples.txt"
CONT_MATRIX="$SWEEP/aadr_continental.matrix.tsv.gz"; CONT_SAMPLES="$SWEEP/aadr_continental.samples.txt"
POPS_COMBINED="$SWEEP/pops.combined.tsv"
for f in "$WIDE_SITES_TSV" "$KGP_MATRIX" "$KGP_SAMPLES" "$CONT_MATRIX" "$CONT_SAMPLES" \
         "$RAW/1kgp.pops.tsv" "$SWEEP/pops.continental.tsv"; do
  [[ -s "$f" ]] || die "missing $f"
done

# Combined pop map (1000G fine codes + continental codes).
cat "$RAW/1kgp.pops.tsv" "$SWEEP/pops.continental.tsv" > "$POPS_COMBINED"
log "combined pop map: $(wc -l < "$POPS_COMBINED" | tr -d ' ') samples; codes: $(cut -f2 "$POPS_COMBINED" | sort -u | tr '\n' ' ')"

# Top-CAP sites (Fst-ranked). CAP=200000 keeps the whole wide set.
SITESET="$SWEEP/sites.final.top${CAP}.tsv"
awk 'NR>1{print $5"\t"$1"\t"$2}' "$WIDE_SITES_TSV" | sort -k1,1gr | awk -v n="$CAP" 'NR<=n' \
  | awk '{print $2"\t"$3}' | sort -k1,1 -k2,2n -u > "$SITESET"
log "final site set: $(wc -l < "$SITESET" | tr -d ' ') sites (cap $CAP)"

# Filter both matrices to the site set.
KGP_F="$SWEEP/1kgp.final.matrix.tsv.gz"; CONT_F="$SWEEP/cont.final.matrix.tsv.gz"
for pair in "$KGP_MATRIX:$KGP_F" "$CONT_MATRIX:$CONT_F"; do
  src="${pair%%:*}"; dst="${pair##*:}"
  gzip -dc "$src" | awk 'NR==FNR{keep[$1"\t"$2]=1;next} (($1"\t"$2) in keep)' "$SITESET" - | gzip > "$dst"
  log "  filtered $(basename "$src") -> $(gzip -dc "$dst" | wc -l | tr -d ' ') rows"
done

FINE_OUT_C="$SWEEP/ancestry_freq_global.cont.${CAP}.bin"
PCA_OUT_C="$SWEEP/ancestry_pca.cont.${CAP}.bin"
MATRICES="$KGP_F,$CONT_F"; SAMPLES="$KGP_SAMPLES,$CONT_SAMPLES"

log "fine-panel (1000G + continental) -> $(basename "$FINE_OUT_C")"
cargo run --release -q -p navigator-panelbuild -- fine-panel \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPS_COMBINED" \
  --out "$FINE_OUT_C" --min-call-rate "$MIN_CALL_RATE"

log "pca (1000G + continental, k=$PCA_COMPONENTS) -> $(basename "$PCA_OUT_C")"
cargo run --release -q -p navigator-panelbuild -- pca \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPS_COMBINED" \
  --out "$PCA_OUT_C" --components "$PCA_COMPONENTS" --min-call-rate "$MIN_CALL_RATE"

log "done: $(ls -lh "$FINE_OUT_C" | awk '{print $5}') fine, $(ls -lh "$PCA_OUT_C" | awk '{print $5}') pca"
log "score: cargo run --release -q -p navigator-panelbuild --example score_modern_from_tsv -- \\"
log "         $FINE_OUT_C $SWEEP/james.wgs.dosage.tsv $PCA_OUT_C"
