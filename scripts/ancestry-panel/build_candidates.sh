#!/usr/bin/env bash
# Depth sweep — STEP 2 (fast): for each site-count cap, subset the WIDE 1000G matrix (STEP 1) to the
# top-N Fst-ranked sites and build candidate fine-admixture + PCA assets. No re-slicing of the BCF —
# each cap is just a row filter of the one wide matrix. Outputs land in $WORK/sweep as
# ancestry_{freq_global,pca}.<N>.bin, ready for score_modern_from_tsv.
#
#   CAPS="20000 100000 200000" ./build_candidates.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool cargo

CAPS="${CAPS:-20000 100000 200000}"
SWEEP="$WORK/sweep"
WIDE_SITES_TSV="$SWEEP/wide_sites.tsv"       # contig position ref alt fst af_*  (fst = col 5, ranking)
WIDE_MATRIX="$SWEEP/1kgp.wide.matrix.tsv.gz"
WIDE_SAMPLES="$SWEEP/1kgp.wide.samples.txt"
POPS="$SWEEP/pops.modern.tsv"
for f in "$WIDE_SITES_TSV" "$WIDE_MATRIX" "$WIDE_SAMPLES" "$POPS"; do
  [[ -s "$f" ]] || die "missing $f (run sweep_modern_depth.sh first)"
done

for N in $CAPS; do
  log "=== cap N=$N ==="
  SITESET="$SWEEP/sites.top${N}.tsv"          # CHROM<TAB>POS of the top-N Fst sites
  MATRIX_N="$SWEEP/1kgp.top${N}.matrix.tsv.gz"
  FINE_N="$SWEEP/ancestry_freq_global.${N}.bin"
  PCA_N="$SWEEP/ancestry_pca.${N}.bin"

  # Top-N sites by Fst (desc), emitted as CHROM<TAB>POS. Use awk `NR<=n` (not `head`, which closes
  # the pipe early and trips SIGPIPE under `set -o pipefail`) to take the top N.
  awk 'NR>1{print $5"\t"$1"\t"$2}' "$WIDE_SITES_TSV" | sort -k1,1gr | awk -v n="$N" 'NR<=n' \
    | awk '{print $2"\t"$3}' | sort -k1,1 -k2,2n -u > "$SITESET"
  log "  selected $(wc -l < "$SITESET" | tr -d ' ') sites"

  # Filter the wide matrix to those sites (CHROM=col1, POS=col2).
  if [[ -s "$MATRIX_N" ]]; then
    log "  have $(basename "$MATRIX_N")"
  else
    gzip -dc "$WIDE_MATRIX" | awk 'NR==FNR{keep[$1"\t"$2]=1;next} (($1"\t"$2) in keep)' "$SITESET" - \
      | gzip > "$MATRIX_N"
    log "  matrix rows: $(gzip -dc "$MATRIX_N" | wc -l | tr -d ' ')"
  fi

  # Fine-admixture AF panel (26-pop) and modern PCA at this cap.
  log "  fine-panel -> $(basename "$FINE_N")"
  cargo run --release -q -p navigator-panelbuild -- fine-panel \
    --matrix "$MATRIX_N" --samples "$WIDE_SAMPLES" --pops "$POPS" \
    --out "$FINE_N" --min-call-rate "$MIN_CALL_RATE"

  log "  pca (k=$PCA_COMPONENTS) -> $(basename "$PCA_N")"
  cargo run --release -q -p navigator-panelbuild -- pca \
    --matrix "$MATRIX_N" --samples "$WIDE_SAMPLES" --pops "$POPS" \
    --out "$PCA_N" --components "$PCA_COMPONENTS" --min-call-rate "$MIN_CALL_RATE"

  log "  cap $N done: $(ls -lh "$FINE_N" | awk '{print $5}') fine, $(ls -lh "$PCA_N" | awk '{print $5}') pca"
done

log "STEP 2 complete. Candidates in $SWEEP:"
ls -lh "$SWEEP"/ancestry_freq_global.*.bin "$SWEEP"/ancestry_pca.*.bin 2>/dev/null >&2 || true
log "Score with: cargo run --release -q -p navigator-panelbuild --example score_modern_from_tsv -- \\"
log "              $SWEEP/ancestry_freq_global.<N>.bin  $SWEEP/james.wgs.dosage.tsv  $SWEEP/ancestry_pca.<N>.bin"
