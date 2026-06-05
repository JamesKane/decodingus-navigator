#!/usr/bin/env bash
# Stage 5 — build the global ancestry assets from the merged matrices + pop map.
#
# Produces, on CHM13: the global PCA loadings+centroids ($PCA_OUT, feeds estimate_pca_gmm
# and estimate_nmonte) and the global per-population AF panel ($FINE_OUT, feeds the global
# admixture). The AF super-pop panel ($PANEL_OUT) was already built in stage 3. Writes a
# provenance manifest with checksums + source dataset versions.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool cargo

MATRICES="$(cat "$TMP/matrices.list")"; SAMPLES="$(cat "$TMP/samples.list")"
POPMAP="$TMP/pops.${BUILD}.tsv"
[[ -n "$MATRICES" && -s "$POPMAP" ]] || die "missing matrices/pop map (run 04_build_matrices.sh)"

# Global PCA loadings + per-population centroids/variances.
# NOTE: this builds the PCA over ALL labelled samples. Projection-mode (basis = modern,
# project ancient) is a navigator-panelbuild refinement (see AncestryAnalysis.md) — until it
# lands, exclude very-low-coverage ancient samples from the pop map to limit axis distortion.
log "navigator-panelbuild pca (k=$PCA_COMPONENTS) -> $PCA_OUT"
cargo run --release -q -p navigator-panelbuild -- pca \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPMAP" \
  --out "$PCA_OUT" --components "$PCA_COMPONENTS" --min-call-rate "$MIN_CALL_RATE"

# Global per-population allele-frequency panel (fine admixture over all labelled pops).
log "navigator-panelbuild fine-panel -> $FINE_OUT"
cargo run --release -q -p navigator-panelbuild -- fine-panel \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPMAP" \
  --out "$FINE_OUT" --min-call-rate "$MIN_CALL_RATE"

# Provenance manifest (checksums + source versions) — published alongside the assets.
log "writing manifest $MANIFEST"
{
  printf '{\n'
  printf '  "build": "%s",\n' "$BUILD"
  printf '  "assetVersion": %s,\n' "$ASSET_VERSION"
  printf '  "panelParams": { "maxSites": %s, "minFst": %s, "pcaComponents": %s },\n' "$MAX_SITES" "$MIN_FST" "$PCA_COMPONENTS"
  printf '  "sources": { "aadr": "%s/%s", "kgpChm13": "1KGP-CHM13v2.0", "hgdp1kg": "gnomAD-v3", "sgdp": "SGDP" },\n' "$AADR_DATASET" "$AADR_VERSION"
  printf '  "assets": {\n'
  printf '    "panel": { "file": "%s", "sha256": "%s" },\n' "$(basename "$PANEL_OUT")" "$(sha256_of "$PANEL_OUT")"
  printf '    "pca":   { "file": "%s", "sha256": "%s" },\n' "$(basename "$PCA_OUT")"   "$(sha256_of "$PCA_OUT")"
  printf '    "freq":  { "file": "%s", "sha256": "%s" }\n'  "$(basename "$FINE_OUT")"  "$(sha256_of "$FINE_OUT")"
  printf '  }\n}\n'
} > "$MANIFEST"

log "stage 5 complete. Assets in $ASSETS:"
ls -lh "$PANEL_OUT" "$PCA_OUT" "$FINE_OUT" "$MANIFEST" >&2
