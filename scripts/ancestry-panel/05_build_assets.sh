#!/usr/bin/env bash
# Stage 5 — build all global ancestry + IBD assets from the merged matrices + pop map, then the
# integrity manifest. Produces (in $ASSETS):
#   ancestry_pca_<build>.bin          modern PCA (PC1×PC2 scatter reference)
#   ancestry_pca_ancient_<build>.bin  PCA with ancient deep components projected (GMM/nMonte)
#   ancestry_freq_global_<build>.bin  fine per-population AF (fine admixture)
#   genetic_map_<build>.bin           IBD recombination map (bp->cM)
#   ibd_panel_<build>.bin             chip-compatible multi-build IBD SNP panel
#   ancestry_manifest_<build>.json    sha256 of every *_<build>.bin (clients verify against it)
# (The super-pop AF panel ancestry_panel_<build>.bin is built in stage 03.)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool cargo

MATRICES="$(cat "$TMP/matrices.list")"; SAMPLES="$(cat "$TMP/samples.list")"
POPMAP="$TMP/pops.${BUILD}.tsv"
[[ -n "$MATRICES" && -s "$POPMAP" ]] || die "missing matrices/pop map (run 04_build_matrices.sh)"

# Projection basis = the modern reference pops (everything that is NOT an AADR component-map label —
# i.e. the 1000G fine pops + SGDP regional pops). The ancient deep components are projected onto it.
BASIS="$TMP/basis_pops.txt"
awk -F'\t' '$0 !~ /^#/ && NF>=2 { print $2 }' "$AADR_COMPONENT_MAP" | sort -u > "$TMP/ancient_components.txt"
cut -f2 "$POPMAP" | sort -u | grep -vxF -f "$TMP/ancient_components.txt" > "$BASIS" || true
log "projection basis: $(wc -l < "$BASIS") modern pops; $(wc -l < "$TMP/ancient_components.txt") deep components"

# Modern-only pop map (samples whose pop is in the modern basis) — for the clean scatter PCA.
MODERN_POPS="$TMP/modern_pops.tsv"
awk 'NR==FNR{b[$1]=1;next} ($2 in b){print}' "$BASIS" "$POPMAP" > "$MODERN_POPS"

# (1) Modern PCA — basis pops only (the PC1×PC2 scatter reference set).
log "panelbuild pca (modern, k=$PCA_COMPONENTS) -> $PCA_OUT"
cargo run --release -q -p navigator-panelbuild -- pca \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$MODERN_POPS" \
  --out "$PCA_OUT" --components "$PCA_COMPONENTS" --min-call-rate "$MIN_CALL_RATE"

# (2) Ancient PCA — full pop set, deep components PROJECTED onto the modern basis (the app's GMM /
#     nMonte deep-ancestry classification prefers this asset; falls back to the modern PCA).
log "panelbuild pca (ancient projected, k=$PCA_COMPONENTS) -> $PCA_ANCIENT_OUT"
cargo run --release -q -p navigator-panelbuild -- pca \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPMAP" --basis-pops "$BASIS" \
  --out "$PCA_ANCIENT_OUT" --components "$PCA_COMPONENTS" --min-call-rate "$MIN_CALL_RATE"

# (3) Fine-population AF panel (fine admixture over the full labelled pop set).
log "panelbuild fine-panel -> $FINE_OUT"
cargo run --release -q -p navigator-panelbuild -- fine-panel \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPMAP" \
  --out "$FINE_OUT" --min-call-rate "$MIN_CALL_RATE"

# (4) IBD genetic map (recombination map, GRCh38 -> CHM13). Best-effort — IBD falls back to uniform.
build_genetic_map "$GMAP_OUT" || log "WARN: genetic map not built (IBD will use uniform 1 cM/Mb)"

# (5) Chip-compatible multi-build IBD panel. Best-effort.
build_ibd_panel "$IBD_PANEL_OUT" || log "WARN: IBD panel not built (IBD / chip matching unavailable)"

# (6) Asset integrity manifest (sha256 of every *_<build>.bin) — run last so it covers everything.
log "panelbuild manifest -> $MANIFEST"
cargo run --release -q -p navigator-panelbuild -- manifest --dir "$ASSETS" --build "$BUILD" --out "$MANIFEST" \
  || die "panelbuild manifest failed"

log "stage 5 complete. Assets in $ASSETS:"
ls -lh "$PANEL_OUT" "$PCA_OUT" "$PCA_ANCIENT_OUT" "$FINE_OUT" "$GMAP_OUT" "$IBD_PANEL_OUT" "$MANIFEST" 2>/dev/null >&2 || true
