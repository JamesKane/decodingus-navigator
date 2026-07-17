#!/usr/bin/env bash
# Stage 5 — build all global ancestry + IBD assets from the merged matrices + pop map, then the
# integrity manifest. Produces (in $ASSETS):
#   ancestry_pca_<build>.bin          modern PCA (PC1×PC2 scatter reference)
#   ancestry_freq_global_<build>.bin  fine per-population AF (fine admixture)
#   ancestry_freq_ancient_<build>.bin deep-source AF: WHG/ANF/Steppe (deep ancestry) + its gates
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

# (2a) Ascertainment manifest (Option A′): the consumer-array-assayed sites, on CHM13. Allele-freq
#      admixture is only valid where sample and reference share ascertainment; the AADR/1240k universe
#      includes capture sites arrays don't assay, and on those the deep estimate is unstable across
#      data sources. Map the array's rsIDs to CHM13 via the stage-02 1240k liftover BED (name field is
#      "rsid|ref|alt", CHM13 pos is the BED end). See docs/design/ancient-ancestry-rebuild.md §4.
BED_1240K="$TMP/1240k_sites.${BUILD}.bed"
ASCERTAIN_ARG=()
if [[ -n "$CHIP_MANIFEST" ]]; then
  [[ -s "$CHIP_MANIFEST" ]] || die "CHIP_MANIFEST=$CHIP_MANIFEST not found"
  [[ -s "$BED_1240K" ]] || die "missing $BED_1240K (run 02_liftover_panel_sites.sh)"
  case "$CHIP_MANIFEST" in *.gz) rd=(gzip -dc) ;; *) rd=(cat) ;; esac
  "${rd[@]}" "$CHIP_MANIFEST" | awk 'BEGIN{FS="[ \t,]"} $0 !~ /^#/ && $1 ~ /^rs/ {print $1}' | sort -u > "$TMP/chip_rsids.txt"
  awk -F'\t' 'NR==FNR{keep[$1]=1; next} { split($4,a,"|"); if (a[1] in keep) print $1"\t"$3 }' \
    "$TMP/chip_rsids.txt" "$BED_1240K" | sort -k1,1 -k2,2n -u > "$ASCERTAIN_SITES"
  [[ -s "$ASCERTAIN_SITES" ]] || die "no array rsIDs matched the 1240k panel — is CHIP_MANIFEST an rsID list?"
  log "ascertainment: $(wc -l < "$TMP/chip_rsids.txt") array rsIDs -> $(wc -l < "$ASCERTAIN_SITES") CHM13 sites -> $ASCERTAIN_SITES"
  ASCERTAIN_ARG=(--ascertain-sites "$ASCERTAIN_SITES")
else
  log "WARN: CHIP_MANIFEST unset — building the FULL (unascertained) ancient panel. It will NOT pass the"
  log "WARN: §3.4 stability gate; step (6) below is expected to fail. Set CHIP_MANIFEST to ship deep ancestry."
fi

# (2) Ancient deep-source AF panel (WHG/ANF/Steppe) — the deep-ancestry asset.
#     Its own builder, NOT a column subset of the fine panel: the fine-panel builder writes 0.0 for a
#     population with no called sample at a site, which is indistinguishable from a real "alt absent".
#     The 1000G pops are called nearly everywhere so that barely hurts them, but ancient genomes are
#     sparse — most sites would enter the mixture as fake zero-frequency evidence and the fit would
#     track missingness instead of ancestry. Here a site survives only if EVERY source has
#     >= $ANCIENT_MIN_CALLED real calls, and (with A′) only if the site is on the array manifest.
log "panelbuild ancient-panel ($ANCIENT_COMPONENTS, >=$ANCIENT_MIN_CALLED calls/source) -> $ANCIENT_OUT"
cargo run --release -q -p navigator-panelbuild -- ancient-panel \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPMAP" \
  --components "$ANCIENT_COMPONENTS" --min-called "$ANCIENT_MIN_CALLED" \
  "${ASCERTAIN_ARG[@]}" \
  --out "$ANCIENT_OUT" --sites-tsv "$TMP/ancient_sites.tsv"

# (3) Fine-population AF panel (fine admixture over the full labelled pop set).
log "panelbuild fine-panel -> $FINE_OUT"
cargo run --release -q -p navigator-panelbuild -- fine-panel \
  --matrix "$MATRICES" --samples "$SAMPLES" --pops "$POPMAP" \
  --out "$FINE_OUT" --min-call-rate "$MIN_CALL_RATE"

# (4) IBD genetic map (recombination map, GRCh38 -> CHM13). Best-effort — IBD falls back to uniform.
build_genetic_map "$GMAP_OUT" || log "WARN: genetic map not built (IBD will use uniform 1 cM/Mb)"

# (5) Chip-compatible multi-build IBD panel. Best-effort.
build_ibd_panel "$IBD_PANEL_OUT" || log "WARN: IBD panel not built (IBD / chip matching unavailable)"

# (6) VALIDATION GATE for the deep-ancestry asset. Do NOT publish an ancient panel that fails this.
#     It simulates individuals from populations whose ancestry is known and checks what comes back:
#     the estimator must round-trip mixtures it was given, put a NW-European near Steppe 40-55 /
#     ANF 25-40 / WHG 10-25, and REJECT samples the three sources cannot express. The previous
#     ancient asset shipped fabricated numbers precisely because nobody ran this.
log "validating $ANCIENT_OUT (deep-ancestry gates)"
cargo run --release -q -p navigator-panelbuild -- validate-ancient \
  --ancient "$ANCIENT_OUT" --reference "$FINE_OUT" \
  || die "deep-ancestry validation FAILED — do not publish this asset"

# (7) Asset integrity manifest (sha256 of every *_<build>.bin) — run last so it covers everything.
log "panelbuild manifest -> $MANIFEST"
cargo run --release -q -p navigator-panelbuild -- manifest --dir "$ASSETS" --build "$BUILD" --out "$MANIFEST" \
  || die "panelbuild manifest failed"

log "stage 5 complete. Assets in $ASSETS:"
ls -lh "$PANEL_OUT" "$PCA_OUT" "$FINE_OUT" "$ANCIENT_OUT" "$GMAP_OUT" "$IBD_PANEL_OUT" "$MANIFEST" 2>/dev/null >&2 || true
