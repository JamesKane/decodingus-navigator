#!/usr/bin/env bash
# Stage the ancestry/IBD reference assets that the installer bundles (the full offline installer —
# see docs/design/packaging-and-release.md, decision A). cargo-packager's `resources` points at
# packaging/staging/ancestry; this populates it before packaging (via before-packaging-command).
#
#   Local builds:  copy from the developer's ~/.decodingus/ancestry/ (already built/downloaded).
#   CI builds:     fetch from the asset CDN by manifest (TODO: wire $NAVIGATOR_ASSET_CDN once the
#                  CDN base is final — the manifest-verified download path already exists in-app).
#
# Idempotent: only copies files that are missing or changed. Never fails the build on a missing
# source (a lean dev package without the big assets still builds; first run just has no bundle).
set -euo pipefail

# Resolve the repo root from this script's location so it works regardless of cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$SCRIPT_DIR/staging/ancestry"
mkdir -p "$STAGE"

SRC="${NAVIGATOR_ASSET_SRC:-$HOME/.decodingus/ancestry}"

# The full Option-A bundle: ancestry panels/PCA/freqs + manifest + genetic map + IBD panel.
PATTERNS=(
  "ancestry_panel_"*.bin
  "ancestry_pca_"*.bin
  "ancestry_pca_ancient_"*.bin
  "ancestry_freq_global_"*.bin
  "ancestry_manifest_"*.json
  "genetic_map_"*.bin
  "ibd_panel_"*.bin
)

if [ -d "$SRC" ]; then
  echo "stage-assets: copying bundled assets from $SRC"
  copied=0
  for pat in "${PATTERNS[@]}"; do
    for f in "$SRC"/$pat; do
      [ -e "$f" ] || continue
      name="$(basename "$f")"
      # Copy when missing or size-differs (cheap freshness check; the app re-verifies via sha256).
      if [ ! -f "$STAGE/$name" ] || [ "$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f")" != "$(stat -f%z "$STAGE/$name" 2>/dev/null || stat -c%s "$STAGE/$name")" ]; then
        cp "$f" "$STAGE/$name"
        copied=$((copied + 1))
      fi
    done
  done
  echo "stage-assets: staged $copied file(s) into $STAGE"
else
  # TODO: CI path — fetch the manifest-listed assets from $NAVIGATOR_ASSET_CDN here.
  echo "stage-assets: WARNING — no asset source at $SRC and no CDN configured; bundling an empty asset set." >&2
fi

# cargo-packager requires the resource dir to exist even if empty.
touch "$STAGE/.staged"
