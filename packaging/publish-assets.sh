#!/usr/bin/env bash
# Publish the bundled ancestry/IBD reference assets to the GitHub release that the installer staging
# pulls from (see packaging/stage-assets.sh + docs/design/packaging-and-release.md). Run this once
# (and again whenever the assets are regenerated) from a machine with the built assets in
# ~/.decodingus/ancestry and an authenticated `gh`.
#
#   ./packaging/publish-assets.sh [build]      # build defaults to chm13v2.0 → release tag assets-<build>
#
# The release is a data store, not a source release: it holds the 10 `*_<build>.bin` files plus
# `ancestry_manifest_<build>.json`. `stage-assets.sh` downloads them by manifest at package time.
set -euo pipefail

BUILD="${1:-chm13v2.0}"
REPO="${NAVIGATOR_ASSET_REPO:-JamesKane/decodingus-navigator}"
SRC="${NAVIGATOR_ASSET_SRC:-$HOME/.decodingus/ancestry}"
TAG="assets-${BUILD}"

# The bundle: every per-build .bin + the manifest. Match stage-assets.sh's expectations.
files=()
for f in "$SRC"/*"${BUILD}".bin "$SRC/ancestry_manifest_${BUILD}.json"; do
  [ -e "$f" ] || { echo "publish-assets: missing $f" >&2; exit 1; }
  files+=("$f")
done

echo "publish-assets: ${#files[@]} file(s) → $REPO release '$TAG'"
if gh release view "$TAG" --repo "$REPO" >/dev/null 2>&1; then
  echo "publish-assets: release exists — uploading (clobber)"
  gh release upload "$TAG" "${files[@]}" --repo "$REPO" --clobber
else
  echo "publish-assets: creating release"
  gh release create "$TAG" "${files[@]}" --repo "$REPO" \
    --title "Ancestry/IBD assets ($BUILD)" \
    --notes "Reference assets bundled into the offline installer; consumed by \`packaging/stage-assets.sh\`. Regenerable — not source."
fi
echo "publish-assets: done. Set NAVIGATOR_ASSET_RELEASE=$TAG in the release workflow."
