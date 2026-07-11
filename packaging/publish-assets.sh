#!/usr/bin/env bash
# Publish the bundled ancestry/IBD reference assets to the GitHub release that the installer staging
# pulls from (see packaging/stage-assets.sh + docs/design/packaging-and-release.md). Run this once
# (and again whenever the assets are regenerated) from a machine with the built assets in
# ~/.decodingus/ancestry and an authenticated `gh`.
#
#   ./packaging/publish-assets.sh [build]      # build defaults to chm13v2.0 → release tag assets-<build>
#   ./packaging/publish-assets.sh ysnp         # the full Y-SNP dictionary → release tag assets-ysnp
#
# The release is a data store, not a source release: it holds the 10 `*_<build>.bin` files,
# `ancestry_manifest_<build>.json`, and the HipSTR references (`<build>.hipstr_reference.bed.gz` +
# `GRCh38.hipstr_reference.bed.gz`). `stage-assets.sh` downloads them at package time.
set -euo pipefail

REPO="${NAVIGATOR_ASSET_REPO:-JamesKane/decodingus-navigator}"

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
  else shasum -a 256 "$1" | awk '{print $1}'; fi
}

# ── ysnp mode: publish the full Y-SNP dictionary (~200 MB, too big to bundle) + a small sha256
# manifest the app verifies against (App::ensure_ysnp_dictionary downloads it on first import). ──
if [ "${1:-}" = "ysnp" ]; then
  SRC="${NAVIGATOR_YSNP_SRC:-$HOME/.decodingus/ysnp}"
  TAG="assets-ysnp"
  DICT="$SRC/dictionary.tsv"
  [ -e "$DICT" ] || { echo "publish-assets: missing $DICT (build it with scripts/ysnp-dictionary)" >&2; exit 1; }
  sha="$(sha256_of "$DICT")"
  bytes="$(stat -f%z "$DICT" 2>/dev/null || stat -c%s "$DICT")"
  MANIFEST="$SRC/ysnp_manifest.json"
  cat > "$MANIFEST" <<JSON
{
  "build": "all",
  "generated_at": "",
  "assets": {
    "dictionary.tsv": { "sha256": "$sha", "bytes": $bytes }
  }
}
JSON
  echo "publish-assets: ysnp dictionary ($bytes bytes, sha256 $sha) → $REPO release '$TAG'"
  if gh release view "$TAG" --repo "$REPO" >/dev/null 2>&1; then
    gh release upload "$TAG" "$DICT" "$MANIFEST" --repo "$REPO" --clobber
  else
    gh release create "$TAG" "$DICT" "$MANIFEST" --repo "$REPO" \
      --title "Y-SNP dictionary" \
      --notes "Full YBrowse-derived Y-SNP catalog + sha256 manifest; downloaded on first Y-SNP import by App::ensure_ysnp_dictionary. Regenerable — not source."
  fi
  echo "publish-assets: done ($TAG)."
  exit 0
fi

BUILD="${1:-chm13v2.0}"
SRC="${NAVIGATOR_ASSET_SRC:-$HOME/.decodingus/ancestry}"
STR_SRC="${NAVIGATOR_STR_SRC:-$HOME/.decodingus/str}"
TAG="assets-${BUILD}"

# The bundle: every per-build .bin + the manifest. Match stage-assets.sh's expectations.
files=()
for f in "$SRC"/*"${BUILD}".bin "$SRC/ancestry_manifest_${BUILD}.json"; do
  [ -e "$f" ] || { echo "publish-assets: missing $f" >&2; exit 1; }
  files+=("$f")
done

# STR HipSTR references: the build-specific one + the shared GRCh38 one, by the exact names
# stage-assets.sh downloads from this release ("${BUILD}.hipstr_reference.bed.gz" and
# "GRCh38.hipstr_reference.bed.gz"). Without these in the release, a clean CI runner (no local
# ~/.decodingus/str) silently bundles no STR reference and fresh installs can't do STR calling.
for f in "$STR_SRC/${BUILD}.hipstr_reference.bed.gz" "$STR_SRC/GRCh38.hipstr_reference.bed.gz"; do
  [ -e "$f" ] || { echo "publish-assets: missing $f (build it into $STR_SRC, or set NAVIGATOR_STR_SRC)" >&2; exit 1; }
  files+=("$f")
done

echo "publish-assets: ${#files[@]} file(s) → $REPO release '$TAG'"
if gh release view "$TAG" --repo "$REPO" >/dev/null 2>&1; then
  echo "publish-assets: release exists — uploading (clobber)"
  gh release upload "$TAG" "${files[@]}" --repo "$REPO" --clobber
else
  echo "publish-assets: creating release"
  gh release create "$TAG" "${files[@]}" --repo "$REPO" \
    --title "Ancestry/IBD + STR assets ($BUILD)" \
    --notes "Reference assets bundled into the offline installer; consumed by \`packaging/stage-assets.sh\` (ancestry/IBD .bin + manifest, plus HipSTR references). Regenerable — not source."
fi
echo "publish-assets: done. Set NAVIGATOR_ASSET_RELEASE=$TAG in the release workflow."
