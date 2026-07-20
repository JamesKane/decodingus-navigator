#!/usr/bin/env bash
# Stage the ancestry/IBD reference assets that the installer bundles (the full offline installer —
# see docs/design/packaging-and-release.md, decision A). cargo-packager's `resources` points at
# packaging/staging/ancestry; this populates it before packaging (via before-packaging-command).
#
#   Local builds:  copy from the developer's ~/.decodingus/ancestry/ (already built/downloaded).
#   CI builds:     download from a GitHub release (NAVIGATOR_ASSET_RELEASE tag) by manifest.
#
# Idempotent: only fetches files that are missing. Never fails the build when *no* source is
# configured (a lean dev package without the big assets still builds; first run just has no bundle),
# but a configured release that fails to download / verify IS fatal (don't ship a broken bundle).
set -euo pipefail

# Resolve the repo root from this script's location so it works regardless of cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$SCRIPT_DIR/staging/ancestry"
mkdir -p "$STAGE"

SRC="${NAVIGATOR_ASSET_SRC:-$HOME/.decodingus/ancestry}"

# Copy files matching any of the given glob patterns from a source dir into a stage dir, when missing
# or size-differs (a cheap freshness check; the app re-verifies via sha256). Echoes the copied count.
stage_from_dir() {  # src_dir  stage_dir  pattern...
  local src="$1" stage="$2"; shift 2
  local copied=0 pat f name
  for pat in "$@"; do
    for f in "$src"/$pat; do
      [ -e "$f" ] || continue
      name="$(basename "$f")"
      if [ ! -f "$stage/$name" ] || [ "$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f")" != "$(stat -f%z "$stage/$name" 2>/dev/null || stat -c%s "$stage/$name")" ]; then
        cp "$f" "$stage/$name"
        copied=$((copied + 1))
      fi
    done
  done
  echo "$copied"
}

# GitHub-release source (CI). The repo is public, so release assets download over plain HTTPS with
# no token — works on every runner and inside the manylinux build container alike.
ASSET_RELEASE="${NAVIGATOR_ASSET_RELEASE:-}"
ASSET_REPO="${NAVIGATOR_ASSET_REPO:-JamesKane/decodingus-navigator}"
ASSET_BUILD="${NAVIGATOR_ASSET_BUILD:-chm13v2.0}"

# The full Option-A bundle: ancestry panels/PCA/freqs + manifest + genetic map + IBD panel.
PATTERNS=(
  "ancestry_panel_"*.bin
  "ancestry_pca_"*.bin
  "ancestry_pca_ancient_"*.bin
  "ancestry_freq_global_"*.bin
  "ancestry_qpadm_"*.bin
  "ancestry_manifest_"*.json
  "genetic_map_"*.bin
  "ibd_panel_"*.bin
)

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
  else shasum -a 256 "$1" | awk '{print $1}'; fi
}

if [ -d "$SRC" ]; then
  echo "stage-assets: copying bundled assets from $SRC"
  copied="$(stage_from_dir "$SRC" "$STAGE" "${PATTERNS[@]}")"
  echo "stage-assets: staged $copied file(s) into $STAGE"
elif [ -n "$ASSET_RELEASE" ]; then
  base="https://github.com/$ASSET_REPO/releases/download/$ASSET_RELEASE"
  manifest="ancestry_manifest_${ASSET_BUILD}.json"
  echo "stage-assets: fetching assets from release '$ASSET_RELEASE' ($base)"
  # The manifest lists every data asset with its sha256; fetch it first, then the assets it names.
  curl -fSL --retry 3 --retry-delay 2 -o "$STAGE/$manifest" "$base/$manifest"
  names="$(grep -oE '"[A-Za-z0-9_.-]+\.bin"' "$STAGE/$manifest" | tr -d '"' | sort -u)"
  [ -n "$names" ] || { echo "stage-assets: ERROR — no .bin assets listed in $manifest" >&2; exit 1; }
  fetched=0
  for name in $names; do
    if [ ! -f "$STAGE/$name" ]; then
      echo "  downloading $name"
      curl -fSL --retry 3 --retry-delay 2 -o "$STAGE/$name" "$base/$name"
    fi
    # Corruption guard: the downloaded file's hash must appear in the manifest (the app then does the
    # authoritative per-file verification at first run via AssetManifest).
    sha="$(sha256_of "$STAGE/$name")"
    grep -q "$sha" "$STAGE/$manifest" || {
      echo "stage-assets: ERROR — checksum mismatch for $name (got $sha; not in $manifest)" >&2
      exit 1
    }
    fetched=$((fetched + 1))
  done
  echo "stage-assets: fetched manifest + $fetched verified asset(s) into $STAGE"
else
  echo "stage-assets: WARNING — no asset source (NAVIGATOR_ASSET_SRC dir or NAVIGATOR_ASSET_RELEASE tag); bundling an empty asset set." >&2
fi

# --- chrY private-Y filtering masks (CHM13) -------------------------------------------------------
# Bundled like the ancestry panels, but the SOURCE is the checked-in repo `assets/masks/` (gzipped
# BEDs, small enough to live in git) rather than ~/.decodingus — so packaging is reproducible from a
# clean checkout. Seeded to ~/.decodingus/masks/ on first run (navigator_app::seed_bundled_masks).
MASK_STAGE="$SCRIPT_DIR/staging/masks"
mkdir -p "$MASK_STAGE"
MASK_SRC="${NAVIGATOR_MASK_SRC:-$SCRIPT_DIR/../assets/masks}"
if [ -d "$MASK_SRC" ]; then
  mask_copied="$(stage_from_dir "$MASK_SRC" "$MASK_STAGE" "chrY_"*.bed.gz "chrY_"*.bed)"
  echo "stage-assets: staged $mask_copied mask file(s) into $MASK_STAGE"
else
  echo "stage-assets: WARNING — no mask source at $MASK_SRC; bundling no chrY masks." >&2
fi

# --- STR reference (HipSTR BEDs) ------------------------------------------------------------------
# Bundled like the ancestry panels (too big for git): sourced from the dev's ~/.decodingus/str (or
# NAVIGATOR_STR_SRC), else fetched from the asset release by known name. Seeded to ~/.decodingus/str/
# on first run (navigator_app::seed_bundled_str).
STR_STAGE="$SCRIPT_DIR/staging/str"
mkdir -p "$STR_STAGE"
STR_SRC="${NAVIGATOR_STR_SRC:-$HOME/.decodingus/str}"
if [ -d "$STR_SRC" ]; then
  str_copied="$(stage_from_dir "$STR_SRC" "$STR_STAGE" "*.hipstr_reference.bed.gz")"
  echo "stage-assets: staged $str_copied STR reference(s) into $STR_STAGE"
elif [ -n "$ASSET_RELEASE" ]; then
  base="https://github.com/$ASSET_REPO/releases/download/$ASSET_RELEASE"
  for name in "${ASSET_BUILD}.hipstr_reference.bed.gz" "GRCh38.hipstr_reference.bed.gz"; do
    if [ ! -f "$STR_STAGE/$name" ]; then
      curl -fSL --retry 3 --retry-delay 2 -o "$STR_STAGE/$name" "$base/$name" \
        && echo "  downloaded $name" || echo "  (no $name in release; skipping)"
    fi
  done
else
  echo "stage-assets: WARNING — no STR source at $STR_SRC; bundling no STR reference." >&2
fi

# cargo-packager requires the resource dirs to exist even if empty.
touch "$STAGE/.staged" "$MASK_STAGE/.staged" "$STR_STAGE/.staged"
