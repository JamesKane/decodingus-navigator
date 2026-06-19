#!/usr/bin/env bash
# Stage 6 — publish the released assets + manifest to the CDN origin.
#
# Uploads the panel/PCA/freq assets and the manifest under $CDN_REMOTE/$CDN_PREFIX. The CDN
# (e.g. CloudFront in front of the S3 bucket) serves them to Navigator clients, which verify
# each download against the manifest's sha256. Versioned by $ASSET_VERSION so a new build
# doesn't clobber the one clients are pinned to. Dry-run by default; pass --apply to upload.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"

APPLY=0; [[ "${1:-}" == "--apply" ]] && APPLY=1
DEST="$CDN_REMOTE/$CDN_PREFIX/v$ASSET_VERSION"

# Pick the uploader: aws cli for s3://, rclone for any remote.
if [[ "$CDN_REMOTE" == s3://* ]]; then require_tool aws; CP=(aws s3 cp); LS_NOTE="aws s3";
else require_tool rclone; CP=(rclone copyto); LS_NOTE="rclone"; fi

for f in "$PANEL_OUT" "$PCA_OUT" "$FINE_OUT" "$MANIFEST"; do
  [[ -s "$f" ]] || die "missing asset $f (run 05_build_assets.sh)"
  log "verify $(basename "$f") sha256=$(sha256_of "$f")"
done

log "publish target: $DEST  ($LS_NOTE)$([[ $APPLY -eq 0 ]] && echo '  [dry-run]')"
for f in "$PANEL_OUT" "$PCA_OUT" "$FINE_OUT" "$MANIFEST"; do
  dst="$DEST/$(basename "$f")"
  if [[ $APPLY -eq 1 ]]; then
    "${CP[@]}" "$f" "$dst"
  else
    log "  would upload $(basename "$f") -> $dst"
  fi
done

# Update a stable "latest" pointer so clients can discover the newest version.
if [[ $APPLY -eq 1 ]]; then
  "${CP[@]}" "$MANIFEST" "$CDN_REMOTE/$CDN_PREFIX/latest.json"
  log "stage 6 complete: published v$ASSET_VERSION + refreshed latest.json"
else
  log "dry-run complete. Re-run with --apply to upload."
fi
