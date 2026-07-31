#!/usr/bin/env bash
set -uo pipefail
S="$1"; SC="$2"; BIN="$3"
[ -s "$SC/runs/$S.json" ] && { echo "SKIP $S (done)"; exit 0; }
t0=$SECONDS
for c in chr21 chr22; do
  "$BIN" call --subject "$S" --contig $c --out /dev/null >/dev/null 2>>"$SC/runs/$S.log" || { echo "FAIL-call $S $c"; exit 1; }
done
"$BIN" archaic-segments --subject "$S" --json > "$SC/runs/$S.json" 2>>"$SC/runs/$S.log" || { echo "FAIL-seg $S"; exit 1; }
echo "OK $S $((SECONDS-t0))s"
