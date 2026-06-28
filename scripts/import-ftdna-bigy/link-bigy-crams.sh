#!/usr/bin/env bash
# EXAMPLE: bulk-link FTDNA Big Y CRAMs to the workspace subjects they belong to, matched by FTDNA
# kit number. Adapt the paths to your own Big Y export and workspace — nothing here is specific to
# a particular machine; pass --root (and --db if yours isn't the default location).
#
# This doubles as a worked example of `navigator ingest --external-id` — resolving an existing
# subject by a vendor id instead of a donor name, skipping unknown ids, and forcing the test type.
#
# Expected layout — one directory per FTDNA sample id, named by the kit number:
#
#     <ROOT>/<KIT>/Big_Y-700/CP086569.2/chrYM.cram      (preferred)
#     <ROOT>/<KIT>/Big_Y-500/CP086569.2/chrYM.cram      (fallback when there is no 700)
#
# We link each CRAM to the workspace subject carrying that <KIT> as an `external_id` (source FTDNA).
# Directories whose name is NOT a known FTDNA kit — e.g. anonymous Big Y UUID folders, or kits not
# loaded into this workspace — are skipped. When both Big_Y-700 and Big_Y-500 exist we take only the
# 700 (the 500 data is merged into it).
#
# Idempotent: a CRAM already recorded as an alignment is skipped (both here, by a DB pre-check, and
# as a backstop inside `navigator ingest`). Re-running only links what is new.
#
# Usage:
#   link-bigy-crams.sh --root DIR [--apply] [--db FILE] [--navigator BIN] [--project NAME]
#
# --root is required. Default is a DRY RUN: it reports what it would link without writing; pass
# --apply to ingest. (Portable to the bash 3.2 that ships with macOS — no associative arrays.)
set -euo pipefail

# The Big Y export root — required, no default (this is an example; point it at your own tree).
ROOT=""
# The Navigator workspace SQLite DB. Defaults to the standard per-user location the GUI/CLI use.
DB="${HOME}/.decodingus/navigator-rs.db"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAVIGATOR=""
PROJECT=""
APPLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply)        APPLY=1; shift ;;
    --root)         ROOT="$2"; shift 2 ;;
    --db)           DB="$2"; shift 2 ;;
    --navigator)    NAVIGATOR="$2"; shift 2 ;;
    --project)      PROJECT="$2"; shift 2 ;;
    -h|--help)      sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *)              echo "error: unknown argument: $1" >&2; exit 2 ;;
  esac
done

die() { echo "error: $*" >&2; exit 1; }

command -v sqlite3 >/dev/null || die "sqlite3 not found on PATH"
[[ -n "$ROOT" ]] || die "pass --root <DIR> (the Big Y export root, one directory per kit number)"
[[ -d "$ROOT" ]] || die "Big Y root not found: $ROOT (is the volume mounted?)"
[[ -f "$DB"   ]] || die "workspace DB not found: $DB"

# Locate the navigator binary: explicit flag, then a release/debug build in the repo, then PATH.
if [[ -z "$NAVIGATOR" ]]; then
  for cand in "$HERE/../../target/release/navigator" "$HERE/../../target/debug/navigator" "$(command -v navigator || true)"; do
    if [[ -n "$cand" && -x "$cand" ]]; then NAVIGATOR="$cand"; break; fi
  done
fi
[[ -n "$NAVIGATOR" && -x "$NAVIGATOR" ]] || die "navigator binary not found — build it (cargo build -p navigator-ui) or pass --navigator"

echo "root      : $ROOT"
echo "db        : $DB"
echo "navigator : $NAVIGATOR"
echo "mode      : $([[ $APPLY == 1 ]] && echo APPLY || echo 'DRY RUN (pass --apply to write)')"
[[ -n "$PROJECT" ]] && echo "project   : $PROJECT"
echo

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Known FTDNA kit numbers in the workspace, the CRAM paths already linked as alignments, and the
# on-disk directory names — all sorted so we can set-match with comm/grep (bash-3.2 friendly).
sqlite3 "$DB" "SELECT external_id FROM external_id WHERE source='FTDNA';" | sort -u > "$WORK/known_kits.txt"
sqlite3 "$DB" "SELECT bam_path FROM alignment WHERE bam_path IS NOT NULL;" | sort -u > "$WORK/linked.txt"
for d in "$ROOT"/*/; do basename "$d"; done | sort -u > "$WORK/dirs.txt"

# Matched kits = directory names that are also known FTDNA kits. Everything else is skipped.
comm -12 "$WORK/known_kits.txt" "$WORK/dirs.txt" > "$WORK/matched.txt"

n_dirs=$(wc -l < "$WORK/dirs.txt" | tr -d ' ')
n_known=$(wc -l < "$WORK/known_kits.txt" | tr -d ' ')
n_matched=$(wc -l < "$WORK/matched.txt" | tr -d ' ')
echo "directories on disk           : $n_dirs"
echo "known FTDNA kits in workspace : $n_known"
echo "matched (dir name == kit)     : $n_matched"
echo "unmatched dirs (skipped)      : $((n_dirs - n_matched))"
echo

no_cram=0 already=0 to_link=0 linked=0 failed=0

while IFS= read -r kit; do
  [[ -n "$kit" ]] || continue

  # Prefer Big_Y-700; fall back to Big_Y-500 only when there is no 700.
  cram=""
  for ver in "Big_Y-700" "Big_Y-500"; do
    c="$ROOT/$kit/$ver/CP086569.2/chrYM.cram"
    if [[ -f "$c" ]]; then cram="$c"; break; fi
  done
  if [[ -z "$cram" ]]; then
    echo "NO-CRAM   $kit (matched, but no Big_Y CRAM on disk)"
    no_cram=$((no_cram + 1))
    continue
  fi

  if grep -Fxq "$cram" "$WORK/linked.txt"; then
    echo "SKIP      $kit (already linked)"
    already=$((already + 1))
    continue
  fi

  to_link=$((to_link + 1))
  if [[ $APPLY == 0 ]]; then
    echo "WOULD     $kit -> $cram"
    continue
  fi

  echo "LINK      $kit -> $cram"
  # The Big_Y-700/500 layout names the test definitively; force it (the CRAM has no .bai for the
  # coverage-shape detector, so inference would fall back to WGS).
  args=(ingest --external-id "$kit" --id-source FTDNA --skip-unmatched --test-type "Big Y")
  [[ -n "$PROJECT" ]] && args+=(--project "$PROJECT")
  args+=("$cram")
  if "$NAVIGATOR" "${args[@]}"; then
    linked=$((linked + 1))
  else
    echo "FAIL      $kit ($cram)" >&2
    failed=$((failed + 1))
  fi
done < "$WORK/matched.txt"

echo
echo "── summary ──────────────────────────────"
echo "matched kits with a CRAM to link : $to_link"
echo "already linked (skipped)         : $already"
echo "matched but no CRAM on disk      : $no_cram"
echo "directories with no kit match    : $((n_dirs - n_matched))"
if [[ $APPLY == 1 ]]; then
  echo "linked this run                  : $linked"
  echo "failed                           : $failed"
else
  echo "(dry run — re-run with --apply to link the $to_link above)"
fi
[[ $failed -eq 0 ]]
