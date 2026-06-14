#!/usr/bin/env bash
# Helper: generate the modern-source sample->population maps the pipeline needs, from public
# metadata, so they're reproducible instead of hand-curated:
#
#   $RAW/1kgp.pops.tsv     1000 Genomes (2504 unrelated) -> 26 fine populations   [REQUIRED]
#   $RAW/sgdp.pops.tsv     SGDP -> Population_ID                                   [SGDP_ENABLE=1]
#   $RAW/hgdp1kg.pops.tsv  gnomAD HGDP+1KG -> fine population                      [HGDP_1KG_ENABLE=1]
#
# stage 04 (add_popmap) reads these to label each modern source's samples. The AADR *ancient*
# labels are produced separately by derive_aadr_pops.sh (driven from the curated component map).
#
# Columns are resolved by header NAME (not fixed index), so upstream reordering won't silently
# mis-map. Idempotent: existing maps are kept unless --force. URLs override via env.
# (URLs verified 2026-06-13.)
#
#   ./derive_modern_pops.sh [--force]
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool curl
ensure_dirs

FORCE=0; [[ "${1:-}" == "--force" ]] && FORCE=1

KGP_PANEL_URL="${KGP_PANEL_URL:-http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel}"
SGDP_META_URL="${SGDP_META_URL:-https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.samples.txt}"
HGDP_META_URL="${HGDP_META_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv}"

# 1-based index of the tab-column whose header == $2 (leading '#' stripped), else 0.
col_idx() {  # <file> <colname>
  LC_ALL=C awk -F'\t' -v want="$2" '
    NR==1 { for (i=1;i<=NF;i++){ h=$i; sub(/^#/,"",h); if (h==want){print i; exit} } print 0; exit }
  ' "$1"
}

# Build a sample<TAB>pop map by header-resolved columns. Skips blanks, "NA", and #-comment rows
# (so a leading "#..." header line is ignored). LC_ALL=C tolerates non-UTF8 bytes (SGDP place names).
build_map() {  # <src-file> <sample-col-name> <pop-col-name> <out>
  local src="$1" sname="$2" pname="$3" out="$4" sc pc
  sc="$(col_idx "$src" "$sname")"; pc="$(col_idx "$src" "$pname")"
  [[ "$sc" -gt 0 && "$pc" -gt 0 ]] || die "$(basename "$src"): missing column '$sname' or '$pname'"
  LC_ALL=C awk -F'\t' -v s="$sc" -v p="$pc" \
    'NR>1 && $0 !~ /^#/ && $s!="" && $p!="" && $p!="NA" { print $s"\t"$p }' "$src" > "$out"
  [[ -s "$out" ]] || die "$(basename "$out"): produced no rows (check the source format)"
  log "wrote $(basename "$out"): $(wc -l < "$out" | tr -d ' ') samples, $(cut -f2 "$out" | LC_ALL=C sort -u | wc -l | tr -d ' ') populations"
}

have() { [[ -s "$1" && "$FORCE" == 0 ]]; }

# ── 1000G: sample -> fine population ('pop' column) ──────────────────────────────
if have "$RAW/1kgp.pops.tsv"; then log "have 1kgp.pops.tsv (skip; --force to rebuild)"; else
  fetch "$KGP_PANEL_URL" "1kgp_call_samples.panel" || die "fetch 1000G panel failed"
  build_map "$RAW/1kgp_call_samples.panel" sample pop "$RAW/1kgp.pops.tsv"
  log "  note: 2504 unrelated samples (the right PCA basis); the 698 related samples in the BCF stay unlabelled."
fi

# ── SGDP: SGDP_ID -> Population_ID ───────────────────────────────────────────────
if have "$RAW/sgdp.pops.tsv"; then log "have sgdp.pops.tsv (skip)"; else
  fetch "$SGDP_META_URL" "sgdp_metadata.samples.txt" || die "fetch SGDP meta failed"
  build_map "$RAW/sgdp_metadata.samples.txt" SGDP_ID Population_ID "$RAW/sgdp.pops.tsv"
  log "  note: keyed on SGDP_ID — verify it matches the PLINK .fam IID after the SGDP fetch."
fi

# ── gnomAD HGDP+1KG: sample 's' -> hgdp_tgp_meta.Population ───────────────────────
if have "$RAW/hgdp1kg.pops.tsv"; then log "have hgdp1kg.pops.tsv (skip)"; else
  fetch "$HGDP_META_URL" "gnomad_hgdp_tgp_meta.tsv" || die "fetch gnomAD meta failed"
  build_map "$RAW/gnomad_hgdp_tgp_meta.tsv" s "hgdp_tgp_meta.Population" "$RAW/hgdp1kg.pops.tsv"
  log "  note: this callset re-includes the 1KGP samples (overlap with the 1000G source) — label only"
  log "        HGDP samples, or accept the double-count, when wiring it into the basis."
fi

log "modern pop maps ready in $RAW (1kgp required; sgdp/hgdp1kg consumed only when those sources are enabled)."
