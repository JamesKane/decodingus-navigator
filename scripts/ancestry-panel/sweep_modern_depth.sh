#!/usr/bin/env bash
# Depth sweep — STEP 1 (heavy): build a WIDE (default 200k) Fst-ranked modern site set and the
# 1000G genotype matrix cut to it, so the fine-admixture + PCA assets can be rebuilt at several
# site-count caps (20k / 100k / 200k) WITHOUT re-slicing the 12 GB BCF per cap.
#
# Self-contained: reads only cached stage-1/2/3 intermediates already under $WORK (the restricted
# 1000G AF VCFs, the local 3202 BCF, the lifted 1240k BEDs) and writes everything under $WORK/sweep,
# so the shipping build state in $TMP is untouched. Modern-source only (1000G): the British fine
# structure the sweep validates lives entirely in the 1000G 26-pop set; AADR/SGDP feed the
# ancient/qpAdm assets, which this sweep leaves alone. Run STEP 2 with build_candidates.sh afterward.
#
# Override the wide cap / Fst floor via env: WIDE_MAX_SITES=200000 WIDE_MIN_FST=0.03 ./sweep_modern_depth.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool bcftools
require_tool tabix
require_tool cargo

WIDE_MAX_SITES="${WIDE_MAX_SITES:-200000}"
WIDE_MIN_FST="${WIDE_MIN_FST:-0.03}"

SWEEP="$WORK/sweep"; mkdir -p "$SWEEP"
FILTERED="$TMP/1kgp-chm13-1240k"   # stage-3 output: 1000G AF VCFs restricted to 1240k∩CHM13
[[ -d "$FILTERED" ]] && ls "$FILTERED"/*.vcf.gz >/dev/null 2>&1 \
  || die "missing restricted 1000G AF VCFs in $FILTERED (run 03_select_panel.sh once)"

# Local mirror of the phased biallelic 3202 genotype BCF (its .csi sits next to it).
KGP_BCF="${KGP_BCF:-$RAW/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf.gz}"
[[ -s "$KGP_BCF" ]] || die "missing local 1000G genotype BCF $KGP_BCF"

WIDE_PANEL="$SWEEP/wide_panel.bin"
WIDE_SITES_TSV="$SWEEP/wide_sites.tsv"        # contig pos ref alt fst per-pop-AF... (Fst column = the ranking)
WIDE_REGIONS="$SWEEP/wide_regions.tsv"        # sorted-unique CHROM<TAB>POS
KGP_WIDE_VCF="$SWEEP/1kgp.wide.vcf.gz"
KGP_WIDE_MATRIX="$SWEEP/1kgp.wide.matrix.tsv.gz"
KGP_WIDE_SAMPLES="$SWEEP/1kgp.wide.samples.txt"

# ── (A) select the WIDE Fst-ranked panel from the cached restricted 1000G AF VCFs ───────────────
if [[ -s "$WIDE_SITES_TSV" ]]; then
  log "have $(basename "$WIDE_SITES_TSV") (skip panel selection)"
else
  log "panelbuild panel WIDE (max_sites=$WIDE_MAX_SITES min_fst=$WIDE_MIN_FST) -> $WIDE_PANEL"
  cargo run --release -q -p navigator-panelbuild -- panel \
    --vcf-dir "$FILTERED" --out "$WIDE_PANEL" \
    --max-sites "$WIDE_MAX_SITES" --min-fst "$WIDE_MIN_FST" \
    --sites-tsv "$WIDE_SITES_TSV"
fi

# ── (B) wide regions (CHROM POS), sorted-unique — the target list for the BCF slice ─────────────
awk 'NR>1 { printf "%s\t%s\n", $1, $2 }' "$WIDE_SITES_TSV" | sort -k1,1 -k2,2n -u > "$WIDE_REGIONS"
log "wide site universe: $(wc -l < "$WIDE_REGIONS" | tr -d ' ') sites"

# ── (C) 1000G genotype matrix cut to the wide regions (slice local 12 GB BCF once) ──────────────
slice_at "$KGP_BCF" "$WIDE_REGIONS" "$KGP_WIDE_VCF" || die "1000G BCF slice failed"
if [[ -s "$KGP_WIDE_MATRIX" ]]; then
  log "have $(basename "$KGP_WIDE_MATRIX") (skip matrix)"
else
  log "matrix_from_vcf -> $(basename "$KGP_WIDE_MATRIX")"
  matrix_from_vcf "$KGP_WIDE_VCF" "$WIDE_REGIONS" "$KGP_WIDE_MATRIX" "$KGP_WIDE_SAMPLES"
fi

# Pop map for the modern (1000G 26-pop) fine set — reused by STEP 2.
cp -f "$RAW/1kgp.pops.tsv" "$SWEEP/pops.modern.tsv"

log "STEP 1 complete."
log "  wide sites : $WIDE_SITES_TSV ($(wc -l < "$WIDE_REGIONS" | tr -d ' ') sites, Fst-ranked)"
log "  1000G matrix: $KGP_WIDE_MATRIX ($(zcat "$KGP_WIDE_MATRIX" 2>/dev/null | wc -l | tr -d ' ') rows)"
log "  samples    : $KGP_WIDE_SAMPLES ($(wc -l < "$KGP_WIDE_SAMPLES" | tr -d ' ') samples)"
log "Next: ./build_candidates.sh   (STEP 2 — subset to caps + build pca/fine candidates)"
