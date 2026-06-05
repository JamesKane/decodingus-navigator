#!/usr/bin/env bash
# Stage 4 — assemble per-source genotype matrices at the panel sites, on CHM13.
#
# For each reference source: convert to VCF (AADR is EIGENSTRAT), liftover to CHM13 +
# align alleles to the CHM13 reference, cut to the panel sites, and emit a panelbuild
# matrix (CHROM POS REF ALT GT...) + parallel sample list. Also derive the unified
# sample->population map (modern fine pops + ancient deep components).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool bcftools
require_tool CrossMap "pip install CrossMap"
require_tool tabix

REGIONS="$TMP/panel_regions.${BUILD}.tsv"
[[ -s "$REGIONS" ]] || die "missing $REGIONS (run 03_select_panel.sh)"
MATRICES=(); SAMPLES=(); POPMAP="$TMP/pops.${BUILD}.tsv"; : > "$POPMAP"

emit_source() {  # <tag> <vcf-on-CHM13.gz>
  local tag="$1" vcf="$2"
  local mat="$TMP/${tag}.matrix.tsv.gz" samp="$TMP/${tag}.samples.txt"
  log "matrix: $tag"
  matrix_from_vcf "$vcf" "$REGIONS" "$mat" "$samp"
  MATRICES+=("$mat"); SAMPLES+=("$samp")
}

# ── AADR (ancient deep components) ──────────────────────────────────────────────
# EIGENSTRAT/packed -> VCF (convertf to PED, then plink2 to VCF), liftover hg19 -> CHM13.
AADR_PREFIX="$(ls "$RAW/${AADR_DATASET}"*.geno 2>/dev/null | head -1 | sed 's/\.geno$//' || true)"
if [[ -n "$AADR_PREFIX" ]]; then
  if [[ ! -s "$TMP/aadr.chm13.vcf.gz" ]]; then
    require_tool convertf "EIGENSOFT/ADMIXTOOLS"
    require_tool plink2
    log "AADR EIGENSTRAT -> VCF"
    cat > "$TMP/convertf.par" <<EOF
genotypename:    ${AADR_PREFIX}.geno
snpname:         ${AADR_PREFIX}.snp
indivname:       ${AADR_PREFIX}.ind
outputformat:    PACKEDPED
genotypeoutname: $TMP/aadr.bed
snpoutname:      $TMP/aadr.bim
indivoutname:    $TMP/aadr.fam
EOF
    convertf -p "$TMP/convertf.par"
    plink2 --bfile "$TMP/aadr" --recode vcf bgz --out "$TMP/aadr.hg19" --output-chr chrM
    liftover_vcf "$TMP/aadr.hg19.vcf.gz" hg19 "$TMP/aadr.chm13.vcf.gz"
  fi
  emit_source aadr "$TMP/aadr.chm13.vcf.gz"
  # Ancient labels: AADR .anno "Group ID" -> deep component, via the curated map.
  log "deriving ancient sample->component map from .anno + $AADR_COMPONENT_MAP"
  "$HERE/derive_aadr_pops.sh" "$RAW" "$AADR_COMPONENT_MAP" "$TMP/aadr.samples.txt" >> "$POPMAP"
else
  log "AADR genotypes not found — skipping ancient sources (download/unpack AADR first)."
fi

# ── HGDP+1KG and SGDP (modern global) ───────────────────────────────────────────
# Per-chromosome GRCh38 VCFs: liftover -> CHM13, then matrix. Wire the actual fetched
# files here once 01_fetch.sh's HGDP/SGDP lines are filled in (filenames # VERIFY).
for src_tag in hgdp1kg sgdp; do
  if compgen -G "$RAW/${src_tag}*.vcf.gz" >/dev/null; then
    merged="$TMP/${src_tag}.chm13.vcf.gz"
    if [[ ! -s "$merged" ]]; then
      parts=(); for v in "$RAW/${src_tag}"*.vcf.gz; do
        o="$TMP/${src_tag}.$(basename "$v" .vcf.gz).chm13.vcf.gz"
        liftover_vcf "$v" grch38 "$o"; parts+=("$o")
      done
      bcftools concat -Oz -o "$merged" "${parts[@]}" && tabix -f -p vcf "$merged"
    fi
    emit_source "$src_tag" "$merged"
    # Modern pops: map sample -> reported population label (provide $RAW/${src_tag}.pops.tsv).
    [[ -s "$RAW/${src_tag}.pops.tsv" ]] && cat "$RAW/${src_tag}.pops.tsv" >> "$POPMAP" \
      || log "NOTE: provide $RAW/${src_tag}.pops.tsv (sample<TAB>population) for $src_tag"
  else
    log "$src_tag: no VCFs in $RAW (skip — add fetch lines in 01_fetch.sh)"
  fi
done

# ── 1000G on CHM13 (modern backbone, native build, no liftover) ─────────────────
if compgen -G "$KGP_CHM13_DIR/*.vcf.gz" >/dev/null; then
  merged="$TMP/1kgp.chm13.vcf.gz"
  [[ -s "$merged" ]] || { bcftools concat -Oz -o "$merged" "$KGP_CHM13_DIR"/*.vcf.gz && tabix -f -p vcf "$merged"; }
  emit_source 1kgp "$merged"
  [[ -s "$RAW/1kgp.pops.tsv" ]] && cat "$RAW/1kgp.pops.tsv" >> "$POPMAP" \
    || log "NOTE: provide $RAW/1kgp.pops.tsv (sample<TAB>fine-population) for 1000G"
fi

# Record the comma-joined argument lists the asset build consumes.
( IFS=,; echo "${MATRICES[*]}" ) > "$TMP/matrices.list"
( IFS=,; echo "${SAMPLES[*]}"  ) > "$TMP/samples.list"
log "stage 4 complete: $(wc -l < "$POPMAP") labelled samples; matrices in $TMP/matrices.list"
