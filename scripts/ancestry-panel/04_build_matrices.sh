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
AADR_PREFIX="$(ls "$RAW/"*"${AADR_DATASET}"*.geno 2>/dev/null | head -1 | sed 's/\.geno$//' || true)"
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

# Append a source's sample->population map (or note that it's missing).
add_popmap() {  # <tag>
  local t="$1"
  [[ -s "$RAW/$t.pops.tsv" ]] && cat "$RAW/$t.pops.tsv" >> "$POPMAP" \
    || log "NOTE: provide $RAW/$t.pops.tsv (sample<TAB>population) for $t"
}

# ── 1000G on CHM13 (PCA-basis genotypes, native build, no liftover) ──────────────
# Slice the phased biallelic BCF at the panel sites — remote-streamed from $KGP_GT_BCF_URL
# (point it at a local mirror to avoid streaming). The AF files in $KGP_CHM13_DIR are sites-only
# (stage 03's panel) and carry no genotypes, so they are NOT used here.
KGP_GT="$TMP/1kgp.chm13.vcf.gz"
slice_at "$KGP_GT_BCF_URL" "$REGIONS" "$KGP_GT" \
  || log "1000G genotype slice failed (point KGP_GT_BCF_URL at a local mirror?)"
if [[ -s "$KGP_GT" ]]; then emit_source 1kgp "$KGP_GT"; add_popmap 1kgp; fi

# ── HGDP+1KG (gnomAD, GRCh38) — OPTIONAL, remote-sliced at 1240k-in-hg38 ─────────
if [[ "$HGDP_1KG_ENABLE" == 1 ]]; then
  SITES_HG38="$TMP/1240k_sites.hg38.tsv"; merged="$TMP/hgdp1kg.chm13.vcf.gz"
  if [[ -s "$merged" ]]; then
    emit_source hgdp1kg "$merged"; add_popmap hgdp1kg
  elif [[ -s "$SITES_HG38" ]]; then
    export GCS_REQUESTER_PAYS_PROJECT="$HGDP_1KG_GCP_PROJECT"
    parts=()
    for chr in $(seq 1 22) X; do
      # shellcheck disable=SC2059
      src="$HGDP_1KG_BASE_URL/$(printf "$HGDP_1KG_PATTERN" "$chr")"
      sl="$TMP/hgdp1kg.chr${chr}.hg38.vcf.gz"
      slice_at "$src" "$SITES_HG38" "$sl" || { log "  hgdp1kg chr${chr} slice failed (requester-pays project / gs access?)"; continue; }
      o="$TMP/hgdp1kg.chr${chr}.chm13.vcf.gz"; liftover_vcf "$sl" grch38 "$o"; parts+=("$o")
    done
    if ((${#parts[@]})); then
      bcftools concat -Oz -o "$merged" "${parts[@]}" && tabix -f -p vcf "$merged"
      emit_source hgdp1kg "$merged"; add_popmap hgdp1kg
    else
      log "hgdp1kg: no chromosomes sliced — skip"
    fi
  else
    log "hgdp1kg: missing $SITES_HG38 (run stage 02 with the hg19->hg38 chain) — skip"
  fi
fi

# ── SGDP (GRCh38 PLINK) — OPTIONAL ──────────────────────────────────────────────
if [[ "$SGDP_ENABLE" == 1 ]]; then
  merged="$TMP/sgdp.chm13.vcf.gz"
  if [[ -s "$merged" ]]; then
    emit_source sgdp "$merged"; add_popmap sgdp
  elif [[ -s "$RAW/${SGDP_PLINK_PREFIX}.bed" ]]; then
    require_tool plink2
    log "SGDP PLINK -> VCF"
    plink2 --bfile "$RAW/$SGDP_PLINK_PREFIX" --recode vcf bgz --out "$TMP/sgdp.hg38" --output-chr chrM
    liftover_vcf "$TMP/sgdp.hg38.vcf.gz" grch38 "$merged"
    emit_source sgdp "$merged"; add_popmap sgdp
  else
    log "sgdp: PLINK $RAW/${SGDP_PLINK_PREFIX}.bed not found (run stage 01 with SGDP_ENABLE=1) — skip"
  fi
fi

# Record the comma-joined argument lists the asset build consumes.
( IFS=,; echo "${MATRICES[*]}" ) > "$TMP/matrices.list"
( IFS=,; echo "${SAMPLES[*]}"  ) > "$TMP/samples.list"
log "stage 4 complete: $(wc -l < "$POPMAP") labelled samples; matrices in $TMP/matrices.list"
