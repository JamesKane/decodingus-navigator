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
# packed/EIGENSTRAT .geno -> PACKEDPED (convertf) -> VCF (plink2) -> liftover hg19 -> CHM13.
# NOTE: recent AADR releases (v66.p1) ship the .geno as TGENO (transposed packed). Use AdmixTools'
# convertf, which reads TGENO; EIGENSOFT's convertf reads only GENO/EIGENSTRAT and aborts on TGENO.
# We just try convertf and skip gracefully on failure — so the build is convertf-implementation
# agnostic (works if convertf can read the format, builds modern-only if not).
AADR_PREFIX="$(ls "$RAW/"*"${AADR_DATASET}"*.geno 2>/dev/null | head -1 | sed 's/\.geno$//' || true)"
aadr_ready=0
if [[ -z "$AADR_PREFIX" ]]; then
  log "AADR genotypes not found — skipping ancient sources (download/unpack AADR first)."
elif [[ -s "$TMP/aadr.chm13.vcf.gz" ]]; then
  aadr_ready=1
elif ! { command -v convertf >/dev/null 2>&1 && command -v plink2 >/dev/null 2>&1; }; then
  log "WARN: convertf/plink2 not found — skipping ancient sources (install AdmixTools + plink2)."
else
  log "AADR packed/EIGENSTRAT -> VCF (convertf: $(command -v convertf))"
  cat > "$TMP/convertf.par" <<EOF
genotypename:    ${AADR_PREFIX}.geno
snpname:         ${AADR_PREFIX}.snp
indivname:       ${AADR_PREFIX}.ind
outputformat:    PACKEDPED
genotypeoutname: $TMP/aadr.bed
snpoutname:      $TMP/aadr.bim
indivoutname:    $TMP/aadr.fam
EOF
  # Map the panel's CHM13 sites back to AADR SNP ids (via the stage-2 lifted BED) so we recode +
  # lift only the ~20k panel SNPs, not all ~1.23M — converting/lifting the full 23k-sample matrix
  # would be tens of GB and hours of CrossMap, and the matrix is cut to the panel anyway.
  awk 'NR==FNR { p[$1"\t"$2]=1; next } { split($4,a,"|"); if (($1"\t"$3) in p) print a[1] }' \
    "$REGIONS" "$TMP/1240k_sites.${BUILD}.bed" > "$TMP/aadr_panel_ids.txt"
  log "AADR: $(wc -l < "$TMP/aadr_panel_ids.txt" | tr -d ' ') panel SNP ids to extract"
  # id-paste=iid: emit only the IID (the AADR genetic ID) as the VCF sample name. Default plink2
  # pastes FID_IID (here "1_Loschbour.AG"), which then won't match the .anno Genetic ID used by
  # derive_aadr_pops.sh — so the ancient labels would silently come out empty.
  if convertf -p "$TMP/convertf.par" \
     && plink2 --bfile "$TMP/aadr" --extract "$TMP/aadr_panel_ids.txt" \
               --export vcf bgz id-paste=iid --out "$TMP/aadr.hg19" --output-chr chrM; then
    liftover_vcf "$TMP/aadr.hg19.vcf.gz" hg19 "$TMP/aadr.chm13.vcf.gz"
    aadr_ready=1
  else
    log "WARN: AADR convert/recode failed — skipping ancient sources."
    [[ "$(head -c5 "${AADR_PREFIX}.geno" 2>/dev/null)" == "TGENO" ]] && \
      log "      (.geno is TGENO transposed packed — needs AdmixTools convertf; EIGENSOFT's can't read it)"
  fi
fi
if [[ "$aadr_ready" == 1 ]]; then
  emit_source aadr "$TMP/aadr.chm13.vcf.gz"
  # Ancient labels: AADR .anno "Group ID" -> deep component, via the curated map.
  log "deriving ancient sample->component map from .anno + $AADR_COMPONENT_MAP"
  "$HERE/derive_aadr_pops.sh" "$RAW" "$AADR_COMPONENT_MAP" "$TMP/aadr.samples.txt" >> "$POPMAP"
else
  log "building WITHOUT ancient deep components (1000G modern sources only)."
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

# ── SGDP (cteam_extended PLINK, GRCh37/hg19) — OPTIONAL ──────────────────────────
# NB: the cteam set is GRCh37/hg19 (verified by 1240k position overlap), NOT GRCh38, and is
# genome-wide (~34M SNPs). Restrict to the panel SNPs (as hg19 `chrom_pos` .bim ids — chrom without
# the 'chr' prefix to match the cteam .bim) and emit IID-only sample names before lifting hg19->CHM13.
if [[ "$SGDP_ENABLE" == 1 ]]; then
  merged="$TMP/sgdp.chm13.vcf.gz"
  if [[ -s "$merged" ]]; then
    emit_source sgdp "$merged"; add_popmap sgdp
  elif [[ -s "$RAW/${SGDP_PLINK_PREFIX}.bed" ]]; then
    require_tool plink2
    awk '
      FILENAME ~ /panel_regions/ { p[$1"\t"$2]=1; next }
      FILENAME ~ /[.]chm13.*[.]bed$/ { split($4,a,"|"); if (($1"\t"$3) in p) keep[a[1]]=1; next }
      FILENAME ~ /hg19[.]bed$/ { split($4,a,"|"); if (a[1] in keep) { c=$1; sub(/^chr/,"",c); print c"_"$3 } }
    ' "$REGIONS" "$TMP/1240k_sites.${BUILD}.bed" "$TMP/1240k_sites.hg19.bed" | sort -u > "$TMP/sgdp_panel_varids.txt"
    log "SGDP: $(wc -l < "$TMP/sgdp_panel_varids.txt" | tr -d ' ') panel SNP ids to extract"
    if plink2 --bfile "$RAW/$SGDP_PLINK_PREFIX" --extract "$TMP/sgdp_panel_varids.txt" \
              --export vcf bgz id-paste=iid --out "$TMP/sgdp.hg19" --output-chr chrM; then
      liftover_vcf "$TMP/sgdp.hg19.vcf.gz" hg19 "$merged"
      emit_source sgdp "$merged"; add_popmap sgdp
    else
      log "WARN: SGDP plink2 recode failed — skipping SGDP."
    fi
  else
    log "sgdp: PLINK $RAW/${SGDP_PLINK_PREFIX}.bed not found (run stage 01 with SGDP_ENABLE=1) — skip"
  fi
fi

# Record the comma-joined argument lists the asset build consumes.
# The ${arr[*]+…} guard keeps an empty array safe under `set -u` on bash 3.2 (macOS default).
( IFS=,; echo "${MATRICES[*]+"${MATRICES[*]}"}" ) > "$TMP/matrices.list"
( IFS=,; echo "${SAMPLES[*]+"${SAMPLES[*]}"}"  ) > "$TMP/samples.list"
log "stage 4 complete: $(wc -l < "$POPMAP") labelled samples; matrices in $TMP/matrices.list"
