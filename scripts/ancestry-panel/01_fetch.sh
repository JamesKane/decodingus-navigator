#!/usr/bin/env bash
# Stage 1 — retrieve every input dataset + the CHM13 reference and liftover chains.
# Idempotent: existing, non-empty downloads are skipped. URLs come from config.sh
# (several are marked # VERIFY — confirm the current release before a real run).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
ensure_dirs
require_tool curl

# ── reference + chains ──────────────────────────────────────────────────────────
fetch "$CHM13_FASTA_URL"            "chm13v2.0.fa.gz"
fetch "$CHAIN_HG19_TO_CHM13"        "$(basename "$CHAIN_HG19_TO_CHM13")"
fetch "$CHAIN_GRCH38_TO_CHM13"      "$(basename "$CHAIN_GRCH38_TO_CHM13")"
# hg19->hg38 chain only needed to project the 1240k universe into hg38 for GRCh38 source slicing.
if [[ "$HGDP_1KG_ENABLE" == 1 || "$SGDP_ENABLE" == 1 ]]; then
  fetch "$CHAIN_HG19_TO_HG38" "$(basename "$CHAIN_HG19_TO_HG38")"
fi

# Unpack the FASTA + chains CrossMap needs uncompressed.
[[ -s "$RAW/chm13v2.0.fa" ]] || { log "gunzip CHM13 FASTA"; gunzip -kf "$RAW/chm13v2.0.fa.gz"; }
for c in "$CHAIN_HG19_TO_CHM13" "$CHAIN_GRCH38_TO_CHM13" "$CHAIN_HG19_TO_HG38"; do
  b="$RAW/$(basename "$c")"; [[ -s "$b" && "$b" == *.gz ]] && gunzip -kf "$b" || true
done

# ── AADR (ancient deep components) ──────────────────────────────────────────────
# Four EIGENSTRAT files from Harvard Dataverse, fetched by numeric file id and saved under
# the $AADR_FILE_PREFIX stem (.geno/.snp/.ind + .anno annotation). ~7.3 GB for the 1240K set.
aadr_ok=1
for pair in "geno:$AADR_ID_GENO" "snp:$AADR_ID_SNP" "ind:$AADR_ID_IND" "anno:$AADR_ID_ANNO"; do
  ext="${pair%%:*}"; id="${pair##*:}"
  fetch "$AADR_DOWNLOAD_BASE/$id" "${AADR_FILE_PREFIX}.$ext" || { aadr_ok=0; break; }
done
[[ "$aadr_ok" == 1 ]] || \
  log "NOTE: AADR auto-fetch failed — download $AADR_DATASET $AADR_VERSION manually (files ${AADR_FILE_PREFIX}.{geno,snp,ind,anno}) from the Dataverse landing page ($AADR_DATAVERSE_DOI) into $RAW/"

# ── 1000G on CHM13 (modern backbone, native build) ──────────────────────────────
# (a) AF files for the panel (stage 03): per-chrom `withafinfo` VCFs (+ .tbi). ~9.9 GB total —
#     sites+INFO only, carrying the AC_<POP>_unrel/AN_<POP>_unrel the panel needs.
# (b) Genotypes for PCA (stage 04): the phased biallelic BCF is NOT downloaded here — stage 04
#     remote-slices it at the panel sites from $KGP_GT_BCF_URL.
if compgen -G "$KGP_CHM13_DIR/*.vcf.gz" >/dev/null; then
  log "1000G-CHM13 AF: local files present in $KGP_CHM13_DIR (skip)"
else
  log "1000G-CHM13 AF: fetching per-chromosome withafinfo VCFs (~9.9 GB) from $KGP_AF_BASE_URL"
  for chr in $(seq 1 22) X Y; do
    # shellcheck disable=SC2059
    f="$(printf "$KGP_AF_PATTERN" "$chr")"
    if fetch "$KGP_AF_BASE_URL/$f" "$f"; then
      fetch "$KGP_AF_BASE_URL/$f.tbi" "$f.tbi" || log "  no .tbi for chr${chr} (stage 03 will index)"
      mv -f "$RAW/$f" "$KGP_CHM13_DIR/"; [[ -s "$RAW/$f.tbi" ]] && mv -f "$RAW/$f.tbi" "$KGP_CHM13_DIR/"
    else
      log "  skip chr${chr} (adjust KGP_AF_PATTERN/KGP_AF_BASE_URL)"
    fi
  done
fi

# ── HGDP + 1KG (gnomAD, modern global) — OPTIONAL, sliced in stage 04 ────────────
# Not fetched whole (~3.6 TB). Stage 04 remote-slices each chromosome at the 1240k-in-hg38 sites.
if [[ "$HGDP_1KG_ENABLE" == 1 ]]; then
  log "HGDP+1KG: enabled — stage 04 will slice $HGDP_1KG_BASE_URL at the 1240k-in-hg38 sites."
  log "  (anonymous gs:// reads of the public gnomAD bucket work; HGDP_1KG_GCP_PROJECT only needed if it ever flips to requester-pays)"
else
  log "HGDP+1KG: disabled (set HGDP_1KG_ENABLE=1 to include)."
fi

# ── SGDP (modern deep diversity) — OPTIONAL, PLINK fetched whole (~3 GB) ─────────
if [[ "$SGDP_ENABLE" == 1 ]]; then
  require_tool unzip
  sgdp_ok=1
  fetch "$SGDP_BASE_URL/${SGDP_PLINK_PREFIX}.bed" "${SGDP_PLINK_PREFIX}.bed" || sgdp_ok=0
  fetch "$SGDP_BASE_URL/${SGDP_PLINK_PREFIX}.fam" "${SGDP_PLINK_PREFIX}.fam" || sgdp_ok=0
  # The .bim ships zipped at sharehost — fetch + unzip to the plain .bim plink2 expects.
  if [[ ! -s "$RAW/${SGDP_PLINK_PREFIX}.bim" ]]; then
    if fetch "$SGDP_BASE_URL/${SGDP_PLINK_PREFIX}.bim.zip" "${SGDP_PLINK_PREFIX}.bim.zip"; then
      ( cd "$RAW" && unzip -o "${SGDP_PLINK_PREFIX}.bim.zip" >/dev/null ) || sgdp_ok=0
    else sgdp_ok=0; fi
  fi
  [[ "$sgdp_ok" == 1 ]] || log "NOTE: SGDP fetch failed — VERIFY SGDP_BASE_URL/SGDP_PLINK_PREFIX (host paths roll forward; .bim is zipped)."
else
  log "SGDP: disabled (set SGDP_ENABLE=1 to include)."
fi

# ── IBD assets: recombination map (genetic map) + optional GRCh38 FASTA ──────────
# Genetic map (GRCh38 PLINK recombination map) for the IBD genetic-map asset (stage 05 lifts it).
if [[ -d "$RAW/gmap_grch38" ]]; then
  log "genetic map: present (skip)"
elif fetch "$GMAP_URL" "plink.GRCh38.map.zip"; then
  require_tool unzip
  ( cd "$RAW" && unzip -o -q plink.GRCh38.map.zip -d gmap_grch38 ) || log "NOTE: genetic-map unzip failed"
else
  log "NOTE: genetic-map fetch failed — IBD genetic map will fall back to uniform 1 cM/Mb."
fi

# GRCh38 FASTA (+ .fai/.dict) and the hg19->hg38 chain, so the IBD panel can carry a GRCh38 allele
# column (lets GRCh38 consumer chips resolve). ~1 GB download; set IBD_GRCH38=0 to skip.
if [[ "$IBD_GRCH38" == 1 ]]; then
  fetch "$CHAIN_HG19_TO_HG38" "$(basename "$CHAIN_HG19_TO_HG38")" || log "NOTE: hg19->hg38 chain fetch failed"
  hgc="$RAW/$(basename "$CHAIN_HG19_TO_HG38")"; [[ -s "$hgc" && "$hgc" == *.gz ]] && gunzip -kf "$hgc" || true
  if fetch "$GRCH38_FASTA_URL" "$(basename "$GRCH38_FASTA_URL")"; then
    require_tool samtools
    hg38gz="$RAW/$(basename "$GRCH38_FASTA_URL")"; hg38fa="${hg38gz%.gz}"
    [[ -s "$hg38fa" ]] || { log "gunzip GRCh38 FASTA"; gunzip -kf "$hg38gz"; }
    [[ -s "$hg38fa.fai" ]] || samtools faidx "$hg38fa"
    [[ -s "${hg38fa%.fa}.dict" ]] || samtools dict "$hg38fa" -o "${hg38fa%.fa}.dict"
  else
    log "NOTE: GRCh38 FASTA fetch failed — IBD panel grch38 column will be omitted."
  fi
else
  log "IBD GRCh38 column disabled (IBD_GRCH38=0) — IBD panel will be CHM13 + GRCh37 only."
fi

log "stage 1 complete. Inputs under $RAW (+ 1000G AF in $KGP_CHM13_DIR)."
