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
  [[ -n "$HGDP_1KG_GCP_PROJECT" ]] \
    || log "WARN: HGDP_1KG_ENABLE=1 but HGDP_1KG_GCP_PROJECT is unset — requester-pays slicing in stage 04 will fail."
  log "HGDP+1KG: enabled — stage 04 will slice $HGDP_1KG_BASE_URL at panel sites (no bulk download)."
else
  log "HGDP+1KG: disabled (set HGDP_1KG_ENABLE=1 + HGDP_1KG_GCP_PROJECT to include)."
fi

# ── SGDP (modern deep diversity) — OPTIONAL, PLINK fetched whole (~3 GB) ─────────
if [[ "$SGDP_ENABLE" == 1 ]]; then
  sgdp_ok=1
  for ext in bed bim fam; do
    fetch "$SGDP_BASE_URL/${SGDP_PLINK_PREFIX}.$ext" "${SGDP_PLINK_PREFIX}.$ext" || { sgdp_ok=0; break; }
  done
  [[ "$sgdp_ok" == 1 ]] || log "NOTE: SGDP fetch failed — VERIFY SGDP_PLINK_PREFIX/SGDP_BASE_URL (Reich host PLINK names roll forward)."
else
  log "SGDP: disabled (set SGDP_ENABLE=1 to include)."
fi

log "stage 1 complete. Inputs under $RAW (+ 1000G AF in $KGP_CHM13_DIR)."
