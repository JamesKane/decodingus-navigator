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

# Unpack the FASTA + chains CrossMap needs uncompressed.
[[ -s "$RAW/chm13v2.0.fa" ]] || { log "gunzip CHM13 FASTA"; gunzip -kf "$RAW/chm13v2.0.fa.gz"; }
for c in "$CHAIN_HG19_TO_CHM13" "$CHAIN_GRCH38_TO_CHM13"; do
  b="$RAW/$(basename "$c")"; [[ "$b" == *.gz ]] && gunzip -kf "$b" || true
done

# ── AADR (ancient deep components) ──────────────────────────────────────────────
# The AADR ships as a packed-ancestrymap triple (.geno/.snp/.ind) + .anno annotation.
# VERIFY the exact archive name on the landing page; this assembles the conventional path.
AADR_TAR="${AADR_DATASET}_${AADR_VERSION}.tar"
fetch "$AADR_BASE_URL/$AADR_VERSION/$AADR_TAR" "$AADR_TAR" || \
  log "NOTE: AADR auto-fetch failed — download $AADR_DATASET $AADR_VERSION manually from the Reich Lab landing page into $RAW/"
[[ -s "$RAW/$AADR_TAR" ]] && { log "unpack AADR"; tar -xf "$RAW/$AADR_TAR" -C "$RAW"; } || true

# ── HGDP + 1KG (modern global) ──────────────────────────────────────────────────
# gnomAD HGDP+1KG dense subset, per-chromosome GRCh38 VCFs. VERIFY filenames.
log "HGDP+1KG: confirm subset VCF names at $HGDP_1KG_BASE_URL and add fetch lines (per-chrom)."

# ── SGDP (modern deep diversity) ────────────────────────────────────────────────
log "SGDP: confirm VCF/eigenstrat names at $SGDP_BASE_URL and add fetch lines."

# ── 1000G on CHM13 (modern backbone, native build) ──────────────────────────────
if compgen -G "$KGP_CHM13_DIR/*.vcf.gz" >/dev/null; then
  log "1000G-CHM13: local mirror present in $KGP_CHM13_DIR (skip)"
else
  log "1000G-CHM13: fetching per-chromosome VCFs from $KGP_CHM13_BASE_URL (VERIFY names)"
  for chr in $(seq 1 22) X; do
    f="1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.withafinfo.vcf.gz"  # VERIFY
    fetch "$KGP_CHM13_BASE_URL/$f" "$f" && mv -f "$RAW/$f" "$KGP_CHM13_DIR/" || \
      log "  skip chr${chr} (adjust filename pattern in config/this script)"
  done
fi

log "stage 1 complete. Inputs under $RAW (+ 1000G in $KGP_CHM13_DIR)."
