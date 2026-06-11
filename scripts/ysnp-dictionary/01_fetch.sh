#!/usr/bin/env bash
# Stage 1 — download the YBrowse master Y-SNP CSVs, the GRCh38→CHM13 chain, and the CHM13
# FASTA (CrossMap's target reference). All land under $RAW. Resumable.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool curl
require_tool gunzip
ensure_dirs

fetch "$YBROWSE_HG38_URL" "snps_hg38.csv" || die "YBrowse hg38 CSV download failed"
# hg19 is optional — only needed to populate the GRCh37 coordinate. Don't fail the run on it.
fetch "$YBROWSE_HG19_URL" "snps_hg19.csv" || log "WARN: YBrowse hg19 CSV unavailable — GRCh37 coords will be omitted"

fetch "$CHAIN_GRCH38_TO_CHM13" || die "GRCh38→CHM13 chain download failed"

# The CHM13 FASTA is only needed if you later add reference-allele validation (a `CrossMap vcf`
# / bcftools-norm style pass). The coordinate liftover in 02_build uses `CrossMap bed`, which
# needs only the chain — so the ~1 GB FASTA is OPTIONAL. Fetch it only when FETCH_FASTA=1.
if [[ "${FETCH_FASTA:-0}" == "1" ]]; then
  fetch "$CHM13_FASTA_URL" || die "CHM13 FASTA download failed"
  FA_GZ="$RAW/$(basename "$CHM13_FASTA_URL")"; FA="$RAW/chm13v2.0.fa"
  if [[ ! -s "$FA" ]]; then log "unpacking $(basename "$FA_GZ")"; gunzip -c "$FA_GZ" > "$FA" || die "gunzip CHM13 FASTA failed"; fi
  if command -v samtools >/dev/null 2>&1 && [[ ! -s "$FA.fai" ]]; then
    log "indexing CHM13 FASTA"; samtools faidx "$FA" || log "WARN: samtools faidx failed"
  fi
fi

log "stage 1 complete -> $RAW"
