#!/usr/bin/env bash
# Fetch + panel-slice + CHM13-lift the PHASED HGDP WGS (Bergström 2020 "statphase") for the
# copying-LAI haplotype reference. HGDP ships a statistically-phased release, so it enters the
# ancestry_haps panel directly (no phasing pass) alongside the phased 1000G source — its per-population
# depth (French/Sardinian/Basque/Russian/Orcadian…) gives the copying LAI the sub-continental European
# resolution 1000G alone (GBR/CEU/FIN/IBS/TSI) can't.
#
# Run AFTER 04_build_matrices.sh (needs panel_regions + the 1240k CHM13/hg38 site beds) and BEFORE
# 05_build_assets.sh (which picks up $TMP/hgdp.matrix.tsv.gz and appends HGDP into the hap panel).
# Produces: $TMP/hgdp.matrix.tsv.gz + $TMP/hgdp.samples.txt, and appends HGDP labels to $POPMAP.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool bcftools tabix gatk curl

REGIONS="$TMP/panel_regions.${BUILD}.tsv"
CHM13_BED="$TMP/1240k_sites.${BUILD}.bed"
HG38_BED="$TMP/1240k_sites.hg38.bed"
POPMAP="$TMP/pops.${BUILD}.tsv"
HGDP_POPS="${HGDP_POPS:-$RAW/hgdp1kg.pops.tsv}"   # sample<TAB>population (HGDP##### rows used)
VCF="$RAW/hgdp.statphase.autosomes.vcf.gz"
for f in "$REGIONS" "$CHM13_BED" "$HG38_BED"; do [[ -s "$f" ]] || die "missing $f — run 02/03/04 first"; done
[[ -s "$HGDP_POPS" ]] || die "missing HGDP pop map $HGDP_POPS (sample<TAB>population)"

# (1) Download the ~6 GB phased VCF + tabix index (once).
if [[ ! -s "$VCF" || ! -s "$VCF.tbi" ]]; then
  log "downloading HGDP statphase VCF -> $VCF"
  curl -sS --fail -o "$VCF" "$HGDP_STATPHASE_URL"
  curl -sS --fail -o "$VCF.tbi" "$HGDP_STATPHASE_URL.tbi"
fi

# (2) The panel sites in hg38 (autosomes): panel_regions(CHM13 pos) -> 1240k rsID -> hg38 pos, via the
#     lifted 1240k beds whose name field is `rsID|ref|alt`.
SITES_HG38="$TMP/hgdp_panel_sites.hg38.tsv"
awk '
  FNR==NR { key[$1"\t"$2]=1; next }
  FILENAME ~ /chm13/ { split($4,a,"|"); if (($1"\t"$3) in key) rs[a[1]]=1; next }
  FILENAME ~ /hg38/  { split($4,a,"|"); if ((a[1] in rs) && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print $1"\t"$3 }
' "$REGIONS" "$CHM13_BED" "$HG38_BED" | sort -k1,1V -k2,2n > "$SITES_HG38"
log "HGDP: $(wc -l < "$SITES_HG38" | tr -d ' ') autosomal panel sites in hg38"

# (3) Local slice (biallelic SNPs) -> liftover hg38->CHM13 -> phased matrix.
bcftools view -R "$SITES_HG38" -v snps -m2 -M2 --threads 4 -Oz -o "$TMP/hgdp.panel.hg38.vcf.gz" "$VCF"
tabix -f -p vcf "$TMP/hgdp.panel.hg38.vcf.gz"
liftover_vcf "$TMP/hgdp.panel.hg38.vcf.gz" grch38 "$TMP/hgdp.chm13.vcf.gz"
matrix_from_vcf "$TMP/hgdp.chm13.vcf.gz" "$REGIONS" "$TMP/hgdp.matrix.tsv.gz" "$TMP/hgdp.samples.txt"
log "HGDP matrix: $(gzcat "$TMP/hgdp.matrix.tsv.gz" | wc -l | tr -d ' ') sites, $(wc -l < "$TMP/hgdp.samples.txt" | tr -d ' ') samples"

# (4) Append HGDP sample->population labels to the unified pop map (idempotent-ish: 05 dedups by sample).
grep '^HGDP' "$HGDP_POPS" >> "$POPMAP" || log "NOTE: no HGDP rows in $HGDP_POPS"
log "HGDP labels appended to $POPMAP — now run 05_build_assets.sh"
