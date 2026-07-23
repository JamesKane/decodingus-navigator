#!/usr/bin/env bash
# Extract present-day CONTINENTAL-European reference samples from the local AADR 1240k PLINK
# (aadr.bed, already converted in the June build — no convertf needed) at the WIDE sweep sites,
# lift hg19->CHM13, and emit a genotype matrix + pop map keyed by the new fine-population CODES.
# These fill the 1000G fine panel's gap (it has no French/German/Dutch/Swiss reference), so a
# continental European like the validation subject stops smearing into CEU + spurious Iberian/Tuscan.
#
# Present-day source groups (AADR Group ID -> new fine code), all high-coverage diploid (.DG/.SG):
#   French->FRN  Orcadian->ORC  Sardinian->SRD  Basque->BSQ  Italian_North->ITN  Russian->RUS
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool plink2
require_tool bcftools
require_tool tabix

SWEEP="$WORK/sweep"
WIDE_REGIONS="$SWEEP/wide_regions.tsv"
ANNO="$RAW/v66.p1.1240K.aadr.PUB.anno"
AADR_BFILE="$TMP/aadr"                      # aadr.bed/.bim/.fam from the June convertf run
BED_CHM13="$TMP/1240k_sites.${BUILD}.bed"   # col4 = rsid|ref|alt  (maps CHM13 locus -> AADR rsid)
for f in "$WIDE_REGIONS" "$ANNO" "$AADR_BFILE.bed" "$AADR_BFILE.fam" "$BED_CHM13"; do
  [[ -s "$f" ]] || die "missing $f"
done

# AADR present-day Group ID -> fine code (edit the here-doc to add/drop continental pops).
IDS_CODES="$SWEEP/continental.ids_codes.tsv"   # geneticID<TAB>code
: > "$IDS_CODES"
while read -r grp code; do
  [[ -z "$grp" || "$grp" == \#* ]] && continue
  awk -F'\t' -v g="$grp" -v c="$code" 'NR>1 && $15==g && $1 ~ /\.(DG|SG)$/{print $1"\t"c}' "$ANNO" >> "$IDS_CODES"
  n=$(awk -F'\t' -v g="$grp" 'NR>1 && $15==g && $1 ~ /\.(DG|SG)$/{c++}END{print c+0}' "$ANNO")
  log "  $grp -> $code : $n present-day samples"
done <<'PAIRS'
French        FRN
Orcadian      ORC
Sardinian     SRD
Basque        BSQ
Italian_North ITN
Russian       RUS
PAIRS
log "continental samples total: $(wc -l < "$IDS_CODES" | tr -d ' ')"

# plink2 --keep file (FID<TAB>IID) by joining wanted IIDs against aadr.fam.
KEEP="$SWEEP/continental.keep"
awk 'NR==FNR{want[$1]=1; next} ($2 in want){print $1"\t"$2}' "$IDS_CODES" "$AADR_BFILE.fam" > "$KEEP"
log "keep rows (FID IID matched in fam): $(wc -l < "$KEEP" | tr -d ' ')"

# Wide-panel AADR SNP ids: CHM13 (contig,pos) in wide_regions -> rsid via the 1240k CHM13 BED.
SNPS="$SWEEP/continental.snps.txt"
awk 'NR==FNR{p[$1"\t"$2]=1; next} {split($4,a,"|"); if (($1"\t"$3) in p) print a[1]}' \
  "$WIDE_REGIONS" "$BED_CHM13" | sort -u > "$SNPS"
log "wide-panel AADR SNP ids: $(wc -l < "$SNPS" | tr -d ' ')"

# Recode the selected samples+SNPs to VCF (hg19), lift to CHM13, build the matrix.
HG19_VCF="$SWEEP/aadr_continental.hg19.vcf.gz"
CHM13_VCF="$SWEEP/aadr_continental.chm13.vcf.gz"
MATRIX="$SWEEP/aadr_continental.matrix.tsv.gz"
SAMPLES="$SWEEP/aadr_continental.samples.txt"
POPS="$SWEEP/pops.continental.tsv"

log "plink2 recode (keep continental, extract wide SNPs) -> $(basename "$HG19_VCF")"
plink2 --bfile "$AADR_BFILE" --keep "$KEEP" --extract "$SNPS" \
       --export vcf bgz id-paste=iid --output-chr chrM --out "${HG19_VCF%.vcf.gz}"
[[ -s "$HG19_VCF" ]] || die "plink2 recode produced no VCF"

liftover_vcf "$HG19_VCF" hg19 "$CHM13_VCF"
matrix_from_vcf "$CHM13_VCF" "$WIDE_REGIONS" "$MATRIX" "$SAMPLES"

# Pop map for the matrix: sample geneticID -> fine code (only samples that survived recode/lift).
awk 'NR==FNR{code[$1]=$2; next} ($1 in code){print $1"\t"code[$1]}' "$IDS_CODES" "$SAMPLES" > "$POPS"
log "continental matrix: $(gzip -dc "$MATRIX" | wc -l | tr -d ' ') sites; $(wc -l < "$SAMPLES" | tr -d ' ') samples; pop map $(wc -l < "$POPS" | tr -d ' ') rows"
log "done. matrix=$MATRIX samples=$SAMPLES pops=$POPS"
