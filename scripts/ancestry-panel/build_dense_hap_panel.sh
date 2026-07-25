#!/usr/bin/env bash
# Build a DENSE phased-haplotype reference (`ancestry_haps`) for the copying-LAI chromosome painter.
#
# WHY a separate builder. The shipped hap panel inherits its sites from the AIM selection built for
# the *frequency* panel: ~15.6k high-Fst sites genome-wide, ≈0.5 markers/Mb. A frequency estimator
# wants few, maximally-differentiated sites; a haplotype-COPYING model wants the opposite — many
# common sites, so a shared tract carries enough markers to identify whose haplotype it is. At the
# current density a 10 cM tract holds ~5 markers, and `navigator-panelbuild validate-lai` measures
# the consequence: held-out reference individuals are called their own sub-population 24.6% of the
# time (1000G) and *below chance* for the smaller HGDP populations. This builder selects sites for
# the copying model instead: common (MAF >= $MIN_MAF), biallelic SNVs, thinned to one per $SPACING bp.
#
# SITE SPACE. Sites are chosen from the 1000G CHM13-native phased BCF (no liftover, and it carries
# per-site MAF in INFO), then projected to hg38 to slice the HGDP statphase release, then lifted back
# so both sources land on the same CHM13 coordinates. The round-trip drops sites that don't lift
# cleanly in both directions — which is the intended filter: a site that can't be located in both
# references can't be shared by both panels anyway.
#
# COST. Genome-wide this streams a 13 GB BCF and a 6.3 GB VCF and runs two GATK liftovers; budget
# a couple of hours and ~20 GB of scratch. Pilot one chromosome first (CONTIGS=chr1) — accuracy per
# unit of genome is what the density question turns on, so a single chromosome answers it.
#
# USAGE
#   CONTIGS=chr1 SPACING=1000 ./build_dense_hap_panel.sh        # pilot
#   ./build_dense_hap_panel.sh                                   # genome-wide, SPACING=$SPACING
# Output: $OUT (default $TMP/ancestry_haps_dense_<tag>.bin) — it does NOT overwrite the live asset;
# validate it with `navigator-panelbuild validate-lai --haps <out>` before installing.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool bcftools tabix gatk curl

# ── parameters ──────────────────────────────────────────────────────────────────
CONTIGS="${CONTIGS:-$(seq 1 22 | sed 's/^/chr/' | paste -sd, -)}"   # autosomes
SPACING="${SPACING:-2000}"        # keep at most one site per this many bp
MIN_MAF="${MIN_MAF:-0.05}"        # global MAF floor (1000G INFO/MAF)
# Name the run after its contig scope, not a truncated contig list ("chr1chr2chr3chr4chr5chr6" for a
# genome-wide run is actively misleading when the outputs sit next to a real chr1-only pilot).
DEFAULT_TAG="$(echo "$CONTIGS" | tr ',' '\n' | wc -l | tr -d ' ')contigs"
[[ "$(echo "$CONTIGS" | tr ',' '\n' | wc -l | tr -d ' ')" == 22 ]] && DEFAULT_TAG=autosomes
[[ "$CONTIGS" != *,* ]] && DEFAULT_TAG="$CONTIGS"
TAG="${TAG:-${DEFAULT_TAG}_s${SPACING}}"
OUT="${OUT:-$TMP/ancestry_haps_dense_${TAG}.bin}"
WORKDIR="$TMP/dense_$TAG"; mkdir -p "$WORKDIR"

KGP_BCF="${KGP_BCF:-$RAW/$(basename "$KGP_GT_BCF_URL")}"
HGDP_VCF="${HGDP_VCF:-$RAW/hgdp.statphase.autosomes.vcf.gz}"
# sample<TAB>population for BOTH sources. The gnomAD HGDP+1KG metadata map covers all 4150 samples;
# $TMP/pops.<build>.tsv only carries the ~130 SGDP-derived HGDP rows, which would drop most of the
# HGDP panel on the floor.
POPMAP="${POPMAP:-$RAW/hgdp1kg.pops.tsv}"
# CHM13 -> hg38, the reverse of $CHAIN_GRCH38_TO_CHM13; needed to place CHM13-chosen sites in the
# hg38-coordinate HGDP release. Same UCSC hub directory as the forward chain.
CHAIN_CHM13_TO_HG38_URL="${CHAIN_CHM13_TO_HG38_URL:-https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/chm13v2-hg38.over.chain.gz}"

for f in "$KGP_BCF" "$HGDP_VCF" "$POPMAP" "$RAW/chm13v2.0.fa" "$RAW/hg38.fa"; do
  [[ -s "$f" ]] || die "missing input: $f"
done
log "dense hap panel: contigs=$CONTIGS spacing=${SPACING}bp minMAF=$MIN_MAF -> $OUT"

# ── (1) candidate sites from the 1000G CHM13 BCF: common biallelic SNVs, thinned by spacing ─────
SITES_CHM13="$WORKDIR/sites.chm13.tsv"
if [[ ! -s "$SITES_CHM13" ]]; then
  log "scanning $(basename "$KGP_BCF") for MAF>=$MIN_MAF SNVs"
  : > "$SITES_CHM13.part"
  for c in ${CONTIGS//,/ }; do
    bcftools query -r "$c" -i "MAF>=$MIN_MAF && TYPE=\"snp\"" -f '%CHROM\t%POS\n' "$KGP_BCF" \
      | awk -v s="$SPACING" 'BEGIN{last=-1e9} $2-last >= s {print; last=$2}' >> "$SITES_CHM13.part"
    log "  $c: $(wc -l < "$SITES_CHM13.part" | tr -d ' ') sites so far"
  done
  mv "$SITES_CHM13.part" "$SITES_CHM13"
fi
log "candidate sites (CHM13): $(wc -l < "$SITES_CHM13" | tr -d ' ')"

# ── (2) project the sites into hg38 so the HGDP release can be sliced at them ────────────────────
# A sites-only VCF is the liftover unit; genotypes aren't needed to move coordinates.
CHAIN_REV="$RAW/$(basename "$CHAIN_CHM13_TO_HG38_URL" .gz)"
if [[ ! -s "$CHAIN_REV" ]]; then
  fetch "$CHAIN_CHM13_TO_HG38_URL" "$(basename "$CHAIN_CHM13_TO_HG38_URL")" || die "cannot fetch the CHM13->hg38 chain"
  gunzip -kf "$RAW/$(basename "$CHAIN_CHM13_TO_HG38_URL")"
fi
SITES_HG38="$WORKDIR/sites.hg38.tsv"
if [[ ! -s "$SITES_HG38" ]]; then
  log "projecting sites CHM13 -> hg38"
  bcftools view -G -T "$SITES_CHM13" -Oz -o "$WORKDIR/sites.chm13.vcf.gz" "$KGP_BCF"
  tabix -f -p vcf "$WORKDIR/sites.chm13.vcf.gz"
  gatk_lift "$WORKDIR/sites.chm13.vcf.gz" "$CHAIN_REV" "$RAW/hg38.fa" "$WORKDIR/sites.hg38.vcf.gz"
  bcftools query -f '%CHROM\t%POS\n' "$WORKDIR/sites.hg38.vcf.gz" | sort -k1,1V -k2,2n > "$SITES_HG38"
fi
log "sites placed in hg38: $(wc -l < "$SITES_HG38" | tr -d ' ')"

# ── (3) HGDP: slice at those hg38 sites, lift back to CHM13, emit the phased matrix ──────────────
if [[ ! -s "$WORKDIR/hgdp.matrix.tsv.gz" ]]; then
  log "slicing HGDP at the dense sites (this streams the whole release once)"
  bcftools view -T "$SITES_HG38" -v snps -m2 -M2 --threads 4 -Oz -o "$WORKDIR/hgdp.hg38.vcf.gz" "$HGDP_VCF"
  tabix -f -p vcf "$WORKDIR/hgdp.hg38.vcf.gz"
  liftover_vcf "$WORKDIR/hgdp.hg38.vcf.gz" grch38 "$WORKDIR/hgdp.chm13.vcf.gz"
  bcftools query -f '%CHROM\t%POS\n' "$WORKDIR/hgdp.chm13.vcf.gz" | sort -k1,1V -k2,2n > "$WORKDIR/sites.final.tsv"
  matrix_from_vcf "$WORKDIR/hgdp.chm13.vcf.gz" "$WORKDIR/sites.final.tsv" \
    "$WORKDIR/hgdp.matrix.tsv.gz" "$WORKDIR/hgdp.samples.txt"
fi
log "HGDP matrix: $(gzcat "$WORKDIR/hgdp.matrix.tsv.gz" | wc -l | tr -d ' ') sites × $(wc -l < "$WORKDIR/hgdp.samples.txt" | tr -d ' ') samples"

# ── (4) 1000G: the same final sites, straight out of the CHM13-native BCF ────────────────────────
if [[ ! -s "$WORKDIR/1kgp.matrix.tsv.gz" ]]; then
  log "extracting the 1000G matrix at the final sites"
  bcftools query -T "$WORKDIR/sites.final.tsv" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$KGP_BCF" \
    | gzip > "$WORKDIR/1kgp.matrix.tsv.gz.part" || die "bcftools query failed on $KGP_BCF"
  mv "$WORKDIR/1kgp.matrix.tsv.gz.part" "$WORKDIR/1kgp.matrix.tsv.gz"
  bcftools query -l "$KGP_BCF" > "$WORKDIR/1kgp.samples.txt"
fi
log "1000G matrix: $(gzcat "$WORKDIR/1kgp.matrix.tsv.gz" | wc -l | tr -d ' ') sites × $(wc -l < "$WORKDIR/1kgp.samples.txt" | tr -d ' ') samples"

# ── (5) union into the bit-packed asset ─────────────────────────────────────────────────────────
log "panelbuild hap-panel -> $OUT"
cargo run --release -q -p navigator-panelbuild --manifest-path "$HERE/../../Cargo.toml" -- hap-panel \
  --matrix "$WORKDIR/1kgp.matrix.tsv.gz,$WORKDIR/hgdp.matrix.tsv.gz" \
  --samples "$WORKDIR/1kgp.samples.txt,$WORKDIR/hgdp.samples.txt" \
  --pops "$POPMAP" --out "$OUT"

log "done: $OUT"
log "validate before installing:"
log "  navigator-panelbuild validate-lai --haps $OUT --contigs ${CONTIGS%%,*} --pops GBR,CEU,FIN,TSI,IBS,French,Sardinian,Orcadian --replicates 3"
