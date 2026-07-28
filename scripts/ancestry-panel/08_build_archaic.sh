#!/usr/bin/env bash
# Stage 8 — archaic (Neanderthal / Denisovan) marker panel.
#
# Design: documents/design/ArchaicAncestry_Design.md §4 (assets) and §10 (M1).
#
# Produces $ASSETS/archaic_markers_${BUILD}.bin — the Tier A asset the 23andMe-equivalent marker
# count reads. Independent of stages 03-06 (it does not use the 1240k AIM universe); it needs only
# stage 01's 1000G-on-CHM13 allele-frequency VCFs, for the African-outgroup filter.
#
# LICENSING — read before changing the fetch logic. The EVA archaic genomes are governed by the
# Ft. Lauderdale principles, not an explicit open licence (design §2). We therefore fetch them here
# at BUILD time and redistribute ONLY our derived sites. Do not add the raw archaic VCFs, their
# FilterBed masks, or the Ensembl ancestral sequence to any bundled/published asset.
#
# The pipeline straddles the liftover boundary, because polarity must be assigned in GRCh37 (the
# build both the archaic VCFs and the EPO ancestral sequence use) while the asset ships in CHM13:
#
#   1  bcftools  per-genome biallelic-SNV genotype tables, masked by each genome's own FilterBed
#   2  panelbuild archaic-candidates   (GRCh37) polarize + require an archaic hom-derived call
#   3  CrossMap  GRCh37 -> CHM13, as stage 02 does
#   4  bcftools  African / non-African allele frequencies at the lifted positions
#   5  panelbuild archaic-panel        (CHM13) orient, filter, classify, write the asset
#
# Step 5's orientation against the CHM13 FASTA is mandatory: CrossMap is not allele-aware, so a
# large share of sites arrive ref/alt-swapped, and nothing downstream fails loudly if they are left
# that way (the ancient-ancestry §7.16 defect).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool bcftools "brew install bcftools"
require_tool CrossMap "pip install CrossMap"
require_tool tar
ensure_dirs

ANC_DIR="$RAW/ancestral_GRCh37"
mkdir -p "$ARCHAIC_RAW" "$ANC_DIR" "$TMP"

CHROMS="${ARCHAIC_CHROMS:-$(seq 1 22 | tr '\n' ' ')}"

# ── 0. inputs ───────────────────────────────────────────────────────────────────
# Ensembl release-75 EPO ancestral sequence (GRCh37) — the polarity resource. ~766 MB tarball
# expanding to one FASTA per chromosome, which is exactly what archaic-candidates --ancestral wants.
if [[ ! -s "$ANC_DIR/homo_sapiens_ancestor_1.fa" ]]; then
  fetch "$ANCESTRAL_GRCH37_URL" "homo_sapiens_ancestor_GRCh37_e71.tar.bz2" \
    || die "could not fetch the EPO ancestral sequence"
  log "extracting ancestral sequence -> $ANC_DIR"
  tar -xjf "$RAW/homo_sapiens_ancestor_GRCh37_e71.tar.bz2" -C "$ANC_DIR" --strip-components=1 \
    || die "ancestral extract failed"
fi

# Archaic genomes. Order matters: it is the order navigator_analysis::archaic::ARCHAIC_GENOMES
# fixes, and the per-site call arrays are indexed by it.
#   name | per-chromosome VCF URL template (%s = chromosome) | FilterBed URL template
ARCHAIC_NAMES=(AltaiNeandertal Vindija33.19 Chagyrskaya8 Denisova3)

archaic_vcf_url() {  # <name> <chrom>
  case "$1" in
    AltaiNeandertal) printf "$ALTAI_VCF_PATTERN" "$2";;
    Vindija33.19)    printf "$VINDIJA_VCF_PATTERN" "$2";;
    Chagyrskaya8)    printf "$CHAGYRSKAYA_VCF_PATTERN" "$2";;
    Denisova3)       printf "$DENISOVA_VCF_PATTERN" "$2";;
    *) die "unknown archaic genome: $1";;
  esac
}
# Is this file BGZF (block-gzip), i.e. capable of indexed random access?
#
# Not academic: EVA publishes **Chagyrskaya as plain gzip** while Altai / Vindija / Denisova are
# BGZF (verified 2026-07-27 from the gzip FEXTRA "BC" subfield). Plain gzip cannot be indexed at all
# — `bcftools index` reports "a format that cannot be usefully indexed" — yet the server ships a
# .tbi next to it anyway, which is an upstream inconsistency: that index is unusable and any -R
# query against the file fails outright. Pass B therefore selects its access mode per file.
is_bgzf() {
  od -A n -t x1 -N 16 "$1" 2>/dev/null | tr -d ' \n' | grep -q '4243'
}

archaic_mask_url() {  # <name> <chrom>
  case "$1" in
    AltaiNeandertal) printf "$ALTAI_MASK_PATTERN" "$2";;
    Vindija33.19)    printf "$VINDIJA_MASK_PATTERN" "$2";;
    Chagyrskaya8)    printf "$CHAGYRSKAYA_MASK_PATTERN" "$2";;
    Denisova3)       printf "$DENISOVA_MASK_PATTERN" "$2";;
    *) die "unknown archaic genome: $1";;
  esac
}

# ── 1. fetch, then genotype in TWO passes ───────────────────────────────────────
# The EVA archaic VCFs are ALL-SITES files: where a genome matches hg19 it emits a
# reference-confident record (ALT=".", GT=0/0) rather than nothing. That makes a single
# variants-only extraction wrong in a way that is silent and one-sided:
#
#   * `HomAncestral` would be indistinguishable from "masked out", and
#   * where the EPO ancestral sequence says the REFERENCE allele is the derived one, a hom-ref
#     archaic genome IS homozygous-derived — the donor state this panel selects on. Keeping only
#     ALT-bearing records would drop every such site, biasing the panel against exactly the sites
#     where hg19 itself carries the archaic allele.
#
# So: pass A discovers the candidate universe from variant records (cheap, small), and pass B
# re-genotypes every genome at that universe WITHOUT the variants-only filter, so reference-confident
# calls come through. Each genome is still masked by its own FilterBed.
for name in "${ARCHAIC_NAMES[@]}"; do
  for c in $CHROMS; do
    vcf="$ARCHAIC_RAW/${name}.chr${c}.vcf.gz"
    mask="$ARCHAIC_RAW/${name}.chr${c}.mask.bed.gz"
    [[ -s "$vcf"  ]] || fetch "$(archaic_vcf_url "$name" "$c")"  "${name}.chr${c}.vcf.gz"      "$ARCHAIC_RAW" || log "WARN: no VCF for $name chr$c"
    [[ -s "$vcf.tbi" ]] || fetch "$(archaic_vcf_url "$name" "$c").tbi" "${name}.chr${c}.vcf.gz.tbi" "$ARCHAIC_RAW" || true
    [[ -s "$mask" ]] || fetch "$(archaic_mask_url "$name" "$c")" "${name}.chr${c}.mask.bed.gz" "$ARCHAIC_RAW" || log "WARN: no mask for $name chr$c (unfiltered)"
    # Build a missing index only where one is possible — see is_bgzf below.
    if [[ -s "$vcf" ]] && is_bgzf "$vcf" && [[ ! -s "$vcf.tbi" && ! -s "$vcf.csi" ]]; then
      log "indexing $(basename "$vcf") (no index published)"
      bcftools index -f -t "$vcf" 2>/dev/null || log "WARN: indexing failed for $(basename "$vcf")"
    fi
  done
done

# Pass A — the candidate universe: positions where ANY archaic genome carries a biallelic SNV.
UNIVERSE="$TMP/archaic_universe.grch37.tsv"
if [[ ! -s "$UNIVERSE" ]]; then
  log "pass A: discovering the candidate universe (variant records only)"
  : > "$UNIVERSE.part"
  for name in "${ARCHAIC_NAMES[@]}"; do
    for c in $CHROMS; do
      vcf="$ARCHAIC_RAW/${name}.chr${c}.vcf.gz"; mask="$ARCHAIC_RAW/${name}.chr${c}.mask.bed.gz"
      [[ -s "$vcf" ]] || continue
      region_args=(); [[ -s "$mask" ]] && region_args=(-T "$mask")
      bcftools view -m2 -M2 -v snps ${region_args[@]+"${region_args[@]}"} -Ou "$vcf" \
        | bcftools query -f '%CHROM\t%POS\n' >> "$UNIVERSE.part" || log "WARN: pass A failed for $name chr$c"
    done
  done
  sort -k1,1 -k2,2n -u "$UNIVERSE.part" > "$UNIVERSE" && rm -f "$UNIVERSE.part"
  log "candidate universe: $(wc -l < "$UNIVERSE") positions"
fi

# Pass B — every genome's state at the universe, reference-confident records included.
for name in "${ARCHAIC_NAMES[@]}"; do
  out="$TMP/archaic.${name}.tsv.gz"
  [[ -s "$out" ]] && { log "have $(basename "$out") (skip)"; continue; }
  log "pass B: genotyping $name at the universe"
  : > "$TMP/.archaic.$name.part"
  for c in $CHROMS; do
    vcf="$ARCHAIC_RAW/${name}.chr${c}.vcf.gz"; mask="$ARCHAIC_RAW/${name}.chr${c}.mask.bed.gz"
    [[ -s "$vcf" ]] || continue
    # Access mode per file: -R index-jumps straight to the universe positions (fast, but needs BGZF
    # plus an index); -T streams the file and filters (no index, works on plain gzip). Chagyrskaya
    # is plain gzip, so it takes the streaming path.
    if is_bgzf "$vcf" && [[ -s "$vcf.tbi" || -s "$vcf.csi" ]]; then
      sel=(-R "$UNIVERSE")
    else
      sel=(-T "$UNIVERSE")
    fi
    # The genome's own quality mask is applied as a SECOND view in the pipe: bcftools takes at most
    # one -T per invocation, and a BED is accepted directly (expanding it to per-position rows would
    # be hundreds of millions of lines for a large chromosome).
    if [[ -s "$mask" ]]; then
      bcftools view "${sel[@]}" -Ou "$vcf" \
        | bcftools view -T "$mask" -Ou \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
        >> "$TMP/.archaic.$name.part" || log "WARN: pass B failed for $name chr$c"
    else
      bcftools view "${sel[@]}" -Ou "$vcf" \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
        >> "$TMP/.archaic.$name.part" || log "WARN: pass B failed for $name chr$c"
    fi
  done
  [[ -s "$TMP/.archaic.$name.part" ]] || die "no calls extracted for $name"
  gzip -c "$TMP/.archaic.$name.part" > "$out" && rm -f "$TMP/.archaic.$name.part"
  log "$name: $(gzip -dc "$out" | wc -l) rows"
done

# ── 2. candidates (GRCh37): polarity + archaic hom-derived ───────────────────────
CAND="$TMP/archaic_candidates.grch37.tsv"
CAND_BED="$TMP/archaic_candidates.grch37.bed"
if [[ ! -s "$CAND" ]]; then
  log "selecting candidates (polarity from EPO ancestral, GRCh37)"
  cargo run --release -q -p navigator-panelbuild -- archaic-candidates \
    --archaic "$TMP/archaic.AltaiNeandertal.tsv.gz,$TMP/archaic.Vindija33.19.tsv.gz,$TMP/archaic.Chagyrskaya8.tsv.gz,$TMP/archaic.Denisova3.tsv.gz" \
    --ancestral "$ANC_DIR" \
    --min-archaic-called "${ARCHAIC_MIN_CALLED:-1}" \
    --out "$CAND" --out-bed "$CAND_BED" || die "archaic-candidates failed"
fi

# ── 3. lift GRCh37 -> CHM13 (same mechanism as stage 02) ────────────────────────
LIFTED="$TMP/archaic_candidates.${BUILD}.bed"
if [[ ! -s "$LIFTED" ]]; then
  CHAIN="$(chain_for hg19)"
  [[ -s "$CHAIN" ]] || die "hg19->CHM13 chain missing ($CHAIN) — run 01_fetch.sh"
  log "CrossMap bed -> $BUILD"
  CrossMap bed "$CHAIN" "$CAND_BED" "$LIFTED" || die "CrossMap bed failed"
  log "lifted $(wc -l < "$LIFTED") of $(wc -l < "$CAND_BED") candidates"
fi

# ── 3b. also project the candidates into GRCh38 ─────────────────────────────────
# So the panel carries hg38 coordinates and a GRCh38 alignment can be genotyped without a runtime
# liftover. GRCh37 needs no lift at all — those are the archaic VCFs' own coordinates, carried
# straight through. The hg38 lift is NOT allele-aware, so archaic-panel orients it against the hg38
# reference exactly as it does for CHM13.
LIFTED38="$TMP/archaic_candidates.hg38.bed"
HG38_FA="${HG38_FASTA:-$RAW/hg38.fa}"
if [[ ! -s "$LIFTED38" ]]; then
  HG38_CHAIN="$RAW/$(basename "$CHAIN_HG19_TO_HG38" .gz)"
  if [[ -s "$HG38_CHAIN" ]]; then
    log "CrossMap bed -> hg38"
    CrossMap bed "$HG38_CHAIN" "$CAND_BED" "$LIFTED38" || log "WARN: hg38 lift failed — GRCh38 alignments will fall back to consensus coverage"
  else
    log "NOTE: no hg19->hg38 chain at $HG38_CHAIN — skipping GRCh38 loci"
  fi
fi

# ── 4. African-outgroup frequencies at the lifted positions (CHM13) ─────────────
# From the per-super-pop AC_<POP>_unrel / AN_<POP>_unrel INFO fields the stage-01 `withafinfo`
# VCFs already carry (design §9 Q2 — no new data source). Emitted as
# CHROM POS REF ALT AF_AFR AF_NONAFR, with AF stated for the VCF's own ALT; the builder
# re-expresses it against the derived base.
OG_AF="$TMP/archaic_outgroup_af.${BUILD}.tsv"
if [[ ! -s "$OG_AF" ]]; then
  awk '{ printf "%s\t%d\n", $1, $3 }' "$LIFTED" | sort -k1,1 -k2,2n -u > "$TMP/archaic_sites.${BUILD}.tsv"
  : > "$OG_AF.part"
  for c in $CHROMS; do
    # shellcheck disable=SC2059
    af_vcf="$KGP_CHM13_DIR/$(printf "$KGP_AF_PATTERN" "$c")"
    [[ -s "$af_vcf" ]] || { log "WARN: missing $(basename "$af_vcf") — run 01_fetch.sh"; continue; }
    bcftools query -R "$TMP/archaic_sites.${BUILD}.tsv" \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC_AFR_unrel\t%INFO/AN_AFR_unrel\t%INFO/AC_AMR_unrel\t%INFO/AN_AMR_unrel\t%INFO/AC_EAS_unrel\t%INFO/AN_EAS_unrel\t%INFO/AC_EUR_unrel\t%INFO/AN_EUR_unrel\t%INFO/AC_SAS_unrel\t%INFO/AN_SAS_unrel\n' \
      "$af_vcf" >> "$OG_AF.part" || log "WARN: AF query failed for chr$c"
  done
  # AFR frequency, and a pooled non-African frequency over AMR+EAS+EUR+SAS. Rows with a zero or
  # missing denominator are dropped rather than defaulted — a fabricated 0 would sail through the
  # "rare in Africa" filter and silently inflate the panel.
  awk -F'\t' 'BEGIN{OFS="\t"}
    { afr_ac=$5+0; afr_an=$6+0;
      non_ac=$7+$9+$11+$13; non_an=$8+$10+$12+$14;
      if (afr_an<=0 || non_an<=0) next;
      printf "%s\t%s\t%s\t%s\t%.6f\t%.6f\n", $1,$2,$3,$4, afr_ac/afr_an, non_ac/non_an }' \
    "$OG_AF.part" > "$OG_AF" && rm -f "$OG_AF.part"
  log "outgroup AF rows: $(wc -l < "$OG_AF")"
fi

# ── 5. the asset (CHM13): orient, filter, classify ──────────────────────────────
REF_FA="${CHM13_FASTA:-$RAW/chm13v2.0.fa}"
[[ -s "$REF_FA" ]] || die "CHM13 FASTA not found at $REF_FA (set CHM13_FASTA); orientation cannot be skipped"
# The orientation step uses an INDEXED FASTA reader, and a missing .fai surfaces as a confusing
# "No such file or directory" naming the .fa itself — so check for it explicitly.
[[ -s "$REF_FA.fai" ]] || { command -v samtools >/dev/null 2>&1 && samtools faidx "$REF_FA"; } \
  || die "no index at $REF_FA.fai (run: samtools faidx $REF_FA)"
OUT="$ASSETS/archaic_markers_${BUILD}.bin"
log "building $OUT"
cargo run --release -q -p navigator-panelbuild -- archaic-panel \
  --candidates "$CAND" \
  --lifted "$LIFTED" \
  --outgroup-af "$OG_AF" \
  --reference "$REF_FA" \
  --max-afr-freq "${ARCHAIC_MAX_AFR_FREQ:-0.01}" \
  --min-non-afr-freq "${ARCHAIC_MIN_NON_AFR_FREQ:-0.05}" \
  --sites-tsv "$TMP/archaic_markers.${BUILD}.tsv" \
  ${LIFTED38:+$([[ -s "$LIFTED38" && -s "$HG38_FA" ]] && echo "--lifted-hg38 $LIFTED38 --reference-hg38 $HG38_FA")} \
  --out "$OUT" || die "archaic-panel failed"

log "stage 8 complete: $OUT"
log "NEXT (design §10, M1 checkpoints — neither is automated):"
log "  A. calibrate --max-afr-freq / --min-non-afr-freq against the hmmix Zenodo set before shipping"
log "  B. intersect $TMP/archaic_markers.${BUILD}.tsv with a 23andMe v5 raw file to measure the"
log "     real chip call rate — that number is the honest ceiling for the reported count"
log "Then run 09 (manifest) so the asset is verifiable, and do NOT publish the raw archaic inputs."
