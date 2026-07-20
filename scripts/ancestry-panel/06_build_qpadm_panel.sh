#!/usr/bin/env bash
# Stage 6 — build the SHIPPING deep-ancestry asset: the full-1240k qpAdm panel (Patterson et al. 2022
# config, docs/design/ancient-ancestry-rebuild.md §7.14). Unlike the 20k AIM assets (stages 03-05),
# this is built over the WHOLE 1240k SNP set — the precision the qpAdm f4 model needs (±1-2% vs ±40%
# at 20k). It reads the full-1240k AADR PLINK from stage 04 (convertf → PACKEDPED), NOT the 20k
# matrices, so its own pipeline: subset the AADR to the source+outgroup individuals, lift hg19→CHM13
# by rsID, and build the per-population AF panel.
#
# Produces $QPADM_OUT (ancestry_qpadm_<build>.bin) — an AncestryPanel with populations
#   [<sources…> <outgroups…>]  (sources first: the app reads the first len(ANCIENT_COMPONENTS) as the
#   qpAdm "left" set, the rest as the "right"/outgroup set).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"; source "$HERE/lib.sh"
require_tool cargo plink2

AADR_PLINK="$TMP/aadr"                              # stage 04 convertf output (.bed/.bim/.fam, hg19 1240k)
BED_1240K="$TMP/1240k_sites.${BUILD}.bed"           # stage 02 liftover (name field: rsid|ref|alt, CHM13 pos = BED end)
ANNO="$RAW/${AADR_FILE_PREFIX}.anno"                # AADR .anno (Genetic ID col1, Group ID col15)
[[ -f "$AADR_PLINK.bed" ]] || die "missing AADR PLINK ($AADR_PLINK.bed) — run 04_build_matrices.sh (needs convertf)"
[[ -s "$BED_1240K" ]] || die "missing $BED_1240K — run 02_liftover_panel_sites.sh"
[[ -s "$ANNO" ]] || die "missing $ANNO — run 01_fetch.sh"
[[ -n "$ANCIENT_OUTGROUPS" ]] || die "ANCIENT_OUTGROUPS is empty — qpAdm needs sister outgroups (qpadm_rightpops.txt)"

# The qpAdm populations = sources ∪ outgroups (the codes the config carries, comma-separated).
WANT="$ANCIENT_COMPONENTS,$ANCIENT_OUTGROUPS"
log "qpAdm populations: $WANT"

# 1. sample(Genetic ID) -> component label, from the anno Group ID (col15) via the component map,
#    restricted to the qpAdm population codes.
awk -F'\t' -v WANT="$WANT" '
  BEGIN{ n=split(WANT,w,","); for(i=1;i<=n;i++){ gsub(/^ +| +$/,"",w[i]); if(w[i]!="") keep[w[i]]=1 } }
  NR==FNR{ if($2 in keep) lab[$1]=$2; next }
  FNR>1 && ($15 in lab){ print $1"\t"lab[$15] }
' "$QPADM_COMPONENT_MAP" "$ANNO" > "$TMP/qpadm.pops.tsv"
[[ -s "$TMP/qpadm.pops.tsv" ]] || die "no AADR samples mapped to the qpAdm populations — check aadr_component_map.tsv"
log "mapped $(wc -l < "$TMP/qpadm.pops.tsv") individuals: $(cut -f2 "$TMP/qpadm.pops.tsv" | sort | uniq -c | sort -rn | awk '{printf "%s=%s ",$2,$1}')"

# 2. Subset the AADR PLINK to those individuals, relabel FID = component, keep autosomes, export .traw.
#    (--update-ids runs at load, so keep by the ORIGINAL numeric FID first, then relabel.)
awk 'NR==FNR{fid[$2]=$1;next} ($1 in fid){print fid[$1]"\t"$1"\t"$2}' "$AADR_PLINK.fam" "$TMP/qpadm.pops.tsv" > "$TMP/qpadm.selfid.tsv"
awk -F'\t' '{print $1"\t"$2}'          "$TMP/qpadm.selfid.tsv" > "$TMP/qpadm.keep.txt"
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$2}' "$TMP/qpadm.selfid.tsv" > "$TMP/qpadm.updateids.txt"
plink2 --bfile "$AADR_PLINK" --keep "$TMP/qpadm.keep.txt" --make-bed --out "$TMP/qpadm.sub" >/dev/null 2>&1 \
  || die "plink2 subset failed"
plink2 --bfile "$TMP/qpadm.sub" --update-ids "$TMP/qpadm.updateids.txt" --chr 1-22 --output-chr 26 \
  --export A-transpose --out "$TMP/qpadm" >/dev/null 2>&1 || die "plink2 export failed"

# 3. .traw (hg19, per-individual COUNTED-allele counts) -> CHM13 genotype matrix
#    (CHROM POS bed_ref bed_alt GT…), lifting hg19→CHM13 by rsID and orienting every GT to count the
#    1240k ALT allele (so it matches how the app genotypes the target). Sites whose CHM13 position or
#    alleles can't be reconciled are dropped.
awk -F'\t' -v BED="$BED_1240K" '
  function comp(x){ return (x=="A"?"T":x=="T"?"A":x=="C"?"G":x=="G"?"C":x) }
  function gt(v){ return (v=="NA"?"./.":(v==0?"0/0":(v==1?"0/1":"1/1"))) }
  BEGIN{ while((getline l < BED)>0){ split(l,f,"\t"); split(f[4],a,"|");
           rs=a[1]; cc[rs]=f[1]; cp[rs]=f[3]; rf[rs]=a[2]; al[rs]=a[3] } }
  NR==1{ next }
  { rs=$2; if(!(rs in cc)) next;
    counted=$5; alt=$6; ba=al[rs];
    flip=(ba==counted||comp(ba)==counted)?0:((ba==alt||comp(ba)==alt)?1:-1);
    if(flip<0) next;
    line=cc[rs]"\t"cp[rs]"\t"rf[rs]"\t"al[rs];
    for(i=7;i<=NF;i++){ v=$i; if(flip && v!="NA") v=2-v; line=line"\t"gt(v) }
    print line;
  }' "$TMP/qpadm.traw" | gzip > "$TMP/qpadm_chm13.matrix.tsv.gz"
head -1 "$TMP/qpadm.traw" | cut -f7- | tr '\t' '\n' > "$TMP/qpadm.samples.txt"
# The plink2 --export sample names are `<FID>_<IID>` = `<component>_<IID>` (after --update-ids), so the
# pop map keys on the full name and the label is the prefix before the first `_` (components have none).
awk -F'_' '{print $0"\t"$1}' "$TMP/qpadm.samples.txt" > "$TMP/qpadm.popmap.tsv"
log "CHM13 matrix: $(zcat < "$TMP/qpadm_chm13.matrix.tsv.gz" | wc -l) sites × $(wc -l < "$TMP/qpadm.samples.txt") samples"

# 4. Build the AncestryPanel (per-population AF; floor 1 — outgroups are legitimately small, and the
#    f4 jackknife accounts for the noise). Sources first, then outgroups, per --components/--outgroups.
log "panelbuild ancient-panel (qpAdm: src $ANCIENT_COMPONENTS ; out $ANCIENT_OUTGROUPS) -> $QPADM_OUT"
# --reference orients every site so REF = the actual CHM13 base (the matrix inherits hg19 allele
# labels from the liftover BED, which are ~30% swapped relative to CHM13; docs §7.16). This makes the
# asset CHM13-canonical and joinable with the other ancestry assets (super/fine/IBD).
CHM13_FASTA="${CHM13_FASTA:-$RAW/chm13v2.0.fa}"
[[ -s "$CHM13_FASTA" && -s "$CHM13_FASTA.fai" ]] || die "missing indexed CHM13 FASTA ($CHM13_FASTA[.fai]) — needed to orient the panel to CHM13"
cargo run --release -q -p navigator-panelbuild -- ancient-panel \
  --matrix "$TMP/qpadm_chm13.matrix.tsv.gz" --samples "$TMP/qpadm.samples.txt" --pops "$TMP/qpadm.popmap.tsv" \
  --components "$ANCIENT_COMPONENTS" --outgroups "$ANCIENT_OUTGROUPS" \
  --min-called 1 --outgroup-min-called 1 \
  --reference "$CHM13_FASTA" \
  --out "$QPADM_OUT" || die "qpAdm panel build failed"

# Refresh the integrity manifest so the qpAdm asset gets a checksum (only when writing to the live
# $ASSETS dir — a QPADM_OUT override, e.g. a test build, must not touch the shipped manifest).
if [[ "$QPADM_OUT" == "$ASSETS/"* ]]; then
  log "refreshing manifest -> $MANIFEST"
  cargo run --release -q -p navigator-panelbuild -- manifest --dir "$ASSETS" --build "$BUILD" --out "$MANIFEST" \
    || log "WARN: manifest refresh failed (qpAdm asset will fail verification until rebuilt)"
fi

log "stage 6 complete: $QPADM_OUT"
ls -lh "$QPADM_OUT" 2>/dev/null >&2 || true
