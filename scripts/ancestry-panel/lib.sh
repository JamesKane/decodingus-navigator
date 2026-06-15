# shellcheck shell=bash
# Shared helpers for the ancestry-panel pipeline. Source AFTER config.sh.

# Timestamped logging to stderr.
log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

# Ensure a command exists, with an install hint.
require_tool() {
  local tool="$1" hint="${2:-}"
  command -v "$tool" >/dev/null 2>&1 || die "missing tool: $tool ${hint:+($hint)}"
}

# Create the working tree (idempotent).
ensure_dirs() {
  mkdir -p "$RAW" "$TMP" "$ASSETS" "$LOG" "$KGP_CHM13_DIR"
}

# Resumable download to $RAW/<name> (skips if present and non-empty).
# Returns non-zero (instead of die) on failure so callers can make a source optional.
fetch() {
  local url="$1" name="${2:-$(basename "$1")}" dest="$RAW/${2:-$(basename "$1")}"
  if [[ -s "$dest" ]]; then log "have $name (skip)"; return 0; fi
  log "fetch $name <- $url"
  # reichdata.hms.harvard.edu serves a broken TLS chain — allow insecure for that host only.
  # Expand with the ${arr[@]+…} guard so an empty array is safe under `set -u` on bash 3.2 (macOS default).
  local insecure=(); [[ "$url" == *reichdata.hms.harvard.edu* ]] && insecure=(-k)
  if ! curl -fL --no-progress-meter ${insecure[@]+"${insecure[@]}"} --retry 3 --retry-delay 5 -C - -o "$dest.part" "$url"; then
    log "download failed: $url"; rm -f "$dest.part"; return 1
  fi
  mv "$dest.part" "$dest"
}

# sha256 of a file (portable: shasum on macOS, sha256sum on Linux).
sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}';
  else shasum -a 256 "$1" | awk '{print $1}'; fi
}

# CHM13 chain selector by source build (grch38|hg19) → chain file path under $RAW.
chain_for() {
  case "$1" in
    hg19|GRCh37|grch37)  echo "$RAW/$(basename "$CHAIN_HG19_TO_CHM13" .gz)";;
    grch38|GRCh38|hg38)  echo "$RAW/$(basename "$CHAIN_GRCH38_TO_CHM13" .gz)";;
    *) die "no CHM13 chain configured for source build: $1";;
  esac
}

# Allele-aware VCF liftover to an arbitrary target reference via GATK/Picard LiftoverVcf — NOT
# CrossMap. CrossMap's `vcf` mode blanks the ALT allele whenever the target reference base differs
# from the source REF (most polymorphic sites, and inverted/rearranged segments), silently
# destroying the genotypes. LiftoverVcf reverse-complements alleles on minus-strand chain blocks
# and, with RECOVER_SWAPPED_REF_ALT, swaps REF<->ALT + flips genotypes when the target ref equals
# the ALT; irreconcilable sites go to a REJECT file. Output is target-ref-correct and sorted.
# MAX_RECORDS_IN_RAM bounds the in-memory sort (spills to disk) so wide many-sample VCFs don't OOM.
#   gatk_lift <in.vcf.gz> <chain> <target.fa> <out.vcf.gz>
gatk_lift() {
  local in="$1" chain="$2" fa="$3" out="$4" dict="${3%.fa}.dict"
  [[ -s "$fa" ]] || die "target FASTA not found at $fa"
  [[ -s "$dict" ]] || samtools dict "$fa" -o "$dict" || die "samtools dict failed for $fa"
  local lift="$TMP/$(basename "$out" .gz).lift.vcf.gz"
  gatk --java-options "${GATK_JAVA_OPTS:--Xmx16g}" LiftoverVcf -I "$in" -O "$lift" -C "$chain" -R "$fa" \
    --REJECT "$TMP/$(basename "$out" .gz).reject.vcf.gz" \
    --RECOVER_SWAPPED_REF_ALT true --WARN_ON_MISSING_CONTIG true \
    --MAX_RECORDS_IN_RAM "${GATK_MAX_RECORDS_IN_RAM:-2000}" --TMP_DIR "$TMP" \
    --CREATE_INDEX false > "$LOG/liftover_$(basename "$out").log" 2>&1 \
    || die "GATK LiftoverVcf failed ($in -> $out); see $LOG/liftover_$(basename "$out").log"
  bcftools view "$lift" -Oz -o "$out" && tabix -f -p vcf "$out" || die "bcftools index failed on $out"
  rm -f "$lift"
}

# Liftover a VCF to CHM13 (the project default target). Wrapper over gatk_lift.
#   liftover_vcf <in.vcf[.gz]> <source-build> <out.vcf.gz>
liftover_vcf() {
  local in="$1" src="$2" out="$3"
  log "liftover $(basename "$in") ($src -> $BUILD) via GATK LiftoverVcf"
  gatk_lift "$in" "$(chain_for "$src")" "$RAW/chm13v2.0.fa" "$out"
}

# Slice a (possibly remote) VCF/BCF down to <regions> (a CHROM<TAB>POS tsv or BED), writing a
# bgzipped, tabix-indexed VCF. htslib streams only the indexed byte ranges, so a multi-TB remote
# callset costs only the panel-site records. The source index (.tbi/.csi) must sit next to <src>.
# For gs:// requester-pays buckets, export GCS_REQUESTER_PAYS_PROJECT before calling.
#   slice_at <src-url-or-path> <regions> <out.vcf.gz>
slice_at() {
  local src="$1" regions="$2" out="$3"
  [[ -s "$out" ]] && { log "have $(basename "$out") (skip)"; return 0; }
  [[ -s "$regions" ]] || { log "slice skipped: regions $regions missing"; return 1; }
  log "slice $(basename "$src") @ $(basename "$regions") -> $(basename "$out")"
  if ! bcftools view -R "$regions" -Oz -o "$out.part" "$src"; then
    log "slice failed: $src"; rm -f "$out.part"; return 1
  fi
  mv "$out.part" "$out"; tabix -f -p vcf "$out"
}

# Extract a panelbuild genotype matrix (CHROM POS REF ALT GT...) at the panel sites,
# plus the parallel sample list.  bcftools query is the format panelbuild ingests.
#   matrix_from_vcf <in.vcf.gz> <sites.tsv> <out-matrix.tsv.gz> <out-samples.txt>
matrix_from_vcf() {
  local in="$1" sites="$2" mat="$3" samp="$4"
  bcftools query -l "$in" > "$samp"
  bcftools query -R "$sites" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$in" \
    | gzip > "$mat" || die "bcftools query failed on $in"
}

# Emit GRCh37/hg19 ##contig header lines (chr-prefixed, standard lengths) for a sites-only VCF.
hg19_contigs() {
  while read -r c l; do echo "##contig=<ID=$c,length=$l>"; done <<'EOF'
chr1 249250621
chr2 243199373
chr3 198022430
chr4 191154276
chr5 180915260
chr6 171115067
chr7 159138663
chr8 146364022
chr9 141213431
chr10 135534747
chr11 135006516
chr12 133851895
chr13 115169878
chr14 107349540
chr15 102531392
chr16 90354753
chr17 81195210
chr18 78077248
chr19 59128983
chr20 63025520
chr21 48129895
chr22 51304566
chrX 155270560
chrY 59373566
EOF
}

# Build the IBD genetic-map asset: lift the GRCh38 recombination map to CHM13 (coordinate-only —
# no alleles, so CrossMap is fine here) and serialize via panelbuild. cM is made monotonic per
# chromosome after the lift (a lifted position can precede an earlier one).
#   build_genetic_map <out.bin>
build_genetic_map() {
  local out="$1" gmdir="$RAW/gmap_grch38/chr_in_chrom_field"
  [[ -d "$gmdir" ]] || { log "WARN: genetic-map input $gmdir missing (run 01_fetch.sh) — skipping genetic map"; return 1; }
  require_tool CrossMap; require_tool cargo
  local chain="$RAW/$(basename "$CHAIN_GRCH38_TO_CHM13" .gz)"
  log "genetic map: GRCh38 -> CHM13 liftover"
  awk '{printf "%s\t%d\t%d\t%s\n",$1,$4-1,$4,$3}' "$gmdir"/plink.chrchr*.GRCh38.map > "$TMP/gmap.hg38.bed"
  CrossMap bed "$chain" "$TMP/gmap.hg38.bed" "$TMP/gmap.chm13.bed" >/dev/null 2>&1 || die "CrossMap genetic-map failed"
  { printf 'chrom\tpos\tcM\n'
    awk '{print $1"\t"$3"\t"$4}' "$TMP/gmap.chm13.bed" | sort -k1,1 -k2,2n \
      | awk -F'\t' '{if($1!=c){c=$1;mx=-1} v=$3+0; if(v<mx)v=mx; else mx=v; print $1"\t"$2"\t"v}'
  } > "$TMP/genetic_map.chm13.txt"
  cargo run --release -q -p navigator-panelbuild -- genetic-map --input "$TMP/genetic_map.chm13.txt" --out "$out" \
    || die "panelbuild genetic-map failed"
}

# Build the chip-compatible, multi-build IBD panel from the AADR 1240k site set (a strong consumer-
# array backbone). GRCh37 is native (.snp); CHM13 via GATK allele-aware lift; GRCh38 too if the hg38
# FASTA is present (else those columns are omitted). panelbuild drops strand-ambiguous palindromes.
#   build_ibd_panel <out.bin>
build_ibd_panel() {
  local out="$1" snp; snp="$(ls "$RAW/"*"${AADR_DATASET}"*.snp 2>/dev/null | head -1 || true)"
  [[ -s "$snp" ]] || { log "WARN: AADR .snp not found — skipping IBD panel"; return 1; }
  require_tool gatk; require_tool samtools; require_tool cargo
  # acgt_snp: keep ACGT bi-allelic autosome+XY rows; emit "chrN<TAB>pos<TAB>rsid<TAB>ref<TAB>alt".
  local acgt='function ok(c){return c=="X"||c=="Y"||(c~/^[0-9]+$/&&c>=1&&c<=22)}
    {c=$2; if(c==23)c="X"; else if(c==24)c="Y"; r=$5;a=$6; if(ok(c)&&r~/^[ACGT]$/&&a~/^[ACGT]$/) print "chr"c"\t"$4"\t"$1"\t"r"\t"a}'
  log "IBD panel: hg19 sites VCF from $(basename "$snp")"
  { echo '##fileformat=VCFv4.2'; hg19_contigs; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    LC_ALL=C awk "$acgt" "$snp" | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t.\t.\t.\n",$1,$2,$3,$4,$5}'
  } > "$TMP/ibd.hg19.vcf"
  bgzip -f "$TMP/ibd.hg19.vcf"
  bcftools sort "$TMP/ibd.hg19.vcf.gz" -Oz -o "$TMP/ibd.hg19.sorted.vcf.gz" 2>/dev/null
  log "IBD panel: GATK lift hg19 -> CHM13"
  gatk_lift "$TMP/ibd.hg19.sorted.vcf.gz" "$RAW/$(basename "$CHAIN_HG19_TO_CHM13" .gz)" "$RAW/chm13v2.0.fa" "$TMP/ibd.chm13.vcf.gz"
  bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' "$TMP/ibd.chm13.vcf.gz" > "$TMP/ibd.chm13.tsv"
  LC_ALL=C awk "$acgt" "$snp" | awk -F'\t' '{print $3"\t"$1"\t"$2"\t"$4"\t"$5}' > "$TMP/ibd.grch37.tsv"  # rsid c p r a
  local hg38fa="$RAW/$(basename "$GRCH38_FASTA_URL" .gz)" hg38chain="$RAW/$(basename "$CHAIN_HG19_TO_HG38" .gz)" have38=0
  if [[ -s "$hg38fa" && -s "$hg38chain" ]]; then
    log "IBD panel: GATK lift hg19 -> GRCh38"
    gatk_lift "$TMP/ibd.hg19.sorted.vcf.gz" "$hg38chain" "$hg38fa" "$TMP/ibd.hg38.vcf.gz"
    bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' "$TMP/ibd.hg38.vcf.gz" > "$TMP/ibd.grch38.tsv"; have38=1
  else
    log "IBD panel: no GRCh38 FASTA — grch38 columns omitted (GRCh37 chips + WGS still covered)"
  fi
  log "IBD panel: join multi-build sites table"
  { if [[ "$have38" == 1 ]]; then
      printf 'rsid\tchm13_contig\tchm13_pos\tchm13_ref\tchm13_alt\tgrch37_contig\tgrch37_pos\tgrch37_ref\tgrch37_alt\tgrch38_contig\tgrch38_pos\tgrch38_ref\tgrch38_alt\n'
      LC_ALL=C awk -F'\t' '
        FILENAME ~ /grch37/ {g37[$1]=$2"\t"$3"\t"$4"\t"$5; next}
        FILENAME ~ /grch38/ {g38[$1]=$2"\t"$3"\t"$4"\t"$5; next}
        ($1 in g37){ e=($1 in g38)?g38[$1]:"\t\t\t"; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,g37[$1],e }
      ' "$TMP/ibd.grch37.tsv" "$TMP/ibd.grch38.tsv" "$TMP/ibd.chm13.tsv"
    else
      printf 'rsid\tchm13_contig\tchm13_pos\tchm13_ref\tchm13_alt\tgrch37_contig\tgrch37_pos\tgrch37_ref\tgrch37_alt\n'
      LC_ALL=C awk -F'\t' 'NR==FNR{g[$1]=$2"\t"$3"\t"$4"\t"$5;next} ($1 in g){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"g[$1]}' "$TMP/ibd.grch37.tsv" "$TMP/ibd.chm13.tsv"
    fi
  } > "$TMP/ibd_sites_multibuild.tsv"
  cargo run --release -q -p navigator-panelbuild -- ibd-panel --input "$TMP/ibd_sites_multibuild.tsv" --out "$out" --build "$BUILD" \
    || die "panelbuild ibd-panel failed"
}
