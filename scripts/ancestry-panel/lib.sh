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

# Liftover a VCF to CHM13 with CrossMap, then normalise alleles against the CHM13
# reference (set REF to the assembly allele, swap/flip REF<->ALT as needed, drop
# irreconcilable sites). Produces a bgzipped, tabix-indexed VCF aligned to CHM13.
#   liftover_vcf <in.vcf[.gz]> <source-build> <out.vcf.gz>
liftover_vcf() {
  local in="$1" src="$2" out="$3" chain; chain="$(chain_for "$src")"
  local fa="$RAW/chm13v2.0.fa" dict="$RAW/chm13v2.0.dict"
  [[ -s "$fa" ]] || die "CHM13 FASTA not unpacked at $fa (run 01_fetch.sh)"
  [[ -s "$dict" ]] || samtools dict "$fa" -o "$dict" || die "samtools dict failed for $fa"
  # Allele-aware liftover via GATK/Picard LiftoverVcf — NOT CrossMap. CrossMap's `vcf` mode blanks
  # the ALT allele whenever the target (CHM13) reference base differs from the source REF (i.e. most
  # polymorphic sites, and CHM13's inverted/rearranged segments), silently destroying ~3/4 of the
  # genotypes. LiftoverVcf reverse-complements alleles on minus-strand chain blocks and, with
  # RECOVER_SWAPPED_REF_ALT, swaps REF<->ALT and flips the genotypes when the target ref equals the
  # ALT (rather than dropping it); irreconcilable sites go to a REJECT file. Output is already
  # target-ref-correct and coordinate-sorted.
  log "liftover $(basename "$in") ($src -> $BUILD) via GATK LiftoverVcf"
  local lift="$TMP/$(basename "$out" .gz).lift.vcf.gz"
  # MAX_RECORDS_IN_RAM bounds the in-memory sort so it spills to disk — without it LiftoverVcf holds
  # every record (here ~20k x 23k genotypes) in RAM and OOMs. Heap + spill dir are overridable.
  gatk --java-options "${GATK_JAVA_OPTS:--Xmx16g}" LiftoverVcf -I "$in" -O "$lift" -C "$chain" -R "$fa" \
    --REJECT "$TMP/$(basename "$out" .gz).reject.vcf.gz" \
    --RECOVER_SWAPPED_REF_ALT true --WARN_ON_MISSING_CONTIG true \
    --MAX_RECORDS_IN_RAM "${GATK_MAX_RECORDS_IN_RAM:-2000}" --TMP_DIR "$TMP" \
    --CREATE_INDEX false > "$LOG/liftover_$(basename "$in").log" 2>&1 \
    || die "GATK LiftoverVcf failed on $in (see $LOG/liftover_$(basename "$in").log)"
  bcftools view "$lift" -Oz -o "$out" && tabix -f -p vcf "$out" \
    || die "bcftools index failed on $in"
  rm -f "$lift"
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
