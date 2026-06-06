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
  local insecure=(); [[ "$url" == *reichdata.hms.harvard.edu* ]] && insecure=(-k)
  if ! curl -fL "${insecure[@]}" --retry 3 --retry-delay 5 -C - -o "$dest.part" "$url"; then
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
  local fa="$RAW/chm13v2.0.fa"
  [[ -s "$fa" ]] || die "CHM13 FASTA not unpacked at $fa (run 01_fetch.sh)"
  log "liftover $(basename "$in") ($src -> $BUILD)"
  CrossMap vcf "$chain" "$in" "$fa" "$TMP/$(basename "$out" .gz)" \
    || die "CrossMap failed on $in"
  # Align alleles to the CHM13 reference; -c s swaps/flips ref-mismatched records,
  # drops what can't be reconciled. Ancient pseudo-haploid GTs pass through.
  bcftools norm -c s -f "$fa" "$TMP/$(basename "$out" .gz)" -Oz -o "$out" \
    || die "bcftools norm failed on $in"
  tabix -f -p vcf "$out"
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
