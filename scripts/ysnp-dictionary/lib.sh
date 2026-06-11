# shellcheck shell=bash
# Shared helpers for the Y-SNP dictionary pipeline. Source AFTER config.sh.

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

require_tool() {
  local tool="$1" hint="${2:-}"
  command -v "$tool" >/dev/null 2>&1 || die "missing tool: $tool ${hint:+($hint)}"
}

ensure_dirs() { mkdir -p "$RAW" "$TMP" "$ASSETS" "$LOG"; }

# Resumable download to $RAW/<name> (skips if present and non-empty).
fetch() {
  local url="$1" name="${2:-$(basename "$1")}" dest="$RAW/${2:-$(basename "$1")}"
  if [[ -s "$dest" ]]; then log "have $name (skip)"; return 0; fi
  log "fetch $name <- $url"
  if ! curl -fL --retry 3 --retry-delay 5 -C - -o "$dest.part" "$url"; then
    log "download failed: $url"; rm -f "$dest.part"; return 1
  fi
  mv "$dest.part" "$dest"
}

# Strip surrounding double-quotes from a field value (YBrowse CSV is fully quoted).
# Used inside awk via gsub; provided here for shell-side use if needed.
unquote() { sed -e 's/^"//' -e 's/"$//'; }

# Parse a YBrowse quoted CSV ($1) into dictionary rows for build key ($2):
#   name<TAB>build<TAB>chrom<TAB>position<TAB>strand<TAB>ancestral<TAB>derived
# Keeps only single-base A/C/G/T SNP alleles on the Y. The seqid varies by extract
# (`chrY` in hg38, `hg19ChrY` in hg19, sometimes `Y`) — any Y spelling is accepted and the
# emitted chrom is normalized to $Y_CONTIG. Alleles are the + strand of the source build.
ybrowse_to_rows() {
  local csv="$1" build="$2"
  awk -F',' -v build="$build" -v ycontig="$Y_CONTIG" '
    NR==1 { next }                                   # header
    {
      for (i=1; i<=12; i++) { gsub(/^"|"$/, "", $i) }
      seqid=$1; start=$4; strand=$7; name=$9; anc=$11; der=$12;
      if (seqid !~ /(^|[^A-Za-z])[Cc]hr[Yy]$/ && seqid != "Y" && seqid != "chrY") next;
      if (name == "" || start !~ /^[0-9]+$/) next;
      if (anc !~ /^[ACGTacgt]$/ || der !~ /^[ACGTacgt]$/) next;
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", name, build, ycontig, start, strand, toupper(anc), toupper(der);
    }' "$csv"
}
