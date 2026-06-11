# Y-SNP name → locus dictionary build

Builds the **SNP-name → locus** asset that gives a BISDNA (or any name-only Y panel)
export its missing coordinates. A name like `CTS10003` resolves to a position plus
ancestral/derived alleles, **per reference build** — so the codebase stays
build-agnostic. Consumed by `navigator-domain::ysnp_dict` (the loader) and the BISDNA
importer. See `docs/design/bisdna-import.md`.

## Source

[YBrowse](https://ybrowse.org/) (Thomas Krahn) — the canonical Y-SNP catalog, ~2M
named variants with positions, ancestral/derived alleles, and strand. Native extracts
cover **GRCh38** and **GRCh37**; **CHM13v2 (`hs1`)** is added by lifting the GRCh38
coordinate with the same chain the ancestry pipeline uses. ~95% of Y-SNPs lift
cleanly; the rest are simply absent on `hs1` (placeable on GRCh38, dropped for CHM13).

## Asset format

`$ASSETS/dictionary.tsv` (default `~/.decodingus/ysnp/dictionary.tsv`) — one row per
(name, build):

```
name	build	chrom	position	strand	ancestral	derived
CTS10003	GRCh38	chrY	15311491	+	C	T
CTS10003	hs1	chrY	14800123	+	C	T
```

- **Alleles are the + strand of the named build.** A GRCh38→hs1 lift that inverts
  strand has its alleles complemented at build time, so the importer/placement always
  compares against the sample's reference + strand base.
- `aliases.tsv` is header-only by default: YBrowse lists synonymous SNP names as
  separate rows at the same locus, so every name is already a direct entry. Add manual
  `alias<TAB>canonical` overrides there if ever needed.

## Stages

```bash
cd scripts/ysnp-dictionary
./01_fetch.sh                          # YBrowse CSVs + GRCh38→CHM13 chain + CHM13 FASTA
./02_build.sh                          # → $ASSETS/dictionary.tsv (+ aliases.tsv)
./03_restrict_panel.sh <results.txt>   # optional: → $ASSETS/chromo2-panel.tsv (checked-in manifest)
```

`config.sh` holds every URL/path (override via env). URLs marked `# VERIFY` are
best-known as of writing — confirm before a production run (YBrowse refreshes ~weekly;
UCSC chain filenames can roll).

### Stage 3 — the checked-in panel manifest

The chromo2 chip probes a fixed ~14k SNP set. `03_restrict_panel.sh` filters the full
dictionary to just the names a given `results.txt` carries, producing a small
`chromo2-panel.tsv` (same format as `dictionary.tsv`) that can be **committed to the
repo** so BISDNA import works offline. It also reports any chip names the dictionary
couldn't resolve (these can never be placed — investigate before shipping the manifest).

## Requirements

`curl`, `awk`, `sort`, [`CrossMap`](https://github.com/liguowang/CrossMap)
(`pip install CrossMap`), and `gunzip`; `samtools` optional (FASTA `.fai`).

## Licensing

Confirm YBrowse redistribution terms before committing a derived `chromo2-panel.tsv`
into the repo (see `docs/design/bisdna-import.md` §8).
