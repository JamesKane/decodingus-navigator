# BISDNA Chromo2 Y-SNP Import — Design

**Status:** Design only (no implementation). Drafted 2026-06-10.
**Branch context:** `rust-rewrite`. Targets `navigator-domain`, `navigator-store`,
`navigator-app`, `navigator-ui`, plus a new asset-build script under `scripts/`.

## 1. Goal

Import a **BISDNA chromo2** Y-chromosome chip export (`results.txt`) as real
**Y-SNP variant calls** that can drive haplogroup placement and cross-source
reconciliation — not as a throwaway QC summary.

The file gives a **SNP name**, an Illumina genotype, and a positive/negative
verdict per marker. It does **not** give positions or ancestral/derived alleles.
The crux of this feature is the missing piece: a **SNP-name → locus dictionary**
that resolves each name to a coordinate (in whatever reference build we're placing
against) and its ancestral/derived alleles.

### Why this matters now

The CLI ingest work (commit `8e26317`) surfaced that this file currently
**misroutes to `import_chip_profile_from_csv`**, which computes a
`ChipSummary` (call/no-call counts) and **discards every per-SNP call**. For a
Y-only panel whose entire value *is* the named SNP calls, that throws away the
signal. BISDNA must instead produce a `VariantSet` of Y-SNP `VariantCall`s.

## 2. The source file

`/Volumes/nas/Genomics/mine/results.txt` — BISDNA chromo2 Y raw data. Real tallies
from the actual file (14,219 markers):

```
header (line 1): a paragraph of prose describing the format
header (line 2): SNPID<TAB>genotype<TAB>result
data:            CTS10003<TAB>CC<TAB>negative
```

| `result` value | count  | meaning                                              |
|----------------|--------|------------------------------------------------------|
| `negative`     | 13,858 | ancestral allele carried                             |
| `positive`     | 341    | derived allele carried                               |
| `(positive)`   | 8      | positive on a back-mutation-prone marker (still derived) |
| `no_call`      | 12     | genotype `00`, undetermined — drop                   |

Genotypes are **two characters** (`AA`, `GG`, `AG`, `AT`, `00`). Key facts from the
file's own preamble:

- Alleles are called on the **Illumina TOP strand**, *not* necessarily the
  reference + strand. (This is the strand trap — see §5.)
- `AA`/`GG` are haploid calls written doubled (Illumina autosomal clustering).
- Apparent heterozygotes (`AG`, `AT`) on a positive marker are still **positive
  calls** — clustering artifacts, not real hets.
- `(positive)` in parentheses = a back-mutation marker (e.g. S163) called derived.
  A genuinely back-mutated call would read `back-mutated`; none appear in this file
  but the parser must tolerate the label.

## 3. Target data model

BISDNA is a genotyping array, but a **Y-only** one. Map it to the existing
variant-set model, *not* the chip-summary model:

- **Run:** `SequenceRun.test_type = "ARRAY_BISDNA"` — already in the catalog
  (`testtype.rs`, target `YChromosome`). No new test type needed.
- **Data:** a `VariantSet` (`variants.rs`) of Y-SNP `VariantCall`s:
  - `contig` = the Y contig **in the emitted build** (`chrY` for GRCh38/CHM13).
  - `position` = the SNP's locus **in the emitted build** (see §4 build-agnosticism).
  - `reference` = ancestral allele, `alternate` = derived allele (both + strand of
    the emitted build).
  - `genotype` = `"1"` for a derived (positive) call, `"0"` for ancestral
    (negative). Haploid Y — single allele, no `0/1` diploid notation.
  - `rs_id` = `None` (BISDNA carries no rsIDs).
- **Source type:** `SourceType::Chip` (`snp_weight = 0.5`). It *is* a chip; the 0.5
  reconciliation weight is correct relative to Sanger (1.0) and WGS (0.85–0.95).

This drops straight into `haplo::score`, which consumes a `HashMap<i64 position,
char base>`: for each call, `position → derived_char` (positive) or
`ancestral_char` (negative). It also drops into `reconcile_variants` for
cross-source SNP concordance.

## 4. The SNP-name dictionary (the core deliverable)

A lookup from **SNP name → locus**, build-agnostic by construction.

### 4.1 Schema (build-agnostic)

Mirror the shape the DecodingUs Y-tree already uses (`coordinates` keyed by build
label — `hs1`, `GRCh38`, `GRCh37`), so we never bake a single build in:

```
SnpEntry {
  name:        String,                 // canonical (e.g. "CTS10003")
  aliases:     Vec<String>,            // synonyms ("M269" == "PF6517", ...)
  coordinates: Map<BuildKey, Coord>,   // "GRCh38" | "GRCh37" | "hs1" | ...
}
Coord {
  chrom:     String,    // "chrY"
  position:  i64,
  strand:    char,      // '+' / '-' on this build (for the strand check, §5)
  ancestral: String,    // + strand of THIS build
  derived:   String,    // + strand of THIS build
}
```

Per-build alleles matter: a `-` strand liftover flips the allele, so ancestral/
derived are stored **per coordinate**, not once. The importer is handed a **target
build** (the build the subject is being placed against — from the subject's
alignment, the active Y-tree provider's build, or an explicit arg) and reads
`coordinates[target_build]`. CHM13v2 (`hs1`) is just one entry among GRCh38/GRCh37;
nothing is CHM13-specific.

### 4.2 Source — YBrowse master index (recommended)

[YBrowse](https://ybrowse.org/) (Thomas Krahn) is the canonical Y-SNP name catalog
— **2M+ names** with positions, ancestral/derived alleles, strand, and aliases. The
master extract `snps_hg38.csv` is GFF-shaped:

```
chrom, source(point), class(snp), start, end, score, strand, ..., name, name, ancestral, derived, haplogroup
```

It ships **hg38 and hg19** (`snps_hg38.csv`, `snps_hg19.csv`) natively → populates
the `GRCh38` and `GRCh37` coordinates directly. **CHM13 (`hs1`)** is added by
**liftover of the GRCh38 coordinate** through the existing Y chain
(`navigator-refgenome` already does GRCh38→CHM13 Y liftover with rev-comp handling —
the same path validated for GFX0457637). ~95% of Y-SNPs lift cleanly; the rest fall
in structurally divergent regions and are simply absent from the `hs1` map (dropped
for CHM13 placement, still usable on GRCh38).

### 4.3 Packaging — versioned asset, not hand-edited config

The full catalog is ~2M rows — a **data asset**, not a config file you edit by hand.
Follow the established asset pattern (`scripts/ancestry-panel/`, assets cached under
`~/.decodingus/`):

- `scripts/ysnp-dictionary/` — staged build: fetch `snps_hg38.csv` + `snps_hg19.csv`
  → parse → liftover hg38→hs1 → emit a compact, sorted asset
  `~/.decodingus/ysnp/dictionary.tsv` (name, build, chrom, pos, strand, anc, der)
  plus an `aliases.tsv`. Regenerable, version-stamped (YBrowse is updated ~weekly).
- Loader in `navigator-domain` (or `navigator-refgenome`) reads the asset into an
  in-memory `HashMap<name → SnpEntry>` (case-folded; aliases resolved). ~2M entries
  is fine in memory; or restrict to the chromo2 panel (§4.4) to shrink it.

**Optional config knob** (this is the genuinely user-facing "configuration file"
the original note gestured at): a small TOML/conf entry pointing at the asset
location + pinned YBrowse version, alongside `test_types.conf` / `feature_toggles.conf`.
The bulk data stays in the generated asset; the conf just pins which asset.

### 4.4 Chip-panel manifest (optional optimization)

The chromo2 chip probes a **fixed ~14k SNP set**. We can pre-filter the dictionary
to just those names and ship a ~14k-row `chromo2-panel.tsv` asset — small enough to
**check into the repo**. This makes BISDNA import work offline with no 2M-row asset,
and doubles as the chip's design manifest. The full YBrowse asset remains the
fallback for names outside the captured panel. *Recommended* as the default path;
the full asset is the superset for other Y chips later.

### 4.5 Fallback — tree-derived dictionary

Every loaded Y-tree (FTDNA or DecodingUs) already carries name + position +
ancestral/derived for its **branch-defining** SNPs (`haplo::Locus`). When the
dictionary asset is absent, we can build a name→locus map from the active tree
covering exactly the SNPs that matter for placement (a subset of the 14k, but the
load-bearing subset). Useful as a zero-asset degraded mode; not a full replacement
(it can't resolve private/equivalent SNPs the tree doesn't define).

## 5. Strand & result handling (correctness-critical)

The Illumina TOP strand is the trap: a `GG` genotype does **not** mean the + strand
base is G. Two ways to resolve it, in order of robustness:

**Primary — trust the `result` column.** BISDNA has already done the
positive/negative determination against its own probe design. So:

- `positive` / `(positive)` → emit the dictionary's **derived** allele, genotype `1`.
- `negative` → emit the dictionary's **ancestral** allele, genotype `0`.
- `no_call` (`00`) → drop.
- `back-mutated` → the lineage is derived but the base reads ancestral. Flag and
  **exclude from the placement `calls` map** (a position→base model can't represent
  "derived lineage showing ancestral base"); retain in the `VariantSet` with a note
  so reconciliation can see it. Rare (0 in this file).

This sidesteps strand entirely — we emit the build-correct + strand allele straight
from the dictionary, keyed only on the boolean verdict.

**Secondary — cross-check the genotype** (a QC signal, not the source of truth).
Take the called base (`GG`→`G`; for a het positive like `AG`, the non-ancestral
base), compare it (and its complement, since strand is unknown) against the
dictionary derived/ancestral. A mismatch on both strands flags a probable
name-collision or stale dictionary entry. Surface as an import warning; don't fail
the row.

**Unknown names.** A SNP name absent from the dictionary can't be placed — count
and report (`N markers unresolved`), don't silently drop. (The CLI memory note's
"2 bogus markers" problem came from exactly this kind of silent mishandling.)

## 6. Detection & routing

Today `filetype::detect` scores STR vs chip and BISDNA lands as `ChipData` (or
`Unknown`). Add a **dedicated detector** ahead of the STR/chip scorer — the BISDNA
fingerprint is unambiguous:

- header line 2 is `SNPID\tgenotype\tresult` (case-insensitive), **and/or**
- the third column is dominated by `positive`/`negative`/`no_call`/`back-mutated`,
  **and** the first column holds named Y-SNPs (CTS/M/L/S/PF/Z/FGC/Y… + digits),
  **and** no `rs` IDs / chromosome-position columns.

Add `DetectedData::YSnpPanel` and route it to the new BISDNA importer. (BISDNA's
multi-line prose preamble also needs the head reader to skip non-tabular lead text —
the parser keys off the `SNPID/genotype/result` header row, not line 1.)

## 7. Module & wiring plan

- **`navigator-domain/src/bisdna.rs`** (new): pure parser.
  `parse(text) -> Result<Vec<BisdnaCall>>` where `BisdnaCall { name, genotype,
  verdict }`; `verdict ∈ {Positive, Negative, NoCall, BackMutated}`. No IO, no
  dictionary — just the file. Unit-tested against the real header + sampled rows.
- **`navigator-domain/src/ysnp_dict.rs`** (new): dictionary types (§4.1) + asset
  loader + `resolve(name, build) -> Option<Coord>` (alias-aware). Pure over an
  already-loaded asset; the fetch/build lives in the script + app layer.
- **`navigator-app`**: `import_bisdna_from_file(biosample_guid, path, build)` —
  parse → resolve each name against the dictionary for `build` → build
  `VariantCall`s → persist as a `VariantSet` (`SourceType::Chip`,
  `source_label = file name`), create/attach the `ARRAY_BISDNA` run. Report
  resolved / unresolved / no-call / strand-mismatch counts.
- **`add_data`**: route `DetectedData::YSnpPanel → import_bisdna_from_file`. Build
  comes from the subject's existing Y alignment if present, else the active Y-tree
  provider's build, else a default (configurable).
- **`filetype::detect`**: add the BISDNA detector + `YSnpPanel` variant (§6).
- **CLI `ingest`**: inherits the routing fix automatically; add coverage to the
  headless path that surfaced the bug.
- **`navigator-ui`**: a Data Sources card for the Y-SNP panel (resolved-marker
  count, unresolved warnings, derived-call count) — mirrors the chip card.

## 8. Open decisions

1. **Dictionary source / packaging** — *recommend* §4.4: ship a checked-in ~14k-row
   chromo2 panel manifest (built once from YBrowse + liftover), with the full
   regenerable YBrowse asset (§4.2–4.3) as the superset/fallback. Alternative: full
   asset only (heavier, no repo bloat, needs the build step before first import).
2. **Back-mutated handling** — flag-and-exclude from placement (§5) vs. model a
   richer call state. Recommend flag-and-exclude for v1 (0 cases in the sample).
3. **Emitted build** — auto-derive from the subject's alignment / tree provider vs.
   always emit all available builds into the VariantSet. Recommend single
   target-build emission (reconciliation compares within one build); revisit if
   subjects mix builds across sources.
4. **YBrowse licensing/attribution** for a bundled manifest — confirm redistribution
   terms before checking the panel into the repo.

## 9. Phased plan

1. **Parser** — `bisdna.rs` + tests against the real file. No dictionary yet.
2. **Dictionary asset** — `scripts/ysnp-dictionary/` build + `ysnp_dict.rs` loader;
   produce the chromo2 panel manifest; validate it resolves ~all 14k names.
3. **Importer + routing** — `import_bisdna_from_file`, `YSnpPanel` detection,
   `add_data` wiring; validate against `results.txt` end-to-end.
4. **Placement validation** — feed the emitted calls to `haplo::score` and confirm
   a sane terminal haplogroup (sanity-check against the file's positive SNPs and any
   known result for this kit).
5. **UI + reconciliation** — Data Sources card; confirm the calls participate in
   `reconcile_variants` alongside any WGS/Sanger sources for the same subject.

---

**Sources:** [YBrowse](https://ybrowse.org/) · [ISOGG Y-SNP indexes](https://isogg.org/wiki/Y-SNP_indexes) · YBrowse master extracts `snps_hg38.csv` / `snps_hg19.csv`.
