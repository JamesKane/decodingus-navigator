# Decoding-Us Navigator User Guide

> **Alpha release.** Navigator is under active development. The analyses below are usable today, but
> outputs, file layouts, and the UI may still change between releases. Where a module's output has not
> yet been independently validated, the guide says so.

Welcome to the **Decoding-Us Navigator**, your private, local companion for advanced genomic analysis. This application lets you explore your DNA data with professional-grade bioinformatics directly on your own computer, keeping your genetic privacy intact while contributing to citizen science.

## Table of Contents
1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Installation & Setup](#installation--setup)
4. [Getting Started](#getting-started)
   - [First-Time Setup: Bringing Your Own Reference Genomes](#first-time-setup-bringing-your-own-reference-genomes)
5. [Core Features](#core-features)
   - [Workspace Management](#workspace-management)
   - [Importing Data](#importing-data)
   - [Project Import (batch)](#project-import-batch-with-the-sidecar-fast-path)
   - [Batch import strategies for existing data collections](#batch-import-strategies-for-existing-data-collections)
   - [Importing an FTDNA group project](#importing-an-ftdna-group-project)
   - [Running Analyses](#running-analyses)
   - [The Branch Report tool](#the-branch-report-tool)
   - [Exporting & Sharing Results](#exporting--sharing-results)
6. [The Command Line](#the-command-line)
7. [Data Management & Privacy](#data-management--privacy)
8. [Settings](#settings)
9. [Advanced Usage](#advanced-usage)
10. [Troubleshooting](#troubleshooting)

---

## Introduction

Decoding-Us Navigator runs a complete bioinformatics stack on your desktop. Unlike cloud services where you must upload your raw DNA, Navigator does all the heavy lifting locally. This "edge-computing" approach means:

- **Privacy First:** Your raw genomic files (BAM/CRAM, chip raw data, etc.) never leave your machine.
- **Data Sovereignty:** You own your data. Only optional, anonymized summaries are shared if you choose to connect to the Decoding-Us Federation.
- **No external tooling:** Navigator is a single self-contained Rust application. There is **no Java runtime, no GATK, no samtools/bcftools** to install — the analysis engine ([noodles](https://github.com/zaeleus/noodles)) is built in. That keeps the download tiny: each installer is one file of roughly **60–85 MB** (Windows ≈ 62 MB, the Linux `.deb`/AppImage packages ≈ 74–78 MB, the universal macOS `.dmg` ≈ 83 MB because it bundles both Apple Silicon and Intel), and that single file *is* the whole application. A conventional stack has to install a Java runtime — which by itself is larger than any of these — and then GATK, samtools, and bcftools on top of it: a full JDK runs ~150–300 MB, the GATK distribution ~300–400 MB, and samtools/bcftools/HTSlib another few tens of MB, so the traditional toolchain typically lands somewhere between **500 MB and well over 1 GB** installed. Navigator does the same work in a single sub-100 MB download.
- **Accessibility:** Complex command-line bioinformatics is wrapped in an easy-to-use desktop interface, with an optional scriptable CLI for power users.

## System Requirements

- **Operating System:** macOS, Windows, or Linux.
- **Runtime:** None required. Navigator ships as a single native binary — no Java, no Python, no external bioinformatics tools.
- **Memory (RAM):**
  - Minimum: 4 GB
  - Recommended: 8 GB or more for large Whole Genome Sequencing files.
- **CPU:** Analysis is parallelized across cores; more cores means faster coverage and haplogroup calling.
- **Disk Space:** Room for your sequencing files plus cached reference genomes (roughly 5–10 GB per reference build).

## Installation & Setup

### Prebuilt installers (recommended)

For most people the simplest path is to grab a prebuilt installer from the [GitHub Releases page](https://github.com/JamesKane/decodingus-navigator/releases/latest) — download, install, and launch. There is a package for every common desktop:

| Platform | Package |
|----------|---------|
| **macOS** (Apple Silicon + Intel, universal) | signed, notarized `.dmg` |
| **Windows** (x64) | `.exe` setup installer |
| **Linux** (Debian/Ubuntu family, x86_64 / ARM64) | `.deb` |
| **Linux** (any distro, x86_64 / ARM64) | self-contained `.AppImage` |

On Linux, take the `.deb` if you are on a Debian/Ubuntu-family distribution and the AppImage if you would rather have a single self-contained executable. Each release also ships a `SHA256SUMS` file so you can verify your download. Because these are Alpha builds, newer tags land as bugs are fixed; the [latest release page](https://github.com/JamesKane/decodingus-navigator/releases/latest) always points at the freshest packages.

### Building from source

Because Navigator is one self-contained Rust binary with no external tools, building from source is genuinely easy — this is the path if you are on a platform without a prebuilt installer (FreeBSD, or a less common Linux setup), or if you simply prefer to build your own. Install [Cargo](https://www.rust-lang.org/tools/install) (the Rust toolchain) first if you don't have it, then:

```bash
git clone https://github.com/JamesKane/decodingus-navigator
cd decodingus-navigator

# Build the whole workspace (use --release for an optimized build)
cargo build --release

# Launch the desktop app
cargo run -p navigator-ui
```

The optimized binary is named `navigator` and lands at `target/release/navigator`. Once built, you can launch it directly:

```bash
./target/release/navigator
```

### Running it

Running `navigator` with no arguments opens the graphical Workbench. Running it with a subcommand (`ingest`, `subjects`, `show`, `projects`, `call`, `branch-report`, `liftvcf`) runs in headless mode against the same workspace — see [The Command Line](#the-command-line).

## Getting Started

### The Workbench
When you launch Navigator you land in the **Workbench**, organized around three top-level tabs:

- **Dashboard** — A high-level overview of your projects and subjects.
- **Subjects** — The master table of every research subject (biosample). Select a row to open its detail panel on the right.
- **Projects** — Your project groupings and their member counts.

The Subjects table shows each subject's ID, name, Y-DNA and mtDNA haplogroups, sex, originating center, and analysis status at a glance.

### First Launch
On first launch, Navigator creates its local workspace database automatically at `~/.decodingus/navigator-rs.db`. No manual configuration is required.

### First-Time Setup: Bringing Your Own Reference Genomes
By default Navigator downloads and caches the reference builds it needs (GRCh38, GRCh37, CHM13v2) on first use, so most people never touch a reference file. But if you already run a bioinformatics toolchain, you almost certainly have the **exact** reference FASTAs your alignments were built against. Pointing Navigator at those files instead of letting it download its own copy has three benefits: it guarantees the coordinate space matches your data bit-for-bit (same contig names, same sequence), it saves the download and the several GB of duplicate cache per build, and it lets you work fully offline.

Register your references once, before you start importing, from **⚙ Settings → Reference genomes**. That panel shows one row per build with these columns:

| Column | What it does |
|--------|--------------|
| **Build** | The build key Navigator resolves against: `GRCh38`, `GRCh37`, `chm13v2.0`, or `chm13v2.0_maskedY_rCRS`. |
| **Status** | Whether that build is currently cached, overridden, or absent. |
| **Local FASTA** | Path to *your* reference FASTA. Type it or use 📂 to browse. When set, Navigator uses this file as-is and never downloads that build. |
| **Auto-download** | Untick to forbid Navigator from ever fetching that build — useful when you want to guarantee only your file is used, or you are offline. |
| **Integrity** | **Verify** hashes the file (and checks it against a pinned SHA-256, if you set one). |

Requirements for a reference you supply:

- **It must be decompressed** (a plain `.fa` / `.fasta` / `.fna`), and it **must be `faidx`-indexed** — a `.fai` file has to sit next to it (e.g. `chm13v2.0.fa` + `chm13v2.0.fa.fai`). If you have the FASTA but not the index, create it with any faidx tool from your existing toolchain (`samtools faidx chm13v2.0.fa`).
- **Match the row to the build your alignments declare.** Point the `chm13v2.0` row at your CHM13v2 FASTA, the `GRCh38` row at your GRCh38 FASTA, and so on. Navigator picks the reference per alignment from the build it detects in each file's header, then uses the FASTA you registered for that build. (`chm13`, `chm13v2`, and `hs1` are all treated as the same `chm13v2.0` build.)
- **Contig names must agree** with your alignments. This is automatic when it is literally the same FASTA your aligner used — which is the whole point of bringing your own.

Overrides are stored in `~/.decodingus/config/reference_sources.json`, which you can also hand-edit (the Settings panel just writes this file). Each key is a build; per build you may set `local_path` (use this exact FASTA, never download), `url` (an alternate download mirror), `sha256` (a pinned integrity hash — a download that doesn't match is rejected), and `auto_download` (`false` to hard-forbid fetching that build):

```json
{
  "references": {
    "chm13v2.0": {
      "local_path": "/Volumes/Refs/chm13v2.0.fa",
      "auto_download": false
    },
    "GRCh38": {
      "url": "https://my-mirror.example/GRCh38.fa",
      "sha256": "…",
      "auto_download": true
    }
  }
}
```

These overrides are global (per build, applied to every alignment and analysis on that build), so registering them once at first-time setup covers every subject you import afterward.

## Core Features

### Workspace Management
Organize your research:
- **Subjects (Biosamples):** Create an entry for each individual you study. The subject detail panel has sub-tabs for:
  - **Overview** — identity, summary status, and consensus haplogroup assignments.
  - **Y-DNA** — split into **Haplogroup** (placement and supporting branch evidence), **SNP** (the full genotyped-variant table, including **Private** off-backbone calls and **Imported** vendor Y-SNPs), and **STR** (Y-STR panel reports).
  - **mtDNA** — **Summary** (maternal haplogroup consensus) and **Variants** (rCRS-relative mutation list and heteroplasmy).
  - **Autosomal** — **Summary** plus a **Profile** diploid genotype table from the SNV/indel caller.
  - **Ancestry** — admixture, PCA, fine-population breakdown, and DNA painting.
  - **IBD Matches** — shared-segment detection and network match suggestions.
  - **Sources** — the per-result hub where you add files and see every run, alignment, and profile attached to the subject.
- **Projects:** Group related subjects (e.g. "Family Study", "Ancient DNA") and assign an administrator.

### Importing Data
Navigator auto-detects the type of any file you import and routes it appropriately. Supported sources:

| Source | What it is |
|--------|-----------|
| **BAM / CRAM** | Aligned sequencing reads (attached to a sequencing run). |
| **VCF / GVCF** | Variant calls from any caller; GVCF additionally carries callable-region context for a fast haplogroup path. |
| **mtDNA FASTA** | Mitochondrial sequence (`.fasta`/`.fa`/`.fna`/`.fas`, plain or `.gz`) for maternal-lineage assignment. |
| **Chip / array raw data** | Consumer genotype files from 23andMe, AncestryDNA, MyHeritage, and Living DNA. Y and mtDNA haplogroups (and autosomal ancestry) are placed on import. |
| **Y-STR profiles** | Short-tandem-repeat CSV/TSV exports (e.g. FTDNA / YSEQ), marker name + repeat count. |
| **Y-SNP panels** | BISDNA chromo2 genotyped Y-SNP exports, imported as real variant calls. |

To import in the desktop app: select a subject, open the **Sources** tab, and add a file. Navigator computes a checksum, detects the platform/test type, and files the data under the right run, alignment, or profile.

### Project Import (batch, with the sidecar fast path)
When you have many samples to load — for example a whole sequencing project staged on a NAS — use **Project Import** in the desktop app to ingest an entire directory tree in one pass. Navigator scans the folder, creates the project and one subject per sample, and attaches each sample's files.

#### Expected directory layout
Project Import expects a **two-level** layout: a project folder whose immediate subfolders are each one sample.

```
MyProject/                              ← project (named after this folder)
├── HG00096/                            ← one subject (named after this folder)
│   ├── HG00096.chm13.cram              ← alignment (+ HG00096.chm13.cram.crai)
│   ├── HG00096.chm13.chrY.g.vcf.gz     ← Y sidecar (+ .tbi)
│   ├── HG00096.chm13.chrM.g.vcf.gz     ← mtDNA sidecar (+ .tbi)
│   ├── HG00096.chm13.chrYM.callable.summary.txt
│   ├── HG00096.chm13.sex
│   ├── coverage.txt
│   └── stats.txt
├── HG00097/
│   └── ...
```

- The **project name** is the top folder's name; **each immediate subfolder is one subject**, named after the folder.
- Files inside a sample folder are found up to two levels deep. Hidden (dot) folders are skipped, and a subfolder with no alignment or variant file is ignored.

#### The sidecar "hot path"
Walking a 10–12 GB CRAM to place a haplogroup takes many minutes. If the pipeline that produced the alignment also left its per-sample intermediate ("sidecar") files **next to the CRAM**, Navigator reads those instead of touching the CRAM — turning per-sample placement from minutes into seconds (HG00096 places to R1b1a1b1a1a in ~5 s versus a ~22-minute CRAM walk). The fast path is **on by default**; it runs during import and returns quickly.

Recognized sidecars (matched by file-name suffix, case-insensitive — the sample-name prefix can be anything):

| Sidecar file | What it provides | Completeness |
|--------------|------------------|--------------|
| `*.chrY.g.vcf.gz` (+ `.tbi`) | Y-DNA haplogroup | Full |
| `*.chrM.g.vcf.gz` (+ `.tbi`) | mtDNA haplogroup | Full |
| `*.sex` (contains `male`/`female`) | Genetic sex | Full |
| `stats.txt` (`samtools stats`) / `*.flagstat` | Read metrics (counts, mean read length, insert size) | Full |
| `coverage.txt` (`samtools coverage`) + `*.callable.summary.txt` | Coverage roll-up (genome-wide mean depth, per-contig stats, callable bases) | Partial ("lite") |
| `*wgs*metric*` (Picard `CollectWgsMetrics`) | Genome-wide depth distribution | Supplemental |
| `*alignment_summary*` (Picard `CollectAlignmentSummaryMetrics`) | Read metrics | Supplemental |

Notes and requirements:
- **GVCFs must be ploidy-1 (haploid) chrY/chrM GVCFs**, and the matching `.tbi` tabix index must sit beside each one so Navigator can read just the needed positions.
- **The build must match.** Navigator reads the build token from the GVCF file name (e.g. the `chm13` in `HG00096.chm13.chrY.g.vcf.gz`) and only takes the fast path when it matches the alignment's reference build. `chm13`, `chm13v2`, and `hs1` are treated as the same build. If the builds differ, Navigator falls back to walking the CRAM (it will not lift GVCF coordinates).
- **A reference genome is still required.** Even on the fast path, Navigator reads the reference FASTA at the relevant positions. Let Navigator resolve/download the reference from the detected build, or point it at an explicit FASTA — which must have its `.fai` index alongside.
- `coverage.txt` and `stats.txt` are matched by exact name; the GVCF/`.sex`/`.callable.summary` files are matched by suffix.

The lite coverage roll-up is the only **partial** result: median depth, the `pct_Nx` thresholds, and the full depth histogram are not in `coverage.txt` and are filled in later by deep analysis.

#### What the fast path does *not* cover
Some analyses always need the CRAM and are **not** produced from sidecars: autosomal **ancestry**, the **full coverage histogram** (median, `pct_10x`/`pct_20x`, depth distribution), **structural variants**, the **diploid SNV/indel caller**, and **IBD** panel genotyping. These run only when you trigger **deep analysis** — use **Analyze All** on the project (or run analysis on a subject). Deep analysis is additive: haplogroups, sex, and read metrics already placed by the fast path are **not** recomputed, and the lite coverage is upgraded in place to the full result.

> **Where to find it:** Project Import and the sidecar fast path are available in the **desktop app**. The headless `navigator ingest` command imports files and directories via auto-detection; a directory argument is treated as one staged sample, so the sidecar fast path applies to it too.

### Batch import strategies for existing data collections
Real-world collections come in two shapes, and each has its own best path in. The dividing question is simple: **is the on-disk layout already `{project}/{sample}/files…`, with folder names you're happy to use as subject names?** If yes, use the desktop **Project Import** directly. If the layout is deeper, uses opaque identifiers, or keeps its human-readable names in a separate manifest, script the CLI instead.

#### Strategy A — a clean project tree (use Project Import as-is)
This is the layout Project Import was built for. For example, a PGP-style collection where each sample is a top-level folder named for the donor:

```
PGP_Harvard/                                     ← project
├── hu46DD40/                                    ← subject (named "hu46DD40")
│   ├── hu46DD40.chm13_HG002Y.cram (+ .crai)     ← alignment
│   ├── hu46DD40.chm13.chrY.g.vcf.gz (+ .tbi)    ← Y sidecar
│   ├── hu46DD40.chm13.chrM.g.vcf.gz (+ .tbi)    ← mtDNA sidecar
│   ├── hu46DD40.chm13.chrYM.callable.summary.txt
│   ├── hu46DD40.chm13.sex
│   ├── coverage.txt
│   └── stats.txt
├── hu0F18A8/
└── …
```

Point Project Import at `/Volumes/Genomics/PGP_Harvard`, leave the fast path on, and go. Each `hu…` folder becomes a subject named after the folder, the sidecars place Y/mtDNA/sex/read-metrics in seconds, and you run **Analyze All** afterward for the deep results (ancestry, full coverage, SV, diploid calls, IBD). No scripting required.

#### Strategy B — a deep tree with an external map (script the CLI)
Pipelines that key everything by UUID and record the human-readable identity in a side manifest do **not** fit the two-level scanner. A D2C-style repository is the canonical example:

```
D2C/
├── _manifests/
│   └── biosample_map.tsv                         ← subject → name/lab/kit + file paths
├── 0a0e8267-dc23-4be4-b86f-4190e59de02b/         ← biosample (opaque UUID)
│   └── 1aceb711-b601-44f5-8835-b361aa95c6e3/     ← analysis run (UUID)
│       ├── b38/          chrYM.cram, gatk3/, gatk4/, coverage.txt, stats.txt
│       └── CP086569.2/   chrYM.cram, gatk3/, gatk4/, coverage.txt, stats.txt
└── …
```

Handing this tree to Project Import goes wrong in three ways: subjects would be named by **opaque UUIDs** instead of the friendly `Dante-23823` names; the **lab and kit** metadata lives only in `biosample_map.tsv`, which the scanner does not read; and each biosample holds **multiple reference builds** (`b38` and `CP086569.2`) with the callable BEDs a directory deeper (`…/CP086569.2/gatk3/callable_status.bed`) than the scanner descends. So the manifest — not the directory names — is the source of truth, and you drive the import from it.

The map has one row per subject, tab-separated, with the columns Navigator cares about:

```
subject   name             lab    kit           y_tier         y_artifact   cram   callable   coverage   stats
```

`name` is the friendly subject label; `cram`/`callable`/`coverage`/`stats` are absolute paths **as the producing pipeline saw them** (e.g. `/mnt/md0/Repo/…`), so on your machine you remap that prefix onto your local mount. Loop the rows and call `navigator ingest` once per subject, taking `name` for `--subject` and pointing at the one reference directory you want per run:

```bash
#!/usr/bin/env bash
set -euo pipefail
ROOT=/Volumes/Genomics/D2C            # local mount
NAV=./target/release/navigator

# Skip the header row; read only the columns we use.
tail -n +2 "$ROOT/_manifests/biosample_map.tsv" |
while IFS=$'\t' read -r subject name lab kit y_tier y_artifact cram callable coverage stats; do
  # Remap the pipeline path (/mnt/md0/Repo/…) onto the local mount, then take its directory —
  # that folder (…/CP086569.2) holds the CRAM plus its coverage.txt / stats.txt sidecars.
  local_cram="${cram/\/mnt\/md0\/Repo/$ROOT}"
  sample_dir="$(dirname "$local_cram")"
  [ -d "$sample_dir" ] || { echo "skip $name — $sample_dir missing"; continue; }

  "$NAV" ingest \
    --subject "$name" \
    --project "D2C" \
    --test-type "Big Y" \
    "$sample_dir"
done
```

Notes on this pattern:

- **One reference build per run.** Point each `ingest` at a specific reference subdirectory (`CP086569.2` for the Y/T2T build, or `b38` for GRCh38) rather than the biosample root, so you don't fold two builds into a single sequencing run. Run the loop twice against different subdirs if you want both.
- **A directory argument is one staged sample**, so the fast path applies: the CRAM's neighboring `coverage.txt` / `stats.txt` are picked up automatically. If a per-subject artifact lives in a *different* tree (the map's `y_artifact` column often points a Y GVCF at a separate `…/ytree/flat/…` path), add its remapped path as an extra argument to the same `ingest` call — `ingest` accepts multiple files and directories at once.
- **`--test-type`** forces the sequencing-run type when the folder layout tells you what it is (these `chrYM.cram` files are Y-focused), which is more reliable than letting a CRAM without a `.bai` fall back to generic WGS.
- **Idempotent and resumable.** `ingest` finds-or-creates the subject and project and skips already-imported paths, so you can re-run the loop after adding kits, or after fixing a few `skip` lines, without creating duplicates.
- The map's `lab`/`kit` columns aren't consumed by `ingest` directly; sequencing-lab and instrument are inferred from each alignment's header during analysis. Use `name` for the subject label, and keep the map alongside your collection as the record of provenance.

After the loop, run deep analysis (**Analyze All** on the `D2C` project in the desktop app, or `navigator` analysis per subject) to add everything beyond the fast-path haplogroups.

### Importing an FTDNA group project
If you administer a FamilyTreeDNA **group project** (a surname or haplogroup project), Navigator can ingest the project's roster, genealogy, and Y-STR chart in one pass. This is a different importer from [Project Import](#project-import-batch-with-the-sidecar-fast-path) above: Project Import walks a folder of *sequencing files*; this reads the four **CSV exports** that FamilyTreeDNA's Group Administration Pages (GAP) produce. It creates one subject per kit, records each member's paternal/maternal most-distant-known-ancestor and vendor kit number, attaches the Y-STR panel from the results chart, and files everyone into the project — without any BAM/CRAM.

> These CSVs are the **administrator** exports. Only a project's admin or co-admin can download them, from the project's GAP pages. This importer is for running your own project's data; it is not a way to pull another project's members.

#### Recommended structure
Keep one project's four exports together in a single folder named for the project. Downloading all four "Download to Excel" exports from GAP gives you exactly these files, each already prefixed with the project name:

```
R1b-CTS4466Plus/                                       ← one folder per project
├── R1b-CTS4466Plus_Member_Information_20260619.csv    ← roster (kits, names, consent flags)
├── R1b-CTS4466Plus_Paternal_Ancestry_20260619.csv     ← paternal ancestor + Y clade subgroup
├── R1b-CTS4466Plus_Maternal_Ancestry_20260619.csv     ← maternal ancestor + mtDNA subgroup
└── R1b-CTS4466Plus_YDNA_Results_Overview.csv          ← wide Y-STR marker chart (DYS…)
```

| Export | What it contributes |
|--------|---------------------|
| **Member Information** | The roster: kit number, member name, and the FTDNA consent flags (`Access Granted`, `Publicly Share DNA Results`). Provides the kit → identity spine. |
| **Paternal Ancestry** | Each kit's paternal most-distant-known ancestor (name, place, country, map coordinates) plus the paternal-clade **Sub Group** path, which supplies a provisional Y terminal and the project subgroup label. |
| **Maternal Ancestry** | The maternal most-distant-known ancestor and the mtDNA subgroup, in the same layout. |
| **YDNA Results Overview** | The wide Y-STR marker table (DYS-prefixed columns). Attaches a Y-STR panel profile (Y-12 … Y-700, sized to the populated markers) to each kit. |

All four are optional — a roster-only or ancestry-only import is valid — but the full set gives the richest result. Files are recognized by their **header content, not their names**, so a renamed export still routes correctly; the filename's project-name prefix is used only to name/target the project.

#### How to import
In the desktop app's **Projects** area, use **Import FTDNA project** and select the CSVs together in the file picker (pick all four at once). Navigator then:

1. **Classifies** each file, joins all rows **by kit number**, and matches every kit against your existing workspace.
2. Shows a **dry-run plan** — nothing is written yet. Each kit is marked **New** (create a subject), **Auto-merge** (an exact FTDNA kit number already in the workspace — locked, always reused), or **Needs confirm** (a fuzzy candidate matched on shared Y-terminal SNP, near-zero Y-STR genetic distance, or overlapping names — you confirm or reject each).
3. On **commit**, applies your resolutions. For each kit it creates or reuses a subject, attaches the FTDNA kit number as a vendor id, stores the member name and paternal/maternal ancestor (MDKA) rows, adds project membership tagged with the clade subgroup, and — for newly created subjects — saves the Y-STR profile. An unresolved fuzzy row defaults to **New**, so it never silently merges.

The project name comes from the export filename prefix (`R1b-CTS4466Plus`): if a project of that name is already open or exists, the kits go into it; otherwise Navigator creates it on commit. Re-running the import later is safe — kits already imported under their FTDNA kit number auto-merge rather than duplicating.

> **What this does *not* import:** sequencing reads or variant calls. It brings in roster, genealogy, and Y-STR only. To add a member's BAM/CRAM, Big Y variant CSV, or VCF, open that subject's **Sources** tab (or use Project Import) and add the file there; it attaches to the same subject the group import created.

### Running Analyses
Open a subject's detail panel and run any module from the relevant tab, or use **Full Analyze** to run a complete pass over all of a subject's data. Results are cached, so re-running is instant when nothing has changed.

Available analyses:

| Analysis | Status | What it gives you |
|----------|--------|-------------------|
| **Coverage / Callable Loci** | Validated | Mean depth, coverage distribution, per-contig depth histograms, and which bases are callable per contig (1×–100× thresholds). |
| **Read Metrics** | Validated | Read-length and insert-size distributions, platform/instrument detection, library orientation, and sequencing-lab inference. |
| **Sex Inference** | Validated | Inferred genetic sex with a confidence score. |
| **Y-DNA Haplogroup** | Validated | Terminal haplogroup plus ranked candidates and supporting branch evidence. Handles GRCh37/GRCh38/CHM13v2 coordinates automatically, against either the DecodingUs or FTDNA tree. |
| **Y-STR Profiles** | Validated | FTDNA/YSEQ-style panel tables (Y-12 … Y-111, YSEQ tiers) with per-marker consensus and conflict detection across sources. |
| **mtDNA Haplogroup** | Validated | Terminal maternal haplogroup from sequence or alignment, with rCRS↔CHM13 mapping. |
| **mtDNA Variants & Heteroplasmy** | Validated (variants); screening (heteroplasmy) | rCRS-relative mutation list (HVR1/HVR2/coding) plus site-level heteroplasmy. Heteroplasmy is a screening pass, not a clinical caller. |
| **Private Y Variants** | Validated | Off-backbone calls — finer branches and novel candidate variants, reconciled across sources. |
| **Ancestry** | Validated | Admixture across fine populations / continental groups (ADMIXTURE, PCA projection + GMM, and an nMonte/G25-style estimate), a geographic map, fine-population breakdown, and DNA-painting local ancestry. |
| **Diploid Variant Calling** | Validated on test data | De-novo **diploid** SNV + indel calls, exportable as a whole-genome VCF (per subject or per alignment). |
| **IBD Detection** | Validated (detection) | Pairwise shared-segment detection and relationship estimates, using a real recombination map. The match-discovery / network UI is still in progress. |
| **Structural Variants (SV)** | Built, output unvalidated | Deletions, duplications, inversions, and breakends. Reliable output needs ≥10× coverage. |

Navigator also reconciles Y/mtDNA haplogroups across multiple runs and alignments per subject into a single genome-level **consensus** assignment, rather than voting on per-run labels.

### The Branch Report tool
The **Branch Report** answers a narrow, practical question: *for an arbitrary branch of the tree, how does this sample genotype at every marker that defines it and its descendants?* You give it any Y or mtDNA node — not just the one the sample was placed on — and it genotypes that node's whole **descendant subtree** fresh, marker by marker, showing the observed base, the derived/ancestral call, and the supporting read evidence for each.

That "any node, subtree-wide" behavior is what makes it a checking tool rather than a placement view. The normal haplogroup card walks the sample's *assigned* path from root to terminal. The Branch Report instead genotypes the subtree you name, so **sibling branches the sample is *ancestral* for are reported too** — which is exactly what you need to confirm a variant sits where it should. Point two researchers' samples at the same parent node and you can see, side by side, that the SNP defining one sibling branch is derived in the sample that belongs there and ancestral (absent) in the one that doesn't. If a variant were mis-mapped or placed on the wrong branch, the two reports would disagree at that marker, and you would catch it before it propagated into the shared tree.

**Where to find it.** Open a subject's detail panel and go to the **Y-DNA** tab (for the Y tree) or the **mtDNA** tab (for the mtDNA tree). The Branch Report card has a node text box and a **Load** button. Type a node and load it:

- The node can be a **haplogroup name** (`R-M269`, `R-FGC29071`, `H2a`) or a **defining marker** (`FGC29071`) — either resolves to the same subtree.
- Loading a **shallow** node (say `R-M269`, or the tree root) pulls in tens of thousands of markers, so it can take a moment; a terminal or near-terminal branch is near-instant. There is an optional depth limit (see the CLI below) if you only want the top few levels.

**What each row shows.** One row per defining marker in the subtree, columns: `node` / `parent` (where the marker sits on the tree), `marker`, `pos`, `anc>der` (ancestral→derived alleles), `obs` (the observed base), `status` (**derived** = the sample carries it, **ancestral** = it doesn't, **no-call** = no confident base), then `AD` / `DP` / `GQ` read evidence and a `note` (flags like `indel/MNV`, `hom-ref block`, or `no call`). The card header summarizes the tally — *N markers: d derived / a ancestral / n no-call* — and whether the evidence came from a **gVCF** sidecar (rich DP/AD/GQ) or a live **pileup**.

**Reading it — a worked example.** Here is the TSV a subject placed at `R-FGC29071` produces when you query that node (evidence columns shown as `.` here because this run came from a pileup rather than a gVCF sidecar):

```
# DUNavigator Y-DNA branch report — node R-FGC29071 (chrY); 4 derived / 2 ancestral / 2 no-call
node        parent      marker              chrom  pos       ancestral  derived  observed_base  status     GT  AD  DP  GQ  source  note
R-FGC29071  R-FGC29067  FGC29071            chrY   15570629  A          C        C              derived    1   .   .   .   pileup
R-FGC29071  R-FGC29067  FGC29076            chrY   20512639  G          T        T              derived    1   .   .   .   pileup
R-FGC29071  R-FGC29067  chrY:14583465G>T    chrY   14583465  G          T        T              derived    1   .   .   .   pileup
R-FGC29071  R-FGC29067  chrY:3332132A>T     chrY   3332132   A          T        T              derived    1   .   .   .   pileup
R-MF41134   R-FGC29071  BY74966             chrY   8442212   T          G        .              nocall     .   .   .   .   pileup  no call
R-MF41134   R-FGC29071  chrY:12803849C>T    chrY   12803849  C          T        .              nocall     .   .   .   .   pileup  no call
R-MF41134   R-FGC29071  chrY:3464631C>T     chrY   3464631   C          T        C              ancestral  0   .   .   .   pileup
R-Y178014   R-MF41134   chrY:11687241T>C    chrY   11687241  T          C        T              ancestral  0   .   .   .   pileup
```

Read top to bottom it tells a clear story: the four markers that define `R-FGC29071` itself are all **derived** (the sample observes the derived base — `C`, `T`, `T`, `T`), which is what puts the sample on this branch. The rows below belong to the **descendant** subtree — the child branch `R-MF41134` and its child `R-Y178014` — and there the sample is **ancestral** or **no-call**, meaning it does *not* descend into them. That contrast is the whole point: it confirms the placement terminates at `R-FGC29071` and does not belong on the deeper branches. If a collaborator's sample were a true match on a deeper branch, their report at the same node would show those `R-MF41134` markers flipping to **derived** instead — and if a variant were mis-mapped, the two reports would disagree at exactly that row.

**Sharing it.** The **Export** button writes this TSV (the `GT` column is VCF-style: `1` derived, `0` ancestral, `.` no-call), which is the format to hand another researcher when you are cross-checking placements between labs — they load the same node on their own sample and diff the two files marker for marker.

### Exporting & Sharing Results
Result cards carry an **Export** action that writes a shareable file via a save dialog. Available formats:

| Result | Formats |
|--------|---------|
| Coverage | TSV, self-contained HTML |
| Read metrics | TSV |
| Ancestry | TSV, self-contained HTML |
| mtDNA variants | TSV |
| Callable loci | BED4 (0-based, half-open) |
| IBD segments | TSV |
| Diploid variants | VCF (per alignment, or a subject-level consensus across same-build alignments) |

The same diploid VCF export is also available headlessly via the [`call`](#the-command-line) subcommand.

## The Command Line

The same `navigator` binary is fully scriptable. With a subcommand it opens the *same* workspace database as the GUI, does its work, and exits. This is ideal for bulk-loading a directory of files, querying results, or producing VCFs in a pipeline.

```bash
# Import everything in a folder into a subject (creating the subject/project if needed)
navigator ingest --subject "Jane Doe" --project "Family Study" --recursive /Volumes/nas/Genomics/jane/

# List all subjects with their data-source counts
navigator subjects

# Show one subject's runs, alignments, profiles, and haplogroup consensus
navigator show --subject "Jane Doe"

# List projects with subject counts
navigator projects

# Call de-novo diploid SNVs + indels to a VCF (whole genome, or one contig)
navigator call --subject "Jane Doe" --out jane.vcf
navigator call --subject "Jane Doe" --contig chr21 --out jane.chr21.vcf

# Branch report: genotype a subject at every defining marker of a Y (or mtDNA) node's subtree
navigator branch-report --subject "Jane Doe" --node R-FGC29071 --tree y
navigator branch-report --subject "Jane Doe" --node H2a --tree mt --tsv jane.mt.branch.tsv

# Lift a VCF from one reference build to another
navigator liftvcf --in calls.GRCh38.vcf.gz --from GRCh38 --to chm13v2.0 --out calls.chm13.vcf.gz
```

Useful flags:
- `--subject` / `-s` — donor identifier (found by exact match, or created on `ingest`).
- `--project` / `-p` — project to assign the subject to (found or created).
- `--sex` — recorded only when a new subject is created (e.g. `male` / `female`).
- `--recursive` / `-r` — recurse into directories instead of importing only their immediate files.
- `--alignment` / `-a` — (for `call` / `branch-report`) target a specific alignment id from `show --json`; omit to use the subject's sole alignment (`branch-report` prefers a CHM13/HiFi alignment when the subject has several).
- `--contig` / `-c` — (for `call`) restrict to a single contig (e.g. `chrM`, `chr21`); default is every primary chromosome.
- `--node` / `-n`, `--tree` / `-t`, `--depth` — (for `branch-report`) the node to report (a haplogroup name like `R-FGC29071` or a defining marker like `FGC29071`), which tree to read (`y` or `mt`), and an optional cap on how many levels below the node to descend (default: the whole subtree).
- `--tsv` — (for `branch-report`) write the report as TSV to a file instead of printing a table; `--json` emits JSON instead (the two are mutually exclusive).
- `--out` / `-o` — (for `call` / `liftvcf`) write the VCF to a file instead of stdout.
- `--in` / `-i`, `--to` / `-t`, `--from` / `-f`, `--filter-par` — (for `liftvcf`) input VCF, target build, optional source build (inferred from the header when omitted), and whether to drop variants landing in the target chrY PAR.
- `--db` — point at an alternate workspace database (defaults to `~/.decodingus/navigator-rs.db`).
- `--json` — emit machine-readable JSON instead of a table (on `subjects`, `show`, `projects`, `branch-report`).

If you're running from source without an installed binary, prefix with `cargo run -p navigator-ui --`:

```bash
cargo run -p navigator-ui -- subjects --json
```

## Data Management & Privacy

### Where is my data?
All application data lives under your home directory in `~/.decodingus/`:

```
~/.decodingus/
├── navigator-rs.db      # Workspace database (SQLite): subjects, projects, runs, alignments, profiles
├── config/              # settings.json (your saved preferences)
├── references/          # Downloaded reference genomes (indexed FASTA)
├── liftover/            # Chain files for build-to-build coordinate conversion
├── masks/               # Callable-region BED masks
├── trees/               # Cached Y-DNA / mtDNA haplotrees (JSON)
├── ysnp/                # Y-SNP dictionary assets
└── ancestry/            # Pre-built ancestry panels and PCA loadings
```

### Reference Genomes
Navigator manages reference genomes for you. It downloads and caches standard builds (GRCh38, GRCh37, CHM13v2) as needed and converts coordinates between builds automatically — you don't need to hunt for reference files. If you already have the exact FASTAs from your own toolchain, you can register them so Navigator uses yours instead of downloading — see [First-Time Setup: Bringing Your Own Reference Genomes](#first-time-setup-bringing-your-own-reference-genomes).

### Cloud Integration (Optional)
Navigator includes support for the **AT Protocol** — the same federated network behind [Bluesky](https://bsky.app) — to publish summaries to a Personal Data Store (PDS) in the Decoding-Us Federation. Everything else in Navigator works fully offline; contributing your results back is opt-in, and it's how the shared, community-built haplogroup tree grows denser.

To contribute, you sign in with AT Protocol credentials, and Navigator publishes your *results* (haplogroup placements and variant observations, not your raw genome) to your own data store on the network.

- **Privacy:** Even with publishing enabled, your raw BAM/CRAM and chip files are **never** uploaded. Only anonymized summaries (haplogroup assignments, coverage QC statistics, ancestry estimates, IBD attestations) are shared, with your explicit consent.
- Publishing is durable: queued summaries are retried with backoff if the network or PDS is briefly unavailable.
- Configure the AppView endpoint in [Settings](#settings) or via the `DECODINGUS_APPVIEW_URL` environment variable.

Two recommendations for signing in comfortably:

- **Use a dedicated profile, not your main Bluesky account.** Make a separate handle for your genomics contributions and sign Navigator in with that. It keeps your genealogy activity cleanly separated from your personal social account, and if you ever want to hand off or retire the contributing identity, you can do it without touching your everyday presence.
- **A private PDS is nice to have, not required.** In AT Protocol terms your data lives in a Personal Data Store. Running your own PDS gives you the fullest ownership, but self-hosting one is genuinely a homelab project. If that's not your thing, use a hosted PDS (the default Bluesky one is fine) and you still keep control of your records and can move them later. Self-hosting is the enthusiast option, not the price of admission.

If you never sign in at all, Navigator remains a complete local analysis tool — contributing is a choice, not a toggle you have to flip to get value.

## Settings

Open the **⚙ Settings** dialog from the app bar to configure (saved to `~/.decodingus/config/settings.json`; environment variables take precedence over saved settings):

- **Connection** — the Federation **AppView URL** for haplotree updates and publishing.
- **Appearance** — light/dark **theme** and **UI scale**.
- **Reference** — the reference-genome cache directory and whether to **prompt before downloading** large reference files.
- **Advanced** — the **Y-tree provider** (`decodingus` or `ftdna`) and the haplotree cache **TTL** (days before refetch; `0` = always refetch).

## Advanced Usage

### Direct Database Access
The workspace is a standard SQLite database at `~/.decodingus/navigator-rs.db`. Power users can query it with any SQLite tool (the `sqlite3` CLI, DB Browser for SQLite, DBeaver, etc.):

```bash
sqlite3 ~/.decodingus/navigator-rs.db ".tables"
```

Close the Navigator app first to avoid write contention.

### Environment Variables
Tune behavior without changing code:

| Variable | Purpose | Default |
|----------|---------|---------|
| `NAVIGATOR_ANALYSIS_THREADS` | Worker threads for per-contig analysis fan-out. | Auto (based on cores) |
| `NAVIGATOR_BGZF_THREADS` | BGZF decompression workers for BAM/CRAM reading. | Auto |
| `NAVIGATOR_Y_TREE_PROVIDER` | Y-tree source: `decodingus` or `ftdna`. | `decodingus` |
| `NAVIGATOR_TREE_TTL_DAYS` | Days to cache haplotrees before refetching (0 = always refetch). | `7` |
| `NAVIGATOR_REFGENOME_DIR` | Root directory for reference/liftover/mask caches. | `~/.decodingus` |
| `NAVIGATOR_TREE_DIR` | Haplotree cache directory. | `~/.decodingus/trees` |
| `NAVIGATOR_ANCESTRY_PANEL` / `NAVIGATOR_ANCESTRY_PCA` | Override paths to pre-built ancestry assets. | `~/.decodingus/ancestry/...` |
| `DECODINGUS_APPVIEW_URL` | Federation AppView endpoint (haplotree updates + publishing). | `http://localhost:9000` |

## Troubleshooting

**Q: Analysis is slow.**
A: WGS analysis is computationally intensive. Navigator parallelizes across CPU cores automatically; you can cap or raise the worker count with `NAVIGATOR_ANALYSIS_THREADS`. A `--release` build is significantly faster than a debug build.

**Q: I can't find my reference genome.**
A: Navigator downloads references on demand. If you are offline, run an analysis at least once while online to cache the necessary files.

**Q: A haplogroup result looks out of date or under-placed.**
A: Haplotrees are cached for `NAVIGATOR_TREE_TTL_DAYS` (default 7). Lower that value or set it to `0` in [Settings](#settings) (or via the environment) to force a fresh fetch, then re-run the analysis.

**Q: My file wasn't recognized on import.**
A: Navigator auto-detects by extension and content fingerprint. Confirm the file is one of the [supported formats](#importing-data). Consumer chip exports from less common vendors may not be detected; the file is still recorded but won't be analyzed.

**Q: A sample imported with only haplogroups and basic metrics.**
A: That's the project-import [fast path](#the-sidecar-hot-path) using sidecar files. Run **Analyze All** (or analyze the subject) to add ancestry, the full coverage histogram, structural variants, the diploid caller, and IBD genotyping from the alignment itself.
