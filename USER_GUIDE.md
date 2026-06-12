# Decoding-Us Navigator User Guide

Welcome to the **Decoding-Us Navigator**, your private, local companion for advanced genomic analysis. This application lets you explore your DNA data with professional-grade bioinformatics directly on your own computer, keeping your genetic privacy intact while contributing to citizen science.

## Table of Contents
1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Installation & Setup](#installation--setup)
4. [Getting Started](#getting-started)
5. [Core Features](#core-features)
   - [Workspace Management](#workspace-management)
   - [Importing Data](#importing-data)
   - [Running Analyses](#running-analyses)
6. [The Command Line](#the-command-line)
7. [Data Management & Privacy](#data-management--privacy)
8. [Advanced Usage](#advanced-usage)
9. [Troubleshooting](#troubleshooting)

---

## Introduction

Decoding-Us Navigator runs a complete bioinformatics stack on your desktop. Unlike cloud services where you must upload your raw DNA, Navigator does all the heavy lifting locally. This "edge-computing" approach means:

- **Privacy First:** Your raw genomic files (BAM/CRAM, chip raw data, etc.) never leave your machine.
- **Data Sovereignty:** You own your data. Only optional, anonymized summaries are shared if you choose to connect to the Decoding-Us Federation.
- **No external tooling:** Navigator is a single self-contained Rust application. There is **no Java runtime, no GATK, no samtools/bcftools** to install — the analysis engine ([noodles](https://github.com/zaeleus/noodles)) is built in.
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

Navigator is built from source with [Cargo](https://www.rust-lang.org/tools/install) (the Rust toolchain). Install Rust first if you don't have it, then:

### Build and run the desktop app

```bash
# Build the whole workspace (use --release for an optimized build)
cargo build --release

# Launch the desktop app
cargo run -p navigator-ui
```

The optimized binary is named `navigator` and lands at `target/release/navigator`. Once built, you can launch it directly:

```bash
./target/release/navigator
```

Running `navigator` with no arguments opens the graphical Workbench. Running it with a subcommand (`ingest`, `subjects`, `show`, `projects`) runs in headless mode against the same workspace — see [The Command Line](#the-command-line).

## Getting Started

### The Workbench
When you launch Navigator you land in the **Workbench**, organized around three top-level tabs:

- **Dashboard** — A high-level overview of your projects and subjects.
- **Subjects** — The master table of every research subject (biosample). Select a row to open its detail panel on the right.
- **Projects** — Your project groupings and their member counts.

The Subjects table shows each subject's ID, name, Y-DNA and mtDNA haplogroups, sex, originating center, and analysis status at a glance.

### First Launch
On first launch, Navigator creates its local workspace database automatically at `~/.decodingus/navigator-rs.db`. No manual configuration is required.

## Core Features

### Workspace Management
Organize your research:
- **Subjects (Biosamples):** Create an entry for each individual you study. The subject detail panel has sub-tabs for **Overview**, **Y-DNA**, **mtDNA**, **Ancestry**, **IBD Matches**, and **Data Sources**.
- **Projects:** Group related subjects (e.g. "Family Study", "Ancient DNA") and assign an administrator.

### Importing Data
Navigator auto-detects the type of any file you import and routes it appropriately. Supported sources:

| Source | What it is |
|--------|-----------|
| **BAM / CRAM** | Aligned sequencing reads (attached to a sequencing run). |
| **VCF / GVCF** | Variant calls from any caller; GVCF additionally carries callable-region context for a fast haplogroup path. |
| **mtDNA FASTA** | Mitochondrial sequence (`.fasta`/`.fa`/`.fna`, plain or `.gz`) for maternal-lineage assignment. |
| **Chip / array raw data** | Consumer genotype files from 23andMe, AncestryDNA, MyHeritage, Living DNA, and FTDNA. Y and mtDNA haplogroups are placed on import. |
| **Y-STR profiles** | Short-tandem-repeat CSV/TSV exports (e.g. FTDNA/YSEQ), marker name + repeat count. |
| **Y-SNP panels** | BISDNA chromo2 genotyped Y-SNP exports, imported as real variant calls. |

To import in the desktop app: select a subject, open the **Data Sources** tab, and add a file. Navigator computes a checksum, detects the platform/test type, and files the data under the right run, alignment, or profile.

### Running Analyses
Open a subject's detail panel and run any module from the relevant tab, or use **Full Analyze** to run a complete pass over all of a subject's data. Results are cached, so re-running is instant when nothing has changed.

Available analyses:

| Analysis | Status | What it gives you |
|----------|--------|-------------------|
| **Coverage / Callable Loci** | Validated | Mean depth, coverage distribution, and which bases are callable per contig (1×–100× thresholds). |
| **Read Metrics** | Validated | Read-length and insert-size distributions, platform detection, library orientation. |
| **Sex Inference** | Validated | Inferred genetic sex with a confidence score. |
| **Y-DNA Haplogroup** | Validated | Terminal haplogroup plus ranked candidates and supporting branch evidence. Handles GRCh37/GRCh38/CHM13v2 coordinates automatically. |
| **mtDNA Haplogroup** | Validated | Terminal maternal haplogroup from sequence or alignment, with rCRS↔CHM13 mapping. |
| **mtDNA Heteroplasmy** | Validated | Site-level heteroplasmy (depth and allele fraction) from alignments. |
| **Private Y Variants** | Validated | Off-backbone calls — finer branches and novel candidate variants. |
| **Ancestry** | Validated | Admixture across 26 fine populations / 8 continents, PCA projection, a geographic map, and DNA-painting local ancestry. |
| **IBD Detection** | Validated (detection) | Pairwise shared-segment detection and relationship estimates. The match-discovery UI is still in progress. |
| **Structural Variants (SV)** | Built, output unvalidated | Deletions, inversions, and copy-number changes. Reliable output needs ≥10× coverage. |

Navigator also reconciles Y/mtDNA haplogroups across multiple runs and alignments per subject, producing a single consensus assignment.

## The Command Line

The same `navigator` binary is fully scriptable. With a subcommand it opens the *same* workspace database as the GUI, does its work, and exits. This is ideal for bulk-loading a directory of files or querying results.

```bash
# Import everything in a folder into a subject (creating the subject/project if needed)
navigator ingest --subject "Jane Doe" --project "Family Study" --recursive /Volumes/nas/Genomics/jane/

# List all subjects with their data-source counts
navigator subjects

# Show one subject's runs, alignments, profiles, and haplogroup consensus
navigator show --subject "Jane Doe"

# List projects with subject counts
navigator projects
```

Useful flags:
- `--subject` / `-s` — donor identifier (found by exact match, or created on `ingest`).
- `--project` / `-p` — project to assign the subject to (found or created).
- `--sex` — recorded only when a new subject is created (e.g. `male` / `female`).
- `--recursive` / `-r` — recurse into directories instead of importing only their immediate files.
- `--db` — point at an alternate workspace database (defaults to `~/.decodingus/navigator-rs.db`).
- `--json` — emit machine-readable JSON instead of a table (on `subjects`, `show`, `projects`).

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
├── references/          # Downloaded reference genomes (indexed FASTA)
├── liftover/            # Chain files for build-to-build coordinate conversion
├── masks/               # Callable-region BED masks
├── trees/               # Cached Y-DNA / mtDNA haplotrees (JSON)
├── ysnp/                # Y-SNP dictionary assets
├── ancestry/            # Pre-built ancestry panels and PCA loadings
└── navigator-lang       # Your saved UI language choice
```

### Reference Genomes
Navigator manages reference genomes for you. It downloads and caches standard builds (GRCh38, GRCh37, CHM13v2) as needed and converts coordinates between builds automatically — you don't need to hunt for reference files.

### Cloud Integration (Optional)
Navigator includes support for the **AT Protocol** to publish summaries to a Personal Data Store (PDS) in the Decoding-Us Federation.
- **Privacy:** Even with publishing enabled, your raw BAM/CRAM and chip files are **never** uploaded. Only anonymized summaries (haplogroup assignments, coverage QC statistics, ancestry estimates) are shared, with your explicit consent.
- Configure the AppView endpoint with the `DECODINGUS_APPVIEW_URL` environment variable.

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
A: Haplotrees are cached for `NAVIGATOR_TREE_TTL_DAYS` (default 7). Lower that value or set it to `0` to force a fresh fetch, then re-run the analysis.

**Q: My file wasn't recognized on import.**
A: Navigator auto-detects by extension and content fingerprint. Confirm the file is one of the [supported formats](#importing-data). Consumer chip exports from less common vendors may not be detected; the file is still recorded but won't be analyzed.
