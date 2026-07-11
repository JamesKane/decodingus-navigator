# Decoding-Us Navigator

Decoding-Us Navigator is an edge-computing companion application to [decoding-us.com](https://decoding-us.com). It analyzes BAM/CRAM files and consumer DNA data directly on your local machine, empowering citizen scientists with advanced bioinformatics while preserving privacy.

Navigator is a **single self-contained Rust application**: a pure-Rust analysis stack ([noodles](https://github.com/zaeleus/noodles)) and an [egui](https://github.com/emilk/egui) desktop UI — **no JVM, no GATK, no samtools/bcftools** to install. The same binary doubles as a scriptable command-line tool.

> **Repository note:** the Rust rewrite is the trunk on `main`. The original Scala/ScalaFX
> implementation is preserved in git history only. The crate topology and developer setup are in
> **[`crates/README.md`](crates/README.md)**; the design of the rewrite is documented in
> **[`documents/design/RustRewrite_Plan.md`](documents/design/RustRewrite_Plan.md)** with resume
> notes in **[`documents/design/HANDOFF.md`](documents/design/HANDOFF.md)**. End-user documentation
> lives in the **[User Guide](USER_GUIDE.md)**.

## Privacy-Preserving Analysis

All analysis runs locally. Your raw genomic files never leave your machine. Only anonymized summaries are optionally shared, with your explicit consent:
- Haplogroup assignments (Y-DNA and mtDNA)
- General coverage / QC statistics
- Ancestry estimates
- IBD (autosomal match) attestations

Data sharing uses the AT Protocol Personal Data Store (PDS) for user-controlled data ownership. Even with publishing enabled, your raw BAM/CRAM and chip files are **never** uploaded.

## Goal

Decoding-Us Navigator wraps complex bioinformatics in an intuitive interface. It is designed for hobbyists and citizen scientists, making advanced genetic analysis accessible without requiring programming expertise — while remaining fully scriptable for power users.

## Cross-Platform

A single native binary runs on macOS, Windows, and Linux with no runtime dependencies to install.

## The Workbench
![Workbench Screenshot](documents/images/MainWorkbench.png)

## Features

### Workspace Management
- Create and manage multiple projects and subjects (biosamples)
- Subject-centric detail view: Overview, Y-DNA, mtDNA, Autosomal, Ancestry, IBD Matches, and Sources
- Persistent workspace stored locally in SQLite
- Search and filter across projects and subjects

### Data Import (auto-detected)
- **BAM / CRAM** aligned reads
- **VCF / GVCF** variant calls (GVCF carries callable-region context for a fast haplogroup path)
- **mtDNA FASTA** sequences
- **Consumer chip raw data** — 23andMe, AncestryDNA, MyHeritage, Living DNA, FTDNA (Y and mtDNA haplogroups placed on import)
- **CompleteGenomics masterVar** (`var-*-ASM.tsv[.bz2]`) genome-wide variant tables
- **Y-STR profiles** — FTDNA/YSEQ-style CSV/TSV exports
- **Y-SNP panels** — BISDNA chromo2 genotyped exports
- **FTDNA Big Y** — named-variant CSV exports

Imports automatically compute a checksum and detect platform (Illumina, PacBio, Oxford Nanopore, MGI, Ion Torrent, Complete Genomics) and test type (WGS, WES, HiFi, CLR, Nanopore, Targeted Panel).

**Project Import (batch + sidecar fast path):** import a whole `project/sample/…` directory tree at once. When per-sample pipeline "sidecar" files (`*.chrY.g.vcf.gz`, `*.chrM.g.vcf.gz`, `*.sex`, `stats.txt`, `coverage.txt`) sit next to the CRAM, Navigator places haplogroups, sex, read metrics, and a lite coverage roll-up from them — no multi-GB CRAM walk (seconds instead of minutes per sample). See the [User Guide](USER_GUIDE.md#project-import-batch-with-the-sidecar-fast-path) for the required layout and rules.

### Analysis Capabilities
- **Coverage / Callable Loci** — mean depth, coverage distribution, per-contig depth histograms, callable bases per contig (1×–100×)
- **Read Metrics** — read length, insert size, platform/instrument detection, library orientation, sequencing-lab inference
- **Sex Inference** — genetic sex with a confidence score
- **Y-DNA & mtDNA Haplogroups** — terminal assignment with ranked candidates, across the DecodingUs and FTDNA tree providers; multi-source reconciliation into a single genome-level consensus per subject
- **mtDNA Variants & Heteroplasmy** — rCRS-relative mutation list plus site-level depth and allele fraction
- **Private Y Variants** — off-backbone calls (finer branches + novel candidates), reconciled across sources
- **Diploid Variant Calling** — de-novo diploid SNV + indel caller, exportable as a whole-genome VCF (per alignment or per-subject consensus)
- **Ancestry** — admixture (26 fine populations / 8 continents), PCA projection, geographic map, DNA-painting local ancestry
- **IBD Detection** — pairwise shared-segment detection with a per-chromosome segment browser and relationship estimates, using a real recombination map. Federated match suggestions surface candidate relatives from the Federation, and an encrypted, consent-gated channel exchanges IBD segments and signed attestations between edges (end-to-end round-trip validation is gated on a live AppView broker).
- **Structural Variants** — deletions, duplications, inversions, breakends (output unvalidated; needs ≥10× coverage)
- **Liftover** — automatic coordinate conversion between GRCh38, GRCh37, and CHM13v2

### Reference Genome Management
- Automatic reference download and caching (GRCh38, GRCh37, CHM13v2)
- Configurable cache directory and on-demand retrieval with size estimates

### Analysis Caching
- Result caching keyed by input file signature to avoid redundant work
- Subject-organized artifact storage for intermediate analysis files

### Data Storage

All application data is stored under `~/.decodingus/`:

```
~/.decodingus/
├── navigator-rs.db      # Workspace database (SQLite): subjects, projects, runs, alignments, profiles
├── config/              # settings.json (saved preferences)
├── references/          # Downloaded reference genomes (indexed FASTA)
├── liftover/            # Chain files for build-to-build coordinate conversion
├── masks/               # Callable-region BED masks
├── trees/               # Cached Y-DNA / mtDNA haplotrees (JSON)
├── ysnp/                # Y-SNP dictionary assets
└── ancestry/            # Pre-built ancestry panels and PCA loadings
```

The database is created automatically on first launch. No configuration is required. Query it directly with any SQLite tool (close the app first):

```bash
sqlite3 ~/.decodingus/navigator-rs.db ".tables"
```

### Optional Cloud Integration
- AT Protocol OAuth (PKCE/DPoP) and PDS record publishing
- Optional upload of summary data with user consent; durable retry with backoff if the network or PDS is briefly unavailable
- AppView endpoint configured via `DECODINGUS_APPVIEW_URL`

## Requirements

- No runtime dependencies to install (no Java, no external bioinformatics tools)
- 4 GB RAM minimum (8 GB recommended for large BAM files)
- A [Rust toolchain](https://www.rust-lang.org/tools/install) to build from source

## Building & Running

```bash
# Build the whole workspace (release = optimized)
cargo build --release

# Run the desktop app
cargo run -p navigator-ui
# ...or run the built binary directly
./target/release/navigator

# Run tests (some integration tests are #[ignore] — live BAMs / network)
cargo test --workspace

# Lint gate (must be clean per commit)
cargo clippy --all-targets -- -D warnings
```

### Command-line use

The same `navigator` binary runs headless when given a subcommand, against the same workspace database:

```bash
navigator ingest --subject "Jane Doe" --project "Family Study" --recursive /path/to/files
navigator subjects --json
navigator show --subject "Jane Doe"
navigator projects
navigator call --subject "Jane Doe" --out jane.vcf          # de-novo diploid SNV/indel VCF
navigator liftvcf --in calls.GRCh38.vcf.gz --from GRCh38 --to chm13v2.0 --out calls.chm13.vcf.gz
```

Additional subcommands cover PDS sign-in (`login`), diagnostics, and maintenance. See the **[User Guide](USER_GUIDE.md)** for full usage, and **[`crates/README.md`](crates/README.md)** for the crate topology and developer setup.
