# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Run Commands

```bash
# Build the whole workspace
cargo build

# Run the desktop app (egui GUI)
cargo run -p navigator-ui

# Run the headless CLI (same binary, with a subcommand)
cargo run -p navigator-ui -- subjects --json

# Run all tests (some integration tests are #[ignore] — live BAMs / network)
cargo test --workspace

# Lint gate (must be clean per commit)
cargo clippy --all-targets -- -D warnings

# Run a single test
cargo test -p navigator-analysis some_test_name
```

The built binary is named `navigator` (`target/debug/navigator` or `target/release/navigator`). Run with no subcommand to launch the GUI; run with `ingest` / `subjects` / `show` / `projects` for headless mode.

## Architecture Overview

Decoding-Us Navigator is a Rust desktop application for local bioinformatics analysis of BAM/CRAM files (and consumer DNA data). It is a single self-contained binary — **no JVM, and no external bioinformatics tools** (no GATK, samtools, or bcftools). The analysis engine is built on [noodles](https://github.com/zaeleus/noodles); the UI is [egui](https://github.com/emilk/egui).

> This is a ground-up Rust rewrite (branch `rust-rewrite`) of an earlier Scala/ScalaFX app. The
> authoritative crate map is **[`crates/README.md`](crates/README.md)**; the design is in
> **[`documents/design/RustRewrite_Plan.md`](documents/design/RustRewrite_Plan.md)** and resume notes
> in **[`documents/design/HANDOFF.md`](documents/design/HANDOFF.md)**.

### Key Architectural Patterns

**Command/Query API**: The UI is a thin egui shell holding view-state and dispatching to a single application API in `navigator-app` (`App`). Business logic lives behind that API, not in the UI.

**Workspace Persistence**: Workspace state is stored in a SQLite database at `~/.decodingus/navigator-rs.db` via `navigator-store` (sqlx, with versioned migrations). The workspace holds subjects (biosamples), projects, sequence runs, alignments, and profiles. The same database backs both the GUI and the CLI.

**Crate dependency rule**: `ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`. Respect this layering — lower crates must not depend on higher ones.

### Crate Structure (under `crates/`)

- `navigator-domain` — Pure desktop aggregate types; re-exports shared `du-domain`.
- `navigator-analysis` — The htsjdk/GATK replacement: noodles BAM/CRAM/FASTA I/O, coverage/callable, haploid caller, Y & mtDNA haplogroups, IBD, sex, read_metrics, sv, ancestry (admixture/PCA/painting).
- `navigator-store` — SQLite (sqlx) persistence + versioned migrations.
- `navigator-refgenome` — Reference/chain retrieval, on-disk cache, and liftover gateway.
- `navigator-sync` — AT-Proto OAuth (PKCE/DPoP) + PDS record publishing.
- `navigator-app` — The single command/query API the UI dispatches to.
- `navigator-ui` — egui desktop shell + the `navigator` binary (GUI + clap CLI).
- `navigator-panelbuild` — Offline tool (not shipped): builds ancestry panels/PCA assets.

Shared crates (`du-domain`, `du-atproto`, `du-bio`) live in the sibling repo `../decodingus-shared/crates/` and are wired in as path dependencies during co-development.

### Analysis Flow

1. A BAM/CRAM file (or VCF/GVCF, chip raw data, STR/Y-SNP export, mtDNA FASTA) is imported via the UI or the `ingest` CLI; `app.add_data` auto-detects the type.
2. A header probe infers reference build / aligner / platform / test type.
3. `navigator-refgenome` resolves/downloads the appropriate reference genome and chains.
4. Analysis runs: coverage/callable, read metrics, sex, SV, Y/mtDNA haplogroups, ancestry (parallelized per contig via rayon).
5. Results are persisted to SQLite and cached as on-disk artifacts under `~/.decodingus/`.
6. Summary records can optionally be published to a PDS via `navigator-sync`.

### Storage Layout (`~/.decodingus/`)

`navigator-rs.db` (SQLite workspace), `references/`, `liftover/`, `masks/`, `trees/` (cached haplotrees), `ysnp/` (Y-SNP dictionary), `ancestry/` (pre-built panels/PCA).

### Key Dependencies

- **noodles** — Pure-Rust BAM/CRAM/FASTA/VCF I/O (replaces HTSJDK/GATK).
- **egui** — Immediate-mode GUI (replaces ScalaFX).
- **sqlx** — Async, compile-time-checked SQLite access.
- **tokio** — Async runtime.
- **rayon** — Per-contig data parallelism.
- **nalgebra** — Linear algebra (ancestry PCA / mixtures).
- **clap** — Headless CLI.
- **serde / serde_json** — Serialization and federated records.

### Useful Environment Variables

`NAVIGATOR_ANALYSIS_THREADS`, `NAVIGATOR_BGZF_THREADS`, `NAVIGATOR_Y_TREE_PROVIDER` (`decodingus`/`ftdna`), `NAVIGATOR_TREE_TTL_DAYS`, `NAVIGATOR_REFGENOME_DIR`, `NAVIGATOR_TREE_DIR`, `NAVIGATOR_ANCESTRY_PANEL` / `NAVIGATOR_ANCESTRY_PCA`, `DECODINGUS_APPVIEW_URL`.
