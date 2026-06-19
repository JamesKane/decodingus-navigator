# Decoding-Us Navigator - Gemini Context

This document provides a comprehensive overview of the Decoding-Us Navigator project, intended to serve as a contextual reference for the Gemini AI.

## Project Overview

Decoding-Us Navigator is an edge-computing companion application for [decoding-us.com](https://decoding-us.com). Developed entirely in **Rust**, it provides a native graphical user interface built with **[egui](https://github.com/emilk/egui)** for advanced bioinformatics analysis directly on the user's local machine. The analysis engine is built on **[noodles](https://github.com/zaeleus/noodles)** (pure-Rust BAM/CRAM/FASTA/VCF I/O), with a strong emphasis on user privacy through local analysis. It is a single self-contained binary — **no JVM, no GATK/samtools/bcftools** required at build or runtime.

**Key features include:**

*   **Privacy-Preserving Analysis:** All analysis is performed locally; only anonymized summary information is optionally shared via the AT Protocol Personal Data Store (PDS).
*   **User-Friendly Interface:** Simplifies complex bioinformatics for citizen scientists without requiring programming expertise, with a scriptable CLI for power users.
*   **Cross-Platform Compatibility:** A single native binary runs on macOS, Windows, and Linux.
*   **Comprehensive Genomic Analysis:** Coverage/Callable Loci, Read Metrics, Sex Inference, Y-DNA and mtDNA Haplogroup Determination, mtDNA Heteroplasmy, Private Y Variants, Ancestry (admixture/PCA/painting), IBD detection, Structural Variants, and Liftover.
*   **Reference Genome Management:** Automatic download, caching, and conversion between reference builds (GRCh38, GRCh37, CHM13v2).
*   **Analysis Caching:** Avoids redundant analysis with file-hash-based result caching.
*   **Optional Cloud Integration:** AT Protocol authentication (PKCE/DPoP) and PDS publishing for optional summary data upload.

## Technologies Used

*   **Language:** Rust
*   **Build Tool:** Cargo (the workspace is a multi-crate Cargo workspace under `crates/`)
*   **Runtime:** Native binary — no managed runtime
*   **UI Framework:** egui (immediate-mode GUI)
*   **Core Bioinformatics:** noodles (BAM/CRAM/FASTA/VCF)
*   **Persistence:** SQLite via sqlx, with versioned migrations
*   **Async / Parallelism:** tokio (async runtime) + rayon (per-contig data parallelism)
*   **Numerics:** nalgebra (ancestry PCA / mixture models)
*   **CLI:** clap
*   **Serialization:** serde / serde_json (incl. AT-Protocol records)

## Building and Running

The project uses `cargo` for building, testing, and running. The built binary is named `navigator`.

### Build

```bash
cargo build            # or `cargo build --release` for an optimized binary
```

### Run the desktop app

```bash
cargo run -p navigator-ui
```

### Run the headless CLI

The same binary runs headless when given a subcommand (`ingest`, `subjects`, `show`, `projects`):

```bash
cargo run -p navigator-ui -- subjects --json
```

### Run Tests

```bash
cargo test --workspace
```

### Lint (per-commit gate)

```bash
cargo clippy --all-targets -- -D warnings
```

## Development Conventions

*   **Language:** Idiomatic Rust — clear ownership, `Result`-based error handling, small focused modules.
*   **Build System:** `cargo`; the workspace and dependencies are defined in the root `Cargo.toml` and per-crate `Cargo.toml` files.
*   **Crate Layering:** Respect the dependency rule `ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`. Shared crates (`du-domain`, `du-atproto`, `du-bio`) live in the sibling repo `../decodingus-shared/`.
*   **UI Development:** egui in `navigator-ui` is a thin shell — view-state and dispatch only; business logic lives behind the `navigator-app` API.
*   **Error Handling:** Robust handling for file I/O, network operations, and large-file parsing; surface failures rather than panicking.
*   **Concurrency:** Use tokio for async I/O and rayon for CPU-bound per-contig fan-out to keep the UI responsive during long analyses.
*   **Testing & Linting:** `cargo test --workspace`; keep `cargo clippy --all-targets -- -D warnings` clean. Some integration tests are `#[ignore]` and gated on env vars (real BAMs / network).

See **[`crates/README.md`](crates/README.md)** for the crate topology, **[`documents/design/RustRewrite_Plan.md`](documents/design/RustRewrite_Plan.md)** for the design, and **[`documents/design/HANDOFF.md`](documents/design/HANDOFF.md)** for resume notes.
