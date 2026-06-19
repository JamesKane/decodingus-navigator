# DEVELOPER.md: Onboarding New Contributors

Welcome to the Decoding-Us Navigator project! This document gives new contributors an overview of the technology stack and where contributions can make the most impact.

## Project Overview

Decoding-Us Navigator is an edge-computing companion to [decoding-us.com](https://decoding-us.com). It enables local, privacy-preserving analysis of BAM/CRAM files (and consumer DNA data) for citizen scientists. It is a **single self-contained Rust application**: a pure-Rust analysis engine plus a native desktop UI, with **no JVM and no external bioinformatics tools** (no GATK, samtools, or bcftools) at build or runtime.

## Technology Stack

### Language
*   **Rust:** The entire application is Rust, chosen for performance, memory safety, a single dependency-free binary, and a strong type system.

### Build Tool
*   **Cargo:** Standard Rust tooling for building, testing, and running. The repository is a Cargo workspace of several crates.

### Core Libraries & Frameworks
*   **[noodles](https://github.com/zaeleus/noodles):** Pure-Rust I/O and processing for BAM/CRAM/FASTA/VCF — the GATK/HTSJDK replacement.
*   **[egui](https://github.com/emilk/egui):** Immediate-mode GUI for the cross-platform desktop interface (replacing ScalaFX).
*   **[sqlx](https://github.com/launchbadge/sqlx):** Async, compile-time-checked SQLite access for the local workspace database, with versioned migrations.
*   **[tokio](https://tokio.rs):** Async runtime for I/O and background analysis tasks.
*   **[rayon](https://github.com/rayon-rs/rayon):** Data parallelism for per-contig analysis fan-out.
*   **[nalgebra](https://nalgebra.org):** Linear algebra for ancestry PCA / mixture models.
*   **[clap](https://github.com/clap-rs/clap):** The headless CLI (`navigator ingest/subjects/show/projects`).
*   **[serde](https://serde.rs) + serde_json:** Serialization, including federated PDS records.

### Data Formats
*   **BAM / CRAM:** Aligned sequencing reads (input).
*   **VCF / GVCF:** Variant calls (input and intermediate output).
*   **FASTA:** Reference genomes and mtDNA sequences.
*   **SQLite:** The local workspace database.
*   **JSON:** Cached haplotrees, config, and AT-Protocol summary records.

## Project Structure

The repository is a Cargo workspace. Crates live under `crates/`:

```
crates/
├── navigator-domain      # Pure desktop aggregate types (re-exports shared du-domain)
├── navigator-analysis    # The htsjdk/GATK replacement: noodles I/O, coverage, caller,
│                         #   haplogroups, mtDNA, IBD, sex, read_metrics, sv, ancestry
├── navigator-store       # SQLite (sqlx) persistence + versioned migrations
├── navigator-refgenome   # Reference/chain retrieval + on-disk cache + liftover gateway
├── navigator-sync        # AT-Proto OAuth (PKCE/DPoP) + PDS record publishing
├── navigator-app         # The single command/query API the UI dispatches to
├── navigator-ui          # egui desktop shell + the `navigator` binary (GUI + clap CLI)
└── navigator-panelbuild  # Offline tool (not shipped): builds ancestry panels/PCA assets
```

**Dependency rule:** `ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`.

Shared crates (`du-domain`, `du-atproto`, `du-bio`) live in the sibling repo `../decodingus-shared/crates/` and are wired in as path dependencies during co-development.

See **[`crates/README.md`](crates/README.md)** for the authoritative crate topology, the design in **[`documents/design/RustRewrite_Plan.md`](documents/design/RustRewrite_Plan.md)**, and resume notes in **[`documents/design/HANDOFF.md`](documents/design/HANDOFF.md)**.

## Areas for Contribution

### 1. Performance Improvements
*   **Analysis throughput:** The walkers parallelize per-contig (rayon). Look for further fan-out opportunities and ways to reduce passes over large BAM/CRAM files.
*   **BGZF / I/O:** Decompression and point-query patterns for CRAM and tabix-indexed sites.
*   **Memory bounds:** Keep peak memory predictable for WGS-scale inputs (load semaphores, N-masking).

### 2. Maintainability & Code Quality
*   **Idiomatic Rust:** Favor clear ownership, `Result`-based error handling, and small focused modules.
*   **Clippy clean:** `cargo clippy --all-targets -- -D warnings` is a per-commit gate.
*   **Documentation:** Improve `///` doc comments for non-obvious analysis logic and public APIs.

### 3. Test Coverage
*   **Unit tests:** Expand coverage in `navigator-analysis` and `navigator-store`.
*   **Integration tests:** Some are `#[ignore]` and gated on env vars (real BAMs / network) — see HANDOFF.
*   **Parity harness:** The GATK-vs-Rust golden-gate tests exist but are `#[ignore]`; automating them is open work.

### 4. New Features & Enhancements
*   **IBD matching system:** Detection math is done; the consent / match-discovery / chromosome-browser UI is pending.
*   **SV validation:** SV output is unvalidated and needs a ≥10× sample to validate.
*   **Federation sync:** Granular per-record publishing paths exist; broaden coverage.
*   **UI/UX:** New visualizations, better progress feedback, more intuitive workflows in the egui shell.

## Getting Started

1.  **Install Rust:** Use [rustup](https://www.rust-lang.org/tools/install).
2.  **Clone the repository:** `git clone [repository-url]` (the Rust workspace lives on `main`). For co-development of shared crates, also clone `decodingus-shared` as a sibling directory.
3.  **Build the workspace:**
    ```bash
    cargo build
    ```
4.  **Run the desktop app:**
    ```bash
    cargo run -p navigator-ui
    ```
5.  **Run the tests and linter:**
    ```bash
    cargo test --workspace
    cargo clippy --all-targets -- -D warnings
    ```

We look forward to your contributions! If you have any questions, please don't hesitate to reach out.
