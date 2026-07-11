# Decoding-Us Navigator — Rust rewrite

A ground-up Rust port of the Scala/ScalaFX Navigator. Replaces GATK/HTSJDK with a pure-Rust
analysis stack ([noodles](https://github.com/zaeleus/noodles)) and ScalaFX with an
[egui](https://github.com/emilk/egui) desktop UI — **no JVM, no GATK/samtools/bcftools at
runtime**. The Rust rewrite is the trunk on `main`; the legacy Scala implementation was removed at
cutover and survives in git history only.

Design: [`documents/design/RustRewrite_Plan.md`](../documents/design/RustRewrite_Plan.md).
Resume notes: [`documents/design/HANDOFF.md`](../documents/design/HANDOFF.md).

## Workspace topology

Dependency rule: `ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`.

| Crate | Role |
|-------|------|
| `navigator-domain` | Pure desktop-only aggregate types (re-exports shared `du-domain`). |
| `navigator-analysis` | The htsjdk/GATK replacement: noodles BAM/CRAM/FASTA/VCF I/O, coverage/callable, haploid + **diploid SNV/indel** callers, Y & mtDNA haplogroups, IBD, **sex**, **read_metrics**, **sv**, **ancestry** (admixture/PCA/painting), CompleteGenomics masterVar. |
| `navigator-store` | SQLite (`sqlx`) persistence, versioned migrations. |
| `navigator-refgenome` | Reference/chain retrieval + on-disk cache + liftover gateway. |
| `navigator-sync` | AT-Proto OAuth (PKCE/DPoP) + PDS record publishing. |
| `navigator-app` | The single command/query API the UI dispatches to. |
| `navigator-ui` | egui desktop shell (thin: view-state + dispatch only). |
| `navigator-panelbuild` | **Offline tool** (not shipped): builds the ancestry panels/PCA/fine assets from 1000G+SGDP genotype data. |

Shared crates live in the sibling repo `../decodingus-shared/crates/{du-domain,du-atproto,du-bio}`
(path deps during co-development).

## Build / run / test

```bash
cargo build                                   # whole workspace
cargo run -p navigator-ui                      # the desktop app
cargo test --workspace                         # unit + integration (ignored = live/network)
cargo clippy --all-targets -- -D warnings      # per-commit gate (must be clean)
```

Some integration tests are `#[ignore]` and gated on env vars (real BAMs / network) — see HANDOFF.

## Feature status (high level)

**Built + validated:** workspace/import (incl. chip, Y-STR, BISDNA Y-SNP, CompleteGenomics
masterVar, batch project import + sidecar fast path), BAM/CRAM reader, coverage/callable, haploid
caller, **diploid SNV/indel caller** + whole-genome VCF export, Y & mtDNA haplogroups incl. CHM13
liftover (chain rev-comp + rotation-aware rCRS↔chrM map), multi-source reconciliation into a
genome-level consensus, IBD detection + identity + per-chromosome segment browser, panel genotyping,
refgenome gateway, OAuth + durable publishing, **ancestry** (admixture composition, 26 fine pops
across 8 continents, PCA scatter, geographic map, DNA-painting local ancestry), and **sex /
read-metrics / SV** wired into both single-alignment and the bulk `analyze_project` flow. Federated
IBD is live: device-key-signed records, an encrypted edge exchange channel, signed attestations, and
AppView-mined network match suggestions surfaced in the UI.

**Pending:** live-network validation of federated IBD exchange (the encrypted channel, consent
round-trip, signed attestations, and UI inbox are implemented and unit-tested, but the full edge-to-
edge round-trip needs a running AppView broker + a partner edge online, and the AppView's symmetric
counterpart-discovery is still being speced); parity-harness automation (GATK-vs-Rust golden gate —
tests exist but `#[ignore]`). SV *output* is unvalidated (needs a ≥10× sample; the coverage gate
works). `navigator-migrate` (H2→SQLite) is **dropped** — no pre-beta data continuity needed.

## Offline-built assets (not committed; regenerable)

Ancestry panels/PCA live in `~/.decodingus/ancestry/` (built by `navigator-panelbuild` from the
1000G+SGDP genotype matrix archived at `~/Genomics/archive/1kgp_chm13_pca_build/`). The app loads
them via `$NAVIGATOR_ANCESTRY_PANEL` / `$NAVIGATOR_ANCESTRY_PCA` or the cache dir. See HANDOFF for
the build recipes and the EC2 genotype-extraction workflow.
