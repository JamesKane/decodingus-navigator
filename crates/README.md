# Decoding-Us Navigator — Rust rewrite

A ground-up Rust port of the Scala/ScalaFX Navigator. Replaces GATK/HTSJDK with a pure-Rust
analysis stack ([noodles](https://github.com/zaeleus/noodles)) and ScalaFX with an
[egui](https://github.com/emilk/egui) desktop UI — **no JVM, no GATK/samtools/bcftools at
runtime**. Lives on branch `rust-rewrite`; coexists with the legacy Scala build until cutover.

Design: [`documents/design/RustRewrite_Plan.md`](../documents/design/RustRewrite_Plan.md).
Resume notes: [`documents/design/HANDOFF.md`](../documents/design/HANDOFF.md).

## Workspace topology

Dependency rule: `ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`.

| Crate | Role |
|-------|------|
| `navigator-domain` | Pure desktop-only aggregate types (re-exports shared `du-domain`). |
| `navigator-analysis` | The htsjdk/GATK replacement: noodles BAM/CRAM/FASTA I/O, coverage, caller, haplogroups, mtDNA, IBD, **sex**, **read_metrics**, **sv**, **ancestry** (admixture/PCA/painting). |
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

**Built + validated:** workspace/import, BAM/CRAM reader, coverage/callable, haploid caller,
Y & mtDNA haplogroups incl. CHM13 liftover (chain rev-comp + rotation-aware rCRS↔chrM map),
multi-source reconciliation (6 phases), IBD detection + identity, panel genotyping, refgenome
gateway, OAuth + publishing, batch import + project report, **ancestry** (admixture composition,
26 fine pops across 8 continents, PCA scatter, geographic map, DNA-painting local ancestry), and
**sex / read-metrics / SV** wired into both single-alignment and the bulk `analyze_project` flow.

**Pending:** parity-harness automation (GATK-vs-Rust golden gate — tests exist but `#[ignore]`);
IBD **matching** system (consent/match-discovery/chromosome-browser UI; detection math is done);
granular per-record sync (publish paths exist). `navigator-migrate` (H2→SQLite) is **dropped** (no
pre-beta data continuity needed). SV *output* is unvalidated (needs a ≥10× sample; gate works).

## Offline-built assets (not committed; regenerable)

Ancestry panels/PCA live in `~/.decodingus/ancestry/` (built by `navigator-panelbuild` from the
1000G+SGDP genotype matrix archived at `~/Genomics/archive/1kgp_chm13_pca_build/`). The app loads
them via `$NAVIGATOR_ANCESTRY_PANEL` / `$NAVIGATOR_ANCESTRY_PCA` or the cache dir. See HANDOFF for
the build recipes and the EC2 genotype-extraction workflow.
