# DUNavigator → Rust Rewrite Plan

**Status:** Planning
**Date:** 2026-06-01
**Decisions locked:** egui/eframe GUI · purpose-built pure-Rust haploid caller (no GATK/JVM) · shared crates extracted, Navigator in its own Cargo workspace

---

## 1. Why this is more tractable than it looks

The sister repo `/Users/jkane/Development/decodingus` (web rewrite) has already built
the foundation we need. It is a clean Cargo workspace whose lower layers are
**desktop-reusable as-is**:

| decodingus crate | What it gives Navigator |
|:---|:---|
| `du-domain` | Pure domain types, strongly-typed IDs, serde enums — **zero IO**. Directly fixes our type-triplication debt. |
| `du-atproto` | DID/handle parsing, did:key verification, PDS resolution, **OAuth: PKCE, ES256, DPoP, PAR** — the exact work from the AT Proto OAuth plan. |
| `du-bio` | **Coordinate math + text-format parsing only** (pure Rust, IO-light): `callable` (BED interval merge), `liftover` (UCSC chain), `vcf` (text reader for catalog ingest), `ybrowse`. **BAM/CRAM I/O and variant calling are explicitly out of scope** — those are Navigator-side (`navigator-analysis`). Depends only on `du-domain` + `thiserror`; no noodles. |
| `du-db` (pattern) | sqlx query-module-per-aggregate + JSONB→domain mapping pattern to copy for SQLite. |

> **Scope note (resolved at extraction, 2026-06-01):** an earlier draft put the noodles
> BAM/CRAM I/O and the walker ports inside shared `du-bio`. That was narrowed: shared `du-bio`
> stays coordinate-math + text-parsing (so the decodingus server never pulls in noodles), and
> **all noodles raw-read I/O + the walkers live in Navigator's `navigator-analysis` crate.** The
> haploid caller was already Navigator-side. See §2 and §4a for the corrected ownership.

The web UI (`du-web`: Axum + Askama + **HTMX**, server-rendered) is **not** reusable as a
desktop UI — but its domain/business logic and the crates above are. Hence: **egui for
the UI, reuse the lower crates.**

Master architecture reference for the sister repo: `~/.claude/plans/robust-knitting-lampson.md`.

---

## 2. Target topology — shared crates + Navigator workspace

The genuinely shared crates are extracted from the decodingus repo into
**`decodingus-shared`**, a sibling repo at the Development root. All three projects sit
side-by-side so both apps depend on the shared crates and fixes flow both ways:

```
/Users/jkane/Development/
├─ decodingus/                  # web app (Axum + Askama + HTMX)
├─ decodingus-shared/           # extracted shared crates, versioned independently
│    du-domain/                 # types, IDs, enums  (no IO)
│    du-atproto/                # DID + OAuth/DPoP    (identity)
│    du-bio/                    # coordinate math + text parsing: liftover + callable + vcf/ybrowse (NO BAM/CRAM, NO caller)
└─ DUNavigator/                 # desktop app — moved up from the legacy scala/ parent
     crates/
       navigator-domain/        # desktop-only types not in du-domain
                                #   (SequenceRun, Alignment, AnalysisArtifact,
                                #    YProfile, IBD, Workspace/Project aggregate)
       navigator-store/         # local persistence: SQLite via sqlx
       navigator-analysis/      # noodles BAM/CRAM I/O + the walkers + the haploid caller
                                #   (uses du-bio for liftover/callable/coordinate math)
       navigator-sync/          # PDS push/pull + AsyncSync (completed, not stubbed)
       navigator-app/           # application/command layer (no UI, no ScalaFX analog)
       navigator-ui/            # egui/eframe front end (thin: render + dispatch)
     Cargo.toml                 # workspace; depends on ../decodingus-shared/* (path or git)
```

**Dependency rule (the antidote to today's god object):** dependencies point *down only*
— `ui → app → {analysis, store, sync} → domain`. The UI never touches the DB or a
processor directly; it dispatches commands and renders state. No crate below `app`
imports anything UI-related.

---

## 3. Persistence — SQLite via sqlx (replaces H2/Slick)

- **Engine:** SQLite (embedded, file-based, perfect for a desktop app), accessed through
  `sqlx` (already in the workspace; supports the `sqlite` feature). Cross-platform, no
  server, single file under `~/.decodingus/`.
- **Kill the Slick-22-tuple JSONB workaround.** sqlx has no tuple-arity limit, so model
  complex children as **proper rows** (e.g. `file_info`, `contig_metrics`, `ibd_segment`
  tables) instead of `asJson.noSpaces` blobs. Reserve JSON columns for genuinely
  schema-less payloads (AT Proto record snapshots), using `sqlx::types::Json<T>` —
  exactly the `du-db` pattern.
- **Migrations:** `sqlx migrate` (versioned, checked-in SQL), replacing the 12 forward-only
  Flyway-style files. Author down-migrations.
- **One source of truth.** Persisted state is authoritative; the UI reads from a
  projection. No more `_workspace.value = …` imperative mutation racing async H2 writes
  (the current race-condition hotspot in `WorkbenchViewModel`).

---

## 4. GATK / HTSJDK replacement strategy

`noodles` (pure Rust) replaces htsjdk for BAM/CRAM/VCF/FASTA/BGZF/index I/O — and crucially
gives us a **clean Windows build** that `rust-htslib` (C bindings) would not. **This noodles
I/O lives in Navigator's `navigator-analysis` crate, not in shared `du-bio`** (which stays
IO-light coordinate math + text parsing, per §1). Navigator depends on `du-bio` for the
liftover/callable/coordinate primitives and adds the raw-read layer on top.

### 4a. The walkers → mechanical port onto noodles pileups

These already reimplement GATK logic in Scala-over-htsjdk; they are pileup loops with
statistics and port almost directly. **All targets are in `navigator-analysis`**; the
`du-bio` references are the shared coordinate/text primitives it builds on:

| Scala component | Replaces (GATK) | Port target |
|:---|:---|:---|
| `CoverageCallableWalker` | `CollectWgsMetrics` + `CallableLoci` | `navigator-analysis::coverage` over `noodles` pileup, using `du-bio::callable` for the BED merge/summary |
| `UnifiedMetricsWalker` | `CollectAlignmentSummaryMetrics` + `CollectInsertSizeMetrics` | `navigator-analysis::read_metrics` |
| `SvEvidenceWalker` | (custom BreakDancer/Pindel-style) | `navigator-analysis::sv` |
| `SexInference` | (index-based ratio) | `navigator-analysis::sex` (noodles BAI metadata) |
| `LiftoverProcessor` | `LiftoverVcf` | `du-bio::liftover` (UCSC chain parse — pure coordinate math, shared) |
| `ReferenceQuerier` | htsjdk FASTA | `navigator-analysis` FASTA via noodles, contig caching |

### 4b. The genuine gap: variant calling — a purpose-built **haploid caller**

There is no pure-Rust GATK, and your flows need three things from it. **Confirmed
requirement:** de-novo Y/mtDNA calls are needed for **private-variant matching and new
branch creation** — so the caller is not just a force-call genotyper. It has two modes:

1. **Force-call genotyping at known sites** (haplogroup tree sites, ancestry-informative
   SNPs). Pileup at each target position → call ref/alt by depth + base-quality majority.
   Straightforward.
2. **De-novo discovery across the contig** (Y and mtDNA) → **private variants** for branch
   creation. Because Y and mtDNA are **haploid (ploidy = 1)**, this is *pileup-consensus
   calling*, not the diploid local-reassembly that makes `HaplotypeCaller`/`Mutect2`
   expensive to reproduce. Algorithm: walk callable positions, compute the consensus
   non-reference allele where depth ≥ min and base-quality/MAPQ filters pass and allele
   fraction ≥ threshold; subtract known tree positions → private set. mtDNA runs at high
   depth over 16.5 kb (cheap); Y is large but limited to reliably callable regions.

**The risk to validate, stated honestly:** SNPs are tractable; **indels and homopolymer
runs are where a naive pileup caller diverges from GATK** — local reassembly is exactly
what improves indel accuracy, and mtDNA homopolymer/length-heteroplasmy sites
(e.g. ~309, ~16193, the 3107 "N") are notorious. Mitigations, in order of preference:
   - Add a **light local realignment** around candidate indels before calling.
   - Restrict branch-defining/private calls to **SNPs** initially; treat indels as
     advisory until validated.
   - Keep a **pinned external caller fallback** (see §4d) usable for indel-sensitive
     mtDNA work during validation only.

`CollectWgsMetrics`, `BuildBamIndex`, `IndexFeatureFile`: not separate algorithms —
metrics fold into the coverage walker; indexing is a noodles call.

### 4c. Validation — parity is a first-class deliverable

Before the Rust caller is trusted, build a **golden-truth harness**: run the existing
GATK pipeline and the new Rust caller on a panel of samples (varied depth, platforms,
builds; include known-hard mtDNA homopolymer cases) and assert agreement on called
genotypes, private-variant sets, and haplogroup assignment. This harness is also the
regression guard for the rewrite. No flow flips to the Rust caller until it passes.

### 4d. Transition bridge (optional, time-boxed)

If the caller's validation lags the rest of the rewrite, the JVM GATK can run as a
**subprocess sidecar** behind the same `navigator-analysis` interface so UI/domain work
proceeds unblocked — explicitly temporary, removed once §4c passes. (This is *not* the
shipped architecture; the end state is JVM-free.)

---

## 5. GUI — egui / eframe

- **Single static binary** per OS, GPU-rendered, no system webview or JRE — the simplest
  robust cross-platform story for Win/Lin/macOS, which is the stated priority.
- **Architecture:** `navigator-ui` is thin. It renders immutable view-state and dispatches
  **commands** to `navigator-app`; long-running analysis runs on a worker thread (or
  `tokio` runtime) and streams progress back via a channel, which the egui repaint loop
  reads. No business logic, no DB calls, no domain decisions in the UI — directly fixing
  the "dialogs make domain decisions" debt (`FingerprintMatchDialog`, etc., become app-layer
  policy with a UI prompt only when truly needed).
- **What maps over:** the 37 ScalaFX dialogs/panels become egui panels/windows; data
  tables and coverage/haplogroup charts are an egui strength. `SubjectDetailView`
  (3,104 LOC) decomposes into per-tab widgets backed by app-layer queries.
- Immediate-mode is awkward for very large forms — accept that; the workbench is
  table/chart/report-heavy, which suits egui.

---

## 6. Tech-debt remediation — explicit mapping

| Today (Scala) | Rewrite fix |
|:---|:---|
| `WorkbenchViewModel` god object (4,021 LOC) | Split across `navigator-app` (commands), `navigator-store`, `navigator-analysis`; `navigator-ui` holds only view-state + dispatch. |
| `HaplogroupResult` / `ScoredHaplogroup` / `HaplogroupAssignments` triplication | One domain type in `du-domain`/`navigator-domain`; serialization variants via serde, not parallel types. |
| Slick 22-tuple → JSONB blobs | Proper SQLite tables; `Json<T>` only for AT Proto snapshots (§3). |
| `EntityConversions` (549 LOC) fragile mapping | sqlx `FromRow` + small `into_domain` per aggregate (the `du-db` pattern). |
| Mixed error handling (`Either`/`Try`/`Option`/exceptions; swallowed failures) | `thiserror` enums per layer (`DomainError`/`StoreError`/`AppError`), propagate with `?`; no proceed-after-failed-persist. |
| `AsyncSyncService` stubbed | `navigator-sync` completed: explicit retry/backoff, conflict policy, offline indicator; reuse `du-atproto` OAuth session. |
| Config fragmentation (3 sources, deprecated overrides) | One validated config (env > file > default). |
| Cache with no invalidation/versioning | Cache key includes **algorithm version**; explicit invalidation + hit metrics. |
| 7-table Y-profile schema + manual concordance | Simplify to a `YDnaProfile { snps, strs, sources, reconciliation }` aggregate. |

"Port as-is" (logic is sound, just translate): haplogroup **scoring** (tree DP traversal),
external clients (ENA, facility, tree providers), IBD detection/relationship math.
"Redesign": workspace state model, persistence, sync, UI, analysis-processor interface
(define one `Processor` trait with `init/run/cleanup` + standard progress signature).

---

## 7. AT Protocol / OAuth

Reuse `du-atproto` — but **partially**. The primitives reuse directly: DID/handle
resolution, PDS discovery, PKCE (S256), DPoP proofs. The confidential-client pieces it also
implements (`private_key_jwt` ES256 client assertion, served client-metadata/JWKS, cookie
session) are for the **decodingus web** client and do **not** apply here: Navigator is a
desktop app, hence a **public/native client** — PKCE only, its own native
`client-metadata.json`, a **loopback redirect** (`http://127.0.0.1:<port>/callback`), and
tokens in the **OS keychain**. Confirm `du-atproto`'s token-exchange builder runs
PKCE-without-client-assertion (small add if not). See
`documents/atmosphere/11-Auth-and-Permissions.md` §6–7 and the server-side companion
`decodingus/rust/docs/atproto-oauth-findings.md`.

This is the OAuth migration the design review called for: app-password `createSession` is
gone; Navigator authenticates via OAuth and requests the `navigatorCore` **write** scope.

**AppView scope reduction (2026-06) shifts work to Navigator.** The AppView no longer
mirrors the network (see `documents/atmosphere/08-AppView-Lifecycle.md`); per-sample data is
authoritative in the Navigator workspace. This adds two Navigator responsibilities:
- **Publish per-sample coverage summaries** as public PDS records (under `navigatorCore`)
  so the AppView can aggregate them on demand. Coverage reporting for the researcher's own
  cohort is also a Navigator-local feature over the local workspace.
- **A variant-proposal submission client** — Navigator posts variant/branch proposals
  directly to the AppView curation API (decoupled from PDS records), replacing the old
  firehose-harvest of private variants.

These live in `navigator-sync` (publishing) and a small `navigator-app` command (proposal
submission). The desktop app remains a **writer** (direct-to-PDS) plus a reader of its own
AppView records; no broad PDS read scope is needed at this stage.

---

## 8. Data migration (H2 → SQLite)

Follow the `du-migrate` precedent: a one-time CLI (`navigator-migrate`) reads the existing
`~/.decodingus/data/workspace.mv.db` (H2) and writes the new SQLite schema. Map
biosamples/projects/sequence-runs/alignments/haplogroups; re-derive anything cheap rather
than migrating it (analysis caches can be recomputed). Keep the H2 file untouched as a
rollback.

---

## 9. Phased roadmap

1. **Extract shared crates.** Promote `du-domain`/`du-atproto`/`du-bio` to the shared
   location; stand up the `navigator` workspace skeleton with the §2 crates compiling empty.
2. **Raw-read layer + finish du-bio.** Finish shared `du-bio` `liftover`/`callable`/`vcf`
   (coordinate math + text parsing). In **`navigator-analysis`**, implement the `noodles`
   BAM/CRAM/FASTA/BGZF/index I/O layer and port the **walkers** (§4a) on top of it, with
   unit tests.
3. **Caller + parity harness.** Build the haploid caller (§4b) and the golden-truth
   validation harness (§4c). Gate on parity. (Sidecar bridge §4d available if needed.)
4. **Store + app layer.** SQLite schema + migrations; `navigator-app` command/query layer;
   `navigator-migrate` from H2.
5. **UI.** egui shell → dashboard, subject detail, haplogroup/coverage/ancestry views,
   reports. Replace dialogs with app-layer policy + thin prompts.
6. **Sync + OAuth.** Wire `du-atproto` OAuth; complete `navigator-sync` (the part that is
   stubbed today).
7. **IBD + Y-profile.** Port detection/relationship math; simplified Y-profile aggregate.
8. **Cutover.** Feature-parity check against the golden harness; ship.

---

## 10. Risks & open questions

**Risks**
- **Indel/homopolymer calling parity** (§4b) — the single biggest technical risk; owned by
  the §4c harness. May constrain v1 branch-creation to SNP evidence.
- **Numerical drift** in metrics/coverage vs GATK — covered by the same harness.
- **Shared-crate coordination** — extraction adds a release axis between two teams/repos.
- **egui polish ceiling** — acceptable for a workbench; revisit only if a report view needs
  rich native widgets.

**Open questions**
1. ~~Where do the shared crates live?~~ **Resolved:** `decodingus-shared`, a sibling repo at
   the Development root (DUNavigator/decodingus/decodingus-shared as siblings). Remaining
   sub-question: path deps vs git deps vs published registry for CI.
2. Does Navigator share the **same SQLite schema shape** as the decodingus Postgres schema
   where domains overlap (biosample/variant/haplogroup), or only the domain *types*?
3. Long-read specifics (PacBio/ONT) the walkers must preserve (e.g. `COUNT_UNPAIRED`,
   read-length handling) — enumerate so parity tests cover them.
4. mtDNA indel fallback: is a pinned external caller acceptable *temporarily* for
   validation, or must everything be pure-Rust from day one?
