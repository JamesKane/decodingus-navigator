# DecodingUs haplogroup tree provider — plan

Last updated: 2026-06-10. Repos: **DUNavigator** (`rust-rewrite`) + **decodingus** AppView
(`rust-rewrite-foundation`). Both are ours. Source of truth for the DecodingUs tree is the
**Rust AppView**, not the live `decoding-us.com` site.

## Goal

Let Navigator place Y-DNA haplogroups against the **DecodingUs** tree (served by our AppView)
in addition to FTDNA. Today Y/mt placement is hardcoded to FTDNA endpoints
(`lib.rs:1651-1658`), even though the parser/scorer (`haplo::parse_ftdna_json`, `haplo::score`)
are generic. DecodingUs is *our* tree — FTDNA-only is a real limitation for Y work.

## The decisive advantage: native CHM13 coordinates

- The FTDNA Y tree is **GRCh38**; for a CHM13 alignment Navigator must **lift every tree
  position** GRCh38→CHM13 via a chain (`lifted_targets`, `tree_build_for_contig → "GRCh38"`).
  That liftover is the costly, fiddly part of Y placement.
- The AppView stores variant coordinates as `Coordinates(BTreeMap<build, BuildCoordinate>)`
  (shared `du-domain::variant`), keyed `"GRCh38"` / `"GRCh37"` / **`"hs1"` (T2T-CHM13)**.
  So the DecodingUs tree can hand Navigator **native CHM13 (`hs1`) positions** — for a CHM13
  alignment, **no liftover at all**: query the tree positions directly.

## Current AppView state (gap to close)

- `GET /api/v1/y-tree` (`du-web/src/api.rs:527`, public) returns nested `HaplogroupNodeDto`
  `{id, name, haplogroup_type, formed_ybp, tmrca_ybp, children}` — **no variants/coordinates**.
- Variants are per-haplogroup only: `GET /api/v1/haplogroups/{name}/variants` → `Vec<VariantDto>`
  with `coordinates` (multi-build). Per-node fetch across the whole tree = thousands of calls →
  infeasible for building a placement tree.
- **Gap:** no single payload with tree **+** each node's defining variants/coordinates.

## Plan

### A. AppView (decodingus) — add a tree-with-variants endpoint

`GET /api/v1/y-tree/full` (and `/api/v1/mt-tree/full`), public, same `RootParams` as today.
Returns the nested tree with each node's defining variants embedded:

- Extend `HaplogroupNodeDto` with `#[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub variants: Vec<VariantDto>` — empty (omitted) on the existing `/y-tree`, populated on
  `/full`. Backward compatible.
- New `du-db` query: fetch all `(haplogroup_id, variant)` for the dna_type (optional root
  subtree) in **one** round trip — join `tree.haplogroup_variant → core.variant`, current
  edges only (`valid_until IS NULL`) — and group by haplogroup_id. Assemble into the nested
  DTO alongside the existing `subtree` nodes.
- Wire route in the public `/api/v1` router; add to the OpenAPI `paths(...)`.

Payload size: the catalogued Y variants (not 113 MB like FTDNA's raw tree) — acceptable as one
download, cached on the Navigator side exactly like the FTDNA JSON is today.

### B. Navigator — provider abstraction + DecodingUs parser + config

1. **Parser** `haplo::parse_decodingus_json(data, target_build) -> Result<HaploTree>`: maps the
   AppView's `y-tree/full` JSON into the existing `HaploTree`/`HaploNode`/`Locus`. Per variant,
   pick the `target_build` coordinate (`"hs1"` for CHM13 alignments, else `"GRCh38"`); drop
   variants lacking that build. Nested `children` → the existing id-keyed `children: Vec<i64>`
   (assign ids by traversal). Same `Locus{position, ancestral, derived, name}` the scorer uses.
2. **Provider abstraction**: a small `enum YTreeProvider { Ftdna, DecodingUs }` (or trait) that
   yields (fetch URL, cache key, parse fn, **native build per contig**). FTDNA path unchanged.
3. **Native build → skip liftover**: `tree_build_for_contig` becomes provider-aware. DecodingUs
   on a CHM13 alignment returns the build the tree was parsed in (`hs1`/CHM13) == the alignment
   build → `lifted_targets` returns `None` → direct query, no chain. (FTDNA stays GRCh38 + lift.)
4. **Configurable host** (per user: local for testing, switchable): an AppView base URL resolved
   from env `DECODINGUS_APPVIEW_URL` (default `http://localhost:8080`), used to build the
   `/api/v1/y-tree/full` URL. The fetch+cache reuses the existing `fetch_tree` machinery
   (on-disk cache under `~/.decodingus/trees`, keyed `decodingus-ytree`).
5. **Provider selection**: env `NAVIGATOR_Y_TREE_PROVIDER = decodingus | ftdna`. Default —
   **see decision Q2**. `assign_y_haplogroup` / `private_y_variants` pick fetch+parse by provider.

### Scope

- **Y-DNA only.** DecodingUs has no mtDNA tree (the Scala provider threw on MTDNA); mtDNA stays
  FTDNA. The provider enum reflects that (DecodingUs = Y only).
- mt placement, the FTDNA path, and all scoring/placement logic are untouched.

## Validation

- AppView: a `du-web` test hitting `/api/v1/y-tree/full` returns nodes with non-empty
  `variants[].coordinates` incl. an `hs1` entry for a known CHM13-mapped SNP; the plain
  `/y-tree` still omits variants.
- Navigator: unit test `parse_decodingus_json` on a small fixture → expected `HaploTree`.
- Live: run the AppView locally (`DECODINGUS_APPVIEW_URL=http://localhost:<port>`), place the
  GFX0457637 CHM13 sample via the DecodingUs provider, assert it still yields **R-FGC29071**
  (matching the FTDNA result) — and confirm it took the **no-liftover** direct-query path.
- `cargo test`/`cargo build` green + warning-free in both repos.

## Decisions to confirm

- **Q1 — AppView endpoint shape**: new `/api/v1/y-tree/full` (recommended; clean, explicit) vs
  extend `/api/v1/y-tree?variants=true`.
- **Q2 — Default Y provider**: make **DecodingUs the default** (native CHM13, no liftover, our
  tree) with FTDNA as fallback, vs keep **FTDNA default** and DecodingUs opt-in.
- **Q3 — Config mechanism**: env vars now (MVP; no config system exists yet — gap #4) vs a small
  persisted setting + later a Settings UI.
