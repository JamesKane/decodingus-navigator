# Pangenome / GAM as Informative Data Sources — Forward-Looking Design

**Status:** Design only, forward-looking (v0, 2026-06-06). Horizon: the "eventual
landing" of the human pangenome as a working reference, not a near-term build.
**Scope:** the **Edge/Navigator** pipeline that turns graph alignments (GAM/GAF)
into the informative records Navigator and AppView already model. Pairs with the
**existing AppView storage model** (`decodingus/rust/migrations/0004_genomics.sql`:
`pangenome_graph`, `pangenome_path`, `canonical_pangenome_variant`,
`pangenome_alignment_metadata`, `reported_variant_pangenome`) and the graph
coordinates already in `decodingus/documents/schema/universal-variant-schema.md`.

> **Framing of the gap (the ask):** the *storage* side is modeled; the *producer*
> isn't. Today a sample's reads become informative data only through the **linear**
> path (BAM/CRAM on CHM13/GRCh38 + liftover). Nothing turns a **graph** alignment
> into Navigator's variants / haplogroups / ancestry / coverage / STR. This doc is
> about that conversion — and, deliberately, about *which signals the graph makes
> available that the linear path cannot*.

## 1. Why the pangenome matters for *this* product

Navigator's whole current stack assumes **one linear reference + liftover** (the
`refgenome`/liftover phase, CHM13↔GRCh38, Y rev-comp tracts, mtDNA rotation). The
pangenome is a different paradigm: a **graph** of many haplotypes where linear
builds are just *paths* through it. For a genetic-genealogy product specifically,
the payoffs are concentrated exactly where we work:

- **Reference-bias elimination** — reads align to the *best-matching haplotype*,
  not a single ref. Wins in **divergent single-copy sequence and novel insertions**
  (SNV-dense regions absent or mismatched in a single linear ref). **NOT** in large
  near-identical repeats — see §1.1.
- **Native structural & novel sequence** — insertions/CNV/segmental duplications
  and **sequence absent from any linear ref** are first-class graph objects, not
  liftover casualties. Directly feeds **private-variant / haplogroup discovery**.
- **Haplotype-resolved signal** — the graph *is* a panel of real haplotypes
  (HPRC's 47+ diploid assemblies, growing). Which paths a sample traverses is a
  direct **ancestry / local-ancestry / phasing** signal, richer than projecting
  onto a single ref's AIMs.
- **Liftover obsolescence (eventually)** — the graph is the multi-reference
  coordinate system. Long-term, graph coordinates *subsume* the liftover machinery
  rather than adding to it.

### 1.1 Where the graph helps vs. where it hurts — Y loops are graph-hostile (empirical)

> **Empirical finding (own experiments, summer 2025):** the **Y-chromosome loops**
> — the palindromes (P1–P8) and ampliconic repeats, megabase-scale inverted/direct
> repeats at >99.9% identity — are **difficult-to-impossible to use directly in a
> GAM.** In the graph these collapse into **cyclic/tangled subgraphs**; a short read
> has no unique anchor, so its traversal is ambiguous and the alignment is
> unusable there. **Surjection to a linear path was required to get anything
> usable.** This is a graph *failure mode*, not a corner case.

Consequences — bake this into the design rather than discovering it twice:

- The graph's value is in **divergent single-copy sequence + novel insertions**.
  **Large near-identical repeats** (Y palindromes/amplicons, and by extension other
  segmental duplications) are *worse* in a graph than on a linear ref, because the
  repeat collapses into loops with no unique placement.
- **For the Y specifically** (our haplogroup core): the genealogically informative
  SNPs largely sit in the **X-degenerate single-copy** male-specific region — those
  are graph-friendly. The **ampliconic/palindromic** regions are **surject-only**.
  So Y processing is *inherently hybrid*: native graph signal where single-copy,
  **mandatory surjection where looped**. (FTDNA's Big Y BED already excludes most of
  the ampliconic mess — a useful prior for which regions to even attempt natively.)
- **Surjection is therefore not merely a convenience** to reuse the linear pipeline
  (as §4 first framed it) — it is **required for correctness** in the repeat/loop
  regions. The design must classify regions into *graph-native-safe* vs
  *surject-only* and route accordingly (§7).

## 2. Tooling reality (verify at implementation; field moves fast)

- **Graph artifacts:** HPRC **Data Release 2** (May 2025), graphs **v1.1**
  (Minigraph-Cactus + pggb), distributed as **GFA / GBZ / Giraffe indexes / VCF /
  odgi**. GBZ (GBWT + graph) is the compact production form.
- **Pure-Rust graph access — yes:** `gbwt-rs` (jltsiren) reads **GBZ/GBWT** in
  Rust → Navigator can hold/traverse the graph and its haplotype paths **without
  htslib/C++**, consistent with the noodles-only stance.
- **Graph alignment — external C++:** the fast short-read graph mapper is
  **vg Giraffe** (C++); there is no pure-Rust equivalent. So *producing* the
  alignment is an **external tool** (vg Giraffe → **GAM**, or GraphAligner →
  **GAF** for long reads). This is the same shape as our existing reliance on
  external alignment for linear BAMs.
- **Reading alignments in Rust:**
  - **GAF** = plain TSV (rGFA path strings) → trivial pure-Rust parse.
  - **GAM** = vg protobuf stream → `prost` + vg's `.proto`, or sidestep via
    `vg surject` (GAM→BAM on a chosen reference path) so the **existing noodles
    linear pipeline just works**.
- **Graph genotyping/coverage (external, today):** `vg pack` (packed coverage over
  nodes/edges → `.pack`) and `vg call` (snarl genotyping → VCF). Pure-Rust
  re-implementation of these is a *later* option, not a prerequisite.

**Net constraint:** Navigator does **not** align to the graph itself. It consumes
graph-aligner *output*. The design question is **how much it extracts natively
(snarl traversal, path support, node coverage) vs. how much it lets vg flatten to
VCF** — that ratio *is* the "informative data sources" question.

## 3. The informative-signal taxonomy (the core)

What a graph alignment yields, what it adds over the linear path, how to extract
it, and where it lands in the **existing** model. Ordered by value-for-effort.

| # | Signal | Beyond linear | Extraction | Lands in |
| --- | --- | --- | --- | --- |
| S1 | **Reference-path variant calls (via surjection)** | fewer FN in divergent single-copy regions; **and the only usable signal in Y loops** (§1.1) | `vg surject`→BAM→existing caller, **or** `vg call`→VCF | existing `VariantSet` / `variant_call` (Navigator); `reported_variant_pangenome` (AppView) — drop-in |
| S2 | **Snarl (bubble) genotypes** | SVs + multi-allelic + novel alleles as first-class sites; no liftover | snarl decomposition + traversal support (`vg call`, or native via `gbwt-rs` + `.pack`) | `reported_variant_pangenome` (`variant_nodes[]`, `variant_edges[]`, `zygosity`, `allele_fraction`) — table already shaped for this |
| S3 | **Node/path coverage** | presence/absence of segments; CNV; true callable footprint per haplotype | `vg pack` → per-node depth; aggregate to PATH/NODE/REGION | `pangenome_alignment_metadata` (`metric_level` GRAPH/PATH/NODE/REGION) — already modeled |
| S4 | **Haplotype-path support → ancestry** | which real HPRC/1000G haplotypes the reads resemble; *graph-native* ancestry + local ancestry | GBWT path support along traversals (`gbwt-rs`) → per-region nearest-haplotype profile | Navigator `ancestry` (new graph-native method alongside ADMIXTURE/PCA/G25); `haplotype_information` JSONB |
| S5 | **Y / mtDNA haplogroup, natively (single-copy only)** | derived/ancestral off node traversal; no liftover — **but only in X-degenerate single-copy Y; ampliconic markers fall back to surject (§1.1)** | encode phylo markers as graph nodes (single-copy regions) → traversal = call; loop-region markers via S1 | Navigator `haplo.rs` → `haplogroup_call` (new `source_key="graph:<graph>"`, distinct from linear) |
| S6 | **Novel-sequence / private variants** | reads supporting off-reference nodes or new insertions = candidate private SNPs/branches | non-reference node support above threshold within callable footprint | private-variant feed → AppView **haplogroup-discovery** system |
| S7 | **STR / VNTR genotypes** | repeat structure is explicit in the graph; count traversals | snarl-local repeat-unit traversal counting | `StrProfile` (WGS-derived), `ystr` modal aggregation |

**Reading of the taxonomy:** S1 is a near-free win (reuse the whole existing
linear pipeline via surject). **S2–S3 are the first genuinely graph-native records**
and the tables already exist for them. **S4–S5 are the differentiated payoff** for a
genealogy product — ancestry and haplogroup *without* liftover, from real
haplotypes. S6–S7 are upside that connects to systems we're already building
(discovery, STR).

## 4. Pipeline architecture — three postures

| | A. Flatten (vg→VCF) | B. Native graph-genotyping | C. Hybrid (recommended) |
| --- | --- | --- | --- |
| Align | vg Giraffe (ext) | vg Giraffe (ext) | vg Giraffe (ext) |
| Genotype/cover | `vg call`/`vg pack` → VCF/coverage | Navigator reads GAM/`.pack` + GBZ, derives natively | vg for align + **`vg pack`** (coverage); Navigator derives S2/S4/S5 from `.pack` + GBZ |
| New Rust | almost none (reuse VCF import) | a lot (snarl logic, path support) | moderate (`gbwt-rs` + `.pack` reader + snarl traversal) |
| Signals | S1 only (linear-projected) | S1–S7 | S1 (surject) + S2–S5 now, S6–S7 later |
| Verdict | quickest, loses the point | maximal, heavy, premature | **staged: cheap S1 first, graph-native next** |

**Recommendation — C, staged *and region-gated*.** Begin by treating the graph as
"a better way to get the *same* records" (surject + `vg call` → existing
`VariantSet`/coverage) — this de-risks ingestion, benefits divergent single-copy
calling, **and is the only usable path through the Y loops (§1.1)**. Then add the
**native extractors** (S2/S3 from `.pack`+GBZ via `gbwt-rs`, then S4/S5) **scoped to
graph-native-safe regions only**, with repeat/loop regions permanently routed
through surjection. This keeps Navigator pure-Rust for everything except the
alignment/`pack` steps (already an external-tool boundary). The native path is an
*enhancement over* a surjection backbone, never a replacement for it.

## 5. Coordinate & identity implications

- **Coexistence, not replacement (at first).** The graph is *another reference*
  alongside CHM13/GRCh38. `universal-variant-schema.md` already carries graph
  coordinates (`{type: pangenome, graph, node}`); a variant can hold *both* a
  linear coord and a node-set, joined by `vg`-deconstruct against the reference
  path. **Don't rip out liftover** — run graph and linear side by side, reconcile
  (existing reconciliation handles multi-source), and let the graph *earn* primacy.
- **The long game.** Once graph calls are trusted, the graph becomes the **primary
  coordinate authority** and linear builds become reference *paths* through it —
  liftover collapses into "project to a path." That is the actual "landing."
- **Graph identity & versioning.** Pin every record to a specific
  `pangenome_graph` (name + GBZ checksum + HPRC version) — graph topology changes
  between releases (v1.0→v1.1 changed even the path-name separator). Node IDs are
  **not stable across graph versions**; store the graph version with any
  `variant_nodes[]`, and re-derive on graph upgrade rather than assuming stability.
- **Federation.** Graph-native *aggregate* records (ancestry, coverage, haplogroup)
  federate fine via `du-domain::fed`. Per-node genotypes are higher-resolution PII
  risk than linear — keep federated outputs aggregate, same as today.

## 6. Staging / roadmap

- **P0 — exists:** AppView storage model + graph coordinates in the schema.
- **P1 — graph ingest + surject (S1):** register a `pangenome_graph` (GBZ via
  `gbwt-rs`); accept a GAM/GAF; `vg surject`→BAM so the **existing** caller/coverage
  path runs; write `reported_variant_pangenome` + `pangenome_alignment_metadata`.
  *Smallest slice; immediate divergent-region benefit.*
- **P2 — native snarl genotypes + node coverage (S2/S3):** `.pack` reader + snarl
  traversal in Rust (or `vg call` to start) → first graph-native records.
- **P3 — graph-native ancestry + haplogroup (S4/S5):** GBWT path-support profiling
  → new ancestry method; phylo-marker-node traversal → `haplogroup_call`. *The
  differentiated payoff; removes liftover from the Y/mt path.*
- **P4 — discovery + STR (S6/S7)** and **graph-as-coordinate-authority** — the
  liftover-subsuming end state.

du-bio (shared) is the natural home for graph coordinate math / snarl logic so both
Navigator and AppView reason about node-sets consistently.

## 7. Open questions / decisions

1. **External `vg` dependency — accept it?** Navigator is pure-Rust by stance, but
   graph alignment (Giraffe) and `pack`/`call` are C++. Confirm `vg` as an external
   tool (like the aligner already is) vs. waiting for a Rust mapper (not on the
   horizon). Gates everything.
2. **GAM vs GAF as the ingest format** — GAF is pure-Rust-trivial but is the
   long-read/GraphAligner format; Giraffe emits GAM (protobuf). Surject sidesteps
   both for S1. Decide per signal.
3. **Which graph(s)** — HPRC v1.1 whole-genome is large; do we ship a **Y+mt
   sub-graph** first (small, directly serves the haplogroup core) before
   whole-genome? Likely yes — **but it must mask/exclude the ampliconic loops**
   (§1.1), which the graph can't serve; those stay surject-only.
3a. **Region classification (gating the native path)** — build a per-graph map of
   *graph-native-safe* (single-copy, well-behaved) vs *surject-only* (palindromes,
   ampliconic, segmental dups) regions. Seed from the Y X-degenerate definition +
   FTDNA Big Y BED. Native extractors (S2/S4/S5) run only inside the safe set;
   everything else falls back to S1 surjection. **This is the load-bearing decision
   the Y-loop finding forces.**
4. **Node-ID instability across versions** — re-derive vs. map. Pick a policy before
   storing any `variant_nodes[]` at scale.
5. **Where does this run** — same Edge box as the current GATK/noodles pipeline;
   confirm `vg` fits the resource/packaging model (it's heavyweight).
6. **Priority vs. current roadmap** — this is a *horizon* design; it should not pull
   focus from FTDNA import, ancestry, or the AppView launch. Sequence it after.

## 8. Next step

Confirm Q1 (accept `vg` as an external tool) and Q3 (Y+mt sub-graph first, ampliconic
loops masked), then the P1 slice (graph register + **surject** → existing pipeline)
is a self-contained, low-risk proof — and surjection is now known to be the
*durable backbone* (§1.1), not just a bootstrap. The native extractors (S2/S4/S5)
are an enhancement layered on top, gated by the region map (Q3a). The signal
taxonomy (§3) is the part to pressure-test against what we actually want to report.
