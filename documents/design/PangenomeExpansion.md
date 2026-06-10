# Pangenome Expansion — Design Notes

Status: **design notes / not yet scheduled**. Drafted 2026-06-09; updated
2026-06-10 with measured build-footprint numbers (see "Empirical resource
footprint").
Source: the personal-pangenome Y-haplogroup research (lab notebook at
`~/Genomics/personal/Y_haplogroup_experiment_notes.md`; blog write-up
"Finding a Y Haplogroup in a Pangenome (and Why It Almost Didn't Work)").

These are notes on what needs to be **designed** before the Navigator gains any
pangenome-aware capability. They build on the existing linear BAM/CRAM stack
(`navigator-analysis`, region-aware Y profile, multi-source reconciliation,
CHM13 liftover) rather than replacing it. Several findings here are pure
upgrades that need no graph at all.

---

## Implementation status against the Rust rewrite (assessed 2026-06-10)

Legend: 🟢 **here** (built, in the current codebase) · 🟡 **planned/partial**
(scaffolding or adjacent capability exists; gap is the specific algorithm/field)
· 🔴 **needs new design** (nothing in the codebase yet).

| Item | Status | Evidence in tree / gap |
|---|---|---|
| Path-supported parsimony placement | 🟢 | **Built + validated 2026-06-10** (`haplo::path_admissible`, wired at `assemble_assignment`). Shipped as a *guard over the Kulczynski ranking*, not a standalone descent: report the best-ranked candidate whose root→node lineage crosses no **contradicted** branch (sample net-ancestral, `a > d`; no-calls/stray errors pass). Fixes the distal-Y tunnel artifact without disturbing the proportional placement. **Validated on GFX0457637/CHM13 → Y `R-FGC29071`, mt `U5a1b1g`** (both correct). *Note:* a descent router ("follow the most-derived subtree") was tried first and **regressed** GFX to a bushier wrong fork (`R-Z17665`) — absolute derived count favours long/bushy paths; reverted. Catching paralog *false-derived* calls (not honest ancestral ones) still needs the allele-balance filter below. |
| Marker-less node collapse (recursive splice of indel-only nodes) | 🟡 | Handled inline by parsimony descent (SNP-less nodes are pass-through via look-ahead). A standalone recursive tree-rewrite/re-parent pass is still net-new but no longer blocks placement. |
| Haploid allele-balance / paralog filter (near-monoallelic gate) | 🟢 | **Built + validated 2026-06-10.** `caller::is_paralogous` — drops a haploid call when the *second* allele has both `>= min_paralog_minor_reads` reads AND a fraction `> max_minor_allele_fraction` (defaults 2 / 0.20); a lone error read doesn't trip it. Applied at all three haploid call sites (`call_bases_at`, `genotype_sites`, `call_denovo`); diploid ancestry path untouched. GFX0457637/CHM13 holds (Y `R-FGC29071`, mt `U5a1b1g`; Y score 0.546→0.548). |
| Callable-region / paralog mask (BED) infrastructure | 🟢 | `mask.rs` `RegionMask::from_bed` — sorted/coalesced interval mask, binary-search `contains`. Built. |
| Mask wired into Y profile / region-aware quality | 🔴 | `RegionMask` has **no consumer** outside its own module (`mask` is exported but unused by `haplo`/`caller`/`unified`). Wiring is net-new. |
| Canonical `chm13v2.0_maskedY_rCRS` analysis reference | 🔴 | `registry.rs` `Build` enum is `Grch38`/`Grch37`/`Chm13v2` only; reference is plain `chm13v2.0.fa`. No masked+rCRS analysis-set variant. |
| Reference-polarity metadata + guard (CHM13 Y = HG002 = J) | 🔴 | No polarity field/guard anywhere in `navigator-refgenome`. Net-new. |
| rCRS mtDNA path (rotation-aware) | 🟢 | `mtvariants.rs`/`heteroplasmy.rs` already operate on rCRS; rCRS↔chrM handled (memory: rotation-aware map validated). mtDNA-from-pangenome correctly parked. |
| Multi-run haplogroup reconciliation records | 🟢 | `navigator-sync/records.rs`: `HaplogroupReconciliationRecord`, `RunHaplogroupCallRecord` (`source_ref`, `call_method`, `confidence`, scores). |
| Callset **source/completeness** dimension (linear / surjected / graph-genotype-only) + reconciliation guard | 🔴 | No completeness or analysis-source-type field on the records; reconciliation cannot down-rank a panel-conditioned "absent". Net-new schema + guard. |
| Surjected-linear-CRAM ingest with provenance label | 🔴 | No provenance/source-type on alignment ingest (`probe.rs` infers build/aligner/platform, not graph-surjection). Net-new. |
| Local `vg` graph mapping (option B) | 🔴 | Decision-gated, out of scope by the doc's own recommendation (option A). |
| SV "graph-coherent vs linear-fragmented" tagging | 🔴 | `sv` exists but output unvalidated per HANDOFF; the distinction is a future field. |

**One-line read:** the *linear-first substrate* the doc depends on is largely
🟢 here (placement scaffolding, mask loader, rCRS path, reconciliation records).
Phase 1's core algorithms are now 🟢 built + validated (parsimony placement guard,
paralog allele-balance filter); only standalone marker-less collapse remains 🟡.
Everything touching the **graph/provenance/reference-set/polarity
schema** (Phases 2–3) is 🔴 net-new design.

---

## Motivation: what the research established

I built a small personal pangenome (T2T-CHM13 plus eight phased assemblies),
aligned a known R1b sample to it with `vg giraffe`, and tried to recover the
terminal Y haplogroup (R1b-FGC29071) from the graph callset. The outcome was
instructive and directly relevant to this app:

1. **A genotype-only graph callset reports only variation already in the graph.**
   `vg call` (without augmentation) can only emit alleles carried by one of the
   panel assemblies. The sample's terminal lineage was not in the panel, so it
   was simply not callable. This is genome-wide, not Y-specific. It bites every
   rare/private variant.
2. **The omission is silent and ambiguous.** In a sites-only graph VCF, "no
   record" conflates three different states: matches-reference, allele-absent-
   from-graph, and no-coverage. A linear gVCF keeps those distinct. So a graph
   callset carries *less* per-site information than a linear one.
3. **Discovery needs surjection or augmentation.** Surjecting the reads back to
   linear CHM13 and calling a pileup recovered the terminal (the marker that had
   no graph record had eleven reads, all derived). Augmenting the graph with the
   reads also recovered it. Genotype-only graph calling did not, and could not.
4. **A naive "deepest derived node" placement is actively misleading** in the
   repeat-rich distal Y, where paralog/mismapping artifacts produce clustered
   false positives. A path-supported parsimony placement is required.
5. **Reference polarity matters.** CHM13v2.0's chrY is HG002's Y, which is
   haplogroup J. So the reference allele on chrY is *not* universally ancestral;
   derived/ancestral must come from the tree, never from "is it REF or ALT."

The honest framing: pangenome mapping reduces reference-allele bias and helps for
known/common and structural variation, but for discovery you pay a silent
panel-representation bias unless you surject or augment. **For an edge/privacy app
that already operates on linear BAM/CRAM, this is mostly good news: the linear
path we already have is the correct place to do discovery.**

---

## Core design principles (derived from the above)

1. **Linear-first stays correct.** The Navigator should keep doing discovery
   (haplogroups, private SNPs, terminal-Y, coverage, IBD, ancestry) on linear
   BAM/CRAM. The pangenome's role is *upstream mapping quality*, consumed as a
   surjected linear BAM/CRAM, not a new analysis-time dependency.
2. **Never treat a graph-derived callset as complete.** If we ever ingest a
   genotype-only graph VCF, it must be tagged as panel-conditioned, and "absent"
   must never be read as "ancestral/reference." This is a schema and a
   reconciliation concern, not just an analysis one.
3. **Harden the placement core** with the path-supported parsimony algorithm
   (below). This is a pure improvement and needs no pangenome.
4. **Reference and region correctness** are first-class: ref-allele polarity,
   rCRS for mtDNA, Y-PAR masking, and callable-region masks for the cyclic
   distal Y.

---

## The hard tension: pure-Rust runtime vs graph alignment

The Rust rewrite's whole premise is "no JVM, no external bioinformatics tooling"
(noodles, not vg/samtools/bcftools). Graph alignment (`vg giraffe`) and
surjection are C++ in `vg`; there is no mature pure-Rust graph aligner. So a
literal "embed a pangenome aligner in the edge app" is off the table without
abandoning that principle.

Options to decide between:

- **(A) Upstream/offline mapping.** Graph alignment + surjection happens once,
  outside the edge app (a lab pipeline, a server, or a one-time tool), and the
  user feeds the Navigator the resulting **surjected linear CRAM**. The runtime
  stays pure Rust. Recommended default.
- **(B) Optional bundled `vg` pre-step.** Ship `vg` as an optional, clearly
  separate pre-processing helper (not the analysis runtime) for users who want
  to map their own FASTQ to a pangenome locally. Compromises purity for one
  opt-in feature; keeps analysis pure.
- **(C) Pure-Rust graph alignment.** Large, speculative effort. Not now.

The research argues we mostly do not need the graph at analysis time at all. So
the cheapest correct posture is **(A)/(C-consume)**: the Navigator should *accept
and prefer* pangenome-surjected CRAMs (label their provenance), and we punt graph
alignment to upstream tooling. Embedding `vg` (B) is a later, decision-gated
convenience, not a requirement.

---

## Empirical resource footprint (measured from the build)

Disk-usage findings from the actual personal-pangenome build on this machine
(`~/Genomics/personal/`, measured 2026-06-10). These put hard numbers on the
"graph is heavy, keep it upstream" decisions above.

- **Total experiment footprint: ~241 GB** for a *single* nine-haplotype personal
  pangenome (T2T-CHM13 + eight phased assemblies). That is the cost of one panel;
  it does not scale down.
- **`vg giraffe` map-time index set** (what the aligner must load):
  - `*.gbz` (graph) ~1.9 GB
  - `*.min` (minimizer) **~11 GB hifi / ~35–44 GB short-read** — the minimizer
    index dominates and is read-type-specific (short-read is ~4× the hifi one).
  - `*.dist` (distance) ~1.0 GB, `*.snarls` ~40 MB, `*.zipcodes` 0.1–0.45 GB
  - So a usable giraffe index is **~15 GB (hifi) to ~47 GB (short-read)** on its
    own — well past anything we would ship to, or build on, an edge device.
- **Augmentation / full-GAM discovery is the heavy path, now quantified.**
  `combined_alignments.gam.gz` is **~110 GB** for one sample — the single largest
  artifact, and concrete evidence for Open Question #3: augment-based discovery is
  decisively *not* an edge operation, and is a non-trivial server one.
- **Each panel rebuild is a fresh full index set, not an incremental delta.** The
  2025 build and the 2026 short-read/hifi "rebuilt" variants are *distinct* files
  — verified: same-sized `.min`/`.dist` across builds hash differently, so there is
  no byte-level reuse between revisions. Iterating the panel multiplies the
  tens-of-GB cost per revision; index storage/versioning is a real upstream concern.

**Implication.** These numbers reinforce **option A (graph mapping strictly
upstream)** and the surjected-linear-CRAM consumption model. The edge app's input —
a surjected CRAM — is single-digit-to-low-tens of GB and shrinks further as CRAM,
while the graph index and GAM that produced it (tens to >100 GB) never touch the
device. Option B (bundled local `vg`) would require shipping or building one of
these index sets locally, which the short-read footprint makes clearly impractical.

---

## Module-by-module design needs (Rust crates)

### `navigator-analysis` — haplogroups (the biggest concrete win)

Port the path-supported parsimony placement validated in the research. It
complements, not replaces, the existing region-aware quality and multi-source
reconciliation:

- **Gate:** only descend into a child branch with net positive derived support
  (derived calls plus reference-match-consistent sites, minus contradictions).
  Forbids tunnelling through unsupported branches, so isolated deep paralog
  artifacts are unreachable.
- **Route:** among traversable children, follow the one whose subtree holds the
  most positive variant evidence (look-ahead). Keeps the ancient backbone on
  track where an intermediate branch is only weakly supported.
- **Report:** the deepest node with genuine variant evidence; trim a
  reference-match-only tail so a lone coincidental match cannot stretch a call.
- **Collapse marker-less nodes (recursive).** Tree nodes defined only by indels
  (e.g. R1b-A353) are invisible to a SNP caller and block a strict descent.
  Splice them out, re-parenting children to the nearest SNP-bearing ancestor.
  Sample-independent, so safe on sparse and dense data alike.
- **Haploid allele-balance filter.** A true Y/mt-haploid site is near-monoallelic;
  mixed VAF signals paralog/mismapping and should be dropped. This is the cheap,
  principled defense against the distal-Y artifacts.
- **Conflict tolerance, opt-in, dense data only.** Allow a lone contradicting
  call when downstream support overwhelms it. Off by default (keeps sparse-input
  placement conservative).

Notes:
- This logic likely belongs in shared **`du-bio`** so the server and the edge app
  place identically; the Navigator's haplogroups module calls into it.
- It should run uniformly over the input modes we already have: a VCF, or a
  pileup over a linear BAM/CRAM. The reference implementation did exactly this.
- Ties into the existing **region-aware quality**: region reliability and the
  allele-balance/paralog filter are two views of the same "is this call trustworthy
  in a repeat" question. Unify them.

### `navigator-analysis` — sv

SV/large-variant representation is where pangenome-surjected inputs and the graph
genuinely help (coherent structural alleles vs fragmented small-variant calls).
SV output is still unvalidated per HANDOFF; when validated, note the
"graph-coherent vs linear-fragmented" distinction for any pangenome-sourced input.

### `navigator-refgenome`

- **Canonical analysis reference.** Per `docs/chm13-reference-resources.md`, prefer
  the **Y-PAR-masked + rCRS** analysis-set FASTA (`chm13v2.0_maskedY_rCRS.fa.gz`)
  for short-read calling. PAR-masking avoids the X/Y multi-mapping that produced
  some of our artifacts; the rCRS mito matches the mtDNA/haplogroup work.
- **Reference-polarity metadata.** Record which assembly's Y/mt sits in a given
  reference (CHM13v2.0 Y = HG002 = J). Anc/der must come from the tree; the
  reference is only the coordinate system. Add a guard/test so this never silently
  flips.
- **Callable-region masks.** Carry BEDs for the Y palindromes/ampliconic/
  heterochromatic regions and the `unique_to_*` non-syntenic regions already
  catalogued in the CHM13 resources doc. Down-weight or exclude calls there. The
  Y profile already annotates regions; this extends it to a hard callability mask.

### `navigator-sync` / `du-domain` (Atmosphere Lexicon / PDS)

- **Callset provenance + completeness.** A result derived from a genotype-only
  graph callset is panel-conditioned and may have silent omissions. The schema
  should record the analysis source (linear / surjected-from-pangenome /
  graph-genotype-only) so consumers and the reconciliation engine treat
  completeness correctly.
- **Reconciliation guard.** The existing multi-run reconciliation must never let a
  graph callset's "absent" be scored as "ancestral/negative." Extend the
  source/quality tracking (already present) with a completeness dimension, so a
  surjected-linear or BigY result outranks a graph-genotype-only one for
  terminal/private calls.

---

## Decoding Us API alignment

- The **Y tree** used in the research is the same source as the DecodingUs tree
  provider here. Its `hs1`/GRCh38/GRCh37 coordinates plus ancestral/derived
  alleles are exactly what the region/marker work needs.
- The **mt-tree** endpoint is **not yet populated**. mtDNA-from-pangenome is parked
  for two reasons: no mt-tree data, and CHM13v2.0's chrM is not rCRS (the Navigator
  already handles the rotation-aware rCRS↔chrM map, so the right move is to keep
  mtDNA on rCRS via that path, not the pangenome).

---

## Open questions / decisions

1. Embed `vg` as an optional local mapping helper (option B), or keep graph
   mapping strictly upstream (option A)? Default to A.
2. Do we ever ingest a GAM/graph VCF directly, or only surjected linear CRAM with
   provenance? Leaning: only surjected linear, at least initially.
3. Is augment-based discovery ever an edge feature, or strictly server-side? It is
   heavy (a full-GAM pass per sample in the research — the GAM was ~110 GB for one
   sample; see "Empirical resource footprint"); not edge, and a non-trivial server
   cost.
4. Where do the callable-region masks live and how are they versioned alongside the
   reference cache?
5. Atmosphere Lexicon: how to represent callset source/completeness so the
   federation reconciles graph-derived vs linear results correctly.

---

## Proposed phasing

- **Phase 1 (no pangenome needed): harden placement.** Port the parsimony
  path-descent + marker-less-node collapse + haploid allele-balance filter into
  `du-bio`/`navigator-analysis`. Immediate accuracy/robustness win on existing
  linear inputs. Validate against the known panel + the research truth set.
- **Phase 2: reference + region correctness.** Adopt the masked+rCRS analysis-set
  FASTA as canonical, add reference-polarity metadata + a polarity guard, and wire
  callable-region masks into the Y profile/quality.
- **Phase 3: pangenome-surjected inputs.** Accept surjected linear CRAMs, label
  their provenance, and add the reconciliation completeness guard so graph-derived
  results never silently win or lose against linear ones.
- **Phase 4 (decision-gated): local graph mapping.** Optional bundled `vg`
  pre-step, only if there is real user demand and we accept the purity compromise.

---

## Backlog seeds (for `documents/BACKLOG.md`)

- 🟢 Harden Y/mt placement with path-supported parsimony (🟢 **done + validated** —
  `haplo::path_admissible` guard over the Kulczynski ranking, wired at
  `assemble_assignment`; GFX0457637/CHM13 → `R-FGC29071` + `U5a1b1g`) + marker-less-node
  collapse (🟡 — handled inline: no-call/marker-less nodes pass the guard; standalone
  recursive collapse still net-new) + haploid allele-balance filter (🟢 **done +
  validated** — `caller::is_paralogous` at all three haploid call sites; the paralog
  *false-derived* defence the guard intentionally left open). *Pure win, no pangenome.*
- 🔴 Adopt `chm13v2.0_maskedY_rCRS` as the canonical short-read analysis reference
  (new `Build` variant + download entry in `navigator-refgenome/registry.rs`).
- 🔴 Reference-polarity metadata + guard (CHM13 Y = HG002 = J; anc/der from tree only).
- 🟡 Callable-region/paralog masks for distal Y: 🟢 `mask.rs` loader exists; 🔴 wire
  it into the Y profile / region-aware quality (no consumer today).
- 🔴 Callset source/completeness in the Atmosphere Lexicon + reconciliation guard for
  graph-derived (panel-conditioned) results (new field on `RunHaplogroupCallRecord`).
