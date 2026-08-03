# Project Y block tree — design

**Status:** **All three phases implemented** (2026-08-02, branch `feat/project-block-tree`) — the
aggregate + collapse (1), the `ProjectTab::Tree` canvas view (2), and private-variant blocks +
shared-private candidate branches + TSV/HTML export (3). 30 unit tests; validated live on two
1,900-member cohorts. **Candidate detection is inert until private-Y is computed for more than one
subject** — see §9 phase 3 and §11 Q4. Drafted 2026-08-02.
**Closes:** `BACKLOG.md` §3.3 remaining scope ("a genuinely zoomable / searchable *whole-tree* view
with the subject's placement highlighted") — reframed as **cohort-scoped**, which is the form that
actually earns its keep.
**Scope:** `navigator-analysis` (the pure induced-subtree logic), `navigator-app` (the aggregate,
builder and collapse), `navigator-ui` (a new `ProjectTab::Tree` + `ui/blocktree.rs`, phase 2). No new
analysis, no new asset, no migration. *(The draft put collapse in `navigator-domain`; the layering
rules that out — see §7.)*

## 1. Problem

Navigator places each subject on the Y tree **individually**: `ui/descent.rs` draws one subject's
root→terminal path, `ui/branch.rs` gives a per-marker report for one node's subtree. Both are
per-subject by construction.

The question a *group project* asks is cross-subject: where do these members sit relative to one
another, which branches do they share, who is closest to whom, and which members share unnamed
(private) variants and therefore constitute a branch that does not exist in the tree yet. FTDNA's
Block Tree is the canonical surface for exactly this, and is the most-used Big Y feature.

Navigator has no cohort-scoped tree view at all — `ProjectTab` is `Members / Report / Ystr`
(`navigator-ui/src/ui/mod.rs:142`).

This became worth building now because the workspace acquired a real cohort: the multi-lab D2C
collection (FTDNA / FGC / YSEQ / Dante / Nebula), which is Y-heavy — Big Y CRAMs, haploid GVCFs and
STR panels — loaded by `scripts/import-d2c/`.

## 2. What a "block" is

- A **block** is a tree node together with the run of SNPs that are phylogenetically **equivalent**
  on that branch: every member below it carries all of them, and nothing observed separates them.
  `HaploNode.loci` already *is* that set.
- The **view** is the **induced subtree** of the full haplotree spanning the project members'
  terminal haplogroups — root→terminal for every member, unioned — with each member hanging off its
  terminal block.
- Below each member sit their **private variants** as an unnamed terminal block. This is the part
  that makes the view worth more than a table: two members sharing a private variant are a candidate
  *new named branch*. Navigator computes private variants itself
  (`App::donor_private_y` → `PrivateBucket`), so we can surface the shared-private case directly
  rather than waiting for a vendor to name the branch.

## 3. Substrate — all of it already exists

| Need | Existing |
|---|---|
| Topology + defining SNPs | `HaploTree` / `HaploNode` / `Locus` (`navigator-analysis/src/haplo.rs:33`) |
| Tree fetch + parse | `App::fetch_decodingus_y_tree` / `fetch_ftdna_y_tree`; `parse_decodingus_json(json, build_key)` / `parse_ftdna_json` |
| Project membership | `biosample::list_members_for_project` |
| Every member's terminal | `App::haplogroup_terminals()` — **one bulk workspace query**, already used by `project_report` and `project_str_overview` |
| Private variants | `App::donor_private_y(guid) -> Option<PrivateBucket>` (reads cached `private_y`; no recompute) |
| Per-subject precedent | `App::descent_report` (`haplogroup.rs:1786`): provider selection, terminal-name → node-id, ancestor walk |
| Project-wide virtualized chart precedent | `project_str_chart` and the virtualized project Y-STR chart |

This is a **new view over built machinery**, not new science. That is the reason to pick it over the
other open backlog threads.

## 4. The aggregate

Built by `App::project_block_tree(project_id: i64, dna: DnaType) -> Result<Option<ProjectBlockTree>, AppError>`.

```rust
pub struct ProjectBlockTree {
    pub dna: DnaType,
    /// Induced-subtree nodes in a stable pre-order (parent always precedes its children).
    pub blocks: Vec<Block>,
    /// Members with no terminal, or whose terminal is absent from this tree.
    pub unplaced: Vec<UnplacedMember>,
    /// The tree this view was drawn on ("decodingus" | "ftdna") — the loci below belong to it.
    pub provider: String,
    /// The coordinate space `Block::loci` positions are in (see §6).
    pub build_key: String,
}

pub struct Block {
    pub node_id: i64,
    pub name: String,
    pub parent: Option<i64>,
    pub depth: usize,
    /// The equivalent SNPs defining this block (the node's own loci).
    pub loci: Vec<Locus>,
    /// Members whose terminal *is* this node.
    pub members: Vec<BlockMember>,
    /// Members at or below this block — the count the UI badges on a collapsed branch.
    pub subtree_members: usize,
}

pub struct BlockMember {
    pub guid: SampleGuid,
    pub name: String,
    /// Unnamed variants below the terminal. `None` when private-Y has never been computed
    /// for this subject — distinct from `Some(0)`, which means "computed, none found".
    pub private_novel: Option<usize>,
    pub private_total: Option<usize>,
}
```

## 5. Building it

1. Fetch + parse the tree **once for the whole project** (the same provider branch `descent_report`
   uses). This is the expensive step — a multi-MB document — and the per-subject path pays it per
   subject today.
2. Invert `children` into a parent map. `HaploTree` stores children only.
3. Build a **name → node-id index** in one pass. `descent_report` does a linear
   `tree.nodes.iter().find(...)` per call, which is fine for one subject and quadratic for N members.
4. For each member with a Y terminal: resolve name → id, walk to root collecting the path, union all
   paths into a node-id set.
5. Emit `blocks` as that induced subtree in pre-order, each node carrying its `loci`, its own
   members, and a `subtree_members` roll-up.
6. A member whose terminal is **absent from the tree** (provider/build skew — `descent_report`
   already returns `None` for this case per subject, `haplogroup.rs:1822`) goes to `unplaced`. It is
   never silently dropped: on a multi-lab cohort, skew is expected and hiding it would misrepresent
   the project's size.

## 6. The build-key question — decide before coding

`parse_decodingus_json` takes a `build_key`, and `descent_report` chooses it **per subject** via
`subject_y_build_key`. A cohort spans subjects on different builds, so there is no single
per-subject answer.

Resolution: **topology and node names are build-independent**; only the displayed `loci` *positions*
are not. So parse once on a single project-level build key (the modal build across members, falling
back to the app default), record it as `ProjectBlockTree.build_key`, and label the view with it.

This cannot introduce a cross-build **placement** error, because this view consumes terminals that
were already placed elsewhere and never re-places anything (§8). The only build-sensitive thing on
screen is a coordinate label, and it is labelled.

Implementation found two cases the "modal build" rule does *not* cover, so `provider` and `build_key`
are now taken from whatever the fetch actually **resolved to**, not from what was configured:

- The **FTDNA Y tree** is published on GRCh38 regardless of what the members are aligned to, so the
  modal build key would have mislabelled it.
- The **mtDNA path falls back to FTDNA at runtime** (`mt_tree_rcrs` returns its own provider tag)
  when the DecodingUs tree can't be remapped for want of a cached CHM13 reference — and either way
  its loci are **rCRS**, not the Y coordinate space.

## 7. Collapsing

An induced subtree over a deep haplotree is mostly single-child chains carrying no members. Options:

- **(a) Show every node** — faithful, but yields long runs of empty blocks per member.
- **(b) Collapse runs of member-less single-child nodes into one "N intermediate branches" block,
  expandable.** This is what FTDNA does and what makes the view readable.

**Chose (b)**, as a **pure function over `blocks`** so it is unit-testable without a tree fetch.

Three things the implementation settled:

- **It lives in `navigator-app`, not `navigator-domain`.** `Block` carries `haplo::Locus`, and
  `navigator-domain` sits *below* `navigator-analysis` in the layering — it cannot name that type.
  The genuinely tree-shaped half went to `navigator-analysis::haplo` instead (`induced_subtree`,
  `name_index`, `InducedNode`), where it is pure and identity-free; `collapse_blocks` and
  `roll_up_subtree_members` stayed in `navigator-app::blocktree` next to the aggregate they shape.
- **Merging is semantic, not cosmetic.** The absorbed loci join the survivor's own, root-most first,
  because within *this cohort* the run genuinely is one undivided block. The absorbed names are kept
  in `Block::collapsed` so the view can still say what it folded.
- **The tree root gets absorbed too**, when it is member-less with a single child — which it usually
  is. That is wanted: the haplotree root is Y-Adam, and a cohort of R-men should open on `R`, not on
  a chain of empty backbone nodes. The topmost survivor simply becomes the new root
  (`collapse_stops_a_run_at_a_placed_branch` pins this).

`subtree_members` needs no recomputation after a collapse: an absorbed node has no members of its own
and exactly one child, so its count already equals the survivor's.

## 8. UI

A new `ProjectTab::Tree` (`project.tab.tree`, en + es), rendered from a new `ui/blocktree.rs` —
following the one-view-per-module split, and keeping `central.rs` from growing again.

- Layout is precomputed **once per `(project, tree)`** into node rects, not rebuilt per frame; draw
  culled to `ui.clip_rect()`.
- Depth → x, leaf order → y (classic block-tree columns). Scroll + zoom. Click a block to expand its
  equivalent SNPs; click a member to open that subject — the `return_to_project` round-trip already
  exists (`central.rs:713`).
- **No `ComboBox` anywhere** in this view, per the roster-picker rule: a project can hold thousands
  of members, and a `ComboBox` builds a widget per entry per frame.
- Loading goes over the worker thread (`Command::LoadProjectBlockTree` / `Event::ProjectBlockTreeReady`)
  because it fetches and parses a multi-MB tree.

## 9. Phasing

- **Phase 1** — ✅ **done.** The aggregate, `App::project_block_tree`, the collapse function, and 14
  unit tests over synthetic trees (no network, no live DB). No UI. Landed as
  `navigator-analysis/src/haplo.rs` (`induced_subtree` / `name_index` / `InducedNode`, 6 tests) +
  `navigator-app/src/blocktree.rs` (builder, `collapse_blocks`, `roll_up_subtree_members`, 8 tests).
- **Phase 2** — ✅ **done.** `ProjectTab::Tree` + `navigator-ui/src/ui/blocktree.rs`: a pure `layout`
  (tidy-tree — depth→x, parent centred on its children's extent) with 6 tests, clip-culled canvas
  painting, expand-on-click, SNP-name hover, zoom, double-click-to-open-subject, en/es keys.
  `Command::LoadProjectBlockTree` / `Event::ProjectBlockTree`, loaded **lazily on first view of the
  tab** rather than on project select — unlike the STR chart, this fetches a multi-MB haplotree.

  **Validated on live data** (`cargo run -p navigator-app --example blocktree_check -- <project_id>`,
  kept as the harness for phase 3):

  | project | members | blocks | placed | no Y placement | terminal absent from tree | build |
  |---|---|---|---|---|---|---|
  | `R1b-CTS4466Plus` | 1881 | 212 | 243 | 1538 | 100 (72 distinct names) | 1.5 s |
  | `R-L21_South_Irish` | 1906 | 348 | 269 | 1637 combined | — | in-GUI |

  Root absorption behaves as intended: the whole Y-Adam→DF13 backbone folds into a single
  `R-DF13 [1769 SNPs] (+27 folded)` block, and the tree opens on `R-CTS4466` where the cohort
  actually lives. Max depth after collapse is 15; 51 branches folded overall.

  The 100 "terminal absent from tree" cases are **genuine provider skew, not a lookup bug** —
  `F-M89` and `R-A1133` really are absent from the (fresh, 57 MB) DecodingUs tree while `R-A9426` is
  present; those terminals came from FTDNA-provenance calls. This is exactly what `unplaced` exists
  to surface.
- **Phase 3** — ✅ **built.** Per-member private counts (`BlockMember.private_novel`/`private_total`,
  from a new bulk `App::private_y_for_biosamples`), **shared-private detection** inserting candidate
  branches, and `export::block_tree_tsv` / `block_tree_html`. 10 further tests.

  **Candidate branches.** Within a block, variants carried by exactly the same set of members are
  equivalent — the same reasoning that makes a named node's SNPs a block — so each distinct carrier
  set becomes one candidate block, `candidate: true` with a synthetic negative `node_id` and an empty
  `name` (the view localizes the label; the exports write `candidate`). They are amber, never green:
  an inference, not a published branch.

  Three rules keep it from manufacturing branches:
  - **Only novel, unique-sequence calls count.** Off-path-*known* variants support an existing finer
    branch (a placement question). Structural-region calls sit in chrY palindromes and amplicons,
    where two men "sharing" a call are far likelier to share a mapping artefact than an ancestor.
  - **Groups must stay laminar** — any two accepted carrier sets are disjoint or nested, so the
    result is a tree. Nested sets nest as branches; a partial overlap is real conflict (recurrence or
    genuine disagreement) and is **counted** in `candidate_conflicts`, not forced into a shape it
    doesn't fit.
  - **A member with no computed private-Y is not grouped**, and shows `None`, not `0`.

  > **Inert on this workspace, and that is the finding.** Exactly **one** subject in the whole
  > database has private-Y computed, so no block can contain two carriers and zero candidates appear.
  > The logic is unit-tested against all of the above cases, but it cannot demonstrate on live data
  > until private-Y exists for more members. This settles §11 Q4: the batch action is **required**,
  > not optional — see the open question below for where it should live.

  **Performance regression found and fixed while validating.** Adding the private-Y load took the
  R1b-CTS4466Plus build from 1.5 s to 25 s. The cause was not the artifact freshness stats (the first
  guess, and wrong) but `artifact::list_for_alignments` selecting `payload` for *every* artifact kind
  — and `tree-genotype` alone is **2.9 GB across 680 rows**. It was reading gigabytes of JSON to find
  one `private_y` row. Fixed with a new `artifact::list_for_alignments_of_kind`; back to **1.7 s**.
  Note `project_report` still uses the unfiltered query and will have the same cost on such a project.

## 10. Out of scope

- **mtDNA block tree.** The aggregate is `DnaType`-generic so allowing it costs nothing, but the
  cohort surface that matters is Y; validate there first.
- **TMRCA / branch ages.** Needs a mutation-rate model we do not have. `ymatch::Tmrca` exists but is
  pairwise, not tree-wide.
- **Re-placement.** This view *reads* placements. It never re-places, and must not.

## 11. Open questions

1. **Collapse threshold** — collapse a run of ≥2 member-less nodes, or ≥5? Affects readability only;
   pick one and make it a constant.
2. **Scope** — project-only, or also a "whole workspace" mode? The D2C cohort is a single project, so
   project-scoped is sufficient for v1.
3. **Terminals query cost** — `haplogroup_terminals()` is workspace-wide. On a 10k-subject workspace
   that is one large query per view open; cache the result on the aggregate rather than re-querying
   on redraw.
4. **Private-Y coverage** — **answered by phase-3 validation: a batch action is required.** Exactly
   one subject in the workspace has private-Y computed, so shared-private detection can never fire.
   The open part is *where it lives*: a bespoke "compute private-Y for this project" button, or
   folded into the existing project-wide analyze / deep-analyze streaming flow — which is the same
   choice `BACKLOG.md` §1.2 already faces for panel genotyping, and probably wants the same answer.

   **The deeper cause, and the prerequisite now built.** Private-Y is keyed on `alignment_id` and
   sourced from a pileup walk or a GVCF sidecar found beside the BAM — so it is offered for
   BAM/CRAM only. Most of this workspace's Y data is *externally processed VCFs*, which live in
   `variant_set`/`variant_call` keyed on `biosample_guid` with no alignment at all: **7,842 subjects
   have chrY calls, against the 1 with private-Y**. An engine over those call sets is the real fix.

   It could not be written usefully first, because the import threw the evidence away: `variant_call`
   stored only contig/position/ref/alt/rsID/genotype, so nothing downstream could tell a 40× hom-alt
   call from a 2-read artefact. **Migration 0042 + `variants::CallEvidence` now capture QUAL, FILTER,
   DP, GQ and AD** (`ad_alt` follows the genotype-selected ALT, so multi-allelic rows read correctly),
   and `variant_set.call_schema` records whether a set has evidence — derived from what was captured,
   never from the importer version, so it cannot promise evidence a sites-only VCF never had. Absent
   fields stay `None`; reading them as `0` would make a good call look unsupported.

   Real FTDNA Big Y (aengine) files confirm the value: ~218k `PASS` against ~44k FILTER-flagged calls
   in a single sample, plus per-call AD/DP/GQ — exactly the gate a VCF-backed private-Y engine needs.
   Existing sets keep `call_schema = 1` and must be re-imported to gain evidence.

   **Routing bug — found, root-caused, fixed.** A chrY-only Big Y VCF was auto-detected as an
   "Autosomal 1240K call set" and landed in `external_panel_dosage` (266 chrY panel loci matched out
   of ~260k records), creating no Y `variant_set` — so no Y placement and no private-Y source. The
   7,834 existing `FTDNA Big Y (aengine)` sets predate that routing, which is why nothing looked
   broken: this detector would have prevented creating them.

   *Why it existed.* `looks_like_genotyped_callset_vcf` arrived with the external-caller-precedence
   work (`5636286`), whose subject was **autosomal** call sets. It classified any VCF emitting
   explicit `0/0` rows as a call set — a sound "call set vs variant list" test for the autosomal
   problem it was solving — and the commit stated its own assumption plainly: *"chrY/chrM GVCFs are
   the sidecar fast path, discovered in a directory, not here."* That held for the directory flow it
   was designed with, but not for a loose `variants.vcf.gz` handed to `ingest`. Vendor Y products
   report reference sites too, so they tripped the one signal it looked at, and nothing looked at
   **contigs**.

   *The fix.* `vcf_known_lineage_only` gates both the `.g.vcf` and genotyped-`.vcf` branches:
   `##contig=<ID=…>` declarations are authoritative when present, else the records' `CHROM` column.
   No contig evidence at all returns `false` — absence of evidence is not evidence of absence, so an
   unreadable `.g.vcf` keeps the claim its extension makes rather than being demoted on a guess.
   5 tests, including the real aengine header and a chrM-only equivalent.

   Verified end to end on the same file that exposed it: now imports as `FTDNA Big Y (aengine)` /
   `TARGETED_NGS` / GRCh38, **4,485 chrY calls** (matching the ~4,184 average of the sets created
   under the old routing) at `call_schema = 2`, with 1,590 FILTER-flagged calls and DP spanning
   1–4,212. The two fixes compose: a `dp=37, qual=1484` call and a `dp=1, gq=1` one-read artefact are
   now distinguishable, which is precisely what the private-Y gate needs.
5. **Suffixed terminal names** (found in phase-2 validation) — some unresolved terminals are a real
   node name plus a suffix, e.g. `R-A9426:n0` where `R-A9426` *is* in the tree. Stripping the suffix
   and matching the parent node would recover those members, but only if the suffix means what it
   looks like; worth confirming against whatever writes it before special-casing anything.
6. **Whole-workspace mode** — resolved for v1: **project-scoped only**.
