# Project Y block tree — design

**Status:** **Phase 1 implemented** (2026-08-02, branch `feat/project-block-tree`) — the aggregate,
`App::project_block_tree`, the collapse function and 14 unit tests. No UI yet (phase 2). Drafted
2026-08-02.
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
- **Phase 2** — the `ProjectTab::Tree` view: layout, culling, expand/collapse, click-through, i18n
  keys at en/es parity.
- **Phase 3** — private-variant blocks, including **shared-private detection** (two or more members
  carrying the same unnamed variant = a candidate new branch), plus TSV/HTML export alongside the
  existing branch/descent exports in `export.rs`.

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
4. **Private-Y coverage** — `donor_private_y` reads *cached* results, so members who have never had
   private-Y computed show `None`, not `0`. Should Phase 3 offer a batch "compute private-Y for this
   project" action, or leave it to the per-subject path?
