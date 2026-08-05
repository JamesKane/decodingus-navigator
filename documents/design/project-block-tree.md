# Project Y block tree — design

**Status:** **All three phases shipped** (branch `feat/project-block-tree`, 13 commits, not pushed) —
the aggregate + collapse (1), the `ProjectTab::Tree` canvas (2), and private-variant blocks +
candidate branches + export (3), plus four things phase 3 turned out to need: the
`private-y --project` batch, a **VCF-backed private-Y engine** for subjects with no alignment, a
**candidate review surface**, and an artefact-filter stack calibrated against R1b-CTS4466Plus.
Live state there: 248/255 placed members carry private-Y (was 1 workspace-wide), **7 candidate
branches** surviving the filters. Suite 797 passed. Open items in §11. Drafted 2026-08-02.

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
- **Phase 3** — ✅ **built**, then calibrated over several rounds against real data (§11 records what
  is still open). The engine, the batch that feeds it, the review surface, and the filter stack are
  described below.

- **Phase 3 (original scope)** — ✅ **built.** Per-member private counts (`BlockMember.private_novel`/`private_total`,
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

  **Two further filters, added after the first live run.** Batch-computing private-Y for
  R1b-CTS4466Plus (`navigator private-y --project`, 228 alignments, 2,157 novel unique-sequence
  variants) produced 8 candidate branches — and nearly all of them were artefacts:

  - **Proximity.** Every candidate variant sat in a tight cluster: six positions inside 32 bp
    (`16342231…16342263`, gaps of 5–8 bp) and three inside 23 bp. Real Y mutations are megabases
    apart; a handful of novel calls within tens of bases is one misaligned read smearing several
    false SNVs — the same reasoning behind the GVCF path's depth floor. `drop_clustered` removes the
    **whole** cluster (`CANDIDATE_MIN_SEPARATION_BP = 100`), because when several calls share one
    mapping event there is no basis for electing one of them the real mutation.
  - **Cross-block recurrence.** `11311865` and `11311870` each defined a candidate under *two*
    different parent blocks. A variant that arose more than once cannot mark a new branch, and the
    laminar check cannot see this — it reasons inside a single block. `recurrent_positions` runs
    across all blocks before any group is accepted, and the count is surfaced as
    `candidate_recurrent`.

  Result on the same cohort: **8 candidates → 3**, all conflicts gone, and no candidate left showing
  a single member. The survivors each rest on one isolated position with 2–4 carriers. Note the
  recurrent count then reads 0 — not because recurrence is absent, but because proximity caught those
  same positions first; the guard still covers an unclustered recurrent call.

  > **Was inert until the batch existed.** Before it, exactly **one** subject in the database had
  > private-Y, so no block could hold two carriers. This settled §11 Q4: the batch action was
  > **required**, not optional. It now exists as `private-y --project` (resumable — cached buckets
  > are skipped; `--force` recomputes), and the CTS4466 run took private-Y from 1 subject to 229.

  **Performance regression found and fixed while validating.** Adding the private-Y load took the
  R1b-CTS4466Plus build from 1.5 s to 25 s. The cause was not the artifact freshness stats (the first
  guess, and wrong) but `artifact::list_for_alignments` selecting `payload` for *every* artifact kind
  — and `tree-genotype` alone is **2.9 GB across 680 rows**. It was reading gigabytes of JSON to find
  one `private_y` row. Fixed with a new `artifact::list_for_alignments_of_kind`; back to **1.7 s**.
  Note `project_report` still uses the unfiltered query and will have the same cost on such a project.

## 9b. What phase 3 grew into

The original phase-3 scope assumed private-Y already existed for a cohort. It did not — one subject
in the whole workspace had it — so four further pieces were needed before candidate branches could
mean anything.

**`private-y --project`** — computes and persists the bucket for every alignment in a project.
Resumable (cached buckets skipped, `--force` recomputes) and it skips alignments whose file is gone
rather than reporting a missing vendor download as a failure. Took CTS4466 from 1 subject to 229.

**A VCF-backed engine** (`App::private_y_from_variant_set`, migration 0044). Private-Y was keyed on
`alignment_id`; ~1,600 of CTS4466's 1,881 members have no alignment at all. The engine classifies a
call set the same way the alignment path classifies a walk, placing from `vset_base_calls` so the
terminal rests on tree-position genotypes. There is no self-callable mask — a VCF has no coverage
track — so the source's own evidence takes its place.

**A review surface.** A candidate is *our inference*, not the tree's assertion, so it carries its
evidence: clicking one shows every carrier's DP, AD(alt), derived fraction and publish-gate verdict,
and the same rows go to the TSV export. This is what found the paralog the automated filters all
missed — two carriers at DP 413/504 against a median of 57. Build the inspection surface *before*
the filters next time; the filters written in the abstract caught clustering and recurrence, but the
one that mattered came from looking at the data.

**The filter stack**, each added because real data demanded it, not by anticipation:

| gate | rejects |
|---|---|
| novel + unique sequence | off-path-known SNPs; palindrome/amplicon calls |
| per-donor proximity (100 bp) | one donor's calls smeared across a misaligned read |
| cross-block recurrence | a position defining branches under two parents — it arose twice |
| cohort frequency (>25%, abstains <20 donors) | population/reference differences, not private variants |
| cross-candidate clustering (1 kb) | separate "branches" inside one sequencing fragment |
| *VCF sources only:* PASS · DP≥4 · GQ≥20 · AF≥0.95 · hemizygous · ≤3× donor median depth | non-deterministic calls; chrY heterozygotes; collapsed-repeat pile-ups |

Two properties worth remembering. **The filters interact** — dropping the 56.83 Mb cluster changed
which groupings passed the laminar check and surfaced a new candidate, so the set is never simply the
previous one minus removals. And **each trades recall for precision**, which is a judgement about a
cohort rather than a fact about the code: 8 candidates became 3, then 20 → 12 → 9 → 7 as the VCF
engine widened the input and the gates tightened around it.

## 10. Out of scope

- **mtDNA block tree.** The aggregate is `DnaType`-generic so allowing it costs nothing, but the
  cohort surface that matters is Y; validate there first.
- **TMRCA / branch ages.** Needs a mutation-rate model we do not have. `ymatch::Tmrca` exists but is
  pairwise, not tree-wide.
- **Re-placement.** This view *reads* placements. It never re-places, and must not.

## 11. Open questions

### Resolved

1. **Collapse threshold** — a run of ≥2 member-less nodes (`COLLAPSE_MIN_RUN = 2`).
2. **Scope / whole-workspace mode** — **project-scoped only** for v1.
3. **Terminals query cost** — moot in practice. The Tree tab loads lazily and caches the aggregate in
   `project_blocktree`, so `haplogroup_terminals()` runs once per project open rather than per
   redraw; builds measure 0.7–1.7 s on a 9,853-subject project.

### Open

4. **Where the private-Y batch lives — half answered.** The batch itself exists as
   `navigator private-y --project` (resumable; `--force` recomputes; skips alignments whose file is
   gone). **There is no GUI trigger**, so a user who never touches the CLI cannot populate private-Y —
   and without it candidate branches cannot fire at all. The original question stands: a bespoke
   button, or folded into the project-wide analyze / deep-analyze streaming flow. `BACKLOG.md` §1.2
   faces the same choice for panel genotyping and probably wants the same answer.

5. **Suffixed terminal names — untouched.** **162** Y consensus labels carry a `:`-suffix
   (`R-A9426:n0`, where `R-A9426` *is* in the tree), so those members land in `unplaced` instead of on
   a branch. Stripping the suffix and matching the parent would recover them — but only if the suffix
   means what it looks like, and nothing has yet established what writes it. Confirm before
   special-casing.

6. **A regional concentration in the surviving candidates.** Three of seven span **88 kb** at
   10.79–10.88 Mb. That is far beyond one sequencing fragment, so the 1 kb cross-candidate rule
   correctly leaves them; widening it until they vanish would fit the filter to the observation rather
   than to a mechanism. Left visible rather than filtered — a judgement about this cohort.

7. **Per-donor novel counts sit above the alignment path's** (median ~73 vs 3–13). This looks like
   instrument difference — GATK HaplotypeCaller at ploidy 1 against a vendor diploid caller — but it
   is not proven. A subject carrying *both* a CRAM and a vendor VCF would settle it directly, by
   classifying the same donor twice.

### Debts this work incurred

- **`project_report` still uses the unfiltered `artifact::list_for_alignments`**, which selects
  `payload` for every artifact kind — the query that read gigabytes of `tree-genotype` JSON and took
  the block-tree build from 1.5 s to 25 s. `artifact::list_for_alignments_of_kind` exists now; that
  caller was never converted.
- **The canvas has no interaction test coverage.** The click bug that made candidate review
  unreachable shipped in phase 2 and the layout tests could not have caught it: they exercise the
  pure `layout()` function, which has no input handling. Any test of click routing would need to
  drive `egui` directly.
