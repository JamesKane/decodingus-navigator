# Project Y block tree — design

**Status:** **All three phases shipped** (branch `feat/project-block-tree`, 15 commits, not pushed) —
the aggregate + collapse (1), the `ProjectTab::Tree` canvas (2), and private-variant blocks +
candidate branches + export (3), plus four things phase 3 turned out to need: the
`private-y --project` batch, a **VCF-backed private-Y engine** for subjects with no alignment, a
**candidate review surface**, and an artefact-filter stack calibrated against R1b-CTS4466Plus.
Live state there: 248/255 placed members carry private-Y (was 1 workspace-wide), **7 candidate
branches** surviving the filters. Suite 797 passed. The canvas was then **redrawn to Alex
Williamson's Big Tree presentation** (§8) — the thing the name "block tree" refers to. Open items in
§11. Drafted 2026-08-02, rendering revised 2026-08-05.

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

## 8. UI — the Big Tree presentation

A new `ProjectTab::Tree` (`project.tab.tree`, en + es), rendered from a new `ui/blocktree.rs` —
following the one-view-per-module split, and keeping `central.rs` from growing again.

The draft specified "depth → x, leaf order → y … click a block to expand its equivalent SNPs". That
was wrong on three counts, and the corrections all come from one principle:

> **Mutations accrue at a roughly steady rate, so a count of SNPs is a measure of elapsed time.**

Everything below follows from taking that seriously, which is what Williamson's Big Tree does and
what FTDNA's Block Tree borrowed.

- **Top-down, not left-right.** Depth → y. Time runs down the page.
- **A block's height is its SNP count, and nothing is elided.** Every equivalent SNP gets its own
  line, so the box's height *is* how long that branch ran unbroken. There is no cap and no
  expand-on-click: a truncated box is a shortened box, and a shortened box misreports the time. (A
  line may carry several *names* — `BY30547 Y43043` is one mutation with two — but never two
  mutations.) This removes the draft's expand interaction entirely, along with its state.
- **Vertical position is cumulative.** A block hangs directly beneath its parent rather than on a row
  shared with everything at its depth. How far down a block sits is therefore the mutations
  accumulated along the path to it: height reads as time *within* a branch, y as time *along* the
  lineage. Depth-aligned rows would pad every short branch out to the tallest box beside it — a lot
  of empty canvas, and a lie about when the branch happened, since two lineages at one depth would
  draw level even when one had accrued thirty more mutations than the other.
- **Men are blocks of their own.** Each biosample is a grey box on a stem below the block it is
  placed on, taking a horizontal slot beside that block's child subtrees so a lineage that both
  splits and holds men makes room for both. Click one to open its subject — the `return_to_project`
  round-trip already exists (`central.rs:713`). A roster side panel stays alongside, because it is
  the only place the private-variant counts are scannable as a list.
- **The backbone above the cohort is a breadcrumb, not a block.** See §8b.
- ASCII markers only in the canvas: `▾` is absent from the bundled font and renders as tofu.
- Layout is precomputed **once per `(tree, zoom)`** into node rects, not rebuilt per frame; draw
  culled to `ui.clip_rect()`. `layout()` is a pure function over `&[Block]`, so it is testable
  without a canvas.
- **No `ComboBox` anywhere** in this view, per the roster-picker rule: a project can hold thousands
  of members, and a `ComboBox` builds a widget per entry per frame.
- Loading goes over the worker thread (`Command::LoadProjectBlockTree` / `Event::ProjectBlockTreeReady`)
  because it fetches and parses a multi-MB tree.

### 8b. Why the backbone can't be a block

Uncapping the SNP list is only affordable because of one exception. R1b-CTS4466Plus induces a root of
`R-Z290`: **24 folded branches, 1,763 SNPs** — twelve times taller than the entire cohort hanging off
it. The next largest block in that project is **24 SNPs**, so there was never a general problem, only
this one.

A collapsed run is by construction *more than one branch*. Its height is a sum across the tree
**above** the cohort, not one branch's elapsed time — the single place the height-as-time reading
does not hold. So it leaves the canvas for a path across the top
(`Y › A0-T › … › R-P312 › R-Z290  (1763 SNPs upstream)`), exactly as the Big Tree's subclade pages
render everything above the clade in view.

The test is **the fold, not member-lessness**. Two shallow kits sit on `R-Z290`; a member-less guard
would never have fired where it mattered. They keep their roster because the breadcrumb is
selectable.

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

5. **Suffixed terminal names — ✅ diagnosed; the proposed fix was wrong.** The draft suggested
   stripping the `:` suffix and matching the parent. **Do not.** Those 162 labels are not a naming
   convention — they are *stale placements against a superseded tree generation*:

   - all are `:n0`/`:n1`, plus three indel-defined (`:6686542 AGT->A`);
   - the current cached DecodingUs tree contains **zero** such node names;
   - nothing in Navigator writes the suffix — it came from the tree;
   - **153 of the 159 fingerprinted ones sit on one tree fingerprint**, `yt:b211464f1bd97ca6`.

   Stripping the suffix would have silently demoted 153 subjects to a *less derived* branch in order
   to work around a stale placement. Closed as won't-fix, subsumed by §12.

6. **A regional concentration in the surviving candidates.** Three of seven span **88 kb** at
   10.79–10.88 Mb. That is far beyond one sequencing fragment, so the 1 kb cross-candidate rule
   correctly leaves them; widening it until they vanish would fit the filter to the observation rather
   than to a mechanism. Left visible rather than filtered — a judgement about this cohort.

7. **Per-donor novel counts sat above the alignment path's** (median ~73 vs 3–13) — ✅ **cause
   found, and it was not the instrument.** The structural-region masks are CHM13-native, and both
   private-Y paths applied them *only* to a CHM13 source. On a GRCh38 set `regions` was `None`, so
   every palindromic and amplicon call counted as unique sequence. The alignment path masked those
   regions; the VCF path, being mostly GRCh38, did not. Hence one branch in R1b-CTS4466Plus
   averaging **661** private variants beside another averaging 4.

   Fixed by lifting the masks — see §13. The GATK-vs-vendor-caller hypothesis is withdrawn; it was
   a plausible story that happened to be wrong, and the give-away (only GRCh38 sources affected)
   was in the data all along.

### Debts this work incurred

- **`project_report` still uses the unfiltered `artifact::list_for_alignments`**, which selects
  `payload` for every artifact kind — the query that read gigabytes of `tree-genotype` JSON and took
  the block-tree build from 1.5 s to 25 s. `artifact::list_for_alignments_of_kind` exists now; that
  caller was never converted.
- **The canvas has no interaction test coverage, and the rendering rework widened the gap.** The
  click bug that made candidate review unreachable shipped in phase 2 and the layout tests could not
  have caught it: they exercise the pure `layout()` function, which has no input handling. The
  Williamson rework then added member-box clicks and a selectable breadcrumb on top of that, still
  untested — all 13 tests are pure layout. Any test of click routing would need to drive `egui`
  directly.
- **`blocktree_recentre` shipped pinned.** It was set on tab open but cleared only on the empty-tree
  early return, so `horizontal_scroll_offset` re-applied every frame and overwrote any sideways drag.
  Fixed, but it is the second one-shot-flag bug in this view and neither was catchable without the
  interaction coverage above.

## 12. Re-evaluating the workspace when a new tree lands — ✅ built

Placement is **demand-driven per alignment**: `assign_y_haplogroup` computes a
`f:<file hash>|yt:<tree hash>` fingerprint and skips re-scoring when it is unchanged, so a new tree
is only noticed the next time that one subject is re-analysed. Nothing swept the workspace, and
**7 distinct tree generations** had accumulated in `haplogroup_call`.

The sweep itself already existed — `rebuild-signatures` re-places a set of subjects, and the
per-alignment skip makes it cheap on anything already current. What was missing was a **selector**.

### Two independent symptoms

Keying on call fingerprints alone is not enough, and the live workspace shows why. `GMWOF5428705`
holds a call placed against *today's* tree (`E-C116698`) sitting under a consensus of
`E-FT400514:n0`, last reconciled four weeks earlier. **A consensus is derived and persisted
separately, with no tree stamp of its own, so it rots while every call beneath it stays current.**
So `--stale-tree` unions two selectors:

1. `haplogroup_call::biosamples_placed_against_another_tree` — a *source call* stamped with a
   different tree hash. Works for Y (`yt:`) and mt (`mt:`).
2. `App::subjects_labelled_off_tree` — a *derived consensus* naming a branch the current tree does
   not carry. Tested against the tree's `name_index`, so it needs no schema change and catches the
   defect directly: a label absent from the tree is stale by definition, and it is exactly the set
   the block tree drops into `unplaced`.

### The unfingerprinted backlog is a separate job

15,648 Y calls (80%) predate the fingerprint field, so which tree they used is unknowable. Folding
them into the default would make the routine sweep 13,183 subjects — mostly BAM re-walks — for a
tree change that provably affects far fewer. They are opt-in behind `--include-unknown`.

Measured on the live workspace (current tree `yt:b6dfde928041fe28`):

| | current | provably superseded | no fingerprint |
|---|---:|---:|---:|
| Y calls | 198 | 3,780 | 15,648 |
| mt calls | 319 | 3,416 | 0 |

`--stale-tree` selects **5,363** subjects; with `--include-unknown`, 13,183. All 162 `:n0` subjects
of §11.5 are inside the default set, and re-placing one turned `E-FT400514:n0` into `E-C116698`.

### Still open

- **The off-tree-label selector is Y-only.** The fingerprint selector covers both arms (3,416 mt
  calls are on a superseded tree), but `subjects_labelled_off_tree` tests Y consensus labels against
  the Y tree and nothing checks the mt equivalent. Whether mt consensus labels rot the same way is
  unmeasured, not established as safe.
- **No GUI trigger.** This is the CLI half. It is the same question as §11.4 (the private-Y batch)
  and should get the same answer.
- **Nothing *notices* a new tree.** `fetch_tree` refreshes on a TTL, but the sweep is still something
  a user has to think to run. A startup comparison of the tree hash against the workspace's
  most-common stamp would turn this into the "option when a new tree lands" it is meant to be.

## 13. Lifting the chrY structural masks off CHM13

The three curated chrY structural BEDs (amplicons, inverted repeats/palindromes, AZF-DYZ) are
CHM13-native, and both private-Y paths simply skipped them for any other build. That is the cause of
§11.7: a GRCh38 source had **no structural mask at all**, so paralogous sequence — precisely the
sequence that generates spurious novel calls — was counted as unique.

`y_structural_regions_for(build)` now lifts them. The chains were already registered
(`chm13v2-grch38.chain`, `chm13v2-hg19.chain`) and `resolve_chain` already cached them; what was
missing was `Gateway::lift_intervals`.

**PAR and heterochromatin are not lifted.** They are taken as per-build constants, because a chain
is least trustworthy in exactly those places — PAR is shared with chrX and Yq12 is satellite — and
because both are precisely documented per assembly. The palindromes and amplicons *are* lifted,
since they sit in male-specific euchromatin where the chain holds.

**A partial lift is recovered rather than dropped.** Endpoints are tried first, as the exact case.
When they fail — common for amplicons, which are exactly where the assemblies disagree — the
interval body is sampled at 64 points and the dominant target contig's extent is taken, provided at
least a quarter of the samples mapped. This matters because the masks *suppress* calls: too small a
mask admits false novels by the hundred, while too large a one costs a handful of true calls in
known-paralogous sequence. The asymmetry says recover.

Measured survival against the real chains:

| | palindrome | amplicon |
|---|---:|---:|
| CHM13 (native) | 6.19 Mb | 4.41 Mb |
| GRCh38 | 5.65 Mb (91%) | 3.59 Mb (81%) |
| GRCh37 | 5.69 Mb (92%) | 3.64 Mb (83%) |

Endpoint-only lifting gave 79% / 64% on GRCh38; the interior fallback is what recovers the rest. The
remaining shortfall is genuinely non-syntenic — CHM13 carries Y sequence the older assemblies lack.

### Measured on the cohort

Recomputed across R1b-CTS4466Plus (337 buckets). Both caches were **version-bumped** rather than
force-cleared — `private_y` artefact `"3"`→`"4"` and the VCF key `pv1:`→`pv2:` — because the
algorithm changed, and because `--force` cannot reach the VCF cache at all. That left both the old
and new buckets in the database for the same 107 variant sets, so the effect is exactly measurable:

| | before | after |
|---|---:|---:|
| novel-in-unique, total | 7,841 | 5,732 |
| median per donor | 73 | 55 |
| max | 102 | 82 |
| donors over `PRIVATE_Y_QC_WARN` | 106 | 64 |

**26.9% of the calls were paralogous sequence.** Real, and not the whole story.

## 14. The residual is recurrent, not private — open

The 5,732 surviving calls fall on **356 distinct positions**. 92% of them sit at a position seen in
ten or more of the 107 donors, and three positions (12,023,935 · 15,285,765 · 15,732,901) appear in
**every** donor. A variant carried by the entire cohort is not private by any definition — it is
either a real branch variant the tree should name, or a systematic artefact. Either way it must be
blocked, and something already exists to block it.

The cause is asset density, not logic. Both cohort masks apply on GRCh38, but they lifted very
differently:

| | CHM13 | GRCh38 | |
|---|---:|---:|---|
| `chrY_callable_mask` | 14.96 Mb | 14.56 Mb | 97% — fine |
| `chrY_cohort_shared_sites` | 323,414 pos | 105,609 pos | **33%** |

The shared-sites blocklist *is* the recurrent-artefact filter ("every position that varies with ≥2
carriers across the cohort, plus homoplasy hotspots"), and two thirds of it did not survive the
CrossMap lift to GRCh38. Those are single-base intervals, so there is no interior to recover the way
§13 recovers a span — both endpoints have to map or the position is lost.

Re-lifting with our own chain would likely yield about the same. The right fix is to **compute the
GRCh38 blocklist natively from the cohort's own GRCh38 calls** rather than lifting a CHM13 product —
an asset-pipeline job under `scripts/`, and a deliberate one, since the masks are manifest-verified
(`asset-manifest-verification`). Until then, GRCh38 private-Y counts stay roughly an order of
magnitude above FTDNA's, and the block tree's QC suppression (§11, amber edge) is what stands between
that and a misleading diagram.
