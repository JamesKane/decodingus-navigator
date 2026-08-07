//! The project **block tree** view (`impl NavigatorApp`) — the cohort haplotree for the open
//! project, drawn the way Alex Williamson's Big Tree draws it (the presentation FTDNA's Block Tree
//! borrowed): **top-down**, depth increasing downward, and each block showing *its equivalent SNPs*
//! rather than a count of them.
//!
//! That last point is the whole idea. A block is the run of phylogenetically equivalent mutations on
//! a branch — the order within it is unknowable — so the SNP list **is** the block, and printing
//! "17 SNPs" withholds exactly what the view exists to show. The members move out to a roster beside
//! the tree, as the Big Tree puts them in a table below it: the diagram carries the phylogeny, the
//! roster carries the men.
//!
//! The aggregate is built off the UI thread (`App::project_block_tree`, see
//! `documents/design/project-block-tree.md`); this module only lays it out and paints it.
//!
//! The backbone above the cohort is a **breadcrumb, not a block**. The Big Tree's subclade pages do
//! the same: `R-P312/S116 > Z46577 > Z290 > L21/S145 > … > CTS4466/S1136` runs as a path across the
//! top, and the diagram starts at the clade in view. Without that, a cohort whose induced root folds
//! a thousand-SNP backbone opens on one absurd box that is all of the canvas and none of the cohort.
//!
//! Two performance rules, because a group project can hold thousands of members:
//!
//! - **Layout is computed once per (tree, expansion, zoom)**, not rebuilt per frame. [`layout`] is a
//!   pure function over `&[Block]`, so it is testable without a canvas.
//! - **Drawing is culled to `clip_rect`.** Only blocks actually on screen are painted, so a tree
//!   with thousands of blocks costs the same per frame as one with a dozen.

use std::collections::HashMap;

use navigator_app::Block;

use super::*;

// Canvas geometry, at zoom 1.0.
// Sized to the content: a SNP name is short (`A2594`), and the longest are position-derived
// (`14405732-C-T`). A wide box multiplied across hundreds of blocks is what turns the canvas into
// empty space.
const BOX_W: f32 = 84.0;
const ROW_H: f32 = 12.0; // one line of SNP text inside a block
const H_GAP: f32 = 6.0; // horizontal gap between sibling subtrees
/// Vertical gap between a block and its children — **zero**. In the Big Tree a parent block spans
/// the full width of its descendants and they sit flush against its underside, so *containment*
/// carries the parent/child relation and no connector is drawn between levels. A gap here would also
/// corrupt the vertical scale, which is meant to read as accumulated mutations and nothing else.
const V_GAP: f32 = 0.0;
/// Stem length from the last block down to the band of biosample boxes.
const STEM: f32 = 26.0;
/// Ruler tick interval, in SNPs.
const TICK_SNPS: usize = 5;
/// Width of the left gutter carrying the SNP ruler.
const GUTTER_W: f32 = 30.0;
const PAD: f32 = 4.0;
/// A man's box, and the width of a private-variant block. Wide enough that "Private variants" sets
/// on one line and a long kit name (`GMWOF5428705`) is not cropped.
const MEMBER_W: f32 = 84.0;
// No SNP cap: a block's height **is** its elapsed time. Mutations accumulate at a roughly steady
// rate, so the number of phylogenetically equivalent SNPs on a branch is how long that branch ran
// unbroken — and eliding any of them shortens the box, which is to say it misreports the time. A
// line may still carry several *names* (synonyms for one mutation, as `BY30547 Y43043` is one SNP
// with two names); it never carries two mutations.
//
// This is affordable because the one pathological case is handled elsewhere: the backbone above the
// cohort is a breadcrumb, not a block (see `upstream_breadcrumb`).

// Muted, close to the Big Tree's tan-on-parchment but keyed for a dark theme.
const BLOCK_BG: egui::Color32 = egui::Color32::from_rgb(44, 46, 51);
const BLOCK_BG_PLACED: egui::Color32 = egui::Color32::from_rgb(48, 61, 52); // carries members
/// A candidate branch reads as *provisional*: amber, not the green of a published branch. It is an
/// inference from shared private variants, and must never be mistaken for a named haplogroup.
const BLOCK_BG_CANDIDATE: egui::Color32 = egui::Color32::from_rgb(66, 57, 38);
const CANDIDATE_STROKE: egui::Color32 = egui::Color32::from_rgb(190, 148, 70);
const BLOCK_STROKE: egui::Color32 = egui::Color32::from_rgb(78, 84, 94);
const SELECTED_STROKE: egui::Color32 = egui::Color32::from_rgb(120, 170, 220);
const EDGE: egui::Color32 = egui::Color32::from_rgb(72, 78, 88);
const SNP_FG: egui::Color32 = egui::Color32::from_rgb(176, 182, 192);
/// Men are grey against the tree's colour, as the Big Tree draws them — they are the evidence the
/// phylogeny is built from, not part of the phylogeny.
const MEMBER_BG: egui::Color32 = egui::Color32::from_rgb(58, 60, 66);
const MEMBER_FG: egui::Color32 = egui::Color32::from_rgb(198, 202, 210);
/// Private variants get their own colour because they are a different *kind* of claim: unnamed
/// mutations this cohort observed, not branches a tree has published. Teal, as the Big Tree has it.
const PRIVATE_BG: egui::Color32 = egui::Color32::from_rgb(30, 74, 72);
const PRIVATE_STROKE: egui::Color32 = egui::Color32::from_rgb(64, 142, 136);
const PRIVATE_FG: egui::Color32 = egui::Color32::from_rgb(150, 208, 200);

/// One laid-out block: where it sits, and how many SNP names its box shows.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Placed {
    /// Index into the `blocks` slice that was laid out.
    pub idx: usize,
    pub rect: egui::Rect,
}

/// One laid-out man, hanging off the bottom of the block he is placed on.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PlacedMember {
    /// Index into the `blocks` slice.
    pub block: usize,
    /// Index into `blocks[block].members`.
    pub member: usize,
    pub rect: egui::Rect,
}

/// Render a mean over a handful of men: whole numbers plain, otherwise one decimal. Rounding 4.5 to
/// 5 would hide that a branch sits half a mutation from its neighbour; printing `4.0` for an exact 4
/// is just noise.
fn fmt_average(v: f32) -> String {
    if (v - v.round()).abs() < 0.05 {
        format!("{}", v.round() as i64)
    } else {
        format!("{v:.1}")
    }
}

/// The private-variant block below a branch: the mutations its men carry that no branch names yet.
///
/// The mean is over the men whose terminal **is** this block — not its subtree. A branch that both
/// splits and holds men counts only the men standing on it, which is what FTDNA reports too
/// (`R-FGC29071` averages over 2 participants while 7 more sit on branches below it).
///
/// It is drawn **on the same vertical scale as the blocks**, because it measures the same thing —
/// mutations accrued since the named branch above it, which is the time between that branch and the
/// present. That is what makes it belong in the diagram rather than in a tooltip: the ruler reads
/// straight through it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PlacedPrivate {
    /// Index into the `blocks` slice.
    pub block: usize,
    /// Mean private-variant count across the men here that have one.
    pub average: f32,
    /// How many men that mean is over — `private_novel` is `None` when never computed, which is not
    /// zero, so an average over 2 of 30 men must not read as the branch's.
    pub counted: usize,
    /// Men dropped as implausible (see [`PRIVATE_Y_QC_WARN`]). Never silently: the block is marked
    /// and the hover says how many, because "we excluded a third of this branch" is a finding.
    pub suppressed: usize,
    pub rect: egui::Rect,
}

/// One ruler graduation: a y on the canvas and the mutations accumulated to it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Tick {
    pub y: f32,
    pub snps: usize,
}

/// The laid-out tree: block rects, member rects, the SNP ruler, and the total canvas size.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Layout {
    pub placed: Vec<Placed>,
    pub members: Vec<PlacedMember>,
    pub privates: Vec<PlacedPrivate>,
    /// Graduations for the left gutter, along the lineage that accrued the most mutations.
    pub ticks: Vec<Tick>,
    pub size: egui::Vec2,
}

/// Graduations for the SNP ruler, walking the **deepest lineage** — the one that accrued the most
/// mutations, and so the one that reaches furthest down the canvas.
///
/// The ticks are computed rather than spaced evenly, because evenly spaced would be wrong: each
/// block spends one row on its name, so a fixed pixels-per-SNP scale drifts by a row per generation.
/// Walking the lineage and placing each graduation inside the block that contains it keeps the axis
/// honest — the ticks come out *nearly* regular, and where they don't, the irregularity is real.
fn ruler_ticks(blocks: &[Block], placed: &[Placed], row_h: f32, pad: f32) -> Vec<Tick> {
    // Cumulative mutations to the bottom of each block, so "deepest" means most mutations, not most
    // generations — a long slow branch outranks several short ones.
    let index: HashMap<i64, usize> = blocks.iter().enumerate().map(|(i, b)| (b.node_id, i)).collect();
    let mut cum = vec![0usize; blocks.len()];
    let mut best = (0usize, 0usize); // (mutations, block)
    for (i, b) in blocks.iter().enumerate() {
        let above = b.parent.and_then(|p| index.get(&p)).map(|&p| cum[p]).unwrap_or(0);
        cum[i] = above + b.loci.len();
        if cum[i] > best.0 {
            best = (cum[i], i);
        }
    }
    // Walk back up from the deepest block, then read the chain root-first.
    let mut chain = Vec::new();
    let mut at = Some(best.1);
    while let Some(i) = at {
        chain.push(i);
        at = blocks[i].parent.and_then(|p| index.get(&p)).copied();
    }
    chain.reverse();

    let mut ticks = Vec::new();
    let mut seen = 0usize;
    for i in chain {
        let n = blocks[i].loci.len();
        // The SNP rows start below the name row.
        let body_top = placed[i].rect.top() + pad + row_h;
        // Every multiple of TICK_SNPS that falls inside this block's run of mutations.
        let mut k = (seen / TICK_SNPS + 1) * TICK_SNPS;
        while k <= seen + n {
            ticks.push(Tick {
                y: body_top + (k - seen) as f32 * row_h,
                snps: k,
            });
            k += TICK_SNPS;
        }
        seen += n;
    }
    ticks
}

/// Lines a block's box needs: the branch name, one line per equivalent SNP, and the member count.
fn lines_for(b: &Block) -> usize {
    // The name's row + one row per SNP. Every SNP, always — see the note above. The old member-count
    // row is gone: the men are boxes in the band below, so counting them here was both redundant and
    // a row of height that no mutation paid for.
    1 + b.loci.len()
}

/// The folded backbone above the cohort, as a path rather than a box.
///
/// Returns `(path, snps)` when the induced root is a *collapsed run* — a chain of branches the
/// cohort descends through, folded into one block because no split within it separates any two
/// members. R1b-CTS4466Plus opens on `R-Z290`, which is 24 folded branches and 1,763 SNPs: 25
/// branch-lengths of backbone that would be twelve times taller than the cohort hanging off it.
///
/// The test is the *fold*, not whether men sit on it. A collapsed run is by construction more than
/// one branch, so its height is a sum across the tree above the cohort rather than one branch's
/// elapsed time — the one place where height-as-time doesn't hold. A root that was never collapsed
/// is a single genuine branch and stays in the canvas at full height like any other.
///
/// Men parked on the backbone (shallow kits, typically) keep their roster: the breadcrumb selects
/// the block, so nothing is lost but the box.
pub(crate) fn upstream_breadcrumb(blocks: &[Block]) -> Option<(String, usize)> {
    let root = blocks.iter().find(|b| b.parent.is_none())?;
    if root.collapsed.is_empty() {
        return None;
    }
    // `collapsed` is root-most first and the surviving block keeps the deepest name, so appending it
    // reads oldest → youngest, the direction the breadcrumb is travelled.
    let mut path = root.collapsed.clone();
    path.push(root.name.clone());
    Some((path.join("  ›  "), root.loci.len()))
}

/// Lay `blocks` (in pre-order, as [`ProjectBlockTree`] delivers them) onto a canvas: **depth → y**,
/// tidy-tree order → **x**, root at the top. A parent is centred over the horizontal extent of its
/// children, so a branch point sits above the lineages it splits into.
///
/// A block hangs **directly beneath its parent**, not on a row shared with everything at its depth.
/// That makes vertical position cumulative: how far down a block sits is the mutations accumulated
/// along the path to it, so the y axis reads as elapsed time the same way a box's height does.
/// Aligning depths into rows would instead pad every short branch out to the tallest box beside it,
/// which is both a lot of empty canvas and a lie about when the branch happened.
///
/// Pure: no `Ui`, no state. Extents bottom-up, then positions top-down.
pub(crate) fn layout(blocks: &[Block], zoom: f32) -> Layout {
    if blocks.is_empty() {
        return Layout::default();
    }
    let zoom = zoom.clamp(0.4, 2.5);
    let (box_w, row_h, h_gap, v_gap, pad) = (BOX_W * zoom, ROW_H * zoom, H_GAP * zoom, V_GAP * zoom, PAD * zoom);
    let member_w = MEMBER_W * zoom;
    let member_h = row_h + 2.0 * pad;
    let (stem, gutter) = (STEM * zoom, GUTTER_W * zoom);

    let index: HashMap<i64, usize> = blocks.iter().enumerate().map(|(i, b)| (b.node_id, i)).collect();
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); blocks.len()];
    let mut roots: Vec<usize> = Vec::new();
    for (i, b) in blocks.iter().enumerate() {
        match b.parent.and_then(|p| index.get(&p)) {
            Some(&p) => children[p].push(i),
            None => roots.push(i),
        }
    }

    let heights: Vec<f32> = blocks.iter().map(|b| lines_for(b) as f32 * row_h + 2.0 * pad).collect();
    let mut member_slots: Vec<(usize, usize, f32)> = Vec::new();

    // Pass 1, bottom-up: the horizontal extent each subtree needs. `blocks` is pre-order, so
    // iterating in reverse visits every child before its parent.
    // A man occupies a slot beside his block's child subtrees: the Big Tree hangs him off the bottom
    // of his terminal on a stem of his own, so he needs horizontal room like a subtree does.
    let slots = |i: usize| children[i].len() + blocks[i].members.len();
    let mut extent = vec![box_w; blocks.len()];
    for i in (0..blocks.len()).rev() {
        let n = slots(i);
        if n == 0 {
            continue;
        }
        let kids: f32 = children[i].iter().map(|&c| extent[c]).sum::<f32>();
        let men = blocks[i].members.len() as f32 * member_w;
        extent[i] = extent[i].max(kids + men + h_gap * (n - 1) as f32);
    }

    // Pass 2, top-down: hand each subtree its band, then let the block *fill* it.
    let mut left = vec![0.0f32; blocks.len()];
    let mut cursor = gutter;
    for &r in &roots {
        left[r] = cursor;
        cursor += extent[r] + h_gap;
    }
    let mut placed = Vec::with_capacity(blocks.len());
    let mut deepest = 0.0f32;
    // Pre-order, so a parent's top is always settled before its children read it.
    let mut top = vec![0.0f32; blocks.len()];
    // Cumulative SNPs down to each block's top — the quantity the ruler measures.
    let mut snps_above = vec![0usize; blocks.len()];
    for i in 0..blocks.len() {
        // Icicle: the block spans its whole subtree. A parent therefore visibly *contains* the
        // lineages it splits into, which is how the Big Tree shows descent — no elbow needed.
        let rect = egui::Rect::from_min_size(egui::pos2(left[i], top[i]), egui::vec2(extent[i], heights[i]));
        for &c in &children[i] {
            top[c] = rect.bottom() + v_gap;
            snps_above[c] = snps_above[i] + blocks[i].loci.len();
        }

        let n = slots(i);
        if n > 0 {
            let kids: f32 = children[i].iter().map(|&c| extent[c]).sum::<f32>();
            let total = kids + blocks[i].members.len() as f32 * member_w + h_gap * (n - 1) as f32;
            let mut cx = left[i] + (extent[i] - total) / 2.0;
            // Men take a slot to the left of the subclades, so a lineage that both splits and holds
            // men makes room for both. Their boxes are positioned later, once the band is known.
            for m in 0..blocks[i].members.len() {
                member_slots.push((i, m, cx));
                cx += member_w + h_gap;
            }
            for &c in &children[i] {
                left[c] = cx;
                cx += extent[c] + h_gap;
            }
        }
        deepest = deepest.max(rect.bottom());
        placed.push(Placed { idx: i, rect });
    }

    // Private variants, between a branch and its men — flush under the block, on the same scale, so
    // the ruler measures straight through. The span is the men's, not the block's: these mutations
    // belong to the men standing here, not to the subclades that branch off elsewhere under it.
    let mut privates: Vec<PlacedPrivate> = Vec::new();
    for (i, b) in blocks.iter().enumerate() {
        // A donor whose raw novel count trips the workspace's own plausibility threshold is dropped
        // whole, not trimmed. `PRIVATE_Y_QC_WARN` already declares such a count "unusually high for
        // one sample — check for contamination, low/uneven coverage, or a reference-build mismatch",
        // which is a statement about the *sample*, so its gated count is not trustworthy either. One
        // donor at 661 would otherwise set a branch's height single-handed.
        let plausible = |m: &&navigator_app::BlockMember| {
            !m.private_novel
                .is_some_and(|n| n >= navigator_domain::results_context::PRIVATE_Y_QC_WARN)
        };
        let counts: Vec<usize> = b
            .members
            .iter()
            .filter(plausible)
            .filter_map(|m| m.private_publishable)
            .collect();
        let suppressed = b.members.iter().filter(|m| !plausible(m)).count();
        if counts.is_empty() {
            continue;
        }
        let average = counts.iter().sum::<usize>() as f32 / counts.len() as f32;
        let xs: Vec<f32> = member_slots
            .iter()
            .filter(|(bi, _, _)| *bi == i)
            .map(|(_, _, x)| *x)
            .collect();
        let (Some(&lo), Some(&hi)) = (
            xs.iter().min_by(|a, b| a.total_cmp(b)),
            xs.iter().max_by(|a, b| a.total_cmp(b)),
        ) else {
            continue;
        };
        // **One column wide, always** — centred over the men it covers. The figure is a single
        // branch-level statistic, so sizing the box to the number of men would imply it is a
        // per-man quantity, and would make an identical average look different on two branches for
        // no reason but headcount.
        let h = (average * row_h).max(row_h) + 2.0 * pad;
        let rect = egui::Rect::from_min_size(
            egui::pos2((lo + hi) / 2.0, placed[i].rect.bottom()),
            egui::vec2(member_w, h),
        );
        deepest = deepest.max(rect.bottom());
        privates.push(PlacedPrivate {
            block: i,
            average,
            counted: counts.len(),
            suppressed,
            rect,
        });
    }

    // The men sit in one band beneath the whole diagram, as the Big Tree tables them below it,
    // reached by a stem from their block. Sharing a baseline is what makes them scannable: hung from
    // their own blocks they would step down the page in lockstep with the phylogeny, which says
    // nothing about the men.
    let band = deepest + stem;
    let members: Vec<PlacedMember> = member_slots
        .into_iter()
        .map(|(block, member, x)| PlacedMember {
            block,
            member,
            rect: egui::Rect::from_min_size(egui::pos2(x, band), egui::vec2(member_w, member_h)),
        })
        .collect();
    if !members.is_empty() {
        deepest = band + member_h;
    }

    Layout {
        ticks: ruler_ticks(blocks, &placed, row_h, pad),
        placed,
        members,
        privates,
        size: egui::vec2((cursor - h_gap).max(0.0) + pad, deepest + pad),
    }
}

impl NavigatorApp {
    /// The project's Tree tab: the cohort Y block tree, loaded lazily on first view.
    pub(crate) fn project_blocktree_section(&mut self, ui: &mut egui::Ui) {
        let Some(pid) = self.selected_project else { return };

        // Lazy load — the aggregate fetches and parses a multi-MB haplotree, so it is not built on
        // project select like the STR chart is.
        if self.project_blocktree.is_none() && !self.project_blocktree_loading {
            self.project_blocktree_loading = true;
            let _ = self.tx.send(Command::LoadProjectBlockTree(pid));
        }

        if self.project_blocktree_loading && self.project_blocktree.is_none() {
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(egui::RichText::new(self.tr("blocktree.building")).weak());
            });
            return;
        }
        let Some(tree) = self.project_blocktree.as_ref() else {
            ui.add_space(8.0);
            ui.label(egui::RichText::new(self.tr("blocktree.empty")).weak());
            return;
        };

        // Summary line: how much of the project the tree actually accounts for. `unplaced` is shown
        // even when zero-length is impossible — a cohort with skew must not look complete.
        let placed: usize = tree.blocks.iter().map(|b| b.members.len()).sum();
        let unplaced = tree.unplaced.len();
        let summary = format!(
            "{} · {} {} · {} {}",
            self.tr("blocktree.summary.prefix"),
            placed,
            self.tr("blocktree.summary.placed"),
            tree.blocks.len(),
            self.tr("blocktree.summary.blocks"),
        );
        // The coordinate space matters: node names are build-independent, the SNP positions are not,
        // so the view says which tree and which build it is showing.
        let coords = if tree.build_key.is_empty() {
            tree.provider.clone()
        } else {
            format!("{} · {}", tree.provider, tree.build_key)
        };
        let unplaced_msg = (unplaced > 0).then(|| {
            format!(
                "{unplaced} {} — {}",
                self.tr("blocktree.summary.unplaced"),
                self.tr("blocktree.unplaced.hint")
            )
        });
        // Candidate branches are the thing a published tree can't tell you, so they get their own
        // line rather than being left for the user to notice in the canvas.
        let candidates = tree.blocks.iter().filter(|b| b.candidate).count();
        let candidate_msg = (candidates > 0).then(|| {
            let mut s = format!("{candidates} {}", self.tr("blocktree.candidates"));
            if tree.candidate_conflicts > 0 {
                s.push_str(&format!(
                    " · {} {}",
                    tree.candidate_conflicts,
                    self.tr("blocktree.conflicts")
                ));
            }
            if tree.candidate_recurrent > 0 {
                s.push_str(&format!(
                    " · {} {}",
                    tree.candidate_recurrent,
                    self.tr("blocktree.recurrent")
                ));
            }
            s
        });
        // Every label is resolved before the closure: `self.tr` borrows `self`, and the zoom slider
        // needs `&mut` — the two can't coexist inside one closure.
        let roster_empty = self.tr("blocktree.roster.empty").to_string();
        let upstream_snps = self.tr("blocktree.upstream.snps").to_string();
        let upstream_hint = self.tr("blocktree.upstream.hint").to_string();
        let private_label = self.tr("blocktree.private.title").to_string();
        let private_average = self.tr("blocktree.private.average").to_string();
        let (zoom_label, no_placement, candidate_label, export_label) = (
            self.tr("blocktree.zoom").to_string(),
            self.tr("blocktree.none.placed").to_string(),
            self.tr("blocktree.candidate").to_string(),
            self.tr("blocktree.export").to_string(),
        );
        let mut zoom = self.blocktree_zoom;
        let mut export = false;

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(summary).weak());
            ui.separator();
            ui.label(egui::RichText::new(coords).weak().small());
            ui.separator();
            ui.label(egui::RichText::new(zoom_label).weak().small());
            ui.add(egui::Slider::new(&mut zoom, 0.5..=2.0).show_value(false));
            ui.separator();
            export = ui.button(export_label).clicked();
        });
        if let Some(msg) = unplaced_msg {
            ui.label(
                egui::RichText::new(msg)
                    .small()
                    .color(egui::Color32::from_rgb(210, 160, 90)),
            );
        }
        if let Some(msg) = candidate_msg {
            ui.label(egui::RichText::new(msg).small().color(CANDIDATE_STROKE));
        }
        if tree.blocks.is_empty() {
            ui.add_space(8.0);
            ui.label(egui::RichText::new(no_placement).weak());
            self.blocktree_zoom = zoom;
            self.blocktree_recentre = false;
            return;
        }
        // The backbone above the cohort becomes a path across the top; the canvas draws what is
        // left, which is the cohort itself.
        let upstream = upstream_breadcrumb(&tree.blocks);
        let drawn: Vec<Block> = if upstream.is_some() {
            let root = tree.blocks.iter().find(|b| b.parent.is_none()).map(|b| b.node_id);
            tree.blocks
                .iter()
                .filter(|b| Some(b.node_id) != root)
                // The lifted root's children become roots; their connector had nowhere to go anyway.
                .map(|b| {
                    let mut b = b.clone();
                    if b.parent == root {
                        b.parent = None;
                    }
                    b
                })
                .collect()
        } else {
            tree.blocks.clone()
        };
        let mut select_upstream = false;
        if let Some((path, snps)) = &upstream {
            // The lineage the cohort descends through, and how many mutations sit on it. It is a
            // path rather than a block because its height would dwarf everything the cohort is.
            // Clickable, so the men parked on the backbone still reach the roster.
            ui.horizontal_wrapped(|ui| {
                ui.spacing_mut().item_spacing.x = 0.0;
                select_upstream |= ui
                    .add(egui::Label::new(egui::RichText::new(path).small().color(SNP_FG)).sense(egui::Sense::click()))
                    .on_hover_text(&upstream_hint)
                    .clicked();
                ui.label(
                    egui::RichText::new(format!("   ({snps} {upstream_snps})"))
                        .small()
                        .weak(),
                );
            });
        }
        ui.add_space(4.0);

        // Roster first, from the right, so the tree takes whatever is left — the Big Tree keeps the
        // men in a table rather than in the diagram, and the diagram needs the room.
        let roster = self
            .blocktree_selected
            .and_then(|id| tree.blocks.iter().find(|b| b.node_id == id))
            .map(|b| {
                let label = if b.candidate {
                    candidate_label.clone()
                } else {
                    b.name.clone()
                };
                let rows: Vec<(String, Option<usize>)> =
                    b.members.iter().map(|m| (m.name.clone(), m.private_novel)).collect();
                (label, rows, b.subtree_members)
            });

        let lay = layout(&drawn, zoom);
        let mut review: Option<i64> = None;
        let mut select: Option<i64> = None;
        let selected = self.blocktree_selected;
        let mut open_subject: Option<SampleGuid> = None;

        if let Some((label, rows, below)) = &roster {
            egui::SidePanel::right("blocktree_roster")
                .resizable(true)
                .default_width(210.0)
                .show_inside(ui, |ui| {
                    ui.label(egui::RichText::new(label).strong());
                    ui.label(
                        egui::RichText::new(format!("{} here · {} below", rows.len(), below))
                            .small()
                            .weak(),
                    );
                    ui.separator();
                    if rows.is_empty() {
                        ui.label(egui::RichText::new(roster_empty.clone()).small().weak());
                        return;
                    }
                    // Virtualized: a terminal block can hold thousands of kits.
                    let line = ui.text_style_height(&egui::TextStyle::Small) + 2.0;
                    egui::ScrollArea::vertical().auto_shrink([false, false]).show_rows(
                        ui,
                        line,
                        rows.len(),
                        |ui, range| {
                            for (name, novel) in &rows[range] {
                                // `None` means private-Y was never computed — not the same as zero,
                                // so it shows nothing rather than "(0)".
                                let text = match novel {
                                    Some(n) if *n > 0 => format!("{name}  ({n})"),
                                    _ => name.clone(),
                                };
                                ui.label(egui::RichText::new(text).small());
                            }
                        },
                    );
                });
        }

        // Centre the first view on the root. The canvas is far wider than any viewport, and its
        // left edge is empty space belonging to subtrees that hang further down.
        let root_x = lay
            .placed
            .iter()
            .find(|p| drawn[p.idx].parent.is_none())
            .map(|p| p.rect.center().x);
        let mut area = egui::ScrollArea::both().auto_shrink([false, false]);
        if let (Some(x), true) = (root_x, self.blocktree_recentre) {
            area = area.horizontal_scroll_offset((x - ui.available_width() / 2.0).max(0.0));
            // One-shot. Left set, `horizontal_scroll_offset` is re-applied on every frame and the
            // canvas snaps back to the root the instant the user drags sideways.
            self.blocktree_recentre = false;
        }
        area.show(ui, |ui| {
            let (canvas, _) = ui.allocate_exact_size(lay.size, egui::Sense::hover());
            let painter = ui.painter_at(canvas);
            let clip = ui.clip_rect();
            let origin = canvas.min.to_vec2();
            let zoom = zoom.clamp(0.4, 2.5);
            let font = egui::FontId::proportional(11.0 * zoom);
            let small = egui::FontId::proportional(9.5 * zoom);

            // The SNP ruler, in the left gutter: the scale that makes a block's height readable as a
            // quantity rather than an impression. Graduated in mutations accumulated from the top of
            // this view — not from the root of the tree, which is above the cohort and in the
            // breadcrumb.
            {
                let g = egui::Rect::from_min_size(
                    egui::pos2(canvas.left(), canvas.top()),
                    egui::vec2(GUTTER_W * zoom, lay.size.y),
                );
                painter.rect_filled(g, 0.0, egui::Color32::from_gray(34));
                for t in &lay.ticks {
                    let y = canvas.top() + t.y;
                    if !clip.y_range().contains(y) {
                        continue;
                    }
                    painter.line_segment(
                        [egui::pos2(g.right() - 4.0 * zoom, y), egui::pos2(g.right(), y)],
                        egui::Stroke::new(1.0_f32, egui::Color32::from_gray(80)),
                    );
                    painter.text(
                        egui::pos2(g.right() - 6.0 * zoom, y),
                        egui::Align2::RIGHT_CENTER,
                        t.snps.to_string(),
                        small.clone(),
                        egui::Color32::from_gray(120),
                    );
                }
            }

            for p in &lay.placed {
                let rect = p.rect.translate(origin);
                let b = &drawn[p.idx];

                // No connector to the parent: the block sits flush under it and inside its span, so
                // containment shows the descent. An elbow here would be drawing what the geometry
                // already says.

                // Cull: everything below is per-block text layout, the expensive part.
                if !clip.intersects(rect) {
                    continue;
                }

                let (bg, stroke) = match (b.candidate, b.members.is_empty()) {
                    (true, _) => (BLOCK_BG_CANDIDATE, egui::Stroke::new(1.5_f32, CANDIDATE_STROKE)),
                    (false, true) => (BLOCK_BG, egui::Stroke::new(1.0_f32, BLOCK_STROKE)),
                    (false, false) => (BLOCK_BG_PLACED, egui::Stroke::new(1.0_f32, BLOCK_STROKE)),
                };
                painter.rect_filled(rect, 3.0, bg);
                let stroke = if selected == Some(b.node_id) {
                    egui::Stroke::new(2.0_f32, SELECTED_STROKE)
                } else {
                    stroke
                };
                painter.rect_stroke(rect, 3.0, stroke);

                let pad = PAD * zoom;
                let row = ROW_H * zoom;
                let mut y = rect.top() + pad;
                let cx = rect.center().x;
                // Clipped to the box: a label that outgrows its block is then cropped rather than
                // spilling across the canvas, whatever the text turns out to be.
                let inner = painter.with_clip_rect(rect.shrink(1.0));
                let put = |text: String, color: egui::Color32, f: &egui::FontId, y: f32| {
                    inner.text(egui::pos2(cx, y), egui::Align2::CENTER_TOP, text, f.clone(), color);
                };

                // A candidate has no published name — the view supplies the label, localized.
                let (title, title_fg) = if b.candidate {
                    (candidate_label.clone(), CANDIDATE_STROKE)
                } else {
                    (b.name.clone(), egui::Color32::from_gray(230))
                };
                put(title, title_fg, &font, y);
                y += row;

                // The equivalent SNPs themselves, one per line — the block's actual content, and the
                // reason its height means something. Printing a count instead withholds both the
                // mutations and the sense of time the box is carrying.
                for l in &b.loci {
                    put(l.name.clone(), SNP_FG, &small, y);
                    y += row;
                }
                // ONE interact per block. Two on the same rect meant the later one sat on top and
                // swallowed the click: every candidate has members, so the double-click handler
                // always existed for them and single-click never fired.
                let resp = ui.interact(rect, egui::Id::new(("blocktree", b.node_id)), egui::Sense::click());
                if resp.double_clicked() {
                    // Jump to a member's subject page.
                    if let Some(m) = b.members.first() {
                        open_subject = Some(m.guid);
                    }
                } else if resp.clicked() {
                    // A named block expands to show its members; a candidate is an inference, so
                    // clicking it opens the evidence instead of just more names.
                    if b.candidate {
                        review = Some(b.node_id);
                    } else {
                        select = Some(b.node_id);
                    }
                }
                if resp.hovered() && !b.loci.is_empty() {
                    let names: Vec<&str> = b.loci.iter().map(|l| l.name.as_str()).take(40).collect();
                    let extra = b.loci.len().saturating_sub(names.len());
                    let mut tip = names.join(", ");
                    if extra > 0 {
                        tip.push_str(&format!(" … +{extra}"));
                    }
                    if !b.collapsed.is_empty() {
                        tip.push_str(&format!(
                            "\n\n{} branch(es) folded in: {}",
                            b.collapsed.len(),
                            b.collapsed.join(" → ")
                        ));
                    }
                    resp.on_hover_text(tip);
                }
            }

            // Private variants, between each branch and its men.
            for pv in &lay.privates {
                let rect = pv.rect.translate(origin);
                if !clip.intersects(rect) {
                    continue;
                }
                painter.rect_filled(rect, 2.0, PRIVATE_BG);
                // A suppressed donor is a fact about the branch, so the box says so rather than
                // quietly reporting a mean over whoever survived.
                let edge = if pv.suppressed > 0 {
                    egui::Stroke::new(1.5_f32, CANDIDATE_STROKE)
                } else {
                    egui::Stroke::new(1.0_f32, PRIVATE_STROKE)
                };
                painter.rect_stroke(rect, 2.0, edge);

                // The box's height is the measurement, so the text has to fit *it* — never the
                // other way round. A one-mutation block is one row tall, which holds one line, and
                // the line that matters is the number: the title is a label, the value is the
                // finding. So the title appears only when both fit, exactly as the Big Tree drops it
                // from its thin blocks.
                let inner = painter.with_clip_rect(rect.shrink(1.0));
                let pad_z = PAD * zoom;
                let avail = rect.width() - 2.0 * pad_z;
                let title = ui.fonts(|f| f.layout(private_label.clone(), small.clone(), PRIVATE_FG, avail));
                let value = ui.fonts(|f| {
                    f.layout(
                        format!("{private_average} {}", fmt_average(pv.average)),
                        small.clone(),
                        PRIVATE_FG,
                        avail,
                    )
                });
                let both = title.size().y + value.size().y <= rect.height() - 2.0 * pad_z;
                let used = if both {
                    title.size().y + value.size().y
                } else {
                    value.size().y
                };
                // Centred in the box, as the reference centres it — with the text top-aligned once
                // the box is shorter than the text, so what survives the clip is the start of it.
                let mut y = rect.top() + ((rect.height() - used) / 2.0).max(pad_z);
                if both {
                    inner.galley(
                        egui::pos2(rect.center().x - title.size().x / 2.0, y),
                        title.clone(),
                        PRIVATE_FG,
                    );
                    y += title.size().y;
                }
                inner.galley(
                    egui::pos2(rect.center().x - value.size().x / 2.0, y),
                    value.clone(),
                    PRIVATE_FG,
                );

                let resp = ui.interact(
                    rect,
                    egui::Id::new(("blocktree_private", pv.block)),
                    egui::Sense::hover(),
                );
                if resp.hovered() {
                    let b = &drawn[pv.block];
                    // Name the denominator. It is the men whose terminal *is* this block — not the
                    // subtree — and among those, only the ones private-Y has actually been computed
                    // for, since `private_novel` is `None` until then and `None` is not zero.
                    let mut tip = format!(
                        "{}\nOn average {} publishable private variant(s) in {} of {} men placed here",
                        b.name,
                        fmt_average(pv.average),
                        pv.counted,
                        b.members.len()
                    );
                    if pv.suppressed > 0 {
                        tip.push_str(&format!(
                            "\n\n{} man/men excluded: over {} novel calls each, which is unusually high \
                             for one sample — check for contamination, low/uneven coverage, or a \
                             reference-build mismatch",
                            pv.suppressed,
                            navigator_domain::results_context::PRIVATE_Y_QC_WARN
                        ));
                    }
                    let unaccounted = b.members.len() - pv.counted - pv.suppressed;
                    if unaccounted > 0 {
                        tip.push_str(&format!("\n{unaccounted} with no private-Y computed"));
                    }
                    resp.on_hover_text(tip);
                }
            }

            // Men, as the Big Tree draws them: a grey box on a stem below the block they sit on.
            // They are the evidence the phylogeny rests on, so they belong in the diagram — the
            // roster beside it stays for the private-variant counts and for scanning a long list.
            for pm in &lay.members {
                let rect = pm.rect.translate(origin);
                let b = &drawn[pm.block];
                // A stem from the block down to the band. This is the one connector the diagram
                // still draws, because a man's box is the one thing not positioned by containment.
                // Start the stem below the private-variant block when there is one, so the man hangs
                // off his own unnamed mutations rather than appearing to hang off the named branch.
                let from = lay
                    .privates
                    .iter()
                    .find(|pv| pv.block == pm.block)
                    .map(|pv| pv.rect)
                    .unwrap_or(lay.placed[pm.block].rect)
                    .translate(origin);
                let a = egui::pos2(rect.center().x.clamp(from.left(), from.right()), from.bottom());
                let z = egui::pos2(rect.center().x, rect.top());
                if clip.intersects(egui::Rect::from_two_pos(a, z)) {
                    let stroke = egui::Stroke::new(1.0_f32, EDGE);
                    let mid = z.y - (z.y - a.y) * 0.35;
                    painter.line_segment([a, egui::pos2(a.x, mid)], stroke);
                    painter.line_segment([egui::pos2(a.x, mid), egui::pos2(z.x, mid)], stroke);
                    painter.line_segment([egui::pos2(z.x, mid), z], stroke);
                }
                if !clip.intersects(rect) {
                    continue;
                }
                let m = &b.members[pm.member];
                painter.rect_filled(rect, 2.0, MEMBER_BG);
                painter.rect_stroke(rect, 2.0, egui::Stroke::new(1.0_f32, BLOCK_STROKE));
                painter.with_clip_rect(rect.shrink(1.0)).text(
                    rect.center(),
                    egui::Align2::CENTER_CENTER,
                    &m.name,
                    small.clone(),
                    MEMBER_FG,
                );

                let resp = ui.interact(
                    rect,
                    egui::Id::new(("blocktree_member", pm.block, pm.member)),
                    egui::Sense::click(),
                );
                if resp.clicked() {
                    open_subject = Some(m.guid);
                }
                if resp.hovered() {
                    // `None` means private-Y was never computed, which is not the same as zero.
                    let tip = match m.private_novel {
                        Some(n) => format!("{}\n{n} private variant(s)", m.name),
                        None => m.name.clone(),
                    };
                    resp.on_hover_text(tip);
                }
            }
        });

        // Formatted while the tree borrow is still live; written below, after it ends.
        let export_bodies = export.then(|| {
            let name = self
                .overview
                .iter()
                .find(|o| o.project.id == pid)
                .map(|o| o.project.name.clone())
                .unwrap_or_else(|| format!("project-{pid}"));
            let safe: String = name
                .chars()
                .map(|c| if c.is_ascii_alphanumeric() || c == '-' { c } else { '_' })
                .collect();
            (
                navigator_app::export::block_tree_tsv(tree),
                navigator_app::export::block_tree_html(tree, &name),
                safe,
            )
        });

        // Deferred dispatch: mutate state after the closure, never inside it.
        self.blocktree_zoom = zoom;
        if let Some((tsv, html, safe)) = export_bodies {
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name(format!("blocktree_{safe}.tsv"))
                .add_filter("TSV", &["tsv"])
                .add_filter("HTML", &["html"])
                .save_file()
            {
                // One button, either format — chosen by the extension the user typed, as the other
                // two-format exports in this app do.
                let html_wanted = path.extension().is_some_and(|e| e.eq_ignore_ascii_case("html"));
                let body = if html_wanted { html } else { tsv };
                match std::fs::write(&path, body) {
                    Ok(()) => self.status = format!("{} {}", self.tr("blocktree.exported"), path.display()),
                    Err(e) => self.status = format!("write {}: {e}", path.display()),
                }
            }
        }
        if select_upstream {
            self.blocktree_selected = tree.blocks.iter().find(|b| b.parent.is_none()).map(|b| b.node_id);
        }
        if let Some(id) = select {
            self.blocktree_selected = Some(id);
        }
        if let Some(id) = review {
            self.blocktree_review = Some(id);
        }
        if let Some(guid) = open_subject {
            self.return_to_project = Some(pid);
            self.select_sample(guid);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_app::{Block, BlockMember};

    /// Layout depends only on the *number* of members, so the guid can be the nil one — this keeps
    /// `uuid` out of the UI crate's dependencies.
    fn member(name: &str) -> BlockMember {
        BlockMember {
            guid: SampleGuid(Default::default()),
            name: name.into(),
            private_novel: None,
            private_publishable: None,
            private_total: None,
        }
    }

    fn locus(name: &str, position: i64) -> navigator_app::Locus {
        navigator_app::Locus {
            position,
            ancestral: "A".into(),
            derived: "G".into(),
            name: name.into(),
        }
    }

    fn block(id: i64, parent: i64, depth: usize, members: &[&str]) -> Block {
        Block {
            node_id: id,
            name: format!("N{id}"),
            parent: (parent != 0).then_some(parent),
            depth,
            loci: vec![locus("M1", 1)],
            members: members.iter().map(|m| member(m)).collect(),
            subtree_members: members.len(),
            collapsed: Vec::new(),
            candidate: false,
            evidence: Vec::new(),
        }
    }

    /// `R ─┬─> R1(a)`
    ///     `└─> R2(b)`
    fn split() -> Vec<Block> {
        vec![block(1, 0, 0, &[]), block(2, 1, 1, &["a"]), block(3, 1, 1, &["b"])]
    }

    #[test]
    fn layout_places_depth_on_y() {
        // Top-down, as the Big Tree draws it: the root sits above its children.
        let lay = layout(&split(), 1.0);
        assert_eq!(lay.placed.len(), 3);
        assert_eq!(lay.placed[0].rect.top(), 0.0, "the root starts the canvas");
        assert!(lay.placed[1].rect.top() >= lay.placed[0].rect.bottom());
    }

    /// Vertical position is cumulative, so a lineage that accrued more mutations sits lower than its
    /// cousin at the same depth. Aligning depths into rows would flatten exactly that difference.
    #[test]
    fn a_lineage_that_accrued_more_snps_sits_lower() {
        let mut root = block(1, 0, 0, &[]);
        root.subtree_members = 2;
        let mut slow = block(2, 1, 1, &["a"]);
        let mut fast = block(3, 1, 1, &["b"]);
        fast.loci = (0..30).map(|i| locus(&format!("M{i}"), i)).collect();
        let kid_of_slow = block(4, 2, 2, &["c"]);
        let kid_of_fast = block(5, 3, 2, &["d"]);
        slow.subtree_members = 2;
        fast.subtree_members = 2;

        let lay = layout(&[root, slow, fast, kid_of_slow, kid_of_fast], 1.0);
        let (under_slow, under_fast) = (lay.placed[3].rect.top(), lay.placed[4].rect.top());
        assert!(
            under_fast > under_slow + 20.0 * ROW_H,
            "29 extra mutations upstream must push the branch visibly further down"
        );
    }

    /// Icicle: the parent spans its whole subtree, so descent reads as containment.
    #[test]
    fn a_parent_block_spans_its_children() {
        let lay = layout(&split(), 1.0);
        let (root, a, b) = (lay.placed[0].rect, lay.placed[1].rect, lay.placed[2].rect);
        assert!(a.center().x < b.center().x, "siblings run left to right in pre-order");
        assert!(
            root.left() <= a.left() && root.right() >= b.right(),
            "the parent must contain both children horizontally"
        );
        assert!(
            (a.top() - root.bottom()).abs() < 0.01,
            "children sit flush against the parent — no gap to misread as elapsed time"
        );
    }

    #[test]
    fn layout_siblings_do_not_overlap() {
        let lay = layout(&split(), 1.0);
        assert!(
            lay.placed[1].rect.right() <= lay.placed[2].rect.left(),
            "sibling boxes must not overlap"
        );
    }

    /// Height is the Big Tree's proxy for elapsed time, so it must track the SNP count.
    #[test]
    fn a_block_is_as_tall_as_it_has_snps() {
        let mut short = block(1, 0, 0, &[]);
        short.loci = (0..2).map(|i| locus(&format!("M{i}"), i)).collect();
        let mut tall = block(1, 0, 0, &[]);
        tall.loci = (0..20).map(|i| locus(&format!("M{i}"), i)).collect();

        let sh = layout(&[short], 1.0).placed[0].rect.height();
        let th = layout(&[tall], 1.0).placed[0].rect.height();
        assert!(th > sh, "20 equivalent SNPs must stand taller than 2");
        // One line per SNP: the 18 extra mutations are 18 extra rows.
        assert!(
            (th - sh - 18.0 * ROW_H).abs() < 0.5,
            "height tracks the SNP count exactly"
        );
    }

    /// No cap, at any size. A truncated box is a shortened box, and a shortened box is a shorter
    /// span of time than the branch actually ran.
    #[test]
    fn a_large_block_is_never_truncated() {
        let mut b = block(1, 0, 0, &[]);
        b.loci = (0..600).map(|i| locus(&format!("M{i}"), i)).collect();
        let lay = layout(&[b], 1.0);
        // name + 600 SNPs.
        assert!((lay.placed[0].rect.height() - (601.0 * ROW_H + 2.0 * PAD)).abs() < 0.5);
    }

    /// The backbone the cohort merely passed through is upstream context, so it leaves the canvas
    /// for the breadcrumb — otherwise one member-less box is taller than the whole cohort below it.
    #[test]
    fn a_folded_backbone_root_becomes_a_breadcrumb() {
        let mut root = block(1, 0, 0, &[]);
        root.name = "Z290".into();
        root.collapsed = vec!["P312".into(), "L21".into()];
        root.loci = (0..900).map(|i| locus(&format!("M{i}"), i)).collect();

        let (path, snps) = upstream_breadcrumb(&[root.clone(), block(2, 1, 1, &["a"])]).expect("backbone lifts out");
        assert_eq!(path, "P312  ›  L21  ›  Z290", "oldest to youngest");
        assert_eq!(snps, 900);

        // The live cohort's backbone carries two shallow kits, so men on it must not keep it in the
        // canvas — the fold is what makes it upstream.
        root.members = vec![member("a")];
        assert!(upstream_breadcrumb(&[root]).is_some());
    }

    #[test]
    fn an_uncollapsed_root_stays_in_the_canvas() {
        // One genuine branch: its height is one branch's elapsed time, which is the whole point.
        assert!(upstream_breadcrumb(&split()).is_none());
    }

    /// Men are blocks of their own hanging under their terminal, the way the Big Tree stems them.
    #[test]
    fn men_share_one_band_below_the_diagram() {
        let mut blocks = split();
        // Give the two lineages different depths, so a shared baseline is observable.
        blocks[2].loci = (0..20).map(|i| locus(&format!("M{i}"), i)).collect();
        let lay = layout(&blocks, 1.0);
        assert_eq!(lay.members.len(), 2, "one box per man");
        assert!(
            (lay.members[0].rect.top() - lay.members[1].rect.top()).abs() < 0.01,
            "men are tabled on one baseline, not stepped down with the phylogeny"
        );
        for pm in &lay.members {
            let block = lay.placed.iter().find(|p| p.idx == pm.block).unwrap();
            assert!(pm.rect.top() > block.rect.bottom(), "a man hangs below his block");
        }
        assert_ne!(
            lay.members[0].block, lay.members[1].block,
            "each under his own terminal"
        );
    }

    /// The ruler is the scale that makes a block's height a quantity rather than an impression.
    #[test]
    fn the_ruler_graduates_the_deepest_lineage_in_snps() {
        let mut root = block(1, 0, 0, &[]);
        root.loci = (0..7).map(|i| locus(&format!("R{i}"), i)).collect();
        let mut deep = block(2, 1, 1, &["a"]);
        deep.loci = (0..9).map(|i| locus(&format!("D{i}"), i)).collect();
        let shallow = block(3, 1, 1, &["b"]);

        let lay = layout(&[root, deep, shallow], 1.0);
        let snps: Vec<usize> = lay.ticks.iter().map(|t| t.snps).collect();
        assert_eq!(
            snps,
            vec![5, 10, 15],
            "every {TICK_SNPS}th mutation down the 16-SNP lineage"
        );
        for w in lay.ticks.windows(2) {
            assert!(w[1].y > w[0].y, "graduations descend");
        }
        // Tick 5 falls in the root's own run, so it lands inside the root block.
        assert!(lay.ticks[0].y < lay.placed[0].rect.bottom());
        // Tick 10 is past the root's 7, so it lands in the block below.
        assert!(lay.ticks[1].y > lay.placed[0].rect.bottom());
    }

    #[test]
    fn men_do_not_overlap_their_uncles() {
        // A block that both splits and holds men has to make room for both: the men take slots
        // beside the child subtrees rather than sitting on top of them.
        let mut root = block(1, 0, 0, &["m1", "m2", "m3"]);
        root.subtree_members = 4;
        let lay = layout(&[root, block(2, 1, 1, &["a"])], 1.0);

        let child = lay.placed.iter().find(|p| p.idx == 1).unwrap().rect;
        let mut boxes: Vec<egui::Rect> = lay.members.iter().map(|m| m.rect).collect();
        boxes.sort_by(|a, b| a.left().total_cmp(&b.left()));
        for w in boxes.windows(2) {
            assert!(w[0].right() <= w[1].left(), "men must not overlap each other");
        }
        for r in &boxes {
            let man_under_root = r.top() < child.bottom() && r.bottom() > child.top();
            assert!(
                !(man_under_root && r.right() > child.left() && r.left() < child.right()),
                "a man must not sit on a sibling subtree"
            );
        }
    }

    /// The canvas has to grow to hold men that hang past the last block row.
    #[test]
    fn the_canvas_makes_room_for_the_deepest_man() {
        let lay = layout(&split(), 1.0);
        let lowest = lay.members.iter().map(|m| m.rect.bottom()).fold(0.0f32, f32::max);
        assert!(lay.size.y >= lowest, "men must not be cut off the bottom");
    }

    /// Private variants are drawn on the block scale, because they measure the same thing: the
    /// mutations between the named branch and the present.
    #[test]
    fn private_variants_extend_the_time_axis_below_a_branch() {
        let mut b = block(1, 0, 0, &["a", "b"]);
        b.members[0].private_novel = Some(6);
        b.members[0].private_publishable = Some(6);
        b.members[1].private_novel = Some(4);
        b.members[1].private_publishable = Some(4);
        let lay = layout(&[b], 1.0);

        let pv = lay.privates.first().expect("a branch whose men carry private variants");
        assert!((pv.average - 5.0).abs() < 0.001);
        assert_eq!(pv.counted, 2);
        assert!(
            (pv.rect.height() - (5.0 * ROW_H + 2.0 * PAD)).abs() < 0.5,
            "five mutations is five rows, the same scale the ruler graduates"
        );
        assert!(
            (pv.rect.top() - lay.placed[0].rect.bottom()).abs() < 0.01,
            "flush under its branch — the axis must not skip"
        );
        // The men hang below it, and it is centred on them.
        let (lo, hi) = (
            lay.members.iter().map(|m| m.rect.left()).fold(f32::MAX, f32::min),
            lay.members.iter().map(|m| m.rect.right()).fold(f32::MIN, f32::max),
        );
        assert!(
            (pv.rect.center().x - (lo + hi) / 2.0).abs() < 0.01,
            "centred on its men"
        );
        for m in &lay.members {
            assert!(m.rect.top() >= pv.rect.bottom());
        }
    }

    /// One column wide regardless of headcount. The average is a branch-level figure; sizing the box
    /// to the number of men would make an identical average look different on two branches.
    #[test]
    fn private_blocks_are_one_column_wide_whatever_the_headcount() {
        let mut one = block(1, 0, 0, &["a"]);
        one.members[0].private_novel = Some(5);
        one.members[0].private_publishable = Some(5);
        let mut many = block(1, 0, 0, &["a", "b", "c", "d"]);
        for m in &mut many.members {
            m.private_novel = Some(5);
            m.private_publishable = Some(5);
        }

        let w1 = layout(&[one], 1.0).privates[0].rect.width();
        let w4 = layout(&[many], 1.0).privates[0].rect.width();
        assert!((w1 - w4).abs() < 0.01, "same average, same box — {w1} vs {w4}");
        assert!((w1 - MEMBER_W).abs() < 0.01, "one member column wide");
    }

    /// The mean is over the men *placed on* the block, not its subtree — FTDNA reports 4 over 2
    /// participants for R-FGC29071 while 7 more men sit on branches below it.
    #[test]
    fn the_average_covers_the_men_placed_here_not_the_subtree() {
        let mut here = block(1, 0, 0, &["a", "b"]);
        here.members[0].private_novel = Some(4);
        here.members[0].private_publishable = Some(4);
        here.members[1].private_novel = Some(4);
        here.members[1].private_publishable = Some(4);
        here.subtree_members = 5;
        let mut below = block(2, 1, 1, &["c", "d", "e"]);
        for m in &mut below.members {
            m.private_novel = Some(40);
            m.private_publishable = Some(40);
        }

        let lay = layout(&[here, below], 1.0);
        let root = lay.privates.iter().find(|p| p.block == 0).unwrap();
        assert!(
            (root.average - 4.0).abs() < 0.001,
            "the subtree's 40s must not pull it up"
        );
        assert_eq!(root.counted, 2);
    }

    /// One donor at 661 novel calls would otherwise set a branch's height single-handed. The
    /// workspace already declares such a count implausible for one sample; the view honours that.
    #[test]
    fn an_implausible_donor_is_dropped_from_the_average_and_declared() {
        let mut b = block(1, 0, 0, &["ok", "ok2", "junk"]);
        b.members[0].private_novel = Some(4);
        b.members[0].private_publishable = Some(4);
        b.members[1].private_novel = Some(6);
        b.members[1].private_publishable = Some(6);
        b.members[2].private_novel = Some(661);
        b.members[2].private_publishable = Some(661);

        let lay = layout(&[b], 1.0);
        let pv = lay.privates.first().unwrap();
        assert!((pv.average - 5.0).abs() < 0.001, "the 661 must not enter the mean");
        assert_eq!(pv.counted, 2);
        assert_eq!(pv.suppressed, 1, "and the exclusion is reported, never silent");
    }

    /// The average is of the *publishable* count — the one a branch claim could rest on.
    #[test]
    fn the_average_is_of_the_gated_count() {
        let mut b = block(1, 0, 0, &["a"]);
        b.members[0].private_novel = Some(30);
        b.members[0].private_publishable = Some(3);
        let lay = layout(&[b], 1.0);
        assert!((lay.privates[0].average - 3.0).abs() < 0.001);
    }

    #[test]
    fn an_average_reads_as_a_whole_number_when_it_is_one() {
        assert_eq!(fmt_average(4.0), "4");
        assert_eq!(fmt_average(4.5), "4.5", "rounding would hide half a mutation");
    }

    /// `private_novel` is `None` until private-Y has been computed, which is not the same as zero.
    #[test]
    fn a_branch_with_no_private_y_computed_gets_no_block() {
        let lay = layout(&split(), 1.0);
        assert!(
            lay.privates.is_empty(),
            "absent evidence must not be drawn as zero mutations"
        );

        // And a mean is taken only over the men that have one.
        let mut b = block(1, 0, 0, &["a", "b", "c"]);
        b.members[0].private_novel = Some(9);
        b.members[0].private_publishable = Some(9);
        let lay = layout(&[b], 1.0);
        let pv = lay.privates.first().unwrap();
        assert!(
            (pv.average - 9.0).abs() < 0.001,
            "the two unknowns are not counted as 0"
        );
        assert_eq!(pv.counted, 1, "and the view can say the mean is over one of three");
    }

    #[test]
    fn zoom_scales_the_canvas() {
        let small = layout(&split(), 1.0);
        let big = layout(&split(), 2.0);
        assert!(big.size.x > small.size.x && big.size.y > small.size.y);
    }

    #[test]
    fn layout_of_an_empty_tree_is_empty() {
        let lay = layout(&[], 1.0);
        assert!(lay.placed.is_empty());
        assert_eq!(lay.size, egui::Vec2::ZERO);
    }
}
