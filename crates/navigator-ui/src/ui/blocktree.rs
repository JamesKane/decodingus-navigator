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
const V_GAP: f32 = 18.0; // vertical gap between depth levels — room for the connector
const PAD: f32 = 4.0;
/// A man's box: kit names are short (`B169652`, `196206`), so it is narrower than a block.
const MEMBER_W: f32 = 58.0;
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

/// The laid-out tree: block rects, member rects, and the total canvas size.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Layout {
    pub placed: Vec<Placed>,
    pub members: Vec<PlacedMember>,
    pub size: egui::Vec2,
}

/// Lines a block's box needs: the branch name, one line per equivalent SNP, and the member count.
fn lines_for(b: &Block) -> usize {
    // name + one line per SNP + the count line. Every SNP, always — see the note above.
    1 + b.loci.len() + 1
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

    // Pass 2, top-down: hand each subtree its band, centre the box in it, and stack children across.
    let mut left = vec![0.0f32; blocks.len()];
    let mut cursor = 0.0;
    for &r in &roots {
        left[r] = cursor;
        cursor += extent[r] + h_gap;
    }
    let mut placed = Vec::with_capacity(blocks.len());
    let mut members = Vec::new();
    let mut deepest = 0.0f32;
    // Pre-order, so a parent's top is always settled before its children read it.
    let mut top = vec![0.0f32; blocks.len()];
    for i in 0..blocks.len() {
        let x = left[i] + (extent[i] - box_w) / 2.0;
        let rect = egui::Rect::from_min_size(egui::pos2(x, top[i]), egui::vec2(box_w, heights[i]));
        for &c in &children[i] {
            top[c] = rect.bottom() + v_gap;
        }

        let n = slots(i);
        if n > 0 {
            let kids: f32 = children[i].iter().map(|&c| extent[c]).sum::<f32>();
            let total = kids + blocks[i].members.len() as f32 * member_w + h_gap * (n - 1) as f32;
            let mut cx = left[i] + (extent[i] - total) / 2.0;
            // Men to the left of the subclades, so a lineage that both splits and holds men reads
            // left-to-right as "these men here, then the branches below".
            let men_top = rect.bottom() + v_gap;
            for m in 0..blocks[i].members.len() {
                members.push(PlacedMember {
                    block: i,
                    member: m,
                    rect: egui::Rect::from_min_size(egui::pos2(cx, men_top), egui::vec2(member_w, member_h)),
                });
                cx += member_w + h_gap;
            }
            deepest = deepest.max(men_top + member_h);
            for &c in &children[i] {
                left[c] = cx;
                cx += extent[c] + h_gap;
            }
        }
        deepest = deepest.max(rect.bottom());
        placed.push(Placed { idx: i, rect });
    }
    let height = deepest + pad;
    Layout {
        placed,
        members,
        size: egui::vec2((cursor - h_gap).max(0.0) + pad, height),
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
            let index: HashMap<i64, usize> = drawn.iter().enumerate().map(|(i, b)| (b.node_id, i)).collect();

            for p in &lay.placed {
                let rect = p.rect.translate(origin);
                let b = &drawn[p.idx];

                // Connector to the parent, drawn even when the block itself is off-screen on one
                // side, so an edge crossing the viewport still shows.
                if let Some(pi) = b.parent.and_then(|q| index.get(&q)) {
                    let prect = lay.placed[*pi].rect.translate(origin);
                    let a = egui::pos2(prect.center().x, prect.bottom());
                    let z = egui::pos2(rect.center().x, rect.top());
                    let mid = (a.y + z.y) / 2.0;
                    if clip.intersects(egui::Rect::from_two_pos(a, z)) {
                        let stroke = egui::Stroke::new(1.0, EDGE);
                        painter.line_segment([a, egui::pos2(a.x, mid)], stroke);
                        painter.line_segment([egui::pos2(a.x, mid), egui::pos2(z.x, mid)], stroke);
                        painter.line_segment([egui::pos2(z.x, mid), z], stroke);
                    }
                }

                // Cull: everything below is per-block text layout, the expensive part.
                if !clip.intersects(rect) {
                    continue;
                }

                let (bg, stroke) = match (b.candidate, b.members.is_empty()) {
                    (true, _) => (BLOCK_BG_CANDIDATE, egui::Stroke::new(1.5, CANDIDATE_STROKE)),
                    (false, true) => (BLOCK_BG, egui::Stroke::new(1.0, BLOCK_STROKE)),
                    (false, false) => (BLOCK_BG_PLACED, egui::Stroke::new(1.0, BLOCK_STROKE)),
                };
                painter.rect_filled(rect, 3.0, bg);
                let stroke = if selected == Some(b.node_id) {
                    egui::Stroke::new(2.0, SELECTED_STROKE)
                } else {
                    stroke
                };
                painter.rect_stroke(rect, 3.0, stroke);

                let pad = PAD * zoom;
                let row = ROW_H * zoom;
                let mut y = rect.top() + pad;
                let left = rect.left() + pad;
                // Clipped to the box: a label that outgrows its block is then cropped rather than
                // spilling across the canvas, whatever the text turns out to be.
                let inner = painter.with_clip_rect(rect.shrink(1.0));
                let put = |text: String, color: egui::Color32, f: &egui::FontId, y: f32| {
                    inner.text(egui::pos2(left, y), egui::Align2::LEFT_TOP, text, f.clone(), color);
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
                // The men *here* are boxes below the block now, so the line carries only what the
                // canvas can't show at a glance: how many sit anywhere in the subtree beneath.
                // ASCII only: `▾` is absent from the bundled font and renders as a tofu box.
                let sub = match b.subtree_members {
                    below if below > b.members.len() => format!("+{below}"),
                    _ => String::new(),
                };
                if !sub.is_empty() {
                    put(sub, egui::Color32::from_gray(140), &small, y);
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

            // Men, as the Big Tree draws them: a grey box on a stem below the block they sit on.
            // They are the evidence the phylogeny rests on, so they belong in the diagram — the
            // roster beside it stays for the private-variant counts and for scanning a long list.
            for pm in &lay.members {
                let rect = pm.rect.translate(origin);
                let b = &drawn[pm.block];
                let stem = lay.placed[pm.block].rect.translate(origin);
                let a = egui::pos2(stem.center().x, stem.bottom());
                let z = egui::pos2(rect.center().x, rect.top());
                if clip.intersects(egui::Rect::from_two_pos(a, z)) {
                    let stroke = egui::Stroke::new(1.0, EDGE);
                    let mid = (a.y + z.y) / 2.0;
                    painter.line_segment([a, egui::pos2(a.x, mid)], stroke);
                    painter.line_segment([egui::pos2(a.x, mid), egui::pos2(z.x, mid)], stroke);
                    painter.line_segment([egui::pos2(z.x, mid), z], stroke);
                }
                if !clip.intersects(rect) {
                    continue;
                }
                let m = &b.members[pm.member];
                painter.rect_filled(rect, 2.0, MEMBER_BG);
                painter.rect_stroke(rect, 2.0, egui::Stroke::new(1.0, BLOCK_STROKE));
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
        assert!(lay.placed[1].rect.top() > lay.placed[0].rect.bottom());
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

    #[test]
    fn layout_centres_a_parent_over_its_children() {
        let lay = layout(&split(), 1.0);
        let (root, a, b) = (
            lay.placed[0].rect.center().x,
            lay.placed[1].rect.center().x,
            lay.placed[2].rect.center().x,
        );
        assert!(a < b, "siblings run left to right in pre-order");
        assert!((root - (a + b) / 2.0).abs() < 0.5, "root sits over its children");
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
        // name + 600 SNPs + count.
        assert!((lay.placed[0].rect.height() - (602.0 * ROW_H + 2.0 * PAD)).abs() < 0.5);
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
    fn men_hang_below_their_own_block() {
        let lay = layout(&split(), 1.0);
        assert_eq!(lay.members.len(), 2, "one box per man");
        for pm in &lay.members {
            let block = lay.placed.iter().find(|p| p.idx == pm.block).unwrap();
            assert!(pm.rect.top() > block.rect.bottom(), "a man hangs below his block");
        }
        // Each man is under his own terminal, not pooled under the root.
        assert_ne!(lay.members[0].block, lay.members[1].block);
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
