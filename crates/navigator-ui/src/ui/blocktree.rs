//! The project **block tree** view (`impl NavigatorApp`) — the cohort haplotree for the open
//! project, drawn as a canvas of blocks: depth → x, tidy-tree row order → y, members listed inside
//! their own terminal block.
//!
//! The aggregate is built off the UI thread (`App::project_block_tree`, see
//! `documents/design/project-block-tree.md`); this module only lays it out and paints it.
//!
//! Two performance rules, because a group project can hold thousands of members:
//!
//! - **Layout is computed once per (tree, expansion, zoom) and memoized**, not rebuilt per frame.
//!   [`layout`] is a pure function over `&[Block]`, so it is testable without a canvas.
//! - **Drawing is culled to `clip_rect`.** Only blocks actually on screen are painted, so a tree
//!   with thousands of blocks costs the same per frame as one with a dozen.

use std::collections::{HashMap, HashSet};

use navigator_app::Block;

use super::*;

// Canvas geometry, at zoom 1.0.
const COL_W: f32 = 190.0; // horizontal distance between depth levels
const BOX_W: f32 = 168.0;
const ROW_H: f32 = 20.0; // one line of text inside a block
const V_GAP: f32 = 8.0; // vertical gap between sibling subtrees
const PAD: f32 = 6.0;
/// Members listed inside a collapsed block before it shows a "+N more" line. Expanding shows all.
const MEMBERS_PREVIEW: usize = 3;

const BLOCK_BG: egui::Color32 = egui::Color32::from_rgb(46, 52, 60);
const BLOCK_BG_PLACED: egui::Color32 = egui::Color32::from_rgb(50, 74, 58); // carries members
/// A candidate branch reads as *provisional*: amber, not the green of a published branch. It is an
/// inference from shared private variants, and must never be mistaken for a named haplogroup.
const BLOCK_BG_CANDIDATE: egui::Color32 = egui::Color32::from_rgb(74, 62, 38);
const CANDIDATE_STROKE: egui::Color32 = egui::Color32::from_rgb(200, 155, 70);
const BLOCK_STROKE: egui::Color32 = egui::Color32::from_rgb(90, 100, 112);
const EDGE: egui::Color32 = egui::Color32::from_rgb(80, 88, 98);
const MEMBER_FG: egui::Color32 = egui::Color32::from_rgb(150, 200, 230);

/// One laid-out block: where it sits, and how many member lines it shows.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Placed {
    /// Index into the `blocks` slice that was laid out.
    pub idx: usize,
    pub rect: egui::Rect,
    /// Member lines rendered inside the box (may be fewer than the block's members).
    pub shown_members: usize,
}

/// The laid-out tree: block rects plus the total canvas size.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Layout {
    pub placed: Vec<Placed>,
    pub size: egui::Vec2,
}

/// How many text lines a block's box needs: the name, the SNP line, and its member lines.
fn lines_for(b: &Block, expanded: bool) -> (usize, usize) {
    let shown = if expanded {
        b.members.len()
    } else {
        b.members.len().min(MEMBERS_PREVIEW)
    };
    // name + SNP/branch summary line + member lines + a "+N more" line when truncated.
    let more = usize::from(shown < b.members.len());
    (2 + shown + more, shown)
}

/// Lay `blocks` (in pre-order, as [`ProjectBlockTree`] delivers them) onto a canvas: **depth → x**,
/// tidy-tree order → **y**. A parent is centred on the vertical extent of its children, so a branch
/// point sits between the lineages it splits into rather than at the top of them.
///
/// Pure: no `Ui`, no state. Two passes — extents bottom-up, then positions top-down.
pub(crate) fn layout(blocks: &[Block], expanded: &HashSet<i64>, zoom: f32) -> Layout {
    if blocks.is_empty() {
        return Layout::default();
    }
    let zoom = zoom.clamp(0.4, 2.5);
    let (col_w, box_w, row_h, v_gap, pad) = (COL_W * zoom, BOX_W * zoom, ROW_H * zoom, V_GAP * zoom, PAD * zoom);

    let index: HashMap<i64, usize> = blocks.iter().enumerate().map(|(i, b)| (b.node_id, i)).collect();
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); blocks.len()];
    let mut roots: Vec<usize> = Vec::new();
    for (i, b) in blocks.iter().enumerate() {
        match b.parent.and_then(|p| index.get(&p)) {
            Some(&p) => children[p].push(i),
            None => roots.push(i),
        }
    }

    let mut heights = Vec::with_capacity(blocks.len());
    let mut shown = Vec::with_capacity(blocks.len());
    for b in blocks {
        let (n, s) = lines_for(b, expanded.contains(&b.node_id));
        heights.push(n as f32 * row_h + 2.0 * pad);
        shown.push(s);
    }

    // Pass 1, bottom-up: the vertical extent each subtree needs. `blocks` is pre-order, so iterating
    // in reverse visits every child before its parent.
    let mut extent = heights.clone();
    for i in (0..blocks.len()).rev() {
        if children[i].is_empty() {
            continue;
        }
        let kids: f32 = children[i].iter().map(|&c| extent[c]).sum::<f32>() + v_gap * (children[i].len() - 1) as f32;
        extent[i] = extent[i].max(kids);
    }

    // Pass 2, top-down: hand each subtree its band, centre the block within it, and stack children.
    let mut top = vec![0.0f32; blocks.len()];
    let mut cursor = 0.0;
    for &r in &roots {
        top[r] = cursor;
        cursor += extent[r] + v_gap;
    }
    let mut placed = Vec::with_capacity(blocks.len());
    let mut max_x: f32 = 0.0;
    for i in 0..blocks.len() {
        let kids = &children[i];
        if !kids.is_empty() {
            let total: f32 = kids.iter().map(|&c| extent[c]).sum::<f32>() + v_gap * (kids.len() - 1) as f32;
            let mut cy = top[i] + (extent[i] - total) / 2.0;
            for &c in kids {
                top[c] = cy;
                cy += extent[c] + v_gap;
            }
        }
        let x = blocks[i].depth as f32 * col_w;
        let y = top[i] + (extent[i] - heights[i]) / 2.0;
        let rect = egui::Rect::from_min_size(egui::pos2(x, y), egui::vec2(box_w, heights[i]));
        max_x = max_x.max(rect.right());
        placed.push(Placed {
            idx: i,
            rect,
            shown_members: shown[i],
        });
    }
    Layout {
        placed,
        size: egui::vec2(max_x + pad, (cursor - v_gap).max(0.0)),
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
            s
        });
        // Every label is resolved before the closure: `self.tr` borrows `self`, and the zoom slider
        // needs `&mut` — the two can't coexist inside one closure.
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
            return;
        }
        ui.add_space(4.0);

        let lay = layout(&tree.blocks, &self.blocktree_expanded, zoom);
        let mut toggle: Option<i64> = None;
        let mut open_subject: Option<SampleGuid> = None;

        egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
            let (canvas, _) = ui.allocate_exact_size(lay.size, egui::Sense::hover());
            let painter = ui.painter_at(canvas);
            let clip = ui.clip_rect();
            let origin = canvas.min.to_vec2();
            let zoom = zoom.clamp(0.4, 2.5);
            let font = egui::FontId::proportional(11.0 * zoom);
            let small = egui::FontId::proportional(9.5 * zoom);
            let index: HashMap<i64, usize> = tree.blocks.iter().enumerate().map(|(i, b)| (b.node_id, i)).collect();

            for p in &lay.placed {
                let rect = p.rect.translate(origin);
                let b = &tree.blocks[p.idx];

                // Connector to the parent, drawn even when the block itself is off-screen on one
                // side, so an edge crossing the viewport still shows.
                if let Some(pi) = b.parent.and_then(|q| index.get(&q)) {
                    let prect = lay.placed[*pi].rect.translate(origin);
                    let a = egui::pos2(prect.right(), prect.center().y);
                    let z = egui::pos2(rect.left(), rect.center().y);
                    let mid = (a.x + z.x) / 2.0;
                    if clip.intersects(egui::Rect::from_two_pos(a, z)) {
                        let stroke = egui::Stroke::new(1.0, EDGE);
                        painter.line_segment([a, egui::pos2(mid, a.y)], stroke);
                        painter.line_segment([egui::pos2(mid, a.y), egui::pos2(mid, z.y)], stroke);
                        painter.line_segment([egui::pos2(mid, z.y), z], stroke);
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
                painter.rect_stroke(rect, 3.0, stroke);

                let pad = PAD * zoom;
                let row = ROW_H * zoom;
                let mut y = rect.top() + pad;
                let left = rect.left() + pad;
                let put = |text: String, color: egui::Color32, f: &egui::FontId, y: f32| {
                    painter.text(egui::pos2(left, y), egui::Align2::LEFT_TOP, text, f.clone(), color);
                };

                // A candidate has no published name — the view supplies the label, localized.
                let (title, title_fg) = if b.candidate {
                    (candidate_label.clone(), CANDIDATE_STROKE)
                } else {
                    (b.name.clone(), egui::Color32::from_gray(230))
                };
                put(title, title_fg, &font, y);
                y += row;

                // The block's own weight: equivalent SNPs, and how many branches it folded away.
                let mut sub = format!("{} SNP", b.loci.len());
                if b.loci.len() != 1 {
                    sub.push('s');
                }
                if !b.collapsed.is_empty() {
                    sub.push_str(&format!(" · +{} branches", b.collapsed.len()));
                }
                if b.subtree_members > b.members.len() {
                    sub.push_str(&format!(" · {} below", b.subtree_members));
                }
                put(sub, egui::Color32::from_gray(150), &small, y);
                y += row;

                for m in b.members.iter().take(p.shown_members) {
                    // A member's own unshared private count, when private-Y has been computed for
                    // them. `None` means never analyzed, which is not the same as zero — so it shows
                    // nothing rather than "0".
                    let label = match m.private_novel {
                        Some(n) if n > 0 => format!("{}  ({n})", m.name),
                        _ => m.name.clone(),
                    };
                    put(label, MEMBER_FG, &small, y);
                    y += row;
                }
                if p.shown_members < b.members.len() {
                    put(
                        format!("+{} more", b.members.len() - p.shown_members),
                        egui::Color32::from_gray(140),
                        &small,
                        y,
                    );
                }

                // Click a block to expand it (full member list + the SNP names on hover).
                let resp = ui.interact(rect, egui::Id::new(("blocktree", b.node_id)), egui::Sense::click());
                if resp.clicked() {
                    toggle = Some(b.node_id);
                }
                if resp.hovered() && !b.loci.is_empty() {
                    let names: Vec<&str> = b.loci.iter().map(|l| l.name.as_str()).take(40).collect();
                    let extra = b.loci.len().saturating_sub(names.len());
                    let mut tip = names.join(", ");
                    if extra > 0 {
                        tip.push_str(&format!(" … +{extra}"));
                    }
                    if !b.collapsed.is_empty() {
                        tip.push_str(&format!("\n\nfolded: {}", b.collapsed.join(" → ")));
                    }
                    resp.on_hover_text(tip);
                }
                // Double-click a placed block to jump to its first member's subject page.
                if let Some(m) = b.members.first() {
                    if ui
                        .interact(rect, egui::Id::new(("blocktree-open", b.node_id)), egui::Sense::click())
                        .double_clicked()
                    {
                        open_subject = Some(m.guid);
                    }
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
        if let Some(id) = toggle {
            if !self.blocktree_expanded.remove(&id) {
                self.blocktree_expanded.insert(id);
            }
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

    fn block(id: i64, parent: i64, depth: usize, members: &[&str]) -> Block {
        Block {
            node_id: id,
            name: format!("N{id}"),
            parent: (parent != 0).then_some(parent),
            depth,
            loci: Vec::new(),
            members: members.iter().map(|m| member(m)).collect(),
            subtree_members: members.len(),
            collapsed: Vec::new(),
            candidate: false,
        }
    }

    /// `R ─┬─> R1(a)`
    ///     `└─> R2(b)`
    fn split() -> Vec<Block> {
        vec![block(1, 0, 0, &[]), block(2, 1, 1, &["a"]), block(3, 1, 1, &["b"])]
    }

    #[test]
    fn layout_places_depth_on_x() {
        let lay = layout(&split(), &HashSet::new(), 1.0);
        assert_eq!(lay.placed.len(), 3);
        assert_eq!(lay.placed[0].rect.left(), 0.0, "root sits in the first column");
        assert_eq!(lay.placed[1].rect.left(), COL_W);
        assert_eq!(lay.placed[2].rect.left(), COL_W);
    }

    #[test]
    fn layout_centres_a_parent_between_its_children() {
        let lay = layout(&split(), &HashSet::new(), 1.0);
        let (root, a, b) = (
            lay.placed[0].rect.center().y,
            lay.placed[1].rect.center().y,
            lay.placed[2].rect.center().y,
        );
        assert!(a < b, "siblings stack in pre-order");
        assert!(
            (root - (a + b) / 2.0).abs() < 0.5,
            "root should sit between its children"
        );
    }

    #[test]
    fn layout_siblings_do_not_overlap() {
        let lay = layout(&split(), &HashSet::new(), 1.0);
        assert!(
            lay.placed[1].rect.bottom() <= lay.placed[2].rect.top(),
            "sibling boxes must not overlap"
        );
    }

    #[test]
    fn expanding_a_block_shows_every_member() {
        let many: Vec<&str> = vec!["a", "b", "c", "d", "e"];
        let blocks = vec![block(1, 0, 0, &many)];
        let collapsed = layout(&blocks, &HashSet::new(), 1.0);
        assert_eq!(collapsed.placed[0].shown_members, MEMBERS_PREVIEW);

        let expanded: HashSet<i64> = [1].into_iter().collect();
        let open = layout(&blocks, &expanded, 1.0);
        assert_eq!(open.placed[0].shown_members, many.len());
        assert!(
            open.placed[0].rect.height() > collapsed.placed[0].rect.height(),
            "the box must grow to fit the extra member lines"
        );
    }

    #[test]
    fn zoom_scales_the_canvas() {
        let small = layout(&split(), &HashSet::new(), 1.0);
        let big = layout(&split(), &HashSet::new(), 2.0);
        assert!(big.size.x > small.size.x && big.size.y > small.size.y);
    }

    #[test]
    fn layout_of_an_empty_tree_is_empty() {
        let lay = layout(&[], &HashSet::new(), 1.0);
        assert!(lay.placed.is_empty());
        assert_eq!(lay.size, egui::Vec2::ZERO);
    }
}
