//! Pure ancestry / genome-region visualization helpers, extracted from the UI shell. Leaf drawing
//! functions over `egui` with no `App`/`self` state — easy to read, test, and reuse independently of
//! the view code that calls them.

use eframe::egui;
use navigator_app::{AncestryResult, AncestrySegment, AssetStatus, GenomeRegions, IbdSegment, SuperPopulationSummary};
use navigator_domain::ancestry::{population_color, population_name, population_super};

/// Sort key for chromosome names: autosomes 1–22, then X, Y, M, then anything else.
fn chrom_sort_key(chr: &str) -> (u8, i64) {
    let bare = chr.trim_start_matches("chr").to_ascii_uppercase();
    if let Ok(n) = bare.parse::<i64>() {
        (0, n)
    } else {
        match bare.as_str() {
            "X" => (1, 0),
            "Y" => (2, 0),
            "M" | "MT" => (3, 0),
            _ => (4, 0),
        }
    }
}

/// Draw a per-chromosome **IBD-segment ideogram** (gap §8): one horizontal bar per chromosome that
/// carries a shared segment, scaled to the chromosome's true length when `regions` is available (else
/// to the segments' own span), each IBD segment painted as a teal block (brighter = longer in cM) with
/// per-segment hover details. Mirrors [`draw_chromosome_painting`]'s painter approach.
pub(crate) fn draw_ibd_segments(ui: &mut egui::Ui, segments: &[IbdSegment], regions: Option<&GenomeRegions>) {
    use std::collections::BTreeMap;
    let mut by_chr: BTreeMap<String, Vec<&IbdSegment>> = BTreeMap::new();
    for s in segments {
        by_chr.entry(s.chromosome.clone()).or_default().push(s);
    }
    if by_chr.is_empty() {
        ui.label(egui::RichText::new("No shared IBD segments to paint.").weak());
        return;
    }
    let mut chroms: Vec<String> = by_chr.keys().cloned().collect();
    chroms.sort_by_key(|c| chrom_sort_key(c));

    let (label_w, bar_w, bar_h, gap) = (48.0f32, 320.0f32, 12.0f32, 4.0f32);
    for chr in &chroms {
        let segs = &by_chr[chr];
        // Scale to the true chromosome length when known, else to the max segment end.
        let chr_len = regions
            .and_then(|r| r.chromosomes.get(chr))
            .map(|c| c.length)
            .filter(|&l| l > 0)
            .unwrap_or_else(|| segs.iter().map(|s| s.end_position).max().unwrap_or(1))
            .max(1) as f32;
        ui.horizontal(|ui| {
            ui.allocate_ui(egui::vec2(label_w, bar_h), |ui| {
                ui.label(egui::RichText::new(chr).small())
            });
            let (rect, resp) = ui.allocate_exact_size(egui::vec2(bar_w, bar_h), egui::Sense::hover());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 2.0, egui::Color32::from_gray(30));
            let mut hover: Option<String> = None;
            let hover_x = resp.hover_pos().map(|p| p.x);
            for &s in segs.iter() {
                let x0 = rect.left() + (s.start_position.max(0) as f32 / chr_len) * rect.width();
                let x1 = rect.left() + (s.end_position.max(0) as f32 / chr_len) * rect.width();
                let seg =
                    egui::Rect::from_min_max(egui::pos2(x0, rect.top()), egui::pos2(x1.max(x0 + 1.5), rect.bottom()));
                let t = (s.length_cm / 30.0).clamp(0.3, 1.0) as f32;
                let col = egui::Color32::from_rgb(30, (90.0 + 120.0 * t) as u8, (90.0 + 100.0 * t) as u8);
                painter.rect_filled(seg, 0.0, col);
                if let Some(hx) = hover_x {
                    if hx >= seg.left() && hx <= seg.right() {
                        hover = Some(format!(
                            "{}:{}–{} · {:.1} cM{}",
                            s.chromosome,
                            s.start_position,
                            s.end_position,
                            s.length_cm,
                            s.snp_count.map(|n| format!(" · {n} SNPs")).unwrap_or_default()
                        ));
                    }
                }
            }
            if let Some(t) = hover {
                resp.on_hover_text(t);
            }
        });
        ui.add_space(gap);
    }
}

/// Draw the per-chromosome local-ancestry painting: one horizontal bar per autosome (each
/// normalized to full width), segments colored by ancestry, plus a legend of the ancestries shown.
pub(crate) fn draw_chromosome_painting(ui: &mut egui::Ui, segments: &[AncestrySegment]) {
    use std::collections::BTreeMap;
    // Group by autosome number → the two copies' segments. Non-autosomes (X/Y/M / the chr99 fallback)
    // are skipped — this is autosomal local ancestry.
    let mut by_chr: BTreeMap<i64, [Vec<&AncestrySegment>; 2]> = BTreeMap::new();
    for s in segments {
        let Ok(n) = s.contig.trim_start_matches("chr").parse::<i64>() else {
            continue;
        };
        if !(1..=22).contains(&n) {
            continue;
        }
        by_chr.entry(n).or_default()[(s.copy as usize).min(1)].push(s);
    }
    let label_w = 42.0;
    let bar_w = 300.0;
    let copy_h = 7.0; // each of the two copy tracks
    let gap = 2.0;
    for (n, copies) in by_chr {
        // Shared bp span across both copies so the two tracks align.
        let lo = copies.iter().flatten().map(|s| s.start).min().unwrap_or(1);
        let hi = copies
            .iter()
            .flatten()
            .map(|s| s.end)
            .max()
            .unwrap_or(lo + 1)
            .max(lo + 1);
        let span = (hi - lo).max(1) as f32;
        ui.horizontal(|ui| {
            ui.allocate_ui(egui::vec2(label_w, copy_h * 2.0 + gap), |ui| {
                ui.label(format!("chr{n}"))
            });
            let (rect, _) = ui.allocate_exact_size(egui::vec2(bar_w, copy_h * 2.0 + gap), egui::Sense::hover());
            let painter = ui.painter_at(rect);
            for (c, segs) in copies.iter().enumerate() {
                let top = rect.top() + c as f32 * (copy_h + gap);
                let track = egui::Rect::from_min_size(egui::pos2(rect.left(), top), egui::vec2(rect.width(), copy_h));
                painter.rect_filled(track, 2.0, egui::Color32::from_gray(30));
                for s in segs {
                    let x0 = track.left() + (s.start - lo) as f32 / span * track.width();
                    let x1 = track.left() + (s.end - lo) as f32 / span * track.width();
                    let seg_rect = egui::Rect::from_min_max(
                        egui::pos2(x0, track.top()),
                        egui::pos2(x1.max(x0 + 1.0), track.bottom()),
                    );
                    painter.rect_filled(seg_rect, 0.0, parse_hex_color(&population_color(&s.population_code)));
                }
            }
        });
    }
    // Legend: distinct ancestries present.
    let mut seen: Vec<&str> = Vec::new();
    for s in segments {
        if !seen.contains(&s.population_code.as_str()) {
            seen.push(&s.population_code);
        }
    }
    ui.add_space(2.0);
    ui.horizontal_wrapped(|ui| {
        for code in seen {
            let (r, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
            ui.painter()
                .circle_filled(r.center(), 4.0, parse_hex_color(&population_color(code)));
            ui.label(egui::RichText::new(population_name(code)).small());
            ui.add_space(6.0);
        }
    });
}

/// Points along a circle arc from angle `a0` to `a1` (radians), `steps`+1 samples.
fn arc_points(c: egui::Pos2, r: f32, a0: f32, a1: f32, steps: usize) -> Vec<egui::Pos2> {
    (0..=steps)
        .map(|i| {
            let t = a0 + (a1 - a0) * (i as f32 / steps as f32);
            egui::pos2(c.x + r * t.cos(), c.y + r * t.sin())
        })
        .collect()
}

/// Draw a solid **pie** chart (`size`×`size`) from `(percentage, color)` slices. Each slice is a fan
/// from the centre to the outer arc — a simple convex polygon egui tessellates cleanly. (The previous
/// *donut* built an annular sector as one concave path, whose tessellation left a stray wedge in the
/// middle.) Slices below 0.5 % are skipped; the first slice starts at 12 o'clock.
pub(crate) fn draw_pie(ui: &mut egui::Ui, size: f32, slices: &[(f64, egui::Color32)]) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(size, size), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let c = rect.center();
    let r = size * 0.47;
    let total: f32 = slices.iter().map(|(p, _)| *p as f32).sum::<f32>().max(1.0);
    let mut a0 = -std::f32::consts::FRAC_PI_2;
    for (pct, color) in slices {
        if *pct < 0.5 {
            continue;
        }
        let a1 = a0 + (*pct as f32 / total) * std::f32::consts::TAU;
        let mut pts = vec![c];
        pts.extend(arc_points(c, r, a0, a1, 40));
        painter.add(egui::epaint::PathShape {
            points: pts,
            closed: true,
            fill: *color,
            stroke: egui::epaint::PathStroke::NONE,
        });
        a0 = a1;
    }
}

/// Pie chart of the super-population proportions (one slice per super-population, colored by
/// continent).
pub(crate) fn draw_ancestry_donut(ui: &mut egui::Ui, summary: &[SuperPopulationSummary]) {
    let slices: Vec<(f64, egui::Color32)> = summary
        .iter()
        .map(|s| {
            let code = s.populations.first().and_then(|c| population_super(c)).unwrap_or("");
            (s.percentage, parse_hex_color(&population_color(code)))
        })
        .collect();
    draw_pie(ui, 120.0, &slices);
}

/// A generic donut from pre-colored `(percentage, color)` slices (used by the Simple-mode brief for
/// the ancient-ancestry pie, whose components carry their own palette colors). Optionally labels the
/// hole with the largest slice's share.
/// Pie chart from explicit `(percentage, color)` slices (the ancient-component report, which carries
/// its own colors). `_center_pct` is kept for call-site compatibility but no longer rendered — a
/// solid pie has no centre to label.
pub(crate) fn draw_color_donut(ui: &mut egui::Ui, slices: &[(f64, egui::Color32)], _center_pct: Option<f64>) {
    draw_pie(ui, 120.0, slices);
}

/// Draw a detailed ancestry breakdown (the fine-population or ancient-component report): the
/// estimate's `components`, sorted by share, as a name/percentage grid with a proportion bar, plus a
/// provenance line (method + SNP count). `id_salt` keeps each report's grid distinct.
pub(crate) fn draw_population_components(ui: &mut egui::Ui, result: &AncestryResult, _id_salt: &str, top_n: usize) {
    // (friendly name, code, percentage) — the component carries `population_name` (e.g. EEF → "Early
    // European Farmer"), so the legend reads in plain language rather than codes.
    let mut comps: Vec<(&str, &str, f64)> = result
        .components
        .iter()
        .filter(|c| c.percentage >= 0.05)
        .map(|c| (c.population_name.as_str(), c.population_code.as_str(), c.percentage))
        .collect();
    comps.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
    if comps.is_empty() {
        ui.label(egui::RichText::new("No components above 0.05%.").weak());
        return;
    }
    let shown = &comps[..comps.len().min(top_n)];
    let slices: Vec<(f64, egui::Color32)> = shown
        .iter()
        .map(|(_, code, pct)| (*pct, parse_hex_color(&population_color(code))))
        .collect();

    // Horizontal: pie on the left, a colour-swatch legend (friendly name + %) on the right — compact
    // vertically, so the modern + ancient panels can sit side by side in the Advanced view.
    ui.horizontal_top(|ui| {
        draw_pie(ui, 108.0, &slices);
        ui.add_space(12.0);
        ui.vertical(|ui| {
            for (name, code, pct) in shown {
                ui.horizontal(|ui| {
                    let (sw, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                    ui.painter().rect_filled(sw, 2.0, parse_hex_color(&population_color(code)));
                    ui.add_space(2.0);
                    ui.label(egui::RichText::new(format!("{pct:.1}%")).strong());
                    ui.label(*name);
                });
            }
            ui.add_space(2.0);
            ui.label(
                egui::RichText::new(format!(
                    "{} · {}/{} SNPs · confidence {:.0}%",
                    result.method,
                    result.snps_with_genotype,
                    result.snps_analyzed,
                    result.confidence_level * 100.0
                ))
                .weak()
                .small(),
            );
        });
    });
}

/// Draw the super-population composition as a single stacked horizontal bar (segment widths =
/// proportions, colored by continent).
pub(crate) fn draw_composition_bar(ui: &mut egui::Ui, summary: &[SuperPopulationSummary]) {
    let w = ui.available_width().min(360.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(w, 16.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 3.0, egui::Color32::from_gray(30));
    let mut x = rect.left();
    for s in summary {
        let seg_w = rect.width() * (s.percentage as f32 / 100.0).clamp(0.0, 1.0);
        let seg = egui::Rect::from_min_size(egui::pos2(x, rect.top()), egui::vec2(seg_w, rect.height()));
        let code = s.populations.first().and_then(|c| population_super(c)).unwrap_or("");
        painter.rect_filled(seg, 0.0, parse_hex_color(&population_color(code)));
        x += seg_w;
    }
}

/// Parse a `#RRGGBB` hex color, falling back to grey on a malformed string.
pub(crate) fn parse_hex_color(hex: &str) -> egui::Color32 {
    let h = hex.trim_start_matches('#');
    if h.len() == 6 {
        if let (Ok(r), Ok(g), Ok(b)) = (
            u8::from_str_radix(&h[0..2], 16),
            u8::from_str_radix(&h[2..4], 16),
            u8::from_str_radix(&h[4..6], 16),
        ) {
            return egui::Color32::from_rgb(r, g, b);
        }
    }
    egui::Color32::from_gray(128)
}

/// A compact "data sources" line: which ancestry/IBD reference assets are present and
/// integrity-verified (✓ verified · • present-but-unverified · ✗ absent).
pub(crate) fn asset_status_line(ui: &mut egui::Ui, assets: &[AssetStatus]) {
    if assets.is_empty() {
        return;
    }
    ui.horizontal_wrapped(|ui| {
        ui.label(egui::RichText::new("Data sources:").small().weak());
        for a in assets {
            let (mark, col, hover) = if a.verified {
                ("✓", egui::Color32::from_rgb(80, 170, 90), "present, integrity-verified")
            } else if a.present {
                ("•", egui::Color32::from_gray(150), "present (no manifest to verify)")
            } else {
                ("✗", egui::Color32::from_rgb(170, 90, 90), "not installed")
            };
            ui.colored_label(col, egui::RichText::new(format!("{} {mark}", a.name)).small())
                .on_hover_text(hover);
            ui.add_space(4.0);
        }
    });
}

/// A single variant tick on a [`draw_variant_track`] chromosome bar.
pub(crate) struct VariantMark {
    pub name: String,
    pub position: i64,
    pub color: egui::Color32,
    /// Human-readable state for the hover tooltip (e.g. "in-tree derived", "novel").
    pub state: &'static str,
}

/// A shaded background region on a variant track (chrY PAR/heterochromatin, chrM HVR/coding).
pub(crate) struct TrackRegion {
    pub start: i64,
    pub end: i64,
    pub color: egui::Color32,
    pub label: String,
}

/// Draw a single-chromosome **variant track**: one horizontal bar scaled to `length`, optional
/// shaded background regions, and a vertical tick per variant colored by its consensus state. Hover
/// over the bar surfaces the nearest variant (`name · pos · state`) and any region under the cursor.
/// Replaces the genome-wide karyotype ideogram for the Y/mt variant views. Mirrors the
/// [`draw_ibd_segments`] painter approach (no `egui_plot`).
pub(crate) fn draw_variant_track(
    ui: &mut egui::Ui,
    chrom_label: &str,
    length: i64,
    regions: &[TrackRegion],
    variants: &[VariantMark],
) {
    if length <= 0 {
        ui.label(egui::RichText::new("Chromosome length unavailable.").weak());
        return;
    }
    let len = length as f32;
    let bar_h = 22.0f32;
    let full_w = ui.available_width().min(720.0);

    let (rect, resp) = ui.allocate_exact_size(egui::vec2(full_w, bar_h), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 3.0, egui::Color32::from_gray(28));

    // Background region shading.
    for r in regions {
        let x0 = rect.left() + (r.start.max(0) as f32 / len).clamp(0.0, 1.0) * rect.width();
        let x1 = rect.left() + (r.end.max(0) as f32 / len).clamp(0.0, 1.0) * rect.width();
        let seg = egui::Rect::from_min_max(egui::pos2(x0, rect.top()), egui::pos2(x1.max(x0 + 1.0), rect.bottom()));
        painter.rect_filled(seg, 0.0, r.color);
    }
    painter.rect_stroke(rect, 3.0, egui::Stroke::new(1.0_f32, egui::Color32::from_gray(70)));

    // Variant ticks.
    let hover_x = resp.hover_pos().map(|p| p.x);
    let mut nearest: Option<(f32, String)> = None;
    for v in variants {
        let x = rect.left() + (v.position.max(0) as f32 / len).clamp(0.0, 1.0) * rect.width();
        painter.line_segment(
            [egui::pos2(x, rect.top() + 1.0), egui::pos2(x, rect.bottom() - 1.0)],
            egui::Stroke::new(1.5_f32, v.color),
        );
        if let Some(hx) = hover_x {
            let d = (hx - x).abs();
            if d <= 6.0 && nearest.as_ref().map_or(true, |(bd, _)| d < *bd) {
                nearest = Some((d, format!("{} · {} · {}", v.name, v.position, v.state)));
            }
        }
    }

    if let Some((_, text)) = nearest {
        resp.on_hover_text(text);
    } else if let (Some(hx), false) = (hover_x, regions.is_empty()) {
        // No tick nearby → report the region under the cursor, if any.
        let g = (((hx - rect.left()) / rect.width()).clamp(0.0, 1.0) * len) as i64;
        if let Some(r) = regions.iter().find(|r| r.start <= g && g < r.end) {
            resp.on_hover_text(format!("{} · {:.0} bp", r.label, g as f64));
        }
    }

    ui.add_space(2.0);
    ui.label(
        egui::RichText::new(format!(
            "{chrom_label} · {} variants · {:.0} kb",
            variants.len(),
            len / 1000.0
        ))
        .small()
        .weak(),
    );

    // Region legend.
    if !regions.is_empty() {
        let mut seen: Vec<&str> = Vec::new();
        ui.horizontal_wrapped(|ui| {
            for r in regions {
                if seen.contains(&r.label.as_str()) {
                    continue;
                }
                seen.push(&r.label);
                let (rr, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                ui.painter().rect_filled(rr, 2.0, r.color);
                ui.label(egui::RichText::new(r.label.as_str()).small().weak());
                ui.add_space(6.0);
            }
        });
    }
}

/// PC1×PC2 ancestry scatter (gap §8): reference population centroids (colored by population) as the
/// backdrop, the donor's projected coordinate marked as a white diamond. `sample` is the donor's
/// (PC1, PC2); `reference` is `(population_code, pc1, pc2)`. Hover a point for its population name.
pub(crate) fn draw_pca_scatter(ui: &mut egui::Ui, sample: Option<(f64, f64)>, reference: &[(String, f64, f64)]) {
    use egui_plot::{MarkerShape, Plot, PlotPoints, Points};
    if sample.is_none() && reference.is_empty() {
        ui.label(egui::RichText::new("PCA coordinates not available for this estimate.").weak());
        return;
    }
    Plot::new("pca_scatter")
        .height(320.0)
        .data_aspect(1.0)
        .x_axis_label("PC1")
        .y_axis_label("PC2")
        .show(ui, |plot_ui| {
            for (code, x, y) in reference {
                let color = parse_hex_color(&population_color(code));
                plot_ui.points(
                    Points::new(PlotPoints::new(vec![[*x, *y]]))
                        .radius(4.0_f32)
                        .color(color)
                        .name(population_name(code)),
                );
            }
            if let Some((x, y)) = sample {
                plot_ui.points(
                    Points::new(PlotPoints::new(vec![[x, y]]))
                        .radius(7.0_f32)
                        .shape(MarkerShape::Diamond)
                        .color(egui::Color32::WHITE)
                        .name("you"),
                );
            }
        });
    // Compact super-population legend (matches the dot colors via a representative member).
    let mut seen: Vec<&str> = Vec::new();
    let mut legend: Vec<(String, egui::Color32)> = Vec::new();
    for (code, _, _) in reference {
        let Some(sup) = population_super(code) else { continue };
        if !seen.contains(&sup) {
            seen.push(sup);
            legend.push((population_name(sup), parse_hex_color(&population_color(code))));
        }
    }
    if !legend.is_empty() {
        ui.horizontal_wrapped(|ui| {
            for (name, color) in &legend {
                let (r, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                ui.painter().circle_filled(r.center(), 4.0, *color);
                ui.label(egui::RichText::new(name.as_str()).small());
                ui.add_space(6.0);
            }
        });
    }
}

/// egui_plot bar chart. Shared by the whole-genome and per-contig coverage views.
pub(crate) fn coverage_histogram_chart(ui: &mut egui::Ui, hist: &[u64], title: &str) {
    use egui_plot::{Bar, BarChart, Plot};
    ui.label(format!("Depth histogram — {title}  (depth ≥1; x = depth, y = bases)"));
    // Skip depth 0 (uncovered + reference-N): it typically dwarfs the coverage peak and would
    // flatten the rest of the distribution. That count is the table's NoCov / callable breakdown.
    let bars: Vec<Bar> = hist
        .iter()
        .enumerate()
        .skip(1)
        .map(|(depth, &count)| Bar::new(depth as f64, count as f64).width(0.9))
        .collect();
    let max_depth = hist.len().max(2) as f64;
    let max_count = hist.iter().skip(1).copied().max().unwrap_or(1) as f64;
    let chart = BarChart::new(bars).name("bases");
    // Fixed, non-interactive view: lock pan/zoom/scroll and pin the bounds to the data so the
    // axes can't drift into negative space or be dragged off-screen.
    Plot::new(format!("coverage_histogram_{title}"))
        .height(180.0)
        .allow_drag(false)
        .allow_zoom(false)
        .allow_scroll(false)
        .allow_boxed_zoom(false)
        .clamp_grid(true)
        .set_margin_fraction(egui::vec2(0.0, 0.05))
        .include_x(0.0)
        .include_x(max_depth)
        .include_y(0.0)
        .include_y(max_count)
        .show(ui, |plot_ui| plot_ui.bar_chart(chart));
}
