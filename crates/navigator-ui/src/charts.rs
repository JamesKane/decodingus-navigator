//! Pure ancestry / genome-region visualization helpers, extracted from the UI shell. Leaf drawing
//! functions over `egui` with no `App`/`self` state — easy to read, test, and reuse independently of
//! the view code that calls them.

use eframe::egui;
use navigator_app::{AncestryResult, AncestrySegment, AssetStatus, ChromosomeRegions, GenomeRegions, IbdSegment, SuperPopulationSummary};
use navigator_domain::ancestry::{population_color, population_name, population_super};

use crate::ui::ACCENT;

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
            ui.allocate_ui(egui::vec2(label_w, bar_h), |ui| ui.label(egui::RichText::new(chr).small()));
            let (rect, resp) = ui.allocate_exact_size(egui::vec2(bar_w, bar_h), egui::Sense::hover());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 2.0, egui::Color32::from_gray(30));
            let mut hover: Option<String> = None;
            let hover_x = resp.hover_pos().map(|p| p.x);
            for &s in segs.iter() {
                let x0 = rect.left() + (s.start_position.max(0) as f32 / chr_len) * rect.width();
                let x1 = rect.left() + (s.end_position.max(0) as f32 / chr_len) * rect.width();
                let seg = egui::Rect::from_min_max(egui::pos2(x0, rect.top()), egui::pos2(x1.max(x0 + 1.5), rect.bottom()));
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
        let Ok(n) = s.contig.trim_start_matches("chr").parse::<i64>() else { continue };
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
        let hi = copies.iter().flatten().map(|s| s.end).max().unwrap_or(lo + 1).max(lo + 1);
        let span = (hi - lo).max(1) as f32;
        ui.horizontal(|ui| {
            ui.allocate_ui(egui::vec2(label_w, copy_h * 2.0 + gap), |ui| ui.label(format!("chr{n}")));
            let (rect, _) = ui.allocate_exact_size(egui::vec2(bar_w, copy_h * 2.0 + gap), egui::Sense::hover());
            let painter = ui.painter_at(rect);
            for (c, segs) in copies.iter().enumerate() {
                let top = rect.top() + c as f32 * (copy_h + gap);
                let track = egui::Rect::from_min_size(egui::pos2(rect.left(), top), egui::vec2(rect.width(), copy_h));
                painter.rect_filled(track, 2.0, egui::Color32::from_gray(30));
                for s in segs {
                    let x0 = track.left() + (s.start - lo) as f32 / span * track.width();
                    let x1 = track.left() + (s.end - lo) as f32 / span * track.width();
                    let seg_rect = egui::Rect::from_min_max(egui::pos2(x0, track.top()), egui::pos2(x1.max(x0 + 1.0), track.bottom()));
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
            ui.painter().circle_filled(r.center(), 4.0, parse_hex_color(&population_color(code)));
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

/// Draw a donut chart of the super-population proportions (one wedge per super-population,
/// colored by continent), with the dominant share in the centre.
pub(crate) fn draw_ancestry_donut(ui: &mut egui::Ui, summary: &[SuperPopulationSummary]) {
    let size = 120.0;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(size, size), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let c = rect.center();
    let (r_out, r_in) = (size * 0.46, size * 0.28);
    let total: f32 = summary.iter().map(|s| s.percentage as f32).sum::<f32>().max(1.0);
    let mut a0 = -std::f32::consts::FRAC_PI_2; // start at 12 o'clock
    for s in summary {
        if s.percentage < 0.5 {
            continue;
        }
        let a1 = a0 + (s.percentage as f32 / total) * std::f32::consts::TAU;
        let code = s.populations.first().and_then(|c| population_super(c)).unwrap_or("");
        let mut pts = arc_points(c, r_out, a0, a1, 32);
        pts.extend(arc_points(c, r_in, a1, a0, 32)); // inner arc, reversed → closed ring sector
        painter.add(egui::epaint::PathShape {
            points: pts,
            closed: true,
            fill: parse_hex_color(&population_color(code)),
            stroke: egui::epaint::PathStroke::NONE,
        });
        a0 = a1;
    }
    if let Some(top) = summary.first() {
        painter.text(
            c,
            egui::Align2::CENTER_CENTER,
            format!("{:.0}%", top.percentage),
            egui::FontId::proportional(18.0),
            egui::Color32::WHITE,
        );
    }
}

/// Draw a detailed ancestry breakdown (the fine-population or ancient-component report): the
/// estimate's `components`, sorted by share, as a name/percentage grid with a proportion bar, plus a
/// provenance line (method + SNP count). `id_salt` keeps each report's grid distinct.
pub(crate) fn draw_population_components(ui: &mut egui::Ui, result: &AncestryResult, id_salt: &str, top_n: usize) {
    let mut comps: Vec<(&str, f64)> =
        result.components.iter().filter(|c| c.percentage >= 0.05).map(|c| (c.population_code.as_str(), c.percentage)).collect();
    comps.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    if comps.is_empty() {
        ui.label(egui::RichText::new("No components above 0.05%.").weak());
        return;
    }
    let max = comps.first().map(|c| c.1).unwrap_or(1.0).max(1e-6);
    egui::Grid::new(format!("{id_salt}_grid")).striped(true).num_columns(3).show(ui, |ui| {
        for (name, pct) in comps.iter().take(top_n) {
            ui.label(*name);
            ui.label(egui::RichText::new(format!("{pct:.1}%")).strong());
            // A small proportion bar (relative to the top component) for at-a-glance ranking.
            let (rect, _) = ui.allocate_exact_size(egui::vec2(120.0, 10.0), egui::Sense::hover());
            let p = ui.painter_at(rect);
            p.rect_filled(rect, 2.0, ui.visuals().faint_bg_color);
            let mut fill = rect;
            fill.set_width(rect.width() * (pct / max) as f32);
            p.rect_filled(fill, 2.0, ACCENT.gamma_multiply(0.7));
            ui.end_row();
        }
    });
    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(format!("{} · {}/{} SNPs · confidence {:.0}%", result.method, result.snps_with_genotype, result.snps_analyzed, result.confidence_level * 100.0))
            .weak()
            .small(),
    );
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
fn parse_hex_color(hex: &str) -> egui::Color32 {
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

/// Giemsa-stain → color for the cytoband ideogram (the standard UCSC palette, tuned for the dark
/// theme): `gneg` light → `gpos100` near-black, `acen` (centromere) red, `gvar`/`stalk` tinted.
fn stain_color(stain: &str) -> egui::Color32 {
    match stain {
        "gneg" => egui::Color32::from_gray(225),
        "gpos25" => egui::Color32::from_gray(170),
        "gpos50" => egui::Color32::from_gray(120),
        "gpos75" => egui::Color32::from_gray(80),
        "gpos100" => egui::Color32::from_gray(45),
        "acen" => egui::Color32::from_rgb(200, 70, 70),
        "gvar" => egui::Color32::from_rgb(120, 140, 185),
        "stalk" => egui::Color32::from_rgb(110, 165, 160),
        _ => egui::Color32::from_gray(140),
    }
}

/// The chromosomes to draw, in karyotype order (1–22, X, Y); non-nuclear / alt / random contigs
/// are skipped. Tolerates a `chr` prefix on the names.
fn karyotype_order(regions: &GenomeRegions) -> Vec<(&String, &ChromosomeRegions)> {
    fn rank(name: &str) -> Option<u32> {
        let s = name.strip_prefix("chr").unwrap_or(name);
        match s {
            "X" => Some(23),
            "Y" => Some(24),
            _ => s.parse::<u32>().ok().filter(|n| (1..=22).contains(n)),
        }
    }
    let mut v: Vec<(u32, &String, &ChromosomeRegions)> =
        regions.chromosomes.iter().filter_map(|(n, c)| rank(n).map(|r| (r, n, c))).collect();
    v.sort_by_key(|(r, _, _)| *r);
    v.into_iter().map(|(_, n, c)| (n, c)).collect()
}

/// A compact legend mapping Giemsa stains to their ideogram colors.
pub(crate) fn ideogram_legend(ui: &mut egui::Ui) {
    ui.horizontal_wrapped(|ui| {
        for (label, stain) in
            [("gneg", "gneg"), ("gpos50", "gpos50"), ("gpos100", "gpos100"), ("centromere", "acen"), ("gvar", "gvar")]
        {
            let (r, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
            ui.painter().rect_filled(r, 2.0, stain_color(stain));
            ui.label(egui::RichText::new(label).small().weak());
            ui.add_space(8.0);
        }
    });
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
            ui.colored_label(col, egui::RichText::new(format!("{} {mark}", a.name)).small()).on_hover_text(hover);
            ui.add_space(4.0);
        }
    });
}

/// Draw the chromosome ideogram: one horizontal bar per chromosome (scaled to the longest), its
/// cytobands as Giemsa-stained segments, with a hover tooltip naming the band under the cursor.
pub(crate) fn draw_ideogram(ui: &mut egui::Ui, regions: &GenomeRegions) {
    let order = karyotype_order(regions);
    if order.is_empty() {
        ui.label(egui::RichText::new("No chromosome data.").weak());
        return;
    }
    let max_len = order.iter().map(|(_, c)| c.length).max().unwrap_or(1).max(1) as f32;
    let label_w = 30.0;
    let row_h = 16.0;
    let full_w = ui.available_width().min(760.0);
    let bar_area = (full_w - label_w - 6.0).max(60.0);
    let text_color = ui.visuals().text_color();

    for (name, c) in &order {
        if c.length <= 0 {
            continue;
        }
        let (rect, resp) = ui.allocate_exact_size(egui::vec2(full_w, row_h), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        let short = name.strip_prefix("chr").unwrap_or(name);
        painter.text(
            egui::pos2(rect.left() + 2.0, rect.center().y),
            egui::Align2::LEFT_CENTER,
            short,
            egui::FontId::proportional(11.0),
            text_color,
        );

        let x0 = rect.left() + label_w;
        let bar_w = bar_area * (c.length as f32 / max_len);
        let bar = egui::Rect::from_min_size(egui::pos2(x0, rect.top() + 2.0), egui::vec2(bar_w, row_h - 4.0));
        painter.rect_filled(bar, 3.0, egui::Color32::from_gray(28));
        let len = c.length as f32;
        for b in &c.cytobands {
            let bx0 = x0 + bar_w * (b.start as f32 / len);
            let bx1 = x0 + bar_w * (b.end as f32 / len);
            let seg = egui::Rect::from_min_max(egui::pos2(bx0, bar.top()), egui::pos2(bx1.max(bx0 + 0.5), bar.bottom()));
            painter.rect_filled(seg, 0.0, stain_color(&b.stain));
        }
        painter.rect_stroke(bar, 3.0, egui::Stroke::new(1.0, egui::Color32::from_gray(70)));

        // Hover → the band (and Mb position) under the cursor.
        if let Some(pos) = resp.hover_pos() {
            if pos.x >= x0 && pos.x <= x0 + bar_w && bar_w > 0.0 {
                let g = (((pos.x - x0) / bar_w).clamp(0.0, 1.0) * len) as i64;
                let band = c.cytobands.iter().find(|b| b.start <= g && g < b.end);
                let name = band.map(|b| format!("{short}{}", b.name)).unwrap_or_else(|| short.to_string());
                let stain = band.map(|b| b.stain.as_str()).unwrap_or("");
                resp.on_hover_text(format!("{name}  ·  {stain}  ·  {:.1} Mb", g as f64 / 1e6));
            }
        }
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
