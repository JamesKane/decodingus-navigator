//! Generic, reusable UI widgets and small formatters extracted from the view shell — leaf functions
//! over `egui` with no `App`/`self` state (tables, cards, chips, stat tiles, dropdowns, and the
//! number/guid formatters they share). The view code in `ui.rs` composes these.

use eframe::egui;
use navigator_app::{CallState, HaploAssignment};
use navigator_domain::variants::VariantCall;
use navigator_domain::workspace::Biosample;

use crate::ui::ACCENT;

/// Short, stable subject id for the table's ID column (first 8 chars of the guid + ellipsis).
pub(crate) fn short_guid(b: &Biosample) -> String {
    let s = b.guid.0.to_string();
    if s.len() > 9 {
        format!("{}…", &s[..9])
    } else {
        s
    }
}

/// Paint a table header row at the column offsets used by [`table_row`].
pub(crate) fn table_header(ui: &mut egui::Ui, cols: &[(&str, f32)]) {
    let total_w: f32 = cols.iter().map(|c| c.1).sum();
    let (rect, _) = ui.allocate_exact_size(egui::vec2(total_w, 24.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let mut x = rect.left() + 8.0;
    let color = ui.visuals().weak_text_color();
    for (name, w) in cols {
        painter.text(
            egui::pos2(x, rect.center().y),
            egui::Align2::LEFT_CENTER,
            name,
            egui::FontId::proportional(12.5),
            color,
        );
        x += w;
    }
    painter.hline(rect.x_range(), rect.bottom(), egui::Stroke::new(1.0, ui.visuals().widgets.noninteractive.bg_stroke.color));
}

/// Paint one clickable table row; returns true when clicked. `status_col` (if any) is rendered
/// as a small accent-coloured badge (the "Status" cell).
pub(crate) fn table_row(ui: &mut egui::Ui, cols: &[(&str, f32)], cells: &[String], selected: bool, status_col: Option<usize>) -> bool {
    let total_w: f32 = cols.iter().map(|c| c.1).sum();
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(total_w, 28.0), egui::Sense::click());
    let painter = ui.painter_at(rect);
    if selected {
        painter.rect_filled(rect, 4.0, ACCENT.gamma_multiply(0.6));
    } else if resp.hovered() {
        painter.rect_filled(rect, 4.0, ui.visuals().faint_bg_color);
    }
    let text_color = if selected { egui::Color32::WHITE } else { ui.visuals().text_color() };
    let mut x = rect.left() + 8.0;
    for (i, ((_, w), val)) in cols.iter().zip(cells).enumerate() {
        let cy = rect.center().y;
        if Some(i) == status_col && val != "-" {
            // a muted pill behind the status text
            let galley = painter.layout_no_wrap(val.clone(), egui::FontId::proportional(11.5), egui::Color32::from_rgb(225, 190, 90));
            let pad = egui::vec2(7.0, 3.0);
            let pill = egui::Rect::from_min_size(egui::pos2(x, cy - galley.size().y / 2.0 - pad.y), galley.size() + pad * 2.0);
            painter.rect_filled(pill, 8.0, egui::Color32::from_rgb(70, 58, 28));
            painter.galley(pill.min + pad, galley, egui::Color32::PLACEHOLDER);
        } else {
            // Elide to the column width (single line, trailing …) so long values — e.g. an
            // ISOGG-longhand Y haplogroup — can't spill into the next column.
            let mut job = egui::text::LayoutJob::single_section(
                val.clone(),
                egui::TextFormat { font_id: egui::FontId::proportional(13.0), color: text_color, ..Default::default() },
            );
            job.wrap = egui::text::TextWrapping {
                max_width: (w - 12.0).max(0.0),
                max_rows: 1,
                overflow_character: Some('…'),
                ..Default::default()
            };
            let galley = ui.fonts(|f| f.layout_job(job));
            painter.galley(egui::pos2(x, cy - galley.size().y / 2.0), galley, text_color);
        }
        x += w;
    }
    resp.clicked()
}

/// A rounded section card with an optional bold title (the Data Sources look).
pub(crate) fn card(ui: &mut egui::Ui, title: &str, body: impl FnOnce(&mut egui::Ui)) {
    egui::Frame::group(ui.style())
        .fill(ui.visuals().faint_bg_color)
        .rounding(egui::Rounding::same(8.0))
        .inner_margin(egui::Margin::same(14.0))
        .show(ui, |ui| {
            ui.set_width(ui.available_width());
            if !title.is_empty() {
                ui.label(egui::RichText::new(title).strong().size(15.0));
                ui.add_space(8.0);
            }
            body(ui);
        });
}

/// A small rounded chip/badge (provider tag, Y/mt badge).
pub(crate) fn chip(ui: &mut egui::Ui, text: &str, bg: egui::Color32, fg: egui::Color32) -> egui::Response {
    let font = egui::FontId::proportional(11.5);
    let galley = ui.painter().layout_no_wrap(text.to_string(), font, fg);
    let pad = egui::vec2(7.0, 3.0);
    let (rect, response) = ui.allocate_exact_size(galley.size() + pad * 2.0, egui::Sense::hover());
    ui.painter().rect_filled(rect, 6.0, bg);
    ui.painter().galley(rect.min + pad, galley, egui::Color32::PLACEHOLDER);
    response
}

/// 3-letter provider abbreviation for the run chip (PACBIO → PAC).
pub(crate) fn provider_abbrev(platform: &str) -> String {
    let p = platform.trim();
    if p.is_empty() || p.eq_ignore_ascii_case("unknown") {
        "SEQ".into()
    } else {
        p.chars().take(3).collect::<String>().to_uppercase()
    }
}

/// Compact read count: 9_900 → "9.9K", 1_200_000 → "1.2M".
pub(crate) fn fmt_reads(n: Option<i64>) -> String {
    match n {
        None => "—".into(),
        Some(v) if v >= 1_000_000 => format!("{:.1}M", v as f64 / 1e6),
        Some(v) if v >= 1_000 => format!("{:.1}K", v as f64 / 1e3),
        Some(v) => v.to_string(),
    }
}

/// A centered empty-state placeholder for a work area with no selection.
pub(crate) fn empty_state(ui: &mut egui::Ui, title: &str, hint: &str) {
    ui.add_space(56.0);
    ui.vertical_centered(|ui| {
        ui.heading(title);
        ui.label(egui::RichText::new(hint).weak());
    });
}

/// A dashboard stat tile: a big number over a muted label, in a rounded card.
pub(crate) fn stat_card(ui: &mut egui::Ui, label: &str, value: usize) {
    egui::Frame::group(ui.style())
        .fill(ui.visuals().faint_bg_color)
        .inner_margin(egui::Margin::symmetric(18.0, 14.0))
        .show(ui, |ui| {
            ui.set_min_width(120.0);
            ui.vertical(|ui| {
                ui.label(egui::RichText::new(value.to_string()).size(28.0).strong());
                ui.label(egui::RichText::new(label).weak());
            });
        });
}

/// Format an optional mean/median depth (one decimal), "—" when not computed.
pub(crate) fn fmt_depth(o: Option<f64>) -> String {
    o.map(|v| format!("{v:.1}")).unwrap_or_else(|| "—".into())
}

/// Format an optional fraction (0–1) as a percentage, "—" when not computed.
pub(crate) fn fmt_pct(o: Option<f64>) -> String {
    o.map(|v| format!("{:.1}%", v * 100.0)).unwrap_or_else(|| "—".into())
}

/// Render a haplogroup assignment: terminal + lineage + alternatives, then the child
/// branches with per-SNP evidence that explains why descent stopped.
pub(crate) fn show_assignment(ui: &mut egui::Ui, a: &HaploAssignment) {
    let Some(top) = a.ranked.first() else {
        ui.label("No match."); // free helper (no `self`); i18n when it takes a `lang` param
        return;
    };
    ui.label(format!("Haplogroup: {}   ({}/{} mutations, score {:.3})", top.name, top.matched, top.expected, top.score));
    ui.label(format!("Lineage: {}", top.lineage.join(" › ")));
    let alts: Vec<String> = a.ranked.iter().skip(1).take(3).map(|r| format!("{} ({:.3})", r.name, r.score)).collect();
    if !alts.is_empty() {
        ui.label(format!("Alternatives: {}", alts.join(", ")));
    }
    for b in &a.branches {
        egui::CollapsingHeader::new(format!("child {} — {}/{} SNPs derived", b.name, b.derived, b.snps.len()))
            .id_salt(("branch", &b.name))
            .show(ui, |ui| {
                egui::Grid::new(("branch_snps", &b.name)).striped(true).num_columns(3).show(ui, |ui| {
                    for s in &b.snps {
                        ui.label(&s.name);
                        ui.label(format!("{}{}>{}", s.position, s.ancestral, s.derived));
                        let (txt, col) = match s.state {
                            CallState::Derived => ("derived", egui::Color32::from_rgb(60, 160, 60)),
                            CallState::Ancestral => ("ancestral", egui::Color32::from_rgb(170, 120, 40)),
                            CallState::NoCall => ("no-call", egui::Color32::GRAY),
                        };
                        ui.colored_label(col, txt);
                        ui.end_row();
                    }
                });
            });
    }
}

/// A readable "change" string for a variant call, covering the indel forms the mtDNA
/// derivation stores (one allele empty).
pub(crate) fn variant_change(c: &VariantCall) -> String {
    if c.alternate.is_empty() {
        format!("{}del", c.reference) // deletion
    } else if c.reference.is_empty() {
        format!("ins{}", c.alternate) // insertion
    } else {
        format!("{}>{}", c.reference, c.alternate) // substitution
    }
}

/// Trim a string, returning `None` when empty.
pub(crate) fn opt(s: &str) -> Option<String> {
    let t = s.trim();
    (!t.is_empty()).then(|| t.to_string())
}

/// A labeled dropdown that sets `value` to one of `options` (string codes).
pub(crate) fn combo(ui: &mut egui::Ui, label: &str, id: &str, value: &mut String, options: &[&str]) {
    ui.horizontal(|ui| {
        ui.label(label);
        egui::ComboBox::from_id_salt(id).selected_text(value.clone()).show_ui(ui, |ui| {
            for opt in options {
                ui.selectable_value(value, opt.to_string(), *opt);
            }
        });
    });
}
