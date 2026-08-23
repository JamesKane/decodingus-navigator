//! Generic, reusable UI widgets and small formatters, taken out of the view shell. They are leaf
//! functions over `egui` with no `App` or `self` state. They cover tables, cards, chips, stat
//! tiles, dropdowns, and the number and guid formatters they share. The view code in `ui.rs`
//! composes these.

use eframe::egui;
use navigator_app::{CallState, HaploAssignment};
use navigator_domain::variants::VariantCall;

/// Click-to-sort state, and the inline filter state of each column, for one sticky-header table.
/// The `App` holds it, so that the chosen sort column, the direction, and the typed filters survive
/// across frames.
#[derive(Default)]
pub(crate) struct TableControls {
    sort_col: Option<usize>,
    ascending: bool,
    /// The filter text of each column, with the column as the index. It grows on demand, and it
    /// drops an empty entry.
    filters: Vec<String>,
}

impl TableControls {
    /// Start with a sort by `col`, from small to large (natural order), before the user clicks any
    /// header. The table then opens in a sensible order: numbered FTDNA kits flow 1, 2, 10, 100.
    pub(crate) fn sorted_by(col: usize) -> Self {
        Self {
            sort_col: Some(col),
            ascending: true,
            filters: Vec::new(),
        }
    }

    /// A mutable handle to the filter text of column `col`. It grows the backing store as needed.
    pub(crate) fn filter_mut(&mut self, col: usize) -> &mut String {
        if self.filters.len() <= col {
            self.filters.resize(col + 1, String::new());
        }
        &mut self.filters[col]
    }

    /// The raw filter text of each column, so that the row caches can detect a change. A compare of
    /// this makes no allocation, and a call to [`Self::filter_norm`] for each column on every frame
    /// does.
    pub(crate) fn filters_raw(&self) -> &[String] {
        &self.filters
    }

    /// Trimmed, lower-case filter text for `col`, empty when nothing sets it. It is ready for
    /// `contains`.
    pub(crate) fn filter_norm(&self, col: usize) -> String {
        self.filters
            .get(col)
            .map(|s| s.trim().to_lowercase())
            .unwrap_or_default()
    }

    pub(crate) fn sort_col(&self) -> Option<usize> {
        self.sort_col
    }

    pub(crate) fn ascending(&self) -> bool {
        self.ascending
    }

    /// The first click on a column sorts it from small to large. A click on the active column
    /// flips the direction.
    pub(crate) fn toggle_sort(&mut self, col: usize) {
        if self.sort_col == Some(col) {
            self.ascending = !self.ascending;
        } else {
            self.sort_col = Some(col);
            self.ascending = true;
        }
    }

    fn arrow(&self, col: usize) -> &'static str {
        match self.sort_col {
            Some(c) if c == col => {
                // ⏶/⏷ (U+23F6/U+23F7), not ▲/▼ (U+25B2/U+25BC). Only the first pair has glyphs in
                // the Proportional family of egui. See
                // `every_source_string_literal_is_renderable`.
                if self.ascending {
                    " ⏶"
                } else {
                    " ⏷"
                }
            }
            _ => "",
        }
    }
}

/// Draw one sticky-header cell: a clickable sort label (click to sort, click again to flip) over an
/// inline filter input for that column. Pass `filterable = false` for a column that takes no filter
/// (for example an actions column), to draw the label only.
pub(crate) fn sortable_header(ui: &mut egui::Ui, ctl: &mut TableControls, col: usize, label: &str, filterable: bool) {
    ui.vertical(|ui| {
        let title = format!("{label}{}", ctl.arrow(col));
        if ui
            .add(egui::Label::new(egui::RichText::new(title).strong()).sense(egui::Sense::click()))
            .on_hover_text("Click to sort")
            .clicked()
        {
            ctl.toggle_sort(col);
        }
        if filterable {
            ui.add(
                egui::TextEdit::singleline(ctl.filter_mut(col))
                    .hint_text("filter")
                    .desired_width(f32::INFINITY)
                    .font(egui::TextStyle::Small),
            );
        }
    });
}

/// Natural ("human") ordering. It compares strings so that a run of digits sorts by numeric value,
/// and not lexically: `Kit-2` < `Kit-10` < `Kit-100`. A run that is not digits compares without
/// regard to case. For a run of digits it compares the trimmed length first, then the digits, so a
/// number of any length is safe.
pub(crate) fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let mut ai = a.chars().peekable();
    let mut bi = b.chars().peekable();
    loop {
        match (ai.peek().copied(), bi.peek().copied()) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(ca), Some(cb)) if ca.is_ascii_digit() && cb.is_ascii_digit() => {
                let na = take_digits(&mut ai);
                let nb = take_digits(&mut bi);
                let ta = na.trim_start_matches('0');
                let tb = nb.trim_start_matches('0');
                let ord = ta.len().cmp(&tb.len()).then_with(|| ta.cmp(tb));
                if ord != Ordering::Equal {
                    return ord;
                }
                // Equal value. Keep the order stable by the count of zeros at the start, so `01`
                // comes before `1`.
                let ord = (na.len() - ta.len()).cmp(&(nb.len() - tb.len()));
                if ord != Ordering::Equal {
                    return ord;
                }
            }
            (Some(ca), Some(cb)) => {
                let ord = ca.to_ascii_lowercase().cmp(&cb.to_ascii_lowercase());
                if ord != Ordering::Equal {
                    return ord;
                }
                ai.next();
                bi.next();
            }
        }
    }
}

fn take_digits(it: &mut std::iter::Peekable<std::str::Chars>) -> String {
    let mut s = String::new();
    while let Some(&c) = it.peek() {
        if !c.is_ascii_digit() {
            break;
        }
        s.push(c);
        it.next();
    }
    s
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

/// Capitalize the first character of a string (for sentence-casing a generated lowercase phrase).
pub(crate) fn capitalize_first(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        Some(c) => c.to_uppercase().collect::<String>() + chars.as_str(),
        None => String::new(),
    }
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

/// True when the user activated a button inside a clickable row. Either the user clicked the button
/// directly, or the row took the press while the pointer was over the button. Both lists of
/// selectable rows (runs, and their alignments) need the same rule, so it has one definition.
pub(crate) fn button_hit(button: &Option<egui::Response>, row_clicked: bool) -> bool {
    button
        .as_ref()
        .is_some_and(|r| r.clicked() || (row_clicked && r.contains_pointer()))
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

/// Format an optional mean or median depth (one decimal). `"—"` when nothing computed it.
pub(crate) fn fmt_depth(o: Option<f64>) -> String {
    o.map(|v| format!("{v:.1}")).unwrap_or_else(|| "—".into())
}

/// Format an optional fraction (0–1) as a percentage. `"—"` when nothing computed it.
pub(crate) fn fmt_pct(o: Option<f64>) -> String {
    o.map(|v| format!("{:.1}%", v * 100.0)).unwrap_or_else(|| "—".into())
}

/// Draw a haplogroup assignment: the terminal, the lineage and the alternatives. Then draw the
/// child branches, with the evidence at each SNP that explains why descent stopped.
pub(crate) fn show_assignment(ui: &mut egui::Ui, a: &HaploAssignment) {
    let Some(top) = a.ranked.first() else {
        ui.label("No match."); // free helper (no `self`); i18n when it takes a `lang` param
        return;
    };
    ui.label(format!(
        "Haplogroup: {}   ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    ));
    ui.label(format!("Lineage: {}", top.lineage.join(" › ")));
    let alts: Vec<String> = a
        .ranked
        .iter()
        .skip(1)
        .take(3)
        .map(|r| format!("{} ({:.3})", r.name, r.score))
        .collect();
    if !alts.is_empty() {
        ui.label(format!("Alternatives: {}", alts.join(", ")));
    }
    for b in &a.branches {
        egui::CollapsingHeader::new(format!(
            "child {} — {}/{} SNPs derived",
            b.name,
            b.derived,
            b.snps.len()
        ))
        .id_salt(("branch", &b.name))
        .show(ui, |ui| {
            egui::Grid::new(("branch_snps", &b.name))
                .striped(true)
                .num_columns(3)
                .show(ui, |ui| {
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

/// A readable "change" string for a variant call. It also covers the indel forms that the mtDNA
/// derivation stores, where one allele is empty.
pub(crate) fn variant_change(c: &VariantCall) -> String {
    if c.alternate.is_empty() {
        format!("{}del", c.reference) // deletion
    } else if c.reference.is_empty() {
        format!("ins{}", c.alternate) // insertion
    } else {
        format!("{}>{}", c.reference, c.alternate) // substitution
    }
}

/// Trim a string. Returns `None` when the result is empty.
pub(crate) fn opt(s: &str) -> Option<String> {
    let t = s.trim();
    (!t.is_empty()).then(|| t.to_string())
}

/// A labeled dropdown that sets `value` to one of `options` (string codes).
pub(crate) fn combo(ui: &mut egui::Ui, label: &str, id: &str, value: &mut String, options: &[&str]) {
    ui.horizontal(|ui| {
        ui.label(label);
        egui::ComboBox::from_id_salt(id)
            .selected_text(value.clone())
            .show_ui(ui, |ui| {
                for opt in options {
                    ui.selectable_value(value, opt.to_string(), *opt);
                }
            });
    });
}

#[cfg(test)]
mod tests {
    use super::natural_cmp;
    use std::cmp::Ordering;

    #[test]
    fn numbered_kits_sort_by_value_not_lexically() {
        let mut kits = vec!["Kit-100", "Kit-2", "Kit-10", "Kit-1"];
        kits.sort_by(|a, b| natural_cmp(a, b));
        assert_eq!(kits, vec!["Kit-1", "Kit-2", "Kit-10", "Kit-100"]);
    }

    #[test]
    fn pure_numbers_and_long_numbers() {
        assert_eq!(natural_cmp("9", "10"), Ordering::Less);
        // Longer-than-u64 digit runs still compare by value (length then digits), no overflow.
        assert_eq!(
            natural_cmp("99999999999999999999", "100000000000000000000"),
            Ordering::Less
        );
    }

    #[test]
    fn case_insensitive_and_mixed() {
        assert_eq!(natural_cmp("alpha", "ALPHA"), Ordering::Equal);
        assert_eq!(natural_cmp("file2", "file10"), Ordering::Less);
        assert_eq!(natural_cmp("b", "a10"), Ordering::Greater);
    }

    #[test]
    fn prefix_is_less() {
        assert_eq!(natural_cmp("Kit", "Kit-1"), Ordering::Less);
    }
}
