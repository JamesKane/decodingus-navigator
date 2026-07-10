//! Per-marker branch report card (`impl NavigatorApp`): the sample's genotype at every defining
//! marker of a chosen Y/mtDNA node's descendant subtree, for spot-checking placement accuracy and
//! exchanging observations with other researchers. Node-triggered (a text input + Load button, not
//! lazy), with a TSV export. Mirrors the descent card's worker/state wiring.

use super::*;

const DERIVED: egui::Color32 = egui::Color32::from_rgb(70, 150, 80); // carries the derived allele
const ANCESTRAL: egui::Color32 = egui::Color32::from_rgb(200, 120, 40); // ancestral (off this branch)
const NOCALL: egui::Color32 = egui::Color32::from_rgb(110, 110, 110); // no confident base

impl NavigatorApp {
    /// Render the branch-report card for `dna`: a node text input + Load, then the subtree table +
    /// TSV export. Self-contained (spinner while loading, hint when empty), safe to drop into a tab.
    pub(crate) fn branch_card(&mut self, ui: &mut egui::Ui, guid: SampleGuid, dna: DnaType) {
        let loading = self.branch_loading.iter().any(|(g, d)| *g == guid && *d == dna);

        // Node input + Load (Enter also loads).
        let mut do_load = false;
        ui.horizontal(|ui| {
            ui.label(self.tr("branch.node"));
            let resp = match dna {
                DnaType::Y => ui.add(
                    egui::TextEdit::singleline(&mut self.branch_node_y)
                        .desired_width(170.0)
                        .hint_text("R-M269"),
                ),
                DnaType::Mt => ui.add(
                    egui::TextEdit::singleline(&mut self.branch_node_mt)
                        .desired_width(170.0)
                        .hint_text("H2a"),
                ),
            };
            let node_empty = match dna {
                DnaType::Y => self.branch_node_y.trim().is_empty(),
                DnaType::Mt => self.branch_node_mt.trim().is_empty(),
            };
            let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            let clicked = ui
                .add_enabled(!loading && !node_empty, egui::Button::new(self.tr("branch.load")))
                .clicked();
            if (clicked || enter) && !loading && !node_empty {
                do_load = true;
            }
        });

        if do_load {
            let node = match dna {
                DnaType::Y => self.branch_node_y.trim().to_string(),
                DnaType::Mt => self.branch_node_mt.trim().to_string(),
            };
            self.branch_reports.retain(|(g, d, _)| !(*g == guid && *d == dna));
            self.branch_loading.push((guid, dna));
            let _ = self.tx.send(Command::LoadBranchReport { guid, dna, node, depth: None });
        }

        if self.branch_loading.iter().any(|(g, d)| *g == guid && *d == dna) {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(egui::RichText::new(self.tr("branch.loading")).weak());
            });
            return;
        }

        let present = self
            .branch_reports
            .iter()
            .any(|(g, d, r)| *g == guid && *d == dna && r.is_some());
        if !present {
            let no_alignment = self
                .branch_reports
                .iter()
                .any(|(g, d, r)| *g == guid && *d == dna && r.is_none());
            ui.add_space(4.0);
            let key = if no_alignment { "branch.noAlignment" } else { "branch.hint" };
            ui.label(egui::RichText::new(self.tr(key)).weak());
            return;
        }

        // Render inside a block so the report borrow ends before the (mutating) export button; the
        // block yields the formatted TSV + filename the button needs.
        let (tsv, fname) = {
            let report = self
                .branch_reports
                .iter()
                .find(|(g, d, _)| *g == guid && *d == dna)
                .and_then(|(_, _, r)| r.as_ref())
                .unwrap();
            self.render_branch(ui, report);
            let safe: String = report
                .root
                .chars()
                .map(|c| if c.is_ascii_alphanumeric() || c == '-' { c } else { '_' })
                .collect();
            (
                navigator_app::export::branch_report_tsv(report),
                format!("branch_{}_{}.tsv", safe, report.contig),
            )
        };
        ui.add_space(6.0);
        if ui.button(self.tr("branch.export")).clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name(&fname)
                .add_filter("TSV", &["tsv"])
                .save_file()
            {
                match std::fs::write(&path, tsv) {
                    Ok(()) => self.status = format!("{} {}", self.tr("branch.exported"), path.display()),
                    Err(e) => self.status = format!("write {}: {e}", path.display()),
                }
            }
        }
    }

    fn render_branch(&self, ui: &mut egui::Ui, report: &navigator_app::BranchReport) {
        let (d, a, n) = report.counts();
        ui.label(
            egui::RichText::new(format!(
                "{} · {} markers: {d} derived / {a} ancestral / {n} no-call · {}",
                report.root,
                report.rows.len(),
                if report.gvcf_backed { "gVCF" } else { "pileup" },
            ))
            .small()
            .weak(),
        );
        ui.add_space(4.0);

        // A subtree rooted at a shallow node (R-M269, or the tree root) carries tens of thousands of
        // markers — a Grid lays out every row per frame and beach-balls. Fixed-width columns through
        // ScrollArea::show_rows build only the visible slice (same idiom as the consensus-panel table).
        const W_NODE: f32 = 110.0;
        const W_MARKER: f32 = 90.0;
        const W_POS: f32 = 80.0;
        const W_ALLELES: f32 = 60.0;
        const W_OBS: f32 = 30.0;
        const W_STATUS: f32 = 70.0;
        const W_AD: f32 = 56.0;
        const W_NUM: f32 = 40.0;
        let row_h = ui.text_style_height(&egui::TextStyle::Small) + 4.0;
        let head = |ui: &mut egui::Ui, w: f32, t: &str| {
            ui.add_sized([w, row_h], egui::Label::new(egui::RichText::new(t).strong().small()));
        };
        ui.horizontal(|ui| {
            head(ui, W_NODE, "node");
            head(ui, W_MARKER, "marker");
            head(ui, W_POS, "pos");
            head(ui, W_ALLELES, "anc>der");
            head(ui, W_OBS, "obs");
            head(ui, W_STATUS, "status");
            head(ui, W_AD, "AD");
            head(ui, W_NUM, "DP");
            head(ui, W_NUM, "GQ");
            ui.label(egui::RichText::new("note").strong().small());
        });
        egui::ScrollArea::vertical()
            .max_height(320.0)
            .auto_shrink([false, true])
            .show_rows(ui, row_h, report.rows.len(), |ui, range| {
                let dot = |s: String| if s.is_empty() { ".".to_string() } else { s };
                let opt_u = |o: Option<u32>| o.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string());
                let cell = |ui: &mut egui::Ui, w: f32, t: egui::RichText| {
                    ui.add_sized([w, row_h], egui::Label::new(t.small()).truncate());
                };
                for i in range {
                    let r = &report.rows[i];
                    let (color, label) = match r.state {
                        navigator_app::CallState::Derived => (DERIVED, "derived"),
                        navigator_app::CallState::Ancestral => (ANCESTRAL, "ancestral"),
                        navigator_app::CallState::NoCall => (NOCALL, "no-call"),
                    };
                    ui.horizontal(|ui| {
                        cell(ui, W_NODE, egui::RichText::new(&r.node));
                        cell(ui, W_MARKER, egui::RichText::new(&r.marker));
                        cell(ui, W_POS, egui::RichText::new(r.position.to_string()));
                        cell(ui, W_ALLELES, egui::RichText::new(format!("{}>{}", r.ancestral, r.derived)));
                        cell(
                            ui,
                            W_OBS,
                            egui::RichText::new(dot(r.observed_base.map(|c| c.to_string()).unwrap_or_default())),
                        );
                        cell(ui, W_STATUS, egui::RichText::new(label).color(color));
                        cell(
                            ui,
                            W_AD,
                            egui::RichText::new(dot(r.ad.map(|(rf, al)| format!("{rf},{al}")).unwrap_or_default())),
                        );
                        cell(ui, W_NUM, egui::RichText::new(opt_u(r.dp)));
                        cell(ui, W_NUM, egui::RichText::new(opt_u(r.gq)));
                        ui.label(egui::RichText::new(&r.note).small().weak());
                    });
                }
            });
    }
}
