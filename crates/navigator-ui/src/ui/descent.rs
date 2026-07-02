//! YFull-YReport-style descent visualization (`impl NavigatorApp`): the subject's root→terminal
//! Y/mtDNA path drawn as labelled node badges, each followed by its defining SNPs as call-colored
//! chips. One generic renderer (built on the `DescentReport` aggregate) serves both lineages and
//! two densities — a compact path chain for the Simple view, the full chip wall for Advanced.

use super::*;

// Chip / badge palette (kept local; these are report-specific, not theme tokens).
const NODE_BG: egui::Color32 = egui::Color32::from_rgb(90, 150, 90); // backbone node badge (YFull green)
const TERMINAL_BG: egui::Color32 = egui::Color32::from_rgb(60, 120, 175); // the reported terminal
const DERIVED: egui::Color32 = egui::Color32::from_rgb(70, 150, 80); // sample carries the derived allele
const ANCESTRAL: egui::Color32 = egui::Color32::from_rgb(200, 120, 40); // ancestral (split not supported)
const NOCALL: egui::Color32 = egui::Color32::from_rgb(110, 110, 110); // no confident base

impl NavigatorApp {
    /// Render the descent report for `dna`, loading it lazily off the worker thread on first view and
    /// caching the result. `compact` = the Simple-view path chain; otherwise the full Advanced report.
    /// Additive and self-contained: shows a spinner while loading and a plain note when there's no
    /// placement, so it's safe to drop into any tab.
    pub(crate) fn descent_card(&mut self, ui: &mut egui::Ui, guid: SampleGuid, dna: DnaType, compact: bool) {
        self.ensure_descent(guid, dna);
        let entry = self
            .descent_reports
            .iter()
            .find(|(g, d, _)| *g == guid && *d == dna)
            .map(|(_, _, r)| r.is_some());
        match entry {
            Some(true) => {
                // Render inside a block so the report borrow ends before the (mutating) export button.
                let has_snps = {
                    let report = self
                        .descent_reports
                        .iter()
                        .find(|(g, d, _)| *g == guid && *d == dna)
                        .and_then(|(_, _, r)| r.as_ref())
                        .unwrap();
                    self.render_descent(ui, report, compact);
                    report.nodes.iter().any(|n| !n.snps.is_empty())
                };
                // Export the descent report (mirrors the on-screen grid) to TSV — full view only.
                if !compact && has_snps {
                    ui.add_space(6.0);
                    if ui.button(self.tr("descent.export")).clicked() {
                        let req = navigator_app::ExportRequest::DescentTsv(guid, dna);
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name(req.default_filename())
                            .add_filter("TSV", &["tsv"])
                            .save_file()
                        {
                            let _ = self.tx.send(Command::Export { request: req, path });
                            self.status = format!("Exporting {}…", req.label());
                        }
                    }
                }
            }
            // Loaded but empty → the variant profile isn't built yet (or has no placement). Offer the
            // one-time build, which persists and then feeds this report instantly.
            Some(false) => self.descent_build_prompt(ui, guid, dna),
            None => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(egui::RichText::new(self.tr("descent.loading")).weak());
                });
            }
        }
    }

    /// Shown when there's no cached report: a one-time "Build" affordance that runs (and persists)
    /// the variant profile this report is drawn from, or a plain note if it's built but unplaced.
    fn descent_build_prompt(&mut self, ui: &mut egui::Ui, guid: SampleGuid, dna: DnaType) {
        let (built, loading) = match dna {
            DnaType::Y => (self.y_profile.is_some(), self.y_profile_loading),
            DnaType::Mt => (self.mt_profile.is_some(), self.mt_profile_loading),
        };
        if loading {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(egui::RichText::new(self.tr("descent.building")).weak());
            });
            return;
        }
        if built {
            // Profile exists but yielded no terminal placement → nothing to draw.
            ui.label(egui::RichText::new(self.tr("descent.none")).weak());
            return;
        }
        ui.label(egui::RichText::new(self.tr("descent.buildHint")).weak());
        if ui.button(self.tr("descent.build")).clicked() {
            match dna {
                DnaType::Y => {
                    self.y_profile_loading = true;
                    let _ = self.tx.send(Command::BuildYProfile { biosample_guid: guid });
                }
                DnaType::Mt => {
                    self.mt_profile_loading = true;
                    let _ = self.tx.send(Command::BuildMtProfile { biosample_guid: guid });
                }
            }
        }
    }

    /// Inline replacement for the Simple-view brief's lineage trail: the compact, call-coloured
    /// descent path when the variant profile is built, otherwise the plain root→tip name trail (so
    /// the brief card is never empty and Simple mode never triggers an expensive build). Render-only;
    /// `subject_brief_view` pre-fires [`ensure_descent`] so the report loads.
    pub(crate) fn brief_descent_trail(&self, ui: &mut egui::Ui, guid: SampleGuid, lb: &LineageBrief) {
        let dna = match lb.kind {
            LineageKind::Paternal => DnaType::Y,
            LineageKind::Maternal => DnaType::Mt,
        };
        let report = self
            .descent_reports
            .iter()
            .find(|(g, d, _)| *g == guid && *d == dna)
            .and_then(|(_, _, r)| r.as_ref());
        if let Some(report) = report {
            ui.add_space(4.0);
            self.render_descent_compact(ui, report);
            return;
        }
        // Fallback until the profile is built: the plain collapsible root→tip trail.
        if lb.lineage_path.len() > 1 {
            ui.add_space(2.0);
            egui::CollapsingHeader::new(self.tr("brief.lineageTrail"))
                .id_salt(("brief_trail", matches!(lb.kind, LineageKind::Paternal)))
                .show(ui, |ui| {
                    ui.label(egui::RichText::new(lb.lineage_path.join("  ·  ")).small());
                });
        }
    }

    /// Fire a `LoadDescentReport` command if this (subject, DNA) report isn't already loaded or in
    /// flight. Idempotent — safe to call every frame.
    pub(crate) fn ensure_descent(&mut self, guid: SampleGuid, dna: DnaType) {
        let loaded = self.descent_reports.iter().any(|(g, d, _)| *g == guid && *d == dna);
        let loading = self.descent_loading.iter().any(|(g, d)| *g == guid && *d == dna);
        if !loaded && !loading {
            self.descent_loading.push((guid, dna));
            let _ = self.tx.send(Command::LoadDescentReport { guid, dna });
        }
    }

    fn render_descent(&self, ui: &mut egui::Ui, report: &DescentReport, compact: bool) {
        // The root carries no defining SNPs; drop empty nodes so the chain starts at the first call.
        if report.nodes.iter().all(|n| n.snps.is_empty()) {
            ui.label(egui::RichText::new(self.tr("descent.none")).weak());
            return;
        }
        if compact {
            self.render_descent_compact(ui, report);
        } else {
            self.render_descent_full(ui, report);
        }
    }

    /// Simple view: the lineage as a wrapped chain of node badges with a derived/total count each,
    /// terminal highlighted. No per-SNP chips.
    fn render_descent_compact(&self, ui: &mut egui::Ui, report: &DescentReport) {
        ui.horizontal_wrapped(|ui| {
            let mut first = true;
            for node in report.nodes.iter().filter(|n| !n.snps.is_empty()) {
                if !first {
                    ui.label(egui::RichText::new("›").weak());
                }
                first = false;
                let bg = if node.is_terminal { TERMINAL_BG } else { NODE_BG };
                ui.label(
                    egui::RichText::new(format!(" {} ", node.name))
                        .background_color(bg)
                        .color(egui::Color32::WHITE)
                        .strong(),
                );
                let (d, t) = node_counts(node);
                ui.label(egui::RichText::new(format!("{d}/{t}")).weak().small());
            }
        });
    }

    /// Advanced view: a colour legend, then per node a badge + a wrapped grid of its defining SNPs
    /// as chips coloured by the sample's call (derived / ancestral / no-call).
    fn render_descent_full(&self, ui: &mut egui::Ui, report: &DescentReport) {
        ui.horizontal_wrapped(|ui| {
            legend_chip(ui, self.tr("descent.derived"), DERIVED);
            legend_chip(ui, self.tr("descent.ancestral"), ANCESTRAL);
            legend_chip(ui, self.tr("descent.nocall"), NOCALL);
        });
        ui.add_space(8.0);

        for node in report.nodes.iter().filter(|n| !n.snps.is_empty()) {
            let (d, t) = node_counts(node);
            let bg = if node.is_terminal { TERMINAL_BG } else { NODE_BG };
            ui.horizontal(|ui| {
                ui.label(
                    egui::RichText::new(format!(" {} ", node.name))
                        .background_color(bg)
                        .color(egui::Color32::WHITE)
                        .strong(),
                );
                ui.label(egui::RichText::new(format!("{d}/{t} {}", self.tr("descent.derivedShort"))).weak().small());
            });
            ui.horizontal_wrapped(|ui| {
                for s in &node.snps {
                    let (chip_bg, state_label) = match s.state {
                        CallState::Derived => (DERIVED, self.tr("descent.derived")),
                        CallState::Ancestral => (ANCESTRAL, self.tr("descent.ancestral")),
                        CallState::NoCall => (NOCALL, self.tr("descent.nocall")),
                    };
                    ui.label(
                        egui::RichText::new(format!(" {} ", s.name))
                            .background_color(chip_bg)
                            .color(egui::Color32::WHITE)
                            .small(),
                    )
                    .on_hover_text(format!(
                        "{} · {}>{} @ {} · {}",
                        s.name, s.ancestral, s.derived, s.position, state_label
                    ));
                }
            });
            ui.add_space(8.0);
        }
    }
}

/// (derived, total) defining-SNP counts for one node.
fn node_counts(node: &navigator_app::NodeEvidence) -> (usize, usize) {
    let derived = node.snps.iter().filter(|s| matches!(s.state, CallState::Derived)).count();
    (derived, node.snps.len())
}

fn legend_chip(ui: &mut egui::Ui, label: &str, color: egui::Color32) {
    ui.label(
        egui::RichText::new(format!(" {label} "))
            .background_color(color)
            .color(egui::Color32::WHITE)
            .small(),
    );
}
