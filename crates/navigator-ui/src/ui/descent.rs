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
        match self.descent_reports.iter().find(|(g, d, _)| *g == guid && *d == dna) {
            Some((_, _, Some(report))) => self.render_descent(ui, report, compact),
            Some((_, _, None)) => {
                ui.label(egui::RichText::new(self.tr("descent.none")).weak());
            }
            None => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(egui::RichText::new(self.tr("descent.loading")).weak());
                });
            }
        }
    }

    /// Fire a `LoadDescentReport` command if this (subject, DNA) report isn't already loaded or in
    /// flight. Idempotent — safe to call every frame.
    fn ensure_descent(&mut self, guid: SampleGuid, dna: DnaType) {
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
