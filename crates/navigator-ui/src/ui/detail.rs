//! `impl NavigatorApp` methods extracted from `ui.rs` (the `detail` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    /// Y-STR profiles for the selected subject + an import form (CSV/TSV marker table).
    /// Donor-level Y-STR consensus across all of the subject's panels (Phase 2 rollup): the modal
    /// value per marker, with cross-panel disagreements flagged.
    /// Y-STR report (FTDNA/YSEQ style): summary header (provider toggle, tier badges, conflict
    /// count) + a By-Panel / All-Markers / Consensus view, rendered from the already-loaded
    /// `str_profiles`.
    /// Y-STR called from sequence (HipSTR caller → FTDNA convention) vs the imported vendor profile.
    pub(crate) fn ystr_sequence_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.horizontal(|ui| {
            let have = matches!(&self.str_concordance, Some((g, _, _)) if *g == guid);
            let label = if have { self.tr("common.refresh") } else { self.tr("ystr.callFromSequence") };
            if ui.add_enabled(!self.str_running, egui::Button::new(label)).clicked() {
                self.str_running = true;
                self.status = "Calling Y-STRs from sequence (first run scans chrY)…".into();
                let _ = self.tx.send(Command::StrConcordance { biosample_guid: guid });
            }
            if self.str_running {
                ui.spinner();
            }
            ui.label(egui::RichText::new(self.tr("hint.ystrSequence")).weak().small());
        });

        let have = matches!(&self.str_concordance, Some((g, _, _)) if *g == guid);
        if have {
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.str_seq_query).hint_text("filter marker").desired_width(140.0));
                if !self.str_seq_query.is_empty() && ui.small_button("✕").clicked() {
                    self.str_seq_query.clear();
                }
            });
        }
        let q = self.str_seq_query.to_ascii_lowercase();

        let Some((g, aln, rows)) = &self.str_concordance else { return };
        if *g != guid {
            return;
        }
        let called = rows.iter().filter(|r| r.called.is_some()).count();
        let agree = rows.iter().filter(|r| r.agree).count();
        let compared = rows.iter().filter(|r| r.called.is_some() && r.imported.is_some() && r.calibrated).count();
        ui.add_space(4.0);
        ui.label(
            egui::RichText::new(format!("aln #{aln}: {called} markers called · {agree}/{compared} calibrated agree with vendor"))
                .weak()
                .small(),
        );
        ui.add_space(4.0);
        // No inner ScrollArea — the tab is already one vertical scroll; nesting clips + captures the
        // wheel (see str_by_panel_view). Flow the (filtered) grid into the page; the filter narrows it.
        egui::Grid::new(("ystr_seq_grid", guid)).num_columns(4).striped(true).spacing([14.0, 2.0]).show(ui, |ui| {
            for h in ["Marker", "Called", "Vendor", ""] {
                ui.label(egui::RichText::new(h).strong().small());
            }
            ui.end_row();
            for r in rows
                .iter()
                .filter(|r| r.called.is_some() || r.imported.is_some())
                .filter(|r| q.is_empty() || r.marker.to_ascii_lowercase().contains(&q))
            {
                ui.label(&r.marker);
                // Colour the called value by calibration status.
                let (txt, col) = match (r.called, r.status.as_str()) {
                    (Some(v), "Reliable" | "ConventionOffset") => (v.to_string(), None),
                    (Some(v), _) => (v.to_string(), Some(egui::Color32::from_rgb(150, 150, 150))), // excluded/uncalibrated
                    (None, _) => ("—".to_string(), Some(egui::Color32::from_rgb(150, 150, 150))),
                };
                match col {
                    Some(c) => ui.colored_label(c, txt),
                    None => ui.label(txt),
                };
                ui.label(r.imported.clone().unwrap_or_else(|| "—".into()));
                // Agreement marker only for calibrated, comparable rows.
                if r.calibrated && r.called.is_some() && r.imported.is_some() {
                    if r.agree {
                        ui.colored_label(egui::Color32::from_rgb(60, 160, 60), "✓");
                    } else {
                        ui.colored_label(egui::Color32::from_rgb(200, 90, 90), "✗");
                    }
                } else {
                    ui.label("");
                }
                ui.end_row();
            }
        });
    }

    /// Cross-subject Y matches (gap §2): rank every other workspace subject by Y relatedness. Button-
    /// driven (reads cached profiles); flows into the page scroll + filter (no nested ScrollArea).
    pub(crate) fn ymatch_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        // Controls: project filter + find button (mirror the assign-project picker's local-copy idiom).
        let whole = self.tr("ymatch.wholeWorkspace").to_string();
        let projects: Vec<(i64, String)> =
            self.overview.iter().map(|o| (o.project.id, o.project.name.clone())).collect();
        let mut chosen = self.y_match_project;
        let sel_text = match chosen {
            Some(pid) => projects.iter().find(|(id, _)| *id == pid).map(|(_, n)| n.clone()).unwrap_or_else(|| format!("project {pid}")),
            None => whole.clone(),
        };
        let find_label = self.tr("ymatch.find").to_string();
        let hint = self.tr("hint.yMatches").to_string();
        let mut do_find = false;
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt(("ymatch_project", guid))
                .selected_text(sel_text)
                .width(220.0)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut chosen, None, &whole);
                    for (id, name) in &projects {
                        ui.selectable_value(&mut chosen, Some(*id), name);
                    }
                });
            if ui.add_enabled(!self.y_matches_running, egui::Button::new(find_label)).clicked() {
                do_find = true;
            }
            if self.y_matches_running {
                ui.spinner();
            }
            ui.label(egui::RichText::new(hint).weak().small());
        });
        self.y_match_project = chosen;
        if do_find {
            self.y_matches_running = true;
            self.status = "Finding Y matches across the workspace…".into();
            let _ = self.tx.send(Command::YMatches { biosample_guid: guid, project_id: chosen });
        }

        let have = matches!(&self.y_matches, Some((g, _)) if *g == guid);
        if have {
            ui.horizontal(|ui| {
                ui.add(egui::TextEdit::singleline(&mut self.y_match_query).hint_text("filter").desired_width(140.0));
                if !self.y_match_query.is_empty() && ui.small_button("✕").clicked() {
                    self.y_match_query.clear();
                }
            });
        }
        let q = self.y_match_query.to_ascii_lowercase();
        let caveat = self.tr("ymatch.tmrcaCaveat").to_string();

        let Some((g, matches)) = &self.y_matches else { return };
        if *g != guid {
            return;
        }
        if matches.is_empty() {
            ui.add_space(4.0);
            ui.label(egui::RichText::new(self.tr("ymatch.none")).weak());
            return;
        }
        ui.add_space(4.0);
        ui.label(egui::RichText::new(format!("{} matches", matches.len())).weak().small());
        ui.add_space(4.0);
        // No inner ScrollArea — the tab is already one vertical scroll (see str_by_panel_view).
        egui::Grid::new(("ymatch_grid", guid)).num_columns(7).striped(true).spacing([14.0, 2.0]).show(ui, |ui| {
            for h in ["Subject", "Shared", "Novel", "Divergence", "STR-GD", "Signal", "TMRCA"] {
                ui.label(egui::RichText::new(h).strong().small());
            }
            ui.end_row();
            for m in matches.iter().filter(|m| {
                q.is_empty()
                    || m.donor.to_ascii_lowercase().contains(&q)
                    || m.terminal.as_deref().is_some_and(|t| t.to_ascii_lowercase().contains(&q))
            }) {
                ui.label(&m.donor);
                let snp_backed = m.signal != YSignal::Str;
                ui.label(if snp_backed { m.shared_derived.to_string() } else { "—".into() });
                ui.label(if snp_backed { m.shared_novel.to_string() } else { "—".into() });
                ui.label(m.divergence.clone().unwrap_or_else(|| "—".into()));
                ui.label(match m.str_gd {
                    Some(gd) => format!("{gd} / {}", m.str_markers),
                    None => "—".into(),
                });
                let sig = match m.signal {
                    YSignal::SnpStr => "SNP+STR",
                    YSignal::Snp => "SNP",
                    YSignal::Str => "STR",
                    YSignal::None => "—",
                };
                ui.label(egui::RichText::new(sig).small());
                let tmrca = m.snp_tmrca.as_ref().or(m.str_tmrca.as_ref());
                match tmrca {
                    Some(t) => {
                        ui.label(format!("~{:.0} gen / ~{:.0} yr", t.generations, t.years)).on_hover_text(&caveat);
                    }
                    None => {
                        ui.label("—");
                    }
                }
                ui.end_row();
            }
        });
    }

    pub(crate) fn ystr_report_section(&mut self, ui: &mut egui::Ui) {
        if self.str_profiles.is_empty() {
            ui.label(egui::RichText::new("No STR profiles yet — import one under Data Sources.").weak());
            return;
        }

        let comparison = strprofile::compare_profiles(&self.str_profiles);
        let multi_provider = comparison.providers.len() > 1;

        // Working copies (written back after rendering) so row clicks can mutate selection while
        // `self.str_profiles` is borrowed immutably.
        let mut sel_provider = self.str_provider.clone().unwrap_or_else(|| {
            self.str_profiles
                .iter()
                .max_by_key(|p| p.markers.len())
                .and_then(|p| p.provider.clone())
                .unwrap_or_else(|| "FTDNA".to_string())
        });
        let mut view = self.str_report_view;
        let mut filter = std::mem::take(&mut self.str_marker_filter);

        // Provider toggle (only when >1 provider).
        if multi_provider {
            ui.horizontal(|ui| {
                ui.label("Provider:");
                for prov in &comparison.providers {
                    if ui.selectable_label(&sel_provider == prov, prov).clicked() {
                        sel_provider = prov.clone();
                    }
                }
            });
        }

        // The profile to display: the most-complete one for the selected provider.
        let canon = strpanel::canonical_provider(&sel_provider);
        let profile_idx = self
            .str_profiles
            .iter()
            .enumerate()
            .filter(|(_, p)| strpanel::canonical_provider(p.provider.as_deref().unwrap_or("FTDNA")) == canon)
            .max_by_key(|(_, p)| p.markers.len())
            .or_else(|| self.str_profiles.iter().enumerate().max_by_key(|(_, p)| p.markers.len()))
            .map(|(i, _)| i);
        let Some(idx) = profile_idx else { return };
        let marker_count = self.str_profiles[idx].markers.len();

        let amber = egui::Color32::from_rgb(220, 150, 60);
        ui.horizontal(|ui| {
            ui.heading(marker_count.to_string());
            ui.label(egui::RichText::new("markers").weak());
            if !comparison.conflicts.is_empty() {
                ui.colored_label(amber, format!("⚠ {} conflict(s)", comparison.conflicts.len()));
            }
        });

        // Tier "reached" badges.
        ui.horizontal_wrapped(|ui| {
            for (name, filled) in strpanel::tier_badges(&sel_provider, marker_count) {
                let (glyph, color) = if filled {
                    ("●", egui::Color32::from_rgb(120, 180, 120))
                } else {
                    ("○", egui::Color32::from_gray(110))
                };
                ui.colored_label(color, format!("{glyph} {name}"))
                    .on_hover_text(if filled { "reached" } else { "not reached" });
            }
        });

        // View selector.
        ui.horizontal(|ui| {
            ui.selectable_value(&mut view, StrReportView::ByPanel, "By panel");
            ui.selectable_value(&mut view, StrReportView::AllMarkers, "All markers");
            ui.selectable_value(&mut view, StrReportView::Consensus, "Consensus");
        });
        ui.separator();

        match view {
            StrReportView::ByPanel => {
                str_by_panel_view(ui, &self.str_profiles[idx], &sel_provider, &comparison);
            }
            StrReportView::AllMarkers => {
                str_all_markers_view(ui, &self.str_profiles[idx], &sel_provider, &comparison, &mut filter);
            }
            StrReportView::Consensus => self.str_consensus_section(ui),
        }

        self.str_report_view = view;
        self.str_provider = Some(sel_provider);
        self.str_marker_filter = filter;
    }

    fn str_consensus_section(&mut self, ui: &mut egui::Ui) {
        if self.str_profiles.is_empty() {
            ui.label(egui::RichText::new("No STR profiles yet — import one under Data Sources.").weak());
            return;
        }
        let consensus = strprofile::consensus_markers(&self.str_profiles);
        ui.label(
            egui::RichText::new(format!("{} markers from {} panel(s)", consensus.len(), self.str_profiles.len()))
                .weak(),
        );
        let conflicts = consensus.iter().filter(|m| m.conflict).count();
        if conflicts > 0 {
            ui.colored_label(
                egui::Color32::from_rgb(220, 150, 60),
                format!("⚠ {conflicts} marker(s) disagree across panels"),
            );
        }
        egui::Grid::new("str_consensus").striped(true).num_columns(3).show(ui, |ui| {
            ui.strong(self.tr("table.marker"));
            ui.strong(self.tr("table.value"));
            ui.strong(self.tr("table.panels"));
            ui.end_row();
            for m in &consensus {
                ui.label(&m.marker);
                if m.conflict {
                    ui.colored_label(egui::Color32::from_rgb(220, 150, 60), &m.value);
                } else {
                    ui.label(&m.value);
                }
                ui.label(m.panels.to_string());
                ui.end_row();
            }
        });
    }

    /// Donor-level ancestry headline (Phase 3): the best estimate across the subject's sources,
    /// with which source + method it came from.
    /// The donor's projected (PC1, PC2), from whichever loaded estimate carries PCA coordinates
    /// (the PCA / nMonte methods, or ADMIXTURE with PCA attached).
    fn sample_pca(&self) -> Option<(f64, f64)> {
        [
            self.donor_ancestry.as_ref().map(|(_, r)| r),
            self.ancient_ancestry.as_ref(),
            self.nmonte_ancestry.as_ref(),
        ]
        .into_iter()
        .flatten()
        .find_map(|r| {
            let c = r.pca_coordinates.as_ref()?;
            (c.len() >= 2).then(|| (c[0], c[1]))
        })
    }

    /// PCA scatter: the donor's PC1×PC2 against the reference population centroids. The donor's
    /// coordinate is always projected in the CHM13 consensus PCA frame, so the reference centroids are
    /// loaded from that same asset (once, guarded) — not the selected source's build, which would mix
    /// frames (or miss the asset entirely) and collapse the plot onto the lone donor point.
    pub(crate) fn pca_scatter_section(&mut self, ui: &mut egui::Ui) {
        let key = navigator_app::CONSENSUS_SOURCE_ID;
        let loaded = matches!(&self.pca_reference, Some((a, _)) if *a == key);
        if !loaded && self.pca_reference_attempted != Some(key) {
            self.pca_reference_attempted = Some(key);
            let _ = self.tx.send(Command::LoadPcaReference);
        }
        let reference: &[(String, f64, f64)] =
            self.pca_reference.as_ref().filter(|(a, _)| *a == key).map(|(_, r)| r.as_slice()).unwrap_or(&[]);
        // Don't render a degenerate one-point plot: without the reference cloud the scatter
        // auto-zooms onto the donor alone and is meaningless. Surface the missing asset instead.
        if reference.is_empty() {
            ui.label(egui::RichText::new(self.tr("pca.referenceMissing")).weak());
            return;
        }
        draw_pca_scatter(ui, self.sample_pca(), reference);
    }

    pub(crate) fn donor_ancestry_summary(&self, ui: &mut egui::Ui) {
        let Some((aln, r)) = &self.donor_ancestry else {
            ui.label(egui::RichText::new("No ancestry estimate for any source yet.").weak());
            return;
        };
        ui.horizontal(|ui| {
            draw_ancestry_donut(ui, &r.super_population_summary);
            ui.add_space(8.0);
            ui.vertical(|ui| {
                if let Some(top) = r.super_population_summary.first() {
                    ui.heading(format!("{} {:.1}%", top.super_population, top.percentage));
                }
                ui.label(format!(
                    "{}/{} SNPs · confidence {:.0}%",
                    r.snps_with_genotype,
                    r.snps_analyzed,
                    r.confidence_level * 100.0
                ));
                ui.label(
                    egui::RichText::new(format!("best source: alignment #{aln} · {} · {}", r.method, r.reference_version))
                        .small()
                        .weak(),
                );
                ui.add_space(4.0);
                draw_composition_bar(ui, &r.super_population_summary);
            });
        });
    }

    /// The build/refresh control shared by the Y/mt/autosomal consensus cards: a button that
    /// reads "Refresh" once a profile exists (else `build_label_key`), an inline spinner while
    /// `loading`, and a weak cost hint. On click it sets `status` and dispatches `command`;
    /// returns whether it was clicked so the caller can set its own loading flag.
    #[allow(clippy::too_many_arguments)]
    fn profile_build_control(
        &mut self,
        ui: &mut egui::Ui,
        has_profile: bool,
        loading: bool,
        build_label_key: &'static str,
        cost_hint_key: &'static str,
        status: &str,
        command: Command,
    ) -> bool {
        let mut clicked = false;
        ui.horizontal(|ui| {
            let label = if has_profile { self.tr("common.refresh") } else { self.tr(build_label_key) };
            if ui.add_enabled(!loading, egui::Button::new(label)).clicked() {
                self.status = status.into();
                let _ = self.tx.send(command);
                clicked = true;
            }
            if loading {
                ui.spinner();
            }
            ui.label(egui::RichText::new(self.tr(cost_hint_key)).weak().small());
        });
        clicked
    }

    /// Multi-source Y-variant profile: per-SNP concordance across the subject's Y sources, with
    /// status (confirmed/novel/conflict/single) and per-source provenance.
    pub(crate) fn y_variant_profile_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        // Build/refresh control first (it mutates self / dispatches), before borrowing the profile.
        if self.profile_build_control(
            ui,
            self.y_profile.is_some(),
            self.y_profile_loading,
            "btn.buildYProfile",
            "hint.yProfileCost",
            "Building Y variant profile…",
            Command::BuildYProfile { biosample_guid: guid },
        ) {
            self.y_profile_loading = true;
        }

        // Read-only source-audit (per-source provenance + per-conflict evidence) — opens a modal
        // over the cached profile; no re-genotyping.
        if self.y_profile.is_some() && ui.small_button(self.tr("audit.open")).clicked() {
            self.audit_y_profile = true;
        }

        let Some(profile) = &self.y_profile else {
            if !self.y_profile_loading {
                ui.label(egui::RichText::new(self.tr("hint.yProfileBuild")).weak());
            }
            return;
        };
        let mut filter = self.y_profile_filter;
        let mut query = std::mem::take(&mut self.y_profile_query);
        draw_consensus_profile(ui, profile, &mut filter, &mut query, "SNP", "Y variants", "y_variant_profile", &self.y_snp_names);
        self.y_profile_filter = filter;
        self.y_profile_query = query;
    }

    /// Multi-source mtDNA consensus profile: per-mutation concordance across the subject's mt
    /// sources (alignments' chrM placement, imported mtDNA sequences, the chip mt panel). Mirrors
    /// the Y-variant card over the same generic consensus engine.
    pub(crate) fn mt_variant_profile_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.profile_build_control(
            ui,
            self.mt_profile.is_some(),
            self.mt_profile_loading,
            "btn.buildMtProfile",
            "hint.mtProfileCost",
            "Building mtDNA consensus profile…",
            Command::BuildMtProfile { biosample_guid: guid },
        ) {
            self.mt_profile_loading = true;
        }

        let Some(profile) = &self.mt_profile else {
            if !self.mt_profile_loading {
                ui.label(egui::RichText::new(self.tr("hint.mtProfileBuild")).weak());
            }
            return;
        };
        let mut filter = self.mt_profile_filter;
        let mut query = std::mem::take(&mut self.mt_profile_query);
        // mtDNA mutations are already named (rCRS notation) — no Y-SNP catalogue annotation.
        draw_consensus_profile(ui, profile, &mut filter, &mut query, "Mutation", "mtDNA mutations", "mt_variant_profile", &std::collections::HashMap::new());
        self.mt_profile_filter = filter;
        self.mt_profile_query = query;
    }

    /// Multi-source autosomal (diploid 0/1/2) consensus over the canonical IBD-panel sites. Build/
    /// Refresh recomputes (panel-genotypes every WGS + chip source); the cached snapshot loads instantly.
    /// Autosomal-consensus **Summary** sub-tab: the build/refresh control plus a one-line digest
    /// (site count + overall confidence) once the profile exists. The heavy per-site table lives on
    /// the Profile sub-tab ([`autosomal_profile_table`]).
    pub(crate) fn autosomal_summary_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.profile_build_control(
            ui,
            self.auto_profile.is_some(),
            self.auto_profile_loading,
            "btn.buildAutosomalProfile",
            "hint.autoProfileCost",
            "Building autosomal consensus profile…",
            Command::BuildAutosomalProfile { biosample_guid: guid },
        ) {
            self.auto_profile_loading = true;
        }

        let Some(profile) = &self.auto_profile else {
            if !self.auto_profile_loading {
                ui.label(egui::RichText::new(self.tr("hint.autoProfileBuild")).weak());
            }
            return;
        };
        let s = &profile.summary;
        ui.label(format!(
            "{} sites · {:.0}% confidence",
            s.total,
            s.overall_confidence * 100.0
        ));
        ui.label(egui::RichText::new(self.tr("hint.autoProfileTable")).weak().small());
    }

    /// Autosomal-consensus **Profile** sub-tab: the full per-site reconciled table. Renders nothing
    /// until the profile has been built from the Summary sub-tab.
    pub(crate) fn autosomal_profile_table(&mut self, ui: &mut egui::Ui) {
        let Some(profile) = &self.auto_profile else {
            ui.label(egui::RichText::new(self.tr("hint.autoProfileBuild")).weak());
            return;
        };
        let mut filter = self.auto_profile_filter;
        let mut query = std::mem::take(&mut self.auto_profile_query);
        draw_diploid_profile(ui, profile, &mut filter, &mut query);
        self.auto_profile_filter = filter;
        self.auto_profile_query = query;
    }

    /// Consensus dashboard (Overview): the subject's source-of-truth at a glance — consensus Y/mt
    /// haplogroups, top ancestry, the autosomal concordance one-liner, and a source inventory.
    pub(crate) fn overview_dashboard(&mut self, ui: &mut egui::Ui, _guid: SampleGuid) {
        let none = self.tr("hint.noConsensusYet");
        card(ui, self.tr("card.consensusSummary"), |ui| {
            let line = |ui: &mut egui::Ui, label: &str, cons: &Option<navigator_app::Consensus>| {
                ui.horizontal(|ui| {
                    ui.strong(label);
                    match cons {
                        Some(c) => {
                            ui.label(&c.haplogroup);
                            ui.label(egui::RichText::new(format!("({} source(s), conf {:.2})", c.run_count, c.confidence)).weak().small());
                        }
                        None => {
                            ui.label(egui::RichText::new(none).weak());
                        }
                    }
                });
            };
            line(ui, "Y-DNA:", &self.consensus_y);
            line(ui, "mtDNA:", &self.consensus_mt);
            ui.horizontal(|ui| {
                ui.strong("Ancestry:");
                match self.donor_ancestry.as_ref().and_then(|(_, r)| r.super_population_summary.first()) {
                    Some(top) => {
                        ui.label(format!("{} {:.1}%", top.super_population, top.percentage));
                    }
                    None => {
                        ui.label(egui::RichText::new(none).weak());
                    }
                }
            });
            ui.horizontal(|ui| {
                ui.strong("Autosomal:");
                match &self.auto_profile {
                    Some(p) => {
                        ui.label(format!(
                            "{} sites · {} confirmed · {} conflict · {:.0}% confidence",
                            p.summary.total, p.summary.confirmed, p.summary.conflict, p.summary.overall_confidence * 100.0
                        ));
                    }
                    None => {
                        ui.label(egui::RichText::new(none).weak());
                    }
                }
            });
        });
        ui.add_space(10.0);
        card(ui, self.tr("card.sourceInventory"), |ui| {
            ui.label(format!("{} sequencing run(s)", self.runs.len()));
            ui.label(format!("{} chip/array profile(s)", self.chip_profiles.len()));
            ui.label(format!("{} STR profile(s)", self.str_profiles.len()));
            ui.label(format!("{} mtDNA sequence(s)", self.mtdna_sequences.len()));
            ui.add_space(4.0);
            ui.label(egui::RichText::new(self.tr("hint.sourcesTabHint")).weak().small());
        });
    }

    /// The per-sequencing-result hub: the source lists (runs/alignments/chips/STR/mtDNA) plus, for the
    /// selected source, the inherently-per-result views (coverage, sex/metrics/SV, ideogram,
    /// heteroplasmy, chrM de-novo) and that source's own Y/mt haplogroup placement.
    pub(crate) fn sources_tab(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        self.data_sources_tab(ui, guid);
        ui.add_space(10.0);
        ui.separator();
        let Some(id) = self.selected_alignment else {
            ui.add_space(4.0);
            ui.label(egui::RichText::new(self.tr("sources.pickHint")).weak());
            return;
        };
        ui.add_space(6.0);
        ui.strong(self.tr("sources.selectedDetail"));
        ui.add_space(4.0);
        if ui
            .add(egui::Button::new(egui::RichText::new(self.tr("action.runFullAnalysis")).color(egui::Color32::WHITE)).fill(ACCENT))
            .clicked()
        {
            self.start_full_analysis(id);
        }
        ui.add_space(8.0);
        card(ui, self.tr("card.coverage"), |ui| self.coverage_section(ui, id));
        ui.add_space(10.0);
        card(ui, self.tr("card.sexMetrics"), |ui| self.sex_metrics_section(ui, id));
        ui.add_space(10.0);
        card(ui, self.tr("card.yHaplogroup"), |ui| self.y_haplogroup_section(ui, id));
        ui.add_space(10.0);
        card(ui, self.tr("card.mtHaplogroup"), |ui| self.mt_haplogroup_section(ui, id));
        ui.add_space(10.0);
        card(ui, self.tr("card.mtDenovo"), |ui| self.denovo_section(ui, id, "chrM"));
        ui.add_space(10.0);
        card(ui, self.tr("card.mtHeteroplasmy"), |ui| self.heteroplasmy_section(ui, id));
    }

    /// Donor-level private-Y union (Phase 3): off-backbone calls pooled + deduped across the
    /// subject's Y-bearing sources.
    pub(crate) fn donor_private_y_section(&mut self, ui: &mut egui::Ui) {
        if self.donor_private_y.is_none() {
            ui.label(egui::RichText::new("No private-Y calls across sources yet — run \"Find private Y variants\".").weak());
            return;
        }
        let (pos_h, chg_h, dep_h, cls_h) =
            (self.tr("table.position"), self.tr("table.change"), self.tr("table.depth"), self.tr("table.class"));
        if let Some(b) = &self.donor_private_y {
            ui.label(format!("{} novel + {} off-path  (union across sources, terminal {})", b.novel(), b.off_path(), b.terminal));
        }
        ui.horizontal(|ui| {
            ui.add(egui::TextEdit::singleline(&mut self.private_y_query).hint_text("filter pos / name").desired_width(160.0));
            if !self.private_y_query.is_empty() && ui.small_button("✕").clicked() {
                self.private_y_query.clear();
            }
        });
        let q = self.private_y_query.to_ascii_lowercase();
        let bucket = self.donor_private_y.as_ref().unwrap();
        let names = &self.y_snp_names; // catalogued Y-SNP name at a novel call's site, if any
        // Filter to matching variants (position, off-path name, "novel", or the catalogued name); the
        // table is bounded to a fixed-height scroll pane (a WGS bucket runs to thousands of rows). A
        // hard cap keeps a pathological bucket from flooding even the pane.
        const CAP: usize = 1000;
        let matched: Vec<_> = bucket
            .variants
            .iter()
            .filter(|v| {
                q.is_empty()
                    || v.position.to_string().contains(&q)
                    || names.get(&v.position).is_some_and(|n| n.to_ascii_lowercase().contains(&q))
                    || match &v.class {
                        PrivateClass::OffPathKnown(n) => n.to_ascii_lowercase().contains(&q),
                        PrivateClass::Novel => "novel".contains(q.as_str()),
                    }
            })
            .collect();
        ui.label(egui::RichText::new(format!("{} shown", matched.len())).weak().small());
        let pane_h = profile_pane_height(ui, matched.len());
        egui::ScrollArea::vertical().id_salt("donor_privy_scroll").max_height(pane_h).auto_shrink([false, false]).show(ui, |ui| {
            egui::Grid::new("donor_privy").striped(true).num_columns(4).show(ui, |ui| {
                ui.strong(pos_h);
                ui.strong(chg_h);
                ui.strong(dep_h);
                ui.strong(cls_h);
                ui.end_row();
                let teal = egui::Color32::from_rgb(90, 190, 190);
                for v in matched.iter().take(CAP) {
                    ui.label(v.position.to_string());
                    ui.label(format!("{}>{}", v.reference, v.alternate));
                    ui.label(v.depth.to_string());
                    match &v.class {
                        // A "novel" call that lands on a catalogued Y-SNP: surface that name (it's not
                        // on the placed lineage, but it is a known site, not a brand-new variant).
                        PrivateClass::Novel => match names.get(&v.position) {
                            Some(name) => ui
                                .colored_label(teal, format!("novel · {name}"))
                                .on_hover_text("catalogued Y-SNP at this site (off the placed lineage)"),
                            None => ui.colored_label(egui::Color32::from_rgb(60, 160, 60), "novel"),
                        },
                        PrivateClass::OffPathKnown(name) => ui.label(format!("off-path: {name}")),
                    };
                    ui.end_row();
                }
            });
        });
        if matched.len() > CAP {
            ui.label(egui::RichText::new(format!("…and {} more — filter to narrow", matched.len() - CAP)).weak());
        }
    }

    fn str_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let mut want_delete: Option<DataDelete> = None;
        for p in &self.str_profiles {
            let provider = p.provider.as_deref().unwrap_or("—");
            let header = format!("{} — {} markers  ({provider})", p.panel_name, p.markers.len());
            egui::CollapsingHeader::new(header).id_salt(("str", p.id)).show(ui, |ui| {
                if ui.small_button(self.tr("delete.thisProfile")).clicked() {
                    want_delete = Some(DataDelete::Str { id: p.id, guid, label: format!("STR profile “{}”", p.panel_name) });
                }
                egui::Grid::new(("str_markers", p.id)).striped(true).num_columns(2).show(ui, |ui| {
                    ui.strong(self.tr("table.marker"));
                    ui.strong(self.tr("table.value"));
                    ui.end_row();
                    for m in &p.markers {
                        ui.label(&m.marker);
                        ui.label(&m.value);
                        ui.end_row();
                    }
                });
            });
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("str.import"), |ui| {
            // Bind labels first — `self.tr()` (immutable) can't share the statement with the
            // `&mut self.forms.*` below (the i18n borrow gotcha).
            let (panel_lbl, provider_lbl, source_lbl) =
                (self.tr("form.panel"), self.tr("form.provider"), self.tr("form.source"));
            combo(ui, panel_lbl, "str_panel", &mut self.forms.str_panel, strprofile::KNOWN_PANELS);
            combo(ui, provider_lbl, "str_provider", &mut self.forms.str_provider, strprofile::KNOWN_PROVIDERS);
            combo(ui, source_lbl, "str_source", &mut self.forms.str_source, strprofile::KNOWN_SOURCES);
            if ui.button(self.tr("str.chooseCsv")).clicked() {
                if let Some(path) = rfd::FileDialog::new().add_filter("STR table", &["csv", "tsv", "txt"]).pick_file() {
                    let _ = self.tx.send(Command::ImportStrProfile {
                        biosample_guid: guid,
                        panel_name: self.forms.str_panel.clone(),
                        provider: opt(&self.forms.str_provider),
                        source: opt(&self.forms.str_source),
                        path,
                    });
                }
            }
            ui.label(self.tr("str.expectsRows"));
        });
    }

    /// SNP variant sets for the selected subject + an import form (VCF or CSV/TSV).
    pub(crate) fn variants_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.variant_sets.is_empty() {
            ui.label(egui::RichText::new("No variants imported yet.").weak());
        }
        const MAX_ROWS: usize = 500;
        let mut want_delete: Option<DataDelete> = None;
        for s in &self.variant_sets {
            let build = s.reference_build.as_deref().map(|b| format!(" · {b}")).unwrap_or_default();
            let header = format!("{} — {} call(s){build}", s.source_label, s.calls.len());
            egui::CollapsingHeader::new(header).id_salt(("vset", s.id)).show(ui, |ui| {
                ui.horizontal(|ui| {
                    if ui.small_button(self.tr("delete.thisProfile")).clicked() {
                        want_delete = Some(DataDelete::Variant { id: s.id, guid, label: format!("variant set “{}”", s.source_label) });
                    }
                    // A Y-SNP panel (BISDNA): place a Y haplogroup from its derived calls.
                    let is_y_panel = s.source_type == SourceType::Chip
                        && s.calls.iter().any(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"));
                    if is_y_panel && ui.small_button(self.tr("ysnp.placeHaplogroup")).clicked() {
                        let _ = self.tx.send(Command::AssignYBisdna { biosample_guid: guid });
                    }
                });
                let pane_h = profile_pane_height(ui, s.calls.len().min(MAX_ROWS));
                egui::ScrollArea::vertical().id_salt(("vcalls_scroll", s.id)).max_height(pane_h).auto_shrink([false, false]).show(ui, |ui| {
                    egui::Grid::new(("vcalls", s.id)).striped(true).num_columns(4).show(ui, |ui| {
                        for h in ["table.position", "table.change", "table.rsid", "table.genotype"] {
                            ui.strong(self.tr(h));
                        }
                        ui.end_row();
                        for c in s.calls.iter().take(MAX_ROWS) {
                            ui.label(format!("{} {}", c.contig, c.position));
                            ui.label(variant_change(c));
                            ui.label(c.rs_id.as_deref().unwrap_or("—"));
                            ui.label(c.genotype.as_deref().unwrap_or("—"));
                            ui.end_row();
                        }
                    });
                });
                if s.calls.len() > MAX_ROWS {
                    ui.label(format!("…and {} more", s.calls.len() - MAX_ROWS));
                }
            });
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("variants.import"), |ui| {
            let labels: Vec<&str> = SourceType::ALL.iter().map(|t| t.as_str()).collect();
            let source_lbl = self.tr("form.source");
            combo(ui, source_lbl, "variant_source", &mut self.forms.variant_source_type, &labels);
            let source_type = SourceType::from_code(&self.forms.variant_source_type);

            if ui.button(self.tr("chip.import")).clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("variants", &["vcf", "csv", "tsv", "txt"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportVariants { biosample_guid: guid, path, source_type });
                }
            }
            ui.label(self.tr("chip.formatHint"));

            ui.separator();
            ui.label(self.tr("str.pasteCalls"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.variant_manual_label).hint_text("source label (e.g. YSEQ panel)"));
            ui.add(
                egui::TextEdit::multiline(&mut self.forms.variant_manual_text)
                    .hint_text("contig,position,ref,alt per line")
                    .desired_rows(3),
            );
            let ready = !self.forms.variant_manual_text.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new(self.tr("str.addPasted"))).clicked() {
                let label = opt(&self.forms.variant_manual_label).unwrap_or_else(|| source_type.as_str().to_string());
                let _ = self.tx.send(Command::AddVariants {
                    biosample_guid: guid,
                    source_label: label,
                    source_type,
                    text: self.forms.variant_manual_text.clone(),
                });
                self.forms.variant_manual_text.clear();
                self.forms.variant_manual_label.clear();
            }
        });
    }

    /// Genotyping-array (chip) profiles for the selected subject + an import form. The
    /// parser computes the QC summary and guesses the vendor; the dropdown can override it.
    fn chip_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let mut want_delete: Option<DataDelete> = None;
        for p in &self.chip_profiles {
            let s = &p.summary;
            let call_rate = if s.total_markers_possible > 0 {
                100.0 * s.total_markers_called as f64 / s.total_markers_possible as f64
            } else {
                0.0
            };
            let ver = p.chip_version.as_deref().map(|v| format!(" {v}")).unwrap_or_default();
            let header = format!("{}{ver} — {} markers, {:.1}% call rate", p.provider, s.total_markers_possible, call_rate);
            egui::CollapsingHeader::new(header).id_salt(("chip", p.id)).show(ui, |ui| {
                if ui.small_button(self.tr("delete.thisProfile")).clicked() {
                    want_delete = Some(DataDelete::Chip { id: p.id, guid, label: format!("chip profile ({})", p.provider) });
                }
                egui::Grid::new(("chip_qc", p.id)).striped(true).num_columns(2).show(ui, |ui| {
                    let row = |ui: &mut egui::Ui, k: &str, v: String| {
                        ui.label(k);
                        ui.label(v);
                        ui.end_row();
                    };
                    row(ui, "Markers possible", s.total_markers_possible.to_string());
                    row(ui, "Markers called", s.total_markers_called.to_string());
                    row(ui, "No-call rate", format!("{:.2}%", s.no_call_rate * 100.0));
                    row(ui, "Het rate (autosomal)", s.het_rate.map(|h| format!("{:.2}%", h * 100.0)).unwrap_or_else(|| "—".into()));
                    row(ui, "Autosomal called", s.autosomal_markers_called.to_string());
                    row(ui, "Y called", s.y_markers_called.to_string());
                    row(ui, "MT called", s.mt_markers_called.to_string());
                    if let Some(file) = &p.source_file_name {
                        row(ui, "Source file", file.clone());
                    }
                });
            });
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("chip.section"), |ui| {
            ui.horizontal(|ui| {
                ui.label(self.tr("form.provider"));
                egui::ComboBox::from_id_salt("chip_provider")
                    .selected_text(self.forms.chip_provider.clone())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.forms.chip_provider, AUTO_DETECT.to_string(), AUTO_DETECT);
                        for p in chipprofile::KNOWN_PROVIDERS {
                            ui.selectable_value(&mut self.forms.chip_provider, p.to_string(), *p);
                        }
                    });
            });
            if ui.button(self.tr("chip.chooseCsv")).clicked() {
                if let Some(path) = rfd::FileDialog::new().add_filter("array data", &["csv", "txt", "tsv"]).pick_file() {
                    let provider = (self.forms.chip_provider != AUTO_DETECT).then(|| self.forms.chip_provider.clone());
                    let _ = self.tx.send(Command::ImportChipProfile { biosample_guid: guid, provider, path });
                }
            }
            ui.label(self.tr("chip.rawHint"));
        });
    }

    /// mtDNA FASTA sequences for the selected subject + an import form. Per sequence: place the
    /// haplogroup, or show its rCRS-relative mutation list (derived against the bundled rCRS).
    pub(crate) fn mtdna_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.mtdna_sequences.is_empty() {
            ui.label(egui::RichText::new("No mtDNA sequences yet.").weak());
        }

        // Bind before the &self loop borrow — used inside the per-row closures.
        let assign_lbl = self.tr("common.assignHaplogroup");
        let mutations_lbl = self.tr("btn.showMtMutations");
        let delete_lbl = self.tr("common.delete");
        let export_lbl = self.tr("mt.exportVariants");
        let mut want_delete: Option<DataDelete> = None;
        // Status to apply after the &self loop (the loop borrow blocks &mut self).
        let mut export_status: Option<String> = None;
        for m in &self.mtdna_sequences {
            let name = m.source_file_name.as_deref().or(m.defline.as_deref()).unwrap_or("mtDNA");
            ui.horizontal(|ui| {
                ui.label(format!("{name} — {} bp, {} N", m.length(), m.n_count));
                if ui.button(mutations_lbl).clicked() {
                    let _ = self.tx.send(Command::LoadMtdnaVariants { mtdna_id: m.id });
                }
                if ui.button(assign_lbl).clicked() {
                    self.status = "Assigning haplogroup (fetching FTDNA tree)…".into();
                    let _ = self.tx.send(Command::AssignMtdnaHaplogroup { mtdna_id: m.id });
                }
                if ui.button(delete_lbl).clicked() {
                    want_delete = Some(DataDelete::Mtdna { id: m.id, guid, label: format!("mtDNA sequence “{name}”") });
                }
            });
            // Show the haplogroup result for this sequence, if any.
            if let Some((id, assignment)) = &self.mtdna_haplogroup {
                if *id == m.id {
                    show_assignment(ui, assignment);
                }
            }
            // Show the rCRS mutation list for this sequence, if loaded.
            if let Some(variants) = self.mtdna_variants.get(&m.id) {
                mtdna_mutations_view(ui, m.id, variants);
                if !variants.is_empty() && ui.button(export_lbl).clicked() {
                    let req = navigator_app::ExportRequest::MtdnaTsv(m.id);
                    if let Some(path) = rfd::FileDialog::new()
                        .set_file_name(req.default_filename())
                        .add_filter("TSV", &["tsv"])
                        .save_file()
                    {
                        let _ = self.tx.send(Command::Export { request: req, path });
                        export_status = Some(format!("Exporting {}…", req.label()));
                    }
                }
            }
        }
        if let Some(s) = export_status {
            self.status = s;
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("mt.importFasta"), |ui| {
            if ui.button(self.tr("mt.chooseFasta")).clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("FASTA", &["fa", "fasta", "fna", "fas"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportMtdna { biosample_guid: guid, path });
                }
            }
            ui.label(self.tr("mt.fullSeq"));
        });
    }

    pub(crate) fn panels_section(&mut self, ui: &mut egui::Ui) {
        ui.label(self.tr("table.panels"));
        let mut pick = None;
        for info in &self.panels {
            let label = format!("{}  ({} sites)", info.panel.name, info.site_count);
            if ui.selectable_label(self.selected_panel == Some(info.panel.id), label).clicked() {
                pick = Some(info.panel.id);
            }
        }
        if let Some(id) = pick {
            self.select_panel(id);
        }
        ui.add(egui::TextEdit::singleline(&mut self.forms.panel_import_name).hint_text("new panel name"));
        if ui
            .add_enabled(!self.forms.panel_import_name.trim().is_empty(), egui::Button::new(self.tr("mt.importSitesVcf")))
            .clicked()
        {
            if let Some(path) = rfd::FileDialog::new().add_filter("VCF", &["vcf"]).pick_file() {
                let _ = self.tx.send(Command::ImportPanel {
                    name: self.forms.panel_import_name.trim().to_string(),
                    path,
                });
                self.forms.panel_import_name.clear();
            }
        }
    }

    /// When an import is blocked on uncached reference builds, prompt to download them (with
    /// a progress bar); on completion the import auto-retries (see the `ReferenceReady` event).
    pub(crate) fn reference_prompt(&mut self, ui: &mut egui::Ui) {
        if self.reference_needs.is_empty() && self.reference_progress.is_none() {
            return;
        }
        ui.add_space(6.0);
        egui::Frame::group(ui.style()).show(ui, |ui| {
            if !self.reference_needs.is_empty() {
                ui.label(self.tr("refdl.required"));
                for b in &self.reference_needs {
                    ui.label(format!("  • {} (~{} MB)", b.build, b.est_bytes / 1_000_000));
                }
                if ui
                    .add_enabled(self.reference_progress.is_none(), egui::Button::new(self.tr("common.downloadContinue")))
                    .clicked()
                {
                    for build in self.reference_needs.iter().map(|b| b.build.clone()).collect::<Vec<_>>() {
                        let _ = self.tx.send(Command::ResolveReference { build });
                    }
                    self.status = "Downloading reference(s)…".into();
                }
            }
            if let Some((build, received, total)) = self.reference_progress.clone() {
                let text = match total {
                    Some(t) => format!("{build}: {} / {} MB", received / 1_000_000, t / 1_000_000),
                    None => format!("{build}: {} MB", received / 1_000_000),
                };
                match total {
                    Some(t) if t > 0 => {
                        ui.add(egui::ProgressBar::new(received as f32 / t as f32).text(text));
                    }
                    _ => {
                        ui.horizontal(|ui| {
                            ui.spinner();
                            ui.label(text);
                        });
                    }
                }
            }
        });
    }

    /// A per-sample coverage/haplogroup table for the open project, with per-row coverage
    /// recompute and a CSV export. Coverage/haplogroup cells show "—" until computed.
    pub(crate) fn project_report_section(&mut self, ui: &mut egui::Ui) {
        if self.project_report.is_empty() {
            return;
        }
        ui.add_space(12.0);
        ui.separator();
        ui.horizontal(|ui| {
            ui.heading(self.tr("projects.report"));
            let busy = self.analyzing || self.running;
            if ui.add_enabled(!busy, egui::Button::new(self.tr("projects.analyzeAll"))).clicked() {
                if let Some(pid) = self.selected_project {
                    self.analyzing = true;
                    self.deep_progress = None;
                    self.status = "Deep-analyzing project (per sample; runs in the background)…".into();
                    let _ = self.tx.send(Command::DeepAnalyzeProject(pid));
                }
            }
            if self.analyzing && ui.button(self.tr("common.cancel")).clicked() {
                let _ = self.tx.send(Command::CancelAnalysis);
                self.status = "Cancelling deep analysis…".into();
            }
            if ui.button(self.tr("projects.exportCsv")).clicked() {
                let csv = navigator_app::report_csv(&self.project_report);
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("CSV", &["csv"])
                    .set_file_name("project_report.csv")
                    .save_file()
                {
                    self.status = match std::fs::write(&path, csv) {
                        Ok(()) => format!("Wrote {}", path.display()),
                        Err(e) => format!("CSV write failed: {e}"),
                    };
                }
            }
        });

        // Streaming deep-analyze progress (sample N of M, current donor).
        if let Some((done, total, sample, fraction)) = self.deep_progress.clone() {
            ui.add(
                egui::ProgressBar::new(fraction)
                    .text(format!("Analyzing {}/{} — {sample}", done + 1, total)),
            );
        }

        let running = self.running || self.analyzing;
        let mut recompute: Option<i64> = None;
        let mut assign_y: Option<i64> = None;
        egui::Grid::new("project_report_grid").striped(true).num_columns(15).show(ui, |ui| {
            for h in ["report.sample", "report.alns", "report.meanCov", "report.median", "report.cov10x", "report.cov20x", "report.callable", "report.y", "report.mtdna", "report.sex", "report.readLen", "report.pctAln", "report.insert", "report.sv", "report.actions"] {
                ui.strong(self.tr(h));
            }
            ui.end_row();
            for r in &self.project_report {
                ui.label(&r.biosample.donor_identifier);
                ui.label(r.alignment_count.to_string());
                // Mean coverage, with a "lite" badge when it's a partial sidecar estimate that a
                // deep walk (the per-row coverage button) would upgrade.
                if r.coverage_partial {
                    ui.horizontal(|ui| {
                        ui.label(fmt_depth(r.mean_coverage));
                        ui.add(egui::Label::new(
                            egui::RichText::new("lite").small().color(egui::Color32::from_rgb(180, 140, 40)),
                        ))
                        .on_hover_text("Lite coverage from the pipeline sidecar — run coverage to compute the full per-base distribution.");
                    });
                } else {
                    ui.label(fmt_depth(r.mean_coverage));
                }
                ui.label(fmt_depth(r.median_coverage));
                ui.label(fmt_pct(r.pct_10x));
                ui.label(fmt_pct(r.pct_20x));
                ui.label(r.callable_bases.map(|v| v.to_string()).unwrap_or_else(|| "—".into()));
                ui.label(r.y_haplogroup.clone().unwrap_or_else(|| "—".into()));
                ui.label(r.mt_haplogroup.clone().unwrap_or_else(|| "—".into()));
                ui.label(r.sex.clone().unwrap_or_else(|| "—".into()));
                ui.label(fmt_depth(r.mean_read_length));
                ui.label(fmt_pct(r.pct_aligned));
                ui.label(fmt_depth(r.median_insert_size));
                ui.label(r.sv_count.map(|v| v.to_string()).unwrap_or_else(|| "—".into()));
                if let Some(aln) = r.primary_alignment_id {
                    ui.horizontal(|ui| {
                        if ui.add_enabled(!running, egui::Button::new(self.tr("btn.cov"))).clicked() {
                            recompute = Some(aln);
                        }
                        if ui.add_enabled(!running, egui::Button::new(self.tr("report.y"))).clicked() {
                            assign_y = Some(aln);
                        }
                    });
                } else {
                    ui.label("—");
                }
                ui.end_row();
            }
        });
        if let Some(aln) = recompute {
            self.running = true;
            self.status = "Recomputing coverage…".into();
            let _ = self.tx.send(Command::RunCoverage(aln));
        }
        if let Some(aln) = assign_y {
            self.status = "Assigning Y haplogroup…".into();
            let _ = self.tx.send(Command::AssignYHaplogroup { alignment_id: aln });
        }
    }

    pub(crate) fn samples_section(&mut self, ui: &mut egui::Ui) {
        let pid = self.selected_project.unwrap();
        let name = self
            .overview
            .iter()
            .find(|o| o.project.id == pid)
            .map(|o| o.project.name.clone())
            .unwrap_or_else(|| "project".into());
        ui.heading(format!("Samples — {name}"));
        ui.separator();
        if self.samples.is_empty() {
            ui.label(self.tr("projects.noSamples"));
        }
        let mut pick = None;
        for s in &self.samples {
            let label = format!(
                "{}  ({}, {})",
                s.donor_identifier,
                s.sample_accession.as_deref().unwrap_or("—"),
                s.sex.as_deref().unwrap_or("—"),
            );
            if ui.selectable_label(self.selected_sample == Some(s.guid), label).clicked() {
                pick = Some(s.guid);
            }
        }
        if let Some(guid) = pick {
            self.select_sample(guid);
        }
        ui.label(self.tr("projects.addSubjectsHint"));
    }

    /// The Data Sources tab: sequencing runs (cards with expandable alignments), chip/array,
    /// and STR profiles — each in a rounded card.
    fn data_sources_tab(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.add_space(4.0);
        card(ui, self.tr("card.sequencingRuns"), |ui| self.runs_card(ui, guid));
        ui.add_space(10.0);
        card(ui, self.tr("card.chipProfiles"), |ui| {
            if self.chip_profiles.is_empty() {
                ui.label(egui::RichText::new("No chip/array data").weak());
            }
            self.chip_section(ui, guid);
        });
        ui.add_space(10.0);
        card(ui, self.tr("card.strProfiles"), |ui| {
            if self.str_profiles.is_empty() {
                ui.label(egui::RichText::new("No STR profiles").weak());
            }
            self.str_section(ui, guid);
        });
    }

    /// The sequencing-runs body: one card per run (provider chip, title, read meta, Y/mt
    /// badges); the selected run expands to its alignment rows + the add-alignment form.
    fn runs_card(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.runs.is_empty() {
            ui.label(egui::RichText::new("No sequencing runs yet.").weak());
        }
        // Clone the small lists so we can call &mut self methods (add forms) inside the loop.
        let runs = self.runs.clone();
        let alignments = self.alignments.clone();
        let coverage_by_aln = self.coverage_by_aln.clone();
        let mut pick_run = None;
        let mut pick_aln = None;
        let mut want_delete: Option<DataDelete> = None;
        let mut want_edit_run: Option<EditRun> = None;
        let mut want_edit_aln: Option<EditAlignment> = None;
        let mut want_merge: Option<i64> = None;

        for r in &runs {
            let selected = self.selected_run == Some(r.id);
            let frame = egui::Frame::group(ui.style())
                .fill(if selected { ACCENT.gamma_multiply(0.18) } else { ui.visuals().extreme_bg_color })
                .stroke(if selected { egui::Stroke::new(1.0, ACCENT) } else { egui::Stroke::NONE })
                .rounding(egui::Rounding::same(6.0))
                .inner_margin(egui::Margin::same(10.0));
            let inner = frame.show(ui, |ui| {
                let mut edit_btn: Option<egui::Response> = None;
                let mut del_btn: Option<egui::Response> = None;
                let mut merge_btn: Option<egui::Response> = None;
                ui.horizontal(|ui| {
                    chip(ui, &provider_abbrev(&r.platform_name), ACCENT.gamma_multiply(0.3), ACCENT);
                    // Lab chip (FGC/FTDNA/YSEQ/Dante/Nebula…) when the sequencing facility is known.
                    if let Some(lab) = r.sequencing_facility.as_deref().filter(|s| !s.is_empty()) {
                        let abbr = navigator_domain::labs::abbreviation(lab, 6);
                        chip(ui, &abbr, egui::Color32::from_rgb(40, 70, 55), egui::Color32::from_rgb(150, 220, 180))
                            .on_hover_text(format!("Sequencing lab: {}", navigator_domain::labs::display_name(lab)));
                    }
                    ui.add_space(4.0);
                    let tt = testtype::by_code(&r.test_type);
                    ui.vertical(|ui| {
                        let plat = if r.platform_name.is_empty() { "—" } else { r.platform_name.as_str() };
                        let title = format!(
                            "{}  ·  {}  ·  {}",
                            testtype::display_name(&r.test_type),
                            plat,
                            r.instrument_model.as_deref().unwrap_or("—")
                        );
                        ui.label(egui::RichText::new(title).strong());
                        // Instrument serial (the lab crowd-source key) + flowcell, when inferred.
                        let inst = match (r.instrument_id.as_deref(), r.flowcell_id.as_deref()) {
                            (Some(i), Some(f)) => format!("   Instr: {i} · FC: {f}"),
                            (Some(i), None) => format!("   Instr: {i}"),
                            _ => String::new(),
                        };
                        // Library-level metrics: total reads + read/insert length (reads *aligned*
                        // is a per-alignment stat, shown on the alignment row, not here).
                        let read_len = r.mean_read_length.map(|v| format!("{v:.0} bp")).unwrap_or_else(|| "—".into());
                        let insert = r.mean_insert_size.map(|v| format!("{v:.0} bp")).unwrap_or_else(|| "—".into());
                        ui.label(
                            egui::RichText::new(format!(
                                "Reads: {}   Read len: {}   Insert: {}   {}{}",
                                fmt_reads(r.total_reads),
                                read_len,
                                insert,
                                r.library_layout.as_deref().unwrap_or("—"),
                                inst,
                            ))
                            .weak()
                            .small(),
                        );
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        del_btn = Some(ui.small_button("🗑").on_hover_text("Delete run + its alignments"));
                        edit_btn = Some(ui.small_button("✏").on_hover_text("Edit run"));
                        if runs.len() > 1 {
                            merge_btn = Some(ui.small_button("⤵").on_hover_text("Merge this run into another"));
                        }
                        if let Some(t) = tt {
                            let mt = matches!(t.target, testtype::TargetType::WholeGenome | testtype::TargetType::MtDna);
                            let y = matches!(t.target, testtype::TargetType::WholeGenome | testtype::TargetType::YChromosome);
                            if mt {
                                chip(ui, "mt", egui::Color32::from_rgb(70, 60, 90), egui::Color32::from_rgb(200, 180, 230));
                            }
                            if y {
                                chip(ui, "Y", egui::Color32::from_rgb(40, 70, 55), egui::Color32::from_rgb(150, 220, 180));
                            }
                        }
                    });
                });
                (edit_btn, del_btn, merge_btn)
            });
            let (edit_btn, del_btn, merge_btn) = inner.inner;
            // Row selection is sensed on the whole frame, which can swallow the inner buttons'
            // clicks; treat a button as hit when it was clicked OR the row swallowed the click
            // while the pointer was over it.
            let row_clicked = inner.response.interact(egui::Sense::click()).clicked();
            let hit = |b: &Option<egui::Response>| {
                b.as_ref().is_some_and(|r| r.clicked() || (row_clicked && r.contains_pointer()))
            };
            if hit(&edit_btn) {
                want_edit_run = Some(EditRun {
                    id: r.id,
                    guid,
                    test_type: r.test_type.clone(),
                    platform_name: r.platform_name.clone(),
                    instrument_model: r.instrument_model.clone().unwrap_or_default(),
                    library_layout: r.library_layout.clone().unwrap_or_default(),
                    sequencing_facility: r.sequencing_facility.clone().unwrap_or_default(),
                });
            } else if hit(&del_btn) {
                want_delete = Some(DataDelete::Run {
                    id: r.id,
                    guid,
                    label: format!("run “{}”", testtype::display_name(&r.test_type)),
                });
            } else if hit(&merge_btn) {
                want_merge = Some(r.id);
            } else if row_clicked {
                pick_run = Some(r.id);
            }

            // Selected run → its alignment rows + the add-alignment form.
            if selected {
                ui.indent(("alns", r.id), |ui| {
                    for a in &alignments {
                        let asel = self.selected_alignment == Some(a.id);
                        let (cov_s, call_s) = match coverage_by_aln.get(&a.id) {
                            Some(c) => (format!("{:.1}", c.mean_coverage), c.callable_bases.to_string()),
                            None => ("–".to_string(), "–".to_string()),
                        };
                        let row = egui::Frame::group(ui.style())
                            .fill(if asel { ACCENT.gamma_multiply(0.14) } else { ui.visuals().widgets.noninteractive.bg_fill })
                            .rounding(egui::Rounding::same(6.0))
                            .inner_margin(egui::Margin::symmetric(10.0, 8.0))
                            .show(ui, |ui| {
                                let mut edit_btn: Option<egui::Response> = None;
                                let mut del_btn: Option<egui::Response> = None;
                                ui.horizontal(|ui| {
                                    ui.label(egui::RichText::new(&a.reference_build).color(ACCENT).strong());
                                    ui.label(egui::RichText::new(if a.bam_path.is_some() { a.aligner.as_str() } else { "Unknown" }).weak());
                                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                        del_btn = Some(ui.small_button("🗑").on_hover_text("Delete alignment"));
                                        edit_btn = Some(ui.small_button("✏").on_hover_text("Edit alignment"));
                                        ui.add_space(10.0);
                                        ui.label(egui::RichText::new(format!("Callable: {call_s}")).weak().small());
                                        ui.add_space(10.0);
                                        ui.label(egui::RichText::new(format!("Coverage: {cov_s}")).weak().small());
                                    });
                                });
                                (edit_btn, del_btn)
                            });
                        let (edit_btn, del_btn) = row.inner;
                        let row_clicked = row.response.interact(egui::Sense::click()).clicked();
                        let hit = |b: &Option<egui::Response>| {
                            b.as_ref().is_some_and(|r| r.clicked() || (row_clicked && r.contains_pointer()))
                        };
                        if hit(&edit_btn) {
                            want_edit_aln = Some(EditAlignment {
                                id: a.id,
                                run_id: r.id,
                                reference_build: a.reference_build.clone(),
                                aligner: a.aligner.clone(),
                                variant_caller: a.variant_caller.clone().unwrap_or_default(),
                            });
                        } else if hit(&del_btn) {
                            want_delete = Some(DataDelete::Alignment {
                                id: a.id,
                                run_id: r.id,
                                label: format!("alignment {} ({})", a.id, a.reference_build),
                            });
                        } else if row_clicked {
                            pick_aln = Some(a.id);
                        }
                    }
                    ui.add_space(4.0);
                    self.add_alignment_form(ui, r.id);
                });
            }
            ui.add_space(6.0);
        }

        if let Some(id) = pick_run {
            self.select_run(id);
        }
        if let Some(id) = pick_aln {
            self.select_alignment(id);
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }
        if want_edit_run.is_some() {
            self.edit_run = want_edit_run;
        }
        if want_edit_aln.is_some() {
            self.edit_alignment = want_edit_aln;
        }
        if let Some(secondary) = want_merge {
            // Default the target (primary) to the first other run; the modal lets the user change it.
            let primary = runs.iter().map(|r| r.id).find(|&id| id != secondary);
            self.merge_runs = Some(MergeRuns { guid, secondary, primary });
        }
        self.add_test_form(ui, guid);
    }

    /// The "Add test" (sequencing run) form.
    fn add_test_form(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.collapsing(self.tr("run.addTest"), |ui| {
            ui.horizontal(|ui| {
                ui.label(self.tr("form.testType"));
                let current = testtype::display_name(&self.forms.run_test_type).to_string();
                egui::ComboBox::from_id_salt("test_type").selected_text(current).show_ui(ui, |ui| {
                    for t in testtype::CATALOG {
                        ui.selectable_value(
                            &mut self.forms.run_test_type,
                            t.code.to_string(),
                            format!("{}  ·  {}", t.display_name, t.target.label()),
                        );
                    }
                });
            });
            ui.add(egui::TextEdit::singleline(&mut self.forms.run_platform).hint_text("platform (optional, e.g. ILLUMINA)"));
            let ready = testtype::by_code(&self.forms.run_test_type).is_some();
            if ui.add_enabled(ready, egui::Button::new(self.tr("run.addTest"))).clicked() {
                let platform = opt(&self.forms.run_platform).unwrap_or_else(|| "UNKNOWN".into());
                let _ = self.tx.send(Command::AddRun(NewSequenceRun {
                    biosample_guid: guid,
                    platform_name: platform,
                    instrument_model: None,
                    test_type: self.forms.run_test_type.clone(),
                    library_layout: None,
                    total_reads: None,
                    pf_reads_aligned: None,
                    mean_read_length: None,
                    mean_insert_size: None,
                }));
                self.forms.run_platform.clear();
            }
        });
    }

    /// The "Add alignment" form for a run. Picking a BAM/CRAM probes its header to auto-fill the
    /// reference build + aligner; the reference FASTA is never asked for (resolved from the build).
    fn add_alignment_form(&mut self, ui: &mut egui::Ui, run_id: i64) {
        ui.collapsing(self.tr("aln.add"), |ui| {
            ui.horizontal(|ui| {
                if ui.button(self.tr("common.pickBamCram")).clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("alignment", &["bam", "cram"]).pick_file() {
                        self.forms.aln_bam = p.to_string_lossy().into_owned();
                        // Probe the header to auto-fill build + aligner.
                        let _ = self.tx.send(Command::ProbeAlignment { path: p });
                        self.status = "Reading header…".into();
                    }
                }
                ui.label(if self.forms.aln_bam.is_empty() { "—" } else { self.forms.aln_bam.as_str() });
            });
            ui.add(
                egui::TextEdit::singleline(&mut self.forms.aln_reference_build)
                    .hint_text("reference build (auto-detected; editable)"),
            );
            ui.add(
                egui::TextEdit::singleline(&mut self.forms.aln_aligner).hint_text("aligner (auto-detected; editable)"),
            );
            ui.label(egui::RichText::new("Reference FASTA is resolved from the build automatically.").weak().small());
            let ready = !self.forms.aln_reference_build.trim().is_empty()
                && !self.forms.aln_aligner.trim().is_empty()
                && !self.forms.aln_bam.is_empty();
            if ui.add_enabled(ready, egui::Button::new(self.tr("aln.add"))).clicked() {
                let _ = self.tx.send(Command::AddAlignment(NewAlignment {
                    sequence_run_id: run_id,
                    reference_build: self.forms.aln_reference_build.trim().to_string(),
                    aligner: self.forms.aln_aligner.trim().to_string(),
                    variant_caller: None,
                    bam_path: opt(&self.forms.aln_bam),
                    reference_path: None, // resolved on demand from the build
                    content_sha256: None,
                }));
                self.forms.aln_reference_build.clear();
                self.forms.aln_aligner.clear();
                self.forms.aln_bam.clear();
            }
        });
    }

    /// Map a consensus variant's state to a tick color + hover label for the variant track.
    fn variant_mark(v: &navigator_app::YProfileVariant) -> Option<VariantMark> {
        use navigator_domain::consensus::ConsensusState;
        let (color, state) = match (&v.consensus, v.in_tree) {
            (ConsensusState::Derived, true) => (egui::Color32::from_rgb(90, 180, 110), "in-tree derived"),
            (ConsensusState::Derived, false) => (egui::Color32::from_rgb(210, 150, 60), "novel / private"),
            (ConsensusState::Ancestral, _) => (egui::Color32::from_gray(120), "ancestral"),
            (ConsensusState::NoCall, _) => return None,
        };
        Some(VariantMark { name: v.name.clone(), position: v.position, color, state })
    }

    /// Lazily resolve catalogued Y-SNP names for the two Y-SNP tables' position-only / novel calls.
    /// Gathers every variant position from the Y consensus profile + the private-Y union and asks the
    /// worker for `position → name` once per subject (re-armed when either source reloads). No-op until
    /// at least one source is present.
    pub(crate) fn ensure_y_snp_names(&mut self, guid: SampleGuid) {
        if self.y_snp_names_requested {
            return;
        }
        let mut positions: Vec<i64> = Vec::new();
        if let Some(p) = &self.y_profile {
            positions.extend(p.variants.iter().filter(|v| v.name.is_empty()).map(|v| v.position));
        }
        if let Some(b) = &self.donor_private_y {
            positions.extend(b.variants.iter().map(|v| v.position));
        }
        if positions.is_empty() {
            return; // nothing to annotate yet (sources not loaded)
        }
        positions.sort_unstable();
        positions.dedup();
        self.y_snp_names_requested = true;
        let _ = self.tx.send(Command::LoadYSnpNames { biosample_guid: guid, positions });
    }

    /// chrY **variant track**: the Y consensus profile's called variants plotted along chromosome Y,
    /// PAR/centromere regions shaded from the selected alignment's genome regions when loaded.
    /// Replaces the genome-wide karyotype ideogram on the Y-DNA → SNP variants sub-tab.
    pub(crate) fn y_variant_track(&mut self, ui: &mut egui::Ui) {
        let Some(profile) = &self.y_profile else {
            ui.label(egui::RichText::new(self.tr("hint.yProfileBuild")).weak());
            return;
        };
        let marks: Vec<VariantMark> = profile.variants.iter().filter_map(Self::variant_mark).collect();
        if marks.is_empty() {
            ui.label(egui::RichText::new(self.tr("hint.noVariantsTrack")).weak());
            return;
        }

        // chrY length + PAR shading from the selected alignment's genome regions (lazily fetched).
        let (mut length, mut regions): (i64, Vec<TrackRegion>) = (62_460_029, Vec::new()); // CHM13 chrY fallback
        if let Some(id) = self.selected_alignment {
            if let Some(build) = self.alignments.iter().find(|a| a.id == id).map(|a| a.reference_build.clone()) {
                let loaded = matches!(&self.genome_regions, Some((aid, _)) if *aid == id);
                if !loaded && self.regions_attempted != Some(id) {
                    self.regions_attempted = Some(id);
                    self.loading_regions = true;
                    let _ = self.tx.send(Command::LoadGenomeRegions { alignment_id: id, build });
                }
            }
            if let Some((aid, gr)) = &self.genome_regions {
                if *aid == id {
                    if let Some(chr_y) = gr.chromosomes.get("chrY").or_else(|| gr.chromosomes.get("Y")) {
                        if chr_y.length > 0 {
                            length = chr_y.length;
                        }
                        for (s, e) in &chr_y.par {
                            regions.push(TrackRegion { start: *s, end: *e, color: egui::Color32::from_rgb(40, 60, 90), label: "PAR".into() });
                        }
                        if let Some((s, e)) = chr_y.centromere {
                            regions.push(TrackRegion { start: s, end: e, color: egui::Color32::from_rgb(90, 45, 45), label: "centromere".into() });
                        }
                    }
                }
            }
        }
        // Guard against build-mismatched positions overrunning the bar.
        length = length.max(marks.iter().map(|m| m.position).max().unwrap_or(0) + 1);
        draw_variant_track(ui, "chrY", length, &regions, &marks);
    }

    /// chrM **variant track**: the mtDNA consensus profile's mutations plotted along the 16,569 bp
    /// mitochondrial genome, with HVR1/HVR2 control regions shaded. On the mtDNA → Variants sub-tab.
    pub(crate) fn mt_variant_track(&mut self, ui: &mut egui::Ui) {
        let Some(profile) = &self.mt_profile else {
            ui.label(egui::RichText::new(self.tr("hint.mtProfileBuild")).weak());
            return;
        };
        let marks: Vec<VariantMark> = profile.variants.iter().filter_map(Self::variant_mark).collect();
        if marks.is_empty() {
            ui.label(egui::RichText::new(self.tr("hint.noVariantsTrack")).weak());
            return;
        }
        let regions = vec![
            TrackRegion { start: 16_024, end: 16_569, color: egui::Color32::from_rgb(50, 70, 50), label: "HVR1".into() },
            TrackRegion { start: 1, end: 576, color: egui::Color32::from_rgb(50, 70, 50), label: "HVR2".into() },
        ];
        draw_variant_track(ui, "chrM", 16_569, &regions, &marks);
    }

}
