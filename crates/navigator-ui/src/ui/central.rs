//! `impl NavigatorApp` methods extracted from `ui.rs` (the `central` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    /// Route files dropped onto the window through the unified importer, attaching them to
    /// the selected subject (auto-detected). No-op when nothing was dropped.
    pub(crate) fn handle_file_drops(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if dropped.is_empty() {
            return;
        }
        let Some(guid) = self.selected_sample else {
            self.status = "Select a subject before dropping data files.".into();
            return;
        };
        // Route every dropped path (files and/or folders) through the batch importer in one go —
        // folders are walked for data files; the result comes back as a single summary modal.
        let paths: Vec<std::path::PathBuf> = dropped.into_iter().filter_map(|f| f.path).collect();
        if !paths.is_empty() {
            self.status = format!("Importing {} dropped item(s)…", paths.len());
            let _ = self.tx.send(Command::AddDataBatch {
                biosample_guid: guid,
                paths,
            });
        }
    }

    /// While files are being dragged over the window, dim the screen and show whether the
    /// drop will land on a subject.
    pub(crate) fn paint_drop_hint(&self, ctx: &egui::Context) {
        if ctx.input(|i| i.raw.hovered_files.is_empty()) {
            return;
        }
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("drop_hint")));
        let rect = ctx.screen_rect();
        painter.rect_filled(rect, 0.0, egui::Color32::from_black_alpha(160));
        let text = if self.selected_sample.is_some() {
            "Drop to add data to this subject"
        } else {
            "Select a subject first, then drop"
        };
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            text,
            egui::FontId::proportional(24.0),
            egui::Color32::WHITE,
        );
    }

    /// The Projects work area: the open project's samples + coverage/haplogroup report.
    pub(crate) fn projects_central(&mut self, ui: &mut egui::Ui) {
        let Some(pid) = self.selected_project else {
            empty_state(ui, self.tr("empty.projects.title"), self.tr("empty.projects.hint"));
            return;
        };
        if let Some(ov) = self.overview.iter().find(|o| o.project.id == pid) {
            let proj = ov.project.clone();
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.vertical(|ui| {
                    ui.heading(&proj.name);
                    ui.label(egui::RichText::new(format!("Administrator: {}", proj.administrator)).weak());
                    if let Some(d) = &proj.description {
                        ui.label(egui::RichText::new(d).weak().small());
                    }
                });
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .add(
                            egui::Button::new(
                                egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE),
                            )
                            .fill(DANGER),
                        )
                        .clicked()
                    {
                        self.confirm_delete_project = Some((proj.id, proj.name.clone()));
                    }
                    if ui.button(self.tr("common.edit")).clicked() {
                        self.edit_project = Some(EditProject {
                            id: proj.id,
                            name: proj.name.clone(),
                            description: proj.description.clone().unwrap_or_default(),
                            administrator: proj.administrator.clone(),
                        });
                    }
                });
            });
            ui.separator();
        }
        // Members vs Report on tabs — both can run to thousands of rows, so they don't stack.
        ui.add_space(4.0);
        self.project_tab = self.sub_bar(ui, self.project_tab, &ProjectTab::ALL);
        match self.project_tab {
            ProjectTab::Members => self.samples_section(ui),
            ProjectTab::Report => self.project_report_section(ui),
            ProjectTab::Ystr => self.project_ystr_section(ui),
        }
    }

    /// The Subjects work area: the selected subject's detail — header + sub-tabs.
    /// A segmented sub-tab bar (one row of selectable labels + a separator). Takes the current
    /// selection by value and returns the new one, so callers avoid a `&mut self.field` borrow clash
    /// with `self.tr` inside the row.
    pub(crate) fn sub_bar<T: Copy + PartialEq>(&self, ui: &mut egui::Ui, current: T, items: &[(T, &'static str)]) -> T {
        let mut sel = current;
        ui.horizontal(|ui| {
            for (variant, key) in items {
                ui.selectable_value(&mut sel, *variant, self.tr(key));
            }
        });
        ui.separator();
        sel
    }

    /// First-run call-to-action for the Simple-mode empty workspace: an "Import DNA" primary action
    /// (pick file(s) → create a subject + import in one step) and a secondary "Add New Subject" that
    /// reveals the inline name form. Without this the empty state told the user to add a subject but
    /// offered no control to do so (Simple mode hides the left panel until a subject exists).
    fn simple_first_run(&mut self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.add_space(56.0);
            ui.heading(self.tr("firstRun.title"));
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("firstRun.hint")).weak());
            ui.add_space(20.0);

            // Import DNA (primary): create a subject named from the file and import in one step.
            if ui
                .add(
                    egui::Button::new(egui::RichText::new(self.tr("firstRun.importDna")).color(egui::Color32::WHITE))
                        .fill(ACCENT)
                        .min_size(egui::vec2(240.0, 34.0)),
                )
                .on_hover_text(self.tr("firstRun.importDnaHint"))
                .clicked()
            {
                if let Some(paths) = rfd::FileDialog::new()
                    .add_filter(
                        "data files",
                        &[
                            "vcf", "gz", "bgz", "csv", "tsv", "txt", "fa", "fasta", "fna", "fas", "bam", "cram",
                        ],
                    )
                    .pick_files()
                {
                    let donor =
                        first_run_subject_name(&paths).unwrap_or_else(|| self.tr("firstRun.defaultName").to_string());
                    self.status = format!("Importing {} file(s) into a new subject…", paths.len());
                    let _ = self.tx.send(Command::CreateSubjectAndImport {
                        donor_identifier: donor,
                        sex: None,
                        paths,
                    });
                }
            }

            ui.add_space(10.0);
            // Add New Subject (secondary): name a subject now, add data later.
            if ui.button(self.tr("subjects.addNew")).clicked() {
                self.forms.show_add_subject = !self.forms.show_add_subject;
            }
            if self.forms.show_add_subject {
                ui.add_space(8.0);
                self.add_subject_form(ui);
            }
        });
    }

    /// When the selected subject has imported data that hasn't been analyzed yet, show a prominent
    /// call-to-action to run the analysis pipeline — the brief stays empty until it runs. Shown only
    /// in the `Pending` state (has alignments, coverage not yet computed); hidden once analysis
    /// completes (→ `Complete`) or when the subject has nothing to analyze (no status row).
    pub(crate) fn simple_analyze_prompt(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        // Gate on the brief itself (rebuilt whenever the subject's data changes — import, clear,
        // delete+re-add) rather than the separate `subject_status` census map, whose async refresh
        // lagged behind those flows and left the prompt missing.
        let needs_analysis = matches!(&self.subject_brief, Some((g, b)) if *g == guid && b.needs_analysis);
        if !needs_analysis {
            return;
        }
        let running = self.analysis.is_some();
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.label(egui::RichText::new(self.tr("brief.analyzeTitle")).strong());
            ui.add_space(2.0);
            ui.label(egui::RichText::new(self.tr("brief.analyzeHint")).weak().small());
            ui.add_space(8.0);
            if ui
                .add_enabled(
                    !running,
                    egui::Button::new(egui::RichText::new(self.tr("brief.analyzeAction")).color(egui::Color32::WHITE))
                        .fill(ACCENT)
                        .min_size(egui::vec2(200.0, 32.0)),
                )
                .clicked()
            {
                self.start_analysis_for_subject(guid);
            }
        });
        ui.add_space(6.0);
    }

    pub(crate) fn subjects_central(&mut self, ui: &mut egui::Ui) {
        let Some(guid) = self.selected_sample else {
            // First launch (Simple mode, empty workspace) has no left panel and so no reachable
            // Add-Subject button — give it a real call-to-action instead of a dead-end hint.
            if self.ui_mode == UiMode::Simple && self.all_biosamples.is_empty() {
                self.simple_first_run(ui);
            } else {
                empty_state(ui, self.tr("empty.subjects.title"), self.tr("empty.subjects.hint"));
            }
            return;
        };
        self.subject_detail_header(ui, guid);
        ui.add_space(6.0);
        // Simple mode hides the per-DNA-type tabs — the subject view *is* the plain-language brief.
        if self.ui_mode == UiMode::Simple {
            ui.separator();
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.add_space(4.0);
                // A background reference-genome download (kicked off after import) only otherwise
                // renders in the Advanced Projects panel — surface its progress bar here too so a
                // Simple-mode user on a slow connection sees why the app is busy. No-ops when idle.
                self.reference_prompt(ui);
                // When this subject has data but hasn't been analyzed yet, offer a prominent Analyze
                // action right at the top — the pipeline (coverage → haplogroups → ancestry) must run
                // before the brief has anything to show.
                self.simple_analyze_prompt(ui, guid);
                // Optional AI-assisted narration sits above the facts (additive, clearly labelled).
                self.simple_ai_narration(ui, guid);
                // The brief's paternal/maternal cards now render the compact descent path in place of
                // the old lineage trail (see brief_descent_trail).
                self.subject_brief_view(ui, guid);
                // Export the brief as a shareable "DNA Story" once it's built.
                if matches!(&self.subject_brief, Some((g, _)) if *g == guid) {
                    ui.add_space(8.0);
                    self.export_row(ui, &[navigator_app::ExportRequest::SubjectBriefHtml(guid)]);
                }
                ui.add_space(10.0);
                // Relatives are live/online, so they render outside the precomputed brief.
                self.simple_relatives_section(ui);
                // Ask-my-results chat (local AI; enabled-gated).
                ui.add_space(10.0);
                self.simple_chat_section(ui, guid);
                // A discreet bridge to the full power-user view for the curious.
                ui.add_space(14.0);
                ui.separator();
                ui.add_space(4.0);
                if ui
                    .link(self.tr("brief.seeData"))
                    .on_hover_text(self.tr("brief.seeDataHint"))
                    .clicked()
                {
                    self.enter_advanced_mode();
                }
            });
            return;
        }
        // Advanced shows the full tab strip.
        ui.horizontal(|ui| {
            for (tab, key) in DetailTab::ALL {
                if ui
                    .selectable_label(self.detail_tab == tab, egui::RichText::new(self.tr(key)).strong())
                    .clicked()
                {
                    self.detail_tab = tab;
                }
            }
        });
        ui.separator();
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.add_space(4.0);
            match self.detail_tab {
                DetailTab::Overview => self.overview_dashboard(ui, guid),
                DetailTab::YDna => {
                    self.y_sub = self.sub_bar(ui, self.y_sub, &YSub::ALL);
                    match self.y_sub {
                        // Compact landing: the subject's Y consensus (source of truth across sources).
                        YSub::Haplogroup => {
                            card(ui, self.tr("card.yHaplogroup"), |ui| {
                                if self.consensus_y.is_some() {
                                    self.consensus_block(ui, "Y-DNA", DnaType::Y);
                                } else {
                                    ui.label(egui::RichText::new(self.tr("hint.noConsensusYet")).weak());
                                }
                            });
                            if self.consensus_y.is_some() {
                                ui.add_space(10.0);
                                card(ui, self.tr("card.descentReport"), |ui| {
                                    self.descent_card(ui, guid, DnaType::Y, false);
                                });
                                ui.add_space(10.0);
                                card(ui, self.tr("card.branchReport"), |ui| {
                                    self.branch_card(ui, guid, DnaType::Y);
                                });
                            }
                        }
                        // The heavy SNP surface: a compact chrY variant track as shared context, then
                        // the heavy tables one at a time (each runs to thousands of rows on a WGS).
                        YSub::Snp => {
                            self.ensure_y_snp_names(guid);
                            card(ui, self.tr("card.variantTrack"), |ui| self.y_variant_track(ui));
                            ui.add_space(10.0);
                            self.y_snp_sub = self.sub_bar(ui, self.y_snp_sub, &YSnpSub::ALL);
                            match self.y_snp_sub {
                                YSnpSub::Profile => {
                                    card(ui, self.tr("card.yVariantProfile"), |ui| self.y_variant_profile_section(ui, guid));
                                }
                                YSnpSub::Private => {
                                    card(ui, self.tr("card.privateYUnion"), |ui| self.donor_private_y_section(ui));
                                }
                                YSnpSub::Imported => {
                                    card(ui, self.tr("card.snpVariants"), |ui| self.variants_section(ui, guid));
                                }
                            }
                        }
                        // STR analysis, separated from SNP.
                        YSub::Str => {
                            card(ui, self.tr("card.ystrConsensus"), |ui| self.ystr_report_section(ui));
                            ui.add_space(10.0);
                            card(ui, self.tr("card.ystrSequence"), |ui| self.ystr_sequence_section(ui, guid));
                            ui.add_space(10.0);
                            card(ui, self.tr("card.yMatches"), |ui| self.ymatch_section(ui, guid));
                        }
                    }
                }
                DetailTab::MtDna => {
                    self.mt_sub = self.sub_bar(ui, self.mt_sub, &MtSub::ALL);
                    match self.mt_sub {
                        MtSub::Summary => {
                            card(ui, self.tr("card.mtHaplogroup"), |ui| {
                                if self.consensus_mt.is_some() {
                                    self.consensus_block(ui, "mtDNA", DnaType::Mt);
                                } else {
                                    ui.label(egui::RichText::new(self.tr("hint.noConsensusYet")).weak());
                                }
                            });
                            if self.consensus_mt.is_some() {
                                ui.add_space(10.0);
                                card(ui, self.tr("card.descentReport"), |ui| {
                                    self.descent_card(ui, guid, DnaType::Mt, false);
                                });
                                ui.add_space(10.0);
                                card(ui, self.tr("card.branchReport"), |ui| {
                                    self.branch_card(ui, guid, DnaType::Mt);
                                });
                            }
                        }
                        MtSub::Variants => {
                            card(ui, self.tr("card.variantTrack"), |ui| self.mt_variant_track(ui));
                            ui.add_space(10.0);
                            card(ui, self.tr("card.mtVariantProfile"), |ui| self.mt_variant_profile_section(ui, guid));
                            ui.add_space(10.0);
                            card(ui, self.tr("card.mtSequences"), |ui| self.mtdna_section(ui, guid));
                        }
                    }
                }
                DetailTab::Autosomal => {
                    self.auto_sub = self.sub_bar(ui, self.auto_sub, &AutoSub::ALL);
                    match self.auto_sub {
                        AutoSub::Summary => {
                            card(ui, self.tr("card.autosomalConsensus"), |ui| {
                                self.autosomal_summary_section(ui, guid);
                                ui.add_space(6.0);
                                ui.label(
                                    egui::RichText::new(self.tr("hint.consensusVcf")).weak().small(),
                                );
                                self.export_row(ui, &[navigator_app::ExportRequest::ConsensusDiploidVcf(guid)]);
                            });
                        }
                        AutoSub::Profile => {
                            card(ui, self.tr("card.autosomalConsensus"), |ui| self.autosomal_profile_table(ui));
                        }
                    }
                }
                DetailTab::Ancestry => {
                    // Consensus is the source of truth: estimate from the subject's pooled autosomal
                    // consensus (no per-alignment BAM walk), decoupled from any selected source.
                    card(ui, self.tr("card.donorAncestry"), |ui| {
                        ui.horizontal(|ui| {
                            let label = if self.donor_ancestry.is_some() { self.tr("common.refresh") } else { self.tr("btn.estimateAncestry") };
                            if ui.add_enabled(!self.estimating_donor_ancestry, egui::Button::new(label)).clicked() {
                                self.estimating_donor_ancestry = true;
                                self.status = "Estimating ancestry from consensus…".into();
                                let _ = self.tx.send(Command::EstimateAncestryFromConsensus { biosample_guid: guid });
                            }
                            if self.estimating_donor_ancestry {
                                ui.spinner();
                            }
                            ui.label(egui::RichText::new(self.tr("hint.ancestryConsensus")).weak().small());
                        });
                        asset_status_line(ui, &self.asset_status);
                        self.donor_ancestry_summary(ui);
                        // Publish the subject's consensus ancestry breakdown (one record per method)
                        // — available once it's been estimated.
                        if self.donor_ancestry.is_some() {
                            self.publish_row(ui, "Publish ancestry to PDS", Command::PublishAncestry { biosample_guid: guid });
                        }
                    });
                    // PC1×PC2 projection: the donor against the reference populations.
                    if self.donor_ancestry.is_some() {
                        ui.add_space(10.0);
                        card(ui, self.tr("card.pcaScatter"), |ui| self.pca_scatter_section(ui));
                    }
                    // Detailed reports from the same consensus estimate (persisted alongside the
                    // super-population ADMIXTURE): modern fine populations + ancient components.
                    if let Some(r) = &self.fine_ancestry {
                        ui.add_space(10.0);
                        card(ui, self.tr("card.ancestryModern"), |ui| draw_population_components(ui, r, "anc_fine", 18));
                    }
                    // Deep (ancient) ancestry via qpAdm — its own on-demand trigger: it genotypes the
                    // best whole-genome CHM13 alignment at the full 1240k (~1-2 min), separate from
                    // the fast consensus estimate above.
                    if navigator_app::ANCIENT_ANCESTRY_ENABLED {
                        ui.add_space(10.0);
                        card(ui, self.tr("card.ancestryAncient"), |ui| {
                            ui.horizontal(|ui| {
                                let label = if self.ancient_ancestry.is_some() {
                                    self.tr("common.refresh")
                                } else {
                                    self.tr("btn.estimateDeepAncestry")
                                };
                                if ui
                                    .add_enabled(!self.estimating_deep_ancestry, egui::Button::new(label))
                                    .clicked()
                                {
                                    self.estimating_deep_ancestry = true;
                                    self.status = "Estimating deep ancestry (qpAdm, ~1–2 min)…".into();
                                    let _ = self.tx.send(Command::EstimateDeepAncestry { biosample_guid: guid });
                                }
                                if self.estimating_deep_ancestry {
                                    ui.spinner();
                                }
                                ui.label(egui::RichText::new(self.tr("hint.deepAncestry")).weak().small());
                            });
                            if let Some(r) = &self.ancient_ancestry {
                                draw_population_components(ui, r, "anc_ancient", 18);
                            }
                        });
                    }
                    // Chromosome painting (diploid local ancestry) — its own section, from the consensus.
                    ui.add_space(10.0);
                    card(ui, self.tr("card.chromosomePainting"), |ui| {
                        ui.horizontal(|ui| {
                            let painted = matches!(&self.painting, Some((id, s)) if *id == navigator_app::CONSENSUS_SOURCE_ID && !s.is_empty());
                            let label = if painted { self.tr("common.refresh") } else { self.tr("ancestry.paint") };
                            if ui.add_enabled(!self.painting_running, egui::Button::new(label)).clicked() {
                                self.painting_running = true;
                                self.status = "Painting local ancestry from consensus…".into();
                                let _ = self.tx.send(Command::PaintAncestryFromConsensus { biosample_guid: guid });
                            }
                            if self.painting_running {
                                ui.spinner();
                            }
                            ui.label(egui::RichText::new(self.tr("hint.chromosomePainting")).weak().small());
                        });
                        if let Some((id, segs)) = &self.painting {
                            if *id == navigator_app::CONSENSUS_SOURCE_ID && !segs.is_empty() {
                                ui.add_space(8.0);
                                draw_chromosome_painting(ui, segs);
                            }
                        }
                    });
                }
                DetailTab::Sources => self.sources_tab(ui, guid),
                DetailTab::IbdMatches => {
                    // Subject-level IBD over the pooled consensus is the primary path.
                    card(ui, self.tr("card.consensusIbd"), |ui| self.consensus_ibd_section(ui, guid));
                    ui.add_space(10.0);
                    card(ui, self.tr("card.networkSuggestions"), |ui| self.network_suggestions_section(ui));
                    ui.add_space(10.0);
                    card(ui, self.tr("card.encryptedExchange"), |ui| self.exchange_section(ui, guid));
                    // Per-source compare + within-subject identity (the QC gate) — advanced.
                    if self.selected_alignment.is_some() {
                        ui.add_space(10.0);
                        let per_source = self.tr("card.panelGenotypingIbd");
                        egui::CollapsingHeader::new(per_source).id_salt("ibd_per_source").show(ui, |ui| self.genotyping_section(ui));
                    }
                }
            }
        });
    }

    /// The subject-detail header: big name, ID + sex, and Add Data / Edit / Delete actions.
    fn subject_detail_header(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let bio = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|b| b.guid == guid)
            .cloned();
        let Some(bio) = bio else { return };
        ui.add_space(6.0);
        // When the subject was opened from a project's report, offer a way back to that project.
        if let Some(pid) = self.return_to_project {
            if ui.button(self.tr("detail.backToProject")).clicked() {
                self.nav = Nav::Projects;
                self.selected_project = Some(pid);
                self.project_tab = ProjectTab::Report;
                self.return_to_project = None;
            }
            ui.add_space(2.0);
        }
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.heading(&bio.donor_identifier);
                ui.label(egui::RichText::new(bio.sex.as_deref().unwrap_or("Unknown")).weak());
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE))
                            .fill(DANGER),
                    )
                    .clicked()
                {
                    self.confirm_delete = Some(guid);
                }
                ui.menu_button(self.tr("common.clearData"), |ui| {
                    if ui
                        .button(self.tr("resetHaplo.action"))
                        .on_hover_text(self.tr("resetHaplo.hint"))
                        .clicked()
                    {
                        self.confirm_reset_haplo = Some(guid);
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui
                        .button(self.tr("clear.allData"))
                        .on_hover_text(self.tr("clear.hint"))
                        .clicked()
                    {
                        self.confirm_clear = Some(guid);
                        ui.close_menu();
                    }
                });
                if ui.button(self.tr("common.edit")).clicked() {
                    self.edit_subject = Some(EditSubject {
                        guid,
                        donor_identifier: bio.donor_identifier.clone(),
                        sample_accession: bio.sample_accession.clone().unwrap_or_default(),
                        description: bio.description.clone().unwrap_or_default(),
                        center_name: bio.center_name.clone().unwrap_or_default(),
                        sex: bio.sex.clone().unwrap_or_default(),
                    });
                }
                ui.menu_button(self.tr("detail.addData"), |ui| {
                    if ui.button(self.tr("detail.addFiles")).clicked() {
                        ui.close_menu();
                        if let Some(paths) = rfd::FileDialog::new()
                            .add_filter(
                                "data files",
                                &[
                                    "vcf", "gz", "csv", "tsv", "txt", "fa", "fasta", "fna", "fas", "bam", "cram",
                                ],
                            )
                            .pick_files()
                        {
                            self.status = format!("Importing {} file(s)…", paths.len());
                            let _ = self.tx.send(Command::AddDataBatch {
                                biosample_guid: guid,
                                paths,
                            });
                        }
                    }
                    if ui.button(self.tr("detail.addFolder")).clicked() {
                        ui.close_menu();
                        if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                            self.status = format!("Scanning {}…", dir.display());
                            let _ = self.tx.send(Command::AddDataBatch {
                                biosample_guid: guid,
                                paths: vec![dir],
                            });
                        }
                    }
                });
            });
        });
    }

    /// A simple at-a-glance dashboard: counts + account state.
    pub(crate) fn dashboard_central(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);
        ui.heading(self.tr("dash.title"));
        ui.add_space(12.0);
        let (projects, subjects, alignments) = (
            self.overview.len(),
            self.all_biosamples.len(),
            self.all_alignments.len(),
        );
        let (lp, ls, la) = (
            self.tr("dash.projects"),
            self.tr("dash.subjects"),
            self.tr("dash.alignments"),
        );
        ui.horizontal_wrapped(|ui| {
            stat_card(ui, lp, projects);
            stat_card(ui, ls, subjects);
            stat_card(ui, la, alignments);
        });
        ui.add_space(16.0);
        match &self.account {
            Some(did) => {
                ui.label(format!("Signed in as {did}"));
                ui.label(if self.online { "● online" } else { "○ offline" });
            }
            None => {
                ui.label(egui::RichText::new("Not signed in — connect a PDS from the top bar to publish.").weak());
            }
        }
    }

    /// The bottom action bar for the Subjects view: selection count + batch actions.
    pub(crate) fn action_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("actions").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                let selected = self.selected_sample.is_some();
                ui.label(format!(
                    "{} {}",
                    if selected { 1 } else { 0 },
                    self.tr("action.selected")
                ));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .add_enabled(selected, egui::Button::new(self.tr("action.addToProject")))
                        .clicked()
                    {
                        if let Some(guid) = self.selected_sample {
                            let current = self
                                .all_biosamples
                                .iter()
                                .chain(self.samples.iter())
                                .find(|b| b.guid == guid)
                                .and_then(|b| b.project_id);
                            self.assign_project = Some((guid, current));
                        }
                    }
                    if ui
                        .add_enabled(selected, egui::Button::new(self.tr("action.batchAnalyze")))
                        .clicked()
                    {
                        if let Some(id) = self.selected_alignment {
                            self.start_full_analysis(id);
                        } else {
                            self.status = "Select an alignment (Data Sources) to run analysis.".into();
                        }
                    }
                    // Compare needs a second subject (multi-select) — disabled for now.
                    let _ = ui.add_enabled(false, egui::Button::new(self.tr("action.compare")));
                });
            });
            ui.add_space(2.0);
        });
    }
}

/// A friendly default subject name derived from the first picked file's stem, stripping the known
/// data-file extensions (including double extensions like `.vcf.gz`). `None` if nothing usable is
/// left, so the caller can fall back to a localized default. The user can rename later.
fn first_run_subject_name(paths: &[std::path::PathBuf]) -> Option<String> {
    const EXTS: [&str; 18] = [
        ".g.vcf.gz", ".vcf.gz", ".vcf.bgz", ".fasta.gz", ".fa.gz", ".fna.gz", ".bam", ".cram", ".vcf", ".fasta",
        ".fa", ".fna", ".fas", ".csv", ".tsv", ".txt", ".gz", ".bgz",
    ];
    let name = paths.first()?.file_name()?.to_str()?;
    let lower = name.to_ascii_lowercase();
    let stem = EXTS
        .iter()
        .find_map(|ext| lower.strip_suffix(ext).map(|s| &name[..s.len()]))
        .unwrap_or(name);
    let stem = stem.trim();
    (!stem.is_empty()).then(|| stem.to_string())
}

#[cfg(test)]
mod tests {
    use super::first_run_subject_name;
    use std::path::PathBuf;

    fn name(p: &str) -> Option<String> {
        first_run_subject_name(&[PathBuf::from(p)])
    }

    #[test]
    fn derives_subject_name_from_file_stem() {
        assert_eq!(name("/data/HG002.bam").as_deref(), Some("HG002"));
        assert_eq!(name("HG00096.chm13.chrY.g.vcf.gz").as_deref(), Some("HG00096.chm13.chrY"));
        assert_eq!(name("MyKit.vcf.gz").as_deref(), Some("MyKit"));
        assert_eq!(name("genome_Full.CRAM").as_deref(), Some("genome_Full")); // extension match is case-insensitive
        assert_eq!(name("relative.fasta").as_deref(), Some("relative"));
        // No recognized extension → keep the whole name; empty/none → fall back handled by caller.
        assert_eq!(name("plainname").as_deref(), Some("plainname"));
        assert_eq!(first_run_subject_name(&[]), None);
    }
}
