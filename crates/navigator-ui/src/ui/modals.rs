//! `impl NavigatorApp` methods extracted from `ui.rs` (the `modals` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    /// The "Full Analysis" progress modal: a dimmed backdrop + centered card with the current
    /// step, a progress bar + percent, and a Cancel button. Shown while `self.analysis` is set.
    pub(crate) fn analysis_modal(&mut self, ctx: &egui::Context) {
        let Some(p) = self.analysis.clone() else { return };
        // Dim everything behind the dialog.

        modal_frame(ctx, "analysis_modal", 460.0, |ui| {
            ui.label(egui::RichText::new("Full Analysis").strong().size(16.0));
            ui.separator();
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(egui::RichText::new("Analysis in progress…").weak());
                let elapsed = (self.frame_time - p.started).max(0.0) as u64;
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(egui::RichText::new(format!("{:02}:{:02}", elapsed / 60, elapsed % 60)).weak());
                });
            });
            ui.add_space(10.0);
            ui.label(format!("Step {}/{}: {} — {}", p.step, p.total, p.label, p.detail));
            ui.add_space(10.0);
            ui.horizontal(|ui| {
                // `animate` shimmers the bar so a long step reads as working, not stalled.
                ui.add(
                    egui::ProgressBar::new(p.fraction)
                        .desired_width(360.0)
                        .rounding(4.0)
                        .animate(true),
                );
                ui.label(format!("{}%", (p.fraction * 100.0).round() as i32));
            });
            ui.add_space(6.0);
            ui.label(
                egui::RichText::new("Whole-genome steps (coverage) can take several minutes on a WGS BAM.")
                    .weak()
                    .small(),
            );
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button(self.tr("common.cancel")).clicked() {
                    let _ = self.tx.send(Command::CancelAnalysis);
                    self.status = "Cancelling…".into();
                }
            });
        });
    }

    /// The Edit-subject modal: editable fields over a dimmed backdrop. Save sends an
    /// `UpdateBiosample` command; the resulting `BiosamplesChanged` event refreshes the lists.
    pub(crate) fn edit_subject_modal(&mut self, ctx: &egui::Context) {
        let Some(mut edit) = self.edit_subject.clone() else {
            return;
        };

        let mut close = false;
        modal_frame(ctx, "edit_subject_modal", 420.0, |ui| {
            ui.label(egui::RichText::new(self.tr("edit.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(6.0);
            let field = |ui: &mut egui::Ui, label: &str, value: &mut String, hint: &str| {
                ui.label(label);
                ui.add(
                    egui::TextEdit::singleline(value)
                        .hint_text(hint)
                        .desired_width(f32::INFINITY),
                );
                ui.add_space(4.0);
            };
            field(
                ui,
                self.tr("edit.identifier"),
                &mut edit.donor_identifier,
                "donor identifier",
            );
            field(
                ui,
                self.tr("edit.accession"),
                &mut edit.sample_accession,
                "accession (optional)",
            );
            field(
                ui,
                self.tr("edit.description"),
                &mut edit.description,
                "description (optional)",
            );
            field(ui, self.tr("edit.center"), &mut edit.center_name, "center (optional)");
            field(ui, self.tr("edit.sex"), &mut edit.sex, "sex (optional)");
            ui.add_space(10.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add_enabled(
                        !edit.donor_identifier.trim().is_empty(),
                        egui::Button::new(self.tr("common.save")).fill(ACCENT),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::UpdateBiosample {
                        guid: edit.guid,
                        donor_identifier: edit.donor_identifier.trim().to_string(),
                        sample_accession: opt(&edit.sample_accession),
                        description: opt(&edit.description),
                        center_name: opt(&edit.center_name),
                        sex: opt(&edit.sex),
                    });
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.edit_subject = None;
        } else {
            self.edit_subject = Some(edit);
        }
    }

    /// The Settings / Preferences modal: connection (AppView URL, Y-tree provider), appearance
    /// (theme, language, tree-cache TTL), reference genomes (local FASTA + auto-download per build),
    /// and a read-only advanced section. Self-mutation/dispatch is deferred until after the closure
    /// so only `self.tr` (immutable) is used inside it.
    pub(crate) fn settings_modal(&mut self, ctx: &egui::Context) {
        if !self.show_settings {
            return;
        }

        let mut form = self.settings_form.clone();
        let mut theme_dark = self.dark_mode;
        let mut lang = self.lang;
        let prev_lang = self.lang;
        let (mut close, mut save) = (false, false);
        // Deferred actions (dispatched after the closure, since only `self.tr` is used inside it).
        let mut verify_build: Option<String> = None;
        let mut lift_request = false;
        let mut test_llm: Option<String> = None;

        modal_frame(ctx, "settings_modal", 580.0, |ui| {
            ui.label(egui::RichText::new(self.tr("settings.title")).strong().size(16.0));
            ui.separator();
            egui::ScrollArea::vertical().max_height(460.0).show(ui, |ui| {
                // --- Connection ---
                ui.label(egui::RichText::new(self.tr("settings.connection")).strong());
                ui.horizontal(|ui| {
                    ui.label(self.tr("settings.appviewUrl"));
                    ui.add(
                        egui::TextEdit::singleline(&mut form.appview_url)
                            .hint_text("http://localhost:9000")
                            .desired_width(320.0),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label(self.tr("settings.yTreeProvider"));
                    let cur = if form.y_tree_provider.eq_ignore_ascii_case("ftdna") {
                        "FTDNA"
                    } else {
                        "Decoding-Us"
                    };
                    egui::ComboBox::from_id_salt("settings_y_provider")
                        .selected_text(cur)
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut form.y_tree_provider, "decodingus".to_string(), "Decoding-Us");
                            ui.selectable_value(&mut form.y_tree_provider, "ftdna".to_string(), "FTDNA");
                        });
                });
                ui.add_space(8.0);

                // --- Appearance ---
                ui.label(egui::RichText::new(self.tr("settings.appearance")).strong());
                ui.horizontal(|ui| {
                    ui.label(self.tr("settings.theme"));
                    ui.selectable_value(&mut theme_dark, true, self.tr("settings.dark"));
                    ui.selectable_value(&mut theme_dark, false, self.tr("settings.light"));
                });
                ui.horizontal(|ui| {
                    ui.label(self.tr("settings.uiScale"));
                    ui.add(
                        egui::Slider::new(&mut form.ui_scale, 0.8..=2.5)
                            .step_by(0.05)
                            .fixed_decimals(2),
                    );
                    if ui.small_button("100%").clicked() {
                        form.ui_scale = 1.0;
                    }
                });
                ui.horizontal(|ui| {
                    ui.label(self.tr("settings.language"));
                    egui::ComboBox::from_id_salt("settings_lang")
                        .selected_text(lang.label())
                        .show_ui(ui, |ui| {
                            for &l in crate::i18n::Lang::all() {
                                ui.selectable_value(&mut lang, l, l.label());
                            }
                        });
                });
                ui.horizontal(|ui| {
                    ui.label(self.tr("settings.treeTtl"));
                    ui.add(egui::TextEdit::singleline(&mut form.tree_ttl_days).desired_width(60.0));
                });
                ui.add_space(8.0);

                // --- AI assistant (local LLM) ---
                ui.label(egui::RichText::new(self.tr("settings.ai")).strong());
                ui.checkbox(&mut form.llm_enabled, self.tr("settings.ai.enable"));
                ui.add_enabled_ui(form.llm_enabled, |ui| {
                    ui.horizontal(|ui| {
                        ui.label(self.tr("settings.ai.baseUrl"));
                        ui.add(
                            egui::TextEdit::singleline(&mut form.llm_base_url)
                                .hint_text(navigator_app::llm::DEFAULT_LLM_BASE_URL)
                                .desired_width(300.0),
                        );
                    });
                    // Quick-pick host ports.
                    ui.horizontal(|ui| {
                        ui.label(self.tr("settings.ai.presets"));
                        if ui.small_button("LM Studio").clicked() {
                            form.llm_base_url = "http://localhost:1234/v1".into();
                        }
                        if ui.small_button("Ollama").clicked() {
                            form.llm_base_url = "http://localhost:11434/v1".into();
                        }
                        if ui.small_button("llama.cpp").clicked() {
                            form.llm_base_url = "http://localhost:8080/v1".into();
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui
                            .add_enabled(!self.llm_testing, egui::Button::new(self.tr("settings.ai.test")))
                            .clicked()
                        {
                            test_llm = Some(form.llm_base_url.trim().to_string());
                        }
                        if self.llm_testing {
                            ui.spinner();
                        }
                        if let Some(msg) = &self.llm_test_msg {
                            ui.label(egui::RichText::new(msg).weak().small());
                        }
                    });
                    // Model picker — populated by a successful Test connection.
                    ui.horizontal(|ui| {
                        ui.label(self.tr("settings.ai.model"));
                        let current = if form.llm_model.is_empty() {
                            self.tr("settings.ai.modelAuto").to_string()
                        } else {
                            form.llm_model.clone()
                        };
                        egui::ComboBox::from_id_salt("settings_llm_model")
                            .selected_text(current)
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut form.llm_model,
                                    String::new(),
                                    self.tr("settings.ai.modelAuto"),
                                );
                                for m in &self.llm_models {
                                    ui.selectable_value(&mut form.llm_model, m.clone(), m);
                                }
                            });
                    });
                    // Privacy line — turns to a warning for a non-loopback URL.
                    if navigator_app::llm::is_loopback_url(&form.llm_base_url) {
                        ui.label(egui::RichText::new(self.tr("settings.ai.local")).weak().small());
                    } else {
                        ui.label(
                            egui::RichText::new(self.tr("settings.ai.remoteWarn"))
                                .small()
                                .color(egui::Color32::from_rgb(230, 170, 80)),
                        );
                    }
                });
                ui.add_space(8.0);

                // --- Reference genomes ---
                ui.label(egui::RichText::new(self.tr("settings.references")).strong());
                ui.checkbox(&mut form.prompt_before_download, self.tr("settings.promptDownload"));
                egui::Grid::new("settings_refs")
                    .striped(true)
                    .num_columns(5)
                    .show(ui, |ui| {
                        for h in [
                            "settings.build",
                            "settings.status",
                            "settings.localFasta",
                            "settings.autoDownload",
                            "settings.integrity",
                        ] {
                            ui.strong(self.tr(h));
                        }
                        ui.end_row();
                        for row in &mut form.references {
                            ui.label(&row.build);
                            ui.label(egui::RichText::new(&row.status).weak());
                            ui.horizontal(|ui| {
                                ui.add(
                                    egui::TextEdit::singleline(&mut row.local_path)
                                        .hint_text("(none)")
                                        .desired_width(180.0),
                                );
                                if ui.button("📂").on_hover_text(self.tr("settings.browse")).clicked() {
                                    if let Some(p) = rfd::FileDialog::new()
                                        .add_filter("FASTA", &["fa", "fasta", "fna", "gz"])
                                        .pick_file()
                                    {
                                        row.local_path = p.display().to_string();
                                    }
                                }
                            });
                            ui.checkbox(&mut row.auto_download, "");
                            ui.horizontal(|ui| {
                                if ui.small_button(self.tr("settings.verify")).clicked() {
                                    verify_build = Some(row.build.clone());
                                }
                                if !row.verify.is_empty() {
                                    ui.label(egui::RichText::new(&row.verify).small().weak());
                                }
                            });
                            ui.end_row();
                        }
                    });
                if form.references.is_empty() {
                    ui.label(egui::RichText::new(self.tr("settings.loadingRefs")).weak());
                }
                ui.add_space(8.0);

                // --- Tools: VCF liftover ---
                ui.label(egui::RichText::new(self.tr("liftvcf.title")).strong());
                ui.label(egui::RichText::new(self.tr("liftvcf.hint")).weak().small());
                ui.horizontal(|ui| {
                    ui.label(self.tr("liftvcf.input"));
                    ui.add(
                        egui::TextEdit::singleline(&mut form.lift_in)
                            .hint_text("input.vcf[.gz]")
                            .desired_width(260.0),
                    );
                    if ui.button("📂").clicked() {
                        if let Some(p) = rfd::FileDialog::new().add_filter("VCF", &["vcf", "gz"]).pick_file() {
                            form.lift_in = p.display().to_string();
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label(self.tr("liftvcf.target"));
                    egui::ComboBox::from_id_salt("liftvcf_target")
                        .selected_text(&form.lift_target)
                        .show_ui(ui, |ui| {
                            for b in ["chm13v2.0", "GRCh38", "GRCh37"] {
                                ui.selectable_value(&mut form.lift_target, b.to_string(), b);
                            }
                        });
                    ui.checkbox(&mut form.lift_filter_par, self.tr("liftvcf.filterPar"));
                });
                ui.horizontal(|ui| {
                    ui.label(self.tr("liftvcf.output"));
                    ui.add(
                        egui::TextEdit::singleline(&mut form.lift_out)
                            .hint_text("lifted.vcf[.gz]")
                            .desired_width(260.0),
                    );
                    if ui.button("📂").clicked() {
                        if let Some(p) = rfd::FileDialog::new()
                            .add_filter("VCF", &["vcf", "gz"])
                            .set_file_name("lifted.vcf")
                            .save_file()
                        {
                            form.lift_out = p.display().to_string();
                        }
                    }
                });
                let lift_ready = !form.lift_in.trim().is_empty() && !form.lift_out.trim().is_empty();
                if ui
                    .add_enabled(lift_ready, egui::Button::new(self.tr("liftvcf.run")))
                    .clicked()
                {
                    lift_request = true;
                }
                ui.add_space(8.0);

                // --- Advanced (read-only) ---
                ui.label(egui::RichText::new(self.tr("settings.advanced")).strong());
                ui.label(
                    egui::RichText::new(format!(
                        "{}: {}",
                        self.tr("settings.cacheDir"),
                        AppSettings::cache_base_dir().display()
                    ))
                    .weak(),
                );
                ui.label(egui::RichText::new(self.tr("settings.advancedEnv")).weak());
            });
            ui.separator();
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .button(egui::RichText::new(self.tr("common.save")).strong())
                    .clicked()
                {
                    save = true;
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });

        // Live-apply theme + UI scale + language (immediate feedback).
        if theme_dark != self.dark_mode {
            self.dark_mode = theme_dark;
            apply_theme(ctx, self.dark_mode);
        }
        if (ctx.zoom_factor() - form.ui_scale).abs() > f32::EPSILON {
            ctx.set_zoom_factor(form.ui_scale.clamp(0.5, 3.0));
        }
        if lang != prev_lang {
            self.lang = lang;
            crate::i18n::save_lang(lang);
        }

        if save {
            let appview = form.appview_url.trim().to_string();
            let settings = AppSettings {
                y_tree_provider: Some(form.y_tree_provider.clone()),
                appview_url: (!appview.is_empty()).then_some(appview),
                tree_ttl_days: form.tree_ttl_days.trim().parse::<u64>().ok(),
                theme: Some(if self.dark_mode {
                    "dark".to_string()
                } else {
                    "light".to_string()
                }),
                prompt_before_download: Some(form.prompt_before_download),
                ui_scale: Some(form.ui_scale),
                // Interface mode is toggled from the app bar, not this dialog — preserve it.
                ui_mode: AppSettings::load().ui_mode,
                llm_enabled: Some(form.llm_enabled),
                llm_base_url: {
                    let u = form.llm_base_url.trim().to_string();
                    (!u.is_empty()).then_some(u)
                },
                llm_model: {
                    let m = form.llm_model.trim().to_string();
                    (!m.is_empty()).then_some(m)
                },
            };
            match settings.save() {
                Ok(()) => self.status = self.tr("settings.saved").to_string(),
                Err(e) => self.status = format!("Could not save settings: {e}"),
            }
            // Reflect the AI toggle immediately (gates the "Polish with AI" affordance).
            self.ai_enabled = form.llm_enabled;
            for row in &form.references {
                let local = row.local_path.trim().to_string();
                let _ = self.tx.send(Command::SetReferenceOverride {
                    build: row.build.clone(),
                    local_path: (!local.is_empty()).then_some(local),
                    auto_download: row.auto_download,
                });
            }
        }

        // Deferred dispatch (only `self.tr` was used inside the closure).
        if let Some(build) = verify_build {
            self.status = format!("Verifying {build}…");
            let _ = self.tx.send(Command::VerifyReference { build });
        }
        if let Some(base_url) = test_llm {
            self.llm_testing = true;
            self.llm_test_msg = Some(self.tr("settings.ai.testing").to_string());
            let _ = self.tx.send(Command::TestLlmConnection { base_url });
        }
        if lift_request {
            self.status = "Lifting VCF…".into();
            let _ = self.tx.send(Command::LiftVcf {
                source: None, // inferred from the VCF header
                target: form.lift_target.clone(),
                in_vcf: std::path::PathBuf::from(form.lift_in.trim()),
                out_vcf: std::path::PathBuf::from(form.lift_out.trim()),
                filter_par: form.lift_filter_par,
            });
        }

        if close {
            self.show_settings = false;
        } else {
            self.settings_form = form;
        }
    }

    /// The Delete-subject confirmation modal. Confirm sends a `DeleteBiosample` command; the app
    /// layer refuses (surfaced via the status bar) when the subject still has dependent data.
    pub(crate) fn delete_subject_modal(&mut self, ctx: &egui::Context) {
        let Some(guid) = self.confirm_delete else { return };
        let name = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|b| b.guid == guid)
            .map(|b| b.donor_identifier.clone())
            .unwrap_or_else(|| guid.0.to_string());

        let mut close = false;
        modal_frame(ctx, "delete_subject_modal", 400.0, |ui| {
            ui.label(egui::RichText::new(self.tr("delete.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(8.0);
            ui.label(format!("{} “{}”?", self.tr("delete.confirm"), name));
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("delete.note")).weak().small());
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE))
                            .fill(DANGER),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::DeleteBiosample(guid));
                    if self.selected_sample == Some(guid) {
                        self.selected_sample = None;
                    }
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.confirm_delete = None;
        }
    }

    /// The Clear-data confirmation modal. Confirm sends a `ClearBiosampleData` command, which resets
    /// the subject's analysis (runs, alignments, haplogroups, ancestry, profiles…) while keeping the
    /// subject itself — the recovery tool for a botched import.
    pub(crate) fn clear_subject_modal(&mut self, ctx: &egui::Context) {
        let Some(guid) = self.confirm_clear else { return };
        let name = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|b| b.guid == guid)
            .map(|b| b.donor_identifier.clone())
            .unwrap_or_else(|| guid.0.to_string());

        let mut close = false;
        modal_frame(ctx, "clear_subject_modal", 420.0, |ui| {
            ui.label(egui::RichText::new(self.tr("clear.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(8.0);
            ui.label(format!("{} “{}”?", self.tr("clear.confirm"), name));
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("clear.note")).weak().small());
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(self.tr("common.clearData")).color(egui::Color32::WHITE))
                            .fill(DANGER),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::ClearBiosampleData(guid));
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.confirm_clear = None;
        }
    }

    /// Confirm resetting only the subject's haplogroup placement (stale-lineage cleanup) — keeps
    /// coverage/ancestry/imported data; the placement re-derives on the next full analysis / re-import.
    pub(crate) fn reset_haplo_modal(&mut self, ctx: &egui::Context) {
        let Some(guid) = self.confirm_reset_haplo else { return };
        let name = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|b| b.guid == guid)
            .map(|b| b.donor_identifier.clone())
            .unwrap_or_else(|| guid.0.to_string());

        let mut close = false;
        modal_frame(ctx, "reset_haplo_modal", 440.0, |ui| {
            ui.label(egui::RichText::new(self.tr("resetHaplo.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(8.0);
            ui.label(format!("{} “{}”?", self.tr("resetHaplo.confirm"), name));
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("resetHaplo.note")).weak().small());
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(
                            egui::RichText::new(self.tr("resetHaplo.action")).color(egui::Color32::WHITE),
                        )
                        .fill(DANGER),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::ClearHaplogroupData(guid));
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.confirm_reset_haplo = None;
        }
    }

    /// Summary modal after a batch Add Data / drag-and-drop: per-file detected type + any skipped
    /// files with the reason. Dismissed with Close.
    pub(crate) fn batch_import_modal(&mut self, ctx: &egui::Context) {
        let Some(summary) = self.batch_import.clone() else {
            return;
        };

        let mut close = false;
        modal_frame(ctx, "batch_import_modal", 460.0, |ui| {
            ui.label(egui::RichText::new(self.tr("import.title")).strong().size(16.0));
            ui.label(
                egui::RichText::new(format!(
                    "{} imported · {} skipped",
                    summary.imported.len(),
                    summary.skipped.len()
                ))
                .weak(),
            );
            ui.separator();
            if summary.imported.is_empty() && summary.skipped.is_empty() {
                ui.label(self.tr("import.none"));
            }
            egui::ScrollArea::vertical().max_height(320.0).show(ui, |ui| {
                egui::Grid::new("batch_import_grid")
                    .num_columns(2)
                    .striped(true)
                    .show(ui, |ui| {
                        for (name, kind) in &summary.imported {
                            ui.colored_label(egui::Color32::from_rgb(60, 160, 60), "✓");
                            ui.label(format!("{name} — {kind}"));
                            ui.end_row();
                        }
                        for (name, reason) in &summary.skipped {
                            ui.colored_label(egui::Color32::from_rgb(190, 140, 40), "•");
                            ui.label(egui::RichText::new(format!("{name} — {reason}")).weak());
                            ui.end_row();
                        }
                    });
            });
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button(self.tr("common.close")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.batch_import = None;
        }
    }

    /// Confirmation modal for deleting a data-source row (run/alignment/profile). Confirm sends
    /// the variant's worker command; the resulting change event refreshes the affected list.
    pub(crate) fn data_delete_modal(&mut self, ctx: &egui::Context) {
        let Some(target) = self.confirm_data_delete.clone() else {
            return;
        };

        let mut close = false;
        modal_frame(ctx, "data_delete_modal", 400.0, |ui| {
            ui.label(egui::RichText::new(self.tr("delete.dataTitle")).strong().size(16.0));
            ui.separator();
            ui.add_space(8.0);
            ui.label(format!("{} {}?", self.tr("delete.confirm"), target.label()));
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("delete.dataNote")).weak().small());
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE))
                            .fill(DANGER),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(target.command());
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.confirm_data_delete = None;
        }
    }

    /// Destructive merge-sequence-runs modal: the `secondary` run's alignments are reparented onto a
    /// chosen `primary` run and the now-empty secondary is deleted. Mirrors the data-delete confirm.
    pub(crate) fn merge_runs_modal(&mut self, ctx: &egui::Context) {
        let Some(mut m) = self.merge_runs.clone() else { return };

        // Run label: "WGS · NovaSeq (#id)".
        let label = |id: i64| -> String {
            self.runs
                .iter()
                .find(|r| r.id == id)
                .map(|r| {
                    format!(
                        "{} · {} (#{})",
                        testtype::display_name(&r.test_type),
                        if r.platform_name.is_empty() {
                            "—"
                        } else {
                            &r.platform_name
                        },
                        r.id
                    )
                })
                .unwrap_or_else(|| format!("run #{id}"))
        };
        let others: Vec<i64> = self.runs.iter().map(|r| r.id).filter(|&id| id != m.secondary).collect();

        let (mut close, mut confirm) = (false, false);
        modal_frame(ctx, "merge_runs_modal", 440.0, |ui| {
            ui.label(egui::RichText::new(self.tr("merge.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(8.0);
            ui.label(format!("{}: {}", self.tr("merge.moveFrom"), label(m.secondary)));
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.label(self.tr("merge.into"));
                let sel = m
                    .primary
                    .map(label)
                    .unwrap_or_else(|| self.tr("merge.pick").to_string());
                egui::ComboBox::from_id_salt("merge_primary")
                    .selected_text(sel)
                    .show_ui(ui, |ui| {
                        for id in &others {
                            ui.selectable_value(&mut m.primary, Some(*id), label(*id));
                        }
                    });
            });
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("merge.note")).weak().small());
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let ready = m.primary.is_some();
                if ui
                    .add_enabled(
                        ready,
                        egui::Button::new(egui::RichText::new(self.tr("merge.run")).color(egui::Color32::WHITE))
                            .fill(DANGER),
                    )
                    .clicked()
                {
                    confirm = true;
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });

        if confirm {
            if let Some(primary) = m.primary {
                self.status = self.tr("merge.running").to_string();
                let _ = self.tx.send(Command::MergeSequenceRuns {
                    biosample_guid: m.guid,
                    primary,
                    secondary: m.secondary,
                });
            }
        }
        if close {
            self.merge_runs = None;
        } else {
            self.merge_runs = Some(m); // keep the picker selection across frames
        }
    }

    /// Read-only Y-profile **source audit**: a per-source provenance table (label · type · method
    /// tier weight · variants contributed) and a per-conflict evidence list (each conflicting variant
    /// with every source's call), so the user can see what drove — or disagreed with — each consensus
    /// call. Pure over the cached `y_profile`; no schema change, no re-genotyping.
    pub(crate) fn y_profile_audit_modal(&mut self, ctx: &egui::Context) {
        if !self.audit_y_profile {
            return;
        }
        let Some(profile) = self.y_profile.clone() else {
            self.audit_y_profile = false;
            return;
        };

        let state_glyph = |s: navigator_domain::consensus::ConsensusState| match s {
            navigator_domain::consensus::ConsensusState::Derived => "derived",
            navigator_domain::consensus::ConsensusState::Ancestral => "ancestral",
            navigator_domain::consensus::ConsensusState::NoCall => "no-call",
        };
        let mut close = false;
        modal_frame(ctx, "y_profile_audit_modal", 560.0, |ui| {
            ui.label(egui::RichText::new(self.tr("audit.title")).strong().size(16.0));
            if let Some(t) = &profile.terminal {
                ui.label(egui::RichText::new(format!("terminal {t}")).weak());
            }
            ui.separator();
            egui::ScrollArea::vertical().max_height(440.0).show(ui, |ui| {
                // --- Per-source provenance ---
                ui.label(egui::RichText::new(self.tr("audit.sources")).strong());
                egui::Grid::new("yaudit_sources")
                    .striped(true)
                    .num_columns(4)
                    .show(ui, |ui| {
                        for h in ["audit.source", "audit.type", "audit.tier", "audit.variants"] {
                            ui.strong(self.tr(h));
                        }
                        ui.end_row();
                        for s in &profile.sources {
                            ui.label(&s.label);
                            ui.label(egui::RichText::new(s.source_type.as_str()).small());
                            ui.label(format!("{:.2}", s.source_type.snp_weight()));
                            ui.label(s.variant_count.to_string());
                            ui.end_row();
                        }
                    });
                ui.label(egui::RichText::new(self.tr("audit.tierNote")).weak().small());
                ui.add_space(10.0);

                // --- Conflicts: who disagrees, at which variant ---
                let conflicts: Vec<_> = profile
                    .variants
                    .iter()
                    .filter(|v| v.status == YVariantStatus::Conflict)
                    .collect();
                ui.label(egui::RichText::new(format!("{} ({})", self.tr("audit.conflicts"), conflicts.len())).strong());
                if conflicts.is_empty() {
                    ui.label(egui::RichText::new(self.tr("audit.noConflicts")).weak());
                } else {
                    let amber = egui::Color32::from_rgb(220, 150, 60);
                    for v in conflicts.iter().take(200) {
                        let name = if v.name.is_empty() {
                            format!("@{}", v.position)
                        } else {
                            v.name.clone()
                        };
                        ui.label(
                            egui::RichText::new(format!("{name}  ·  {}", v.position))
                                .color(amber)
                                .strong(),
                        );
                        ui.indent(("yaudit", v.position), |ui| {
                            for src in &v.sources {
                                ui.label(
                                    egui::RichText::new(format!(
                                        "{} ({}, tier {:.2}) → {}",
                                        src.label,
                                        src.source_type.as_str(),
                                        src.source_type.snp_weight(),
                                        state_glyph(src.state)
                                    ))
                                    .small(),
                                );
                            }
                        });
                    }
                }
            });
            ui.separator();
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button(self.tr("common.close")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.audit_y_profile = false;
        }
    }

    /// The Add-to-Project picker: a dropdown of projects (plus "no project"). Save sends
    /// `AssignBiosampleProject`; the resulting `BiosamplesChanged` event refreshes the lists.
    pub(crate) fn assign_project_modal(&mut self, ctx: &egui::Context) {
        let Some((guid, mut chosen)) = self.assign_project else {
            return;
        };

        let selected_text = match chosen {
            Some(pid) => self
                .overview
                .iter()
                .find(|o| o.project.id == pid)
                .map(|o| o.project.name.clone())
                .unwrap_or_else(|| format!("project {pid}")),
            None => self.tr("projects.noProject").to_string(),
        };
        let mut close = false;
        let mut commit = false;
        modal_frame(ctx, "assign_project_modal", 360.0, |ui| {
            ui.label(egui::RichText::new(self.tr("action.addToProject")).strong().size(16.0));
            ui.separator();
            ui.add_space(8.0);
            egui::ComboBox::from_id_salt("assign_project_combo")
                .selected_text(selected_text)
                .width(300.0)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut chosen, None, self.tr("projects.noProject"));
                    for o in &self.overview {
                        ui.selectable_value(&mut chosen, Some(o.project.id), &o.project.name);
                    }
                });
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.add(egui::Button::new(self.tr("common.save")).fill(ACCENT)).clicked() {
                    commit = true;
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if commit {
            let _ = self.tx.send(Command::AssignBiosampleProject {
                guid,
                project_id: chosen,
            });
        }
        if close {
            self.assign_project = None;
        } else {
            self.assign_project = Some((guid, chosen));
        }
    }

    /// The Edit-project modal: name / administrator / description. Save sends `UpdateProject`.
    pub(crate) fn edit_project_modal(&mut self, ctx: &egui::Context) {
        let Some(mut edit) = self.edit_project.clone() else {
            return;
        };

        let mut close = false;
        modal_frame(ctx, "edit_project_modal", 400.0, |ui| {
            ui.label(egui::RichText::new(self.tr("editProject.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(6.0);
            ui.label(self.tr("editProject.name"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.name)
                    .hint_text("name")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            ui.label(self.tr("editProject.admin"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.administrator)
                    .hint_text("administrator")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            ui.label(self.tr("editProject.description"));
            ui.add(
                egui::TextEdit::multiline(&mut edit.description)
                    .hint_text("description (optional)")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(10.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add_enabled(
                        !edit.name.trim().is_empty(),
                        egui::Button::new(self.tr("common.save")).fill(ACCENT),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::UpdateProject {
                        id: edit.id,
                        name: edit.name.trim().to_string(),
                        description: opt(&edit.description),
                        administrator: opt(&edit.administrator).unwrap_or_else(|| "unknown".into()),
                    });
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.edit_project = None;
        } else {
            self.edit_project = Some(edit);
        }
    }

    /// The Delete-project confirmation modal. Confirm sends `DeleteProject`; the app layer
    /// refuses (surfaced via the status bar) while subjects still belong to the project.
    pub(crate) fn delete_project_modal(&mut self, ctx: &egui::Context) {
        let Some((id, name)) = self.confirm_delete_project.clone() else {
            return;
        };

        let mut close = false;
        modal_frame(ctx, "delete_project_modal", 400.0, |ui| {
            ui.label(
                egui::RichText::new(self.tr("editProject.deleteTitle"))
                    .strong()
                    .size(16.0),
            );
            ui.separator();
            ui.add_space(8.0);
            ui.label(format!("{} “{}”?", self.tr("delete.confirm"), name));
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("editProject.deleteNote")).weak().small());
            ui.add_space(12.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE))
                            .fill(DANGER),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::DeleteProject(id));
                    if self.selected_project == Some(id) {
                        self.selected_project = None;
                    }
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.confirm_delete_project = None;
        }
    }

    /// The Edit-run modal: test type (dropdown) + platform / instrument / library layout. Read
    /// metrics are analysis-derived and not editable here. Save sends `UpdateSequenceRun`.
    pub(crate) fn edit_run_modal(&mut self, ctx: &egui::Context) {
        let Some(mut edit) = self.edit_run.clone() else { return };

        let mut close = false;
        modal_frame(ctx, "edit_run_modal", 400.0, |ui| {
            ui.label(egui::RichText::new(self.tr("editRun.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(6.0);
            ui.label(self.tr("form.testType"));
            let current = testtype::display_name(&edit.test_type).to_string();
            egui::ComboBox::from_id_salt("edit_run_test_type")
                .selected_text(current)
                .width(360.0)
                .show_ui(ui, |ui| {
                    for t in testtype::CATALOG {
                        ui.selectable_value(
                            &mut edit.test_type,
                            t.code.to_string(),
                            format!("{}  ·  {}", t.display_name, t.target.label()),
                        );
                    }
                });
            ui.add_space(4.0);
            ui.label(self.tr("editRun.platform"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.platform_name)
                    .hint_text("platform (e.g. ILLUMINA)")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            ui.label(self.tr("editRun.instrument"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.instrument_model)
                    .hint_text("instrument model (optional)")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            ui.label(self.tr("editRun.layout"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.library_layout)
                    .hint_text("library layout (optional, e.g. PAIRED)")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            // Lab / sequencing facility — a dropdown from the labs catalog ("(none)" clears
            // it). Resolved automatically from the instrument id once the AppView lookup
            // ships (roadmap D8); set manually here meanwhile.
            ui.label(self.tr("editRun.lab"));
            let lab_text = if edit.sequencing_facility.is_empty() {
                "(none)".to_string()
            } else {
                edit.sequencing_facility.clone()
            };
            egui::ComboBox::from_id_salt("edit_run_lab")
                .selected_text(lab_text)
                .width(360.0)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut edit.sequencing_facility, String::new(), "(none)");
                    for name in navigator_domain::labs::sequence_run_lab_names() {
                        ui.selectable_value(&mut edit.sequencing_facility, name.to_string(), name);
                    }
                });
            ui.add_space(10.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let ready = testtype::by_code(&edit.test_type).is_some();
                if ui
                    .add_enabled(ready, egui::Button::new(self.tr("common.save")).fill(ACCENT))
                    .clicked()
                {
                    let _ = self.tx.send(Command::UpdateSequenceRun {
                        id: edit.id,
                        biosample_guid: edit.guid,
                        platform_name: edit.platform_name.trim().to_string(),
                        instrument_model: opt(&edit.instrument_model),
                        test_type: edit.test_type.clone(),
                        library_layout: opt(&edit.library_layout),
                        sequencing_facility: opt(&edit.sequencing_facility),
                    });
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.edit_run = None;
        } else {
            self.edit_run = Some(edit);
        }
    }

    /// The Edit-alignment modal: reference build / aligner / variant caller. File paths are
    /// managed by import/probe. Save sends `UpdateAlignment`.
    pub(crate) fn edit_alignment_modal(&mut self, ctx: &egui::Context) {
        let Some(mut edit) = self.edit_alignment.clone() else {
            return;
        };

        let mut close = false;
        modal_frame(ctx, "edit_alignment_modal", 400.0, |ui| {
            ui.label(egui::RichText::new(self.tr("editAln.title")).strong().size(16.0));
            ui.separator();
            ui.add_space(6.0);
            ui.label(self.tr("editAln.build"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.reference_build)
                    .hint_text("reference build")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            ui.label(self.tr("editAln.aligner"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.aligner)
                    .hint_text("aligner")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(4.0);
            ui.label(self.tr("editAln.caller"));
            ui.add(
                egui::TextEdit::singleline(&mut edit.variant_caller)
                    .hint_text("variant caller (optional)")
                    .desired_width(f32::INFINITY),
            );
            ui.add_space(10.0);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let ready = !edit.reference_build.trim().is_empty() && !edit.aligner.trim().is_empty();
                if ui
                    .add_enabled(ready, egui::Button::new(self.tr("common.save")).fill(ACCENT))
                    .clicked()
                {
                    let _ = self.tx.send(Command::UpdateAlignment {
                        id: edit.id,
                        sequence_run_id: edit.run_id,
                        reference_build: edit.reference_build.trim().to_string(),
                        aligner: edit.aligner.trim().to_string(),
                        variant_caller: opt(&edit.variant_caller),
                    });
                    close = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });
        if close {
            self.edit_alignment = None;
        } else {
            self.edit_alignment = Some(edit);
        }
    }

    /// FTDNA import review: a dry-run plan grouped into Needs-confirmation (per-row Merge/New/Skip),
    /// Auto-merged, and New subjects, with a Commit. Deferred dispatch — only `self.tr` is used inside
    /// the closure; resolution edits + the commit are applied after.
    pub(crate) fn ftdna_review_modal(&mut self, ctx: &egui::Context) {
        let Some(plan) = self.ftdna_plan.clone() else { return };
        let (n_new, n_merge, n_confirm) = plan.counts();
        let resolutions = self.ftdna_resolutions.clone();
        let (mut close, mut commit) = (false, false);
        let mut set_res: Vec<(String, FtdnaResolution)> = Vec::new();

        modal_frame(ctx, "ftdna_review", 660.0, |ui| {
            ui.label(egui::RichText::new(self.tr("ftdna.reviewTitle")).strong().size(16.0));
            ui.label(egui::RichText::new(&plan.project_name).weak());
            ui.label(
                egui::RichText::new(format!(
                    "{n_new} {} · {n_merge} {} · {n_confirm} {}",
                    self.tr("ftdna.new"),
                    self.tr("ftdna.autoMerged"),
                    self.tr("ftdna.needsConfirm"),
                ))
                .weak(),
            );
            // Recognized-input diagnostics: a missing roster (0) explains all-orphan/no-match results.
            let s = &plan.stats;
            let roster_txt = egui::RichText::new(format!(
                "{}: {} · {}: {} · {}: {} · {}: {}",
                self.tr("ftdna.statRoster"),
                s.roster,
                self.tr("ftdna.statAncestry"),
                s.paternal + s.maternal,
                self.tr("ftdna.statYstr"),
                s.ystr,
                self.tr("ftdna.statScanned"),
                s.scanned_subjects,
            ))
            .small();
            // Red when no roster was recognized (likely the Member_Information file was not selected).
            if s.roster == 0 {
                ui.label(roster_txt.color(egui::Color32::from_rgb(210, 120, 70)));
                ui.label(
                    egui::RichText::new(self.tr("ftdna.noRosterHint"))
                        .small()
                        .color(egui::Color32::from_rgb(210, 120, 70)),
                );
            } else {
                ui.label(roster_txt.weak());
            }
            ui.separator();
            egui::ScrollArea::vertical().max_height(420.0).show(ui, |ui| {
                // Needs confirmation — the only group with per-row actions.
                if n_confirm > 0 {
                    ui.label(egui::RichText::new(self.tr("ftdna.needsConfirm")).strong());
                    for row in plan.rows.iter() {
                        let MatchKind::NeedsConfirm { candidates } = &row.kind else {
                            continue;
                        };
                        egui::Frame::group(ui.style()).show(ui, |ui| {
                            ui.label(egui::RichText::new(&row.label).strong());
                            if let Some(y) = &row.y_terminal {
                                ui.label(
                                    egui::RichText::new(format!("Y {y} · {} STR", row.ystr_count))
                                        .weak()
                                        .small(),
                                );
                            }
                            // Mutually-exclusive choice per kit: merge into a candidate, new, or skip.
                            // Radio buttons (not selectable labels) so the active choice is obvious.
                            let cur = resolutions.get(&row.kit_number);
                            for c in candidates {
                                let sel = matches!(cur, Some(FtdnaResolution::Merge(g)) if *g == c.guid);
                                let label = format!(
                                    "{} {} ({:.0}% · {})",
                                    self.tr("ftdna.mergeInto"),
                                    c.donor_identifier,
                                    c.score * 100.0,
                                    c.reasons.join(", ")
                                );
                                if ui.radio(sel, label).clicked() {
                                    set_res.push((row.kit_number.clone(), FtdnaResolution::Merge(c.guid)));
                                }
                            }
                            ui.horizontal(|ui| {
                                if ui
                                    .radio(matches!(cur, Some(FtdnaResolution::New)), self.tr("ftdna.itsNew"))
                                    .clicked()
                                {
                                    set_res.push((row.kit_number.clone(), FtdnaResolution::New));
                                }
                                if ui
                                    .radio(matches!(cur, Some(FtdnaResolution::Skip)), self.tr("ftdna.skip"))
                                    .clicked()
                                {
                                    set_res.push((row.kit_number.clone(), FtdnaResolution::Skip));
                                }
                            });
                        });
                    }
                }

                // Auto-merged (exact kit#) — informational.
                if n_merge > 0 {
                    ui.add_space(4.0);
                    ui.collapsing(format!("{} ({n_merge})", self.tr("ftdna.autoMerged")), |ui| {
                        for row in plan.rows.iter() {
                            if let MatchKind::AutoMerge { donor_identifier, .. } = &row.kind {
                                ui.label(format!("{} → {donor_identifier}", row.kit_number));
                            }
                        }
                    });
                }

                // New subjects — informational (orphans flagged).
                ui.add_space(4.0);
                ui.collapsing(format!("{} ({n_new})", self.tr("ftdna.new")), |ui| {
                    for row in plan.rows.iter() {
                        if matches!(row.kind, MatchKind::New) {
                            let badge = if row.in_roster { "" } else { " · orphan" };
                            ui.label(format!("{}{badge}", row.label));
                        }
                    }
                });
            });
            ui.separator();
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .button(egui::RichText::new(self.tr("ftdna.commit")).strong())
                    .clicked()
                {
                    commit = true;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    close = true;
                }
            });
        });

        for (kit, res) in set_res {
            self.ftdna_resolutions.insert(kit, res);
        }
        if commit {
            let resolutions = self.ftdna_resolutions.clone();
            self.status = self.tr("ftdna.committing").to_string();
            let _ = self.tx.send(Command::CommitFtdnaImport { plan, resolutions });
            self.ftdna_plan = None;
        } else if close {
            self.ftdna_plan = None;
            self.ftdna_resolutions.clear();
        }
    }
}

/// Shared modal scaffold: a dimmed full-screen backdrop + a centered `Frame::window` of `width`.
/// `id` namespaces the dim layer + area; `add_contents` draws the modal body.
fn modal_frame(ctx: &egui::Context, id: &str, width: f32, add_contents: impl FnOnce(&mut egui::Ui)) {
    let painter = ctx.layer_painter(egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(format!("{id}_dim")),
    ));
    painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));
    egui::Area::new(egui::Id::new(format!("{id}_modal")))
        .order(egui::Order::Foreground)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            egui::Frame::window(ui.style())
                .inner_margin(egui::Margin::same(18.0))
                .show(ui, |ui| {
                    ui.set_width(width);
                    add_contents(ui);
                });
        });
}
