//! `impl NavigatorApp` methods extracted from `ui.rs` (the `chrome` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    /// Kick off the full-analysis pipeline for an alignment and show the modal immediately.
    pub(crate) fn start_full_analysis(&mut self, alignment_id: i64) {
        self.analysis = Some(AnalysisModal {
            step: 1,
            total: 8,
            label: "Starting".into(),
            detail: "preparing pipeline".into(),
            fraction: 0.0,
            started: self.frame_time,
        });
        self.status = format!("Running full analysis on alignment #{alignment_id}…");
        let _ = self.tx.send(Command::RunFullAnalysis { alignment_id });
    }

    /// Translate a catalog key for the active language. Returns `&'static str` (catalogs are
    /// embedded), so it never borrows `self` — convenient inside egui closures.
    pub(crate) fn tr(&self, key: &'static str) -> &'static str {
        crate::i18n::tr(self.lang, key)
    }

    /// The top app bar: product title (left), settings + language + account controls (right).
    pub(crate) fn app_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("appbar").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                ui.heading(self.tr("app.name"));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("⚙").on_hover_text(self.tr("settings.title")).clicked() {
                        self.settings_form = SettingsForm::from_settings();
                        self.show_settings = true;
                        let _ = self.tx.send(Command::LoadReferenceSettings);
                    }
                    ui.separator();
                    // Theme toggle lives in Settings → Appearance (no duplicate app-bar control).
                    let prev_lang = self.lang;
                    egui::ComboBox::from_id_salt("lang")
                        .selected_text(self.lang.label())
                        .show_ui(ui, |ui| {
                            for &l in crate::i18n::Lang::all() {
                                ui.selectable_value(&mut self.lang, l, l.label());
                            }
                        });
                    if self.lang != prev_lang {
                        crate::i18n::save_lang(self.lang); // persist across restarts
                    }
                    ui.separator();
                    self.account_controls(ui);
                });
            });
            ui.add_space(2.0);
        });
    }

    /// The primary navigation strip (Dashboard / Subjects / Projects).
    pub(crate) fn nav_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("nav").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                for (nav, icon, key) in [
                    (Nav::Dashboard, "📊", "nav.dashboard"),
                    (Nav::Subjects, "👥", "nav.subjects"),
                    (Nav::Projects, "📁", "nav.projects"),
                ] {
                    let label = format!("{icon}  {}", self.tr(key));
                    if ui
                        .selectable_label(self.nav == nav, egui::RichText::new(label).strong())
                        .clicked()
                    {
                        self.nav = nav;
                    }
                }
            });
            ui.add_space(2.0);
        });
    }

    /// Account / sign-in controls (reused in the app bar).
    fn account_controls(&mut self, ui: &mut egui::Ui) {
        match &self.account {
            Some(did) => {
                let did = did.clone();
                if ui.button(self.tr("account.signOut")).clicked() {
                    let _ = self.tx.send(Command::Logout);
                }
                // A local did:key identity self-certifies (no PDS) — show a "local" chip; a real PDS
                // account shows online/offline.
                if did.starts_with("did:key:") {
                    ui.colored_label(egui::Color32::from_rgb(150, 160, 220), self.tr("account.localIdentity"));
                } else if self.online {
                    ui.colored_label(egui::Color32::from_rgb(80, 190, 120), self.tr("account.online"));
                } else {
                    ui.colored_label(egui::Color32::from_rgb(220, 150, 60), self.tr("account.offline"));
                }
                let short: String = did.chars().take(22).collect();
                ui.label(egui::RichText::new(short).weak());
            }
            None => {
                if self.logging_in {
                    ui.spinner();
                }
                let ready = !self.forms.login_handle.trim().is_empty() && !self.logging_in;
                if ui
                    .add_enabled(ready, egui::Button::new(self.tr("account.signIn")))
                    .clicked()
                {
                    self.logging_in = true;
                    self.status = "Opening browser to authorize…".into();
                    let _ = self.tx.send(Command::Login {
                        handle: self.forms.login_handle.trim().to_string(),
                    });
                }
                let hint = self.tr("account.handleHint");
                ui.add(
                    egui::TextEdit::singleline(&mut self.forms.login_handle)
                        .hint_text(hint)
                        .desired_width(180.0),
                );
                ui.label(self.tr("account.pds"));
                // Self-certifying did:key identity — a dev/local-stack affordance (federation only, no
                // PDS repo). Compiled out of distributed (release) builds; opt in for a release dev
                // build with `--features dev-identity`.
                if cfg!(any(debug_assertions, feature = "dev-identity"))
                    && ui
                        .add_enabled(!self.logging_in, egui::Button::new(self.tr("account.useLocal")))
                        .clicked()
                {
                    let _ = self.tx.send(Command::UseLocalIdentity);
                }
            }
        }
    }

    /// The left panel, routed by the active nav tab. Hidden on the Dashboard.
    pub(crate) fn left_panel(&mut self, ctx: &egui::Context) {
        match self.nav {
            Nav::Dashboard => {}
            Nav::Projects => {
                egui::SidePanel::left("left")
                    .min_width(240.0)
                    .show(ctx, |ui| self.projects_side(ui));
            }
            Nav::Subjects if self.subjects_collapsed => {
                // Collapsed: a thin strip with just an expand button, handing the detail panel
                // the full width for charts/tables.
                egui::SidePanel::left("left")
                    .resizable(false)
                    .exact_width(34.0)
                    .show(ctx, |ui| {
                        ui.add_space(6.0);
                        if ui.button("▶").on_hover_text(self.tr("subjects.expand")).clicked() {
                            self.subjects_collapsed = false;
                        }
                    });
            }
            Nav::Subjects => {
                egui::SidePanel::left("left")
                    .resizable(true)
                    .default_width(680.0)
                    .min_width(420.0)
                    .show(ctx, |ui| self.subjects_side(ui));
            }
        }
    }

    /// Project management: the projects list, new-project form, batch import, and panels.
    fn projects_side(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.heading(self.tr("projects.heading"));
        ui.separator();
        let mut pick = None;
        for ov in &self.overview {
            let label = format!("{}  ({})", ov.project.name, ov.sample_count);
            if ui
                .selectable_label(self.selected_project == Some(ov.project.id), label)
                .clicked()
            {
                pick = Some(ov.project.id);
            }
        }
        if let Some(id) = pick {
            self.select_project(id);
        }

        ui.add_space(12.0);
        ui.separator();
        ui.label(self.tr("projects.new"));
        ui.add(egui::TextEdit::singleline(&mut self.forms.project_name).hint_text("name"));
        ui.add(egui::TextEdit::singleline(&mut self.forms.project_admin).hint_text("administrator"));
        if ui
            .add_enabled(
                !self.forms.project_name.trim().is_empty(),
                egui::Button::new(self.tr("projects.create")),
            )
            .clicked()
        {
            let _ = self.tx.send(Command::CreateProject(NewProject {
                name: self.forms.project_name.trim().to_string(),
                description: None,
                administrator: opt(&self.forms.project_admin).unwrap_or_else(|| "unknown".into()),
            }));
            self.forms.project_name.clear();
            self.forms.project_admin.clear();
        }

        ui.add_space(8.0);
        ui.label(self.tr("projects.batchImportHint"));
        if ui
            .add_enabled(!self.importing, egui::Button::new(self.tr("projects.batchImport")))
            .clicked()
        {
            if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                self.importing = true;
                self.reference_needs.clear();
                self.status = format!("Importing {}…", dir.display());
                let _ = self.tx.send(Command::ImportProjectDir { dir, reference: None });
            }
        }
        self.reference_prompt(ui);

        ui.add_space(12.0);
        ui.separator();
        self.panels_section(ui);
    }

    /// Subjects browser: a search box + "Add New Subject" on one row, then the subjects table.
    fn subjects_side(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui.button("◀").on_hover_text(self.tr("subjects.collapse")).clicked() {
                self.subjects_collapsed = true;
            }
            let btn_w = 160.0;
            let hint = self.tr("subjects.search");
            ui.add(
                egui::TextEdit::singleline(&mut self.subject_search)
                    .hint_text(hint)
                    .desired_width((ui.available_width() - btn_w).max(120.0)),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Button::new(egui::RichText::new(self.tr("subjects.addNew")).color(egui::Color32::WHITE))
                            .fill(ACCENT),
                    )
                    .clicked()
                {
                    self.forms.show_add_subject = !self.forms.show_add_subject;
                }
            });
        });
        if self.forms.show_add_subject {
            self.add_subject_form(ui);
        }
        ui.add_space(6.0);
        self.subjects_table(ui);
    }

    /// All biosamples, independent of any project (the project link is optional). Selecting
    /// one drives the runs → alignments → analysis flow in the central panel; adding one
    /// tags it to the open project if there is one, else leaves it project-less.
    /// The subjects table: columns ID / Name / Y-DNA / mtDNA / Sex / Center / Status, with the
    /// selected row highlighted. Clicking a row selects the subject.
    fn subjects_table(&mut self, ui: &mut egui::Ui) {
        if self.all_biosamples.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.none")).weak());
            return;
        }
        let total_w: f32 = SUBJECT_COLS.iter().map(|c| c.1).sum();
        let needle = self.subject_search.trim().to_lowercase();
        let mut pick = None;
        let mut shown = 0;
        egui::ScrollArea::both().show(ui, |ui| {
            ui.set_min_width(total_w);
            table_header(ui, &SUBJECT_COLS);
            for s in &self.all_biosamples {
                if !needle.is_empty() {
                    let hay = format!(
                        "{} {}",
                        s.donor_identifier.to_lowercase(),
                        s.sample_accession.as_deref().unwrap_or("").to_lowercase()
                    );
                    if !hay.contains(&needle) {
                        continue;
                    }
                }
                shown += 1;
                // Y/mt from the bulk per-subject summary; the selected row prefers the freshly
                // loaded consensus (reflects a just-run assignment before the summary reloads).
                let sel = self.selected_sample == Some(s.guid);
                let summary = self.haplo_summary.get(&s.guid);
                let y = sel
                    .then(|| self.consensus_y.as_ref().map(|c| c.haplogroup.clone()))
                    .flatten()
                    .or_else(|| summary.and_then(|(y, _)| y.clone()))
                    .unwrap_or_else(|| "-".into());
                let mt = sel
                    .then(|| self.consensus_mt.as_ref().map(|c| c.haplogroup.clone()))
                    .flatten()
                    .or_else(|| summary.and_then(|(_, m)| m.clone()))
                    .unwrap_or_else(|| "-".into());
                let cells = [
                    short_guid(s),
                    s.donor_identifier.clone(),
                    y,
                    mt,
                    s.sex.clone().unwrap_or_else(|| "-".into()),
                    s.center_name.clone().unwrap_or_else(|| "-".into()),
                    "Pending".to_string(),
                ];
                if table_row(ui, &SUBJECT_COLS, &cells, self.selected_sample == Some(s.guid), Some(6)) {
                    pick = Some(s.guid);
                }
            }
            if shown == 0 {
                ui.label(egui::RichText::new(self.tr("subjects.noMatch")).weak());
            }
        });
        if let Some(guid) = pick {
            self.select_sample(guid);
        }
    }

    /// The inline new-subject form (toggled by the "Add New Subject" button).
    fn add_subject_form(&mut self, ui: &mut egui::Ui) {
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_donor).hint_text("donor identifier"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_accession).hint_text("accession (optional)"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_sex).hint_text("sex (optional)"));
            if let Some(pid) = self.selected_project {
                let name = self
                    .overview
                    .iter()
                    .find(|o| o.project.id == pid)
                    .map(|o| o.project.name.as_str());
                ui.label(format!("→ project: {}", name.unwrap_or("(open)")));
            } else {
                ui.label(self.tr("projects.noProject"));
            }
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(
                        !self.forms.sample_donor.trim().is_empty(),
                        egui::Button::new(self.tr("projects.addSubject")).fill(ACCENT),
                    )
                    .clicked()
                {
                    let _ = self.tx.send(Command::AddBiosample(NewBiosample {
                        project_id: self.selected_project,
                        donor_identifier: self.forms.sample_donor.trim().to_string(),
                        sample_accession: opt(&self.forms.sample_accession),
                        sex: opt(&self.forms.sample_sex),
                    }));
                    self.forms.sample_donor.clear();
                    self.forms.sample_accession.clear();
                    self.forms.sample_sex.clear();
                    self.forms.show_add_subject = false;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    self.forms.show_add_subject = false;
                }
            });
        });
    }
}
