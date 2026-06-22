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

    /// Flip the interface mode (Simple ⇄ Advanced), pin it (so the first-run heuristic stops
    /// overriding), and persist the choice. Leaving Advanced for Simple snaps the nav back to a
    /// Simple-visible tab.
    pub(crate) fn toggle_ui_mode(&mut self) {
        self.ui_mode = match self.ui_mode {
            UiMode::Simple => UiMode::Advanced,
            UiMode::Advanced => UiMode::Simple,
        };
        self.ui_mode_pinned = true;
        if let Err(e) = navigator_app::persist_ui_mode(self.ui_mode) {
            self.status = format!("Could not save interface mode: {e}");
        }
        self.normalize_for_mode();
    }

    /// First-run default: until the user pins a mode, derive it from the workspace — a casual user
    /// (no projects, at most one subject) gets Simple; anyone with projects or multiple subjects
    /// gets Advanced. Re-evaluated as subjects/projects load; a no-op once pinned.
    pub(crate) fn apply_ui_mode_heuristic(&mut self) {
        if self.ui_mode_pinned {
            return;
        }
        self.ui_mode = if self.overview.is_empty() && self.all_biosamples.len() <= 1 {
            UiMode::Simple
        } else {
            UiMode::Advanced
        };
        self.normalize_for_mode();
    }

    /// Keep nav consistent with the mode: Simple hides Projects/Community, so snap off them.
    pub(crate) fn normalize_for_mode(&mut self) {
        if self.ui_mode == UiMode::Simple && matches!(self.nav, Nav::Projects | Nav::Community) {
            self.nav = Nav::Subjects;
        }
    }

    /// In Simple mode, the casual user has (usually) one subject — select it automatically so the
    /// brief/overview appears without a list interaction. Fires once per load (guarded by the
    /// existing selection).
    pub(crate) fn auto_select_single_subject(&mut self) {
        if self.ui_mode == UiMode::Simple && self.selected_sample.is_none() && self.all_biosamples.len() == 1 {
            let guid = self.all_biosamples[0].guid;
            self.select_sample(guid);
        }
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
                    // Interface-mode toggle: Simple (casual briefs) ⇄ Advanced (full UI).
                    let (label, hover) = match self.ui_mode {
                        UiMode::Simple => (self.tr("mode.simple"), self.tr("mode.toAdvanced")),
                        UiMode::Advanced => (self.tr("mode.advanced"), self.tr("mode.toSimple")),
                    };
                    if ui.button(label).on_hover_text(hover).clicked() {
                        self.toggle_ui_mode();
                    }
                    ui.separator();
                    // Notification bell → Community / Notifications. Only meaningful when signed in.
                    if self.account.is_some() {
                        let bell = if self.notif_unread > 0 {
                            format!("🔔 {}", self.notif_unread)
                        } else {
                            "🔔".to_string()
                        };
                        let resp = ui.button(bell).on_hover_text(self.tr("community.notifications"));
                        if resp.clicked() {
                            self.nav = Nav::Community;
                            self.community_tab = CommunityTab::Notifications;
                            let _ = self.tx.send(Command::LoadNotifications);
                        }
                        ui.separator();
                    }
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
                // Simple mode is a single-person experience: only Dashboard + the subject view
                // ("My DNA"). Advanced exposes the full workspace.
                let tabs: &[(Nav, &str, &str)] = match self.ui_mode {
                    UiMode::Simple => &[
                        (Nav::Dashboard, "📊", "nav.dashboard"),
                        (Nav::Subjects, "🧬", "nav.myDna"),
                    ],
                    UiMode::Advanced => &[
                        (Nav::Dashboard, "📊", "nav.dashboard"),
                        (Nav::Subjects, "👥", "nav.subjects"),
                        (Nav::Projects, "📁", "nav.projects"),
                        (Nav::Community, "💬", "nav.community"),
                    ],
                };
                for &(nav, icon, key) in tabs {
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
            // Dashboard + Community are full-width (no side panel).
            Nav::Dashboard | Nav::Community => {}
            Nav::Projects if self.projects_collapsed => {
                // Collapsed: a thin strip with just an expand button, handing the detail panel
                // the full width for the wide Y-STR chart.
                egui::SidePanel::left("left")
                    .resizable(false)
                    .exact_width(34.0)
                    .show(ctx, |ui| {
                        ui.add_space(6.0);
                        if ui.button("▶").on_hover_text(self.tr("projects.expand")).clicked() {
                            self.projects_collapsed = false;
                        }
                    });
            }
            Nav::Projects => {
                egui::SidePanel::left("left")
                    .resizable(true)
                    .default_width(300.0)
                    .min_width(240.0)
                    .show(ctx, |ui| self.projects_side(ui));
            }
            // Simple mode: a single-subject experience. With one subject (the common case) the
            // side panel is hidden entirely — the brief/overview fills the window; with several,
            // a minimal "who am I looking at" selector.
            Nav::Subjects if self.ui_mode == UiMode::Simple => {
                if self.all_biosamples.len() > 1 {
                    egui::SidePanel::left("left")
                        .resizable(true)
                        .default_width(260.0)
                        .min_width(200.0)
                        .show(ctx, |ui| self.simple_subjects_side(ui));
                }
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
        ui.horizontal(|ui| {
            if ui.button("◀").on_hover_text(self.tr("projects.collapse")).clicked() {
                self.projects_collapsed = true;
            }
            ui.heading(self.tr("projects.heading"));
        });
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

        // FTDNA project import — imports into the selected project, or creates one named from the
        // exports if none is selected.
        ui.add_space(8.0);
        ui.label(self.tr("ftdna.importHint"));
        let hover = if self.selected_project.is_some() {
            self.tr("ftdna.intoSelected")
        } else {
            self.tr("ftdna.intoNew")
        };
        if ui
            .add_enabled(!self.importing, egui::Button::new(self.tr("ftdna.import")))
            .on_hover_text(hover)
            .clicked()
        {
            if let Some(paths) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).pick_files() {
                self.start_ftdna_import(paths);
            }
        }

        ui.add_space(12.0);
        ui.separator();
        self.panels_section(ui);
    }

    /// Classify the picked files by header sniff, route each to its FTDNA parser slot, and dispatch a
    /// dry-run plan against the open project. Unrecognized files are ignored (noted in the status).
    fn start_ftdna_import(&mut self, paths: Vec<PathBuf>) {
        use navigator_domain::ftdna::{classify, FtdnaFileKind};
        let (mut member, mut paternal, mut maternal, mut ystr) = (None, None, None, None);
        let mut unrecognized = 0;
        for p in paths {
            match std::fs::read_to_string(&p).ok().as_deref().and_then(classify) {
                Some(FtdnaFileKind::Member) => member = Some(p),
                Some(FtdnaFileKind::PaternalAncestry) => paternal = Some(p),
                Some(FtdnaFileKind::MaternalAncestry) => maternal = Some(p),
                Some(FtdnaFileKind::YdnaOverview) => ystr = Some(p),
                None => unrecognized += 1,
            }
        }
        if member.is_none() && paternal.is_none() && maternal.is_none() && ystr.is_none() {
            self.status = self.tr("ftdna.noneRecognized").to_string();
            return;
        }
        // No project selected → import into a new one named from the FTDNA filenames (FTDNA prefixes
        // each export with the project name).
        let project_name = self.selected_project.is_none().then(|| {
            [&member, &paternal, &maternal, &ystr]
                .into_iter()
                .flatten()
                .find_map(|p| ftdna_project_name(p))
                .unwrap_or_else(|| "FTDNA Project".to_string())
        });
        self.importing = true;
        self.status = self.tr("ftdna.planning").to_string();
        if unrecognized > 0 {
            self.status = format!("{} ({} {})", self.status, unrecognized, self.tr("ftdna.unrecognized"));
        }
        let _ = self.tx.send(Command::PlanFtdnaImport {
            project_id: self.selected_project,
            project_name,
            member,
            paternal,
            maternal,
            ystr,
        });
    }

    /// Simple-mode subject selector: a "who am I looking at" list with a free-text filter and a
    /// scrollable, row-virtualized body (a research surface can hold thousands of subjects), plus the
    /// Add-New affordance. Only shown when more than one subject exists.
    fn simple_subjects_side(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);
        ui.heading(self.tr("nav.myDna"));
        ui.separator();

        // Add New (kept above the scroll body so it's always reachable).
        if ui
            .add(
                egui::Button::new(egui::RichText::new(self.tr("subjects.addNew")).color(egui::Color32::WHITE))
                    .fill(ACCENT),
            )
            .clicked()
        {
            self.forms.show_add_subject = !self.forms.show_add_subject;
        }
        if self.forms.show_add_subject {
            self.add_subject_form(ui);
        }

        ui.add_space(6.0);
        let filter_hint = self.tr("subjects.filter");
        ui.add(
            egui::TextEdit::singleline(&mut self.simple_subject_filter)
                .hint_text(filter_hint)
                .desired_width(f32::INFINITY),
        );

        // Filtered view built from immutable reads first, so the scroll closure borrows only locals.
        let needle = self.simple_subject_filter.trim().to_lowercase();
        let rows: Vec<(SampleGuid, String)> = self
            .all_biosamples
            .iter()
            .filter(|b| needle.is_empty() || b.donor_identifier.to_lowercase().contains(&needle))
            .map(|b| (b.guid, format!("👤  {}", b.donor_identifier)))
            .collect();

        ui.add_space(4.0);
        ui.label(egui::RichText::new(format!("{}", rows.len())).weak().small());
        ui.separator();

        if rows.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.noMatch")).weak());
            return;
        }

        let selected = self.selected_sample;
        let mut pick = None;
        let row_h = ui.spacing().interact_size.y;
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show_rows(ui, row_h, rows.len(), |ui, range| {
                for i in range {
                    let (guid, label) = &rows[i];
                    if ui.selectable_label(selected == Some(*guid), label).clicked() {
                        pick = Some(*guid);
                    }
                }
            });
        if let Some(guid) = pick {
            self.select_sample(guid);
        }
    }

    /// Subjects browser: a collapse toggle + "Add New Subject" on one row, then the subjects table
    /// (sort + filter live in the table's sticky header).
    fn subjects_side(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui.button("◀").on_hover_text(self.tr("subjects.collapse")).clicked() {
                self.subjects_collapsed = true;
            }
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
        use egui_extras::{Column, TableBuilder};

        if self.all_biosamples.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.none")).weak());
            return;
        }

        // Build the display rows from immutable `self` reads first, so the render closures below
        // can borrow only the table-control field (and locals) without conflicting.
        struct Row {
            guid: SampleGuid,
            cells: [String; 6],
        }
        let mut rows: Vec<Row> = self
            .all_biosamples
            .iter()
            .map(|s| {
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
                Row {
                    guid: s.guid,
                    cells: [
                        s.donor_identifier.clone(),
                        y,
                        mt,
                        s.sex.clone().unwrap_or_else(|| "-".into()),
                        s.center_name.clone().unwrap_or_else(|| "-".into()),
                        "Pending".to_string(),
                    ],
                }
            })
            .collect();

        // Inline per-column filters (AND across columns), then natural-sort by the active column.
        for col in 0..SUBJECT_COLS.len() {
            let f = self.subjects_table_ctl.filter_norm(col);
            if !f.is_empty() {
                rows.retain(|r| r.cells[col].to_lowercase().contains(&f));
            }
        }
        if let Some(c) = self.subjects_table_ctl.sort_col() {
            let asc = self.subjects_table_ctl.ascending();
            rows.sort_by(|a, b| {
                let o = natural_cmp(&a.cells[c], &b.cells[c]);
                if asc {
                    o
                } else {
                    o.reverse()
                }
            });
        }

        let selected = self.selected_sample;
        let labels = SUBJECT_COLS.map(|(l, _)| l);
        let ctl = &mut self.subjects_table_ctl;
        let mut pick: Option<SampleGuid> = None;

        let mut tb = TableBuilder::new(ui)
            .striped(true)
            .sense(egui::Sense::click())
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
            .auto_shrink([false, false]);
        for (_, w) in SUBJECT_COLS {
            tb = tb.column(Column::initial(w).at_least(56.0).clip(true).resizable(true));
        }
        tb.header(44.0, |mut header| {
            for (i, label) in labels.into_iter().enumerate() {
                header.col(|ui| sortable_header(ui, ctl, i, label, true));
            }
        })
        .body(|body| {
            body.rows(26.0, rows.len(), |mut row| {
                let r = &rows[row.index()];
                row.set_selected(Some(r.guid) == selected);
                for (ci, cell) in r.cells.iter().enumerate() {
                    row.col(|ui| {
                        // The "Status" column is rendered as a muted accent badge.
                        if ci == 5 && cell != "-" {
                            chip(
                                ui,
                                cell,
                                egui::Color32::from_rgb(70, 58, 28),
                                egui::Color32::from_rgb(225, 190, 90),
                            );
                        } else {
                            ui.label(cell);
                        }
                    });
                }
                if row.response().clicked() {
                    pick = Some(r.guid);
                }
            });
        });

        if rows.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.noMatch")).weak());
        }
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

/// Derive the FTDNA project name from an export filename. FTDNA prefixes each file with the project
/// name, e.g. `R1b-CTS4466Plus_Member_Information_20260606.csv` → `R1b-CTS4466Plus`.
fn ftdna_project_name(path: &std::path::Path) -> Option<String> {
    let stem = path.file_name()?.to_str()?;
    for marker in [
        "_Member_Information",
        "_Paternal_Ancestry",
        "_Maternal_Ancestry",
        "_YDNA",
    ] {
        if let Some(idx) = stem.find(marker) {
            let name = stem[..idx].trim();
            if !name.is_empty() {
                return Some(name.to_string());
            }
        }
    }
    None
}
