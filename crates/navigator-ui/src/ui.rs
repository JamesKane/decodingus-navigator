//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use eframe::egui;
use navigator_app::{Coverage, DenovoCall, IbdComparison, PanelGenotype, ProjectOverview};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::strprofile::{self, StrProfile};
use navigator_domain::testtype;
use navigator_domain::variants::VariantSet;
use navigator_domain::workspace::{Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, SequenceRun};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event, NewBiosample, PanelInfo};

#[derive(Default)]
struct Forms {
    project_name: String,
    project_admin: String,
    sample_donor: String,
    sample_accession: String,
    sample_sex: String,
    run_platform: String,
    run_test_type: String,
    aln_reference_build: String,
    aln_aligner: String,
    aln_bam: String,
    aln_reference: String,
    denovo_contig: String,
    ploidy: String,
    panel_import_name: String,
    login_handle: String,
    str_panel: String,
    str_provider: String,
    str_source: String,
}

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    overview: Vec<ProjectOverview>,
    selected_project: Option<i64>,
    samples: Vec<Biosample>,
    /// Every biosample (the project-independent subjects list).
    all_biosamples: Vec<Biosample>,
    selected_sample: Option<SampleGuid>,
    runs: Vec<SequenceRun>,
    /// STR profiles for the selected subject.
    str_profiles: Vec<StrProfile>,
    /// SNP variant sets for the selected subject.
    variant_sets: Vec<VariantSet>,
    selected_run: Option<i64>,
    alignments: Vec<Alignment>,
    selected_alignment: Option<i64>,
    coverage: Option<Coverage>,
    running: bool,
    denovo: Option<Vec<DenovoCall>>,
    running_denovo: bool,
    panels: Vec<PanelInfo>,
    selected_panel: Option<i64>,
    all_alignments: Vec<Alignment>,
    panel_genotypes: Option<Vec<PanelGenotype>>,
    running_genotype: bool,
    ibd_other: Option<i64>,
    ibd_result: Option<IbdComparison>,
    running_ibd: bool,
    /// Signed-in account DID, or `None`. Gates the "Publish" actions.
    account: Option<String>,
    /// Whether the last PDS write reached the server (offline indicator).
    online: bool,
    logging_in: bool,
    publishing: bool,
    forms: Forms,
    status: String,
}

fn opt(s: &str) -> Option<String> {
    let t = s.trim();
    (!t.is_empty()).then(|| t.to_string())
}

/// A labeled dropdown that sets `value` to one of `options` (string codes).
fn combo(ui: &mut egui::Ui, label: &str, id: &str, value: &mut String, options: &[&str]) {
    ui.horizontal(|ui| {
        ui.label(label);
        egui::ComboBox::from_id_salt(id).selected_text(value.clone()).show_ui(ui, |ui| {
            for opt in options {
                ui.selectable_value(value, opt.to_string(), *opt);
            }
        });
    });
}

impl NavigatorApp {
    pub fn new(cc: &eframe::CreationContext<'_>, db_path: PathBuf) -> Self {
        let ctx = cc.egui_ctx.clone();
        let (tx, rx) = worker::spawn(db_path, move || ctx.request_repaint());
        let _ = tx.send(Command::LoadOverview);
        let _ = tx.send(Command::LoadAllBiosamples);
        let _ = tx.send(Command::LoadPanels);
        let _ = tx.send(Command::LoadAllAlignments);
        let _ = tx.send(Command::AuthStatus);
        let _ = tx.send(Command::SyncStatus);
        NavigatorApp {
            tx,
            rx,
            overview: Vec::new(),
            selected_project: None,
            samples: Vec::new(),
            all_biosamples: Vec::new(),
            selected_sample: None,
            runs: Vec::new(),
            str_profiles: Vec::new(),
            variant_sets: Vec::new(),
            selected_run: None,
            alignments: Vec::new(),
            selected_alignment: None,
            coverage: None,
            running: false,
            denovo: None,
            running_denovo: false,
            panels: Vec::new(),
            selected_panel: None,
            all_alignments: Vec::new(),
            panel_genotypes: None,
            running_genotype: false,
            ibd_other: None,
            ibd_result: None,
            running_ibd: false,
            account: None,
            online: true,
            logging_in: false,
            publishing: false,
            forms: Forms {
                denovo_contig: "chrM".into(),
                ploidy: "2".into(),
                run_test_type: "WGS".into(),
                str_panel: "Y-37".into(),
                str_provider: "FTDNA".into(),
                str_source: "DIRECT_TEST".into(),
                ..Forms::default()
            },
            status: "Loading…".into(),
        }
    }

    fn drain_events(&mut self) {
        while let Ok(event) = self.rx.try_recv() {
            match event {
                Event::Overview(v) => {
                    self.status = format!("{} project(s)", v.len());
                    self.overview = v;
                }
                Event::ProjectCreated(p) => {
                    self.select_project(p.id);
                    let _ = self.tx.send(Command::LoadOverview);
                }
                Event::Samples { project_id, samples } => {
                    if self.selected_project == Some(project_id) {
                        self.samples = samples;
                    }
                }
                Event::AllBiosamples(v) => self.all_biosamples = v,
                Event::BiosamplesChanged => {
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadSamples(pid));
                    }
                    let _ = self.tx.send(Command::LoadOverview); // project sample counts changed
                }
                Event::Runs { biosample_guid, runs } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.runs = runs;
                    }
                }
                Event::RunsChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadRuns(guid));
                    }
                }
                Event::StrProfiles { biosample_guid, profiles } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.str_profiles = profiles;
                    }
                }
                Event::StrProfilesChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadStrProfiles(guid));
                    }
                    self.status = "STR profile imported".into();
                }
                Event::VariantSets { biosample_guid, sets } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.variant_sets = sets;
                    }
                }
                Event::VariantSetsChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadVariantSets(guid));
                    }
                    self.status = "Variants imported".into();
                }
                Event::Alignments { sequence_run_id, alignments } => {
                    if self.selected_run == Some(sequence_run_id) {
                        self.alignments = alignments;
                    }
                }
                Event::AlignmentsChanged(run_id) => {
                    if self.selected_run == Some(run_id) {
                        let _ = self.tx.send(Command::LoadAlignments(run_id));
                    }
                    let _ = self.tx.send(Command::LoadAllAlignments); // keep IBD picker current
                }
                Event::Coverage { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.coverage = result;
                    }
                    self.running = false;
                }
                Event::Denovo { alignment_id, contig, result } => {
                    if self.selected_alignment == Some(alignment_id) && self.forms.denovo_contig == contig {
                        self.denovo = result;
                    }
                    self.running_denovo = false;
                }
                Event::Panels(p) => self.panels = p,
                Event::PanelImported => {
                    self.status = "Panel imported".into();
                    let _ = self.tx.send(Command::LoadPanels);
                }
                Event::AllAlignments(a) => self.all_alignments = a,
                Event::PanelGenotypes { alignment_id, panel_id, ploidy, genotypes } => {
                    if self.selected_alignment == Some(alignment_id)
                        && self.selected_panel == Some(panel_id)
                        && self.ploidy() == ploidy
                    {
                        self.panel_genotypes = (!genotypes.is_empty()).then_some(genotypes);
                    }
                    self.running_genotype = false;
                }
                Event::Ibd(cmp) => {
                    self.ibd_result = Some(cmp);
                    self.running_ibd = false;
                }
                Event::Authenticated(account) => {
                    self.status = match &account {
                        Some(did) => format!("Signed in as {did}"),
                        None => "Signed out".into(),
                    };
                    self.account = account;
                    self.logging_in = false;
                }
                Event::Published { kind, uri } => {
                    self.status = format!("Published {kind}: {uri}");
                    self.publishing = false;
                    let _ = self.tx.send(Command::SyncStatus); // refresh the online dot
                }
                Event::SyncOnline(online) => self.online = online,
                Event::Error(e) => {
                    self.status = format!("Error: {e}");
                    self.running = false;
                    self.running_denovo = false;
                    self.running_genotype = false;
                    self.running_ibd = false;
                    self.logging_in = false;
                    self.publishing = false;
                    let _ = self.tx.send(Command::SyncStatus); // a failed publish may have gone offline
                }
            }
        }
    }

    fn select_project(&mut self, id: i64) {
        self.selected_project = Some(id);
        self.samples.clear();
        self.clear_sample_selection();
        let _ = self.tx.send(Command::LoadSamples(id));
    }

    fn select_sample(&mut self, guid: SampleGuid) {
        self.selected_sample = Some(guid);
        self.clear_run_selection();
        self.runs.clear();
        self.str_profiles.clear();
        self.variant_sets.clear();
        let _ = self.tx.send(Command::LoadRuns(guid));
        let _ = self.tx.send(Command::LoadStrProfiles(guid));
        let _ = self.tx.send(Command::LoadVariantSets(guid));
    }

    fn select_run(&mut self, id: i64) {
        self.selected_run = Some(id);
        self.alignments.clear();
        self.selected_alignment = None;
        self.coverage = None;
        let _ = self.tx.send(Command::LoadAlignments(id));
    }

    fn select_alignment(&mut self, id: i64) {
        self.selected_alignment = Some(id);
        self.coverage = None;
        self.denovo = None;
        self.panel_genotypes = None;
        self.ibd_result = None;
        let _ = self.tx.send(Command::LoadCoverage(id));
        let _ = self.tx.send(Command::LoadDenovo { alignment_id: id, contig: self.forms.denovo_contig.clone() });
        if let Some(panel_id) = self.selected_panel {
            let _ = self.tx.send(Command::LoadPanelGenotypes { alignment_id: id, panel_id, ploidy: self.ploidy() });
        }
    }

    fn select_panel(&mut self, panel_id: i64) {
        self.selected_panel = Some(panel_id);
        self.panel_genotypes = None;
        self.ibd_result = None;
        if let Some(aln) = self.selected_alignment {
            let _ = self.tx.send(Command::LoadPanelGenotypes { alignment_id: aln, panel_id, ploidy: self.ploidy() });
        }
    }

    fn ploidy(&self) -> u8 {
        self.forms.ploidy.trim().parse().unwrap_or(2)
    }

    fn clear_sample_selection(&mut self) {
        self.selected_sample = None;
        self.runs.clear();
        self.clear_run_selection();
    }

    fn clear_run_selection(&mut self) {
        self.selected_run = None;
        self.alignments.clear();
        self.selected_alignment = None;
        self.coverage = None;
    }
}

impl eframe::App for NavigatorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.drain_events();
        self.account_panel(ctx);
        self.projects_panel(ctx);
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Status:");
                ui.label(&self.status);
            });
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            if self.selected_project.is_none() && self.selected_sample.is_none() {
                ui.heading("DUNavigator");
                ui.label("Select a subject, or create/select a project.");
                return;
            }
            egui::ScrollArea::vertical().show(ui, |ui| {
                // A project (when open) shows its own samples as a filtered subview.
                if self.selected_project.is_some() {
                    self.samples_section(ui);
                }
                if let Some(guid) = self.selected_sample {
                    self.subject_header(ui);
                    self.str_section(ui, guid);
                    self.variants_section(ui, guid);
                    self.runs_section(ui);
                }
                if self.selected_run.is_some() {
                    self.alignments_section(ui);
                }
                if let Some(id) = self.selected_alignment {
                    self.coverage_section(ui, id);
                    self.denovo_section(ui, id);
                    self.genotyping_section(ui, id);
                }
            });
        });
    }
}

impl NavigatorApp {
    fn account_panel(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("account").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.strong("DUNavigator");
                ui.separator();
                match &self.account {
                    Some(did) => {
                        ui.label(format!("Signed in: {did}"));
                        if self.online {
                            ui.colored_label(egui::Color32::from_rgb(60, 160, 60), "● online");
                        } else {
                            ui.colored_label(egui::Color32::from_rgb(200, 120, 40), "○ offline");
                        }
                        if ui.button("Sign out").clicked() {
                            let _ = self.tx.send(Command::Logout);
                        }
                    }
                    None => {
                        ui.label("PDS:");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.forms.login_handle)
                                .hint_text("handle or did")
                                .desired_width(200.0),
                        );
                        let ready = !self.forms.login_handle.trim().is_empty() && !self.logging_in;
                        if ui.add_enabled(ready, egui::Button::new("Sign in")).clicked() {
                            self.logging_in = true;
                            self.status = "Opening browser to authorize…".into();
                            let _ = self.tx.send(Command::Login {
                                handle: self.forms.login_handle.trim().to_string(),
                            });
                        }
                        if self.logging_in {
                            ui.spinner();
                        }
                    }
                }
            });
        });
    }

    fn projects_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("projects").min_width(220.0).show(ctx, |ui| {
            ui.heading("Projects");
            ui.separator();
            let mut pick = None;
            for ov in &self.overview {
                let label = format!("{}  ({})", ov.project.name, ov.sample_count);
                if ui.selectable_label(self.selected_project == Some(ov.project.id), label).clicked() {
                    pick = Some(ov.project.id);
                }
            }
            if let Some(id) = pick {
                self.select_project(id);
            }

            ui.add_space(12.0);
            ui.separator();
            ui.label("New project");
            ui.add(egui::TextEdit::singleline(&mut self.forms.project_name).hint_text("name"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.project_admin).hint_text("administrator"));
            if ui
                .add_enabled(!self.forms.project_name.trim().is_empty(), egui::Button::new("Create project"))
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

            ui.add_space(12.0);
            ui.separator();
            self.subjects_section(ui);

            ui.add_space(12.0);
            ui.separator();
            self.panels_section(ui);
        });
    }

    /// All biosamples, independent of any project (the project link is optional). Selecting
    /// one drives the runs → alignments → analysis flow in the central panel; adding one
    /// tags it to the open project if there is one, else leaves it project-less.
    fn subjects_section(&mut self, ui: &mut egui::Ui) {
        ui.label("Subjects");
        if self.all_biosamples.is_empty() {
            ui.label("No subjects yet.");
        }
        let mut pick = None;
        for s in &self.all_biosamples {
            let tag = if s.project_id.is_some() { "" } else { "  ·" }; // mark project-less
            let label = format!("{}{}", s.donor_identifier, tag);
            if ui.selectable_label(self.selected_sample == Some(s.guid), label).clicked() {
                pick = Some(s.guid);
            }
        }
        if let Some(guid) = pick {
            self.select_sample(guid);
        }

        ui.collapsing("Add subject", |ui| {
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_donor).hint_text("donor identifier"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_accession).hint_text("accession (optional)"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_sex).hint_text("sex (optional)"));
            if let Some(pid) = self.selected_project {
                let name = self.overview.iter().find(|o| o.project.id == pid).map(|o| o.project.name.as_str());
                ui.label(format!("→ project: {}", name.unwrap_or("(open)")));
            } else {
                ui.label("→ no project");
            }
            if ui
                .add_enabled(!self.forms.sample_donor.trim().is_empty(), egui::Button::new("Add subject"))
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
            }
        });
    }

    /// Heading naming the selected subject above its runs.
    fn subject_header(&mut self, ui: &mut egui::Ui) {
        let Some(guid) = self.selected_sample else { return };
        let donor = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|s| s.guid == guid)
            .map(|s| s.donor_identifier.clone())
            .unwrap_or_else(|| "subject".into());
        ui.add_space(12.0);
        ui.separator();
        ui.heading(format!("Subject — {donor}"));
    }

    /// Y-STR profiles for the selected subject + an import form (CSV/TSV marker table).
    fn str_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("STR profiles");
        if self.str_profiles.is_empty() {
            ui.label("No STR profiles yet.");
        }
        for p in &self.str_profiles {
            let provider = p.provider.as_deref().unwrap_or("—");
            ui.label(format!("{} — {} markers  ({provider})", p.panel_name, p.markers.len()));
        }

        ui.add_space(6.0);
        ui.collapsing("Import STR profile", |ui| {
            combo(ui, "Panel:", "str_panel", &mut self.forms.str_panel, strprofile::KNOWN_PANELS);
            combo(ui, "Provider:", "str_provider", &mut self.forms.str_provider, strprofile::KNOWN_PROVIDERS);
            combo(ui, "Source:", "str_source", &mut self.forms.str_source, strprofile::KNOWN_SOURCES);
            if ui.button("Choose CSV/TSV…").clicked() {
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
            ui.label("Expects rows of marker,value (e.g. DYS393,13).");
        });
    }

    /// SNP variant sets for the selected subject + an import form (VCF or CSV/TSV).
    fn variants_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("SNP variants");
        if self.variant_sets.is_empty() {
            ui.label("No variants imported yet.");
        }
        for s in &self.variant_sets {
            ui.label(format!("{} — {} SNP(s)", s.source_label, s.calls.len()));
        }

        ui.add_space(6.0);
        ui.collapsing("Import variants", |ui| {
            if ui.button("Choose VCF/CSV/TSV…").clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("variants", &["vcf", "csv", "tsv", "txt"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportVariants { biosample_guid: guid, path });
                }
            }
            ui.label("VCF, or rows of contig,position,ref,alt[,rsid][,genotype]. SNP-only.");
        });
    }

    fn panels_section(&mut self, ui: &mut egui::Ui) {
        ui.label("Panels");
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
            .add_enabled(!self.forms.panel_import_name.trim().is_empty(), egui::Button::new("Import sites VCF…"))
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

    fn samples_section(&mut self, ui: &mut egui::Ui) {
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
            ui.label("No samples yet.");
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
        ui.label("Add subjects from the sidebar (they tag to this open project).");
    }

    fn runs_section(&mut self, ui: &mut egui::Ui) {
        let guid = self.selected_sample.unwrap();
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Tests");
        if self.runs.is_empty() {
            ui.label("No tests yet.");
        }
        let mut pick = None;
        for r in &self.runs {
            let platform = if r.platform_name.is_empty() || r.platform_name == "UNKNOWN" {
                String::new()
            } else {
                format!("  ({})", r.platform_name)
            };
            let label = format!("#{}  {}{}", r.id, testtype::display_name(&r.test_type), platform);
            if ui.selectable_label(self.selected_run == Some(r.id), label).clicked() {
                pick = Some(r.id);
            }
        }
        if let Some(id) = pick {
            self.select_run(id);
        }

        ui.add_space(6.0);
        ui.collapsing("Add test", |ui| {
            ui.horizontal(|ui| {
                ui.label("Test type:");
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
            if ui.add_enabled(ready, egui::Button::new("Add test")).clicked() {
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
                // keep the selected test type for adding several of the same kind
            }
        });
    }

    fn alignments_section(&mut self, ui: &mut egui::Ui) {
        let run_id = self.selected_run.unwrap();
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Alignments");
        if self.alignments.is_empty() {
            ui.label("No alignments yet.");
        }
        let mut pick = None;
        for a in &self.alignments {
            let files = if a.bam_path.is_some() { "" } else { "  (no files)" };
            let label = format!("#{}  {} · {}{}", a.id, a.reference_build, a.aligner, files);
            if ui.selectable_label(self.selected_alignment == Some(a.id), label).clicked() {
                pick = Some(a.id);
            }
        }
        if let Some(id) = pick {
            self.select_alignment(id);
        }

        ui.add_space(6.0);
        ui.collapsing("Add alignment", |ui| {
            ui.add(egui::TextEdit::singleline(&mut self.forms.aln_reference_build).hint_text("reference build (e.g. chm13v2.0)"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.aln_aligner).hint_text("aligner (e.g. bwa-mem)"));
            ui.horizontal(|ui| {
                if ui.button("BAM/CRAM…").clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("alignment", &["bam", "cram"]).pick_file() {
                        self.forms.aln_bam = p.to_string_lossy().into_owned();
                    }
                }
                ui.label(if self.forms.aln_bam.is_empty() { "—" } else { self.forms.aln_bam.as_str() });
            });
            ui.horizontal(|ui| {
                if ui.button("Reference…").clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("fasta", &["fa", "fasta", "fna"]).pick_file() {
                        self.forms.aln_reference = p.to_string_lossy().into_owned();
                    }
                }
                ui.label(if self.forms.aln_reference.is_empty() { "—" } else { self.forms.aln_reference.as_str() });
            });
            let ready = !self.forms.aln_reference_build.trim().is_empty()
                && !self.forms.aln_aligner.trim().is_empty()
                && !self.forms.aln_bam.is_empty()
                && !self.forms.aln_reference.is_empty();
            if ui.add_enabled(ready, egui::Button::new("Add alignment")).clicked() {
                let _ = self.tx.send(Command::AddAlignment(NewAlignment {
                    sequence_run_id: run_id,
                    reference_build: self.forms.aln_reference_build.trim().to_string(),
                    aligner: self.forms.aln_aligner.trim().to_string(),
                    variant_caller: None,
                    bam_path: opt(&self.forms.aln_bam),
                    reference_path: opt(&self.forms.aln_reference),
                }));
                self.forms.aln_reference_build.clear();
                self.forms.aln_aligner.clear();
                self.forms.aln_bam.clear();
                self.forms.aln_reference.clear();
            }
        });
    }

    fn coverage_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Coverage");

        let has_paths = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some() && a.reference_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_paths && !self.running, egui::Button::new("Run coverage")).clicked() {
                self.running = true;
                self.status = format!("Running coverage on alignment #{alignment_id}…");
                let _ = self.tx.send(Command::RunCoverage(alignment_id));
            }
            if self.running {
                ui.spinner();
            }
            if !has_paths {
                ui.label("(no BAM/reference path recorded)");
            }
        });

        match &self.coverage {
            None if !self.running => {
                ui.label("No coverage result yet.");
            }
            None => {}
            Some(c) => {
                egui::Grid::new("coverage_metrics").striped(true).num_columns(2).show(ui, |ui| {
                    let row = |ui: &mut egui::Ui, k: &str, v: String| {
                        ui.label(k);
                        ui.label(v);
                        ui.end_row();
                    };
                    row(ui, "Genome territory", c.genome_territory.to_string());
                    row(ui, "Mean coverage", format!("{:.2}", c.mean_coverage));
                    row(ui, "Median coverage", format!("{:.0}", c.median_coverage));
                    row(ui, "Callable bases", c.callable_bases.to_string());
                    row(ui, "% ≥10x", format!("{:.1}%", c.pct_10x * 100.0));
                    row(ui, "% ≥20x", format!("{:.1}%", c.pct_20x * 100.0));
                    row(ui, "% ≥30x", format!("{:.1}%", c.pct_30x * 100.0));
                });
            }
        }

        if self.coverage.is_some() {
            self.publish_row(ui, "Publish summary to PDS", Command::PublishCoverage(alignment_id));
        }
    }

    /// A "Publish to PDS" button + sign-in hint/spinner, shared by the result sections.
    fn publish_row(&mut self, ui: &mut egui::Ui, label: &str, cmd: Command) {
        ui.horizontal(|ui| {
            let ready = self.account.is_some() && !self.publishing;
            if ui.add_enabled(ready, egui::Button::new(label)).clicked() {
                self.publishing = true;
                self.status = "Publishing to PDS…".into();
                let _ = self.tx.send(cmd);
            }
            if self.account.is_none() {
                ui.label("(sign in to publish)");
            }
            if self.publishing {
                ui.spinner();
            }
        });
    }

    fn genotyping_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Panel genotyping & IBD");

        let Some(panel_id) = self.selected_panel else {
            ui.label("Select a panel in the sidebar.");
            return;
        };
        let panel_name = self
            .panels
            .iter()
            .find(|p| p.panel.id == panel_id)
            .map(|p| p.panel.name.clone())
            .unwrap_or_default();
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            ui.label("Ploidy:");
            ui.add(egui::TextEdit::singleline(&mut self.forms.ploidy).desired_width(32.0));
            if ui
                .add_enabled(has_bam && !self.running_genotype, egui::Button::new(format!("Genotype vs {panel_name}")))
                .clicked()
            {
                self.running_genotype = true;
                let _ = self.tx.send(Command::GenotypePanel { alignment_id, panel_id, ploidy: self.ploidy() });
            }
            if self.running_genotype {
                ui.spinner();
            }
        });

        match &self.panel_genotypes {
            Some(genos) => {
                let (mut hr, mut het, mut ha, mut nc) = (0, 0, 0, 0);
                for g in genos {
                    match g.dosage {
                        0 => hr += 1,
                        1 => het += 1,
                        2 => ha += 1,
                        _ => nc += 1,
                    }
                }
                ui.label(format!(
                    "{} sites — {hr} hom-ref, {het} het, {ha} hom-alt, {nc} no-call",
                    genos.len()
                ));
            }
            None if !self.running_genotype => {
                ui.label("Not genotyped against this panel.");
            }
            None => {}
        }

        // IBD compare against another genotyped alignment.
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            ui.label("IBD vs:");
            let current = self
                .all_alignments
                .iter()
                .find(|a| Some(a.id) == self.ibd_other)
                .map(|a| format!("#{} {}", a.id, a.reference_build))
                .unwrap_or_else(|| "(pick alignment)".into());
            egui::ComboBox::from_id_salt("ibd_other").selected_text(current).show_ui(ui, |ui| {
                for a in &self.all_alignments {
                    if a.id != alignment_id {
                        ui.selectable_value(&mut self.ibd_other, Some(a.id), format!("#{} {}", a.id, a.reference_build));
                    }
                }
            });
            let ready = self.ibd_other.is_some() && !self.running_ibd;
            if ui.add_enabled(ready, egui::Button::new("Compare")).clicked() {
                self.running_ibd = true;
                self.ibd_result = None;
                let _ = self.tx.send(Command::CompareIbd {
                    a: alignment_id,
                    b: self.ibd_other.unwrap(),
                    panel_id,
                    ploidy: self.ploidy(),
                });
            }
            if self.running_ibd {
                ui.spinner();
            }
        });

        if let Some(cmp) = &self.ibd_result {
            ui.label(format!(
                "{:?} — total {:.1} cM, {} segment(s), longest {:.1} cM",
                cmp.summary.relationship,
                cmp.summary.total_shared_cm,
                cmp.summary.segment_count,
                cmp.summary.longest_segment_cm,
            ));
            if !cmp.segments.is_empty() {
                egui::Grid::new("ibd_segments").striped(true).num_columns(4).show(ui, |ui| {
                    ui.strong("Chr");
                    ui.strong("Start");
                    ui.strong("End");
                    ui.strong("cM");
                    ui.end_row();
                    for s in &cmp.segments {
                        ui.label(&s.chromosome);
                        ui.label(s.start_position.to_string());
                        ui.label(s.end_position.to_string());
                        ui.label(format!("{:.1}", s.length_cm));
                        ui.end_row();
                    }
                });
            }
        }
    }

    fn denovo_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("De-novo SNP calls (haploid)");

        let has_paths = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some() && a.reference_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            ui.label("Contig:");
            ui.add(egui::TextEdit::singleline(&mut self.forms.denovo_contig).desired_width(80.0));
            let contig = self.forms.denovo_contig.trim().to_string();
            let ready = has_paths && !self.running_denovo && !contig.is_empty();
            if ui.add_enabled(ready, egui::Button::new("Run de-novo")).clicked() {
                self.running_denovo = true;
                self.denovo = None;
                self.status = format!("Calling {contig} on alignment #{alignment_id}…");
                let _ = self.tx.send(Command::RunDenovo { alignment_id, contig });
            }
            if self.running_denovo {
                ui.spinner();
            }
            if !has_paths {
                ui.label("(no BAM/reference path recorded)");
            }
        });

        match &self.denovo {
            None if !self.running_denovo => {
                ui.label("No calls yet — run for the contig above.");
            }
            None => {}
            Some(calls) if calls.is_empty() => {
                ui.label("0 SNP calls.");
            }
            Some(calls) => {
                ui.label(format!("{} SNP call(s)", calls.len()));
                egui::Grid::new("denovo_calls").striped(true).num_columns(4).show(ui, |ui| {
                    ui.strong("Position");
                    ui.strong("Change");
                    ui.strong("Depth");
                    ui.strong("AF");
                    ui.end_row();
                    for c in calls {
                        ui.label(c.position.to_string());
                        ui.label(format!("{}>{}", c.reference_allele, c.alternate_allele));
                        ui.label(c.depth.to_string());
                        ui.label(format!("{:.2}", c.allele_fraction));
                        ui.end_row();
                    }
                });
            }
        }

        let has_calls = self.denovo.as_ref().map(|c| !c.is_empty()).unwrap_or(false);
        let contig = self.forms.denovo_contig.trim().to_string();
        if has_calls && !contig.is_empty() {
            self.publish_row(ui, "Publish variants to PDS", Command::PublishVariants { alignment_id, contig });
        }
    }
}
