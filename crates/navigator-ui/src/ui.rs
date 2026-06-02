//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use eframe::egui;
use navigator_app::{Coverage, ProjectOverview};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, SequenceRun};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event, NewBiosample};

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
}

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    overview: Vec<ProjectOverview>,
    selected_project: Option<i64>,
    samples: Vec<Biosample>,
    selected_sample: Option<SampleGuid>,
    runs: Vec<SequenceRun>,
    selected_run: Option<i64>,
    alignments: Vec<Alignment>,
    selected_alignment: Option<i64>,
    coverage: Option<Coverage>,
    running: bool,
    forms: Forms,
    status: String,
}

fn opt(s: &str) -> Option<String> {
    let t = s.trim();
    (!t.is_empty()).then(|| t.to_string())
}

impl NavigatorApp {
    pub fn new(cc: &eframe::CreationContext<'_>, db_path: PathBuf) -> Self {
        let ctx = cc.egui_ctx.clone();
        let (tx, rx) = worker::spawn(db_path, move || ctx.request_repaint());
        let _ = tx.send(Command::LoadOverview);
        NavigatorApp {
            tx,
            rx,
            overview: Vec::new(),
            selected_project: None,
            samples: Vec::new(),
            selected_sample: None,
            runs: Vec::new(),
            selected_run: None,
            alignments: Vec::new(),
            selected_alignment: None,
            coverage: None,
            running: false,
            forms: Forms::default(),
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
                Event::SamplesChanged(project_id) => {
                    if self.selected_project == Some(project_id) {
                        let _ = self.tx.send(Command::LoadSamples(project_id));
                    }
                    let _ = self.tx.send(Command::LoadOverview); // sample counts changed
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
                Event::Alignments { sequence_run_id, alignments } => {
                    if self.selected_run == Some(sequence_run_id) {
                        self.alignments = alignments;
                    }
                }
                Event::AlignmentsChanged(run_id) => {
                    if self.selected_run == Some(run_id) {
                        let _ = self.tx.send(Command::LoadAlignments(run_id));
                    }
                }
                Event::Coverage { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.coverage = result;
                    }
                    self.running = false;
                }
                Event::Error(e) => {
                    self.status = format!("Error: {e}");
                    self.running = false;
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
        let _ = self.tx.send(Command::LoadRuns(guid));
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
        let _ = self.tx.send(Command::LoadCoverage(id));
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
        self.projects_panel(ctx);
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Status:");
                ui.label(&self.status);
            });
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            if self.selected_project.is_none() {
                ui.heading("DUNavigator");
                ui.label("Select or create a project.");
                return;
            }
            egui::ScrollArea::vertical().show(ui, |ui| {
                self.samples_section(ui);
                if self.selected_sample.is_some() {
                    self.runs_section(ui);
                }
                if self.selected_run.is_some() {
                    self.alignments_section(ui);
                }
                if let Some(id) = self.selected_alignment {
                    self.coverage_section(ui, id);
                }
            });
        });
    }
}

impl NavigatorApp {
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
        });
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

        ui.add_space(6.0);
        ui.collapsing("Add sample", |ui| {
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_donor).hint_text("donor identifier"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_accession).hint_text("accession (optional)"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_sex).hint_text("sex (optional)"));
            if ui
                .add_enabled(!self.forms.sample_donor.trim().is_empty(), egui::Button::new("Add sample"))
                .clicked()
            {
                let _ = self.tx.send(Command::AddBiosample(NewBiosample {
                    project_id: pid,
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

    fn runs_section(&mut self, ui: &mut egui::Ui) {
        let guid = self.selected_sample.unwrap();
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Sequence runs");
        if self.runs.is_empty() {
            ui.label("No runs yet.");
        }
        let mut pick = None;
        for r in &self.runs {
            let label = format!("#{}  {} · {}", r.id, r.platform_name, r.test_type);
            if ui.selectable_label(self.selected_run == Some(r.id), label).clicked() {
                pick = Some(r.id);
            }
        }
        if let Some(id) = pick {
            self.select_run(id);
        }

        ui.add_space(6.0);
        ui.collapsing("Add run", |ui| {
            ui.add(egui::TextEdit::singleline(&mut self.forms.run_platform).hint_text("platform (e.g. ILLUMINA)"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.run_test_type).hint_text("test type (e.g. WGS)"));
            let ready = !self.forms.run_platform.trim().is_empty() && !self.forms.run_test_type.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new("Add run")).clicked() {
                let _ = self.tx.send(Command::AddRun(NewSequenceRun {
                    biosample_guid: guid,
                    platform_name: self.forms.run_platform.trim().to_string(),
                    instrument_model: None,
                    test_type: self.forms.run_test_type.trim().to_string(),
                    library_layout: None,
                    total_reads: None,
                    pf_reads_aligned: None,
                    mean_read_length: None,
                    mean_insert_size: None,
                }));
                self.forms.run_platform.clear();
                self.forms.run_test_type.clear();
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
    }
}
