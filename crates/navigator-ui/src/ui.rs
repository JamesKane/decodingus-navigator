//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use eframe::egui;
use navigator_app::{Coverage, ProjectOverview};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Alignment, Biosample, NewProject};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event};

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    overview: Vec<ProjectOverview>,
    selected_project: Option<i64>,
    samples: Vec<Biosample>,
    selected_sample: Option<SampleGuid>,
    alignments: Vec<Alignment>,
    selected_alignment: Option<i64>,
    coverage: Option<Coverage>,
    running: bool,
    new_project_name: String,
    new_project_admin: String,
    status: String,
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
            alignments: Vec::new(),
            selected_alignment: None,
            coverage: None,
            running: false,
            new_project_name: String::new(),
            new_project_admin: String::new(),
            status: "Loading…".into(),
        }
    }

    /// Apply any events that arrived since the last frame to the view-state.
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
                Event::Alignments { biosample_guid, alignments } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.alignments = alignments;
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
        self.alignments.clear();
        self.selected_alignment = None;
        self.coverage = None;
        let _ = self.tx.send(Command::LoadAlignments(guid));
    }

    fn clear_sample_selection(&mut self) {
        self.selected_sample = None;
        self.alignments.clear();
        self.selected_alignment = None;
        self.coverage = None;
    }

    fn select_alignment(&mut self, id: i64) {
        self.selected_alignment = Some(id);
        self.coverage = None;
        let _ = self.tx.send(Command::LoadCoverage(id));
    }
}

impl eframe::App for NavigatorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.drain_events();

        egui::SidePanel::left("projects").min_width(220.0).show(ctx, |ui| {
            ui.heading("Projects");
            ui.separator();
            let mut to_select = None;
            for ov in &self.overview {
                let label = format!("{}  ({})", ov.project.name, ov.sample_count);
                if ui.selectable_label(self.selected_project == Some(ov.project.id), label).clicked() {
                    to_select = Some(ov.project.id);
                }
            }
            if let Some(id) = to_select {
                self.select_project(id);
            }

            ui.add_space(12.0);
            ui.separator();
            ui.label("New project");
            ui.add(egui::TextEdit::singleline(&mut self.new_project_name).hint_text("name"));
            ui.add(egui::TextEdit::singleline(&mut self.new_project_admin).hint_text("administrator"));
            let named = !self.new_project_name.trim().is_empty();
            if ui.add_enabled(named, egui::Button::new("Create")).clicked() {
                let admin = self.new_project_admin.trim();
                let _ = self.tx.send(Command::CreateProject(NewProject {
                    name: self.new_project_name.trim().to_string(),
                    description: None,
                    administrator: if admin.is_empty() { "unknown".into() } else { admin.to_string() },
                }));
                self.new_project_name.clear();
                self.new_project_admin.clear();
            }
        });

        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Status:");
                ui.label(&self.status);
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let Some(pid) = self.selected_project else {
                ui.heading("DUNavigator");
                ui.label("Select a project to view its samples.");
                return;
            };
            let project_name = self
                .overview
                .iter()
                .find(|o| o.project.id == pid)
                .map(|o| o.project.name.clone())
                .unwrap_or_else(|| "project".into());

            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading(format!("Samples — {project_name}"));
                ui.separator();
                if self.samples.is_empty() {
                    ui.label("No samples yet.");
                }
                let mut pick_sample = None;
                for s in &self.samples {
                    let label = format!(
                        "{}  ({}, {})",
                        s.donor_identifier,
                        s.sample_accession.as_deref().unwrap_or("—"),
                        s.sex.as_deref().unwrap_or("—"),
                    );
                    if ui.selectable_label(self.selected_sample == Some(s.guid), label).clicked() {
                        pick_sample = Some(s.guid);
                    }
                }
                if let Some(guid) = pick_sample {
                    self.select_sample(guid);
                }

                if self.selected_sample.is_some() {
                    self.alignment_section(ui);
                }
            });
        });
    }
}

impl NavigatorApp {
    fn alignment_section(&mut self, ui: &mut egui::Ui) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Alignments");
        if self.alignments.is_empty() {
            ui.label("No alignments for this sample.");
            return;
        }
        let mut pick = None;
        for a in &self.alignments {
            let label = format!("#{}  {} · {}", a.id, a.reference_build, a.aligner);
            if ui.selectable_label(self.selected_alignment == Some(a.id), label).clicked() {
                pick = Some(a.id);
            }
        }
        if let Some(id) = pick {
            self.select_alignment(id);
        }

        if let Some(id) = self.selected_alignment {
            self.coverage_section(ui, id);
        }
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
            let can_run = has_paths && !self.running;
            if ui.add_enabled(can_run, egui::Button::new("Run coverage")).clicked() {
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
