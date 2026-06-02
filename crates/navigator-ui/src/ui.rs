//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use eframe::egui;
use navigator_app::ProjectOverview;
use navigator_domain::workspace::{Biosample, NewProject};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event};

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    overview: Vec<ProjectOverview>,
    selected: Option<i64>,
    samples: Vec<Biosample>,
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
            selected: None,
            samples: Vec::new(),
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
                    self.select_project(p.id); // jump to the new project
                    let _ = self.tx.send(Command::LoadOverview); // and refresh the list
                }
                Event::Samples { project_id, samples } => {
                    if self.selected == Some(project_id) {
                        self.samples = samples;
                    }
                }
                Event::Error(e) => self.status = format!("Error: {e}"),
            }
        }
    }

    fn select_project(&mut self, id: i64) {
        self.selected = Some(id);
        self.samples.clear();
        let _ = self.tx.send(Command::LoadSamples(id));
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
                if ui.selectable_label(self.selected == Some(ov.project.id), label).clicked() {
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

        egui::CentralPanel::default().show(ctx, |ui| match self.selected {
            None => {
                ui.heading("DUNavigator");
                ui.label("Select a project to view its samples.");
            }
            Some(pid) => {
                let name = self
                    .overview
                    .iter()
                    .find(|o| o.project.id == pid)
                    .map(|o| o.project.name.as_str())
                    .unwrap_or("project");
                ui.heading(format!("Samples — {name}"));
                ui.separator();
                if self.samples.is_empty() {
                    ui.label("No samples yet.");
                } else {
                    egui::Grid::new("samples").striped(true).num_columns(3).show(ui, |ui| {
                        ui.strong("Donor");
                        ui.strong("Accession");
                        ui.strong("Sex");
                        ui.end_row();
                        for s in &self.samples {
                            ui.label(&s.donor_identifier);
                            ui.label(s.sample_accession.as_deref().unwrap_or("—"));
                            ui.label(s.sex.as_deref().unwrap_or("—"));
                            ui.end_row();
                        }
                    });
                }
            }
        });
    }
}
