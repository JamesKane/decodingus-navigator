//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use eframe::egui;
use navigator_app::{
    AncestryResult, AuditEntry, BuildNeed, CallState, CompatibilityLevel, Consensus, Coverage,
    DenovoCall, DnaType, HaploAssignment, HeteroplasmySite, IbdComparison, IdentityVerification,
    PanelGenotype, PrivateBucket, PrivateClass, ProjectOverview, ProjectSampleReport,
    ReconciledVariant, SourceType, VariantStatus, VerificationStatus,
};
use navigator_domain::ancestry::population_color;
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::chipprofile::{self, ChipProfile};
use navigator_domain::mtdna::MtdnaSequence;
use navigator_domain::strprofile::{self, StrProfile};
use navigator_domain::testtype;
use navigator_domain::variants::{VariantCall, VariantSet};
use navigator_domain::workspace::{Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, SequenceRun};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event, NewBiosample, PanelInfo, YMask};

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
    chip_provider: String,
    variant_source_type: String,
    variant_manual_label: String,
    variant_manual_text: String,
    /// Manual-override inputs for the Y / mtDNA consensus (corrected haplogroup + reason).
    override_y_haplogroup: String,
    override_y_reason: String,
    override_mt_haplogroup: String,
    override_mt_reason: String,
}

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    overview: Vec<ProjectOverview>,
    selected_project: Option<i64>,
    /// Per-sample coverage/haplogroup report rows for the selected project.
    project_report: Vec<ProjectSampleReport>,
    samples: Vec<Biosample>,
    /// Every biosample (the project-independent subjects list).
    all_biosamples: Vec<Biosample>,
    selected_sample: Option<SampleGuid>,
    runs: Vec<SequenceRun>,
    /// Donor-level haplogroup consensus for the selected subject (Y, mtDNA).
    consensus_y: Option<Consensus>,
    consensus_mt: Option<Consensus>,
    /// Reconciliation audit log for the selected subject (Y, mtDNA).
    audit_y: Vec<AuditEntry>,
    audit_mt: Vec<AuditEntry>,
    /// Last mtDNA heteroplasmy scan: (alignment id, sites).
    heteroplasmy: Option<(i64, Vec<HeteroplasmySite>)>,
    /// STR profiles for the selected subject.
    str_profiles: Vec<StrProfile>,
    /// SNP variant sets for the selected subject.
    variant_sets: Vec<VariantSet>,
    /// Cross-source variant concordance for the selected subject.
    variant_concordance: Vec<ReconciledVariant>,
    /// Chip/array profiles for the selected subject.
    chip_profiles: Vec<ChipProfile>,
    /// mtDNA sequences for the selected subject.
    mtdna_sequences: Vec<MtdnaSequence>,
    /// Chosen rCRS reference FASTA, reused across mtDNA variant derivations.
    rcrs_path: Option<PathBuf>,
    /// Last mtDNA haplogroup assignment: (sequence id, assignment).
    mtdna_haplogroup: Option<(i64, HaploAssignment)>,
    /// Last Y haplogroup assignment: (alignment id, assignment).
    y_haplogroup: Option<(i64, HaploAssignment)>,
    /// Last ancestry estimate: (alignment id, result). `None` result = computed, no estimate.
    ancestry: Option<(i64, Option<AncestryResult>)>,
    estimating_ancestry: bool,
    /// Live genotyping progress for the in-flight estimate: (alignment id, done, total) contigs.
    ancestry_progress: Option<(i64, usize, usize)>,
    /// Last private Y bucket: (alignment id, bucket).
    private_y: Option<(i64, PrivateBucket)>,
    finding_private_y: bool,
    /// Callable-region BED (external mask), reused across private-Y runs.
    y_mask_path: Option<PathBuf>,
    /// Use the sample's own callable-Y BED (self-referential) rather than an external mask.
    y_self_mask: bool,
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
    /// Identity-verification result for the current IBD pair.
    identity: Option<IdentityVerification>,
    /// Signed-in account DID, or `None`. Gates the "Publish" actions.
    account: Option<String>,
    /// Whether the last PDS write reached the server (offline indicator).
    online: bool,
    logging_in: bool,
    publishing: bool,
    /// A batch project-directory import is in flight (disables the button).
    importing: bool,
    /// The dir to retry importing once needed references are downloaded.
    pending_import_dir: Option<PathBuf>,
    /// Reference builds an import is waiting on (prompt the user to download).
    reference_needs: Vec<BuildNeed>,
    /// In-flight reference download: (build, received, total).
    reference_progress: Option<(String, u64, Option<u64>)>,
    /// A project-wide analyze pass is running (disables the report's analyze button).
    analyzing: bool,
    forms: Forms,
    status: String,
}

/// Sentinel option in the chip-provider dropdown that means "let the parser guess".
const AUTO_DETECT: &str = "(auto-detect)";

/// Format an optional mean/median depth (one decimal), "—" when not computed.
fn fmt_depth(o: Option<f64>) -> String {
    o.map(|v| format!("{v:.1}")).unwrap_or_else(|| "—".into())
}

/// Format an optional fraction (0–1) as a percentage, "—" when not computed.
fn fmt_pct(o: Option<f64>) -> String {
    o.map(|v| format!("{:.1}%", v * 100.0)).unwrap_or_else(|| "—".into())
}

/// Parse a `#RRGGBB` hex color, falling back to grey on a malformed string.
fn parse_hex_color(hex: &str) -> egui::Color32 {
    let h = hex.trim_start_matches('#');
    if h.len() == 6 {
        if let (Ok(r), Ok(g), Ok(b)) = (
            u8::from_str_radix(&h[0..2], 16),
            u8::from_str_radix(&h[2..4], 16),
            u8::from_str_radix(&h[4..6], 16),
        ) {
            return egui::Color32::from_rgb(r, g, b);
        }
    }
    egui::Color32::from_gray(128)
}

/// Render a haplogroup assignment: terminal + lineage + alternatives, then the child
/// branches with per-SNP evidence that explains why descent stopped.
fn show_assignment(ui: &mut egui::Ui, a: &HaploAssignment) {
    let Some(top) = a.ranked.first() else {
        ui.label("No match.");
        return;
    };
    ui.label(format!("Haplogroup: {}   ({}/{} mutations, score {:.3})", top.name, top.matched, top.expected, top.score));
    ui.label(format!("Lineage: {}", top.lineage.join(" › ")));
    let alts: Vec<String> = a.ranked.iter().skip(1).take(3).map(|r| format!("{} ({:.3})", r.name, r.score)).collect();
    if !alts.is_empty() {
        ui.label(format!("Alternatives: {}", alts.join(", ")));
    }
    for b in &a.branches {
        egui::CollapsingHeader::new(format!("child {} — {}/{} SNPs derived", b.name, b.derived, b.snps.len()))
            .id_salt(("branch", &b.name))
            .show(ui, |ui| {
                egui::Grid::new(("branch_snps", &b.name)).striped(true).num_columns(3).show(ui, |ui| {
                    for s in &b.snps {
                        ui.label(&s.name);
                        ui.label(format!("{}{}>{}", s.position, s.ancestral, s.derived));
                        let (txt, col) = match s.state {
                            CallState::Derived => ("derived", egui::Color32::from_rgb(60, 160, 60)),
                            CallState::Ancestral => ("ancestral", egui::Color32::from_rgb(170, 120, 40)),
                            CallState::NoCall => ("no-call", egui::Color32::GRAY),
                        };
                        ui.colored_label(col, txt);
                        ui.end_row();
                    }
                });
            });
    }
}

/// A readable "change" string for a variant call, covering the indel forms the mtDNA
/// derivation stores (one allele empty).
fn variant_change(c: &VariantCall) -> String {
    if c.alternate.is_empty() {
        format!("{}del", c.reference) // deletion
    } else if c.reference.is_empty() {
        format!("ins{}", c.alternate) // insertion
    } else {
        format!("{}>{}", c.reference, c.alternate) // substitution
    }
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
            project_report: Vec::new(),
            samples: Vec::new(),
            all_biosamples: Vec::new(),
            selected_sample: None,
            runs: Vec::new(),
            consensus_y: None,
            consensus_mt: None,
            audit_y: Vec::new(),
            audit_mt: Vec::new(),
            heteroplasmy: None,
            str_profiles: Vec::new(),
            variant_sets: Vec::new(),
            variant_concordance: Vec::new(),
            chip_profiles: Vec::new(),
            mtdna_sequences: Vec::new(),
            rcrs_path: None,
            mtdna_haplogroup: None,
            y_haplogroup: None,
            ancestry: None,
            estimating_ancestry: false,
            ancestry_progress: None,
            private_y: None,
            finding_private_y: false,
            y_mask_path: None,
            y_self_mask: true,
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
            identity: None,
            account: None,
            online: true,
            logging_in: false,
            publishing: false,
            importing: false,
            pending_import_dir: None,
            reference_needs: Vec::new(),
            reference_progress: None,
            analyzing: false,
            forms: Forms {
                denovo_contig: "chrM".into(),
                ploidy: "2".into(),
                run_test_type: "WGS".into(),
                str_panel: "Y-37".into(),
                str_provider: "FTDNA".into(),
                str_source: "DIRECT_TEST".into(),
                chip_provider: AUTO_DETECT.into(),
                variant_source_type: "IMPORTED".into(),
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
                Event::ProjectImported(summary) => {
                    let mut msg = format!(
                        "Imported {}: {} sample(s), {} alignment(s)",
                        summary.project.name, summary.samples_total, summary.alignments_created
                    );
                    if summary.alignments_skipped > 0 {
                        msg.push_str(&format!(" ({} already present)", summary.alignments_skipped));
                    }
                    if !summary.missing_index.is_empty() {
                        msg.push_str(&format!("; {} sample(s) missing an index", summary.missing_index.len()));
                    }
                    self.status = msg;
                    self.importing = false;
                    self.pending_import_dir = None;
                    self.reference_needs.clear();
                    self.select_project(summary.project.id);
                    let _ = self.tx.send(Command::LoadOverview);
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                }
                Event::ReferenceNeeded { dir, builds } => {
                    self.importing = false;
                    self.pending_import_dir = Some(dir);
                    self.status = format!("{} reference build(s) need downloading", builds.len());
                    self.reference_needs = builds;
                }
                Event::ReferenceProgress { build, received, total } => {
                    self.reference_progress = Some((build, received, total));
                }
                Event::ReferenceReady { build, path } => {
                    self.status = format!("Reference {build} ready ({})", path.display());
                    self.reference_progress = None;
                    self.reference_needs.retain(|b| b.build != build);
                    // When every needed build is in, retry the import automatically.
                    if self.reference_needs.is_empty() {
                        if let Some(dir) = self.pending_import_dir.take() {
                            self.importing = true;
                            self.status = format!("Importing {}…", dir.display());
                            let _ = self.tx.send(Command::ImportProjectDir { dir, reference: None });
                        }
                    }
                }
                Event::Samples { project_id, samples } => {
                    if self.selected_project == Some(project_id) {
                        self.samples = samples;
                    }
                }
                Event::ProjectReport { project_id, rows } => {
                    if self.selected_project == Some(project_id) {
                        self.project_report = rows;
                    }
                }
                Event::ProjectAnalyzed { project_id, samples, coverage_done, y_done, errors } => {
                    self.analyzing = false;
                    self.status = format!(
                        "Analyzed {samples} sample(s): {coverage_done} coverage, {y_done} Y{}",
                        if errors > 0 { format!(", {errors} error(s)") } else { String::new() }
                    );
                    if self.selected_project == Some(project_id) {
                        let _ = self.tx.send(Command::LoadProjectReport(project_id));
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
                        let _ = self.tx.send(Command::LoadVariantConcordance(guid)); // sources changed
                    }
                    self.status = "Variants imported".into();
                }
                Event::VariantConcordance { biosample_guid, variants } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.variant_concordance = variants;
                    }
                }
                Event::ChipProfiles { biosample_guid, profiles } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.chip_profiles = profiles;
                    }
                }
                Event::ChipProfilesChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadChipProfiles(guid));
                    }
                    self.status = "Chip data imported".into();
                }
                Event::MtdnaSequences { biosample_guid, sequences } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.mtdna_sequences = sequences;
                    }
                }
                Event::MtdnaChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadMtdna(guid));
                    }
                    self.status = "mtDNA sequence imported".into();
                }
                Event::Haplogroup { mtdna_id, assignment } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("mtDNA haplogroup: {} (score {:.3})", top.name, top.score),
                        None => "No haplogroup match".into(),
                    };
                    self.mtdna_haplogroup = Some((mtdna_id, assignment));
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensus(guid)); // a call was recorded
                    }
                }
                Event::YHaplogroup { alignment_id, assignment } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("Y haplogroup: {} (score {:.3})", top.name, top.score),
                        None => "No Y haplogroup match".into(),
                    };
                    self.y_haplogroup = Some((alignment_id, assignment));
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensus(guid));
                    }
                    // A per-row "Assign Y" from the project report just recorded a call —
                    // refresh the report so its Y column fills in.
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::AncestryProgress { alignment_id, done, total } => {
                    self.ancestry_progress = Some((alignment_id, done, total));
                    self.status = format!("Genotyping ancestry panel: {done}/{total} contigs…");
                }
                Event::Ancestry { alignment_id, result } => {
                    self.estimating_ancestry = false;
                    self.ancestry_progress = None;
                    self.status = match &result {
                        Some(r) => match r.primary() {
                            Some(top) => format!(
                                "Ancestry: {} {:.1}% ({}/{} SNPs)",
                                top.population_name, top.percentage, r.snps_with_genotype, r.snps_analyzed
                            ),
                            None => "Ancestry: no estimate".into(),
                        },
                        None => "Ancestry: not yet computed".into(),
                    };
                    self.ancestry = Some((alignment_id, result));
                }
                Event::Consensus { biosample_guid, y, mt } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.consensus_y = y;
                        self.consensus_mt = mt;
                    }
                }
                Event::Audit { biosample_guid, dna_type, entries } => {
                    if self.selected_sample == Some(biosample_guid) {
                        match dna_type {
                            DnaType::Y => self.audit_y = entries,
                            DnaType::Mt => self.audit_mt = entries,
                        }
                    }
                }
                Event::Heteroplasmy { alignment_id, sites } => {
                    self.status = format!("mtDNA heteroplasmy: {} site(s)", sites.len());
                    self.heteroplasmy = Some((alignment_id, sites));
                }
                Event::ReconciliationChanged { biosample_guid, dna_type } => {
                    if self.selected_sample == Some(biosample_guid) {
                        let _ = self.tx.send(Command::LoadConsensus(biosample_guid));
                        let _ = self.tx.send(Command::LoadAudit { biosample_guid, dna_type });
                    }
                }
                Event::PrivateY { alignment_id, bucket } => {
                    self.status = format!("Private Y: {} novel, {} off-path", bucket.novel(), bucket.off_path());
                    self.private_y = Some((alignment_id, bucket));
                    self.finding_private_y = false;
                }
                Event::DataImported { biosample_guid, label } => {
                    self.status = format!("Imported {label}");
                    if self.selected_sample == Some(biosample_guid) {
                        // Reload every data section — detection picked one of them.
                        let _ = self.tx.send(Command::LoadStrProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadVariantSets(biosample_guid));
                        let _ = self.tx.send(Command::LoadChipProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadMtdna(biosample_guid));
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
                    let _ = self.tx.send(Command::LoadAllAlignments); // keep IBD picker current
                }
                Event::Coverage { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.coverage = result;
                    }
                    self.running = false;
                    // A recompute (possibly from the project report) may have filled a cell.
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
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
                Event::Identity(v) => {
                    self.status = format!("Identity: {:?} ({} sites)", v.status, v.sites_compared);
                    self.identity = Some(v);
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
                    self.finding_private_y = false;
                    self.estimating_ancestry = false;
                    self.ancestry_progress = None;
                    let _ = self.tx.send(Command::SyncStatus); // a failed publish may have gone offline
                }
            }
        }
    }

    fn select_project(&mut self, id: i64) {
        self.selected_project = Some(id);
        self.samples.clear();
        self.project_report.clear();
        self.clear_sample_selection();
        let _ = self.tx.send(Command::LoadSamples(id));
        let _ = self.tx.send(Command::LoadProjectReport(id));
    }

    fn select_sample(&mut self, guid: SampleGuid) {
        self.selected_sample = Some(guid);
        self.clear_run_selection();
        self.runs.clear();
        self.str_profiles.clear();
        self.variant_sets.clear();
        self.variant_concordance.clear();
        self.chip_profiles.clear();
        self.mtdna_sequences.clear();
        self.mtdna_haplogroup = None;
        self.consensus_y = None;
        self.consensus_mt = None;
        self.audit_y.clear();
        self.audit_mt.clear();
        self.heteroplasmy = None;
        let _ = self.tx.send(Command::LoadConsensus(guid));
        let _ = self.tx.send(Command::LoadAudit { biosample_guid: guid, dna_type: DnaType::Y });
        let _ = self.tx.send(Command::LoadAudit { biosample_guid: guid, dna_type: DnaType::Mt });
        let _ = self.tx.send(Command::LoadVariantConcordance(guid));
        let _ = self.tx.send(Command::LoadRuns(guid));
        let _ = self.tx.send(Command::LoadStrProfiles(guid));
        let _ = self.tx.send(Command::LoadVariantSets(guid));
        let _ = self.tx.send(Command::LoadChipProfiles(guid));
        let _ = self.tx.send(Command::LoadMtdna(guid));
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
        self.identity = None;
        self.y_haplogroup = None;
        self.private_y = None;
        self.ancestry = None;
        self.estimating_ancestry = false;
        self.ancestry_progress = None;
        let _ = self.tx.send(Command::LoadCoverage(id));
        let _ = self.tx.send(Command::LoadAncestry { alignment_id: id });
        let _ = self.tx.send(Command::LoadDenovo { alignment_id: id, contig: self.forms.denovo_contig.clone() });
        if let Some(panel_id) = self.selected_panel {
            let _ = self.tx.send(Command::LoadPanelGenotypes { alignment_id: id, panel_id, ploidy: self.ploidy() });
        }
    }

    fn select_panel(&mut self, panel_id: i64) {
        self.selected_panel = Some(panel_id);
        self.panel_genotypes = None;
        self.ibd_result = None;
        self.identity = None;
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
        self.handle_file_drops(ctx);
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
                    self.project_report_section(ui);
                }
                if let Some(guid) = self.selected_sample {
                    self.subject_header(ui);
                    self.consensus_section(ui);
                    self.str_section(ui, guid);
                    self.variants_section(ui, guid);
                    self.chip_section(ui, guid);
                    self.mtdna_section(ui, guid);
                    self.runs_section(ui);
                }
                if self.selected_run.is_some() {
                    self.alignments_section(ui);
                }
                if let Some(id) = self.selected_alignment {
                    self.coverage_section(ui, id);
                    self.denovo_section(ui, id);
                    self.y_haplogroup_section(ui, id);
                    self.ancestry_section(ui, id);
                    self.heteroplasmy_section(ui, id);
                    self.genotyping_section(ui, id);
                }
            });
        });
        self.paint_drop_hint(ctx);
    }
}

impl NavigatorApp {
    /// Route files dropped onto the window through the unified importer, attaching them to
    /// the selected subject (auto-detected). No-op when nothing was dropped.
    fn handle_file_drops(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if dropped.is_empty() {
            return;
        }
        let Some(guid) = self.selected_sample else {
            self.status = "Select a subject before dropping data files.".into();
            return;
        };
        let mut sent = 0;
        for f in dropped {
            if let Some(path) = f.path {
                let _ = self.tx.send(Command::AddData { biosample_guid: guid, path });
                sent += 1;
            }
        }
        if sent > 0 {
            self.status = format!("Importing {sent} dropped file(s)…");
        }
    }

    /// While files are being dragged over the window, dim the screen and show whether the
    /// drop will land on a subject.
    fn paint_drop_hint(&self, ctx: &egui::Context) {
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

            ui.add_space(8.0);
            ui.label("Batch import a project folder (per-sample subdirs of BAM/CRAM):");
            if ui
                .add_enabled(!self.importing, egui::Button::new("Batch import project…"))
                .clicked()
            {
                if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                    self.importing = true;
                    self.reference_needs.clear();
                    self.status = format!("Importing {}…", dir.display());
                    // No reference up front: the gateway resolves each build from the cache,
                    // and reports ReferenceNeeded if a download is required.
                    let _ = self.tx.send(Command::ImportProjectDir { dir, reference: None });
                }
            }
            self.reference_prompt(ui);

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
        ui.horizontal(|ui| {
            ui.heading(format!("Subject — {donor}"));
            if ui.button("➕ Add data…").clicked() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("data files", &["vcf", "csv", "tsv", "txt", "fa", "fasta", "fna", "fas"])
                    .pick_file()
                {
                    let _ = self.tx.send(Command::AddData { biosample_guid: guid, path });
                }
            }
        });
        ui.label("Auto-detects STR / SNP variants / chip array / mtDNA from the file.");
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
            let header = format!("{} — {} markers  ({provider})", p.panel_name, p.markers.len());
            egui::CollapsingHeader::new(header).id_salt(("str", p.id)).show(ui, |ui| {
                egui::Grid::new(("str_markers", p.id)).striped(true).num_columns(2).show(ui, |ui| {
                    ui.strong("Marker");
                    ui.strong("Value");
                    ui.end_row();
                    for m in &p.markers {
                        ui.label(&m.marker);
                        ui.label(&m.value);
                        ui.end_row();
                    }
                });
            });
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
        const MAX_ROWS: usize = 500;
        for s in &self.variant_sets {
            let header = format!("{} — {} call(s)", s.source_label, s.calls.len());
            egui::CollapsingHeader::new(header).id_salt(("vset", s.id)).show(ui, |ui| {
                egui::Grid::new(("vcalls", s.id)).striped(true).num_columns(4).show(ui, |ui| {
                    for h in ["Position", "Change", "rsID", "Genotype"] {
                        ui.strong(h);
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
                if s.calls.len() > MAX_ROWS {
                    ui.label(format!("…and {} more", s.calls.len() - MAX_ROWS));
                }
            });
        }

        // Cross-source concordance (shown once ≥2 sources exist).
        if self.variant_sets.len() >= 2 && !self.variant_concordance.is_empty() {
            let confirmed = self.variant_concordance.iter().filter(|v| v.status == VariantStatus::Confirmed).count();
            let conflict = self.variant_concordance.iter().filter(|v| v.status == VariantStatus::Conflict).count();
            let single = self.variant_concordance.iter().filter(|v| v.status == VariantStatus::SingleSource).count();
            egui::CollapsingHeader::new(format!("Cross-source: {confirmed} confirmed · {conflict} conflict · {single} single"))
                .id_salt(("vconc", guid.0))
                .show(ui, |ui| {
                    egui::Grid::new(("vconc_grid", guid.0)).striped(true).num_columns(4).show(ui, |ui| {
                        for h in ["Position", "Allele", "Status", "Sources"] {
                            ui.strong(h);
                        }
                        ui.end_row();
                        for v in self.variant_concordance.iter().take(MAX_ROWS) {
                            ui.label(format!("{} {}", v.contig, v.position));
                            ui.label(&v.allele);
                            let (txt, col) = match v.status {
                                VariantStatus::Confirmed => ("confirmed", egui::Color32::from_rgb(60, 160, 60)),
                                VariantStatus::Conflict => ("conflict", egui::Color32::from_rgb(200, 60, 60)),
                                VariantStatus::SingleSource => ("single", egui::Color32::from_rgb(170, 150, 40)),
                            };
                            ui.colored_label(col, txt);
                            ui.label(format!("{}/{}", v.support, v.total));
                            ui.end_row();
                        }
                    });
                });
        }

        ui.add_space(6.0);
        ui.collapsing("Add / import variants", |ui| {
            let labels: Vec<&str> = SourceType::ALL.iter().map(|t| t.as_str()).collect();
            combo(ui, "Source:", "variant_source", &mut self.forms.variant_source_type, &labels);
            let source_type = SourceType::from_code(&self.forms.variant_source_type);

            if ui.button("Import VCF/CSV/TSV…").clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("variants", &["vcf", "csv", "tsv", "txt"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportVariants { biosample_guid: guid, path, source_type });
                }
            }
            ui.label("VCF, or rows of contig,position,ref,alt[,rsid][,genotype]. SNP-only.");

            ui.separator();
            ui.label("Or paste calls (e.g. Sanger / YSEQ confirmation):");
            ui.add(egui::TextEdit::singleline(&mut self.forms.variant_manual_label).hint_text("source label (e.g. YSEQ panel)"));
            ui.add(
                egui::TextEdit::multiline(&mut self.forms.variant_manual_text)
                    .hint_text("contig,position,ref,alt per line")
                    .desired_rows(3),
            );
            let ready = !self.forms.variant_manual_text.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new("Add pasted calls")).clicked() {
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
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Chip / array data");
        if self.chip_profiles.is_empty() {
            ui.label("No chip data imported yet.");
        }
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

        ui.add_space(6.0);
        ui.collapsing("Import chip data", |ui| {
            ui.horizontal(|ui| {
                ui.label("Provider:");
                egui::ComboBox::from_id_salt("chip_provider")
                    .selected_text(self.forms.chip_provider.clone())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.forms.chip_provider, AUTO_DETECT.to_string(), AUTO_DETECT);
                        for p in chipprofile::KNOWN_PROVIDERS {
                            ui.selectable_value(&mut self.forms.chip_provider, p.to_string(), *p);
                        }
                    });
            });
            if ui.button("Choose CSV/TXT…").clicked() {
                if let Some(path) = rfd::FileDialog::new().add_filter("array data", &["csv", "txt", "tsv"]).pick_file() {
                    let provider = (self.forms.chip_provider != AUTO_DETECT).then(|| self.forms.chip_provider.clone());
                    let _ = self.tx.send(Command::ImportChipProfile { biosample_guid: guid, provider, path });
                }
            }
            ui.label("23andMe / AncestryDNA / MyHeritage raw-data export.");
        });
    }

    /// mtDNA FASTA sequences for the selected subject + an import form, and a
    /// derive-variants-vs-rCRS action per sequence.
    fn mtdna_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("mtDNA sequences");
        if self.mtdna_sequences.is_empty() {
            ui.label("No mtDNA sequences yet.");
        }

        // rCRS reference picker (reused for every derivation this session).
        ui.horizontal(|ui| {
            ui.label("rCRS reference:");
            let label = self
                .rcrs_path
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "not set".into());
            ui.label(label);
            if ui.button("Choose rCRS…").clicked() {
                if let Some(p) = rfd::FileDialog::new().add_filter("FASTA", &["fa", "fasta", "fna", "fas"]).pick_file() {
                    self.rcrs_path = Some(p);
                }
            }
        });

        let rcrs = self.rcrs_path.clone();
        for m in &self.mtdna_sequences {
            let name = m.source_file_name.as_deref().or(m.defline.as_deref()).unwrap_or("mtDNA");
            ui.horizontal(|ui| {
                ui.label(format!("{name} — {} bp, {} N", m.length(), m.n_count));
                if ui
                    .add_enabled(rcrs.is_some(), egui::Button::new("Derive variants vs rCRS"))
                    .clicked()
                {
                    if let Some(path) = rcrs.clone() {
                        let _ = self.tx.send(Command::DeriveMtdnaVariants { mtdna_id: m.id, rcrs_path: path });
                    }
                }
                if ui.button("Assign haplogroup").clicked() {
                    self.status = "Assigning haplogroup (fetching FTDNA tree)…".into();
                    let _ = self.tx.send(Command::AssignMtdnaHaplogroup { mtdna_id: m.id });
                }
            });
            // Show the haplogroup result for this sequence, if any.
            if let Some((id, assignment)) = &self.mtdna_haplogroup {
                if *id == m.id {
                    show_assignment(ui, assignment);
                }
            }
        }

        ui.add_space(6.0);
        ui.collapsing("Import mtDNA FASTA", |ui| {
            if ui.button("Choose FASTA…").clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("FASTA", &["fa", "fasta", "fna", "fas"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportMtdna { biosample_guid: guid, path });
                }
            }
            ui.label("Full mtDNA sequence (~16,569 bp) aligned to rCRS.");
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

    /// When an import is blocked on uncached reference builds, prompt to download them (with
    /// a progress bar); on completion the import auto-retries (see the `ReferenceReady` event).
    fn reference_prompt(&mut self, ui: &mut egui::Ui) {
        if self.reference_needs.is_empty() && self.reference_progress.is_none() {
            return;
        }
        ui.add_space(6.0);
        egui::Frame::group(ui.style()).show(ui, |ui| {
            if !self.reference_needs.is_empty() {
                ui.label("Reference download required before import:");
                for b in &self.reference_needs {
                    ui.label(format!("  • {} (~{} MB)", b.build, b.est_bytes / 1_000_000));
                }
                if ui
                    .add_enabled(self.reference_progress.is_none(), egui::Button::new("Download & continue"))
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
    fn project_report_section(&mut self, ui: &mut egui::Ui) {
        if self.project_report.is_empty() {
            return;
        }
        ui.add_space(12.0);
        ui.separator();
        ui.horizontal(|ui| {
            ui.heading("Project report");
            let busy = self.analyzing || self.running;
            if ui.add_enabled(!busy, egui::Button::new("Analyze all (coverage + Y)")).clicked() {
                if let Some(pid) = self.selected_project {
                    self.analyzing = true;
                    self.status = "Analyzing project (coverage + Y per sample)…".into();
                    let _ = self.tx.send(Command::AnalyzeProject(pid));
                }
            }
            if ui.button("Export CSV…").clicked() {
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

        let running = self.running || self.analyzing;
        let mut recompute: Option<i64> = None;
        let mut assign_y: Option<i64> = None;
        egui::Grid::new("project_report_grid").striped(true).num_columns(10).show(ui, |ui| {
            for h in ["Sample", "Alns", "Mean cov", "Median", "≥10x", "≥20x", "Callable", "Y", "mtDNA", "actions"] {
                ui.strong(h);
            }
            ui.end_row();
            for r in &self.project_report {
                ui.label(&r.biosample.donor_identifier);
                ui.label(r.alignment_count.to_string());
                ui.label(fmt_depth(r.mean_coverage));
                ui.label(fmt_depth(r.median_coverage));
                ui.label(fmt_pct(r.pct_10x));
                ui.label(fmt_pct(r.pct_20x));
                ui.label(r.callable_bases.map(|v| v.to_string()).unwrap_or_else(|| "—".into()));
                ui.label(r.y_haplogroup.clone().unwrap_or_else(|| "—".into()));
                ui.label(r.mt_haplogroup.clone().unwrap_or_else(|| "—".into()));
                if let Some(aln) = r.primary_alignment_id {
                    ui.horizontal(|ui| {
                        if ui.add_enabled(!running, egui::Button::new("Cov")).clicked() {
                            recompute = Some(aln);
                        }
                        if ui.add_enabled(!running, egui::Button::new("Y")).clicked() {
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

    /// Donor-level haplogroup consensus across all recorded sources (runs, Sanger, …).
    fn consensus_section(&mut self, ui: &mut egui::Ui) {
        if self.consensus_y.is_none() && self.consensus_mt.is_none() {
            return; // nothing assigned yet for this subject
        }
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Haplogroup consensus");
        self.consensus_block(ui, "Y-DNA", DnaType::Y);
        self.consensus_block(ui, "mtDNA", DnaType::Mt);
    }

    /// One DNA type's consensus row plus its override controls, audit log, and publish
    /// button. The consensus/audit are cloned up front so the form fields can be borrowed
    /// mutably for the override inputs.
    fn consensus_block(&mut self, ui: &mut egui::Ui, label: &str, dna_type: DnaType) {
        let cons = match dna_type {
            DnaType::Y => self.consensus_y.clone(),
            DnaType::Mt => self.consensus_mt.clone(),
        };
        let Some(c) = cons else { return };
        let Some(guid) = self.selected_sample else { return };

        let (compat, col) = match c.compatibility {
            CompatibilityLevel::Compatible => ("compatible", egui::Color32::from_rgb(60, 160, 60)),
            CompatibilityLevel::MinorDivergence => ("minor divergence", egui::Color32::from_rgb(170, 150, 40)),
            CompatibilityLevel::MajorDivergence => ("major divergence", egui::Color32::from_rgb(200, 120, 40)),
            CompatibilityLevel::Incompatible => ("incompatible", egui::Color32::from_rgb(200, 60, 60)),
        };
        ui.horizontal(|ui| {
            ui.strong(format!("{label}: {}", c.haplogroup));
            ui.label(format!("({} source(s), conf {:.3})", c.run_count, c.confidence));
            ui.colored_label(col, compat);
            if c.overridden {
                ui.colored_label(egui::Color32::from_rgb(120, 120, 220), "manual override");
            }
        });
        for w in &c.warnings {
            ui.label(format!("  ⚠ {w}"));
        }

        // Manual override: set a curator-corrected terminal (e.g. Sanger-confirmed), or clear.
        let (hg_field, reason_field) = match dna_type {
            DnaType::Y => (&mut self.forms.override_y_haplogroup, &mut self.forms.override_y_reason),
            DnaType::Mt => (&mut self.forms.override_mt_haplogroup, &mut self.forms.override_mt_reason),
        };
        ui.horizontal(|ui| {
            ui.label("Override:");
            ui.add(egui::TextEdit::singleline(hg_field).hint_text("haplogroup").desired_width(140.0));
            ui.add(egui::TextEdit::singleline(reason_field).hint_text("reason").desired_width(180.0));
            let hg = hg_field.trim().to_string();
            let reason = reason_field.trim().to_string();
            if ui.add_enabled(!hg.is_empty(), egui::Button::new("Set")).clicked() {
                self.status = format!("Overriding {label} consensus → {hg}");
                let _ = self.tx.send(Command::SetHaploOverride {
                    biosample_guid: guid,
                    dna_type,
                    haplogroup: hg,
                    reason: (!reason.is_empty()).then_some(reason),
                });
            }
            if ui.add_enabled(c.overridden, egui::Button::new("Clear")).clicked() {
                self.status = format!("Clearing {label} override");
                let _ = self.tx.send(Command::ClearHaploOverride { biosample_guid: guid, dna_type });
            }
        });

        // mtDNA heteroplasmy observations from the last scan (folded into the published record).
        let het: Vec<HeteroplasmySite> = if dna_type == DnaType::Mt {
            self.heteroplasmy.as_ref().map(|(_, s)| s.clone()).unwrap_or_default()
        } else {
            Vec::new()
        };
        if dna_type == DnaType::Mt && !het.is_empty() {
            egui::CollapsingHeader::new(format!("heteroplasmy — {} site(s)", het.len()))
                .id_salt(("het", label))
                .show(ui, |ui| {
                    for h in &het {
                        ui.label(format!(
                            "  pos {}: {}/{} minor {:.1}% (depth {})",
                            h.position, h.major_base, h.minor_base, h.minor_fraction * 100.0, h.depth
                        ));
                    }
                });
        }

        // Audit log.
        let audit = match dna_type {
            DnaType::Y => &self.audit_y,
            DnaType::Mt => &self.audit_mt,
        };
        if !audit.is_empty() {
            egui::CollapsingHeader::new(format!("audit log — {} entr{}", audit.len(), if audit.len() == 1 { "y" } else { "ies" }))
                .id_salt(("audit", label))
                .show(ui, |ui| {
                    for e in audit {
                        ui.label(format!("  {} · {} — {}", e.timestamp, e.action, e.note));
                    }
                });
        }

        // Publish the donor-level reconciliation record (gated on sign-in).
        ui.horizontal(|ui| {
            let signed_in = self.account.is_some();
            if ui
                .add_enabled(signed_in && !self.publishing, egui::Button::new(format!("Publish {label} reconciliation")))
                .clicked()
            {
                self.publishing = true;
                self.status = format!("Publishing {label} reconciliation…");
                let _ = self.tx.send(Command::PublishReconciliation {
                    biosample_guid: guid,
                    dna_type,
                    heteroplasmy: het,
                    identity: self.identity.clone(),
                });
            }
            if !signed_in {
                ui.label("(sign in to publish)");
            }
        });
        ui.add_space(6.0);
    }

    /// mtDNA heteroplasmy scan for an alignment (chrM pileup → mixed positions). Results
    /// feed the mtDNA reconciliation record's heteroplasmy observations.
    fn heteroplasmy_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.add_space(12.0);
        ui.separator();
        ui.heading("mtDNA heteroplasmy");
        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new("Scan chrM heteroplasmy")).clicked() {
                self.status = "Scanning chrM pileup for heteroplasmy…".into();
                let _ = self.tx.send(Command::LoadHeteroplasmy { alignment_id });
            }
            if !has_bam {
                ui.label("(no BAM/CRAM path recorded)");
            }
        });

        if let Some((id, sites)) = &self.heteroplasmy {
            if *id == alignment_id {
                if sites.is_empty() {
                    ui.label("No heteroplasmic positions above the noise floor.");
                } else {
                    ui.label(format!("{} heteroplasmic position(s):", sites.len()));
                    for h in sites {
                        ui.label(format!(
                            "  pos {}: {} (major) / {} (minor) — minor {:.1}%, depth {}",
                            h.position, h.major_base, h.minor_base, h.minor_fraction * 100.0, h.depth
                        ));
                    }
                }
            }
        }
    }

    /// Ancestry estimate for an alignment: super-population proportions from the AIMs panel.
    fn ancestry_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Ancestry");

        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui
                .add_enabled(has_bam && !self.estimating_ancestry, egui::Button::new("Estimate ancestry"))
                .clicked()
            {
                self.estimating_ancestry = true;
                self.status = "Estimating ancestry (genotyping AIMs panel)…".into();
                let _ = self.tx.send(Command::EstimateAncestry { alignment_id });
            }
            if self.estimating_ancestry {
                ui.spinner();
            }
            if !has_bam {
                ui.label("(no BAM/CRAM path recorded)");
            }
        });

        // Live genotyping progress (the slow per-contig pass over the BAM).
        if self.estimating_ancestry {
            match self.ancestry_progress {
                Some((id, done, total)) if id == alignment_id && total > 0 => {
                    ui.add(
                        egui::ProgressBar::new(done as f32 / total as f32)
                            .text(format!("genotyping contig {done}/{total}")),
                    );
                }
                _ => {
                    ui.add(egui::ProgressBar::new(0.0).text("preparing…"));
                }
            }
        }

        match &self.ancestry {
            Some((id, Some(result))) if *id == alignment_id => {
                ui.label(format!(
                    "{}/{} panel SNPs genotyped · confidence {:.0}%",
                    result.snps_with_genotype,
                    result.snps_analyzed,
                    result.confidence_level * 100.0
                ));
                ui.add_space(4.0);
                egui::Grid::new("ancestry_components").striped(true).num_columns(3).show(ui, |ui| {
                    ui.strong("Super-population");
                    ui.strong("Share");
                    ui.strong("");
                    ui.end_row();
                    for c in &result.super_population_summary {
                        if c.percentage < 0.5 {
                            continue; // hide trace amounts
                        }
                        ui.label(&c.super_population);
                        ui.label(format!("{:.1}%", c.percentage));
                        // A simple proportion bar tinted by the population's catalog color.
                        let color = parse_hex_color(&population_color(
                            c.populations.first().map(String::as_str).unwrap_or(""),
                        ));
                        let (rect, _) =
                            ui.allocate_exact_size(egui::vec2(160.0, 12.0), egui::Sense::hover());
                        let painter = ui.painter();
                        painter.rect_filled(rect, 2.0, egui::Color32::from_gray(40));
                        let mut filled = rect;
                        filled.set_width(rect.width() * (c.percentage as f32 / 100.0).clamp(0.0, 1.0));
                        painter.rect_filled(filled, 2.0, color);
                        ui.end_row();
                    }
                });
            }
            Some((id, None)) if *id == alignment_id && !self.estimating_ancestry => {
                ui.label("No ancestry estimate yet.");
            }
            _ if !self.estimating_ancestry => {
                ui.label("No ancestry estimate yet.");
            }
            _ => {}
        }
    }

    /// Y-haplogroup assignment for an alignment (calls chrY tree positions; FTDNA tree).
    fn y_haplogroup_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        ui.add_space(12.0);
        ui.separator();
        ui.heading("Y haplogroup");

        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new("Assign Y haplogroup")).clicked() {
                self.status = "Assigning Y haplogroup (fetching FTDNA tree)…".into();
                let _ = self.tx.send(Command::AssignYHaplogroup { alignment_id });
            }
            if !has_bam {
                ui.label("(no BAM/CRAM path recorded)");
            }
        });

        if let Some((id, assignment)) = &self.y_haplogroup {
            if *id == alignment_id {
                show_assignment(ui, assignment);
            }
        }

        // Private bucket: de-novo chrY calls off the assigned backbone (branch candidates).
        let has_ref = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some() && a.reference_path.is_some())
            .unwrap_or(false);
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.label("Callable mask:");
            ui.checkbox(&mut self.y_self_mask, "self-referential (this sample)");
        });
        if !self.y_self_mask {
            ui.horizontal(|ui| {
                let label = self
                    .y_mask_path
                    .as_ref()
                    .and_then(|p| p.file_name())
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_else(|| "none (noisy)".into());
                ui.label(format!("External BED: {label}"));
                if ui.button("Choose BED…").clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("BED", &["bed"]).pick_file() {
                        self.y_mask_path = Some(p);
                    }
                }
            });
        }
        ui.horizontal(|ui| {
            if ui.add_enabled(has_ref && !self.finding_private_y, egui::Button::new("Find private Y variants")).clicked() {
                self.finding_private_y = true;
                self.status = "Finding private Y variants (de-novo chrY)…".into();
                let mask = if self.y_self_mask {
                    YMask::SelfReferential
                } else {
                    self.y_mask_path.clone().map(YMask::Bed).unwrap_or(YMask::None)
                };
                let _ = self.tx.send(Command::FindPrivateY { alignment_id, mask });
            }
            if self.finding_private_y {
                ui.spinner();
            }
            if !has_ref {
                ui.label("(needs BAM + reference path)");
            }
        });
        if let Some((id, bucket)) = &self.private_y {
            if *id == alignment_id {
                ui.label(format!("{} novel + {} off-path, below {}", bucket.novel(), bucket.off_path(), bucket.terminal));
                egui::CollapsingHeader::new("Private variants").id_salt(("privy", alignment_id)).show(ui, |ui| {
                    egui::Grid::new(("privy_grid", alignment_id)).striped(true).num_columns(4).show(ui, |ui| {
                        for h in ["Position", "Change", "Depth", "Class"] {
                            ui.strong(h);
                        }
                        ui.end_row();
                        for v in bucket.variants.iter().take(500) {
                            ui.label(v.position.to_string());
                            ui.label(format!("{}>{}", v.reference, v.alternate));
                            ui.label(v.depth.to_string());
                            match &v.class {
                                PrivateClass::Novel => ui.colored_label(egui::Color32::from_rgb(60, 160, 60), "novel"),
                                PrivateClass::OffPathKnown(name) => ui.label(format!("off-path: {name}")),
                            };
                            ui.end_row();
                        }
                    });
                    if bucket.variants.len() > 500 {
                        ui.label(format!("…and {} more", bucket.variants.len() - 500));
                    }
                });
            }
        }
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
        self.identity = None;
                let _ = self.tx.send(Command::CompareIbd {
                    a: alignment_id,
                    b: self.ibd_other.unwrap(),
                    panel_id,
                    ploidy: self.ploidy(),
                });
            }
            if ui.add_enabled(self.ibd_other.is_some(), egui::Button::new("Verify same donor")).clicked() {
                self.identity = None;
                let _ = self.tx.send(Command::VerifyIdentity {
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

        if let Some(v) = &self.identity {
            let (txt, col) = match v.status {
                VerificationStatus::VerifiedSame => ("same individual", egui::Color32::from_rgb(60, 160, 60)),
                VerificationStatus::LikelySame => ("likely same", egui::Color32::from_rgb(120, 160, 60)),
                VerificationStatus::Uncertain => ("uncertain", egui::Color32::from_rgb(170, 150, 40)),
                VerificationStatus::LikelyDifferent => ("likely different", egui::Color32::from_rgb(200, 120, 40)),
                VerificationStatus::VerifiedDifferent => ("different individuals", egui::Color32::from_rgb(200, 60, 60)),
            };
            ui.horizontal(|ui| {
                ui.label("Identity:");
                ui.colored_label(col, txt);
                if let Some(c) = v.snp_concordance {
                    ui.label(format!("SNP concordance {:.3} over {} sites", c, v.sites_compared));
                }
                if v.y_str_markers > 0 {
                    ui.label(format!("· Y-STR {}/{} differ", v.y_str_distance.unwrap_or(0), v.y_str_markers));
                }
            });
        }

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
