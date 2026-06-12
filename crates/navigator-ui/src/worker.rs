//! Sync↔async bridge. egui runs the immediate-mode loop on the main thread; the
//! `App` (tokio + sqlx) runs on a dedicated worker thread with its own runtime. The UI
//! sends [`Command`]s and drains [`Event`]s each frame — no DB calls or domain
//! decisions on the UI thread (plan §5).
//!
//! Each command is handled on its own task so a long analysis run never blocks quick
//! queries. The command→event mapping ([`handle`]) is pure and unit-tested; [`spawn`]
//! is the thread/runtime/channel glue.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender};
use std::sync::{Arc, Mutex};

use navigator_app::{
    AlignmentProbe, AncestryResult, AncestrySegment, App, AppError, AuditEntry, BuildNeed, Consensus, Coverage,
    DenovoCall, DnaType, HaploAssignment, HeteroplasmySite, IbdComparison, IbdDetectorConfig,
    IdentityVerification, PanelGenotype, PrivateBucket, ProjectImportSummary, ProjectOverview,
    ProjectSampleReport, ReadMetrics, ReconciledVariant, SexInferenceResult, SourceType,
    SvAnalysisResult,
};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::chipprofile::ChipProfile;
use navigator_domain::mtdna::MtdnaSequence;
use navigator_domain::strprofile::StrProfile;
use navigator_domain::variants::VariantSet;
use navigator_domain::workspace::{
    Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, Panel, Project, SequenceRun,
};
use tokio::sync::mpsc::{unbounded_channel, UnboundedSender};

/// Callable-region mask choice for the private-Y bucket.
#[derive(Debug, Clone)]
pub enum YMask {
    /// Self-referential: the sample's own callable-Y BED (adapts to depth + read tech).
    SelfReferential,
    /// An external callable BED (e.g. the Poznik/1KG `b38_sites.bed`).
    Bed(PathBuf),
    /// No mask (noisy — every off-backbone de-novo call).
    None,
}

/// Fields for adding a biosample (the app assigns its `SampleGuid`). `project_id` is
/// optional — biosamples are first-class and need not belong to a project.
#[derive(Debug, Clone)]
pub struct NewBiosample {
    pub project_id: Option<i64>,
    pub donor_identifier: String,
    pub sample_accession: Option<String>,
    pub sex: Option<String>,
}

/// A request from the UI to the worker.
#[derive(Debug, Clone)]
pub enum Command {
    LoadOverview,
    CreateProject(NewProject),
    LoadSamples(i64),
    /// Load the per-sample coverage/haplogroup report for a project.
    LoadProjectReport(i64),
    /// Analyze every sample in a project: coverage + Y haplogroup (fills the report).
    AnalyzeProject(i64),
    /// Load every biosample (subjects list), regardless of project.
    LoadAllBiosamples,
    /// Load donor-level Y/mt terminal haplogroups for every subject (fills the list columns).
    LoadHaploSummary,
    AddBiosample(NewBiosample),
    /// Batch-import a NAS project directory (scan → Project/Biosample/Run/Alignment).
    /// `reference` is optional: `None` lets the gateway resolve each build from the cache
    /// (and report `ReferenceNeeded` if a download is required); `Some` pins a FASTA.
    ImportProjectDir { dir: PathBuf, reference: Option<PathBuf> },
    /// Resolve (download + decompress + index) a reference build, streaming progress.
    ResolveReference { build: String },
    LoadRuns(SampleGuid),
    AddRun(NewSequenceRun),
    /// Load the donor-level Y + mtDNA haplogroup consensus for a subject.
    LoadConsensus(SampleGuid),
    /// Reconcile the subject's variant sets across sources.
    LoadVariantConcordance(SampleGuid),
    LoadStrProfiles(SampleGuid),
    ImportStrProfile {
        biosample_guid: SampleGuid,
        panel_name: String,
        provider: Option<String>,
        source: Option<String>,
        path: PathBuf,
    },
    LoadVariantSets(SampleGuid),
    ImportVariants { biosample_guid: SampleGuid, path: PathBuf, source_type: SourceType },
    /// Manually-entered variant calls (e.g. Sanger/YSEQ confirmations): `contig,pos,ref,alt` rows.
    AddVariants { biosample_guid: SampleGuid, source_label: String, source_type: SourceType, text: String },
    LoadChipProfiles(SampleGuid),
    ImportChipProfile { biosample_guid: SampleGuid, provider: Option<String>, path: PathBuf },
    LoadMtdna(SampleGuid),
    ImportMtdna { biosample_guid: SampleGuid, path: PathBuf },
    /// Derive mtDNA variants for a stored sequence vs an rCRS reference FASTA.
    DeriveMtdnaVariants { mtdna_id: i64, rcrs_path: PathBuf },
    /// Assign an mtDNA haplogroup (fetch the FTDNA tree, rank by the sample's base calls).
    AssignMtdnaHaplogroup { mtdna_id: i64 },
    /// Assign a Y haplogroup from an alignment (call chrY tree positions, rank).
    AssignYHaplogroup { alignment_id: i64 },
    /// Assign a Y haplogroup from the subject's imported BISDNA / Y-SNP-panel calls (no
    /// alignment) — records a donor call.
    AssignYBisdna { biosample_guid: SampleGuid },
    /// Assign an mtDNA haplogroup directly from an alignment's chrM (records a donor call).
    AssignMtdnaHaplogroupFromAlignment { alignment_id: i64 },
    /// Estimate ancestry (super-population proportions) from an alignment via the AIMs panel.
    EstimateAncestry { alignment_id: i64 },
    /// Load the persisted ancestry estimate for an alignment, if any.
    LoadAncestry { alignment_id: i64 },
    /// Load the reference population centroids (PC1,PC2) for the PCA scatter backdrop.
    LoadPcaReference { alignment_id: i64 },
    /// Paint each chromosome with local ancestry (genotypes the BAM; streams progress).
    PaintAncestry { alignment_id: i64 },
    /// Find the private bucket: de-novo chrY calls off the assigned Y backbone, restricted
    /// by the chosen callable mask.
    FindPrivateY { alignment_id: i64, mask: YMask },
    /// Load a previously-computed (self-masked) private-Y bucket from cache.
    LoadPrivateY { alignment_id: i64 },
    /// Unified import: detect the file's type and route it to the right importer.
    AddData { biosample_guid: SampleGuid, path: PathBuf },
    LoadAlignments(i64),
    AddAlignment(NewAlignment),
    /// Resolve the subject's default analysis alignment (highest-coverage, else first).
    DefaultAlignment { biosample_guid: SampleGuid },
    /// Load the subject's donor-level ancestry (best estimate across all sources).
    LoadDonorAncestry { biosample_guid: SampleGuid },
    /// Load the subject's donor-level private-Y union across all sources.
    LoadDonorPrivateY { biosample_guid: SampleGuid },
    /// Probe a BAM/CRAM header for build/aligner/platform/test-type (to auto-fill the form).
    ProbeAlignment { path: PathBuf },
    LoadCoverage(i64),
    RunCoverage(i64),
    LoadSex(i64),
    RunSex(i64),
    LoadReadMetrics(i64),
    RunReadMetrics(i64),
    LoadSv(i64),
    RunSv(i64),
    LoadDenovo { alignment_id: i64, contig: String },
    RunDenovo { alignment_id: i64, contig: String },
    LoadPanels,
    ImportPanel { name: String, path: PathBuf },
    LoadAllAlignments,
    GenotypePanel { alignment_id: i64, panel_id: i64, ploidy: u8 },
    LoadPanelGenotypes { alignment_id: i64, panel_id: i64, ploidy: u8 },
    CompareIbd { a: i64, b: i64, panel_id: i64, ploidy: u8 },
    /// Verify two alignments are the same individual (genotype concordance + Y-STR).
    VerifyIdentity { a: i64, b: i64, panel_id: i64, ploidy: u8 },
    /// Report who's signed in (no side effects) — sent on startup.
    AuthStatus,
    /// Report the current online/offline state (no side effects).
    SyncStatus,
    /// Sign in to a PDS via OAuth (opens a browser); `handle` is a handle or DID.
    Login { handle: String },
    Logout,
    PublishCoverage(i64),
    PublishVariants { alignment_id: i64, contig: String },
    /// Manually override the consensus haplogroup for a subject + DNA type.
    SetHaploOverride { biosample_guid: SampleGuid, dna_type: DnaType, haplogroup: String, reason: Option<String> },
    /// Clear a manual override.
    ClearHaploOverride { biosample_guid: SampleGuid, dna_type: DnaType },
    /// Load the reconciliation audit log for a subject + DNA type.
    LoadAudit { biosample_guid: SampleGuid, dna_type: DnaType },
    /// Scan an alignment's chrM pileup for heteroplasmic positions.
    LoadHeteroplasmy { alignment_id: i64 },
    /// Publish the subject's haplogroup reconciliation record to the signed-in PDS.
    PublishReconciliation {
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: Vec<HeteroplasmySite>,
        identity: Option<IdentityVerification>,
    },
    /// Run the full per-alignment analysis pipeline (coverage → sex → metrics → SV → variant
    /// calling → Y haplogroup → ancestry), streaming `AnalysisProgress` per step. Each step's
    /// own result event is forwarded too, so the detail tabs fill in as it runs.
    RunFullAnalysis { alignment_id: i64 },
    /// Request cancellation of the in-flight full analysis (checked between steps).
    CancelAnalysis,
    /// Update a subject's editable fields. Empty optional values clear the column.
    UpdateBiosample {
        guid: SampleGuid,
        donor_identifier: String,
        sample_accession: Option<String>,
        description: Option<String>,
        center_name: Option<String>,
        sex: Option<String>,
    },
    /// Delete a subject. Refused by the app layer if it still has dependent data.
    DeleteBiosample(SampleGuid),
    /// Delete a sequence run (cascades to its alignments + artifacts). `biosample_guid` is the
    /// owner, so the UI can refresh that subject's run list.
    DeleteSequenceRun { id: i64, biosample_guid: SampleGuid },
    /// Delete a single alignment (cascades to its artifacts). `sequence_run_id` is the owner,
    /// so the UI can refresh that run's alignment list.
    DeleteAlignment { id: i64, sequence_run_id: i64 },
    /// Delete an imported STR profile (and its markers).
    DeleteStrProfile { id: i64, biosample_guid: SampleGuid },
    /// Delete an imported variant set (and its calls).
    DeleteVariantSet { id: i64, biosample_guid: SampleGuid },
    /// Delete an imported chip/array profile.
    DeleteChipProfile { id: i64, biosample_guid: SampleGuid },
    /// Delete an imported mtDNA sequence.
    DeleteMtdnaSequence { id: i64, biosample_guid: SampleGuid },
    /// Assign a subject to a project (`None` clears it). The app layer validates the project.
    AssignBiosampleProject { guid: SampleGuid, project_id: Option<i64> },
    /// Update a project's editable fields.
    UpdateProject { id: i64, name: String, description: Option<String>, administrator: String },
    /// Delete a project. Refused by the app layer while subjects still belong to it.
    DeleteProject(i64),
    /// Update a sequence run's descriptive fields (read metrics preserved). `biosample_guid` is
    /// the owner so the UI can refresh that subject's run list.
    UpdateSequenceRun {
        id: i64,
        biosample_guid: SampleGuid,
        platform_name: String,
        instrument_model: Option<String>,
        test_type: String,
        library_layout: Option<String>,
    },
    /// Update an alignment's descriptive fields. `sequence_run_id` is the owner so the UI can
    /// refresh that run's alignment list.
    UpdateAlignment {
        id: i64,
        sequence_run_id: i64,
        reference_build: String,
        aligner: String,
        variant_caller: Option<String>,
    },
}

/// A panel with its site count, for the panel list.
#[derive(Debug, Clone)]
pub struct PanelInfo {
    pub panel: Panel,
    pub site_count: i64,
}

/// A result/notification from the worker to the UI.
#[derive(Debug, Clone)]
pub enum Event {
    /// Nothing to report (e.g. a cache-load that missed) — the UI ignores it.
    Noop,
    Overview(Vec<ProjectOverview>),
    ProjectCreated(Project),
    /// A project was updated or deleted; reload the overview.
    ProjectsChanged,
    /// A batch project-directory import completed.
    ProjectImported(ProjectImportSummary),
    /// Import needs reference build(s) downloaded first; `dir` lets the UI retry the import
    /// after the user approves and the download finishes.
    ReferenceNeeded { dir: PathBuf, builds: Vec<BuildNeed> },
    /// Bytes received during a reference download (`total` from Content-Length, if known).
    ReferenceProgress { build: String, received: u64, total: Option<u64> },
    /// A reference build finished resolving (cached + indexed).
    ReferenceReady { build: String, path: PathBuf },
    /// Per-sample coverage/haplogroup report for a project.
    ProjectReport { project_id: i64, rows: Vec<ProjectSampleReport> },
    /// A project-wide analyze pass finished (coverage + Y per sample).
    ProjectAnalyzed {
        project_id: i64,
        samples: usize,
        coverage_done: usize,
        y_done: usize,
        sex_done: usize,
        metrics_done: usize,
        sv_done: usize,
        errors: usize,
    },
    Samples { project_id: i64, samples: Vec<Biosample> },
    /// All biosamples (the project-independent subjects list).
    AllBiosamples(Vec<Biosample>),
    /// Per-subject Y/mt terminal haplogroups for the subjects list (`guid → (Y, mt)`).
    HaploSummary(std::collections::HashMap<SampleGuid, (Option<String>, Option<String>)>),
    /// A biosample was added/changed; reload the subjects list (and any open project view).
    BiosamplesChanged,
    Runs { biosample_guid: SampleGuid, runs: Vec<SequenceRun> },
    RunsChanged(SampleGuid),
    /// Donor-level haplogroup consensus for a subject (Y, mtDNA).
    Consensus { biosample_guid: SampleGuid, y: Option<Consensus>, mt: Option<Consensus> },
    /// Variant concordance across the subject's sources.
    VariantConcordance { biosample_guid: SampleGuid, variants: Vec<ReconciledVariant> },
    StrProfiles { biosample_guid: SampleGuid, profiles: Vec<StrProfile> },
    StrProfilesChanged(SampleGuid),
    VariantSets { biosample_guid: SampleGuid, sets: Vec<VariantSet> },
    VariantSetsChanged(SampleGuid),
    ChipProfiles { biosample_guid: SampleGuid, profiles: Vec<ChipProfile> },
    ChipProfilesChanged(SampleGuid),
    MtdnaSequences { biosample_guid: SampleGuid, sequences: Vec<MtdnaSequence> },
    MtdnaChanged(SampleGuid),
    /// A unified import succeeded; `label` describes the detected type. The UI should
    /// reload the subject's data sections.
    DataImported { biosample_guid: SampleGuid, label: String },
    /// mtDNA haplogroup assignment for a sequence (ranked + terminal evidence).
    Haplogroup { mtdna_id: i64, assignment: HaploAssignment },
    /// Y haplogroup assignment for an alignment (ranked + terminal evidence).
    YHaplogroup { alignment_id: i64, assignment: HaploAssignment },
    /// Y haplogroup assignment from the subject's BISDNA / Y-SNP panel (records a donor call).
    YBisdnaHaplogroup { biosample_guid: SampleGuid, assignment: HaploAssignment },
    /// mtDNA haplogroup assignment from an alignment (records a donor call → reload consensus).
    MtHaplogroup { alignment_id: i64, assignment: HaploAssignment },
    /// Progress of an ancestry genotyping pass: `done`/`total` contigs scanned.
    AncestryProgress { alignment_id: i64, done: usize, total: usize },
    /// Ancestry estimate for an alignment (`None` = not yet computed, for `LoadAncestry`).
    Ancestry { alignment_id: i64, result: Option<AncestryResult> },
    /// Reference population centroids (code, PC1, PC2) for the PCA scatter; empty if no loadings.
    PcaReference { alignment_id: i64, points: Vec<(String, f64, f64)> },
    /// Local-ancestry segments per chromosome (the "DNA painting").
    AncestryPainting { alignment_id: i64, segments: Vec<AncestrySegment> },
    /// Private Y variants (off-backbone de-novo calls) for an alignment.
    PrivateY { alignment_id: i64, bucket: PrivateBucket },
    Alignments { sequence_run_id: i64, alignments: Vec<Alignment> },
    AlignmentsChanged(i64),
    /// The subject's default analysis alignment, to auto-select on the detail tabs.
    DefaultAlignment { run_id: i64, alignment_id: i64 },
    /// Donor-level ancestry (best across sources) + the source alignment it came from.
    DonorAncestry { alignment_id: i64, result: AncestryResult },
    /// Donor-level private-Y union across the subject's sources.
    DonorPrivateY { bucket: PrivateBucket },
    /// Header-probe result for the add-alignment form (build/aligner/platform/test-type).
    AlignmentProbe(AlignmentProbe),
    Coverage { alignment_id: i64, result: Option<Coverage> },
    Sex { alignment_id: i64, result: Option<SexInferenceResult> },
    ReadMetrics { alignment_id: i64, result: Option<ReadMetrics> },
    Sv { alignment_id: i64, result: Option<SvAnalysisResult> },
    Denovo { alignment_id: i64, contig: String, result: Option<Vec<DenovoCall>> },
    /// Full-analysis pipeline progress: starting `step` of `total` (1-based), with a `label`
    /// + `detail` and the bar `fraction` (0..1).
    AnalysisProgress { step: usize, total: usize, label: String, detail: String, fraction: f32 },
    /// The full-analysis pipeline finished (or was cancelled).
    AnalysisDone { cancelled: bool },
    Panels(Vec<PanelInfo>),
    PanelImported,
    AllAlignments(Vec<Alignment>),
    PanelGenotypes { alignment_id: i64, panel_id: i64, ploidy: u8, genotypes: Vec<PanelGenotype> },
    Ibd(IbdComparison),
    /// Identity-verification result between two alignments.
    Identity(IdentityVerification),
    /// The reconciliation audit log for a subject + DNA type.
    Audit { biosample_guid: SampleGuid, dna_type: DnaType, entries: Vec<AuditEntry> },
    /// mtDNA heteroplasmy sites for an alignment.
    Heteroplasmy { alignment_id: i64, sites: Vec<HeteroplasmySite> },
    /// A reconciliation override/clear succeeded; the UI reloads consensus + audit.
    ReconciliationChanged { biosample_guid: SampleGuid, dna_type: DnaType },
    /// Current signed-in account (DID), or `None` when signed out.
    Authenticated(Option<String>),
    /// A record was published; `kind` is a human label, `uri` the `at://` URI.
    Published { kind: String, uri: String },
    /// Whether the last PDS write reached the server (offline indicator).
    SyncOnline(bool),
    Error(String),
}

/// Execute one command against the app, mapping success/failure to an [`Event`].
pub async fn handle(app: &App, cmd: Command) -> Event {
    match cmd {
        Command::LoadOverview => match app.project_overview().await {
            Ok(v) => Event::Overview(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::CreateProject(new) => match app.create_project(new).await {
            Ok(p) => Event::ProjectCreated(p),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportProjectDir { dir, reference } => {
            match app.import_project_dir(&dir, reference, "unknown".into(), true).await {
                Ok(summary) => Event::ProjectImported(summary),
                Err(AppError::ReferenceNeeded(builds)) => Event::ReferenceNeeded { dir, builds },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        // ResolveReference is handled in the spawn loop (it streams progress events); reaching
        // here would mean a routing bug.
        Command::ResolveReference { build } => Event::Error(format!("internal: unrouted ResolveReference {build}")),
        Command::LoadSamples(project_id) => match app.list_biosamples(project_id).await {
            Ok(samples) => Event::Samples { project_id, samples },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadProjectReport(project_id) => match app.project_report(project_id).await {
            Ok(rows) => Event::ProjectReport { project_id, rows },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AnalyzeProject(project_id) => match app.analyze_project(project_id).await {
            Ok(s) => Event::ProjectAnalyzed {
                project_id,
                samples: s.samples,
                coverage_done: s.coverage_done,
                y_done: s.y_done,
                sex_done: s.sex_done,
                metrics_done: s.metrics_done,
                sv_done: s.sv_done,
                errors: s.errors.len(),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadHaploSummary => match app.haplogroup_terminals().await {
            Ok(map) => Event::HaploSummary(map),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAllBiosamples => match app.list_all_biosamples().await {
            Ok(samples) => Event::AllBiosamples(samples),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddBiosample(b) => {
            match app
                .add_biosample(b.project_id, b.donor_identifier, b.sample_accession, b.sex)
                .await
            {
                Ok(_) => Event::BiosamplesChanged,
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::UpdateBiosample { guid, donor_identifier, sample_accession, description, center_name, sex } => {
            match app
                .update_biosample(guid, donor_identifier, sample_accession, description, center_name, sex)
                .await
            {
                Ok(_) => Event::BiosamplesChanged,
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::DeleteBiosample(guid) => match app.delete_biosample(guid).await {
            Ok(()) => Event::BiosamplesChanged,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteSequenceRun { id, biosample_guid } => match app.delete_sequence_run(id).await {
            Ok(()) => Event::RunsChanged(biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteAlignment { id, sequence_run_id } => match app.delete_alignment(id).await {
            Ok(()) => Event::AlignmentsChanged(sequence_run_id),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteStrProfile { id, biosample_guid } => match app.delete_str_profile(id).await {
            Ok(()) => Event::StrProfilesChanged(biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteVariantSet { id, biosample_guid } => match app.delete_variant_set(id).await {
            Ok(()) => Event::VariantSetsChanged(biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteChipProfile { id, biosample_guid } => match app.delete_chip_profile(id).await {
            Ok(()) => Event::ChipProfilesChanged(biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteMtdnaSequence { id, biosample_guid } => match app.delete_mtdna_sequence(id).await {
            Ok(()) => Event::MtdnaChanged(biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignBiosampleProject { guid, project_id } => {
            match app.add_biosample_to_project(guid, project_id).await {
                Ok(()) => Event::BiosamplesChanged,
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::UpdateProject { id, name, description, administrator } => {
            match app.update_project(id, name, description, administrator).await {
                Ok(_) => Event::ProjectsChanged,
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::DeleteProject(id) => match app.delete_project(id).await {
            Ok(()) => Event::ProjectsChanged,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::UpdateSequenceRun { id, biosample_guid, platform_name, instrument_model, test_type, library_layout } => {
            match app
                .update_sequence_run(id, platform_name, instrument_model, test_type, library_layout)
                .await
            {
                Ok(_) => Event::RunsChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::UpdateAlignment { id, sequence_run_id, reference_build, aligner, variant_caller } => {
            match app.update_alignment(id, reference_build, aligner, variant_caller).await {
                Ok(_) => Event::AlignmentsChanged(sequence_run_id),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadRuns(biosample_guid) => match app.list_sequence_runs(biosample_guid).await {
            Ok(runs) => Event::Runs { biosample_guid, runs },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddRun(new) => match app.record_sequence_run(new).await {
            Ok(run) => Event::RunsChanged(run.biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadConsensus(guid) => {
            let y = app.haplogroup_consensus(guid, DnaType::Y).await.unwrap_or(None);
            let mt = app.haplogroup_consensus(guid, DnaType::Mt).await.unwrap_or(None);
            Event::Consensus { biosample_guid: guid, y, mt }
        }
        Command::LoadVariantConcordance(guid) => match app.reconcile_variants(guid).await {
            Ok(variants) => Event::VariantConcordance { biosample_guid: guid, variants },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadStrProfiles(guid) => match app.list_str_profiles(guid).await {
            Ok(profiles) => Event::StrProfiles { biosample_guid: guid, profiles },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportStrProfile { biosample_guid, panel_name, provider, source, path } => {
            match app
                .import_str_profile_from_csv(biosample_guid, &panel_name, provider, source, &path)
                .await
            {
                Ok(_) => Event::StrProfilesChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadVariantSets(guid) => match app.list_variant_sets(guid).await {
            Ok(sets) => Event::VariantSets { biosample_guid: guid, sets },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportVariants { biosample_guid, path, source_type } => {
            match app.import_variants_from_file(biosample_guid, &path, source_type).await {
                Ok(_) => Event::VariantSetsChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::AddVariants { biosample_guid, source_label, source_type, text } => {
            match app.add_variants(biosample_guid, &source_label, source_type, &text).await {
                Ok(_) => Event::VariantSetsChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadChipProfiles(guid) => match app.list_chip_profiles(guid).await {
            Ok(profiles) => Event::ChipProfiles { biosample_guid: guid, profiles },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportChipProfile { biosample_guid, provider, path } => {
            match app.import_chip_profile_from_csv(biosample_guid, provider, None, &path).await {
                Ok(_) => Event::ChipProfilesChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadMtdna(guid) => match app.list_mtdna_sequences(guid).await {
            Ok(sequences) => Event::MtdnaSequences { biosample_guid: guid, sequences },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportMtdna { biosample_guid, path } => {
            match app.import_mtdna_from_fasta(biosample_guid, &path).await {
                Ok(_) => Event::MtdnaChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::DeriveMtdnaVariants { mtdna_id, rcrs_path } => {
            match app.derive_mtdna_variants(mtdna_id, &rcrs_path).await {
                Ok(set) => Event::VariantSetsChanged(set.biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::AssignMtdnaHaplogroup { mtdna_id } => match app.assign_mtdna_haplogroup(mtdna_id).await {
            Ok(assignment) => Event::Haplogroup { mtdna_id, assignment },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignYBisdna { biosample_guid } => match app.assign_y_bisdna(biosample_guid, None).await {
            Ok(assignment) => Event::YBisdnaHaplogroup { biosample_guid, assignment },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignYHaplogroup { alignment_id } => match app.assign_y_haplogroup(alignment_id).await {
            Ok(assignment) => Event::YHaplogroup { alignment_id, assignment },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignMtdnaHaplogroupFromAlignment { alignment_id } => {
            match app.assign_mtdna_haplogroup_from_alignment(alignment_id).await {
                Ok(assignment) => Event::MtHaplogroup { alignment_id, assignment },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        // EstimateAncestry is handled in the spawn loop (it streams AncestryProgress); reaching
        // here would mean a routing bug.
        Command::EstimateAncestry { alignment_id } => {
            Event::Error(format!("internal: unrouted EstimateAncestry {alignment_id}"))
        }
        Command::LoadAncestry { alignment_id } => match app.ancestry_for_alignment(alignment_id).await {
            Ok(result) => Event::Ancestry { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadPcaReference { alignment_id } => match app.ancestry_pca_reference(alignment_id).await {
            Ok(points) => Event::PcaReference { alignment_id, points },
            Err(e) => Event::Error(e.to_string()),
        },
        // PaintAncestry is handled in the spawn loop (it streams AncestryProgress); reaching here
        // would mean a routing bug.
        Command::PaintAncestry { alignment_id } => {
            Event::Error(format!("internal: unrouted PaintAncestry {alignment_id}"))
        }
        // RunFullAnalysis streams AnalysisProgress from the spawn loop; CancelAnalysis sets the
        // shared cancel flag there. Reaching here would mean a routing bug.
        Command::RunFullAnalysis { alignment_id } => {
            Event::Error(format!("internal: unrouted RunFullAnalysis {alignment_id}"))
        }
        Command::CancelAnalysis => Event::Error("internal: unrouted CancelAnalysis".into()),
        Command::FindPrivateY { alignment_id, mask } => {
            let result = match mask {
                YMask::SelfReferential => app.private_y_variants_self_masked(alignment_id).await,
                YMask::Bed(p) => app.private_y_variants(alignment_id, Some(&p)).await,
                YMask::None => app.private_y_variants(alignment_id, None).await,
            };
            match result {
                Ok(bucket) => Event::PrivateY { alignment_id, bucket },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadPrivateY { alignment_id } => match app.cached_private_y(alignment_id).await {
            Ok(Some(bucket)) => Event::PrivateY { alignment_id, bucket },
            Ok(None) => Event::Noop,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddData { biosample_guid, path } => match app.add_data(biosample_guid, &path).await {
            Ok(detected) => Event::DataImported { biosample_guid, label: detected.description().to_string() },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAlignments(sequence_run_id) => match app.list_alignments(sequence_run_id).await {
            Ok(alignments) => Event::Alignments { sequence_run_id, alignments },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddAlignment(new) => match app.record_alignment(new).await {
            Ok(a) => Event::AlignmentsChanged(a.sequence_run_id),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ProbeAlignment { path } => match app.probe_alignment(path).await {
            Ok(probe) => Event::AlignmentProbe(probe),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DefaultAlignment { biosample_guid } => match app.default_alignment_for_subject(biosample_guid).await {
            Ok(Some((run_id, alignment_id))) => Event::DefaultAlignment { run_id, alignment_id },
            Ok(None) => Event::Noop,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadDonorAncestry { biosample_guid } => match app.donor_ancestry(biosample_guid).await {
            Ok(Some((alignment_id, result))) => Event::DonorAncestry { alignment_id, result },
            Ok(None) => Event::Noop,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadDonorPrivateY { biosample_guid } => match app.donor_private_y(biosample_guid).await {
            Ok(Some(bucket)) => Event::DonorPrivateY { bucket },
            Ok(None) => Event::Noop,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadCoverage(alignment_id) => match app.cached_coverage(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunCoverage(alignment_id) => match app.run_coverage_for_alignment(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSex(alignment_id) => match app.cached_sex(alignment_id).await {
            Ok(result) => Event::Sex { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunSex(alignment_id) => match app.run_sex(alignment_id).await {
            Ok(result) => Event::Sex { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadReadMetrics(alignment_id) => match app.cached_read_metrics(alignment_id).await {
            Ok(result) => Event::ReadMetrics { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunReadMetrics(alignment_id) => match app.run_read_metrics(alignment_id).await {
            Ok(result) => Event::ReadMetrics { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSv(alignment_id) => match app.cached_sv(alignment_id).await {
            Ok(result) => Event::Sv { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunSv(alignment_id) => match app.run_sv(alignment_id).await {
            Ok(result) => Event::Sv { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadDenovo { alignment_id, contig } => match app.cached_denovo(alignment_id, &contig).await {
            Ok(result) => Event::Denovo { alignment_id, contig, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunDenovo { alignment_id, contig } => {
            match app.run_denovo_for_alignment(alignment_id, contig.clone()).await {
                Ok(result) => Event::Denovo { alignment_id, contig, result: Some(result) },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadPanels => match app.list_panels().await {
            Ok(panels) => {
                let mut infos = Vec::with_capacity(panels.len());
                for panel in panels {
                    let site_count = app.panel_site_count(panel.id).await.unwrap_or(0);
                    infos.push(PanelInfo { panel, site_count });
                }
                Event::Panels(infos)
            }
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportPanel { name, path } => match app.import_panel_from_vcf(&name, &path).await {
            Ok(_) => Event::PanelImported,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAllAlignments => match app.list_all_alignments().await {
            Ok(alns) => Event::AllAlignments(alns),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::GenotypePanel { alignment_id, panel_id, ploidy } => {
            match app.genotype_panel(alignment_id, panel_id, ploidy).await {
                Ok(genotypes) => Event::PanelGenotypes { alignment_id, panel_id, ploidy, genotypes },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadPanelGenotypes { alignment_id, panel_id, ploidy } => {
            match app.cached_panel_genotypes(alignment_id, panel_id, ploidy).await {
                Ok(genotypes) => Event::PanelGenotypes {
                    alignment_id,
                    panel_id,
                    ploidy,
                    genotypes: genotypes.unwrap_or_default(),
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::CompareIbd { a, b, panel_id, ploidy } => {
            match app.compare_ibd(a, b, panel_id, ploidy, IbdDetectorConfig::default()).await {
                Ok(cmp) => Event::Ibd(cmp),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::VerifyIdentity { a, b, panel_id, ploidy } => match app.verify_identity(a, b, panel_id, ploidy).await {
            Ok(v) => Event::Identity(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AuthStatus => Event::Authenticated(app.current_account()),
        Command::SyncStatus => Event::SyncOnline(app.is_online()),
        Command::Login { handle } => match app.login(&handle).await {
            Ok(did) => Event::Authenticated(Some(did)),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::Logout => match app.logout().await {
            Ok(()) => Event::Authenticated(None),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PublishCoverage(alignment_id) => match app.publish_coverage(alignment_id).await {
            Ok(r) => Event::Published { kind: "coverage summary".into(), uri: r.uri },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PublishVariants { alignment_id, contig } => {
            match app.publish_variants(alignment_id, &contig).await {
                Ok(r) => Event::Published { kind: format!("{contig} variants"), uri: r.uri },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::SetHaploOverride { biosample_guid, dna_type, haplogroup, reason } => {
            match app.set_manual_override(biosample_guid, dna_type, &haplogroup, reason.as_deref()).await {
                Ok(()) => Event::ReconciliationChanged { biosample_guid, dna_type },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::ClearHaploOverride { biosample_guid, dna_type } => {
            match app.clear_manual_override(biosample_guid, dna_type).await {
                Ok(()) => Event::ReconciliationChanged { biosample_guid, dna_type },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadAudit { biosample_guid, dna_type } => match app.reconciliation_audit(biosample_guid, dna_type).await {
            Ok(entries) => Event::Audit { biosample_guid, dna_type, entries },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadHeteroplasmy { alignment_id } => match app.mtdna_heteroplasmy(alignment_id).await {
            Ok(sites) => Event::Heteroplasmy { alignment_id, sites },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PublishReconciliation { biosample_guid, dna_type, heteroplasmy, identity } => {
            match app.publish_reconciliation(biosample_guid, dna_type, &heteroplasmy, identity.as_ref()).await {
                Ok(r) => Event::Published { kind: format!("{} reconciliation", dna_type.as_str()), uri: r.uri },
                Err(e) => Event::Error(e.to_string()),
            }
        }
    }
}

/// Resolve a reference build, emitting throttled `ReferenceProgress` events (and waking the
/// UI) as bytes arrive, then a final `ReferenceReady`/`Error`. Run from the spawn loop so it
/// can stream — `handle` returns only a single event.
async fn resolve_reference_streaming(app: &App, build: String, evt_tx: &Sender<Event>, wake: &(dyn Fn() + Send + Sync)) {
    // The progress closure must be Send (it runs in a task) — capture an owned Sender clone
    // and a label, not borrows. Throttle to ~every 25 MB so a multi-GB pull doesn't flood.
    let tx = evt_tx.clone();
    let label = build.clone();
    let mut last_sent = 0u64;
    let mut progress = move |received: u64, total: Option<u64>| {
        if received.saturating_sub(last_sent) >= 25_000_000 || total == Some(received) {
            last_sent = received;
            let _ = tx.send(Event::ReferenceProgress { build: label.clone(), received, total });
            wake();
        }
    };
    let event = match app.resolve_reference(&build, &mut progress).await {
        Ok(path) => Event::ReferenceReady { build, path },
        Err(e) => Event::Error(e.to_string()),
    };
    let _ = evt_tx.send(event);
    wake();
}

/// Estimate ancestry, emitting `AncestryProgress` per genotyped contig (and waking the UI),
/// then a final `Ancestry`/`Error`. Run from the spawn loop so it can stream.
async fn estimate_ancestry_streaming(
    app: &App,
    alignment_id: i64,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    // The progress closure runs on the blocking genotyping thread → must be Send + 'static:
    // capture owned clones of the sender and wake, not borrows.
    let tx = evt_tx.clone();
    let wake_cb = wake.clone();
    let progress = move |done: usize, total: usize| {
        let _ = tx.send(Event::AncestryProgress { alignment_id, done, total });
        wake_cb();
    };
    let event = match app.estimate_ancestry_with_progress(alignment_id, progress).await {
        Ok(result) => Event::Ancestry { alignment_id, result: Some(result) },
        Err(e) => Event::Error(e.to_string()),
    };
    let _ = evt_tx.send(event);
    wake();
}

/// Paint local ancestry, emitting `AncestryProgress` per genotyped contig, then a final
/// `AncestryPainting`/`Error`. Run from the spawn loop so it can stream.
async fn paint_ancestry_streaming(
    app: &App,
    alignment_id: i64,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    let tx = evt_tx.clone();
    let wake_cb = wake.clone();
    let progress = move |done: usize, total: usize| {
        let _ = tx.send(Event::AncestryProgress { alignment_id, done, total });
        wake_cb();
    };
    let event = match app.local_ancestry_with_progress(alignment_id, progress).await {
        Ok(segments) => Event::AncestryPainting { alignment_id, segments },
        Err(e) => Event::Error(e.to_string()),
    };
    let _ = evt_tx.send(event);
    wake();
}

/// Run the full per-alignment analysis pipeline, emitting `AnalysisProgress` before each step
/// and forwarding each step's own result event (so the detail tabs fill in live). `cancel` is
/// checked between steps. Per-step errors are forwarded but don't abort the pipeline (best-effort).
async fn run_full_analysis_streaming<W: Fn() + Send + Sync + 'static>(
    app: &App,
    alignment_id: i64,
    cancel: Arc<AtomicBool>,
    evt_tx: &Sender<Event>,
    wake: Arc<W>,
) {
    cancel.store(false, Ordering::Relaxed);
    let total = 6; // unified metrics + 4 command steps + ancestry

    // Step 1: unified quality metrics — coverage + callable, read-level QC, and sex inference in
    // ONE pass over the alignment (was three separate steps reading the file 2–3×). The slow
    // whole-genome read; stream per-contig sub-progress so the bar advances chromosome by
    // chromosome instead of sitting at 0% for minutes.
    if !cancel.load(Ordering::Relaxed) {
        let _ = evt_tx.send(Event::AnalysisProgress {
            step: 1,
            total,
            label: "Quality metrics".into(),
            detail: "scanning contigs…".into(),
            fraction: 0.0,
        });
        wake();
        // Reuse cached sub-results instead of re-scanning the whole genome (minutes) — only when
        // all three are present, since they're persisted together by the unified walker.
        let cached = match (
            app.cached_coverage(alignment_id).await,
            app.cached_read_metrics(alignment_id).await,
            app.cached_sex(alignment_id).await,
        ) {
            (Ok(Some(cov)), Ok(Some(rm)), Ok(Some(sex))) => Some((cov, rm, Some(sex))),
            _ => None,
        };
        let outcome = match cached {
            Some(triple) => Ok(triple),
            None => {
                // The parallel walker invokes progress from worker threads, so the callback must
                // be Fn + Sync; the event Sender is !Sync, so guard it with a Mutex.
                let evt = Arc::new(Mutex::new(evt_tx.clone()));
                let wk = wake.clone();
                app.run_unified_metrics_with_progress(alignment_id, move |done, tot| {
                    let within = if tot > 0 { done as f32 / tot as f32 } else { 0.0 };
                    if let Ok(tx) = evt.lock() {
                        let _ = tx.send(Event::AnalysisProgress {
                            step: 1,
                            total,
                            label: "Quality metrics".into(),
                            detail: format!("contig {done}/{tot}"),
                            fraction: within / total as f32,
                        });
                    }
                    wk();
                })
                .await
                .map(|r| (r.coverage, r.read_metrics, r.sex))
            }
        };
        // Emit the same per-result events the old separate steps did, so the UI updates identically.
        match outcome {
            Ok((cov, rm, sex)) => {
                let _ = evt_tx.send(Event::Coverage { alignment_id, result: Some(cov) });
                let _ = evt_tx.send(Event::ReadMetrics { alignment_id, result: Some(rm) });
                let _ = evt_tx.send(Event::Sex { alignment_id, result: sex });
            }
            Err(e) => {
                let _ = evt_tx.send(Event::Error(e.to_string()));
            }
        }
        wake();
    }

    // Steps 2–5: the cheaper command-driven steps (run via `handle`, which forwards their events).
    // Y variant discovery is the callable-masked "private Y" pass, NOT a raw whole-chrY de-novo
    // (which is enormous + mostly artifacts); chrM de-novo is fine (small, fully callable).
    let steps: [(&str, &str); 4] = [
        ("Structural variants", "CNV + discordant pairs (needs ≥10×)"),
        ("Variant calling", "chrM de-novo (haploid)"),
        ("Y haplogroup", "placing on the Y tree"),
        ("mtDNA haplogroup", "placing on the mt tree"),
    ];
    for (i, (label, detail)) in steps.iter().enumerate() {
        if cancel.load(Ordering::Relaxed) {
            break;
        }
        let step = i + 2; // steps 2..=5
        let _ = evt_tx.send(Event::AnalysisProgress {
            step,
            total,
            label: (*label).to_string(),
            detail: (*detail).to_string(),
            fraction: (step as f32 - 1.0) / total as f32,
        });
        wake();
        let cmd = match i {
            0 => Command::RunSv(alignment_id),
            1 => Command::RunDenovo { alignment_id, contig: "chrM".into() },
            2 => Command::AssignYHaplogroup { alignment_id },
            _ => Command::AssignMtdnaHaplogroupFromAlignment { alignment_id },
        };
        let ev = handle(app, cmd).await; // runs to completion; we may cancel before the next step
        let _ = evt_tx.send(ev);
        wake();
    }

    // Final step: ancestry (run directly — EstimateAncestry's command path streams separately).
    if !cancel.load(Ordering::Relaxed) {
        let _ = evt_tx.send(Event::AnalysisProgress {
            step: total,
            total,
            label: "Ancestry".to_string(),
            detail: "estimating population proportions".to_string(),
            fraction: (total as f32 - 1.0) / total as f32,
        });
        wake();
        // Reuse a cached estimate instead of re-genotyping the whole-genome panel (minutes).
        let result = match app.ancestry_for_alignment(alignment_id).await {
            Ok(Some(r)) => Ok(r),
            _ => app.estimate_ancestry(alignment_id).await,
        };
        let ev = match result {
            Ok(result) => Event::Ancestry { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        };
        let _ = evt_tx.send(ev);
        wake();
    }

    let _ = evt_tx.send(Event::AnalysisDone { cancelled: cancel.load(Ordering::Relaxed) });
    wake();
}

/// Spawn the worker thread: open the workspace at `db_path` inside the worker's runtime
/// (so the connection pool lives there), then serve commands. `wake` is called after
/// each event so the UI can `request_repaint`. Returns the command sender and event
/// receiver the UI holds.
pub fn spawn(
    db_path: PathBuf,
    wake: impl Fn() + Send + Sync + 'static,
) -> (UnboundedSender<Command>, Receiver<Event>) {
    let (cmd_tx, mut cmd_rx) = unbounded_channel::<Command>();
    let (evt_tx, evt_rx) = std::sync::mpsc::channel::<Event>();
    let wake = Arc::new(wake);

    std::thread::Builder::new()
        .name("navigator-worker".into())
        .spawn(move || {
            let rt = match tokio::runtime::Runtime::new() {
                Ok(rt) => rt,
                Err(e) => {
                    let _ = evt_tx.send(Event::Error(format!("runtime: {e}")));
                    wake();
                    return;
                }
            };
            rt.block_on(async move {
                let app = match App::open(&db_path).await {
                    Ok(app) => app,
                    Err(e) => {
                        let _ = evt_tx.send(Event::Error(format!("open workspace: {e}")));
                        wake();
                        return;
                    }
                };
                // Shared cancel flag for the full-analysis pipeline (one runs at a time).
                let cancel = Arc::new(AtomicBool::new(false));
                while let Some(cmd) = cmd_rx.recv().await {
                    let app = app.clone();
                    let evt_tx: Sender<Event> = evt_tx.clone();
                    let wake = wake.clone();
                    let cancel = cancel.clone();
                    tokio::spawn(async move {
                        match cmd {
                            // Streams ReferenceProgress events as bytes arrive, then a final event.
                            Command::ResolveReference { build } => {
                                resolve_reference_streaming(&app, build, &evt_tx, &*wake).await;
                            }
                            // Streams AncestryProgress per contig, then a final Ancestry/Error.
                            Command::EstimateAncestry { alignment_id } => {
                                estimate_ancestry_streaming(&app, alignment_id, &evt_tx, wake.clone()).await;
                            }
                            // Streams AncestryProgress per contig, then a final AncestryPainting.
                            Command::PaintAncestry { alignment_id } => {
                                paint_ancestry_streaming(&app, alignment_id, &evt_tx, wake.clone()).await;
                            }
                            // Streams AnalysisProgress per step (+ each step's result), then AnalysisDone.
                            Command::RunFullAnalysis { alignment_id } => {
                                run_full_analysis_streaming(&app, alignment_id, cancel, &evt_tx, wake.clone()).await;
                            }
                            // Signals the in-flight full analysis to stop between steps.
                            Command::CancelAnalysis => {
                                cancel.store(true, Ordering::Relaxed);
                            }
                            other => {
                                let event = handle(&app, other).await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                        }
                    });
                }
            });
        })
        .expect("spawn worker thread");

    (cmd_tx, evt_rx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_store::Store;

    async fn app() -> App {
        App::new(Store::open_in_memory().await.unwrap())
    }

    #[tokio::test]
    async fn handle_maps_commands_to_events() {
        let app = app().await;

        // create a project
        let created = handle(&app, Command::CreateProject(NewProject {
            name: "Trio".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await;
        let pid = match created {
            Event::ProjectCreated(p) => p.id,
            other => panic!("expected ProjectCreated, got {other:?}"),
        };

        // overview reflects it
        match handle(&app, Command::LoadOverview).await {
            Event::Overview(v) => {
                assert_eq!(v.len(), 1);
                assert_eq!(v[0].sample_count, 0);
            }
            other => panic!("expected Overview, got {other:?}"),
        }

        // samples for the project (empty)
        match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { project_id, samples } => {
                assert_eq!(project_id, pid);
                assert!(samples.is_empty());
            }
            other => panic!("expected Samples, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn handle_runs_coverage_for_a_stored_alignment() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};
        use std::path::PathBuf;

        let app = app().await;
        let b = app.add_biosample(None, "HG002", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                biosample_guid: b.guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await
            .unwrap();
        let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../navigator-analysis/tests/fixtures");
        let aln = app
            .record_alignment(NewAlignment {
                sequence_run_id: run.id,
                reference_build: "chrM".into(),
                aligner: "synthetic".into(),
                variant_caller: None,
                bam_path: Some(fixtures.join("coverage.bam").to_string_lossy().into_owned()),
                reference_path: Some(fixtures.join("ref.fa").to_string_lossy().into_owned()),
                content_sha256: None,
            })
            .await
            .unwrap();

        // alignments query (keyed by run)
        match handle(&app, Command::LoadAlignments(run.id)).await {
            Event::Alignments { alignments, .. } => assert_eq!(alignments, vec![aln.clone()]),
            other => panic!("expected Alignments, got {other:?}"),
        }

        // cold cache
        match handle(&app, Command::LoadCoverage(aln.id)).await {
            Event::Coverage { result, .. } => assert!(result.is_none()),
            other => panic!("expected Coverage(None), got {other:?}"),
        }

        // run + persist (uses the alignment's stored paths, via spawn_blocking)
        match handle(&app, Command::RunCoverage(aln.id)).await {
            Event::Coverage { alignment_id, result } => {
                assert_eq!(alignment_id, aln.id);
                assert_eq!(result.unwrap().genome_territory, 50);
            }
            other => panic!("expected Coverage(Some), got {other:?}"),
        }

        // now cached
        match handle(&app, Command::LoadCoverage(aln.id)).await {
            Event::Coverage { result, .. } => assert_eq!(result.unwrap().callable_bases, 10),
            other => panic!("expected cached Coverage, got {other:?}"),
        }

        // de-novo on the fixture contig (cold -> run -> cached), per-contig keyed
        match handle(&app, Command::LoadDenovo { alignment_id: aln.id, contig: "chrM".into() }).await {
            Event::Denovo { result, .. } => assert!(result.is_none()),
            other => panic!("expected Denovo(None), got {other:?}"),
        }
        match handle(&app, Command::RunDenovo { alignment_id: aln.id, contig: "chrM".into() }).await {
            Event::Denovo { contig, result, .. } => {
                assert_eq!(contig, "chrM");
                let calls = result.unwrap();
                assert_eq!(calls.iter().map(|c| c.position).collect::<Vec<_>>(), vec![2, 3, 4, 6, 7, 8, 10]);
            }
            other => panic!("expected Denovo(Some), got {other:?}"),
        }
        match handle(&app, Command::LoadDenovo { alignment_id: aln.id, contig: "chrM".into() }).await {
            Event::Denovo { result, .. } => assert_eq!(result.unwrap().len(), 7),
            other => panic!("expected cached Denovo, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn add_commands_create_and_signal_reload() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};

        let app = app().await;
        let pid = match handle(&app, Command::CreateProject(NewProject {
            name: "P".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };

        // add a sample (tagged to the project) -> BiosamplesChanged
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: Some(pid),
            donor_identifier: "HG002".into(),
            sample_accession: None,
            sex: Some("male".into()),
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        // it shows up both under the project and in the all-subjects list
        let guid = match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { samples, .. } => samples[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 1),
            other => panic!("got {other:?}"),
        }

        // a project-less subject is also allowed
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: None,
            donor_identifier: "NA12878".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 2),
            other => panic!("got {other:?}"),
        }

        // add a run -> RunsChanged(sample)
        match handle(&app, Command::AddRun(NewSequenceRun {
            biosample_guid: guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        }))
        .await
        {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        let run_id = match handle(&app, Command::LoadRuns(guid)).await {
            Event::Runs { runs, .. } => runs[0].id,
            other => panic!("got {other:?}"),
        };

        // add an alignment -> AlignmentsChanged(run)
        match handle(&app, Command::AddAlignment(NewAlignment {
            sequence_run_id: run_id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        }))
        .await
        {
            Event::AlignmentsChanged(r) => assert_eq!(r, run_id),
            other => panic!("got {other:?}"),
        }
    }

    #[tokio::test]
    async fn edit_and_delete_subject_commands() {
        use navigator_domain::workspace::NewSequenceRun;

        let app = app().await;

        // a project-less subject we can freely edit and delete
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: None,
            donor_identifier: "draft".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => all[0].guid,
            other => panic!("got {other:?}"),
        };

        // edit: set the identifier + a few optional fields
        match handle(&app, Command::UpdateBiosample {
            guid,
            donor_identifier: "HG002".into(),
            sample_accession: Some("SAMN123".into()),
            description: Some("trio son".into()),
            center_name: None,
            sex: Some("male".into()),
        })
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => {
                let b = &all[0];
                assert_eq!(b.donor_identifier, "HG002");
                assert_eq!(b.sample_accession.as_deref(), Some("SAMN123"));
                assert_eq!(b.description.as_deref(), Some("trio son"));
                assert_eq!(b.center_name, None);
                assert_eq!(b.sex.as_deref(), Some("male"));
            }
            other => panic!("got {other:?}"),
        }

        // adding dependent data makes delete refuse with a conflict
        match handle(&app, Command::AddRun(NewSequenceRun {
            biosample_guid: guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        }))
        .await
        {
            Event::RunsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::DeleteBiosample(guid)).await {
            Event::Error(msg) => assert!(msg.contains("sequencing run"), "unexpected message: {msg}"),
            other => panic!("expected conflict Error, got {other:?}"),
        }
        // still present
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 1),
            other => panic!("got {other:?}"),
        }

        // removing the run clears the conflict, so the subject can then be deleted (the
        // end-to-end 'remove data first' path)
        let run_id = match handle(&app, Command::LoadRuns(guid)).await {
            Event::Runs { runs, .. } => runs[0].id,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::DeleteSequenceRun { id: run_id, biosample_guid: guid }).await {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::DeleteBiosample(guid)).await {
            Event::BiosamplesChanged => {}
            other => panic!("expected clean delete, got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert!(all.iter().all(|b| b.guid != guid)),
            other => panic!("got {other:?}"),
        }

        // a subject with no dependents deletes cleanly
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: None,
            donor_identifier: "spare".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let spare = match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => all.iter().find(|b| b.donor_identifier == "spare").unwrap().guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::DeleteBiosample(spare)).await {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert!(all.iter().all(|b| b.guid != spare)),
            other => panic!("got {other:?}"),
        }
    }

    #[tokio::test]
    async fn assign_biosample_project_command() {
        let app = app().await;
        let pid = match handle(&app, Command::CreateProject(NewProject {
            name: "P".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: None,
            donor_identifier: "loose".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => all[0].guid,
            other => panic!("got {other:?}"),
        };

        // assign into the project
        match handle(&app, Command::AssignBiosampleProject { guid, project_id: Some(pid) }).await {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { samples, .. } => assert_eq!(samples.len(), 1),
            other => panic!("got {other:?}"),
        }

        // assigning to a non-existent project is refused
        match handle(&app, Command::AssignBiosampleProject { guid, project_id: Some(9999) }).await {
            Event::Error(_) => {}
            other => panic!("expected Error, got {other:?}"),
        }

        // clearing the project (None) removes it from the project list
        match handle(&app, Command::AssignBiosampleProject { guid, project_id: None }).await {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { samples, .. } => assert!(samples.is_empty()),
            other => panic!("got {other:?}"),
        }
    }

    #[tokio::test]
    async fn update_and_delete_project_commands() {
        let app = app().await;
        let pid = match handle(&app, Command::CreateProject(NewProject {
            name: "Old".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };

        // edit name/admin/description
        match handle(&app, Command::UpdateProject {
            id: pid,
            name: "Renamed".into(),
            description: Some("a study".into()),
            administrator: "curator".into(),
        })
        .await
        {
            Event::ProjectsChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadOverview).await {
            Event::Overview(v) => {
                let p = &v.iter().find(|o| o.project.id == pid).unwrap().project;
                assert_eq!(p.name, "Renamed");
                assert_eq!(p.description.as_deref(), Some("a study"));
                assert_eq!(p.administrator, "curator");
            }
            other => panic!("got {other:?}"),
        }

        // a project with members refuses delete
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: Some(pid),
            donor_identifier: "member".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::DeleteProject(pid)).await {
            Event::Error(msg) => assert!(msg.contains("subject"), "unexpected message: {msg}"),
            other => panic!("expected conflict Error, got {other:?}"),
        }

        // reassigning the member away clears the conflict, so delete succeeds
        let guid = match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { samples, .. } => samples[0].guid,
            other => panic!("got {other:?}"),
        };
        let _ = handle(&app, Command::AssignBiosampleProject { guid, project_id: None }).await;
        match handle(&app, Command::DeleteProject(pid)).await {
            Event::ProjectsChanged => {}
            other => panic!("expected clean delete, got {other:?}"),
        }
        match handle(&app, Command::LoadOverview).await {
            Event::Overview(v) => assert!(v.iter().all(|o| o.project.id != pid)),
            other => panic!("got {other:?}"),
        }
    }

    #[tokio::test]
    async fn update_run_and_alignment_commands() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};

        let app = app().await;
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: None,
            donor_identifier: "subj".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => all[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::AddRun(NewSequenceRun {
            biosample_guid: guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: None,
            total_reads: Some(1_000),
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        }))
        .await
        {
            Event::RunsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        let run = match handle(&app, Command::LoadRuns(guid)).await {
            Event::Runs { runs, .. } => runs[0].clone(),
            other => panic!("got {other:?}"),
        };

        // edit the run's descriptive fields; the read metric is preserved
        match handle(&app, Command::UpdateSequenceRun {
            id: run.id,
            biosample_guid: guid,
            platform_name: "MGI".into(),
            instrument_model: Some("DNBSEQ-T7".into()),
            test_type: "WGS".into(),
            library_layout: Some("PAIRED".into()),
        })
        .await
        {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadRuns(guid)).await {
            Event::Runs { runs, .. } => {
                let r = &runs[0];
                assert_eq!(r.platform_name, "MGI");
                assert_eq!(r.instrument_model.as_deref(), Some("DNBSEQ-T7"));
                assert_eq!(r.library_layout.as_deref(), Some("PAIRED"));
                assert_eq!(r.total_reads, Some(1_000)); // metric untouched
            }
            other => panic!("got {other:?}"),
        }

        match handle(&app, Command::AddAlignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "grch38".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        }))
        .await
        {
            Event::AlignmentsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        let aln_id = match handle(&app, Command::LoadAlignments(run.id)).await {
            Event::Alignments { alignments, .. } => alignments[0].id,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::UpdateAlignment {
            id: aln_id,
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "minimap2".into(),
            variant_caller: Some("deepvariant".into()),
        })
        .await
        {
            Event::AlignmentsChanged(r) => assert_eq!(r, run.id),
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAlignments(run.id)).await {
            Event::Alignments { alignments, .. } => {
                let a = &alignments[0];
                assert_eq!(a.reference_build, "chm13v2.0");
                assert_eq!(a.aligner, "minimap2");
                assert_eq!(a.variant_caller.as_deref(), Some("deepvariant"));
            }
            other => panic!("got {other:?}"),
        }
    }
}
