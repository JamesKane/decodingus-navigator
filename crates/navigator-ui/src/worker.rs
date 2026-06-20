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
    AlignmentProbe, AncestryResult, AncestrySegment, App, AppError, AuditEntry, BatchImportSummary, BuildNeed,
    Consensus, Coverage, DenovoCall, DnaType, ExchangeSessionInfo, FtdnaGenealogy, FtdnaImportOptions, FtdnaImportPlan,
    FtdnaImportSummary, FtdnaResolution, HaploAssignment, HeteroplasmySite, IbdComparison, IbdDetectorConfig,
    IbdSuggestion, IdentityVerification, IncomingRequest, PanelGenotype, PrivateBucket, ProjectImportSummary,
    ProjectOverview, ProjectSampleReport, ReadMetrics, RefBuildStatus, SexInferenceResult, SourceType,
    StoredIbdExchange, StrConcordanceRow, SvAnalysisResult, YMatch, YstrClustering,
};
use navigator_domain::chipprofile::ChipProfile;
use navigator_domain::du_domain::ids::SampleGuid;
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
    /// Load the ancestry/IBD asset presence + integrity status (the "data sources" line).
    LoadAssetStatus,
    CreateProject(NewProject),
    LoadSamples(i64),
    /// Load the per-sample coverage/haplogroup report for a project.
    LoadProjectReport(i64),
    /// Deep-analyze every sample in a project as a cancellable background job, streaming
    /// per-sample `DeepAnalyzeProgress` and yielding between samples so the UI stays responsive.
    /// Skips what the fast path already filled; cancelled via [`Command::CancelAnalysis`].
    /// (The one-shot `App::analyze_project` is still used headless/by tests.)
    DeepAnalyzeProject(i64),
    /// Load every biosample (subjects list), regardless of project.
    LoadAllBiosamples,
    /// Load donor-level Y/mt terminal haplogroups for every subject (fills the list columns).
    LoadHaploSummary,
    AddBiosample(NewBiosample),
    /// Batch-import a NAS project directory (scan → Project/Biosample/Run/Alignment).
    /// `reference` is optional: `None` lets the gateway resolve each build from the cache
    /// (and report `ReferenceNeeded` if a download is required); `Some` pins a FASTA.
    ImportProjectDir {
        dir: PathBuf,
        reference: Option<PathBuf>,
    },
    /// Dry-run an FTDNA project import: parse + match the (already classified) batch files into a
    /// reviewable plan. Any path may be absent.
    PlanFtdnaImport {
        /// Target project, or `None` to import into a new project named `project_name`.
        project_id: Option<i64>,
        project_name: Option<String>,
        member: Option<PathBuf>,
        paternal: Option<PathBuf>,
        maternal: Option<PathBuf>,
        ystr: Option<PathBuf>,
    },
    /// Commit a reviewed FTDNA import plan with the admin's per-kit resolutions.
    CommitFtdnaImport {
        plan: FtdnaImportPlan,
        resolutions: std::collections::BTreeMap<String, FtdnaResolution>,
    },
    /// Load a subject's imported genealogy (vendor ids + FTDNA member + MDKA) for the detail card.
    LoadGenealogy(SampleGuid),
    /// Autocluster a project's members by Y-STR (suggest SNP branches for STR-only members).
    ClusterProject(i64),
    /// Resolve (download + decompress + index) a reference build, streaming progress.
    ResolveReference {
        build: String,
    },
    LoadRuns(SampleGuid),
    AddRun(NewSequenceRun),
    /// Load the donor-level Y + mtDNA haplogroup consensus for a subject.
    LoadConsensus(SampleGuid),
    LoadStrProfiles(SampleGuid),
    /// Call Y-STRs from the subject's best sequence alignment and compare to the imported vendor
    /// profile (the By-Panel concordance). Heavy on first call (a chrY pass); cached after.
    StrConcordance {
        biosample_guid: SampleGuid,
    },
    /// Rank every other workspace subject against this one by Y relatedness (gap §2). One-vs-all
    /// over the workspace, or one project when `project_id` is set. Consumes cached profiles.
    YMatches {
        biosample_guid: SampleGuid,
        project_id: Option<i64>,
    },
    ImportStrProfile {
        biosample_guid: SampleGuid,
        panel_name: String,
        provider: Option<String>,
        source: Option<String>,
        path: PathBuf,
    },
    LoadVariantSets(SampleGuid),
    ImportVariants {
        biosample_guid: SampleGuid,
        path: PathBuf,
        source_type: SourceType,
    },
    /// Manually-entered variant calls (e.g. Sanger/YSEQ confirmations): `contig,pos,ref,alt` rows.
    AddVariants {
        biosample_guid: SampleGuid,
        source_label: String,
        source_type: SourceType,
        text: String,
    },
    LoadChipProfiles(SampleGuid),
    ImportChipProfile {
        biosample_guid: SampleGuid,
        provider: Option<String>,
        path: PathBuf,
    },
    LoadMtdna(SampleGuid),
    ImportMtdna {
        biosample_guid: SampleGuid,
        path: PathBuf,
    },
    /// Derive mtDNA variants for a stored sequence vs an rCRS reference FASTA.
    /// Derive the mtDNA mutation list (vs the bundled rCRS) for display.
    LoadMtdnaVariants {
        mtdna_id: i64,
    },
    /// Assign an mtDNA haplogroup (fetch the FTDNA tree, rank by the sample's base calls).
    AssignMtdnaHaplogroup {
        mtdna_id: i64,
    },
    /// Assign a Y haplogroup from an alignment (call chrY tree positions, rank).
    AssignYHaplogroup {
        alignment_id: i64,
    },
    /// Full Y placement report: ranked candidates + lineage SNP evidence (gap §8 haplogroup report).
    YHaploReport {
        alignment_id: i64,
    },
    /// Assign a Y haplogroup from the subject's imported BISDNA / Y-SNP-panel calls (no
    /// alignment) — records a donor call.
    AssignYBisdna {
        biosample_guid: SampleGuid,
    },
    /// Assign an mtDNA haplogroup directly from an alignment's chrM (records a donor call).
    AssignMtdnaHaplogroupFromAlignment {
        alignment_id: i64,
    },
    /// Estimate autosomal ancestry from the subject's CONSENSUS (no BAM walk) — the default path.
    EstimateAncestryFromConsensus {
        biosample_guid: SampleGuid,
    },
    /// Paint local ancestry from the subject's CONSENSUS (no BAM walk).
    PaintAncestryFromConsensus {
        biosample_guid: SampleGuid,
    },
    /// Load the cached chromosome painting (if current for the consensus signature) — cheap.
    LoadPainting {
        biosample_guid: SampleGuid,
    },
    /// Load the cached detailed consensus ancestry reports (modern fine + ancient components).
    LoadConsensusAncestryDetail {
        biosample_guid: SampleGuid,
    },
    /// Find the private bucket: de-novo chrY calls off the assigned Y backbone, restricted
    /// by the chosen callable mask.
    FindPrivateY {
        alignment_id: i64,
        mask: YMask,
    },
    /// Load a previously-computed (self-masked) private-Y bucket from cache.
    LoadPrivateY {
        alignment_id: i64,
    },
    /// Unified import: multiple files and/or folders (folders walked for data files), each
    /// auto-detected + routed; returns one [`Event::DataBatchImported`] summary.
    AddDataBatch {
        biosample_guid: SampleGuid,
        paths: Vec<PathBuf>,
    },
    LoadAlignments(i64),
    AddAlignment(NewAlignment),
    /// Resolve the subject's default analysis alignment (highest-coverage, else first).
    DefaultAlignment {
        biosample_guid: SampleGuid,
    },
    /// Load the subject's donor-level ancestry (best estimate across all sources).
    LoadDonorAncestry {
        biosample_guid: SampleGuid,
    },
    /// Load the subject's donor-level private-Y union across all sources.
    LoadDonorPrivateY {
        biosample_guid: SampleGuid,
    },
    /// Load the subject's multi-source Y-variant profile (concordance across all Y sources).
    /// Load the persisted Y-profile snapshot (cheap; no genotyping).
    LoadYProfile {
        biosample_guid: SampleGuid,
    },
    /// Recompute the Y-profile from all sources and persist the snapshot (expensive — re-genotypes).
    BuildYProfile {
        biosample_guid: SampleGuid,
    },
    /// Resolve catalogued Y-SNP names at the given positions (annotates the Y-SNP tables).
    LoadYSnpNames {
        biosample_guid: SampleGuid,
        positions: Vec<i64>,
    },
    /// Load the persisted mtDNA consensus-profile snapshot (cheap; no genotyping).
    LoadMtProfile {
        biosample_guid: SampleGuid,
    },
    /// Recompute the mtDNA consensus profile from all sources and persist (expensive — re-places).
    BuildMtProfile {
        biosample_guid: SampleGuid,
    },
    /// Load the persisted autosomal consensus-profile snapshot (cheap; no genotyping).
    LoadAutosomalProfile {
        biosample_guid: SampleGuid,
    },
    /// Recompute the autosomal consensus from all sources and persist (expensive — panel-genotypes).
    BuildAutosomalProfile {
        biosample_guid: SampleGuid,
    },
    /// Probe a BAM/CRAM header for build/aligner/platform/test-type (to auto-fill the form).
    ProbeAlignment {
        path: PathBuf,
    },
    LoadCoverage(i64),
    /// Cached coverage for several alignments at once (Data Sources alignment rows).
    LoadCoverageBulk(Vec<i64>),
    /// Genome-region metadata (cytoband ideogram) for an alignment's build (Ideogram tab).
    LoadGenomeRegions {
        alignment_id: i64,
        build: String,
    },
    RunCoverage(i64),
    LoadSex(i64),
    RunSex(i64),
    LoadReadMetrics(i64),
    RunReadMetrics(i64),
    LoadSv(i64),
    RunSv(i64),
    LoadDenovo {
        alignment_id: i64,
        contig: String,
    },
    RunDenovo {
        alignment_id: i64,
        contig: String,
    },
    LoadPanels,
    ImportPanel {
        name: String,
        path: PathBuf,
    },
    LoadAllAlignments,
    GenotypePanel {
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
    },
    LoadPanelGenotypes {
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
    },
    CompareIbd {
        a: i64,
        b: i64,
        panel_id: i64,
        ploidy: u8,
    },
    /// Compare two samples (each a WGS alignment or an imported chip) over the chip-compatible IBD
    /// panel — the volume-case path (chip↔chip / chip↔WGS).
    CompareIbdSources {
        a: navigator_app::IbdSource,
        b: navigator_app::IbdSource,
    },
    /// Compare two SUBJECTS over their autosomal consensuses (the subject-level IBD path).
    CompareIbdConsensus {
        a: SampleGuid,
        b: SampleGuid,
    },
    /// Verify two alignments are the same individual (genotype concordance + Y-STR).
    VerifyIdentity {
        a: i64,
        b: i64,
        panel_id: i64,
        ploidy: u8,
    },
    /// Verify two SUBJECTS are the same individual over their pooled autosomal consensus (no panel).
    VerifyIdentityConsensus {
        a: SampleGuid,
        b: SampleGuid,
    },
    /// Federated IBD step 1: fetch the AppView's pseudonymous match suggestions for the
    /// signed-in account (registers the device key on first use).
    LoadIbdSuggestions,
    /// Federated IBD step 2: request an introduction to a suggested candidate.
    IbdIntroduce {
        suggested_sample_guid: String,
    },
    /// Adopt a local self-certifying did:key identity (desktop bootstrap — no PDS/OAuth).
    UseLocalIdentity,
    /// Poll the AppView for inbound exchange requests + consent-ready sessions (the exchange inbox).
    ExchangeInbox,
    /// Consent to (or decline) an inbound exchange request.
    ExchangeConsent {
        request_uri: String,
        given: bool,
    },
    /// Run a full IBD exchange for a subject over a consent-ready session (handshake → dosage
    /// exchange → signed attestations → persist). Long-running; needs the peer online.
    RunIbdExchange {
        info: ExchangeSessionInfo,
        biosample_guid: SampleGuid,
    },
    /// Load the subject's persisted IBD exchange results.
    LoadIbdExchanges {
        biosample_guid: SampleGuid,
    },
    /// Resolve the sequencing lab for runs that have an inferred instrument id but no facility,
    /// via the AppView instrument→lab map (best-effort, cached). Sent on startup + after imports.
    BackfillLabs,
    /// Report who's signed in (no side effects) — sent on startup.
    AuthStatus,
    /// Report the current online/offline state (no side effects).
    SyncStatus,
    /// Sign in to a PDS via OAuth (opens a browser); `handle` is a handle or DID.
    Login {
        handle: String,
    },
    Logout,
    PublishCoverage(i64),
    PublishVariants {
        alignment_id: i64,
        contig: String,
    },
    /// Publish the subject's consensus ancestry breakdown (one record per method) to the signed-in PDS.
    PublishAncestry {
        biosample_guid: SampleGuid,
    },
    /// Attempt to push the ready outbox rows now (also runs periodically + after a publish).
    DrainOutbox,
    /// PULL reconcile: fetch the account's PDS records and reconcile against local (gap §5-p2).
    PullSync,
    /// Re-check tracked source files' accessibility (moved/missing).
    VerifySourceFiles,
    /// Write the IBD match's segments to a CSV/TSV at `path` (the match-browser export).
    ExportIbdSegments {
        segments: Vec<navigator_app::IbdSegment>,
        path: PathBuf,
    },
    /// Load the reference population PC1/PC2 centroids for an alignment's build (PCA scatter backdrop).
    LoadPcaReference,
    /// Export a cached result to `path` (TSV/HTML/BED). `request` carries the kind + source id.
    Export {
        request: navigator_app::ExportRequest,
        path: PathBuf,
    },
    /// Manually override the consensus haplogroup for a subject + DNA type.
    SetHaploOverride {
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        haplogroup: String,
        reason: Option<String>,
    },
    /// Clear a manual override.
    ClearHaploOverride {
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    },
    /// Load the reconciliation audit log for a subject + DNA type.
    LoadAudit {
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    },
    /// Scan an alignment's chrM pileup for heteroplasmic positions.
    LoadHeteroplasmy {
        alignment_id: i64,
    },
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
    RunFullAnalysis {
        alignment_id: i64,
    },
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
    DeleteSequenceRun {
        id: i64,
        biosample_guid: SampleGuid,
    },
    /// Merge `secondary` sequence run into `primary` (reparent alignments, delete the empty run).
    MergeSequenceRuns {
        biosample_guid: SampleGuid,
        primary: i64,
        secondary: i64,
    },
    /// Delete a single alignment (cascades to its artifacts). `sequence_run_id` is the owner,
    /// so the UI can refresh that run's alignment list.
    DeleteAlignment {
        id: i64,
        sequence_run_id: i64,
    },
    /// Delete an imported STR profile (and its markers).
    DeleteStrProfile {
        id: i64,
        biosample_guid: SampleGuid,
    },
    /// Delete an imported variant set (and its calls).
    DeleteVariantSet {
        id: i64,
        biosample_guid: SampleGuid,
    },
    /// Delete an imported chip/array profile.
    DeleteChipProfile {
        id: i64,
        biosample_guid: SampleGuid,
    },
    /// Delete an imported mtDNA sequence.
    DeleteMtdnaSequence {
        id: i64,
        biosample_guid: SampleGuid,
    },
    /// Assign a subject to a project (`None` clears it). The app layer validates the project.
    AssignBiosampleProject {
        guid: SampleGuid,
        project_id: Option<i64>,
    },
    /// Update a project's editable fields.
    UpdateProject {
        id: i64,
        name: String,
        description: Option<String>,
        administrator: String,
    },
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
        /// The sequencing lab (a `labs` display name); `None` clears it.
        sequencing_facility: Option<String>,
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
    /// Load per-build reference-genome settings + cache status for the Settings dialog.
    LoadReferenceSettings,
    /// Set a build's local-FASTA override + auto-download flag (persists reference_sources.json).
    SetReferenceOverride {
        build: String,
        local_path: Option<String>,
        auto_download: bool,
    },
    /// Re-hash a cached reference against its integrity sidecar (Settings "Verify").
    VerifyReference {
        build: String,
    },
    /// Lift a VCF from `source` (inferred when `None`) to `target` build, writing `out`.
    LiftVcf {
        source: Option<String>,
        target: String,
        in_vcf: PathBuf,
        out_vcf: PathBuf,
        filter_par: bool,
    },
    // ---- social (Community tab — signed AppView Edge API) -------------------
    /// List the signed-in account's support threads (team↔tester).
    LoadSupportThreads,
    /// Read one support thread's messages (marks it read server-side).
    LoadSupportThread {
        conversation_id: String,
    },
    /// Open a new support thread to the team.
    OpenSupportThread {
        subject: String,
        body: String,
    },
    /// Reply to an existing support thread.
    ReplySupportThread {
        conversation_id: String,
        body: String,
    },
    /// Read the community feed (announcements + community + federated).
    LoadCommunityFeed,
    /// Post to the community feed (optionally tagged with a topic).
    PostCommunity {
        content: String,
        topic: Option<String>,
    },
    /// Fetch notifications + unread count (also drives the app-bar bell).
    LoadNotifications,
    /// Mark one notification read (`Some`) or all (`None`).
    MarkNotificationRead {
        id: Option<String>,
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
    /// Ancestry/IBD asset presence + integrity (the "data sources" transparency line).
    AssetStatus(Vec<navigator_app::AssetStatus>),
    ProjectCreated(Project),
    /// A project was updated or deleted; reload the overview.
    ProjectsChanged,
    /// A batch project-directory import completed.
    ProjectImported(ProjectImportSummary),
    /// A dry-run FTDNA import plan, ready for the review modal.
    FtdnaPlan(FtdnaImportPlan),
    /// The result of committing an FTDNA import.
    FtdnaImported(FtdnaImportSummary),
    /// A subject's imported genealogy bundle for the detail card.
    Genealogy {
        guid: SampleGuid,
        data: FtdnaGenealogy,
    },
    /// A project's Y-STR clustering (members grouped by branch, STR-only suggestions).
    ProjectClustering {
        project_id: i64,
        clustering: YstrClustering,
    },
    /// Import needs reference build(s) downloaded first; `dir` lets the UI retry the import
    /// after the user approves and the download finishes.
    ReferenceNeeded {
        dir: PathBuf,
        builds: Vec<BuildNeed>,
    },
    /// Bytes received during a reference download (`total` from Content-Length, if known).
    ReferenceProgress {
        build: String,
        received: u64,
        total: Option<u64>,
    },
    /// A reference build finished resolving (cached + indexed).
    ReferenceReady {
        build: String,
        path: PathBuf,
    },
    /// Per-sample coverage/haplogroup report for a project.
    ProjectReport {
        project_id: i64,
        rows: Vec<ProjectSampleReport>,
    },
    /// A project-wide analyze pass finished (coverage + Y per sample). `cancelled` is true when a
    /// streaming deep-analyze was stopped early (counts reflect what completed before the stop).
    ProjectAnalyzed {
        project_id: i64,
        samples: usize,
        coverage_done: usize,
        y_done: usize,
        sex_done: usize,
        metrics_done: usize,
        sv_done: usize,
        errors: usize,
        cancelled: bool,
    },
    /// Per-sample progress of a streaming deep-analyze pass: `done` of `total` samples processed,
    /// `sample` is the donor id currently being analyzed, `fraction` drives the bar (0..1).
    DeepAnalyzeProgress {
        project_id: i64,
        done: usize,
        total: usize,
        sample: String,
        fraction: f32,
    },
    Samples {
        project_id: i64,
        samples: Vec<Biosample>,
    },
    /// All biosamples (the project-independent subjects list).
    AllBiosamples(Vec<Biosample>),
    /// Per-subject Y/mt terminal haplogroups for the subjects list (`guid → (Y, mt)`).
    HaploSummary(std::collections::HashMap<SampleGuid, (Option<String>, Option<String>)>),
    /// A biosample was added/changed; reload the subjects list (and any open project view).
    BiosamplesChanged,
    Runs {
        biosample_guid: SampleGuid,
        runs: Vec<SequenceRun>,
    },
    RunsChanged(SampleGuid),
    /// Donor-level haplogroup consensus for a subject (Y, mtDNA).
    Consensus {
        biosample_guid: SampleGuid,
        y: Option<Consensus>,
        mt: Option<Consensus>,
    },
    StrProfiles {
        biosample_guid: SampleGuid,
        profiles: Vec<StrProfile>,
    },
    StrProfilesChanged(SampleGuid),
    VariantSets {
        biosample_guid: SampleGuid,
        sets: Vec<VariantSet>,
    },
    VariantSetsChanged(SampleGuid),
    ChipProfiles {
        biosample_guid: SampleGuid,
        profiles: Vec<ChipProfile>,
    },
    ChipProfilesChanged(SampleGuid),
    MtdnaSequences {
        biosample_guid: SampleGuid,
        sequences: Vec<MtdnaSequence>,
    },
    MtdnaChanged(SampleGuid),
    /// The rCRS-relative mutation list for an mtDNA sequence.
    MtdnaVariants {
        mtdna_id: i64,
        variants: Vec<navigator_app::MtVariant>,
    },
    /// Y-STR concordance: markers called from sequence (FTDNA-convention) vs the imported vendor
    /// profile, for the By-Panel view. `alignment_id` is the source alignment chosen.
    StrConcordance {
        biosample_guid: SampleGuid,
        alignment_id: i64,
        rows: Vec<StrConcordanceRow>,
    },
    /// Cross-subject Y matches for a subject, ranked best-first (gap §2).
    YMatches {
        biosample_guid: SampleGuid,
        matches: Vec<YMatch>,
    },
    /// A batch import finished; the summary lists per-file imported/skipped outcomes. The UI
    /// shows it in a modal and reloads the subject's data sections.
    DataBatchImported {
        biosample_guid: SampleGuid,
        summary: BatchImportSummary,
    },
    /// mtDNA haplogroup assignment for a sequence (ranked + terminal evidence).
    Haplogroup {
        mtdna_id: i64,
        assignment: HaploAssignment,
    },
    /// Y haplogroup assignment for an alignment (ranked + terminal evidence).
    YHaplogroup {
        alignment_id: i64,
        assignment: HaploAssignment,
    },
    /// Full Y placement report: ranked candidates + lineage SNP evidence.
    YHaploReport {
        alignment_id: i64,
        assignment: HaploAssignment,
        lineage: Vec<navigator_app::SnpEvidence>,
    },
    /// Y haplogroup assignment from the subject's BISDNA / Y-SNP panel (records a donor call).
    YBisdnaHaplogroup {
        biosample_guid: SampleGuid,
        assignment: HaploAssignment,
    },
    /// mtDNA haplogroup assignment from an alignment (records a donor call → reload consensus).
    MtHaplogroup {
        alignment_id: i64,
        assignment: HaploAssignment,
    },
    /// Local-ancestry segments per chromosome (the "DNA painting").
    AncestryPainting {
        alignment_id: i64,
        segments: Vec<AncestrySegment>,
    },
    /// Private Y variants (off-backbone de-novo calls) for an alignment.
    PrivateY {
        alignment_id: i64,
        bucket: PrivateBucket,
    },
    Alignments {
        sequence_run_id: i64,
        alignments: Vec<Alignment>,
    },
    AlignmentsChanged(i64),
    /// The subject's default analysis alignment, to auto-select on the detail tabs.
    DefaultAlignment {
        run_id: i64,
        alignment_id: i64,
    },
    /// Donor-level ancestry (best across sources) + the source alignment it came from.
    DonorAncestry {
        alignment_id: i64,
        result: AncestryResult,
    },
    /// Donor-level private-Y union across the subject's sources.
    DonorPrivateY {
        bucket: PrivateBucket,
    },
    /// The subject's multi-source Y-variant profile.
    YProfile {
        biosample_guid: SampleGuid,
        profile: Option<navigator_app::YProfile>,
    },
    /// Catalogued Y-SNP names at requested positions (`position → name`) for the Y-SNP tables.
    YSnpNames {
        names: std::collections::HashMap<i64, String>,
    },
    /// The subject's multi-source mtDNA consensus profile.
    MtProfile {
        biosample_guid: SampleGuid,
        profile: Option<navigator_app::ConsensusProfile>,
    },
    /// The subject's multi-source autosomal consensus profile (diploid 0/1/2).
    AutosomalProfile {
        biosample_guid: SampleGuid,
        profile: Option<navigator_app::DiploidProfile>,
    },
    /// Detailed consensus ancestry reports: modern fine-population + ancient-component breakdowns.
    ConsensusAncestryDetail {
        biosample_guid: SampleGuid,
        // Boxed: AncestryResult is large, and three of them would bloat the Event enum's size.
        fine: Option<Box<navigator_app::AncestryResult>>,
        ancient: Option<Box<navigator_app::AncestryResult>>,
        nmonte: Option<Box<navigator_app::AncestryResult>>,
    },
    /// Header-probe result for the add-alignment form (build/aligner/platform/test-type).
    AlignmentProbe(AlignmentProbe),
    Coverage {
        alignment_id: i64,
        result: Option<Coverage>,
    },
    /// Cached coverage for several alignments (Data Sources rows): `(alignment_id, result)`.
    CoverageBulk(Vec<(i64, Option<Coverage>)>),
    /// Genome-region metadata (cytoband ideogram) for an alignment's build.
    GenomeRegions {
        alignment_id: i64,
        regions: Option<std::sync::Arc<navigator_app::GenomeRegions>>,
    },
    Sex {
        alignment_id: i64,
        result: Option<SexInferenceResult>,
    },
    ReadMetrics {
        alignment_id: i64,
        result: Option<ReadMetrics>,
    },
    Sv {
        alignment_id: i64,
        result: Option<SvAnalysisResult>,
    },
    Denovo {
        alignment_id: i64,
        contig: String,
        result: Option<Vec<DenovoCall>>,
    },
    /// Full-analysis pipeline progress: starting `step` of `total` (1-based), with a `label`
    /// + `detail` and the bar `fraction` (0..1).
    AnalysisProgress {
        step: usize,
        total: usize,
        label: String,
        detail: String,
        fraction: f32,
    },
    /// The full-analysis pipeline finished (or was cancelled).
    AnalysisDone {
        cancelled: bool,
    },
    Panels(Vec<PanelInfo>),
    PanelImported,
    AllAlignments(Vec<Alignment>),
    PanelGenotypes {
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
        genotypes: Vec<PanelGenotype>,
    },
    Ibd(IbdComparison),
    /// Federated IBD match suggestions from the AppView (may be empty in a single-user dev AppView).
    IbdSuggestions(Vec<IbdSuggestion>),
    /// An introduction request was opened for a candidate (status initially `PENDING`).
    IbdIntroduced {
        suggested_sample_guid: String,
        request_uri: String,
        status: String,
    },
    /// The exchange inbox: inbound requests awaiting our consent + consent-ready sessions.
    ExchangeInbox {
        incoming: Vec<IncomingRequest>,
        ready: Vec<ExchangeSessionInfo>,
    },
    /// A consent was recorded (the UI refreshes the inbox).
    ExchangeConsented,
    /// A full IBD exchange completed for a subject (the UI reloads its results).
    IbdExchangeDone {
        biosample_guid: SampleGuid,
        total_shared_cm: f64,
        segment_count: usize,
        relationship: String,
        agreed: bool,
    },
    /// The subject's persisted IBD exchange results.
    IbdExchanges {
        biosample_guid: SampleGuid,
        rows: Vec<StoredIbdExchange>,
    },
    /// Identity-verification result between two alignments.
    Identity(IdentityVerification),
    /// The reconciliation audit log for a subject + DNA type.
    Audit {
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        entries: Vec<AuditEntry>,
    },
    /// mtDNA heteroplasmy sites for an alignment.
    Heteroplasmy {
        alignment_id: i64,
        sites: Vec<HeteroplasmySite>,
    },
    /// A reconciliation override/clear succeeded; the UI reloads consensus + audit.
    ReconciliationChanged {
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    },
    /// Current signed-in account (DID), or `None` when signed out.
    Authenticated(Option<String>),
    /// A record was published; `kind` is a human label, `uri` the `at://` URI.
    Published {
        kind: String,
        uri: String,
    },
    /// A publish was enqueued to the durable outbox (it'll send now if online, else on reconnect).
    Queued {
        kind: String,
    },
    /// Outbox rows still awaiting a successful push (the "N pending" indicator).
    SyncPending(i64),
    /// A result was exported to `path`; `label` is the human kind (e.g. "coverage (TSV)").
    Exported {
        label: String,
        path: PathBuf,
    },
    /// Whether the last PDS write reached the server (offline indicator).
    SyncOnline(bool),
    /// A PULL reconcile finished (gap §5-p2): the per-action tallies.
    PullDone {
        in_sync: usize,
        applied: usize,
        adopted: usize,
        repushed: usize,
        conflicts: usize,
    },
    /// Source-file accessibility re-check finished; `missing` files are moved/deleted.
    SourceFilesVerified {
        missing: usize,
    },
    /// Reference population PC1/PC2 centroids for the PCA scatter: `(population_code, pc1, pc2)`.
    PcaReference {
        alignment_id: i64,
        points: Vec<(String, f64, f64)>,
    },
    /// How many runs had their sequencing lab filled in by the AppView backfill (`0` ⇒ quiet).
    LabsResolved(usize),
    /// Per-build reference-genome settings + cache status for the Settings dialog.
    ReferenceSettings(Vec<RefBuildStatus>),
    /// A reference override was saved; the UI may reload the settings rows.
    ReferenceSettingsChanged,
    /// Result of a reference integrity check (a short human-readable status per build).
    ReferenceVerified {
        build: String,
        status: String,
    },
    /// A VCF liftover finished (a human-readable stats summary), for the status line.
    VcfLifted {
        summary: String,
    },
    // ---- social (Community tab) --------------------------------------------
    /// The signed-in account's support threads.
    SupportThreads(Vec<navigator_app::SocialThreadSummary>),
    /// One support thread's messages (with its conversation id).
    SupportThread {
        conversation_id: String,
        messages: Vec<navigator_app::SocialMessage>,
    },
    /// A support thread was opened/replied; reload the list (+ the open thread).
    SupportThreadPosted {
        conversation_id: String,
    },
    /// The community feed.
    CommunityFeed(navigator_app::FeedView),
    /// A community post succeeded; reload the feed.
    CommunityPosted,
    /// Notifications + unread count (drives the app-bar bell badge).
    Notifications {
        items: Vec<navigator_app::SocialNotification>,
        unread: i64,
    },
    /// Notifications were marked read; reload them.
    NotificationsMarked,
    Error(String),
}

/// Execute one command against the app, mapping success/failure to an [`Event`].
pub async fn handle(app: &App, cmd: Command) -> Event {
    match cmd {
        Command::LoadOverview => match app.project_overview().await {
            Ok(v) => Event::Overview(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAssetStatus => Event::AssetStatus(navigator_app::ancestry_asset_status()),
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
        Command::PlanFtdnaImport {
            project_id,
            project_name,
            member,
            paternal,
            maternal,
            ystr,
        } => match app
            .plan_ftdna_import(
                project_id,
                project_name,
                member,
                paternal,
                maternal,
                ystr,
                FtdnaImportOptions::default(),
            )
            .await
        {
            Ok(plan) => Event::FtdnaPlan(plan),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::CommitFtdnaImport { plan, resolutions } => match app.commit_ftdna_import(&plan, &resolutions).await {
            Ok(summary) => Event::FtdnaImported(summary),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadGenealogy(guid) => match app.subject_genealogy(guid).await {
            Ok(data) => Event::Genealogy { guid, data },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ClusterProject(project_id) => match app.cluster_project_ystr(project_id).await {
            Ok(clustering) => Event::ProjectClustering { project_id, clustering },
            Err(e) => Event::Error(e.to_string()),
        },
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
        // DeepAnalyzeProject streams DeepAnalyzeProgress from the spawn loop; reaching here is a bug.
        Command::DeepAnalyzeProject(project_id) => {
            Event::Error(format!("internal: unrouted DeepAnalyzeProject {project_id}"))
        }
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
        Command::UpdateBiosample {
            guid,
            donor_identifier,
            sample_accession,
            description,
            center_name,
            sex,
        } => {
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
        Command::MergeSequenceRuns {
            biosample_guid,
            primary,
            secondary,
        } => match app.merge_sequence_runs(biosample_guid, primary, secondary).await {
            Ok(_) => Event::RunsChanged(biosample_guid),
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
        Command::LoadReferenceSettings => Event::ReferenceSettings(app.reference_settings()),
        Command::SetReferenceOverride {
            build,
            local_path,
            auto_download,
        } => match app.set_reference_override(&build, local_path, auto_download) {
            Ok(()) => Event::ReferenceSettingsChanged,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::VerifyReference { build } => match app.verify_reference(&build).await {
            Ok(outcome) => {
                let status = match outcome {
                    navigator_app::VerifyOutcome::Verified => "✓ verified".to_string(),
                    navigator_app::VerifyOutcome::Mismatch { .. } => "✗ mismatch (corrupted?)".to_string(),
                    navigator_app::VerifyOutcome::NoSidecar => "• no checksum on record".to_string(),
                    navigator_app::VerifyOutcome::NotCached => "not cached".to_string(),
                };
                Event::ReferenceVerified { build, status }
            }
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LiftVcf {
            source,
            target,
            in_vcf,
            out_vcf,
            filter_par,
        } => {
            let source = source.or_else(|| navigator_app::infer_vcf_source_build(&in_vcf));
            match source {
                None => Event::Error(format!(
                    "could not infer the source build of {} — set it explicitly",
                    in_vcf.display()
                )),
                Some(src) => {
                    let opts = navigator_app::VcfLiftOpts { filter_par };
                    match app
                        .lift_vcf(&src, &target, in_vcf, out_vcf.clone(), opts, &mut |_, _| {})
                        .await
                    {
                        Ok(s) => Event::VcfLifted {
                            summary: format!(
                                "Lifted {}/{} variants ({} unmapped, {} ref-mismatch) → {}",
                                s.lifted,
                                s.total,
                                s.unmapped,
                                s.ref_mismatch,
                                out_vcf.display()
                            ),
                        },
                        Err(e) => Event::Error(e.to_string()),
                    }
                }
            }
        }
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
        Command::UpdateProject {
            id,
            name,
            description,
            administrator,
        } => match app.update_project(id, name, description, administrator).await {
            Ok(_) => Event::ProjectsChanged,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteProject(id) => match app.delete_project(id).await {
            Ok(()) => Event::ProjectsChanged,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::UpdateSequenceRun {
            id,
            biosample_guid,
            platform_name,
            instrument_model,
            test_type,
            library_layout,
            sequencing_facility,
        } => {
            match app
                .update_sequence_run(
                    id,
                    platform_name,
                    instrument_model,
                    test_type,
                    library_layout,
                    sequencing_facility,
                )
                .await
            {
                Ok(_) => Event::RunsChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::UpdateAlignment {
            id,
            sequence_run_id,
            reference_build,
            aligner,
            variant_caller,
        } => match app.update_alignment(id, reference_build, aligner, variant_caller).await {
            Ok(_) => Event::AlignmentsChanged(sequence_run_id),
            Err(e) => Event::Error(e.to_string()),
        },
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
            Event::Consensus {
                biosample_guid: guid,
                y,
                mt,
            }
        }
        Command::StrConcordance { biosample_guid } => match app.str_concordance_for_subject(biosample_guid).await {
            Ok((alignment_id, rows)) => Event::StrConcordance {
                biosample_guid,
                alignment_id,
                rows,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::YMatches {
            biosample_guid,
            project_id,
        } => match app.y_matches(biosample_guid, project_id).await {
            Ok(matches) => Event::YMatches {
                biosample_guid,
                matches,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadStrProfiles(guid) => match app.list_str_profiles(guid).await {
            Ok(profiles) => Event::StrProfiles {
                biosample_guid: guid,
                profiles,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportStrProfile {
            biosample_guid,
            panel_name,
            provider,
            source,
            path,
        } => {
            match app
                .import_str_profile_from_csv(biosample_guid, &panel_name, provider, source, &path)
                .await
            {
                Ok(_) => Event::StrProfilesChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadVariantSets(guid) => match app.list_variant_sets(guid).await {
            Ok(sets) => Event::VariantSets {
                biosample_guid: guid,
                sets,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportVariants {
            biosample_guid,
            path,
            source_type,
        } => match app.import_variants_from_file(biosample_guid, &path, source_type).await {
            Ok(_) => Event::VariantSetsChanged(biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddVariants {
            biosample_guid,
            source_label,
            source_type,
            text,
        } => {
            match app
                .add_variants(biosample_guid, &source_label, source_type, &text)
                .await
            {
                Ok(_) => Event::VariantSetsChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadChipProfiles(guid) => match app.list_chip_profiles(guid).await {
            Ok(profiles) => Event::ChipProfiles {
                biosample_guid: guid,
                profiles,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportChipProfile {
            biosample_guid,
            provider,
            path,
        } => {
            match app
                .import_chip_profile_from_csv(biosample_guid, provider, None, &path)
                .await
            {
                Ok(_) => Event::ChipProfilesChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadMtdna(guid) => match app.list_mtdna_sequences(guid).await {
            Ok(sequences) => Event::MtdnaSequences {
                biosample_guid: guid,
                sequences,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportMtdna { biosample_guid, path } => {
            match app.import_mtdna_from_fasta(biosample_guid, &path).await {
                Ok(_) => Event::MtdnaChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadMtdnaVariants { mtdna_id } => match app.mtdna_variants(mtdna_id).await {
            Ok(variants) => Event::MtdnaVariants { mtdna_id, variants },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignMtdnaHaplogroup { mtdna_id } => match app.assign_mtdna_haplogroup(mtdna_id).await {
            Ok(assignment) => Event::Haplogroup { mtdna_id, assignment },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignYBisdna { biosample_guid } => match app.assign_y_bisdna(biosample_guid, None).await {
            Ok(assignment) => Event::YBisdnaHaplogroup {
                biosample_guid,
                assignment,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::YHaploReport { alignment_id } => match app.y_haplogroup_report(alignment_id).await {
            Ok((assignment, lineage)) => Event::YHaploReport {
                alignment_id,
                assignment,
                lineage,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignYHaplogroup { alignment_id } => match app.assign_y_haplogroup(alignment_id).await {
            Ok(assignment) => Event::YHaplogroup {
                alignment_id,
                assignment,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AssignMtdnaHaplogroupFromAlignment { alignment_id } => {
            match app.assign_mtdna_haplogroup_from_alignment(alignment_id).await {
                Ok(assignment) => Event::MtHaplogroup {
                    alignment_id,
                    assignment,
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::EstimateAncestryFromConsensus { biosample_guid } => {
            // Estimate from the pooled consensus, then surface it as the donor-level result.
            match app.estimate_ancestry_from_consensus(biosample_guid).await {
                Ok(result) => Event::DonorAncestry {
                    alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                    result,
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::PaintAncestryFromConsensus { biosample_guid } => {
            // Painting from the consensus needs no genotyping pass — fast, no progress stream.
            match app.paint_local_ancestry_from_consensus(biosample_guid).await {
                Ok(segments) => Event::AncestryPainting {
                    alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                    segments,
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadPainting { biosample_guid } => match app.cached_painting(biosample_guid).await {
            Ok(segments) => Event::AncestryPainting {
                alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                segments: segments.unwrap_or_default(),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadConsensusAncestryDetail { biosample_guid } => {
            let fine = app
                .consensus_ancestry(biosample_guid, "FINE_ADMIXTURE")
                .await
                .unwrap_or(None)
                .map(Box::new);
            let ancient = app
                .consensus_ancestry(biosample_guid, "PCA_PROJECTION_GMM")
                .await
                .unwrap_or(None)
                .map(Box::new);
            let nmonte = app
                .consensus_ancestry(biosample_guid, "G25_NMONTE")
                .await
                .unwrap_or(None)
                .map(Box::new);
            Event::ConsensusAncestryDetail {
                biosample_guid,
                fine,
                ancient,
                nmonte,
            }
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
        Command::AddDataBatch { biosample_guid, paths } => {
            match app.add_data_batch(biosample_guid, paths, |_, _| {}).await {
                Ok(summary) => Event::DataBatchImported {
                    biosample_guid,
                    summary,
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadAlignments(sequence_run_id) => match app.list_alignments(sequence_run_id).await {
            Ok(alignments) => Event::Alignments {
                sequence_run_id,
                alignments,
            },
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
        Command::LoadYProfile { biosample_guid } => match app.cached_y_profile(biosample_guid).await {
            Ok(profile) => Event::YProfile {
                biosample_guid,
                profile,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::BuildYProfile { biosample_guid } => match app.build_y_profile(biosample_guid).await {
            Ok(profile) => Event::YProfile {
                biosample_guid,
                profile: Some(profile),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadYSnpNames {
            biosample_guid,
            positions,
        } => match app.y_snp_names_at(biosample_guid, &positions).await {
            Ok(names) => Event::YSnpNames { names },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadMtProfile { biosample_guid } => match app.cached_mt_profile(biosample_guid).await {
            Ok(profile) => Event::MtProfile {
                biosample_guid,
                profile,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::BuildMtProfile { biosample_guid } => match app.build_mt_profile(biosample_guid).await {
            Ok(profile) => Event::MtProfile {
                biosample_guid,
                profile: Some(profile),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAutosomalProfile { biosample_guid } => match app.cached_autosomal_profile(biosample_guid).await {
            Ok(profile) => Event::AutosomalProfile {
                biosample_guid,
                profile,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::BuildAutosomalProfile { biosample_guid } => match app.build_autosomal_profile(biosample_guid).await {
            Ok(profile) => Event::AutosomalProfile {
                biosample_guid,
                profile: Some(profile),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadCoverage(alignment_id) => match app.cached_coverage(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadCoverageBulk(ids) => match app.cached_coverage_bulk(&ids).await {
            Ok(results) => Event::CoverageBulk(results),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadGenomeRegions { alignment_id, build } => match app.genome_regions(&build).await {
            Ok(regions) => Event::GenomeRegions {
                alignment_id,
                regions: Some(regions),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunCoverage(alignment_id) => match app.run_coverage_for_alignment(alignment_id).await {
            Ok(result) => Event::Coverage {
                alignment_id,
                result: Some(result),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSex(alignment_id) => match app.cached_sex(alignment_id).await {
            Ok(result) => Event::Sex { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunSex(alignment_id) => match app.run_sex(alignment_id).await {
            Ok(result) => Event::Sex {
                alignment_id,
                result: Some(result),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadReadMetrics(alignment_id) => match app.cached_read_metrics(alignment_id).await {
            Ok(result) => Event::ReadMetrics { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunReadMetrics(alignment_id) => match app.run_read_metrics(alignment_id).await {
            Ok(result) => Event::ReadMetrics {
                alignment_id,
                result: Some(result),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSv(alignment_id) => match app.cached_sv(alignment_id).await {
            Ok(result) => Event::Sv { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunSv(alignment_id) => match app.run_sv(alignment_id).await {
            Ok(result) => Event::Sv {
                alignment_id,
                result: Some(result),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadDenovo { alignment_id, contig } => match app.cached_denovo(alignment_id, &contig).await {
            Ok(result) => Event::Denovo {
                alignment_id,
                contig,
                result,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunDenovo { alignment_id, contig } => {
            match app.run_denovo_for_alignment(alignment_id, contig.clone()).await {
                Ok(result) => Event::Denovo {
                    alignment_id,
                    contig,
                    result: Some(result),
                },
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
        Command::GenotypePanel {
            alignment_id,
            panel_id,
            ploidy,
        } => match app.genotype_panel(alignment_id, panel_id, ploidy).await {
            Ok(genotypes) => Event::PanelGenotypes {
                alignment_id,
                panel_id,
                ploidy,
                genotypes,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadPanelGenotypes {
            alignment_id,
            panel_id,
            ploidy,
        } => match app.cached_panel_genotypes(alignment_id, panel_id, ploidy).await {
            Ok(genotypes) => Event::PanelGenotypes {
                alignment_id,
                panel_id,
                ploidy,
                genotypes: genotypes.unwrap_or_default(),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::CompareIbd { a, b, panel_id, ploidy } => {
            match app
                .compare_ibd(a, b, panel_id, ploidy, IbdDetectorConfig::default())
                .await
            {
                Ok(cmp) => Event::Ibd(cmp),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::CompareIbdConsensus { a, b } => {
            match app.compare_ibd_consensus(a, b, IbdDetectorConfig::default()).await {
                Ok(cmp) => Event::Ibd(cmp),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::CompareIbdSources { a, b } => {
            match app.compare_ibd_sources(a, b, IbdDetectorConfig::default()).await {
                Ok(cmp) => Event::Ibd(cmp),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::VerifyIdentity { a, b, panel_id, ploidy } => match app.verify_identity(a, b, panel_id, ploidy).await {
            Ok(v) => Event::Identity(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::VerifyIdentityConsensus { a, b } => match app.verify_identity_consensus(a, b).await {
            Ok(v) => Event::Identity(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadIbdSuggestions => match app.ibd_suggestions().await {
            Ok(items) => Event::IbdSuggestions(items),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::IbdIntroduce { suggested_sample_guid } => match app.ibd_introduce(&suggested_sample_guid).await {
            Ok(r) => Event::IbdIntroduced {
                suggested_sample_guid,
                request_uri: r.request_uri,
                status: r.status,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::UseLocalIdentity => match app.use_local_identity() {
            Ok(did) => Event::Authenticated(Some(did)),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ExchangeInbox => match (app.exchange_incoming().await, app.exchange_pending().await) {
            (Ok(incoming), Ok(ready)) => Event::ExchangeInbox { incoming, ready },
            (Err(e), _) | (_, Err(e)) => Event::Error(e.to_string()),
        },
        Command::ExchangeConsent { request_uri, given } => match app.exchange_consent(&request_uri, given).await {
            Ok(_) => Event::ExchangeConsented,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunIbdExchange { info, biosample_guid } => {
            let cfg = IbdDetectorConfig::default();
            match app.open_exchange_session(&info).await {
                Ok(session) => match app
                    .exchange_ibd_for_subject(&session, biosample_guid, &info.request_uri, None, cfg)
                    .await
                {
                    Ok(r) => Event::IbdExchangeDone {
                        biosample_guid,
                        total_shared_cm: r.summary.total_shared_cm,
                        segment_count: r.summary.segment_count,
                        relationship: format!("{:?}", r.summary.relationship),
                        agreed: r.agreed,
                    },
                    Err(e) => Event::Error(e.to_string()),
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadIbdExchanges { biosample_guid } => {
            match app.list_ibd_exchanges_for_subject(biosample_guid).await {
                Ok(rows) => Event::IbdExchanges { biosample_guid, rows },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::BackfillLabs => match app.backfill_run_labs().await {
            Ok(count) => Event::LabsResolved(count),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AuthStatus => Event::Authenticated(app.current_account()),
        Command::SyncStatus => Event::SyncOnline(app.is_online()),
        Command::PullSync => match app.pull_sync().await {
            Ok(o) => Event::PullDone {
                in_sync: o.in_sync,
                applied: o.applied,
                adopted: o.adopted,
                repushed: o.repushed,
                conflicts: o.conflicts,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::VerifySourceFiles => match app.verify_source_files().await {
            Ok(missing) => Event::SourceFilesVerified { missing },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::Login { handle } => match app.login(&handle).await {
            Ok(did) => Event::Authenticated(Some(did)),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::Logout => match app.logout().await {
            Ok(()) => Event::Authenticated(None),
            Err(e) => Event::Error(e.to_string()),
        },
        // Publishes enqueue to the durable outbox then drain — handled in the spawn loop (they emit
        // multiple events: Queued + per-row Published + SyncPending). Reaching here is a routing bug.
        Command::PublishCoverage(id) => Event::Error(format!("internal: unrouted PublishCoverage {id}")),
        Command::PublishVariants { alignment_id, .. } => {
            Event::Error(format!("internal: unrouted PublishVariants {alignment_id}"))
        }
        Command::PublishAncestry { biosample_guid } => {
            Event::Error(format!("internal: unrouted PublishAncestry {biosample_guid}"))
        }
        Command::DrainOutbox => Event::Error("internal: unrouted DrainOutbox".into()),
        Command::Export { request, path } => match app.export_content(&request).await {
            Ok(content) => match std::fs::write(&path, content) {
                Ok(()) => Event::Exported {
                    label: request.label().to_string(),
                    path,
                },
                Err(e) => Event::Error(format!("write {}: {e}", path.display())),
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ExportIbdSegments { segments, path } => {
            match std::fs::write(&path, navigator_app::export::ibd_segments_tsv(&segments)) {
                Ok(()) => Event::Exported {
                    label: "IBD segments (CSV)".into(),
                    path,
                },
                Err(e) => Event::Error(format!("write {}: {e}", path.display())),
            }
        }
        Command::LoadPcaReference => match app.ancestry_pca_reference().await {
            Ok(points) => Event::PcaReference {
                alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                points,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::SetHaploOverride {
            biosample_guid,
            dna_type,
            haplogroup,
            reason,
        } => {
            match app
                .set_manual_override(biosample_guid, dna_type, &haplogroup, reason.as_deref())
                .await
            {
                Ok(()) => Event::ReconciliationChanged {
                    biosample_guid,
                    dna_type,
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::ClearHaploOverride {
            biosample_guid,
            dna_type,
        } => match app.clear_manual_override(biosample_guid, dna_type).await {
            Ok(()) => Event::ReconciliationChanged {
                biosample_guid,
                dna_type,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAudit {
            biosample_guid,
            dna_type,
        } => match app.reconciliation_audit(biosample_guid, dna_type).await {
            Ok(entries) => Event::Audit {
                biosample_guid,
                dna_type,
                entries,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadHeteroplasmy { alignment_id } => match app.mtdna_heteroplasmy(alignment_id).await {
            Ok(sites) => Event::Heteroplasmy { alignment_id, sites },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PublishReconciliation { biosample_guid, .. } => {
            Event::Error(format!("internal: unrouted PublishReconciliation {biosample_guid:?}"))
        }
        // ---- social (Community tab) ----------------------------------------
        Command::LoadSupportThreads => match app.support_threads().await {
            Ok(items) => Event::SupportThreads(items),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSupportThread { conversation_id } => match app.support_thread(&conversation_id).await {
            Ok(messages) => Event::SupportThread {
                conversation_id,
                messages,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::OpenSupportThread { subject, body } => match app.open_support_thread(&subject, &body).await {
            Ok(conversation_id) => Event::SupportThreadPosted { conversation_id },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ReplySupportThread { conversation_id, body } => {
            match app.reply_support_thread(&conversation_id, &body).await {
                Ok(conversation_id) => Event::SupportThreadPosted { conversation_id },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadCommunityFeed => match app.community_feed().await {
            Ok(feed) => Event::CommunityFeed(feed),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PostCommunity { content, topic } => match app.post_community(&content, topic.as_deref(), None).await {
            Ok(_) => Event::CommunityPosted,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadNotifications => match app.notifications().await {
            Ok(n) => Event::Notifications {
                items: n.items,
                unread: n.unread,
            },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::MarkNotificationRead { id } => match app.mark_notification_read(id.as_deref()).await {
            Ok(_) => Event::NotificationsMarked,
            Err(e) => Event::Error(e.to_string()),
        },
    }
}

/// Resolve a reference build, emitting throttled `ReferenceProgress` events (and waking the
/// UI) as bytes arrive, then a final `ReferenceReady`/`Error`. Run from the spawn loop so it
/// can stream — `handle` returns only a single event.
async fn resolve_reference_streaming(
    app: &App,
    build: String,
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    // The progress closure must be Send (it runs in a task) — capture an owned Sender clone
    // and a label, not borrows. Throttle to ~every 25 MB so a multi-GB pull doesn't flood.
    let tx = evt_tx.clone();
    let label = build.clone();
    let mut last_sent = 0u64;
    let mut progress = move |received: u64, total: Option<u64>| {
        if received.saturating_sub(last_sent) >= 25_000_000 || total == Some(received) {
            last_sent = received;
            let _ = tx.send(Event::ReferenceProgress {
                build: label.clone(),
                received,
                total,
            });
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
    let total = 5; // unified metrics + 4 command steps

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
                            detail: format!("scanning genome — {:.0}%", within * 100.0),
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
                let _ = evt_tx.send(Event::Coverage {
                    alignment_id,
                    result: Some(cov),
                });
                let _ = evt_tx.send(Event::ReadMetrics {
                    alignment_id,
                    result: Some(rm),
                });
                let _ = evt_tx.send(Event::Sex {
                    alignment_id,
                    result: sex,
                });
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
            1 => Command::RunDenovo {
                alignment_id,
                contig: "chrM".into(),
            },
            2 => Command::AssignYHaplogroup { alignment_id },
            _ => Command::AssignMtdnaHaplogroupFromAlignment { alignment_id },
        };
        let ev = handle(app, cmd).await; // runs to completion; we may cancel before the next step
        let _ = evt_tx.send(ev);
        wake();
    }

    let _ = evt_tx.send(Event::AnalysisDone {
        cancelled: cancel.load(Ordering::Relaxed),
    });
    wake();
}

/// Deep-analyze every sample in a project one at a time, emitting `DeepAnalyzeProgress` before
/// each sample (so the bar advances sample by sample) and a final `ProjectAnalyzed`. `cancel` is
/// checked before each sample — a stop leaves the already-computed artifacts in place (the pass is
/// additive and idempotent). Each `analyze_biosample` awaits internally, so the worker runtime
/// stays free for quick UI queries between samples.
async fn deep_analyze_project_streaming(
    app: &App,
    project_id: i64,
    cancel: Arc<AtomicBool>,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    cancel.store(false, Ordering::Relaxed);
    let biosamples = match app.list_biosamples(project_id).await {
        Ok(v) => v,
        Err(e) => {
            let _ = evt_tx.send(Event::Error(e.to_string()));
            wake();
            return;
        }
    };
    let total = biosamples.len();
    let (mut samples, mut coverage_done, mut y_done, mut sex_done, mut metrics_done, mut sv_done, mut errors) =
        (0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize);

    for (i, biosample) in biosamples.iter().enumerate() {
        if cancel.load(Ordering::Relaxed) {
            break;
        }
        let _ = evt_tx.send(Event::DeepAnalyzeProgress {
            project_id,
            done: i,
            total,
            sample: biosample.donor_identifier.clone(),
            fraction: if total > 0 { i as f32 / total as f32 } else { 0.0 },
        });
        wake();
        match app.analyze_biosample(biosample).await {
            Ok(o) if o.had_alignment => {
                samples += 1;
                coverage_done += o.coverage_done as usize;
                y_done += o.y_done as usize;
                sex_done += o.sex_done as usize;
                metrics_done += o.metrics_done as usize;
                sv_done += o.sv_done as usize;
                errors += o.errors.len();
            }
            Ok(_) => {} // no BAM-bearing alignment — not counted
            Err(e) => {
                // A structural (DB/IO) failure on one sample: count it, surface it, keep going.
                errors += 1;
                let _ = evt_tx.send(Event::Error(format!("{}: {e}", biosample.donor_identifier)));
                wake();
            }
        }
    }

    let _ = evt_tx.send(Event::ProjectAnalyzed {
        project_id,
        samples,
        coverage_done,
        y_done,
        sex_done,
        metrics_done,
        sv_done,
        errors,
        cancelled: cancel.load(Ordering::Relaxed),
    });
    wake();
}

/// Report an enqueue result (`Queued`/`Error`), then drain the outbox so an online publish sends
/// immediately. `kind` is the human label for the queued feedback.
async fn publish_then_drain(
    app: &App,
    enqueue: Result<(), navigator_app::AppError>,
    kind: &str,
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    match enqueue {
        Ok(()) => {
            let _ = evt_tx.send(Event::Queued { kind: kind.to_string() });
            wake();
            emit_drain(app, evt_tx, wake).await;
        }
        Err(e) => {
            let _ = evt_tx.send(Event::Error(e.to_string()));
            wake();
        }
    }
}

/// Drain the outbox once and emit the outcome: a `Published` per sent row, the online flag, and the
/// remaining pending count.
async fn emit_drain(app: &App, evt_tx: &Sender<Event>, wake: &(dyn Fn() + Send + Sync)) {
    match app.drain_outbox().await {
        Ok(outcome) => {
            for (kind, uri) in outcome.published {
                let _ = evt_tx.send(Event::Published { kind, uri });
            }
            let _ = evt_tx.send(Event::SyncOnline(app.is_online()));
            let _ = evt_tx.send(Event::SyncPending(outcome.pending));
        }
        Err(e) => {
            let _ = evt_tx.send(Event::Error(e.to_string()));
        }
    }
    wake();
}

/// Spawn the worker thread: open the workspace at `db_path` inside the worker's runtime
/// (so the connection pool lives there), then serve commands. `wake` is called after
/// each event so the UI can `request_repaint`. Returns the command sender and event
/// receiver the UI holds.
pub fn spawn(db_path: PathBuf, wake: impl Fn() + Send + Sync + 'static) -> (UnboundedSender<Command>, Receiver<Event>) {
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

                // Background outbox drain: retry pending PDS publishes every 30s (catches
                // offline→online without a user action). Skips work when the queue is empty.
                {
                    let app = app.clone();
                    let evt_tx = evt_tx.clone();
                    let wake = wake.clone();
                    tokio::spawn(async move {
                        let mut tick = tokio::time::interval(std::time::Duration::from_secs(30));
                        tick.set_missed_tick_behavior(tokio::time::MissedTickBehavior::Skip);
                        loop {
                            tick.tick().await;
                            if app.outbox_pending_count().await.unwrap_or(0) > 0 {
                                emit_drain(&app, &evt_tx, &*wake).await;
                            }
                        }
                    });
                }

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
                            // Streams AnalysisProgress per step (+ each step's result), then AnalysisDone.
                            Command::RunFullAnalysis { alignment_id } => {
                                run_full_analysis_streaming(&app, alignment_id, cancel, &evt_tx, wake.clone()).await;
                            }
                            // Streams DeepAnalyzeProgress per sample, then a final ProjectAnalyzed.
                            Command::DeepAnalyzeProject(project_id) => {
                                deep_analyze_project_streaming(&app, project_id, cancel, &evt_tx, wake.clone()).await;
                            }
                            // Signals the in-flight full analysis / deep-analyze to stop between steps.
                            Command::CancelAnalysis => {
                                cancel.store(true, Ordering::Relaxed);
                            }
                            // Publishes enqueue durably, then drain (send-now-if-online). The drain
                            // emits Published per row + SyncPending; we emit Queued for instant feedback.
                            Command::PublishCoverage(id) => {
                                publish_then_drain(
                                    &app,
                                    app.publish_coverage(id).await,
                                    "coverage summary",
                                    &evt_tx,
                                    &*wake,
                                )
                                .await;
                            }
                            Command::PublishVariants { alignment_id, contig } => {
                                let r = app.publish_variants(alignment_id, &contig).await;
                                publish_then_drain(&app, r, &format!("{contig} variants"), &evt_tx, &*wake).await;
                            }
                            Command::PublishAncestry { biosample_guid } => {
                                let r = app.publish_ancestry(biosample_guid).await;
                                publish_then_drain(&app, r, "ancestry breakdown", &evt_tx, &*wake).await;
                            }
                            Command::PublishReconciliation {
                                biosample_guid,
                                dna_type,
                                heteroplasmy,
                                identity,
                            } => {
                                let r = app
                                    .publish_reconciliation(biosample_guid, dna_type, &heteroplasmy, identity.as_ref())
                                    .await;
                                publish_then_drain(
                                    &app,
                                    r,
                                    &format!("{} reconciliation", dna_type.as_str()),
                                    &evt_tx,
                                    &*wake,
                                )
                                .await;
                            }
                            // Periodic / on-reconnect drain of the outbox.
                            Command::DrainOutbox => {
                                emit_drain(&app, &evt_tx, &*wake).await;
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
        let created = handle(
            &app,
            Command::CreateProject(NewProject {
                name: "Trio".into(),
                description: None,
                administrator: "jk".into(),
            }),
        )
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
        match handle(
            &app,
            Command::LoadDenovo {
                alignment_id: aln.id,
                contig: "chrM".into(),
            },
        )
        .await
        {
            Event::Denovo { result, .. } => assert!(result.is_none()),
            other => panic!("expected Denovo(None), got {other:?}"),
        }
        match handle(
            &app,
            Command::RunDenovo {
                alignment_id: aln.id,
                contig: "chrM".into(),
            },
        )
        .await
        {
            Event::Denovo { contig, result, .. } => {
                assert_eq!(contig, "chrM");
                let calls = result.unwrap();
                assert_eq!(
                    calls.iter().map(|c| c.position).collect::<Vec<_>>(),
                    vec![2, 3, 4, 6, 7, 8, 10]
                );
            }
            other => panic!("expected Denovo(Some), got {other:?}"),
        }
        match handle(
            &app,
            Command::LoadDenovo {
                alignment_id: aln.id,
                contig: "chrM".into(),
            },
        )
        .await
        {
            Event::Denovo { result, .. } => assert_eq!(result.unwrap().len(), 7),
            other => panic!("expected cached Denovo, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn add_commands_create_and_signal_reload() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};

        let app = app().await;
        let pid = match handle(
            &app,
            Command::CreateProject(NewProject {
                name: "P".into(),
                description: None,
                administrator: "jk".into(),
            }),
        )
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };

        // add a sample (tagged to the project) -> BiosamplesChanged
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: Some(pid),
                donor_identifier: "HG002".into(),
                sample_accession: None,
                sex: Some("male".into()),
            }),
        )
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
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: None,
                donor_identifier: "NA12878".into(),
                sample_accession: None,
                sex: None,
            }),
        )
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
        match handle(
            &app,
            Command::AddRun(NewSequenceRun {
                biosample_guid: guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            }),
        )
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
        match handle(
            &app,
            Command::AddAlignment(NewAlignment {
                sequence_run_id: run_id,
                reference_build: "chm13v2.0".into(),
                aligner: "bwa".into(),
                variant_caller: None,
                bam_path: None,
                reference_path: None,
                content_sha256: None,
            }),
        )
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
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: None,
                donor_identifier: "draft".into(),
                sample_accession: None,
                sex: None,
            }),
        )
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
        match handle(
            &app,
            Command::UpdateBiosample {
                guid,
                donor_identifier: "HG002".into(),
                sample_accession: Some("SAMN123".into()),
                description: Some("trio son".into()),
                center_name: None,
                sex: Some("male".into()),
            },
        )
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
        match handle(
            &app,
            Command::AddRun(NewSequenceRun {
                biosample_guid: guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            }),
        )
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
        match handle(
            &app,
            Command::DeleteSequenceRun {
                id: run_id,
                biosample_guid: guid,
            },
        )
        .await
        {
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
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: None,
                donor_identifier: "spare".into(),
                sample_accession: None,
                sex: None,
            }),
        )
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
        let pid = match handle(
            &app,
            Command::CreateProject(NewProject {
                name: "P".into(),
                description: None,
                administrator: "jk".into(),
            }),
        )
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: None,
                donor_identifier: "loose".into(),
                sample_accession: None,
                sex: None,
            }),
        )
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
        match handle(
            &app,
            Command::AssignBiosampleProject {
                guid,
                project_id: Some(pid),
            },
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { samples, .. } => assert_eq!(samples.len(), 1),
            other => panic!("got {other:?}"),
        }

        // assigning to a non-existent project is refused
        match handle(
            &app,
            Command::AssignBiosampleProject {
                guid,
                project_id: Some(9999),
            },
        )
        .await
        {
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
        let pid = match handle(
            &app,
            Command::CreateProject(NewProject {
                name: "Old".into(),
                description: None,
                administrator: "jk".into(),
            }),
        )
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };

        // edit name/admin/description
        match handle(
            &app,
            Command::UpdateProject {
                id: pid,
                name: "Renamed".into(),
                description: Some("a study".into()),
                administrator: "curator".into(),
            },
        )
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
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: Some(pid),
                donor_identifier: "member".into(),
                sample_accession: None,
                sex: None,
            }),
        )
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
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: None,
                donor_identifier: "subj".into(),
                sample_accession: None,
                sex: None,
            }),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => all[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(
            &app,
            Command::AddRun(NewSequenceRun {
                biosample_guid: guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: Some(1_000),
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            }),
        )
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
        match handle(
            &app,
            Command::UpdateSequenceRun {
                id: run.id,
                biosample_guid: guid,
                platform_name: "MGI".into(),
                instrument_model: Some("DNBSEQ-T7".into()),
                test_type: "WGS".into(),
                library_layout: Some("PAIRED".into()),
                sequencing_facility: Some("Dante Labs".into()),
            },
        )
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
                assert_eq!(r.sequencing_facility.as_deref(), Some("Dante Labs")); // lab persisted
                assert_eq!(r.total_reads, Some(1_000)); // metric untouched
            }
            other => panic!("got {other:?}"),
        }

        match handle(
            &app,
            Command::AddAlignment(NewAlignment {
                sequence_run_id: run.id,
                reference_build: "grch38".into(),
                aligner: "bwa".into(),
                variant_caller: None,
                bam_path: None,
                reference_path: None,
                content_sha256: None,
            }),
        )
        .await
        {
            Event::AlignmentsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        let aln_id = match handle(&app, Command::LoadAlignments(run.id)).await {
            Event::Alignments { alignments, .. } => alignments[0].id,
            other => panic!("got {other:?}"),
        };
        match handle(
            &app,
            Command::UpdateAlignment {
                id: aln_id,
                sequence_run_id: run.id,
                reference_build: "chm13v2.0".into(),
                aligner: "minimap2".into(),
                variant_caller: Some("deepvariant".into()),
            },
        )
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

    /// The streaming deep-analyze emits one progress event per sample and a final `ProjectAnalyzed`.
    /// Samples without a BAM-bearing alignment are walked (so the bar advances) but not counted —
    /// keeping the test free of any reference/network/tree dependency.
    #[tokio::test]
    async fn deep_analyze_streams_progress_then_a_final_summary() {
        let app = app().await;
        let p = app
            .create_project(NewProject {
                name: "P".into(),
                description: None,
                administrator: "jk".into(),
            })
            .await
            .unwrap();
        app.add_biosample(Some(p.id), "S1", None, None).await.unwrap();
        app.add_biosample(Some(p.id), "S2", None, None).await.unwrap();

        let (tx, rx) = std::sync::mpsc::channel::<Event>();
        let cancel = Arc::new(AtomicBool::new(false));
        let wake: Arc<dyn Fn() + Send + Sync> = Arc::new(|| {});
        deep_analyze_project_streaming(&app, p.id, cancel, &tx, wake).await;

        let events: Vec<Event> = rx.try_iter().collect();
        let progress = events
            .iter()
            .filter(|e| matches!(e, Event::DeepAnalyzeProgress { .. }))
            .count();
        assert_eq!(progress, 2, "one progress event per sample");
        match events.last() {
            Some(Event::ProjectAnalyzed {
                project_id,
                samples,
                cancelled,
                errors,
                ..
            }) => {
                assert_eq!(*project_id, p.id);
                assert_eq!(*samples, 0, "no BAM-bearing alignments → nothing counted");
                assert_eq!(*errors, 0);
                assert!(!*cancelled);
            }
            other => panic!("expected a final ProjectAnalyzed, got {other:?}"),
        }
    }

    /// A cancel raised mid-run stops the loop before the next sample and reports `cancelled`.
    #[tokio::test]
    async fn deep_analyze_honors_a_mid_run_cancel() {
        let app = app().await;
        let p = app
            .create_project(NewProject {
                name: "P".into(),
                description: None,
                administrator: "jk".into(),
            })
            .await
            .unwrap();
        app.add_biosample(Some(p.id), "S1", None, None).await.unwrap();
        app.add_biosample(Some(p.id), "S2", None, None).await.unwrap();

        // The function clears the flag at entry, so re-arm it via a wake hook fired on the first
        // progress emission — simulating the user hitting Cancel after the first sample starts.
        let (tx, rx) = std::sync::mpsc::channel::<Event>();
        let cancel = Arc::new(AtomicBool::new(false));
        let armed = cancel.clone();
        let wake: Arc<dyn Fn() + Send + Sync> = Arc::new(move || armed.store(true, Ordering::Relaxed));
        deep_analyze_project_streaming(&app, p.id, cancel, &tx, wake).await;

        let events: Vec<Event> = rx.try_iter().collect();
        let progress = events
            .iter()
            .filter(|e| matches!(e, Event::DeepAnalyzeProgress { .. }))
            .count();
        assert_eq!(progress, 1, "cancel after S1 skips S2's progress");
        match events.last() {
            Some(Event::ProjectAnalyzed { cancelled, .. }) => assert!(*cancelled),
            other => panic!("expected a cancelled ProjectAnalyzed, got {other:?}"),
        }
    }
}
