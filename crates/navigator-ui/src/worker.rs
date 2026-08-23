//! The bridge between sync and async. egui runs the immediate-mode loop on the main thread. The
//! `App`, which is tokio and sqlx, runs on its own worker thread, with its own runtime. The UI
//! sends [`Command`] values and drains [`Event`] values on each frame. No DB call, and no domain
//! decision, happens on the UI thread (plan §5).
//!
//! Each command runs on its own task, so a long analysis never blocks a quick query. The map from
//! command to event ([`handle`]) is pure, and it has unit tests. [`spawn`] is the glue for the
//! thread, the runtime and the channels.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, Sender};
use std::sync::{Arc, Mutex};

use navigator_app::CancelToken;
use navigator_app::{
    AlignmentProbe, AnalysisStep, AncestryResult, App, AppError, ArchaicMarkerResult, ArchaicSegmentResult, AuditEntry,
    BatchImportSummary, BuildNeed, ChatTurn, Consensus, Coverage, DenovoCall, DescentReport, DmConversationSummary,
    DmMessage, DnaType, ExchangeSessionInfo, FtdnaGenealogy, FtdnaImportOptions, FtdnaImportPlan, FtdnaImportSummary,
    FtdnaResolution, HaploAssignment, HeteroplasmySite, IbdComparison, IbdDetectorConfig, IbdSuggestion,
    IdentityVerification, IncomingRequest, MatchingEntry, NarratedBrief, PaintingResult, PrivateBucket,
    ProjectBlockTree, ProjectImportSummary, ProjectOverview, ProjectSampleReport, ProjectStrChart, ReadMetrics,
    RecruitmentInvitation, RefBuildStatus, RohResult, SexInferenceResult, SignalKind, SourceType, StoredIbdExchange,
    StrConcordanceRow, SubjectAnalysisStatus, SubjectBrief, SvAnalysisResult, YMatch, YstrClustering,
};
use navigator_domain::chipprofile::ChipProfile;
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::identity::NewMdka;
use navigator_domain::mtdna::MtdnaSequence;
use navigator_domain::strprofile::StrProfile;
use navigator_domain::variants::VariantSet;
use navigator_domain::workspace::{
    Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, Project, SequenceRun,
};
use tokio::sync::mpsc::{unbounded_channel, UnboundedSender};

/// Callable-region mask choice for the private-Y bucket.
#[derive(Debug, Clone)]
pub enum YMask {
    /// Self-referential: the sample's own callable-Y BED (adapts to depth + read tech).
    SelfReferential,
    /// An external callable BED (e.g. the Poznik/1KG `b38_sites.bed`).
    Bed(PathBuf),
    /// No mask. It is noisy, because it takes every off-backbone de-novo call.
    None,
}

/// The fields to add a biosample. The app assigns its `SampleGuid`. `project_id` is optional,
/// because a biosample stands on its own, and it does not have to belong to a project.
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
    /// Survey the workspace chores (what each would do if run). Deliberately on demand: two of the
    /// three cost real work to measure, one of them a multi-MB tree fetch.
    SurveyMaintenance,
    /// Run one workspace chore. It streams `ChoreProgress`, and ends with `ChoreDone`.
    RunChore {
        chore: navigator_app::Chore,
        /// Recompute what is already cached (private-Y only; the others have nothing to force).
        force: bool,
    },
    /// Check GitHub Releases for a newer installer and notify (no auto-update).
    CheckForUpdate,
    CreateProject(NewProject),
    LoadSamples(i64),
    /// Load the coverage and haplogroup report of each sample, for a project.
    LoadProjectReport(i64),
    /// Compute and load the Y-STR overview of each member (the FTDNA-style chart) for a project.
    LoadProjectStrChart(i64),
    /// Build the cohort Y **block tree** for a project. It reads and parses a multi-MB haplotree.
    /// So the UI sends it lazily, on the first view of the Tree tab, and not on a project select.
    LoadProjectBlockTree(i64),
    /// Build (off the UI thread) the plain-language Subject Brief for a subject (Simple mode).
    LoadSubjectBrief(SampleGuid),
    /// Build (off the UI thread) a YFull-style Y/mtDNA descent report for a subject.
    LoadDescentReport {
        guid: SampleGuid,
        dna: DnaType,
    },
    /// Build a branch report for each marker, over the subtree of `node`, for a subject. It runs off
    /// the UI thread.
    LoadBranchReport {
        guid: SampleGuid,
        dna: DnaType,
        node: String,
        depth: Option<usize>,
    },
    /// Narrate the brief of a subject through the local LLM ("Polish with AI"). It falls back on any
    /// failure.
    NarrateBrief(SampleGuid),
    /// Ask the local LLM a question about a subject's results (grounded in the brief).
    AskQuestion {
        guid: SampleGuid,
        history: Vec<ChatTurn>,
        question: String,
    },
    /// Explain one result signal in plain language (the "Explain this" of a tab, M5).
    NarrateSignal {
        guid: SampleGuid,
        kind: SignalKind,
    },
    /// Deep-analyze every sample in a project, as a background job the user can cancel. It streams
    /// a `DeepAnalyzeProgress` for each sample, and gives up the thread between samples, so that
    /// the UI stays responsive. It steps over what the fast path already filled.
    /// [`Command::CancelAnalysis`] cancels it. The headless path, and the tests, still use the
    /// one-shot `App::analyze_project`.
    DeepAnalyzeProject(i64),
    /// Load every biosample (subjects list), regardless of project.
    LoadAllBiosamples,
    /// Load donor-level Y/mt terminal haplogroups for every subject (fills the list columns).
    LoadHaploSummary,
    /// Load the analysis status of each subject (Pending or Complete), for the Status column of the
    /// subjects list.
    LoadSubjectStatus,
    AddBiosample(NewBiosample),
    /// Batch-import a NAS project directory: scan → Project, Biosample, Run, Alignment.
    /// `reference` is optional. With `None`, the gateway resolves each build from the cache, and
    /// reports `ReferenceNeeded` when a download is necessary. `Some` pins a FASTA.
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
    /// Commit a reviewed FTDNA import plan, with the resolution the admin chose for each kit.
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
    /// Rank every other subject in the workspace against this one by Y relatedness (gap §2). It is
    /// one against all over the workspace, or over one project when `project_id` has a value. It
    /// reads cached profiles.
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
    /// Assign a Y haplogroup from the imported BISDNA or Y-SNP-panel calls of the subject, with no
    /// alignment. It records a donor call.
    AssignYBisdna {
        biosample_guid: SampleGuid,
    },
    /// Assign an mtDNA haplogroup directly from an alignment's chrM (records a donor call).
    AssignMtdnaHaplogroupFromAlignment {
        alignment_id: i64,
    },
    /// Estimate autosomal ancestry from the CONSENSUS of the subject, with no BAM walk. This is the
    /// default path.
    EstimateAncestryFromConsensus {
        biosample_guid: SampleGuid,
    },
    /// Estimate **deep (ancient) ancestry** through qpAdm. It genotypes the best CHM13 alignment of
    /// the subject at the full 1240k. It is heavy, at about 1 to 2 minutes, and it runs only on an
    /// explicit request.
    EstimateDeepAncestry {
        biosample_guid: SampleGuid,
    },
    /// Paint local ancestry from the subject's CONSENSUS (no BAM walk).
    PaintAncestryFromConsensus {
        biosample_guid: SampleGuid,
    },
    /// Load the cached chromosome painting, when it is current for the consensus signature. The
    /// cost is low.
    LoadPainting {
        biosample_guid: SampleGuid,
    },
    /// Detect runs of homozygosity (ROH) from the subject's CONSENSUS (no BAM walk).
    ComputeRohFromConsensus {
        biosample_guid: SampleGuid,
    },
    /// Load the cached ROH result, when it is current for the consensus signature. The cost is low.
    LoadRoh {
        biosample_guid: SampleGuid,
    },
    /// Count archaic (Neanderthal / Denisovan) marker copies from the subject's CONSENSUS.
    ComputeArchaicFromConsensus {
        biosample_guid: SampleGuid,
    },
    /// Load the cached archaic marker count, when it is current for the consensus signature. The
    /// cost is low.
    LoadArchaic {
        biosample_guid: SampleGuid,
    },
    /// Call Tier B archaic SEGMENTS from cached genome-wide diploid calls.
    CallArchaicSegments {
        biosample_guid: SampleGuid,
    },
    /// Load the cached Tier B segment result. The cost is low.
    LoadArchaicSegments {
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
    /// Unified import: more than one file, or folder, or both. It walks a folder for data files. It
    /// detects the type of each one and sends it to the right place, then returns one
    /// [`Event::DataBatchImported`] summary.
    AddDataBatch {
        biosample_guid: SampleGuid,
        paths: Vec<PathBuf>,
    },
    /// A convenience for a first run: make a subject, and import `paths` into it, in one step. The
    /// empty state of Simple mode can then go from nothing to a full brief with one file pick. It
    /// returns [`Event::SubjectCreatedAndImported`], with the new guid and the import summary.
    CreateSubjectAndImport {
        donor_identifier: String,
        sex: Option<String>,
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
    /// Load the multi-source Y-variant profile of the subject, which is the concordance over all Y
    /// sources.
    ///
    /// Load the persisted Y-profile snapshot. The cost is low, and nothing genotypes.
    LoadYProfile {
        biosample_guid: SampleGuid,
    },
    /// Compute the Y-profile again from all sources, and persist the snapshot. The cost is high,
    /// because it genotypes again.
    BuildYProfile {
        biosample_guid: SampleGuid,
    },
    /// Resolve catalogued Y-SNP names at the given positions (annotates the Y-SNP tables).
    LoadYSnpNames {
        biosample_guid: SampleGuid,
        positions: Vec<i64>,
    },
    /// Load the persisted mtDNA consensus-profile snapshot. The cost is low, and nothing
    /// genotypes.
    LoadMtProfile {
        biosample_guid: SampleGuid,
    },
    /// Compute the mtDNA consensus profile again from all sources, and persist it. The cost is
    /// high, because it places again.
    BuildMtProfile {
        biosample_guid: SampleGuid,
    },
    /// Load the persisted autosomal consensus-profile snapshot. The cost is low, and nothing
    /// genotypes.
    LoadAutosomalProfile {
        biosample_guid: SampleGuid,
    },
    /// Compute the autosomal consensus again from all sources, and persist it. The cost is high,
    /// because it genotypes at the panel.
    BuildAutosomalProfile {
        biosample_guid: SampleGuid,
    },
    /// Probe a BAM/CRAM header for build/aligner/platform/test-type (to auto-fill the form).
    ProbeAlignment {
        path: PathBuf,
    },
    LoadCoverage(i64),
    /// Cached coverage for more than one alignment at a time (the alignment rows of Data Sources).
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
    LoadAllAlignments,
    /// How many alignments of one project a realignment would act on.
    ///
    /// This asks the app, and the UI does not count them. The card that shows this number used to
    /// filter `all_alignments`, which is the whole workspace. So a project whose alignments were
    /// every one already on the target build still read that 35 of them "in this project" could
    /// map again. The batch it would have started covered the project, and was correct. Only the
    /// number was wrong, and that is the worse way round.
    LoadRealignableInProject {
        project_id: i64,
        target_build: String,
    },
    /// Compare two samples over the IBD panel that works with a chip. Each sample is a WGS alignment
    /// or an imported chip. This is the volume path, for chip↔chip and chip↔WGS.
    CompareIbdSources {
        a: navigator_app::IbdSource,
        b: navigator_app::IbdSource,
    },
    /// Compare two SUBJECTS over their autosomal consensuses (the subject-level IBD path).
    CompareIbdConsensus {
        a: SampleGuid,
        b: SampleGuid,
    },
    /// Check whether two SUBJECTS are the same individual, over their pooled autosomal consensus,
    /// with no panel.
    VerifyIdentityConsensus {
        a: SampleGuid,
        b: SampleGuid,
    },
    /// Federated IBD step 1: fetch the AppView's pseudonymous match suggestions for the
    /// signed-in account (registers the device key on first use).
    LoadIbdSuggestions,
    /// Ask for an introduction to a candidate, and record the conversation in the matching ledger.
    RequestIntroduction {
        suggestion: IbdSuggestion,
        biosample_guid: Option<SampleGuid>,
    },
    /// Tell the AppView to drop a candidate from its suggestions.
    DismissCandidate {
        suggested_sample_guid: String,
    },
    /// Adopt a local did:key identity that certifies itself. It is the desktop bootstrap, with no
    /// PDS and no OAuth.
    UseLocalIdentity,
    /// Reconcile the matching ledger against the broker (inbound requests + consent-ready sessions)
    /// and return every conversation.
    RefreshMatching,
    /// Consent to an inbound exchange request, or decline it, and record the decision durably.
    MatchingConsent {
        request_uri: String,
        given: bool,
        biosample_guid: Option<SampleGuid>,
    },
    /// Drop a conversation from the local ledger (forget, not cancel).
    ForgetMatchingRequest {
        request_uri: String,
    },
    /// Run a full IBD exchange for a subject, over a session that has consent: handshake → dosage
    /// exchange → signed attestations → persist. It takes a long time, and the peer must be
    /// online.
    RunIbdExchange {
        info: ExchangeSessionInfo,
        biosample_guid: SampleGuid,
    },
    /// Load the subject's persisted IBD exchange results.
    LoadIbdExchanges {
        biosample_guid: SampleGuid,
    },
    /// Peer DMs (social 3a): open a DM request to a partner DID.
    DmInitiate {
        partner_did: String,
    },
    /// Poll the DM inbox: the inbound DM requests that wait for consent, and the sessions that have
    /// consent and no connection yet.
    LoadDmInbox,
    /// Consent to (or decline) an inbound DM request.
    DmConsent {
        request_uri: String,
        given: bool,
    },
    /// Connect a consent-ready DM session (one-time handshake; persists the session key).
    DmConnect {
        info: ExchangeSessionInfo,
    },
    /// Load the persisted DM conversation list.
    LoadDmConversations,
    /// Load (and mark read) one conversation's transcript.
    LoadDmMessages {
        session_id: String,
    },
    /// Encrypt + relay a message on a conversation.
    DmSend {
        session_id: String,
        text: String,
    },
    /// Pull, decrypt and persist any message that waits on a conversation.
    DmSync {
        session_id: String,
    },
    /// Recruitment 3c: poll the signed-in account's open recruitment invitations.
    LoadRecruitmentInvitations,
    /// Accept (`true`) or decline (`false`) a recruitment invitation.
    RespondRecruitment {
        campaign_id: i64,
        accept: bool,
    },
    /// Resolve the sequencing lab for a run that has an inferred instrument id and no facility. It
    /// goes through the instrument→lab map of the AppView, which is best-effort and cached. The UI
    /// sends it on startup, and after an import.
    BackfillLabs,
    /// Report who signed in. It has no side effect, and the UI sends it on startup.
    AuthStatus,
    /// Report the current online/offline state (no side effects).
    SyncStatus,
    /// Sign in to a PDS through OAuth, which opens a browser. `handle` is a handle or a DID.
    Login {
        handle: String,
    },
    Logout,
    PublishCoverage(i64),
    PublishVariants {
        alignment_id: i64,
        contig: String,
    },
    /// Publish the consensus ancestry breakdown of the subject, one record for each method, to the
    /// PDS that signed in.
    PublishAncestry {
        biosample_guid: SampleGuid,
    },
    /// Publish the subject anchor to the PDS that signed in. That anchor is the anonymized biosample
    /// summary, plus its sequence runs. Every derived record, for coverage or ancestry, ties back to
    /// it.
    PublishBiosample {
        biosample_guid: SampleGuid,
    },
    /// Try to push the ready outbox rows now. It also runs at intervals, and after a publish.
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
    /// Run the full analysis pipeline of one alignment: coverage → sex → metrics → variant calling
    /// → Y haplogroup → ancestry. It streams an `AnalysisProgress` for each step. It also forwards
    /// the result event of each step, so that the detail tabs fill in while it runs. Structural
    /// variants are **not** part of this. See [`Command::RunSv`], which the Sources tab dispatches
    /// on request.
    RunFullAnalysis {
        alignment_id: i64,
    },
    /// Run the full analysis on the representative alignment of a subject. It resolves that
    /// alignment from the guid (see `default_alignment_for_subject`). The Simple "My DNA" view uses
    /// it. A casual user can then analyze with no visit to the Advanced sources table to find an
    /// alignment id.
    AnalyzeSubject {
        biosample_guid: SampleGuid,
    },
    /// Request cancellation of the in-flight full analysis (checked between steps).
    CancelAnalysis,
    /// Realign an off-build alignment onto another reference (design/realignment-module.md). It
    /// takes hours. It streams a `RealignProgress` for each stage, then `RealignDone`.
    StartRealign {
        alignment_id: i64,
        target_build: String,
    },
    /// Stop a realignment in progress. It shares the cancellation registry with the analysis
    /// pipeline, because only one long job runs at a time, by design.
    CancelRealign,
    /// Realign every eligible alignment in a project, one after another.
    ///
    /// One after another, and not in parallel. Each job already fills the cores of the machine, and
    /// it wants ~12 GB. Two at a time would be slower than two in turn, and they could exhaust the
    /// memory. A cancel stops the current job, and abandons the rest of the queue.
    StartProjectRealign {
        project_id: i64,
        target_build: String,
    },
    /// Update a subject's editable fields. Empty optional values clear the column.
    UpdateBiosample {
        guid: SampleGuid,
        donor_identifier: String,
        sample_accession: Option<String>,
        description: Option<String>,
        center_name: Option<String>,
        sex: Option<String>,
    },
    /// Attach a vendor id (kit number) to a subject from the subject editor. `(source, external_id)`
    /// is the dedup anchor; the app layer refuses an id already bound to another subject.
    AddExternalId {
        guid: SampleGuid,
        source: String,
        external_id: String,
    },
    /// Detach a vendor id (by row id) from a subject. `guid` is the owner, so the UI can refresh
    /// that subject's genealogy card.
    DeleteExternalId {
        guid: SampleGuid,
        id: i64,
    },
    /// Insert or update a subject's MDKA (most distant known ancestor) for one lineage.
    UpsertMdka {
        guid: SampleGuid,
        mdka: NewMdka,
    },
    /// Remove a subject's MDKA for one lineage (`Y` | `Mt` | `Auto`).
    DeleteMdka {
        guid: SampleGuid,
        lineage: String,
    },
    /// Delete a subject. Refused by the app layer if it still has dependent data.
    DeleteBiosample(SampleGuid),
    /// Clear all the sequencing data of a subject, and all the analysis data, whether derived or
    /// imported. The subject itself stays.
    ClearBiosampleData(SampleGuid),
    /// Reset the haplogroup placement of the subject, and nothing else. It is the cleanup for a
    /// stale lineage, and the other data stays.
    ClearHaplogroupData(SampleGuid),
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
    /// Load the reference-genome settings of each build, and the cache status, for the Settings
    /// dialog.
    LoadReferenceSettings,
    /// Force a refresh of the cached haplotrees. It clears the session memo and the on-disk cache,
    /// so that a corrected AppView tree arrives with no restart of the app. A profile then
    /// interprets against it again on the next load.
    RefreshTrees,
    /// Health-check a local-LLM server at `base_url` (Settings "Test connection"): lists its models.
    TestLlmConnection {
        base_url: String,
    },
    /// Persist **all** the reference-source overrides at one time (the "References" table of
    /// Settings). One command → one atomic write. One command for each row raced the config file
    /// into corruption (#26).
    SetReferenceOverrides(Vec<navigator_app::ReferenceOverrideInput>),
    /// Hash a cached reference again, and compare it with its integrity sidecar. The integrity-check
    /// button of Settings sends this.
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
    // ---- social (the Community tab, over the signed AppView Edge API) -------
    /// List the support threads of the account that signed in (team↔tester).
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
    /// Post to the community feed, with an optional topic tag. When `publish_pds` has a value, the
    /// post *also* goes to the PDS that signed in, as a federated `feed.post` record (roadmap 3b).
    PostCommunity {
        content: String,
        topic: Option<String>,
        publish_pds: bool,
    },
    /// Fetch notifications + unread count (also drives the app-bar bell).
    LoadNotifications,
    /// Mark one notification read (`Some`) or all (`None`).
    MarkNotificationRead {
        id: Option<String>,
    },
}

/// A result/notification from the worker to the UI.
#[derive(Debug, Clone)]
pub enum Event {
    /// Nothing to report, for example a cache load that missed. The UI ignores it.
    Noop,
    Overview(Vec<ProjectOverview>),
    /// Ancestry/IBD asset presence + integrity (the "data sources" transparency line).
    AssetStatus(Vec<navigator_app::AssetStatus>),
    /// Haplotrees were force-refreshed (N cached files cleared); the UI re-loads open profiles.
    TreesRefreshed(usize),
    ProjectCreated(Project),
    /// Something changed or deleted a project, so load the overview again.
    ProjectsChanged,
    /// A batch project-directory import completed.
    ProjectImported(ProjectImportSummary),
    /// A dry-run FTDNA import plan, ready for the review modal.
    FtdnaPlan(FtdnaImportPlan),
    /// The result of a commit of an FTDNA import.
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
    /// A reference build resolved. It is now in the cache, and it has an index.
    ReferenceReady {
        build: String,
        path: PathBuf,
    },
    /// A build of a BAM or CRAM coordinate index (`.bai`/`.crai`) is in progress, so that a region
    /// query works. `total` is the compressed file size for a BAM, which gives a byte fraction. It
    /// is `None` for a CRAM, which gives a spinner with no fraction.
    IndexProgress {
        file: String,
        done: u64,
        total: Option<u64>,
    },
    /// A coordinate index build ended, or nothing needed one. `built` is the path the code wrote,
    /// when there is one.
    IndexReady {
        built: Option<PathBuf>,
    },
    /// A newer installer exists on GitHub Releases. The app tells the user, and it never updates
    /// itself.
    UpdateAvailable(Box<navigator_app::UpdateInfo>),
    /// The installer-update check ran, and the app is already current. This also covers a check
    /// that did not run.
    UpToDate,
    /// The coverage and haplogroup report of each sample, for a project.
    ProjectReport {
        project_id: i64,
        rows: Vec<ProjectSampleReport>,
    },
    /// The Y-STR overview of each member (the FTDNA-style chart) for a project, computed first.
    ProjectStrChart {
        project_id: i64,
        chart: ProjectStrChart,
    },
    /// The cohort Y block tree for a project. `tree` is `None` when the project has no member at
    /// all. That differs from a tree where placement reached nobody, which comes back with everybody
    /// `unplaced`.
    ProjectBlockTree {
        project_id: i64,
        tree: Box<Option<ProjectBlockTree>>,
    },
    /// The plain-language Subject Brief for a subject (Simple mode).
    SubjectBrief {
        guid: SampleGuid,
        brief: Box<SubjectBrief>,
    },
    /// A YFull-style descent report for a subject's Y or mtDNA lineage (`None` = placed but empty /
    /// not placed; `Err` = a load failure surfaced to the status line).
    DescentReportLoaded {
        guid: SampleGuid,
        dna: DnaType,
        result: Result<Option<DescentReport>, String>,
    },
    /// A branch report for each marker, over the Y or mtDNA subtree of a subject. `None` means there
    /// is no alignment. `Err` is a failure of the load or the lookup, for example a node nothing
    /// found, and it goes to the status line.
    BranchReportLoaded {
        guid: SampleGuid,
        dna: DnaType,
        result: Result<Option<navigator_app::BranchReport>, String>,
    },
    /// A slice of narration text from the stream, while the model writes it. It is a live preview,
    /// and the final BriefNarration is the authoritative one.
    BriefNarrationChunk {
        guid: SampleGuid,
        text: String,
    },
    /// AI-assisted narration of a subject's brief (or a plain-language reason it is unavailable).
    BriefNarration {
        guid: SampleGuid,
        result: Result<NarratedBrief, String>,
    },
    /// A slice of a chat answer from the stream, while the model writes it. It is a live preview.
    ChatAnswerChunk {
        guid: SampleGuid,
        text: String,
    },
    /// Answer to an "ask my results" question (or a plain-language reason it is unavailable).
    ChatAnswer {
        guid: SampleGuid,
        result: Result<String, String>,
    },
    /// A slice of an "Explain this" narration for one signal, from the stream. It is a live
    /// preview.
    SignalNarrationChunk {
        guid: SampleGuid,
        kind: SignalKind,
        text: String,
    },
    /// AI-assisted explanation of one result signal (or a plain-language reason it is unavailable).
    SignalNarration {
        guid: SampleGuid,
        kind: SignalKind,
        result: Result<NarratedBrief, String>,
    },
    /// An analyze pass over the whole project ended, with coverage and Y for each sample.
    /// `cancelled` is true when something stopped a streaming deep-analyze early, and the counts
    /// then cover only what ended before the stop.
    ProjectAnalyzed {
        project_id: i64,
        samples: usize,
        coverage_done: usize,
        y_done: usize,
        sex_done: usize,
        metrics_done: usize,
        errors: usize,
        cancelled: bool,
    },
    /// The progress of a streaming deep-analyze pass, one sample at a time: `done` of `total`
    /// samples so far. `sample` is the donor id under analysis now, and `fraction` drives the bar
    /// (0..1).
    DeepAnalyzeProgress {
        project_id: i64,
        done: usize,
        total: usize,
        sample: String,
        fraction: f32,
    },
    /// The workspace-chore survey: what each chore would do if run now.
    MaintenanceSurvey(Vec<navigator_app::ChoreSurvey>),
    /// The progress of a chore in progress, one item at a time. `label` is the subject it works on
    /// now.
    ChoreProgress {
        chore: navigator_app::Chore,
        done: usize,
        total: usize,
        label: String,
        fraction: f32,
    },
    /// A chore ended, or a cancel stopped it after `outcome.done` items.
    ChoreDone {
        chore: navigator_app::Chore,
        outcome: navigator_app::ChoreOutcome,
    },
    /// The progress of a streaming project-directory import, one sample at a time: `done` of
    /// `total` samples written. `sample` is the sample id the import writes now, and `fraction`
    /// drives the bar.
    ImportProgress {
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
    /// The terminal Y and mt haplogroups of each subject, for the subjects list
    /// (`guid → (Y, mt)`).
    HaploSummary(std::collections::HashMap<SampleGuid, (Option<String>, Option<String>)>),
    /// The analysis status of each subject (Pending or Complete), for the Status column of the
    /// subjects list.
    SubjectStatus(std::collections::HashMap<SampleGuid, SubjectAnalysisStatus>),
    /// Something added or changed a biosample. Load the subjects list again, and any open project
    /// view.
    BiosamplesChanged,
    Runs {
        biosample_guid: SampleGuid,
        runs: Vec<SequenceRun>,
    },
    RunsChanged(SampleGuid),
    /// Something cleared the analysis data of a subject. The UI then loads that subject again in
    /// full, and the list columns too.
    BiosampleDataCleared(SampleGuid),
    /// A subject's haplogroup placement was reset (the UI reloads that subject; other data stays).
    HaplogroupDataReset(SampleGuid),
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
    /// A batch import ended. The summary lists, for each file, whether the import took it or
    /// stepped over it. The UI shows that in a modal, and loads the data sections of the subject
    /// again.
    DataBatchImported {
        biosample_guid: SampleGuid,
        summary: BatchImportSummary,
    },
    /// A first-run "Import DNA" created a subject and imported into it. The UI selects the new
    /// subject (so its brief appears) and shows the import summary.
    SubjectCreatedAndImported {
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
    /// The local-ancestry painting of each chromosome (the "DNA painting"): the segments of each
    /// side, and the side labels.
    AncestryPainting {
        alignment_id: i64,
        result: PaintingResult,
    },
    /// Runs-of-homozygosity result for a subject (segments + F_ROH summary). `None` on a cache-miss load.
    RohResultReady {
        biosample_guid: SampleGuid,
        result: Option<Box<RohResult>>,
    },
    /// Archaic (Tier A) marker count for a subject. `None` on a cache-miss load.
    ArchaicResultReady {
        biosample_guid: SampleGuid,
        result: Option<Box<ArchaicMarkerResult>>,
    },
    /// Tier B archaic segment result. `None` on a cache-miss load.
    ArchaicSegmentsReady {
        biosample_guid: SampleGuid,
        result: Option<Box<ArchaicSegmentResult>>,
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
    /// Result of an on-demand deep (ancient) ancestry estimate. `result` is `None` when the deep model
    /// does not apply to this subject (non-European / no CHM13 alignment / rejected fit).
    DeepAncestryEstimated {
        biosample_guid: SampleGuid,
        // Boxed: AncestryResult is large; keep the Event enum small.
        result: Option<Box<AncestryResult>>,
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
    /// `ancient` is absent when the deep sources can't express this sample's ancestry.
    ConsensusAncestryDetail {
        biosample_guid: SampleGuid,
        // Boxed: AncestryResult is large, and two of them would bloat the Event enum's size.
        fine: Option<Box<navigator_app::AncestryResult>>,
        ancient: Option<Box<navigator_app::AncestryResult>>,
    },
    /// Header-probe result for the add-alignment form (build/aligner/platform/test-type).
    AlignmentProbe(AlignmentProbe),
    Coverage {
        alignment_id: i64,
        result: Option<Coverage>,
    },
    /// Cached coverage for more than one alignment (the Data Sources rows):
    /// `(alignment_id, result)`.
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
    /// The progress of the full-analysis pipeline: `step` of `total` begins, and both count from 1.
    /// It carries a `label`, a `detail`, and the bar `fraction` (0..1).
    AnalysisProgress {
        step: usize,
        total: usize,
        label: String,
        detail: String,
        fraction: f32,
    },
    /// The full-analysis pipeline ended, or a cancel stopped it.
    AnalysisDone {
        cancelled: bool,
    },
    /// A realignment stage began.
    RealignProgress {
        alignment_id: i64,
        /// The subject this job belongs to, so a card can tell whether the run is *theirs*.
        /// `None` only if the lookup failed, in which case no card claims it.
        biosample_guid: Option<SampleGuid>,
        step: usize,
        total: usize,
        label: String,
        detail: String,
    },
    /// A realignment over the whole project ended. It is separate from `RealignDone`, which fires
    /// for each sample, because a batch has its own outcome. That outcome is how many of the queue
    /// ended, and whether the rest went. An empty queue reports `queued: 0`, and never nothing.
    RealignBatchDone {
        queued: usize,
        completed: usize,
        cancelled: bool,
    },
    /// A realignment ended, a cancel stopped it, or it failed. `new_alignment_id` has a value only
    /// on success. The insert of that row comes last, so an absence means nothing registered.
    RealignDone {
        alignment_id: i64,
        biosample_guid: Option<SampleGuid>,
        new_alignment_id: Option<i64>,
        cancelled: bool,
        summary: String,
    },
    AllAlignments(Vec<Alignment>),
    /// The alignments in `project_id` that a realignment to the target build would act on. It
    /// carries the project id. A reply that comes back after the user moves on then goes away, and
    /// it does not appear against the wrong project.
    RealignableInProject {
        project_id: i64,
        ids: Vec<i64>,
    },
    Ibd(IbdComparison),
    /// Federated IBD match suggestions from the AppView (may be empty in a single-user dev AppView).
    IbdSuggestions(Vec<IbdSuggestion>),
    /// The matching ledger: every conversation, with its result. The refresh emits it, and so does
    /// every mutation, so the panel never has to poll the broker again to see its own action.
    Matching(Vec<MatchingEntry>),
    /// Something dismissed a candidate, so the UI drops its row.
    CandidateDismissed {
        suggested_sample_guid: String,
    },
    /// A DM request went out to a partner DID. The UI refreshes the inbox.
    DmInitiated,
    /// The DM inbox: inbound DM requests + consent-ready sessions to connect.
    DmInbox {
        incoming: Vec<IncomingRequest>,
        ready: Vec<ExchangeSessionInfo>,
    },
    /// The store took a DM consent. The UI refreshes the inbox.
    DmConsented,
    /// A DM session connected, and the key is on disk. The UI refreshes the conversation list.
    DmConnected,
    /// The persisted DM conversation list.
    DmConversations(Vec<DmConversationSummary>),
    /// One conversation's transcript.
    DmMessages {
        session_id: String,
        rows: Vec<DmMessage>,
    },
    /// A DM went out. The UI loads the open transcript again.
    DmSent {
        session_id: String,
    },
    /// A DM sync finished; `new_count` newly-stored messages (the UI reloads if > 0).
    DmSynced {
        session_id: String,
        new_count: usize,
    },
    /// The signed-in account's open recruitment invitations.
    RecruitmentInvitations(Vec<RecruitmentInvitation>),
    /// The store took a response to a recruitment invitation. The UI refreshes the invitations and
    /// the notifications.
    RecruitmentResponded,
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
    /// A record went out. `kind` is a human label, and `uri` is the `at://` URI.
    Published {
        kind: String,
        uri: String,
    },
    /// A publish went into the durable outbox. It sends now if the app is online, and on the next
    /// connection if it is not.
    Queued {
        kind: String,
    },
    /// The outbox rows that still wait for a push to succeed (the "N pending" indicator).
    SyncPending(i64),
    /// A result went to `path`. `label` is the human kind, for example "coverage (TSV)".
    Exported {
        label: String,
        path: PathBuf,
    },
    /// Whether the last PDS write reached the server (offline indicator).
    SyncOnline(bool),
    /// A PULL reconcile ended (gap §5-p2). It carries the tally of each action.
    PullDone {
        in_sync: usize,
        applied: usize,
        adopted: usize,
        repushed: usize,
        conflicts: usize,
    },
    /// A second check on whether the source files are reachable ended. A `missing` file moved, or
    /// something deleted it.
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
    /// Local-LLM "Test connection" result: the models the server reports, or a plain-language error.
    LlmConnection(Result<Vec<String>, String>),
    /// The reference-genome settings of each build, and the cache status, for the Settings dialog.
    ReferenceSettings(Vec<RefBuildStatus>),
    /// The store took a reference override. The UI can load the settings rows again.
    ReferenceSettingsChanged,
    /// The result of a reference integrity check: one short human-readable status for each build.
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
    /// Somebody opened a support thread, or replied to one. Load the list again, and the open
    /// thread.
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
    /// Something marked the notifications read. Load them again.
    NotificationsMarked,
    Error(String),
    /// A command failed, **and** a preflight at file level found a concrete cause. `message` is the
    /// original error, and the status bar still shows it. `report` is the diagnosis to paste, and it
    /// names the file at fault.
    ///
    /// It is separate from [`Event::Error`], so that the UI can offer the report. The UI does not
    /// have to guess from a string whether an error has one. It appears only when the preflight
    /// failed a check. A tree download, or a network error, must not raise a file report.
    Diagnosed {
        message: String,
        report: String,
    },
    /// A run stopped because the user cancelled it.
    ///
    /// It differs from `Error`, because this is not a failure. It also differs from `Noop`, which
    /// would leave the control that asked for it in a spin for ever. A standalone SV or de-novo run
    /// has no `AnalysisDone` to clear its flag.
    Cancelled,
}

/// How a cancellation reads after an event flattens it to a string.
const CANCELLED_MESSAGE: &str = "cancelled";

/// The token of whichever cancellable run is in progress, so that `CancelAnalysis` can reach it.
///
/// It replaces one shared `AtomicBool` that each run reset to `false` at its own start. Every
/// command goes through `tokio::spawn`, so that reset raced the click. A cancel that arrived between
/// the spawn and the reset went away with no message. A second run that started at the same time
/// also wiped the cancel the first one still needed. A [`CancelToken`] comes into existence one time
/// for each run, and nothing ever un-cancels it, so there is no window where a cancel can go.
///
/// The generation counter stops a run that ends from a clear of the registration of a *newer* run.
/// That is the same stale-write fault in different clothes.
#[derive(Clone, Default)]
struct CancelRegistry {
    current: Arc<Mutex<Option<(u64, CancelToken)>>>,
    next_gen: Arc<AtomicU64>,
}

impl CancelRegistry {
    /// Register a fresh token for a run that starts. Returns its generation and the token.
    fn begin(&self) -> (u64, CancelToken) {
        let gen = self.next_gen.fetch_add(1, Ordering::Relaxed);
        let token = CancelToken::new();
        *self.current.lock().unwrap() = Some((gen, token.clone()));
        (gen, token)
    }

    /// Retire the registration of this run, but only when no newer run already replaced it.
    fn end(&self, gen: u64) {
        let mut slot = self.current.lock().unwrap();
        if slot.as_ref().is_some_and(|(g, _)| *g == gen) {
            *slot = None;
        }
    }

    /// Cancel whatever runs now. It does nothing when nothing runs, and that is what makes a stray
    /// click harmless, and not something that poisons the next run.
    fn cancel_current(&self) {
        if let Some((_, token)) = self.current.lock().unwrap().as_ref() {
            token.cancel();
        }
    }
}

/// Settle a finished alignment command into the event the UI must see.
///
/// Two things have to happen between a command that fails and a user who reads about it. Both are
/// there to keep the wrong thing off the screen.
///
/// 1. This drops a **cancellation**. It travels as an error, so that it can unwind the walk from
///    deep inside a walker. But the user asked for it, and "Error: cancelled" would report their
///    own click back to them as a failure. The `AnalysisDone { cancelled }` of the run already says
///    it.
/// 2. A **genuine** failure gets a diagnosis at file level.
///
/// The errors this upgrades are the ones that name a path, and not the *right* path. The reader
/// helpers report whatever path the call that failed received. So a bad index, a reference nothing
/// can read, and a file behind a privacy denial all come out as `io error on <the alignment>`. A
/// run of [`App::diagnose_alignment`] probes each of those files on its own, and says which one it
/// is.
///
/// An error with no cause at file level passes through untouched. If every preflight check passes,
/// the failure is truly somewhere else, in a tree fetch, a liftover, or the appview. A clean bill of
/// health would then be worse than nothing.
async fn settle_alignment_command(app: &App, alignment_id: i64, event: Event) -> Event {
    let Event::Error(message) = event else {
        return event;
    };
    // A walk the user cancelled is not a file problem, and not a failure. A diagnosis of it would
    // be a slow lie, and a report of it would contradict the action of the user.
    if message == CANCELLED_MESSAGE {
        return Event::Cancelled;
    }
    match app.diagnose_alignment(alignment_id).await {
        Ok(report) if report.failed() => Event::Diagnosed {
            message,
            report: report.to_string(),
        },
        _ => Event::Error(message),
    }
}

/// Do one command against the app, and map its success or failure to an [`Event`].
///
/// Read the genealogy bundle of a subject again (the vendor ids, the FTDNA member, and the MDKA),
/// and wrap it as an [`Event::Genealogy`]. This is the refresh that follows any genealogy mutation,
/// so that the detail card shows the new state, with no separate "changed" round-trip.
async fn reload_genealogy(app: &App, guid: SampleGuid) -> Event {
    ev(app.subject_genealogy(guid).await, |data| Event::Genealogy {
        guid,
        data,
    })
}

/// Map an app call that can fail to an [`Event`]. `ok` names the success event, and **any** error
/// becomes `Event::Error`, with the `Display` text of that error. This is the one place that writes
/// the policy down.
///
/// Almost every arm of [`handle`] has this shape, so one definition here keeps ~120 call sites from
/// each writing it out. `ok` is usually the event constructor alone
/// (`ev(app.refresh_trees().await, Event::TreesRefreshed)`), and a closure covers a struct variant,
/// and a case with extra fields. It is generic over the error type, so that a store error, an
/// analysis error and an io error all take the same route. It is `FnOnce`, so that the closure can
/// move a captured value into the event.
///
/// An arm that needs more than this keeps its explicit `match`. That covers a three-way split of
/// `Ok(Some)` and `Ok(None)`, an `.await` on the success path, and a message that is not
/// `Display`.
fn ev<T, E: std::fmt::Display>(result: Result<T, E>, ok: impl FnOnce(T) -> Event) -> Event {
    match result {
        Ok(value) => ok(value),
        Err(e) => Event::Error(e.to_string()),
    }
}

pub async fn handle(app: &App, cmd: Command, cancel: &CancelToken) -> Event {
    match cmd {
        Command::LoadOverview => ev(app.project_overview().await, Event::Overview),
        Command::LoadAssetStatus => Event::AssetStatus(navigator_app::ancestry_asset_status()),
        Command::CheckForUpdate => match app.check_for_update().await {
            Ok(Some(info)) => Event::UpdateAvailable(Box::new(info)),
            Ok(None) => Event::UpToDate,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RefreshTrees => ev(app.refresh_trees().await, Event::TreesRefreshed),
        Command::CreateProject(new) => ev(app.create_project(new).await, Event::ProjectCreated),
        // ImportProjectDir streams ImportProgress from the spawn loop, so a path that arrives here
        // is a fault.
        Command::ImportProjectDir { .. } => Event::Error("internal: unrouted ImportProjectDir".into()),
        Command::PlanFtdnaImport {
            project_id,
            project_name,
            member,
            paternal,
            maternal,
            ystr,
        } => ev(
            app.plan_ftdna_import(
                project_id,
                project_name,
                member,
                paternal,
                maternal,
                ystr,
                FtdnaImportOptions::default(),
            )
            .await,
            Event::FtdnaPlan,
        ),
        Command::CommitFtdnaImport { plan, resolutions } => {
            ev(app.commit_ftdna_import(&plan, &resolutions).await, Event::FtdnaImported)
        }
        Command::LoadGenealogy(guid) => ev(app.subject_genealogy(guid).await, |data| Event::Genealogy {
            guid,
            data,
        }),
        Command::ClusterProject(project_id) => ev(app.cluster_project_ystr(project_id).await, |clustering| {
            Event::ProjectClustering { project_id, clustering }
        }),
        // The spawn loop controls ResolveReference, because it streams progress events. A path that
        // arrives here would be a fault in the routes.
        Command::ResolveReference { build } => Event::Error(format!("internal: unrouted ResolveReference {build}")),
        Command::LoadSamples(project_id) => ev(app.list_biosamples(project_id).await, |samples| Event::Samples {
            project_id,
            samples,
        }),
        Command::LoadProjectReport(project_id) => ev(app.project_report(project_id).await, |rows| {
            Event::ProjectReport { project_id, rows }
        }),
        Command::LoadProjectStrChart(project_id) => ev(app.project_str_chart(project_id).await, |chart| {
            Event::ProjectStrChart { project_id, chart }
        }),
        Command::LoadProjectBlockTree(project_id) => ev(app.project_block_tree(project_id, DnaType::Y).await, |tree| {
            Event::ProjectBlockTree {
                project_id,
                tree: Box::new(tree),
            }
        }),
        Command::LoadSubjectBrief(guid) => ev(app.subject_brief(guid).await, |brief| Event::SubjectBrief {
            guid,
            brief: Box::new(brief),
        }),
        Command::LoadDescentReport { guid, dna } => Event::DescentReportLoaded {
            guid,
            dna,
            result: app.descent_report(guid, dna).await.map_err(|e| e.to_string()),
        },
        Command::LoadBranchReport { guid, dna, node, depth } => Event::BranchReportLoaded {
            guid,
            dna,
            result: app
                .branch_report_for_subject(guid, dna, &node, depth)
                .await
                .map_err(|e| e.to_string()),
        },
        // NarrateBrief and AskQuestion stream from the spawn loop, so a path that arrives here is
        // a fault.
        Command::NarrateBrief(guid) => Event::Error(format!("internal: unrouted NarrateBrief {guid}")),
        Command::AskQuestion { guid, .. } => Event::Error(format!("internal: unrouted AskQuestion {guid}")),
        Command::NarrateSignal { guid, .. } => Event::Error(format!("internal: unrouted NarrateSignal {guid}")),
        // DeepAnalyzeProject streams DeepAnalyzeProgress from the spawn loop, so a path that
        // arrives here is a fault.
        Command::DeepAnalyzeProject(project_id) => {
            Event::Error(format!("internal: unrouted DeepAnalyzeProject {project_id}"))
        }
        Command::SurveyMaintenance => ev(app.maintenance_survey().await, Event::MaintenanceSurvey),
        // RunChore streams ChoreProgress from the spawn loop, so a path that arrives here is a
        // fault.
        Command::RunChore { chore, .. } => Event::Error(format!("internal: unrouted RunChore {}", chore.key())),
        Command::LoadSubjectStatus => ev(app.subject_analysis_status().await, Event::SubjectStatus),
        Command::LoadHaploSummary => ev(app.haplogroup_terminals().await, Event::HaploSummary),
        Command::LoadAllBiosamples => ev(app.list_all_biosamples().await, Event::AllBiosamples),
        Command::AddBiosample(b) => ev(
            app.add_biosample(b.project_id, b.donor_identifier, b.sample_accession, b.sex)
                .await,
            |_| Event::BiosamplesChanged,
        ),
        Command::UpdateBiosample {
            guid,
            donor_identifier,
            sample_accession,
            description,
            center_name,
            sex,
        } => ev(
            app.update_biosample(guid, donor_identifier, sample_accession, description, center_name, sex)
                .await,
            |_| Event::BiosamplesChanged,
        ),
        Command::AddExternalId {
            guid,
            source,
            external_id,
        } => match app.add_external_id(guid, &source, &external_id).await {
            Ok(_) => reload_genealogy(app, guid).await,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteExternalId { guid, id } => match app.delete_external_id(id).await {
            Ok(()) => reload_genealogy(app, guid).await,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::UpsertMdka { guid, mdka } => match app.upsert_mdka(guid, mdka).await {
            Ok(()) => reload_genealogy(app, guid).await,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteMdka { guid, lineage } => match app.delete_mdka(guid, &lineage).await {
            Ok(()) => reload_genealogy(app, guid).await,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DeleteBiosample(guid) => ev(app.delete_biosample(guid).await, |_| Event::BiosamplesChanged),
        Command::ClearBiosampleData(guid) => ev(app.clear_biosample_data(guid).await, |_| {
            Event::BiosampleDataCleared(guid)
        }),
        Command::ClearHaplogroupData(guid) => ev(app.clear_haplogroup_data(guid).await, |_| {
            Event::HaplogroupDataReset(guid)
        }),
        Command::DeleteSequenceRun { id, biosample_guid } => ev(app.delete_sequence_run(id).await, |_| {
            Event::RunsChanged(biosample_guid)
        }),
        Command::MergeSequenceRuns {
            biosample_guid,
            primary,
            secondary,
        } => ev(
            app.merge_sequence_runs(biosample_guid, primary, secondary).await,
            |_| Event::RunsChanged(biosample_guid),
        ),
        Command::DeleteAlignment { id, sequence_run_id } => ev(app.delete_alignment(id).await, |_| {
            Event::AlignmentsChanged(sequence_run_id)
        }),
        Command::DeleteStrProfile { id, biosample_guid } => ev(app.delete_str_profile(id).await, |_| {
            Event::StrProfilesChanged(biosample_guid)
        }),
        Command::LoadReferenceSettings => Event::ReferenceSettings(app.reference_settings()),
        Command::TestLlmConnection { base_url } => {
            Event::LlmConnection(app.llm_models_at(&base_url).await.map_err(|e| e.to_string()))
        }
        Command::SetReferenceOverrides(rows) => {
            ev(app.set_reference_overrides(&rows), |_| Event::ReferenceSettingsChanged)
        }
        Command::VerifyReference { build } => ev(app.verify_reference(&build).await, |outcome| {
            let status = match outcome {
                navigator_app::VerifyOutcome::Verified => "✔ verified".to_string(),
                navigator_app::VerifyOutcome::Mismatch { .. } => "✖ mismatch (corrupted?)".to_string(),
                navigator_app::VerifyOutcome::NoSidecar => "• no checksum on record".to_string(),
                navigator_app::VerifyOutcome::NotCached => "not cached".to_string(),
            };
            Event::ReferenceVerified { build, status }
        }),
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
                    ev(
                        app.lift_vcf(&src, &target, in_vcf, out_vcf.clone(), opts, &mut |_, _| {})
                            .await,
                        |s| Event::VcfLifted {
                            summary: format!(
                                "Lifted {}/{} variants ({} unmapped, {} ref-mismatch) › {}",
                                s.lifted,
                                s.total,
                                s.unmapped,
                                s.ref_mismatch,
                                out_vcf.display()
                            ),
                        },
                    )
                }
            }
        }
        Command::DeleteVariantSet { id, biosample_guid } => ev(app.delete_variant_set(id).await, |_| {
            Event::VariantSetsChanged(biosample_guid)
        }),
        Command::DeleteChipProfile { id, biosample_guid } => ev(app.delete_chip_profile(id).await, |_| {
            Event::ChipProfilesChanged(biosample_guid)
        }),
        Command::DeleteMtdnaSequence { id, biosample_guid } => ev(app.delete_mtdna_sequence(id).await, |_| {
            Event::MtdnaChanged(biosample_guid)
        }),
        Command::AssignBiosampleProject { guid, project_id } => {
            ev(app.add_biosample_to_project(guid, project_id).await, |_| {
                Event::BiosamplesChanged
            })
        }
        Command::UpdateProject {
            id,
            name,
            description,
            administrator,
        } => ev(app.update_project(id, name, description, administrator).await, |_| {
            Event::ProjectsChanged
        }),
        Command::DeleteProject(id) => ev(app.delete_project(id).await, |_| Event::ProjectsChanged),
        Command::UpdateSequenceRun {
            id,
            biosample_guid,
            platform_name,
            instrument_model,
            test_type,
            library_layout,
            sequencing_facility,
        } => ev(
            app.update_sequence_run(
                id,
                platform_name,
                instrument_model,
                test_type,
                library_layout,
                sequencing_facility,
            )
            .await,
            |_| Event::RunsChanged(biosample_guid),
        ),
        Command::UpdateAlignment {
            id,
            sequence_run_id,
            reference_build,
            aligner,
            variant_caller,
        } => ev(
            app.update_alignment(id, reference_build, aligner, variant_caller).await,
            |_| Event::AlignmentsChanged(sequence_run_id),
        ),
        Command::LoadRuns(biosample_guid) => ev(app.list_sequence_runs(biosample_guid).await, |runs| Event::Runs {
            biosample_guid,
            runs,
        }),
        Command::AddRun(new) => ev(app.record_sequence_run(new).await, |run| {
            Event::RunsChanged(run.biosample_guid)
        }),
        Command::LoadConsensus(guid) => {
            let y = app.haplogroup_consensus(guid, DnaType::Y).await.unwrap_or(None);
            let mt = app.haplogroup_consensus(guid, DnaType::Mt).await.unwrap_or(None);
            Event::Consensus {
                biosample_guid: guid,
                y,
                mt,
            }
        }
        Command::StrConcordance { biosample_guid } => ev(
            app.str_concordance_for_subject(biosample_guid).await,
            |(alignment_id, rows)| Event::StrConcordance {
                biosample_guid,
                alignment_id,
                rows,
            },
        ),
        Command::YMatches {
            biosample_guid,
            project_id,
        } => ev(app.y_matches(biosample_guid, project_id).await, |matches| {
            Event::YMatches {
                biosample_guid,
                matches,
            }
        }),
        Command::LoadStrProfiles(guid) => ev(app.list_str_profiles(guid).await, |profiles| Event::StrProfiles {
            biosample_guid: guid,
            profiles,
        }),
        Command::ImportStrProfile {
            biosample_guid,
            panel_name,
            provider,
            source,
            path,
        } => ev(
            app.import_str_profile_from_csv(biosample_guid, &panel_name, provider, source, &path)
                .await,
            |_| Event::StrProfilesChanged(biosample_guid),
        ),
        Command::LoadVariantSets(guid) => ev(app.list_variant_sets(guid).await, |sets| Event::VariantSets {
            biosample_guid: guid,
            sets,
        }),
        Command::ImportVariants {
            biosample_guid,
            path,
            source_type,
        } => ev(
            app.import_variants_from_file(biosample_guid, &path, source_type).await,
            |_| Event::VariantSetsChanged(biosample_guid),
        ),
        Command::AddVariants {
            biosample_guid,
            source_label,
            source_type,
            text,
        } => ev(
            app.add_variants(biosample_guid, &source_label, source_type, &text)
                .await,
            |_| Event::VariantSetsChanged(biosample_guid),
        ),
        Command::LoadChipProfiles(guid) => ev(app.list_chip_profiles(guid).await, |profiles| Event::ChipProfiles {
            biosample_guid: guid,
            profiles,
        }),
        Command::ImportChipProfile {
            biosample_guid,
            provider,
            path,
        } => ev(
            app.import_chip_profile_from_csv(biosample_guid, provider, None, &path)
                .await,
            |_| Event::ChipProfilesChanged(biosample_guid),
        ),
        Command::LoadMtdna(guid) => ev(app.list_mtdna_sequences(guid).await, |sequences| {
            Event::MtdnaSequences {
                biosample_guid: guid,
                sequences,
            }
        }),
        Command::ImportMtdna { biosample_guid, path } => {
            ev(app.import_mtdna_from_fasta(biosample_guid, &path).await, |_| {
                Event::MtdnaChanged(biosample_guid)
            })
        }
        Command::LoadMtdnaVariants { mtdna_id } => ev(app.mtdna_variants(mtdna_id).await, |variants| {
            Event::MtdnaVariants { mtdna_id, variants }
        }),
        Command::AssignMtdnaHaplogroup { mtdna_id } => ev(app.assign_mtdna_haplogroup(mtdna_id).await, |assignment| {
            Event::Haplogroup { mtdna_id, assignment }
        }),
        Command::AssignYBisdna { biosample_guid } => {
            ev(app.assign_y_bisdna(biosample_guid, None).await, |assignment| {
                Event::YBisdnaHaplogroup {
                    biosample_guid,
                    assignment,
                }
            })
        }
        Command::YHaploReport { alignment_id } => {
            ev(app.y_haplogroup_report(alignment_id).await, |(assignment, lineage)| {
                Event::YHaploReport {
                    alignment_id,
                    assignment,
                    lineage,
                }
            })
        }
        Command::AssignYHaplogroup { alignment_id } => ev(app.assign_y_haplogroup(alignment_id).await, |assignment| {
            Event::YHaplogroup {
                alignment_id,
                assignment,
            }
        }),
        Command::AssignMtdnaHaplogroupFromAlignment { alignment_id } => ev(
            app.assign_mtdna_haplogroup_from_alignment(alignment_id).await,
            |assignment| Event::MtHaplogroup {
                alignment_id,
                assignment,
            },
        ),
        Command::EstimateAncestryFromConsensus { biosample_guid } => {
            // Estimate from the pooled consensus, then surface it as the donor-level result.
            ev(app.estimate_ancestry_from_consensus(biosample_guid).await, |result| {
                Event::DonorAncestry {
                    alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                    result,
                }
            })
        }
        Command::EstimateDeepAncestry { biosample_guid } => {
            // Heavy: genotypes the best CHM13 alignment at ~1.15M sites, then fits qpAdm f4. Persists
            // the ANCIENT_ADMIXTURE result (or nothing, when the model does not apply).
            ev(app.estimate_deep_ancestry(biosample_guid).await, |result| {
                Event::DeepAncestryEstimated {
                    biosample_guid,
                    result: result.map(Box::new),
                }
            })
        }
        Command::PaintAncestryFromConsensus { biosample_guid } => {
            // A painting from the consensus needs no genotyping pass. It is fast, and it streams no
            // progress.
            ev(
                app.paint_local_ancestry_from_consensus(biosample_guid).await,
                |result| Event::AncestryPainting {
                    alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                    result,
                },
            )
        }
        Command::LoadPainting { biosample_guid } => ev(app.cached_painting(biosample_guid).await, |result| {
            Event::AncestryPainting {
                alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
                result: result.unwrap_or_default(),
            }
        }),
        Command::ComputeRohFromConsensus { biosample_guid } => {
            // ROH from the consensus needs no genotyping pass. It is fast, and it streams no
            // progress.
            ev(app.compute_roh_from_consensus(biosample_guid).await, |result| {
                Event::RohResultReady {
                    biosample_guid,
                    result: Some(Box::new(result)),
                }
            })
        }
        Command::LoadRoh { biosample_guid } => {
            ev(app.cached_roh(biosample_guid).await, |result| Event::RohResultReady {
                biosample_guid,
                result: result.map(Box::new),
            })
        }
        Command::ComputeArchaicFromConsensus { biosample_guid } => {
            // A pure read over the cached consensus and the marker panel, with no genotyping pass.
            ev(app.estimate_archaic_from_consensus(biosample_guid).await, |result| {
                Event::ArchaicResultReady {
                    biosample_guid,
                    result: Some(Box::new(result)),
                }
            })
        }
        Command::LoadArchaic { biosample_guid } => ev(app.cached_archaic(biosample_guid).await, |result| {
            Event::ArchaicResultReady {
                biosample_guid,
                result: result.map(Box::new),
            }
        }),
        Command::CallArchaicSegments { biosample_guid } => {
            ev(app.call_archaic_segments_for_subject(biosample_guid).await, |result| {
                Event::ArchaicSegmentsReady {
                    biosample_guid,
                    result: Some(Box::new(result)),
                }
            })
        }
        Command::LoadArchaicSegments { biosample_guid } => {
            ev(app.cached_archaic_segments(biosample_guid).await, |result| {
                Event::ArchaicSegmentsReady {
                    biosample_guid,
                    result: result.map(Box::new),
                }
            })
        }
        Command::LoadConsensusAncestryDetail { biosample_guid } => {
            let fine = app
                .consensus_ancestry(biosample_guid, "FINE_ADMIXTURE")
                .await
                .unwrap_or(None)
                .map(Box::new);
            // This reads ANCIENT_ADMIXTURE only. The retired PCA_PROJECTION_GMM and G25_NMONTE
            // rows can still sit in a database from before the rebuild, and nothing must show them
            // again.
            let ancient = if navigator_app::ANCIENT_ANCESTRY_ENABLED {
                app.consensus_ancestry(biosample_guid, navigator_app::ANCIENT_ADMIXTURE)
                    .await
                    .unwrap_or(None)
                    .map(Box::new)
            } else {
                None
            };
            Event::ConsensusAncestryDetail {
                biosample_guid,
                fine,
                ancient,
            }
        }
        // RunFullAnalysis streams AnalysisProgress from the spawn loop, and CancelAnalysis sets the
        // shared cancel flag there. A path that arrives here would be a fault in the routes.
        Command::RunFullAnalysis { alignment_id } => {
            Event::Error(format!("internal: unrouted RunFullAnalysis {alignment_id}"))
        }
        Command::AnalyzeSubject { biosample_guid } => {
            Event::Error(format!("internal: unrouted AnalyzeSubject {biosample_guid}"))
        }
        Command::CancelAnalysis => Event::Error("internal: unrouted CancelAnalysis".into()),
        Command::StartRealign { alignment_id, .. } => {
            Event::Error(format!("internal: unrouted StartRealign {alignment_id}"))
        }
        Command::CancelRealign => Event::Error("internal: unrouted CancelRealign".into()),
        Command::StartProjectRealign { project_id, .. } => {
            Event::Error(format!("internal: unrouted StartProjectRealign {project_id}"))
        }
        Command::FindPrivateY { alignment_id, mask } => {
            let result = match mask {
                YMask::SelfReferential => app.private_y_variants_self_masked(alignment_id).await,
                YMask::Bed(p) => app.private_y_variants(alignment_id, Some(&p)).await,
                YMask::None => app.private_y_variants(alignment_id, None).await,
            };
            ev(result, |bucket| Event::PrivateY { alignment_id, bucket })
        }
        Command::LoadPrivateY { alignment_id } => match app.cached_private_y(alignment_id).await {
            Ok(Some(bucket)) => Event::PrivateY { alignment_id, bucket },
            Ok(None) => Event::Noop,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddDataBatch { biosample_guid, paths } => {
            ev(app.add_data_batch(biosample_guid, paths, |_, _| {}).await, |summary| {
                Event::DataBatchImported {
                    biosample_guid,
                    summary,
                }
            })
        }
        Command::CreateSubjectAndImport {
            donor_identifier,
            sex,
            paths,
        } => match app.add_biosample(None, donor_identifier, None, sex).await {
            Ok(bio) => ev(app.add_data_batch(bio.guid, paths, |_, _| {}).await, |summary| {
                Event::SubjectCreatedAndImported {
                    biosample_guid: bio.guid,
                    summary,
                }
            }),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAlignments(sequence_run_id) => ev(app.list_alignments(sequence_run_id).await, |alignments| {
            Event::Alignments {
                sequence_run_id,
                alignments,
            }
        }),
        Command::AddAlignment(new) => ev(app.record_alignment(new).await, |a| {
            Event::AlignmentsChanged(a.sequence_run_id)
        }),
        Command::ProbeAlignment { path } => ev(app.probe_alignment(path).await, Event::AlignmentProbe),
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
        Command::LoadYProfile { biosample_guid } => {
            ev(app.cached_y_profile(biosample_guid).await, |profile| Event::YProfile {
                biosample_guid,
                profile,
            })
        }
        Command::BuildYProfile { biosample_guid } => {
            ev(app.build_y_profile(biosample_guid).await, |profile| Event::YProfile {
                biosample_guid,
                profile: Some(profile),
            })
        }
        Command::LoadYSnpNames {
            biosample_guid,
            positions,
        } => ev(app.y_snp_names_at(biosample_guid, &positions).await, |names| {
            Event::YSnpNames { names }
        }),
        Command::LoadMtProfile { biosample_guid } => ev(app.cached_mt_profile(biosample_guid).await, |profile| {
            Event::MtProfile {
                biosample_guid,
                profile,
            }
        }),
        Command::BuildMtProfile { biosample_guid } => {
            ev(app.build_mt_profile(biosample_guid).await, |profile| Event::MtProfile {
                biosample_guid,
                profile: Some(profile),
            })
        }
        Command::LoadAutosomalProfile { biosample_guid } => {
            ev(app.cached_autosomal_profile(biosample_guid).await, |profile| {
                Event::AutosomalProfile {
                    biosample_guid,
                    profile,
                }
            })
        }
        Command::BuildAutosomalProfile { biosample_guid } => {
            ev(app.build_autosomal_profile(biosample_guid).await, |profile| {
                Event::AutosomalProfile {
                    biosample_guid,
                    profile: Some(profile),
                }
            })
        }
        Command::LoadCoverage(alignment_id) => ev(app.cached_coverage(alignment_id).await, |result| Event::Coverage {
            alignment_id,
            result,
        }),
        Command::LoadCoverageBulk(ids) => ev(app.cached_coverage_bulk(&ids).await, Event::CoverageBulk),
        Command::LoadGenomeRegions { alignment_id, build } => {
            ev(app.genome_regions(&build).await, |regions| Event::GenomeRegions {
                alignment_id,
                regions: Some(regions),
            })
        }
        Command::RunCoverage(alignment_id) => ev(app.run_coverage_for_alignment(alignment_id).await, |result| {
            Event::Coverage {
                alignment_id,
                result: Some(result),
            }
        }),
        Command::LoadSex(alignment_id) => ev(app.cached_sex(alignment_id).await, |result| Event::Sex {
            alignment_id,
            result,
        }),
        Command::RunSex(alignment_id) => ev(app.run_sex(alignment_id).await, |result| Event::Sex {
            alignment_id,
            result: Some(result),
        }),
        Command::LoadReadMetrics(alignment_id) => ev(app.cached_read_metrics(alignment_id).await, |result| {
            Event::ReadMetrics { alignment_id, result }
        }),
        Command::RunReadMetrics(alignment_id) => {
            ev(app.run_read_metrics(alignment_id).await, |result| Event::ReadMetrics {
                alignment_id,
                result: Some(result),
            })
        }
        Command::LoadSv(alignment_id) => ev(app.cached_sv(alignment_id).await, |result| Event::Sv {
            alignment_id,
            result,
        }),
        Command::RunSv(alignment_id) => ev(app.run_sv(alignment_id, cancel.clone()).await, |result| Event::Sv {
            alignment_id,
            result: Some(result),
        }),
        Command::LoadDenovo { alignment_id, contig } => {
            ev(app.cached_denovo(alignment_id, &contig).await, |result| Event::Denovo {
                alignment_id,
                contig,
                result,
            })
        }
        Command::RunDenovo { alignment_id, contig } => ev(
            app.run_denovo_for_alignment(alignment_id, contig.clone()).await,
            |result| Event::Denovo {
                alignment_id,
                contig,
                result: Some(result),
            },
        ),
        Command::LoadAllAlignments => ev(app.list_all_alignments().await, Event::AllAlignments),
        Command::LoadRealignableInProject {
            project_id,
            target_build,
        } => match app.realignable_in_project(project_id, &target_build).await {
            Ok(ids) => Event::RealignableInProject { project_id, ids },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::CompareIbdConsensus { a, b } => ev(
            app.compare_ibd_consensus(a, b, IbdDetectorConfig::default()).await,
            Event::Ibd,
        ),
        Command::CompareIbdSources { a, b } => ev(
            app.compare_ibd_sources(a, b, IbdDetectorConfig::default()).await,
            Event::Ibd,
        ),
        Command::VerifyIdentityConsensus { a, b } => ev(app.verify_identity_consensus(a, b).await, Event::Identity),
        Command::LoadIbdSuggestions => ev(app.ibd_suggestions().await, Event::IbdSuggestions),
        Command::RequestIntroduction {
            suggestion,
            biosample_guid,
        } => match app.request_introduction(&suggestion, biosample_guid).await {
            Ok(_) => ev(app.matching_entries().await, Event::Matching),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::DismissCandidate { suggested_sample_guid } => {
            ev(app.ibd_dismiss(&suggested_sample_guid).await, |_| {
                Event::CandidateDismissed { suggested_sample_guid }
            })
        }
        Command::UseLocalIdentity => ev(app.use_local_identity(), |did| Event::Authenticated(Some(did))),
        Command::RefreshMatching => ev(app.refresh_matching().await, Event::Matching),
        Command::MatchingConsent {
            request_uri,
            given,
            biosample_guid,
        } => match app.matching_consent(&request_uri, given, biosample_guid).await {
            Ok(_) => ev(app.matching_entries().await, Event::Matching),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ForgetMatchingRequest { request_uri } => match app.forget_matching_request(&request_uri).await {
            Ok(()) => ev(app.matching_entries().await, Event::Matching),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunIbdExchange { info, biosample_guid } => {
            let cfg = IbdDetectorConfig::default();
            // The store records a failure on the conversation, and it is not only a toast that
            // goes away. If not, the request sits at READY, and the user can not tell that anything
            // tried it.
            let outcome = match app.open_exchange_session(&info).await {
                Ok(session) => {
                    app.exchange_ibd_for_subject(&session, biosample_guid, &info.request_uri, None, cfg)
                        .await
                }
                Err(e) => Err(e),
            };
            match outcome {
                Ok(r) => Event::IbdExchangeDone {
                    biosample_guid,
                    total_shared_cm: r.summary.total_shared_cm,
                    segment_count: r.summary.segment_count,
                    relationship: format!("{:?}", r.summary.relationship),
                    agreed: r.agreed,
                },
                Err(e) => {
                    let msg = e.to_string();
                    let _ = app.record_matching_failure(&info.request_uri, &msg).await;
                    Event::Error(msg)
                }
            }
        }
        Command::LoadIbdExchanges { biosample_guid } => {
            ev(app.list_ibd_exchanges_for_subject(biosample_guid).await, |rows| {
                Event::IbdExchanges { biosample_guid, rows }
            })
        }
        Command::DmInitiate { partner_did } => ev(app.dm_initiate(&partner_did).await, |_| Event::DmInitiated),
        Command::LoadDmInbox => match (app.dm_incoming().await, app.dm_ready().await) {
            (Ok(incoming), Ok(ready)) => Event::DmInbox { incoming, ready },
            (Err(e), _) | (_, Err(e)) => Event::Error(e.to_string()),
        },
        Command::DmConsent { request_uri, given } => {
            ev(app.dm_consent(&request_uri, given).await, |_| Event::DmConsented)
        }
        Command::DmConnect { info } => ev(app.dm_connect(&info).await, |_| Event::DmConnected),
        Command::LoadDmConversations => ev(app.dm_conversations().await, Event::DmConversations),
        Command::LoadDmMessages { session_id } => ev(app.dm_messages(&session_id).await, |rows| Event::DmMessages {
            session_id,
            rows,
        }),
        Command::DmSend { session_id, text } => {
            ev(app.dm_send(&session_id, &text).await, |_| Event::DmSent { session_id })
        }
        Command::DmSync { session_id } => ev(app.dm_sync(&session_id).await, |new_count| Event::DmSynced {
            session_id,
            new_count,
        }),
        Command::LoadRecruitmentInvitations => ev(app.recruitment_invitations().await, Event::RecruitmentInvitations),
        Command::RespondRecruitment { campaign_id, accept } => {
            ev(app.recruitment_respond(campaign_id, accept).await, |_| {
                Event::RecruitmentResponded
            })
        }
        Command::BackfillLabs => ev(app.backfill_run_labs().await, Event::LabsResolved),
        Command::AuthStatus => Event::Authenticated(app.current_account()),
        Command::SyncStatus => Event::SyncOnline(app.is_online()),
        Command::PullSync => ev(app.pull_sync().await, |o| Event::PullDone {
            in_sync: o.in_sync,
            applied: o.applied,
            adopted: o.adopted,
            repushed: o.repushed,
            conflicts: o.conflicts,
        }),
        Command::VerifySourceFiles => ev(app.verify_source_files().await, |missing| Event::SourceFilesVerified {
            missing,
        }),
        Command::Login { handle } => ev(app.login(&handle).await, |did| Event::Authenticated(Some(did))),
        Command::Logout => ev(app.logout().await, |_| Event::Authenticated(None)),
        // A publish goes into the durable outbox, and then drains. The spawn loop controls that,
        // because it emits more than one event: Queued, then a Published for each row, then a
        // pending-count update. A path that arrives here is a fault in the routes.
        Command::PublishCoverage(id) => Event::Error(format!("internal: unrouted PublishCoverage {id}")),
        Command::PublishVariants { alignment_id, .. } => {
            Event::Error(format!("internal: unrouted PublishVariants {alignment_id}"))
        }
        Command::PublishAncestry { biosample_guid } => {
            Event::Error(format!("internal: unrouted PublishAncestry {biosample_guid}"))
        }
        Command::PublishBiosample { biosample_guid } => {
            Event::Error(format!("internal: unrouted PublishBiosample {biosample_guid}"))
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
        Command::LoadPcaReference => ev(app.ancestry_pca_reference().await, |points| Event::PcaReference {
            alignment_id: navigator_app::CONSENSUS_SOURCE_ID,
            points,
        }),
        Command::SetHaploOverride {
            biosample_guid,
            dna_type,
            haplogroup,
            reason,
        } => ev(
            app.set_manual_override(biosample_guid, dna_type, &haplogroup, reason.as_deref())
                .await,
            |_| Event::ReconciliationChanged {
                biosample_guid,
                dna_type,
            },
        ),
        Command::ClearHaploOverride {
            biosample_guid,
            dna_type,
        } => ev(app.clear_manual_override(biosample_guid, dna_type).await, |_| {
            Event::ReconciliationChanged {
                biosample_guid,
                dna_type,
            }
        }),
        Command::LoadAudit {
            biosample_guid,
            dna_type,
        } => ev(app.reconciliation_audit(biosample_guid, dna_type).await, |entries| {
            Event::Audit {
                biosample_guid,
                dna_type,
                entries,
            }
        }),
        Command::LoadHeteroplasmy { alignment_id } => ev(app.mtdna_heteroplasmy(alignment_id).await, |sites| {
            Event::Heteroplasmy { alignment_id, sites }
        }),
        Command::PublishReconciliation { biosample_guid, .. } => {
            Event::Error(format!("internal: unrouted PublishReconciliation {biosample_guid:?}"))
        }
        // ---- social (Community tab) ----------------------------------------
        Command::LoadSupportThreads => ev(app.support_threads().await, Event::SupportThreads),
        Command::LoadSupportThread { conversation_id } => ev(app.support_thread(&conversation_id).await, |messages| {
            Event::SupportThread {
                conversation_id,
                messages,
            }
        }),
        Command::OpenSupportThread { subject, body } => {
            ev(app.open_support_thread(&subject, &body).await, |conversation_id| {
                Event::SupportThreadPosted { conversation_id }
            })
        }
        Command::ReplySupportThread { conversation_id, body } => ev(
            app.reply_support_thread(&conversation_id, &body).await,
            |conversation_id| Event::SupportThreadPosted { conversation_id },
        ),
        Command::LoadCommunityFeed => ev(app.community_feed().await, Event::CommunityFeed),
        Command::PostCommunity {
            content,
            topic,
            publish_pds,
        } => match app.post_community(&content, topic.as_deref(), None).await {
            // The native post landed. When the user opted into federation, also publish the
            // durable `feed.post` record. A publish that fails reaches the screen, and the post
            // itself stays.
            Ok(_) => {
                if publish_pds {
                    ev(app.publish_feed_post(&content, topic.as_deref()).await, |_| {
                        Event::CommunityPosted
                    })
                } else {
                    Event::CommunityPosted
                }
            }
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadNotifications => ev(app.notifications().await, |n| Event::Notifications {
            items: n.items,
            unread: n.unread,
        }),
        Command::MarkNotificationRead { id } => ev(app.mark_notification_read(id.as_deref()).await, |_| {
            Event::NotificationsMarked
        }),
    }
}

/// Resolve a reference build. It emits a throttled `ReferenceProgress` event as bytes arrive, and
/// wakes the UI, then a final `ReferenceReady` or `Error`. It runs from the spawn loop, so that it
/// can stream, because `handle` returns one event only.
async fn resolve_reference_streaming(
    app: &App,
    build: String,
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    // The progress closure must be Send, because it runs in a task. So capture an owned Sender
    // clone and a label, and not a borrow. Throttle it to about every 25 MB, so that a multi-GB
    // download does not flood the channel.
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
    let event = ev(app.resolve_reference(&build, &mut progress).await, |path| {
        Event::ReferenceReady { build, path }
    });
    let _ = evt_tx.send(event);
    wake();
}

/// Resolve each build that the cache does not hold, with a visible `ReferenceProgress` bar, through
/// [`resolve_reference_streaming`]. A build the cache holds goes by with no message.
///
/// Any BAM or CRAM analysis needs the reference FASTA. The code reads it on demand: the cache
/// first, and a multi-GB download if that misses. Without this, the download runs deep inside a
/// pure
/// `App` method, behind a callback that does nothing. The UI then shows nothing, and a first import
/// looks as though it "did not register". Call this from the worker after an import, and before an
/// analysis that needs a reference, so that the download is visible.
async fn ensure_references_streaming(
    app: &App,
    builds: &[String],
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    for build in builds {
        if app.reference_cached(build) {
            continue;
        }
        resolve_reference_streaming(app, build.clone(), evt_tx, wake).await;
    }
}

/// Build the coordinate index (`.bai`/`.crai`) of the alignment when it is missing. It emits a
/// throttled `IndexProgress`, then a final `IndexReady`. An analysis that runs queries needs the
/// index to seek by region. Without it, that analysis errors, or it falls back to a scan of the
/// whole file. A build up front, with a visible bar, keeps a file that just came in from a look of
/// being stuck on its first analysis. A file that already has an index returns at once with
/// `built: None`, and no progress noise.
async fn ensure_index_streaming(app: &App, alignment_id: i64, evt_tx: &Sender<Event>, wake: &(dyn Fn() + Send + Sync)) {
    // Progress runs on a thread that blocks, so the callback must be Send. Capture owned clones,
    // and not borrows. The analysis layer already throttles it, to about every 32 MB, so forward
    // each tick.
    let tx = evt_tx.clone();
    let label = app
        .reference_build_of_alignment(alignment_id)
        .await
        .ok()
        .flatten()
        .map(|b| format!("alignment {alignment_id} ({b})"))
        .unwrap_or_else(|| format!("alignment {alignment_id}"));
    let progress = move |done: u64, total: Option<u64>| {
        let _ = tx.send(Event::IndexProgress {
            file: label.clone(),
            done,
            total,
        });
    };
    // This pre-flight is the first thing on the path from a button to touch the file. So it is
    // where an alignment or index that nothing can read appears first, and where the raw message
    // says the least. An *index build* that fails reports the path of the alignment. Diagnose it
    // before the report.
    let event = match app.ensure_alignment_index(alignment_id, progress).await {
        Ok(built) => Event::IndexReady { built },
        Err(e) => settle_alignment_command(app, alignment_id, Event::Error(e.to_string())).await,
    };
    // Wake once at the start (the first progress tick may lag on a small file) and after completion.
    wake();
    let _ = evt_tx.send(event);
    wake();
}

/// Build the missing coordinate index for each of a subject's alignments (see
/// [`ensure_index_streaming`]). Used after an import and before a subject-level analysis so every
/// query path is ready.
async fn ensure_indexes_for_subject_streaming(
    app: &App,
    biosample_guid: SampleGuid,
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    if let Ok(ids) = app.alignment_ids_for_subject(biosample_guid).await {
        for id in ids {
            ensure_index_streaming(app, id, evt_tx, wake).await;
        }
    }
}

/// Run a realignment. It emits a `RealignProgress` as each stage begins, and `RealignDone` at the
/// end.
///
/// The reference must resolve before the job starts. A realignment to a build whose FASTA is not in
/// the cache would otherwise stop at the index stage, with no message, while gigabytes download.
async fn run_realign_streaming(
    app: &App,
    alignment_id: i64,
    target_build: String,
    cancel: CancelToken,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    // Whose job this is. The code resolves it one time at the start, and not on each event.
    let biosample_guid = app.subject_of_alignment(alignment_id).await.ok().flatten();

    // Resolve the reference we map to, with a download when that is necessary. Stream its progress
    // the same way as every other command that needs a reference.
    ensure_references_streaming(app, std::slice::from_ref(&target_build), evt_tx, &*wake).await;

    let reference = match app.cached_reference_path(&target_build) {
        Some(path) => path,
        None => {
            let _ = evt_tx.send(Event::RealignDone {
                alignment_id,
                biosample_guid,
                new_alignment_id: None,
                cancelled: false,
                summary: format!("the {target_build} reference is not available"),
            });
            wake();
            return;
        }
    };

    let tx = evt_tx.clone();
    let waker = wake.clone();
    let progress = move |p: navigator_app::realign_job::RealignProgress| {
        let _ = tx.send(Event::RealignProgress {
            alignment_id,
            biosample_guid,
            step: p.stage.step(),
            total: p.total_stages,
            label: p.stage.label().to_string(),
            detail: p.detail,
        });
        waker();
    };

    let params = navigator_app::realign_job::RealignParams {
        target_build: target_build.clone(),
        target_reference: reference,
        preset: None,
        scratch_root: None,
        // It is safe to opt in here, because the scratch path comes from this source alignment and
        // this target build. So anything in it belongs to the job about to run. An intermediate
        // survives at all only when something killed an earlier try outright. The machine went
        // down, the session went down, or somebody force-quit the process. In that case a user who
        // presses Realign again means "carry on", and not "spend four hours on the same file".
        resume: true,
    };

    let event = match app.realign_alignment(alignment_id, params, cancel, progress).await {
        Ok(outcome) => Event::RealignDone {
            alignment_id,
            biosample_guid,
            new_alignment_id: Some(outcome.alignment.id),
            cancelled: false,
            // A job that resumes steps over the stages that count these, so a figure can be truly
            // unknown. To say so is better than a zero the user would read as a result.
            summary: {
                let count = |n: Option<u64>| {
                    n.map(|n| n.to_string())
                        .unwrap_or_else(|| "an unrecorded number of".into())
                };
                format!(
                    "{} reads on {}; {} had been unplaced, {} duplicates marked",
                    count(outcome.reads_written),
                    target_build,
                    count(outcome.source_unmapped_reads),
                    count(outcome.duplicates_marked),
                )
            },
        },
        Err(e) => {
            let cancelled = e.is_cancelled();
            Event::RealignDone {
                alignment_id,
                biosample_guid,
                new_alignment_id: None,
                cancelled,
                summary: if cancelled { String::new() } else { e.to_string() },
            }
        }
    };
    let _ = evt_tx.send(event);
    wake();
}

/// Run the full analysis pipeline of one alignment. It emits an `AnalysisProgress` before each
/// step, and forwards the result event of each step, so that the detail tabs fill in live. It checks
/// `cancel` between steps. It also forwards the error of a step, and that error does not stop the
/// pipeline, which is best-effort.
async fn run_full_analysis_streaming<W: Fn() + Send + Sync + 'static>(
    app: &App,
    alignment_id: i64,
    include_ancestry: bool,
    cancel: CancelToken,
    evt_tx: &Sender<Event>,
    wake: Arc<W>,
) {
    // `App::plan_full_analysis` decides which steps run, and which ones go. A step goes when a
    // trusted external caller already placed this alignment, or when it has no chrM reads. That
    // plan is the one definition, and the CLI shares it. This function only turns those steps into
    // progress and result events.
    //
    // It plans twice, because the decision on the mitochondrion is a guess until step 1 gives
    // coverage. `include_sv = false`, because SV is experimental and costs hours for one
    // whole-genome sample, so no Full Analysis ever takes it in. The "Call SV" button on the
    // Sources tab runs it on request.
    let mut steps = app
        .plan_full_analysis(alignment_id, include_ancestry, false, None)
        .await
        .unwrap_or_else(|_| vec![AnalysisStep::QualityMetrics]);
    let mut total = steps.len();

    // Step 1: unified quality metrics. It does coverage and callable, read-level QC, and sex
    // inference in ONE pass over the alignment. Three separate steps used to read the file 2 to 3
    // times. This is the slow whole-genome read, so stream the sub-progress of each contig. The bar
    // then advances chromosome by chromosome, and does not stay at 0% for minutes.
    if !cancel.is_cancelled() {
        let _ = evt_tx.send(Event::AnalysisProgress {
            step: 1,
            total,
            label: "Quality metrics".into(),
            detail: "scanning contigs…".into(),
            fraction: 0.0,
        });
        wake();
        // Reuse the cached sub-results, instead of a second scan of the whole genome, which costs
        // minutes. Do it only when all three are there, because the unified walker persists them
        // together. The coverage must also cover the right scope: a stale whole-genome result for a
        // targeted-Y test reads as a miss. So the code computes it again, over the target contigs
        // only.
        let cached = match (
            app.cached_coverage_for_analysis(alignment_id).await,
            app.cached_read_metrics(alignment_id).await,
            app.cached_sex(alignment_id).await,
        ) {
            (Ok(Some(cov)), Ok(Some(rm)), Ok(Some(sex))) => Some((cov, rm, Some(sex))),
            _ => None,
        };
        let outcome = match cached {
            Some(triple) => Ok(triple),
            None => {
                // The parallel walker calls progress from worker threads, so the callback must be
                // Fn + Sync. The event Sender is !Sync, so guard it with a Mutex.
                let evt = Arc::new(Mutex::new(evt_tx.clone()));
                let wk = wake.clone();
                app.run_unified_metrics_with_progress(
                    alignment_id,
                    move |done, tot| {
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
                    },
                    cancel.clone(),
                )
                .await
                .map(|r| (r.coverage, r.read_metrics, r.sex))
            }
        };
        // Emit the same result events as the old separate steps, so that the UI updates the same
        // way.
        match outcome {
            Ok((cov, rm, sex)) => {
                // Plan again from the coverage that just ran. A Big Y with zero chrM reads then
                // correctly drops the chrM de-novo step and the mt-placement step. The pre-flight
                // plan above ran before this coverage existed. This also adjusts the total of the
                // steps that remain.
                steps = app
                    .plan_full_analysis(alignment_id, include_ancestry, false, Some(&cov))
                    .await
                    .unwrap_or_else(|_| std::mem::take(&mut steps));
                total = steps.len();
                // Pin a generic FTDNA Targeted-Y to Big Y-500 or Big Y-700, from its callable-chrY
                // footprint. It happens here, and not inside the metrics walker alone. It then
                // also fires on the cached fast path above, which steps over the walker and its
                // internal refine. When the generation changed, load the run card again, so that
                // it
                // shows the new label.
                if let Ok(Some(_)) = app.refine_big_y_generation_for_alignment(alignment_id, &cov).await {
                    if let Ok(guid) = app.biosample_of_alignment(alignment_id).await {
                        let _ = evt_tx.send(Event::RunsChanged(guid));
                    }
                }
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
                // Persist the failure, from a file that is corrupt or that the decoder can not
                // read. The project report then shows "Failed", and not a blank with no message.
                // This matches the CLI path and the batch path.
                app.record_analysis_error(alignment_id, "metrics", &e.to_string()).await;
                let _ = evt_tx.send(Event::Error(e.to_string()));
            }
        }
        wake();
    }

    // The steps that remain run through `handle`, which forwards the result events of each one. Y
    // variant discovery is the "private Y" pass behind the callable mask. It is NOT a raw
    // whole-chrY de-novo, which is enormous and mostly artifacts. chrM de-novo is acceptable,
    // because it is small and fully callable.
    let command_for = |step: &AnalysisStep| match step {
        // Step 1 ran above, with sub-progress; it is never dispatched as a command.
        AnalysisStep::QualityMetrics => None,
        AnalysisStep::StructuralVariants => Some(Command::RunSv(alignment_id)),
        AnalysisStep::MitoDenovo { contig } => Some(Command::RunDenovo {
            alignment_id,
            contig: contig.clone(),
        }),
        AnalysisStep::YHaplogroup => Some(Command::AssignYHaplogroup { alignment_id }),
        AnalysisStep::MtHaplogroup => Some(Command::AssignMtdnaHaplogroupFromAlignment { alignment_id }),
        AnalysisStep::YSignature { biosample_guid } => Some(Command::BuildYProfile {
            biosample_guid: *biosample_guid,
        }),
        AnalysisStep::AutosomalProfile { biosample_guid } => Some(Command::BuildAutosomalProfile {
            biosample_guid: *biosample_guid,
        }),
        AnalysisStep::Ancestry { biosample_guid } => Some(Command::EstimateAncestryFromConsensus {
            biosample_guid: *biosample_guid,
        }),
    };
    // Carry the position of each step in the plan, counted from 1. The progress numbers then stay
    // right, whatever steps the plan took in, and whichever ones this place dispatches instead of
    // the code above.
    let steps: Vec<(usize, String, String, Command)> = steps
        .iter()
        .enumerate()
        .filter_map(|(i, s)| command_for(s).map(|c| (i + 1, s.label().to_string(), s.detail().to_string(), c)))
        .collect();
    for (step, label, detail, cmd) in steps {
        if cancel.is_cancelled() {
            break;
        }
        let _ = evt_tx.send(Event::AnalysisProgress {
            step,
            total,
            label: label.to_string(),
            detail: detail.to_string(),
            fraction: (step as f32 - 1.0) / total as f32,
        });
        wake();
        // This runs to the end, and we can cancel before the next step. A step here goes around the
        // pre-flight that the outer match does for each command. So this is also the only place
        // where its failure can take up a diagnosis at file level.
        let ev = settle_alignment_command(app, alignment_id, handle(app, cmd, &cancel).await).await;
        let _ = evt_tx.send(ev);
        wake();
    }

    let _ = evt_tx.send(Event::AnalysisDone {
        cancelled: cancel.is_cancelled(),
    });
    wake();
}

/// Deep-analyze every sample in a project, one at a time. It emits a `DeepAnalyzeProgress` before
/// each sample, so that the bar advances sample by sample, then a final `ProjectAnalyzed`. It checks
/// `cancel` before each sample, and a stop leaves the artifacts that already exist in place, because
/// the pass is additive and idempotent. Each `analyze_biosample` awaits inside, so the worker
/// runtime stays free for a quick UI query between samples.
///
/// Run one workspace chore, and emit progress for each item.
///
/// The loop lives here, and not in `navigator-app`, for the same reason as the loop of
/// `deep_analyze_project`. The worker owns the event channel and the cancel token, and `App` stays
/// free of UI glue. What each item *does* is an `App` method, and the CLI shares it.
async fn run_chore_streaming(
    app: &App,
    chore: navigator_app::Chore,
    force: bool,
    cancel: CancelToken,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    use navigator_app::Chore;
    let mut outcome = navigator_app::ChoreOutcome::default();

    match chore {
        Chore::PrivateY => {
            let subjects = match app.list_all_biosamples().await {
                Ok(v) => v,
                Err(e) => return fail_chore(chore, e.to_string(), evt_tx, &wake),
            };
            let total = subjects.len();
            let mut novel = 0usize;
            for (i, b) in subjects.iter().enumerate() {
                if cancel.is_cancelled() {
                    break;
                }
                progress(chore, i, total, &b.donor_identifier, evt_tx, &wake);
                match app.refresh_private_y(b.guid, force).await {
                    Ok(r) => {
                        outcome.done += r.computed;
                        outcome.skipped += r.skipped + r.missing_file;
                        outcome.failed += r.failed;
                        novel += r.novel;
                    }
                    Err(e) => {
                        outcome.failed += 1;
                        let _ = evt_tx.send(Event::Error(format!("{}: {e}", b.donor_identifier)));
                        wake();
                    }
                }
            }
            outcome.summary = format!("{} computed · {} novel unique variants", outcome.done, novel);
        }

        Chore::StaleTree => {
            let targets = match app.stale_tree_targets(false).await {
                Ok(v) => v,
                Err(e) => return fail_chore(chore, e.to_string(), evt_tx, &wake),
            };
            let total = targets.len();
            let (mut calls_replaced, mut calls_failed, mut calls_skipped) = (0usize, 0usize, 0usize);
            // One lookup for the whole batch. It does not run a query for each subject only to
            // label a progress line.
            let names: std::collections::HashMap<_, _> = app
                .list_all_biosamples()
                .await
                .unwrap_or_default()
                .into_iter()
                .map(|b| (b.guid, b.donor_identifier))
                .collect();
            for (i, guid) in targets.iter().enumerate() {
                if cancel.is_cancelled() {
                    break;
                }
                let label = names.get(guid).cloned().unwrap_or_else(|| guid.0.to_string());
                progress(chore, i, total, &label, evt_tx, &wake);
                // Place the calls of each alignment again, *and* build the pooled profiles again.
                // See `App::replace_against_current_tree`. A build of the profiles alone, which is
                // what this used to do, left every `haplogroup_call` row on its old tree. That is
                // both the "sources diverge" conflicts on the Y card, and the reason a subject the
                // sweep touched stayed due for ever.
                match app.replace_against_current_tree(*guid).await {
                    Ok(r) => {
                        outcome.done += 1;
                        calls_replaced += r.calls_replaced;
                        calls_failed += r.calls_failed;
                        calls_skipped += r.calls_skipped;
                    }
                    Err(e) => {
                        outcome.failed += 1;
                        let _ = evt_tx.send(Event::Error(format!("{label}: {e}")));
                        wake();
                    }
                }
            }
            // A skip is separate from a failure in the report. "file gone" is normal in a
            // workspace where somebody cleaned out the vendor downloads. To put it into the error
            // count makes a healthy run look broken.
            outcome.summary = format!(
                "{} subject(s) re-placed against the current tree · {calls_replaced} call(s) re-placed{}{}",
                outcome.done,
                if calls_skipped > 0 {
                    format!(" · {calls_skipped} skipped (file gone)")
                } else {
                    String::new()
                },
                if calls_failed > 0 {
                    format!(" · {calls_failed} call(s) failed")
                } else {
                    String::new()
                }
            );
        }

        Chore::PublishOrigins => match app.publish_ancestral_origins(navigator_app::Lineage::Y, false).await {
            Ok(r) => {
                outcome.done = r.publishable;
                outcome.skipped = r.refused;
                outcome.summary = format!(
                    "{} queued · {} with a place, {} country only",
                    r.publishable, r.with_place, r.country_only
                );
            }
            Err(e) => return fail_chore(chore, e.to_string(), evt_tx, &wake),
        },
    }

    let _ = evt_tx.send(Event::ChoreDone { chore, outcome });
    wake();
}

fn progress(
    chore: navigator_app::Chore,
    done: usize,
    total: usize,
    label: &str,
    evt_tx: &Sender<Event>,
    wake: &Arc<dyn Fn() + Send + Sync>,
) {
    let _ = evt_tx.send(Event::ChoreProgress {
        chore,
        done,
        total,
        label: label.to_string(),
        fraction: if total > 0 { done as f32 / total as f32 } else { 0.0 },
    });
    wake();
}

/// A chore that could not start at all. Show the reason, and close it out, so that the UI never
/// leaves a spinner on a job that never began.
fn fail_chore(chore: navigator_app::Chore, err: String, evt_tx: &Sender<Event>, wake: &Arc<dyn Fn() + Send + Sync>) {
    let _ = evt_tx.send(Event::Error(err.clone()));
    let _ = evt_tx.send(Event::ChoreDone {
        chore,
        outcome: navigator_app::ChoreOutcome {
            summary: err,
            ..Default::default()
        },
    });
    wake();
}

async fn deep_analyze_project_streaming(
    app: &App,
    project_id: i64,
    cancel: CancelToken,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    let biosamples = match app.list_biosamples(project_id).await {
        Ok(v) => v,
        Err(e) => {
            let _ = evt_tx.send(Event::Error(e.to_string()));
            wake();
            return;
        }
    };
    let total = biosamples.len();
    let (mut samples, mut coverage_done, mut y_done, mut sex_done, mut metrics_done, mut errors) =
        (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);

    for (i, biosample) in biosamples.iter().enumerate() {
        if cancel.is_cancelled() {
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
        match app.analyze_biosample(biosample, cancel.clone()).await {
            Ok(o) if o.had_alignment => {
                samples += 1;
                coverage_done += o.coverage_done as usize;
                y_done += o.y_done as usize;
                sex_done += o.sex_done as usize;
                metrics_done += o.metrics_done as usize;
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
        errors,
        cancelled: cancel.is_cancelled(),
    });
    wake();
}

/// Import a NAS project directory. It emits an `ImportProgress` for each sample, so that an import
/// of 1000 samples shows a live status, and does not look frozen. It ends with a final
/// `ProjectImported`. A reference build that is missing comes back as `ReferenceNeeded`, which
/// prompts a download, and any other failure comes back as `Error`.
async fn import_project_dir_streaming(
    app: &App,
    dir: PathBuf,
    reference: Option<PathBuf>,
    evt_tx: &Sender<Event>,
    wake: Arc<dyn Fn() + Send + Sync>,
) {
    let progress = {
        let evt_tx = evt_tx.clone();
        let wake = wake.clone();
        move |done: usize, total: usize, sample: &str| {
            let _ = evt_tx.send(Event::ImportProgress {
                done,
                total,
                sample: sample.to_string(),
                fraction: if total > 0 { done as f32 / total as f32 } else { 0.0 },
            });
            wake();
        }
    };
    // Not `ev`: a missing reference is not an error but its own event, so this needs three arms.
    let event = match app
        .import_project_dir_with_progress(&dir, reference, "unknown".into(), true, progress)
        .await
    {
        Ok(summary) => Event::ProjectImported(summary),
        Err(AppError::ReferenceNeeded(builds)) => Event::ReferenceNeeded { dir, builds },
        Err(e) => Event::Error(e.to_string()),
    };
    let _ = evt_tx.send(event);
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

/// Narrate a subject's brief, streaming `BriefNarrationChunk`s as text arrives, then a final
/// `BriefNarration` (the authoritative result, or a plain-language failure reason).
async fn narrate_brief_streaming(app: &App, guid: SampleGuid, evt_tx: &Sender<Event>, wake: &(dyn Fn() + Send + Sync)) {
    let result = app
        .narrate_subject_streaming(guid, |text| {
            let _ = evt_tx.send(Event::BriefNarrationChunk {
                guid,
                text: text.to_string(),
            });
            wake();
        })
        .await
        .map_err(|e| e.to_string());
    let _ = evt_tx.send(Event::BriefNarration { guid, result });
    wake();
}

/// Answer a question about a subject's results, streaming `ChatAnswerChunk`s, then a final
/// `ChatAnswer`.
async fn ask_question_streaming(
    app: &App,
    guid: SampleGuid,
    history: Vec<ChatTurn>,
    question: String,
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    let result = app
        .answer_question_streaming(guid, history, question, |text| {
            let _ = evt_tx.send(Event::ChatAnswerChunk {
                guid,
                text: text.to_string(),
            });
            wake();
        })
        .await
        .map_err(|e| e.to_string());
    let _ = evt_tx.send(Event::ChatAnswer { guid, result });
    wake();
}

/// Explain a single result signal, streaming `SignalNarrationChunk`s as text arrives, then a final
/// `SignalNarration` (the authoritative result, or a plain-language failure reason).
async fn narrate_signal_streaming(
    app: &App,
    guid: SampleGuid,
    kind: SignalKind,
    evt_tx: &Sender<Event>,
    wake: &(dyn Fn() + Send + Sync),
) {
    let result = app
        .narrate_signal_streaming(guid, kind, |text| {
            let _ = evt_tx.send(Event::SignalNarrationChunk {
                guid,
                kind,
                text: text.to_string(),
            });
            wake();
        })
        .await
        .map_err(|e| e.to_string());
    let _ = evt_tx.send(Event::SignalNarration { guid, kind, result });
    wake();
}

/// Drain the outbox one time, and emit the outcome: one `Published` for each row it sent, the online
/// flag, and how many rows still wait.
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

/// Spawn the worker thread. It opens the workspace at `db_path` inside the runtime of the worker,
/// so that the connection pool lives there, then it serves commands. It calls `wake` after each
/// event, so that the UI can `request_repaint`. It returns the command sender, and the event
/// receiver, that the UI holds.
pub fn spawn(db_path: PathBuf, wake: impl Fn() + Send + Sync + 'static) -> (UnboundedSender<Command>, Receiver<Event>) {
    let (cmd_tx, mut cmd_rx) = unbounded_channel::<Command>();
    let (evt_tx, evt_rx) = std::sync::mpsc::channel::<Event>();
    let wake = Arc::new(wake);

    std::thread::Builder::new()
        .name("navigator-worker".into())
        .spawn(move || {
            // 64 MiB stacks. Two independent sources of deep recursion overflow the default 2 MiB
            // stack of a tokio thread. They abort the whole app in mid-batch.
            //
            // The first is the Y and mt tree parse, and the placement. Those recurse to the depth of
            // the haplotree (`flatten_du_node`, and the descent traversal), on a deep lineage, or on
            // the large FTDNA tree. The second is the CRAM decoder of noodles. That recurses on the
            // `spawn_blocking` decode paths, and goes deepest on a CRAM 3.1 file, with its new
            // range, fqzcomp and tokenizer codecs.
            //
            // A whole-genome record decode runs on `reader::decode_pool` instead, and this covers
            // the targeted decodes. See `NAVIGATOR_DECODE_STACK_MB`.
            let rt = match tokio::runtime::Builder::new_multi_thread()
                .enable_all()
                .thread_stack_size(64 * 1024 * 1024)
                .build()
            {
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
                // Registry of the in-flight cancellable run's token (see `CancelRegistry`).
                let cancels = CancelRegistry::default();

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
                    let cancels = cancels.clone();
                    tokio::spawn(async move {
                        match cmd {
                            // Streams ReferenceProgress events as bytes arrive, then a final event.
                            Command::ResolveReference { build } => {
                                resolve_reference_streaming(&app, build, &evt_tx, &*wake).await;
                            }
                            // Import, then resolve the references of the imported alignments at
                            // once, with a visible progress bar. A first CRAM or BAM that needs a
                            // multi-GB reference download otherwise looks as though it did not
                            // register. See `ensure_references_streaming`.
                            Command::AddDataBatch { biosample_guid, paths } => {
                                let event = handle(
                                    &app,
                                    Command::AddDataBatch { biosample_guid, paths },
                                    &CancelToken::none(),
                                )
                                .await;
                                let imported = matches!(event, Event::DataBatchImported { .. });
                                let _ = evt_tx.send(event);
                                wake();
                                if imported {
                                    if let Ok(builds) = app.reference_builds_for_subject(biosample_guid).await {
                                        ensure_references_streaming(&app, &builds, &evt_tx, &*wake).await;
                                    }
                                    ensure_indexes_for_subject_streaming(&app, biosample_guid, &evt_tx, &*wake).await;
                                }
                            }
                            Command::CreateSubjectAndImport {
                                donor_identifier,
                                sex,
                                paths,
                            } => {
                                let event = handle(
                                    &app,
                                    Command::CreateSubjectAndImport {
                                        donor_identifier,
                                        sex,
                                        paths,
                                    },
                                    &CancelToken::none(),
                                )
                                .await;
                                let guid = match &event {
                                    Event::SubjectCreatedAndImported { biosample_guid, .. } => Some(*biosample_guid),
                                    _ => None,
                                };
                                let _ = evt_tx.send(event);
                                wake();
                                if let Some(guid) = guid {
                                    if let Ok(builds) = app.reference_builds_for_subject(guid).await {
                                        ensure_references_streaming(&app, &builds, &evt_tx, &*wake).await;
                                    }
                                    ensure_indexes_for_subject_streaming(&app, guid, &evt_tx, &*wake).await;
                                }
                            }
                            // Resolve the reference and the coordinate index of the subject, or of
                            // the alignment, first, with a progress bar. An analysis that runs
                            // queries then does not start a download, or an index build, with no
                            // message, part of the way through.
                            Command::BuildAutosomalProfile { biosample_guid } => {
                                if let Ok(builds) = app.reference_builds_for_subject(biosample_guid).await {
                                    ensure_references_streaming(&app, &builds, &evt_tx, &*wake).await;
                                }
                                ensure_indexes_for_subject_streaming(&app, biosample_guid, &evt_tx, &*wake).await;
                                let event = handle(
                                    &app,
                                    Command::BuildAutosomalProfile { biosample_guid },
                                    &CancelToken::none(),
                                )
                                .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::StrConcordance { biosample_guid } => {
                                if let Ok(builds) = app.reference_builds_for_subject(biosample_guid).await {
                                    ensure_references_streaming(&app, &builds, &evt_tx, &*wake).await;
                                }
                                ensure_indexes_for_subject_streaming(&app, biosample_guid, &evt_tx, &*wake).await;
                                let event =
                                    handle(&app, Command::StrConcordance { biosample_guid }, &CancelToken::none())
                                        .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::RunSv(alignment_id) => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let (gen, token) = cancels.begin();
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(&app, Command::RunSv(alignment_id), &token).await,
                                )
                                .await;
                                cancels.end(gen);
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::LoadHeteroplasmy { alignment_id } => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(&app, Command::LoadHeteroplasmy { alignment_id }, &CancelToken::none())
                                        .await,
                                )
                                .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::AssignYHaplogroup { alignment_id } => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(&app, Command::AssignYHaplogroup { alignment_id }, &CancelToken::none())
                                        .await,
                                )
                                .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::YHaploReport { alignment_id } => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(&app, Command::YHaploReport { alignment_id }, &CancelToken::none()).await,
                                )
                                .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::AssignMtdnaHaplogroupFromAlignment { alignment_id } => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(
                                        &app,
                                        Command::AssignMtdnaHaplogroupFromAlignment { alignment_id },
                                        &CancelToken::none(),
                                    )
                                    .await,
                                )
                                .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::FindPrivateY { alignment_id, mask } => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(&app, Command::FindPrivateY { alignment_id, mask }, &CancelToken::none())
                                        .await,
                                )
                                .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::RunDenovo { alignment_id, contig } => {
                                if let Ok(Some(build)) = app.reference_build_of_alignment(alignment_id).await {
                                    ensure_references_streaming(&app, std::slice::from_ref(&build), &evt_tx, &*wake)
                                        .await;
                                }
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                let (gen, token) = cancels.begin();
                                let event = settle_alignment_command(
                                    &app,
                                    alignment_id,
                                    handle(&app, Command::RunDenovo { alignment_id, contig }, &token).await,
                                )
                                .await;
                                cancels.end(gen);
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::CompareIbdSources { a, b } => {
                                for src in [a, b] {
                                    if let navigator_app::IbdSource::Alignment(aln) = src {
                                        if let Ok(Some(build)) = app.reference_build_of_alignment(aln).await {
                                            ensure_references_streaming(
                                                &app,
                                                std::slice::from_ref(&build),
                                                &evt_tx,
                                                &*wake,
                                            )
                                            .await;
                                        }
                                        ensure_index_streaming(&app, aln, &evt_tx, &*wake).await;
                                    }
                                }
                                let event =
                                    handle(&app, Command::CompareIbdSources { a, b }, &CancelToken::none()).await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::BuildYProfile { biosample_guid } => {
                                if let Ok(builds) = app.reference_builds_for_subject(biosample_guid).await {
                                    ensure_references_streaming(&app, &builds, &evt_tx, &*wake).await;
                                }
                                ensure_indexes_for_subject_streaming(&app, biosample_guid, &evt_tx, &*wake).await;
                                let event =
                                    handle(&app, Command::BuildYProfile { biosample_guid }, &CancelToken::none()).await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            Command::BuildMtProfile { biosample_guid } => {
                                if let Ok(builds) = app.reference_builds_for_subject(biosample_guid).await {
                                    ensure_references_streaming(&app, &builds, &evt_tx, &*wake).await;
                                }
                                ensure_indexes_for_subject_streaming(&app, biosample_guid, &evt_tx, &*wake).await;
                                let event =
                                    handle(&app, Command::BuildMtProfile { biosample_guid }, &CancelToken::none())
                                        .await;
                                let _ = evt_tx.send(event);
                                wake();
                            }
                            // This streams an AnalysisProgress for each step, plus the result of
                            // each step, then AnalysisDone. Make sure the coordinate index exists
                            // first, because the walker of step 1 seeks by region on each contig.
                            Command::RunFullAnalysis { alignment_id } => {
                                ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                // Advanced: haplogroups/coverage only; ancestry is a separate action.
                                let (gen, cancel) = cancels.begin();
                                run_full_analysis_streaming(&app, alignment_id, false, cancel, &evt_tx, wake.clone())
                                    .await;
                                cancels.end(gen);
                            }
                            // Resolve the subject's representative alignment, then run the same pipeline
                            // as RunFullAnalysis. Lets the Simple view analyze from a subject guid alone.
                            Command::AnalyzeSubject { biosample_guid } => {
                                match app.default_alignment_for_subject(biosample_guid).await {
                                    Ok(Some((_run_id, alignment_id))) => {
                                        ensure_index_streaming(&app, alignment_id, &evt_tx, &*wake).await;
                                        // Simple one-click: include the autosomal ancestry step.
                                        let (gen, cancel) = cancels.begin();
                                        run_full_analysis_streaming(
                                            &app,
                                            alignment_id,
                                            true,
                                            cancel,
                                            &evt_tx,
                                            wake.clone(),
                                        )
                                        .await;
                                        cancels.end(gen);
                                    }
                                    Ok(None) => {
                                        let _ = evt_tx.send(Event::Error(
                                            "No sequencing data to analyze — import a BAM/CRAM for this person first."
                                                .into(),
                                        ));
                                        let _ = evt_tx.send(Event::AnalysisDone { cancelled: false });
                                        wake();
                                    }
                                    Err(e) => {
                                        let _ = evt_tx.send(Event::Error(e.to_string()));
                                        let _ = evt_tx.send(Event::AnalysisDone { cancelled: false });
                                        wake();
                                    }
                                }
                            }
                            // This streams a DeepAnalyzeProgress for each sample, then a final
                            // ProjectAnalyzed.
                            Command::DeepAnalyzeProject(project_id) => {
                                let (gen, cancel) = cancels.begin();
                                deep_analyze_project_streaming(&app, project_id, cancel, &evt_tx, wake.clone()).await;
                                cancels.end(gen);
                            }
                            // A workspace chore streams a ChoreProgress for each item, then
                            // ChoreDone. It walks the whole workspace, so a run inline would freeze
                            // the UI for minutes, with no sign of life.
                            Command::RunChore { chore, force } => {
                                let (gen, cancel) = cancels.begin();
                                run_chore_streaming(&app, chore, force, cancel, &evt_tx, wake.clone()).await;
                                cancels.end(gen);
                            }
                            // This streams an ImportProgress for each sample, then a final
                            // ProjectImported, or ReferenceNeeded, or Error. A large NAS import, of
                            // thousands of samples, otherwise looks frozen until the whole batch
                            // ends.
                            Command::ImportProjectDir { dir, reference } => {
                                import_project_dir_streaming(&app, dir, reference, &evt_tx, wake.clone()).await;
                            }
                            // Signals the in-flight full analysis / deep-analyze to stop between steps.
                            Command::StartRealign {
                                alignment_id,
                                target_build,
                            } => {
                                let (gen, cancel) = cancels.begin();
                                run_realign_streaming(&app, alignment_id, target_build, cancel, &evt_tx, wake.clone())
                                    .await;
                                cancels.end(gen);
                            }
                            Command::StartProjectRealign {
                                project_id,
                                target_build,
                            } => {
                                let queue = app
                                    .realignable_in_project(project_id, &target_build)
                                    .await
                                    .unwrap_or_default();
                                let queued = queue.len();
                                let mut completed = 0usize;
                                let mut abandoned = false;
                                for alignment_id in queue {
                                    let (gen, cancel) = cancels.begin();
                                    run_realign_streaming(
                                        &app,
                                        alignment_id,
                                        target_build.clone(),
                                        cancel.clone(),
                                        &evt_tx,
                                        wake.clone(),
                                    )
                                    .await;
                                    let was_cancelled = cancel.is_cancelled();
                                    cancels.end(gen);
                                    // A cancel abandons the queue, and not only the current
                                    // sample. Somebody who stops a multi-day batch means all of
                                    // it.
                                    if was_cancelled {
                                        abandoned = true;
                                        break;
                                    }
                                    completed += 1;
                                }
                                let _ = evt_tx.send(Event::RealignBatchDone {
                                    queued,
                                    completed,
                                    cancelled: abandoned,
                                });
                                wake();
                            }
                            Command::CancelRealign => {
                                cancels.cancel_current();
                                wake();
                            }
                            Command::CancelAnalysis => {
                                cancels.cancel_current();
                            }
                            // A publish goes into a durable queue, then drains, and it sends now
                            // if the app is online. The drain emits one Published for each row, and
                            // a pending count. We emit Queued for immediate feedback.
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
                            Command::PublishBiosample { biosample_guid } => {
                                let r = app.publish_biosample(biosample_guid).await;
                                publish_then_drain(&app, r, "biosample summary", &evt_tx, &*wake).await;
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
                            // This streams the narration text while the model writes it, then a
                            // final BriefNarration.
                            Command::NarrateBrief(guid) => {
                                narrate_brief_streaming(&app, guid, &evt_tx, &*wake).await;
                            }
                            // This streams a chat answer while the model writes it, then a final
                            // ChatAnswer.
                            Command::AskQuestion {
                                guid,
                                history,
                                question,
                            } => {
                                ask_question_streaming(&app, guid, history, question, &evt_tx, &*wake).await;
                            }
                            // This streams the explanation of one signal while the model writes
                            // it, then a final SignalNarration.
                            Command::NarrateSignal { guid, kind } => {
                                narrate_signal_streaming(&app, guid, kind, &evt_tx, &*wake).await;
                            }
                            other => {
                                let event = handle(&app, other, &CancelToken::none()).await;
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
    async fn handle_edits_genealogy_and_refreshes() {
        let app = app().await;
        let b = app.add_biosample(None, "GFX", None, None).await.unwrap();
        let guid = b.guid;

        // Add a vendor kit → the refreshed genealogy carries it.
        let ev = handle(
            &app,
            Command::AddExternalId {
                guid,
                source: "FTDNA".into(),
                external_id: "B5163".into(),
            },
            &CancelToken::none(),
        )
        .await;
        let data = match ev {
            Event::Genealogy { guid: g, data } if g == guid => data,
            other => panic!("expected Genealogy, got {other:?}"),
        };
        assert_eq!(data.external_ids.len(), 1);
        assert_eq!(data.external_ids[0].external_id, "B5163");
        let kit_id = data.external_ids[0].id;

        // Upsert a paternal MDKA → present in the refresh.
        let ev = handle(
            &app,
            Command::UpsertMdka {
                guid,
                mdka: NewMdka {
                    lineage: "Y".into(),
                    ancestor_name: Some("Thomas Kane".into()),
                    birth_year: Some(1830),
                    ..Default::default()
                },
            },
            &CancelToken::none(),
        )
        .await;
        match ev {
            Event::Genealogy { data, .. } => {
                let y = data.mdka.iter().find(|m| m.lineage == "Y").expect("Y mdka");
                assert_eq!(y.ancestor_name.as_deref(), Some("Thomas Kane"));
            }
            other => panic!("expected Genealogy, got {other:?}"),
        }

        // Delete both → empty genealogy.
        let _ = handle(
            &app,
            Command::DeleteMdka {
                guid,
                lineage: "Y".into(),
            },
            &CancelToken::none(),
        )
        .await;
        let ev = handle(
            &app,
            Command::DeleteExternalId { guid, id: kit_id },
            &CancelToken::none(),
        )
        .await;
        match ev {
            Event::Genealogy { data, .. } => assert!(data.is_empty(), "all genealogy removed"),
            other => panic!("expected Genealogy, got {other:?}"),
        }

        // A bind of a kit that another subject already owns comes back as an error.
        let c = app.add_biosample(None, "OTHER", None, None).await.unwrap();
        app.add_external_id(c.guid, "FTDNA", "B9999").await.unwrap();
        let ev = handle(
            &app,
            Command::AddExternalId {
                guid,
                source: "FTDNA".into(),
                external_id: "B9999".into(),
            },
            &CancelToken::none(),
        )
        .await;
        assert!(matches!(ev, Event::Error(_)), "conflicting id → Error, got {ev:?}");
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
            &CancelToken::none(),
        )
        .await;
        let pid = match created {
            Event::ProjectCreated(p) => p.id,
            other => panic!("expected ProjectCreated, got {other:?}"),
        };

        // overview reflects it
        match handle(&app, Command::LoadOverview, &CancelToken::none()).await {
            Event::Overview(v) => {
                assert_eq!(v.len(), 1);
                assert_eq!(v[0].sample_count, 0);
            }
            other => panic!("expected Overview, got {other:?}"),
        }

        // samples for the project (empty)
        match handle(&app, Command::LoadSamples(pid), &CancelToken::none()).await {
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
            .record_sequence_run(NewSequenceRun::new(b.guid, "ILLUMINA", "WGS"))
            .await
            .unwrap();
        let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../navigator-analysis/tests/fixtures");
        let aln = app
            .record_alignment(NewAlignment {
                bam_path: Some(fixtures.join("coverage.bam").to_string_lossy().into_owned()),
                reference_path: Some(fixtures.join("ref.fa").to_string_lossy().into_owned()),
                ..NewAlignment::new(run.id, "chrM", "synthetic")
            })
            .await
            .unwrap();

        // alignments query (keyed by run)
        match handle(&app, Command::LoadAlignments(run.id), &CancelToken::none()).await {
            Event::Alignments { alignments, .. } => assert_eq!(alignments, vec![aln.clone()]),
            other => panic!("expected Alignments, got {other:?}"),
        }

        // cold cache
        match handle(&app, Command::LoadCoverage(aln.id), &CancelToken::none()).await {
            Event::Coverage { result, .. } => assert!(result.is_none()),
            other => panic!("expected Coverage(None), got {other:?}"),
        }

        // run and persist (it uses the stored paths of the alignment, through `spawn_blocking`)
        match handle(&app, Command::RunCoverage(aln.id), &CancelToken::none()).await {
            Event::Coverage { alignment_id, result } => {
                assert_eq!(alignment_id, aln.id);
                assert_eq!(result.unwrap().genome_territory, 50);
            }
            other => panic!("expected Coverage(Some), got {other:?}"),
        }

        // now cached
        match handle(&app, Command::LoadCoverage(aln.id), &CancelToken::none()).await {
            Event::Coverage { result, .. } => assert_eq!(result.unwrap().callable_bases, 10),
            other => panic!("expected cached Coverage, got {other:?}"),
        }

        // de-novo on the fixture contig (cold -> run -> cached), keyed on the contig
        match handle(
            &app,
            Command::LoadDenovo {
                alignment_id: aln.id,
                contig: "chrM".into(),
            },
            &CancelToken::none(),
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
            &CancelToken::none(),
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
            &CancelToken::none(),
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
            &CancelToken::none(),
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        // it shows up both under the project and in the all-subjects list
        let guid = match handle(&app, Command::LoadSamples(pid), &CancelToken::none()).await {
            Event::Samples { samples, .. } => samples[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 2),
            other => panic!("got {other:?}"),
        }

        // add a run -> RunsChanged(sample)
        match handle(
            &app,
            Command::AddRun(NewSequenceRun::new(guid, "ILLUMINA", "WGS")),
            &CancelToken::none(),
        )
        .await
        {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        let run_id = match handle(&app, Command::LoadRuns(guid), &CancelToken::none()).await {
            Event::Runs { runs, .. } => runs[0].id,
            other => panic!("got {other:?}"),
        };

        // add an alignment -> AlignmentsChanged(run)
        match handle(
            &app,
            Command::AddAlignment(NewAlignment::new(run_id, "chm13v2.0", "bwa")),
            &CancelToken::none(),
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
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

        // new dependent data makes the delete refuse, with a conflict
        match handle(
            &app,
            Command::AddRun(NewSequenceRun::new(guid, "ILLUMINA", "WGS")),
            &CancelToken::none(),
        )
        .await
        {
            Event::RunsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::DeleteBiosample(guid), &CancelToken::none()).await {
            Event::Error(msg) => assert!(msg.contains("sequencing run"), "unexpected message: {msg}"),
            other => panic!("expected conflict Error, got {other:?}"),
        }
        // still present
        match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 1),
            other => panic!("got {other:?}"),
        }

        // a remove of the run clears the conflict, so the delete of the subject then works (the
        // end-to-end 'remove data first' path)
        let run_id = match handle(&app, Command::LoadRuns(guid), &CancelToken::none()).await {
            Event::Runs { runs, .. } => runs[0].id,
            other => panic!("got {other:?}"),
        };
        match handle(
            &app,
            Command::DeleteSequenceRun {
                id: run_id,
                biosample_guid: guid,
            },
            &CancelToken::none(),
        )
        .await
        {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::DeleteBiosample(guid), &CancelToken::none()).await {
            Event::BiosamplesChanged => {}
            other => panic!("expected clean delete, got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let spare = match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
            Event::AllBiosamples(all) => all.iter().find(|b| b.donor_identifier == "spare").unwrap().guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::DeleteBiosample(spare), &CancelToken::none()).await {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
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
            &CancelToken::none(),
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadSamples(pid), &CancelToken::none()).await {
            Event::Samples { samples, .. } => assert_eq!(samples.len(), 1),
            other => panic!("got {other:?}"),
        }

        // an assign to a project that does not exist gets a refusal
        match handle(
            &app,
            Command::AssignBiosampleProject {
                guid,
                project_id: Some(9999),
            },
            &CancelToken::none(),
        )
        .await
        {
            Event::Error(_) => {}
            other => panic!("expected Error, got {other:?}"),
        }

        // a clear of the project (None) removes it from the project list
        match handle(
            &app,
            Command::AssignBiosampleProject { guid, project_id: None },
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadSamples(pid), &CancelToken::none()).await {
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
            &CancelToken::none(),
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
            &CancelToken::none(),
        )
        .await
        {
            Event::ProjectsChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadOverview, &CancelToken::none()).await {
            Event::Overview(v) => {
                let p = &v.iter().find(|o| o.project.id == pid).unwrap().project;
                assert_eq!(p.name, "Renamed");
                assert_eq!(p.description.as_deref(), Some("a study"));
                assert_eq!(p.administrator, "curator");
            }
            other => panic!("got {other:?}"),
        }

        // A delete of a project with members now works. Its members detach, and the subjects stay.
        // The delete no longer gets a refusal.
        match handle(
            &app,
            Command::AddBiosample(NewBiosample {
                project_id: Some(pid),
                donor_identifier: "member".into(),
                sample_accession: None,
                sex: None,
            }),
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadSamples(pid), &CancelToken::none()).await {
            Event::Samples { samples, .. } => samples[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::DeleteProject(pid), &CancelToken::none()).await {
            Event::ProjectsChanged => {}
            other => panic!("expected clean delete, got {other:?}"),
        }
        // Project gone…
        match handle(&app, Command::LoadOverview, &CancelToken::none()).await {
            Event::Overview(v) => assert!(v.iter().all(|o| o.project.id != pid)),
            other => panic!("got {other:?}"),
        }
        // …but the detached subject still exists in the workspace.
        assert!(
            app.list_all_biosamples().await.unwrap().iter().any(|b| b.guid == guid),
            "subject should survive its project's deletion"
        );
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
            &CancelToken::none(),
        )
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        let guid = match handle(&app, Command::LoadAllBiosamples, &CancelToken::none()).await {
            Event::AllBiosamples(all) => all[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(
            &app,
            Command::AddRun(NewSequenceRun {
                total_reads: Some(1_000),
                ..NewSequenceRun::new(guid, "ILLUMINA", "WGS")
            }),
            &CancelToken::none(),
        )
        .await
        {
            Event::RunsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        let run = match handle(&app, Command::LoadRuns(guid), &CancelToken::none()).await {
            Event::Runs { runs, .. } => runs[0].clone(),
            other => panic!("got {other:?}"),
        };

        // edit the descriptive fields of the run, and the read metric stays
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
            &CancelToken::none(),
        )
        .await
        {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadRuns(guid), &CancelToken::none()).await {
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
            Command::AddAlignment(NewAlignment::new(run.id, "grch38", "bwa")),
            &CancelToken::none(),
        )
        .await
        {
            Event::AlignmentsChanged(_) => {}
            other => panic!("got {other:?}"),
        }
        let aln_id = match handle(&app, Command::LoadAlignments(run.id), &CancelToken::none()).await {
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
            &CancelToken::none(),
        )
        .await
        {
            Event::AlignmentsChanged(r) => assert_eq!(r, run.id),
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAlignments(run.id), &CancelToken::none()).await {
            Event::Alignments { alignments, .. } => {
                let a = &alignments[0];
                assert_eq!(a.reference_build, "chm13v2.0");
                assert_eq!(a.aligner, "minimap2");
                assert_eq!(a.variant_caller.as_deref(), Some("deepvariant"));
            }
            other => panic!("got {other:?}"),
        }
    }

    /// The streaming deep-analyze emits one progress event for each sample, then a final
    /// `ProjectAnalyzed`. A sample with no alignment that carries a BAM still goes through the walk,
    /// so that the bar advances, and the count leaves it out. That keeps the test free of any
    /// dependency on a reference, a network, or a tree.
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
        let wake: Arc<dyn Fn() + Send + Sync> = Arc::new(|| {});
        deep_analyze_project_streaming(&app, p.id, CancelToken::new(), &tx, wake).await;

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

        // A wake hook raises Cancel on the first progress event. That is what a user does when
        // they click Cancel after the first sample starts. Nothing arms the token again here. The
        // run no longer resets its own token at its start. That reset is exactly the race that used
        // to swallow a cancel between the spawn and the reset.
        let (tx, rx) = std::sync::mpsc::channel::<Event>();
        let cancel = CancelToken::new();
        let armed = cancel.clone();
        let wake: Arc<dyn Fn() + Send + Sync> = Arc::new(move || armed.cancel());
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
