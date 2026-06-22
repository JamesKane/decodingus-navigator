//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use crate::charts::{
    asset_status_line, coverage_histogram_chart, draw_ancestry_donut, draw_chromosome_painting, draw_composition_bar,
    draw_ibd_segments, draw_pca_scatter, draw_population_components, draw_variant_track, TrackRegion, VariantMark,
};
use crate::widgets::{
    capitalize_first, card, chip, combo, empty_state, fmt_depth, fmt_pct, fmt_reads, natural_cmp, opt, provider_abbrev,
    show_assignment, sortable_header, stat_card, variant_change, TableControls,
};
use eframe::egui;
use navigator_app::{
    AncestryResult, AncestrySegment, AppSettings, AuditEntry, BatchImportSummary, BuildNeed, CallState,
    CompatibilityLevel, Consensus, Coverage, DenovoCall, DnaType, FtdnaGenealogy, FtdnaImportPlan, FtdnaResolution,
    HaploAssignment, HeteroplasmySite, IbdComparison, IbdSuggestion, IdentityVerification, LineageBrief, LineageKind,
    MatchKind, MtRegion, MtVariant, PackStatus, PanelGenotype, PrivateBucket, PrivateClass, ProjectOverview,
    ProjectSampleReport, ProjectStrChart, ReadMetrics, RefBuildStatus, SexInferenceResult, SnpEvidence, SourceType,
    StrConcordanceRow, SubjectBrief, SvAnalysisResult, UiMode, VerificationStatus, YMatch, YProfile, YSignal, YState,
    YVariantStatus, YstrClustering,
};
use navigator_domain::chipprofile::{self, ChipProfile};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::mtdna::MtdnaSequence;
use navigator_domain::strpanel;
use navigator_domain::strprofile::{self, StrComparison, StrProfile};
use navigator_domain::testtype;
use navigator_domain::variants::VariantSet;
use navigator_domain::workspace::{Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, SequenceRun};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event, NewBiosample, PanelInfo, YMask};

#[derive(Default)]
struct Forms {
    /// Whether the inline "Add New Subject" form is expanded.
    show_add_subject: bool,
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

/// Primary navigation tabs in the app bar.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Nav {
    Dashboard,
    Subjects,
    Projects,
    Community,
}

/// Sub-tabs of the Community panel (the signed-in account's social surface).
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum CommunityTab {
    #[default]
    Support,
    Feed,
    Messages,
    Notifications,
}
impl CommunityTab {
    const ALL: [(CommunityTab, &'static str); 4] = [
        (CommunityTab::Support, "community.tab.support"),
        (CommunityTab::Feed, "community.tab.feed"),
        (CommunityTab::Messages, "community.tab.messages"),
        (CommunityTab::Notifications, "community.tab.notifications"),
    ];
}

/// Sub-tabs of the project detail panel (the member list vs the per-sample analysis report).
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum ProjectTab {
    #[default]
    Members,
    Report,
    Ystr,
}
impl ProjectTab {
    const ALL: [(ProjectTab, &'static str); 3] = [
        (ProjectTab::Members, "project.tab.members"),
        (ProjectTab::Report, "project.tab.report"),
        (ProjectTab::Ystr, "project.tab.ystr"),
    ];
}

/// Sub-tabs of the subject detail panel.
#[derive(Clone, Copy, PartialEq, Eq)]
enum DetailTab {
    Overview,
    YDna,
    MtDna,
    Autosomal,
    Ancestry,
    Sources,
    IbdMatches,
}

impl DetailTab {
    /// `(tab, i18n key)` in display order. The DNA-type tabs (Y / mt / Autosomal / Ancestry) show the
    /// subject's *consensus* across all sources; `Sources` is the per-sequencing-result hub.
    const ALL: [(DetailTab, &'static str); 7] = [
        (DetailTab::Overview, "detail.overview"),
        (DetailTab::YDna, "detail.ydna"),
        (DetailTab::MtDna, "detail.mtdna"),
        (DetailTab::Autosomal, "detail.autosomal"),
        (DetailTab::Ancestry, "detail.ancestry"),
        (DetailTab::Sources, "detail.sources"),
        (DetailTab::IbdMatches, "detail.ibd"),
    ];
}

/// Y-DNA sub-tabs: compact haplogroup landing, the heavy SNP surface, and STR (separated from SNP).
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum YSub {
    #[default]
    Haplogroup,
    Snp,
    Str,
}
impl YSub {
    const ALL: [(YSub, &'static str); 3] = [
        (YSub::Haplogroup, "detail.sub.haplogroup"),
        (YSub::Snp, "detail.sub.snp"),
        (YSub::Str, "detail.sub.str"),
    ];
}

/// Y-DNA → SNP variants nested sub-tabs: the heavy tables, one at a time (each runs to thousands of
/// rows on a WGS). The compact variant track stays above the bar as shared context.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum YSnpSub {
    #[default]
    Profile,
    Private,
    Imported,
}
impl YSnpSub {
    const ALL: [(YSnpSub, &'static str); 3] = [
        (YSnpSub::Profile, "detail.sub.yProfile"),
        (YSnpSub::Private, "detail.sub.privateY"),
        (YSnpSub::Imported, "detail.sub.imported"),
    ];
}

/// mtDNA sub-tabs.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum MtSub {
    #[default]
    Summary,
    Variants,
}
impl MtSub {
    const ALL: [(MtSub, &'static str); 2] = [
        (MtSub::Summary, "detail.sub.summary"),
        (MtSub::Variants, "detail.sub.variants"),
    ];
}

/// Autosomal sub-tabs: compact summary vs the heavy diploid profile table.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum AutoSub {
    #[default]
    Summary,
    Profile,
}
impl AutoSub {
    const ALL: [(AutoSub, &'static str); 2] = [
        (AutoSub::Summary, "detail.sub.summary"),
        (AutoSub::Profile, "detail.sub.profile"),
    ];
}

/// Which Y-STR report view is shown in the Y-DNA tab.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum StrReportView {
    /// FTDNA/YSEQ-style tier-grouped marker table.
    #[default]
    ByPanel,
    /// Flat, filterable marker table.
    AllMarkers,
    /// Cross-panel consensus value per marker.
    Consensus,
}

/// One editable reference-genome row in the Settings dialog.
#[derive(Clone)]
struct RefRow {
    build: String,
    status: String,
    local_path: String,
    auto_download: bool,
    /// Last integrity-check result for this build (set by `Event::ReferenceVerified`).
    verify: String,
}

/// Editable Settings-dialog state (loaded from `AppSettings`; reference rows arrive via
/// `Event::ReferenceSettings`).
#[derive(Clone)]
struct SettingsForm {
    appview_url: String,
    y_tree_provider: String, // "decodingus" | "ftdna"
    tree_ttl_days: String,
    prompt_before_download: bool,
    /// UI scale (egui zoom factor); 1.0 = native.
    ui_scale: f32,
    references: Vec<RefRow>,
    /// VCF-liftover tool state (input/output paths, target build, PAR filter).
    lift_in: String,
    lift_out: String,
    lift_target: String,
    lift_filter_par: bool,
}

impl SettingsForm {
    /// Scalar fields from the persisted `AppSettings` (reference rows filled later by the worker).
    fn from_settings() -> Self {
        let s = AppSettings::load();
        SettingsForm {
            appview_url: s.appview_url.unwrap_or_default(),
            y_tree_provider: s.y_tree_provider.unwrap_or_else(|| "decodingus".to_string()),
            tree_ttl_days: s
                .tree_ttl_days
                .map(|d| d.to_string())
                .unwrap_or_else(|| "7".to_string()),
            prompt_before_download: s.prompt_before_download.unwrap_or(true),
            ui_scale: s.ui_scale.unwrap_or(1.0),
            references: Vec::new(),
            lift_in: String::new(),
            lift_out: String::new(),
            lift_target: "chm13v2.0".to_string(),
            lift_filter_par: false,
        }
    }
}

/// The persisted UI scale (egui zoom factor), clamped to a sane range; `1.0` when unset.
fn resolved_ui_scale() -> f32 {
    AppSettings::load().ui_scale.unwrap_or(1.0).clamp(0.5, 3.0)
}

/// One-shot auto UI-scale probe (the "behave like a native app" default). On the first frame the
/// monitor size is known, derive a zoom when the OS reports a ~1.0 scale factor on a clearly
/// high-resolution panel (e.g. native-4K, where macOS itself doesn't up-scale). A Retina / scaled
/// display (native ppp > 1) is already handled by egui's native scaling, so it's left at 1.0. Skipped
/// entirely when a manual scale is persisted (`probed` starts `true`). The result fills the Settings
/// slider but is not persisted until the user saves — re-probed each launch otherwise.
fn run_auto_scale(probed: &mut bool, form: &mut SettingsForm, ctx: &egui::Context) {
    if *probed {
        return;
    }
    let Some(monitor_w) = ctx.input(|i| i.viewport().monitor_size).map(|s| s.x) else {
        return; // monitor size not reported yet — try again next frame
    };
    *probed = true;
    let native_ppp = ctx.native_pixels_per_point().unwrap_or(1.0);
    let physical_w = monitor_w * native_ppp; // monitor_size is logical points → ×ppp = physical px
    let auto = if native_ppp > 1.05 || physical_w < 3000.0 {
        1.0 // OS already HiDPI-scales, or the panel isn't hi-res enough to need help
    } else {
        (physical_w / 1920.0).clamp(1.0, 2.0) // ~2.0 on a 3840-wide 4K reported at scale 1.0
    };
    if (auto - 1.0).abs() > 0.01 {
        ctx.set_zoom_factor(auto);
        form.ui_scale = auto;
    }
}

/// In-flight full-analysis state, driving the modal dialog.
#[derive(Clone)]
struct AnalysisModal {
    step: usize,
    total: usize,
    label: String,
    detail: String,
    fraction: f32,
    /// egui time (seconds) when this step began — for the elapsed-time display.
    started: f64,
}

/// Editable copy of a project, driving the project Edit modal (Some ⇒ the dialog is shown).
#[derive(Clone)]
struct EditProject {
    id: i64,
    name: String,
    description: String,
    administrator: String,
}

/// Editable copy of a sequence run, driving the run Edit modal (Some ⇒ the dialog is shown).
/// Read-metric columns are not editable here, so they're not carried.
#[derive(Clone)]
struct EditRun {
    id: i64,
    guid: SampleGuid,
    test_type: String,
    platform_name: String,
    instrument_model: String,
    library_layout: String,
    sequencing_facility: String,
}

/// Drives the destructive merge-sequence-runs modal (Some ⇒ shown). `secondary` is the run that
/// will be emptied + deleted; `primary` is the chosen merge target (its picker default).
#[derive(Clone)]
struct MergeRuns {
    guid: SampleGuid,
    secondary: i64,
    primary: Option<i64>,
}

/// Editable copy of an alignment, driving the alignment Edit modal (Some ⇒ the dialog is shown).
#[derive(Clone)]
struct EditAlignment {
    id: i64,
    run_id: i64,
    reference_build: String,
    aligner: String,
    variant_caller: String,
}

/// Editable copy of a subject, driving the Edit modal (Some ⇒ the dialog is shown).
#[derive(Clone)]
struct EditSubject {
    guid: SampleGuid,
    donor_identifier: String,
    sample_accession: String,
    description: String,
    center_name: String,
    sex: String,
}

/// A data-source row pending delete confirmation. The `label` is shown in the confirm dialog;
/// the variant carries the ids the worker command needs (and the parent id to refresh).
#[derive(Clone)]
enum DataDelete {
    Run { id: i64, guid: SampleGuid, label: String },
    Alignment { id: i64, run_id: i64, label: String },
    Str { id: i64, guid: SampleGuid, label: String },
    Variant { id: i64, guid: SampleGuid, label: String },
    Chip { id: i64, guid: SampleGuid, label: String },
    Mtdna { id: i64, guid: SampleGuid, label: String },
}

impl DataDelete {
    fn label(&self) -> &str {
        match self {
            DataDelete::Run { label, .. }
            | DataDelete::Alignment { label, .. }
            | DataDelete::Str { label, .. }
            | DataDelete::Variant { label, .. }
            | DataDelete::Chip { label, .. }
            | DataDelete::Mtdna { label, .. } => label,
        }
    }

    fn command(&self) -> Command {
        match *self {
            DataDelete::Run { id, guid, .. } => Command::DeleteSequenceRun {
                id,
                biosample_guid: guid,
            },
            DataDelete::Alignment { id, run_id, .. } => Command::DeleteAlignment {
                id,
                sequence_run_id: run_id,
            },
            DataDelete::Str { id, guid, .. } => Command::DeleteStrProfile {
                id,
                biosample_guid: guid,
            },
            DataDelete::Variant { id, guid, .. } => Command::DeleteVariantSet {
                id,
                biosample_guid: guid,
            },
            DataDelete::Chip { id, guid, .. } => Command::DeleteChipProfile {
                id,
                biosample_guid: guid,
            },
            DataDelete::Mtdna { id, guid, .. } => Command::DeleteMtdnaSequence {
                id,
                biosample_guid: guid,
            },
        }
    }
}

/// Reference population PC1/PC2 centroids for the PCA scatter: `(population_code, pc1, pc2)`.
type PcaCentroids = Vec<(String, f64, f64)>;

/// A full Y-haplogroup placement report for one alignment: ranked candidates + lineage SNP evidence.
struct YReport {
    alignment_id: i64,
    assignment: HaploAssignment,
    lineage: Vec<SnpEvidence>,
}

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    /// In-flight full-analysis progress (Some ⇒ the modal dialog is shown).
    analysis: Option<AnalysisModal>,
    /// Subject being edited (Some ⇒ the Edit modal is shown).
    edit_subject: Option<EditSubject>,
    /// Subject pending delete confirmation (Some ⇒ the confirm dialog is shown).
    confirm_delete: Option<SampleGuid>,
    /// Subject pending "clear all data" confirmation (Some ⇒ the confirm dialog is shown).
    confirm_clear: Option<SampleGuid>,
    /// The last batch-import summary, shown in a modal until dismissed.
    batch_import: Option<BatchImportSummary>,
    /// Y-STR-from-sequence concordance for the selected subject: `(guid, source alignment, rows)`.
    str_concordance: Option<(SampleGuid, i64, Vec<StrConcordanceRow>)>,
    /// Whether a Y-STR-from-sequence call is in flight (the heavy first pass).
    str_running: bool,
    /// Cross-subject Y matches for the selected subject: `(guid, ranked matches)`. Gap §2.
    y_matches: Option<(SampleGuid, Vec<YMatch>)>,
    /// Whether a Y-match search is in flight.
    y_matches_running: bool,
    /// Project filter for the Y-match search (None ⇒ whole workspace).
    y_match_project: Option<i64>,
    /// Text filter over the Y-match table (by donor / haplogroup).
    y_match_query: String,
    /// Data-source row pending delete confirmation (Some ⇒ the confirm dialog is shown).
    confirm_data_delete: Option<DataDelete>,
    /// Subject being assigned to a project: (subject, selected project or None). Some ⇒ picker shown.
    assign_project: Option<(SampleGuid, Option<i64>)>,
    /// Project being edited (Some ⇒ the project Edit modal is shown).
    edit_project: Option<EditProject>,
    /// Project pending delete confirmation: (id, name). Some ⇒ the confirm dialog is shown.
    confirm_delete_project: Option<(i64, String)>,
    /// Sequence run being edited (Some ⇒ the run Edit modal is shown).
    edit_run: Option<EditRun>,
    merge_runs: Option<MergeRuns>,
    /// Whether the read-only Y-profile source-audit modal is open (reads the cached `y_profile`).
    audit_y_profile: bool,
    /// Alignment being edited (Some ⇒ the alignment Edit modal is shown).
    edit_alignment: Option<EditAlignment>,
    /// Current frame's egui time (seconds), captured at the top of `update`.
    frame_time: f64,
    /// Selected primary navigation tab.
    nav: Nav,
    /// Interface mode: Simple (casual single-person briefs) vs. Advanced (full power-user UI).
    ui_mode: UiMode,
    /// Whether the mode was explicitly pinned (env / settings / user toggle). When `false` the
    /// first-run workspace heuristic may still adjust it as data loads.
    ui_mode_pinned: bool,
    /// Precomputed Simple-mode brief for the selected subject `(guid, brief)`; `None` until built.
    subject_brief: Option<(SampleGuid, SubjectBrief)>,
    /// Whether a Subject Brief build is in flight.
    subject_brief_loading: bool,
    /// Free-text filter over the Simple-mode "My DNA" subject selector (matches the donor name).
    simple_subject_filter: String,
    /// Selected subject-detail sub-tab.
    detail_tab: DetailTab,
    /// Active UI language.
    lang: crate::i18n::Lang,
    /// Dark (default) vs light theme.
    dark_mode: bool,
    /// Whether the one-shot auto-UI-scale probe has run (skipped when a manual scale is persisted).
    scale_probed: bool,
    /// Settings dialog open + its editable form.
    show_settings: bool,
    settings_form: SettingsForm,
    /// Sort + inline per-column filter state for the subjects table.
    subjects_table_ctl: TableControls,
    /// Collapse the subjects side panel to a thin strip so the detail panel (charts/tables)
    /// gets the full width.
    subjects_collapsed: bool,
    /// Collapse the projects side panel to a thin strip, handing the detail panel the full width.
    projects_collapsed: bool,
    overview: Vec<ProjectOverview>,
    selected_project: Option<i64>,
    /// Per-sample coverage/haplogroup report rows for the selected project.
    project_report: Vec<ProjectSampleReport>,
    /// Precomputed Y-STR overview (FTDNA-style chart) for the selected project; `None` until the
    /// background build returns. A boolean tracks the in-flight build so the UI can show a spinner.
    project_str_chart: Option<ProjectStrChart>,
    project_str_loading: bool,
    samples: Vec<Biosample>,
    /// Every biosample (the project-independent subjects list).
    all_biosamples: Vec<Biosample>,
    /// Per-subject Y/mt terminal haplogroups for the list columns (`guid → (Y, mt)`).
    haplo_summary: std::collections::HashMap<SampleGuid, (Option<String>, Option<String>)>,
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
    /// Y-STR report view-state: which view, which provider (when multiple), and the marker filter.
    str_report_view: StrReportView,
    str_provider: Option<String>,
    str_marker_filter: String,
    /// SNP variant sets for the selected subject.
    variant_sets: Vec<VariantSet>,
    /// Chip/array profiles for the selected subject.
    chip_profiles: Vec<ChipProfile>,
    /// mtDNA sequences for the selected subject.
    mtdna_sequences: Vec<MtdnaSequence>,
    /// rCRS-relative mutation lists per mtDNA sequence id (loaded on demand).
    mtdna_variants: std::collections::HashMap<i64, Vec<MtVariant>>,
    /// Last mtDNA haplogroup assignment: (sequence id, assignment).
    mtdna_haplogroup: Option<(i64, HaploAssignment)>,
    /// Last Y haplogroup assignment: (alignment id, assignment).
    y_haplogroup: Option<(i64, HaploAssignment)>,
    /// Active sub-tab within the Y-DNA / mtDNA / Autosomal detail tabs.
    y_sub: YSub,
    y_snp_sub: YSnpSub,
    mt_sub: MtSub,
    auto_sub: AutoSub,
    /// Full Y placement report (ranked candidates + lineage SNP evidence) for an alignment.
    y_report: Option<YReport>,
    /// True while the haplogroup report is being built.
    y_report_running: bool,
    /// Last mtDNA-from-alignment haplogroup assignment: (alignment id, assignment).
    mt_haplogroup: Option<(i64, HaploAssignment)>,
    /// Ancestry/IBD reference-asset presence + integrity (the "data sources" line). Loaded once.
    asset_status: Vec<navigator_app::AssetStatus>,
    /// Donor-level ancestry (best across the subject's sources): (source alignment id, result).
    donor_ancestry: Option<(i64, AncestryResult)>,
    /// Detailed consensus ancestry reports: modern fine-population + ancient-component breakdowns.
    fine_ancestry: Option<AncestryResult>,
    ancient_ancestry: Option<AncestryResult>,
    nmonte_ancestry: Option<AncestryResult>,
    /// Reference PC1/PC2 centroids for the PCA scatter, keyed by alignment_id (lazy-loaded).
    pca_reference: Option<(i64, PcaCentroids)>,
    /// Which PCA-reference key we've already dispatched a load for (avoids re-sending every frame).
    pca_reference_attempted: Option<i64>,
    /// Donor-level private-Y union across the subject's sources.
    donor_private_y: Option<PrivateBucket>,
    /// The selected subject's multi-source Y-variant profile.
    y_profile: Option<YProfile>,
    /// Y-variant profile status filter (None = all).
    y_profile_filter: Option<YVariantStatus>,
    /// Text search across the variant/SNP tables (by SNP name / site), per table.
    y_profile_query: String,
    mt_profile_query: String,
    auto_profile_query: String,
    private_y_query: String,
    str_seq_query: String,
    /// Catalogued Y-SNP names at variant positions (`position → name`), used to annotate the two
    /// Y-SNP tables' position-only / novel calls. Resolved once per subject from the Y-SNP dictionary.
    y_snp_names: std::collections::HashMap<i64, String>,
    /// True once we've dispatched the Y-SNP-name resolution for the current subject (avoids re-sending).
    y_snp_names_requested: bool,
    /// True while the (expensive) Y-variant profile is being built.
    y_profile_loading: bool,
    /// The selected subject's multi-source mtDNA consensus profile.
    mt_profile: Option<navigator_app::ConsensusProfile>,
    /// mtDNA consensus-profile status filter (None = all).
    mt_profile_filter: Option<YVariantStatus>,
    /// True while the (expensive) mtDNA consensus profile is being built.
    mt_profile_loading: bool,
    /// The selected subject's multi-source autosomal (diploid 0/1/2) consensus profile.
    auto_profile: Option<navigator_app::DiploidProfile>,
    /// Autosomal consensus status filter (None = all).
    auto_profile_filter: Option<YVariantStatus>,
    /// True while the (expensive) autosomal consensus profile is being built.
    auto_profile_loading: bool,
    /// Whether the consensus-driven donor ancestry estimate is in flight.
    estimating_donor_ancestry: bool,
    /// Local-ancestry painting: (alignment id, segments). `painting_running` while genotyping.
    painting: Option<(i64, Vec<AncestrySegment>)>,
    painting_running: bool,
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
    /// An alignment to auto-select once its run's alignments load (subject-centric default).
    pending_alignment: Option<i64>,
    coverage: Option<Coverage>,
    /// Cached coverage per alignment for the selected run's Data Sources rows (so each row shows
    /// coverage/callable without first selecting that alignment). Keyed by alignment id.
    coverage_by_aln: std::collections::HashMap<i64, Coverage>,
    /// Genome-region metadata (cytoband ideogram) for the selected alignment's build, `(alignment_id,
    /// regions)`. Lazily fetched when the Ideogram tab is opened.
    genome_regions: Option<(i64, std::sync::Arc<navigator_app::GenomeRegions>)>,
    /// True while the cytoBand fetch is in flight.
    loading_regions: bool,
    /// The alignment we've already kicked off (or completed) a region load for — avoids re-firing
    /// the fetch every frame, including after a failure.
    regions_attempted: Option<i64>,
    /// Which contig's depth histogram the coverage view charts: `None` = whole-genome histogram,
    /// `Some(i)` = `coverage.contig_coverage_stats[i]`.
    coverage_hist_contig: Option<usize>,
    sex: Option<SexInferenceResult>,
    read_metrics: Option<ReadMetrics>,
    sv: Option<SvAnalysisResult>,
    running_sex: bool,
    running_metrics: bool,
    running_sv: bool,
    running: bool,
    /// De-novo haploid SNP calls keyed by contig (chrY on the Y-DNA tab, chrM on the mtDNA tab).
    denovo: std::collections::HashMap<String, Vec<DenovoCall>>,
    running_denovo: bool,
    panels: Vec<PanelInfo>,
    selected_panel: Option<i64>,
    all_alignments: Vec<Alignment>,
    panel_genotypes: Option<Vec<PanelGenotype>>,
    running_genotype: bool,
    ibd_other: Option<i64>,
    /// Chip-compatible IBD compare: the two picked sources (each a WGS alignment or an imported chip).
    ibd_src_a: Option<navigator_app::IbdSource>,
    ibd_src_b: Option<navigator_app::IbdSource>,
    /// Subject-level (consensus) IBD compare: the other subject picked for comparison.
    ibd_other_subject: Option<SampleGuid>,
    ibd_result: Option<IbdComparison>,
    running_ibd: bool,
    /// Identity-verification result for the current IBD pair.
    identity: Option<IdentityVerification>,
    /// Federated IBD: pseudonymous match suggestions fetched from the AppView.
    ibd_suggestions: Vec<IbdSuggestion>,
    /// Whether a suggestions fetch is in flight (drives the spinner).
    loading_ibd_suggestions: bool,
    /// Per-candidate introduction status, keyed by `suggested_sample_guid` (e.g. "PENDING").
    ibd_intros: std::collections::HashMap<String, String>,
    /// Encrypted-exchange inbox: inbound requests awaiting consent + consent-ready sessions.
    exchange_incoming: Vec<navigator_app::IncomingRequest>,
    exchange_ready: Vec<navigator_app::ExchangeSessionInfo>,
    /// The selected subject's persisted IBD exchange results.
    exchange_results: Vec<navigator_app::StoredIbdExchange>,
    /// True while an inbox refresh / consent / exchange run is in flight.
    exchange_busy: bool,
    /// Signed-in account DID, or `None`. Gates the "Publish" actions.
    account: Option<String>,
    /// Whether the last PDS write reached the server (offline indicator).
    online: bool,
    /// Outbox rows still awaiting a successful push (the "N pending" sync indicator).
    sync_pending: i64,
    /// True while a PULL reconcile is in flight.
    pulling: bool,
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
    /// Streaming deep-analyze progress: `(done, total, current_sample, fraction)` while running.
    deep_progress: Option<(usize, usize, String, f32)>,
    /// The dry-run FTDNA import plan being reviewed (drives the review modal).
    ftdna_plan: Option<FtdnaImportPlan>,
    /// The admin's per-kit resolutions for the fuzzy rows in [`Self::ftdna_plan`].
    ftdna_resolutions: std::collections::BTreeMap<String, FtdnaResolution>,
    /// The selected subject's imported genealogy (vendor ids + FTDNA member + MDKA), for the
    /// Overview card. `(guid, data)` so a stale bundle from a prior subject isn't shown.
    genealogy: Option<(SampleGuid, FtdnaGenealogy)>,
    /// The current project's Y-STR clustering, keyed by project id (so a stale one isn't shown).
    project_clustering: Option<(i64, YstrClustering)>,
    /// True while the project Y-STR clustering is computing.
    clustering_running: bool,
    /// Active project detail sub-tab (Members vs Report).
    project_tab: ProjectTab,
    /// Filter for the project Members list (kit / name / branch substring).
    member_filter: String,
    /// Sort + inline per-column filter state for the project Report table.
    report_table_ctl: TableControls,
    // ---- Community (social) ------------------------------------------------
    /// Active Community sub-tab (Support / Feed / Notifications).
    community_tab: CommunityTab,
    /// The signed-in account's support threads.
    support_threads: Vec<navigator_app::SocialThreadSummary>,
    /// The opened thread's `(conversation_id, messages)` — `None` when viewing the list.
    open_thread: Option<(String, Vec<navigator_app::SocialMessage>)>,
    /// The loaded community feed.
    feed: Option<navigator_app::FeedView>,
    /// Loaded notifications + the server's unread count (the app-bar bell badge).
    notifications: Vec<navigator_app::SocialNotification>,
    notif_unread: i64,
    /// Whether the social tab has fetched at least once this session (lazy first load).
    community_loaded: bool,
    /// Peer DMs (social 3a). Inbox: inbound DM requests + consent-ready sessions to connect.
    dm_incoming: Vec<navigator_app::IncomingRequest>,
    dm_ready: Vec<navigator_app::ExchangeSessionInfo>,
    /// Persisted conversation list + the opened conversation's `(session_id, transcript)`.
    dm_conversations: Vec<navigator_app::DmConversationSummary>,
    open_dm: Option<(String, Vec<navigator_app::DmMessage>)>,
    /// Composer buffers: start-a-DM partner DID + the open conversation's message draft.
    dm_partner_did: String,
    dm_compose: String,
    /// Whether the Messages sub-tab has loaded at least once this session (lazy first load).
    dm_loaded: bool,
    /// Recruitment 3c: the signed-in account's open recruitment invitations (shown in Notifications).
    recruitment_invitations: Vec<navigator_app::RecruitmentInvitation>,
    /// Composer buffers: new-thread subject/body, open-thread reply, feed post + topic.
    new_thread_subject: String,
    new_thread_body: String,
    thread_reply: String,
    feed_content: String,
    feed_topic: String,
    /// Opt-in: also publish the next community post to the signed-in PDS as a federated
    /// `feed.post` record (roadmap 3b). Off by default — publishing to your own repo is an
    /// explicit, portable public act.
    feed_publish_pds: bool,
    forms: Forms,
    status: String,
}

/// Sentinel option in the chip-provider dropdown that means "let the parser guess".
const AUTO_DETECT: &str = "(auto-detect)";

/// The workbench accent (primary buttons, selection, active tabs).
pub(crate) const ACCENT: egui::Color32 = egui::Color32::from_rgb(45, 125, 246);
/// Destructive-action red (Delete) — used by the subject header in Phase 2.
#[allow(dead_code)]
const DANGER: egui::Color32 = egui::Color32::from_rgb(220, 60, 60);

/// Rows shown before the heavy variant/site tables (Y/mt/autosomal consensus profiles, private-Y,
/// de-novo SNPs) scroll internally. On a WGS these run thousands of rows; bounding them keeps the
/// detail page navigable instead of forcing endless scrolling — but the pane must be tall enough to
/// be useful (the user wants 20-30 rows, not a 3-row slot).
const PROFILE_TABLE_ROWS: usize = 26;

/// Explicit height for a scrollable profile table that shows up to [`PROFILE_TABLE_ROWS`] rows of the
/// given `count`, then scrolls. Sized to the content when there are fewer rows (no empty pane).
/// Computed (rather than a fixed constant) because a nested `ScrollArea` with vertical `auto_shrink`
/// collapses to a few rows inside the page scroll — we pair this with `auto_shrink([false, false])`.
fn profile_pane_height(ui: &egui::Ui, count: usize) -> f32 {
    let row_h = ui.text_style_height(&egui::TextStyle::Body) + ui.spacing().item_spacing.y + 3.0;
    let rows = count.clamp(1, PROFILE_TABLE_ROWS) as f32;
    row_h * (rows + 1.0) // +1 for the header row
}

/// Apply the Decoding-Us workbench look: a dark (or light) palette with the accent blue,
/// rounded widgets, and roomier spacing — the visual base that closes most of the gap to the
/// Scala Workbench. Re-applied on theme toggle.
fn apply_theme(ctx: &egui::Context, dark: bool) {
    use egui::{Color32, Rounding, Stroke};
    // Pin the preference so egui stops following the OS theme — otherwise `theme_preference`
    // defaults to `System` and our styled visuals get clobbered by the host (e.g. a light macOS
    // would show light even with Dark selected in Settings).
    ctx.set_theme(if dark {
        egui::ThemePreference::Dark
    } else {
        egui::ThemePreference::Light
    });
    let mut style = (*ctx.style()).clone();
    let mut v = if dark {
        egui::Visuals::dark()
    } else {
        egui::Visuals::light()
    };

    if dark {
        v.panel_fill = Color32::from_gray(27);
        v.window_fill = Color32::from_gray(32);
        v.extreme_bg_color = Color32::from_gray(20); // text-edit / table body
        v.faint_bg_color = Color32::from_gray(38); // striped rows / cards
        v.override_text_color = Some(Color32::from_gray(224));
        v.widgets.noninteractive.bg_fill = Color32::from_gray(32);
        v.widgets.inactive.bg_fill = Color32::from_gray(52);
        v.widgets.inactive.weak_bg_fill = Color32::from_gray(44);
        v.widgets.hovered.bg_fill = Color32::from_gray(64);
        v.widgets.active.bg_fill = ACCENT;
        v.window_stroke = Stroke::new(1.0, Color32::from_gray(48));
    }
    v.hyperlink_color = ACCENT;
    v.selection.bg_fill = ACCENT.gamma_multiply(0.55);
    v.selection.stroke = Stroke::new(1.0, ACCENT);

    let r = Rounding::same(6.0);
    for w in [
        &mut v.widgets.noninteractive,
        &mut v.widgets.inactive,
        &mut v.widgets.hovered,
        &mut v.widgets.active,
        &mut v.widgets.open,
    ] {
        w.rounding = r;
    }
    v.window_rounding = Rounding::same(10.0);
    v.menu_rounding = r;

    style.visuals = v;
    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(10.0, 6.0);
    ctx.set_style(style);
}

/// Subjects-table columns: `(header, width)`.
const SUBJECT_COLS: [(&str, f32); 6] = [
    ("Name", 180.0),
    ("Y-DNA", 150.0),
    ("mtDNA", 110.0),
    ("Sex", 70.0),
    ("Center", 130.0),
    ("Status", 90.0),
];

mod central;
mod chrome;
mod community;
mod detail;
mod events;
mod ibd;
mod modals;
mod sources;

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
        let _ = tx.send(Command::BackfillLabs); // resolve labs for runs imported before D8 landed
        let _ = tx.send(Command::VerifySourceFiles); // flag any imported file that moved/disappeared
        let _ = tx.send(Command::LoadAssetStatus); // ancestry/IBD "data sources" line
                                                   // Persisted theme wins; default dark. (Must match `dark_mode` below.)
        let dark = !matches!(AppSettings::load().theme.as_deref(), Some("light"));
        apply_theme(&cc.egui_ctx, dark);
        // Persisted UI scale (egui zoom) — fixes tiny text on a native-4K display the OS reports at
        // scale factor 1.0. egui's keyboard zoom (Cmd +/-/0) also works but isn't persisted.
        cc.egui_ctx.set_zoom_factor(resolved_ui_scale());
        NavigatorApp {
            tx,
            rx,
            analysis: None,
            edit_subject: None,
            confirm_delete: None,
            confirm_clear: None,
            batch_import: None,
            str_concordance: None,
            str_running: false,
            y_matches: None,
            y_matches_running: false,
            y_match_project: None,
            y_match_query: String::new(),
            confirm_data_delete: None,
            assign_project: None,
            edit_project: None,
            confirm_delete_project: None,
            edit_run: None,
            merge_runs: None,
            audit_y_profile: false,
            edit_alignment: None,
            frame_time: 0.0,
            nav: Nav::Subjects,
            // Pinned mode (env / settings) wins; else default Simple provisionally and let the
            // first-run workspace heuristic adjust once subjects/projects load (see
            // `apply_ui_mode_heuristic`).
            ui_mode: navigator_app::configured_ui_mode().unwrap_or(UiMode::Simple),
            ui_mode_pinned: navigator_app::configured_ui_mode().is_some(),
            subject_brief: None,
            subject_brief_loading: false,
            simple_subject_filter: String::new(),
            detail_tab: DetailTab::Overview,
            // Persisted choice wins; else honor $LANG (e.g. "es_ES.UTF-8") when it names a
            // supported locale; else English.
            lang: crate::i18n::load_lang()
                .or_else(|| std::env::var("LANG").ok().and_then(|l| crate::i18n::Lang::parse(&l)))
                .unwrap_or(crate::i18n::Lang::En),
            // Persisted theme wins; default dark.
            dark_mode: !matches!(AppSettings::load().theme.as_deref(), Some("light")),
            // A persisted manual scale takes precedence; otherwise probe the monitor on frame 1.
            scale_probed: AppSettings::load().ui_scale.is_some(),
            show_settings: false,
            settings_form: SettingsForm::from_settings(),
            subjects_table_ctl: TableControls::sorted_by(0),
            subjects_collapsed: false,
            projects_collapsed: false,
            overview: Vec::new(),
            selected_project: None,
            project_report: Vec::new(),
            project_str_chart: None,
            project_str_loading: false,
            samples: Vec::new(),
            all_biosamples: Vec::new(),
            haplo_summary: std::collections::HashMap::new(),
            selected_sample: None,
            runs: Vec::new(),
            consensus_y: None,
            consensus_mt: None,
            audit_y: Vec::new(),
            audit_mt: Vec::new(),
            heteroplasmy: None,
            str_profiles: Vec::new(),
            str_report_view: StrReportView::default(),
            str_provider: None,
            str_marker_filter: String::new(),
            variant_sets: Vec::new(),
            chip_profiles: Vec::new(),
            asset_status: Vec::new(),
            mtdna_sequences: Vec::new(),
            mtdna_variants: std::collections::HashMap::new(),
            mtdna_haplogroup: None,
            y_haplogroup: None,
            y_sub: YSub::default(),
            y_snp_sub: YSnpSub::default(),
            mt_sub: MtSub::default(),
            auto_sub: AutoSub::default(),
            y_report: None,
            y_report_running: false,
            mt_haplogroup: None,
            donor_ancestry: None,
            fine_ancestry: None,
            ancient_ancestry: None,
            nmonte_ancestry: None,
            pca_reference: None,
            pca_reference_attempted: None,
            donor_private_y: None,
            y_profile: None,
            y_profile_filter: None,
            y_profile_query: String::new(),
            mt_profile_query: String::new(),
            auto_profile_query: String::new(),
            private_y_query: String::new(),
            str_seq_query: String::new(),
            y_snp_names: std::collections::HashMap::new(),
            y_snp_names_requested: false,
            y_profile_loading: false,
            mt_profile: None,
            mt_profile_filter: None,
            mt_profile_loading: false,
            auto_profile: None,
            auto_profile_filter: None,
            auto_profile_loading: false,
            estimating_donor_ancestry: false,
            painting: None,
            painting_running: false,
            private_y: None,
            finding_private_y: false,
            y_mask_path: None,
            y_self_mask: true,
            selected_run: None,
            alignments: Vec::new(),
            selected_alignment: None,
            pending_alignment: None,
            coverage: None,
            coverage_by_aln: std::collections::HashMap::new(),
            genome_regions: None,
            loading_regions: false,
            regions_attempted: None,
            coverage_hist_contig: None,
            sex: None,
            read_metrics: None,
            sv: None,
            running_sex: false,
            running_metrics: false,
            running_sv: false,
            running: false,
            denovo: std::collections::HashMap::new(),
            running_denovo: false,
            panels: Vec::new(),
            selected_panel: None,
            all_alignments: Vec::new(),
            panel_genotypes: None,
            running_genotype: false,
            ibd_other: None,
            ibd_src_a: None,
            ibd_src_b: None,
            ibd_other_subject: None,
            ibd_result: None,
            running_ibd: false,
            identity: None,
            ibd_suggestions: Vec::new(),
            loading_ibd_suggestions: false,
            ibd_intros: std::collections::HashMap::new(),
            exchange_incoming: Vec::new(),
            exchange_ready: Vec::new(),
            exchange_results: Vec::new(),
            exchange_busy: false,
            account: None,
            online: true,
            sync_pending: 0,
            pulling: false,
            logging_in: false,
            publishing: false,
            importing: false,
            pending_import_dir: None,
            reference_needs: Vec::new(),
            reference_progress: None,
            analyzing: false,
            deep_progress: None,
            ftdna_plan: None,
            ftdna_resolutions: std::collections::BTreeMap::new(),
            genealogy: None,
            project_clustering: None,
            clustering_running: false,
            project_tab: ProjectTab::default(),
            member_filter: String::new(),
            report_table_ctl: TableControls::sorted_by(0),
            community_tab: CommunityTab::default(),
            support_threads: Vec::new(),
            open_thread: None,
            feed: None,
            notifications: Vec::new(),
            notif_unread: 0,
            community_loaded: false,
            dm_incoming: Vec::new(),
            dm_ready: Vec::new(),
            dm_conversations: Vec::new(),
            open_dm: None,
            dm_partner_did: String::new(),
            dm_compose: String::new(),
            dm_loaded: false,
            recruitment_invitations: Vec::new(),
            new_thread_subject: String::new(),
            new_thread_body: String::new(),
            thread_reply: String::new(),
            feed_content: String::new(),
            feed_topic: String::new(),
            feed_publish_pds: false,
            forms: Forms {
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
}

impl eframe::App for NavigatorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.frame_time = ctx.input(|i| i.time);
        run_auto_scale(&mut self.scale_probed, &mut self.settings_form, ctx);
        // While an analysis runs, keep repainting so the spinner/elapsed timer animate even
        // during a long step that emits no events (e.g. whole-genome coverage).
        if self.analysis.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(120));
        }
        self.drain_events();
        // Mode upkeep: first-run heuristic (until pinned) + auto-select the sole subject in Simple.
        self.apply_ui_mode_heuristic();
        self.auto_select_single_subject();
        self.handle_file_drops(ctx);
        self.app_bar(ctx);
        self.nav_bar(ctx);
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.online {
                    ui.colored_label(egui::Color32::from_rgb(80, 190, 120), self.tr("status.online"));
                } else {
                    ui.colored_label(egui::Color32::from_rgb(220, 150, 60), self.tr("status.offline"));
                }
                if self.sync_pending > 0 {
                    ui.separator();
                    ui.colored_label(
                        egui::Color32::from_rgb(220, 150, 60),
                        format!("⟳ {} {}", self.sync_pending, self.tr("status.pending")),
                    )
                    .on_hover_text(self.tr("status.pendingHint"));
                }
                // PULL reconcile — a PDS-repo op, so it needs a real PDS (OAuth) account. A local
                // did:key identity (federation/exchange only) has no PDS repo → show it disabled with
                // a reason rather than letting it fail with a confusing "not signed in".
                if let Some(acct) = self.account.clone() {
                    let is_pds = !acct.starts_with("did:key:");
                    ui.separator();
                    let btn = ui.add_enabled(is_pds && !self.pulling, egui::Button::new(self.tr("sync.pull")).small());
                    let btn = btn.on_hover_text(self.tr(if is_pds { "sync.pullHint" } else { "sync.pullNeedsPds" }));
                    if is_pds && btn.clicked() {
                        self.pulling = true;
                        self.status = self.tr("sync.pulling").to_string();
                        let _ = self.tx.send(Command::PullSync);
                    }
                    if self.pulling {
                        ui.spinner();
                    }
                }
                ui.separator();
                ui.label(egui::RichText::new(self.tr("status.label")).weak());
                ui.label(&self.status);
            });
        });
        // The action bar's batch/compare/add-to-project affordances are power-user features —
        // Advanced only.
        if self.nav == Nav::Subjects && self.ui_mode == UiMode::Advanced {
            self.action_bar(ctx);
        }
        self.left_panel(ctx);
        egui::CentralPanel::default().show(ctx, |ui| match self.nav {
            Nav::Dashboard => self.dashboard_central(ui),
            Nav::Subjects => self.subjects_central(ui),
            Nav::Projects => self.projects_central(ui),
            Nav::Community => self.community_central(ui),
        });
        self.analysis_modal(ctx);
        self.edit_subject_modal(ctx);
        self.delete_subject_modal(ctx);
        self.clear_subject_modal(ctx);
        self.data_delete_modal(ctx);
        self.assign_project_modal(ctx);
        self.edit_project_modal(ctx);
        self.delete_project_modal(ctx);
        self.edit_run_modal(ctx);
        self.merge_runs_modal(ctx);
        self.y_profile_audit_modal(ctx);
        self.edit_alignment_modal(ctx);
        self.settings_modal(ctx);
        self.batch_import_modal(ctx);
        self.ftdna_review_modal(ctx);
        self.paint_drop_hint(ctx);
    }
}

/// Render a depth histogram (`bin d` = bases observed at depth `d`, top bin = ≥255) as an
const STR_CONFLICT: egui::Color32 = egui::Color32::from_rgb(220, 150, 60);

/// FTDNA/YSEQ-style By-Panel view: markers grouped into tiers (Y-12 / Y-25 / …), each tier rendered
/// as transposed mini-grids (marker-name row over value row, ≤12 markers wide). Conflicting markers
/// are amber.
fn str_by_panel_view(ui: &mut egui::Ui, profile: &StrProfile, provider: &str, comparison: &StrComparison) {
    let conflicts: std::collections::HashSet<String> = comparison
        .conflicts
        .iter()
        .map(|c| c.marker.trim().to_uppercase())
        .collect();
    let groups = strpanel::assign_markers_to_panels(&profile.markers, provider);
    if groups.is_empty() {
        ui.label(egui::RichText::new("No markers.").weak());
        return;
    }
    let canon = strpanel::canonical_provider(provider);
    // No inner scroll area here — the detail panel is already wrapped in one vertical
    // ScrollArea, and nesting a second (fixed-height) one clips the panel tables and
    // captures the wheel so the page can't scroll. Let the tiers flow into the page scroll.
    for (tier, markers) in &groups {
        ui.add_space(6.0);
        ui.label(egui::RichText::new(format!("{canon} {tier}  ({} markers)", markers.len())).strong());
        for (ci, chunk) in markers.chunks(12).enumerate() {
            egui::Grid::new(format!("str_tier_{tier}_{ci}"))
                .num_columns(chunk.len())
                .show(ui, |ui| {
                    for mk in chunk {
                        let t = egui::RichText::new(&mk.marker).small();
                        let c = conflicts.contains(&mk.marker.trim().to_uppercase());
                        ui.label(if c { t.color(STR_CONFLICT) } else { t.weak() });
                    }
                    ui.end_row();
                    for mk in chunk {
                        let t = egui::RichText::new(&mk.value).monospace().strong();
                        let c = conflicts.contains(&mk.marker.trim().to_uppercase());
                        ui.label(if c { t.color(STR_CONFLICT) } else { t });
                    }
                    ui.end_row();
                });
        }
    }
}

/// Flat, filterable marker table: Marker | Panel | Value, plus ⚠ | Other when >1 provider (the
/// other provider's disagreeing value). Conflicting values are amber.
fn str_all_markers_view(
    ui: &mut egui::Ui,
    profile: &StrProfile,
    provider: &str,
    comparison: &StrComparison,
    filter: &mut String,
) {
    use std::collections::HashMap;
    let conflict_map: HashMap<String, &strprofile::MarkerConflict> = comparison
        .conflicts
        .iter()
        .map(|c| (c.marker.trim().to_uppercase(), c))
        .collect();
    let mut tier_of: HashMap<String, String> = HashMap::new();
    for (tier, ms) in strpanel::assign_markers_to_panels(&profile.markers, provider) {
        for mk in ms {
            tier_of.insert(mk.marker.trim().to_uppercase(), tier.clone());
        }
    }
    let multi = comparison.providers.len() > 1;
    let this_provider = profile.provider.clone().unwrap_or_default();

    ui.horizontal(|ui| {
        ui.label("Filter:");
        ui.add(
            egui::TextEdit::singleline(filter)
                .hint_text("marker")
                .desired_width(120.0),
        );
        if !filter.is_empty() && ui.button("✕").clicked() {
            filter.clear();
        }
    });
    let f = filter.trim().to_uppercase();
    let cols = if multi { 5 } else { 3 };
    // Flow into the detail panel's outer ScrollArea (no nested vertical scroll — it clips
    // the table and steals the wheel).
    egui::Grid::new("str_all_grid")
        .striped(true)
        .num_columns(cols)
        .show(ui, |ui| {
            ui.strong("Marker");
            ui.strong("Panel");
            ui.strong("Value");
            if multi {
                ui.strong("⚠");
                ui.strong("Other");
            }
            ui.end_row();
            for mk in &profile.markers {
                let norm = mk.marker.trim().to_uppercase();
                if !f.is_empty() && !norm.contains(&f) {
                    continue;
                }
                let conflict = conflict_map.get(&norm).copied();
                ui.label(&mk.marker);
                ui.label(egui::RichText::new(tier_of.get(&norm).map(|s| s.as_str()).unwrap_or("—")).weak());
                let v = egui::RichText::new(&mk.value).monospace();
                ui.label(if conflict.is_some() { v.color(STR_CONFLICT) } else { v });
                if multi {
                    match conflict {
                        Some(c) => {
                            ui.colored_label(STR_CONFLICT, "⚠");
                            let others: Vec<String> = c
                                .by_provider
                                .iter()
                                .filter(|(p, _)| p != &this_provider)
                                .map(|(p, val)| format!("{p}:{val}"))
                                .collect();
                            ui.label(egui::RichText::new(others.join(", ")).monospace().weak());
                        }
                        None => {
                            ui.label("");
                            ui.label("");
                        }
                    }
                }
                ui.end_row();
            }
        });
}

/// Canonical short label + badge color for a consensus status, shared by the Y/mt consensus card
/// and the autosomal diploid card. `Novel` can't arise on the diploid (panel-site) path — see
/// [`navigator_domain::consensus::reconcile_diploid`] — but is mapped here for completeness.
fn consensus_status_badge(status: YVariantStatus) -> (&'static str, egui::Color32) {
    let amber = egui::Color32::from_rgb(220, 150, 60);
    match status {
        YVariantStatus::Confirmed => ("confirmed", egui::Color32::from_rgb(120, 180, 120)),
        YVariantStatus::Novel => ("novel", egui::Color32::from_rgb(120, 150, 220)),
        YVariantStatus::Conflict => ("conflict", amber),
        YVariantStatus::SingleSource => ("single", egui::Color32::from_gray(150)),
        YVariantStatus::Pending => ("pending", egui::Color32::from_gray(150)),
        YVariantStatus::NoCoverage => ("no-cov", egui::Color32::from_gray(110)),
    }
}

/// Shared renderer for a multi-source consensus profile (Y or mtDNA — same generic engine): header
/// (counts + lineage label + provenance), a status filter, and the per-variant grid. `variant_col`
/// names the identity column ("SNP" / "Mutation"); `kind` labels the empty state; `id_salt` keeps the
/// two cards' scroll/grid ids distinct. `snp_names` annotates a position-only/novel row with the
/// catalogued Y-SNP name at that site (empty for mtDNA, whose mutations are already named).
#[allow(clippy::too_many_arguments)]
fn draw_consensus_profile(
    ui: &mut egui::Ui,
    profile: &navigator_app::ConsensusProfile,
    filter: &mut Option<YVariantStatus>,
    query: &mut String,
    variant_col: &str,
    kind: &str,
    id_salt: &str,
    snp_names: &std::collections::HashMap<i64, String>,
) {
    if profile.variants.is_empty() {
        ui.label(egui::RichText::new(format!("No {kind} across sources.")).weak());
        return;
    }
    let s = &profile.summary;
    let mut header = format!(
        "{} confirmed · {} novel · {} conflict · {} single-source · confidence {:.0}%",
        s.confirmed,
        s.novel,
        s.conflict,
        s.single_source,
        s.overall_confidence * 100.0
    );
    if let Some(t) = &profile.terminal {
        header = format!("terminal {t}   —   {header}");
    }
    ui.label(egui::RichText::new(header).weak());
    // Provenance: which tests contributed (label · count).
    if !profile.sources.is_empty() {
        let prov = profile
            .sources
            .iter()
            .map(|src| format!("{} ({})", src.label, src.variant_count))
            .collect::<Vec<_>>()
            .join(" · ");
        ui.label(egui::RichText::new(format!("sources: {prov}")).weak().small());
    }

    let amber = egui::Color32::from_rgb(220, 150, 60);
    ui.horizontal(|ui| {
        ui.label("Show:");
        ui.selectable_value(filter, None, "All");
        ui.selectable_value(filter, Some(YVariantStatus::Conflict), "Conflicts");
        ui.selectable_value(filter, Some(YVariantStatus::Novel), "Novel");
        ui.selectable_value(filter, Some(YVariantStatus::Confirmed), "Confirmed");
        ui.add(
            egui::TextEdit::singleline(query)
                .hint_text("filter SNP / pos")
                .desired_width(140.0),
        );
        if !query.is_empty() && ui.small_button("✕").clicked() {
            query.clear();
        }
    });
    let q = query.to_ascii_lowercase();

    let state_label = |st: YState| match st {
        YState::Derived => "derived",
        YState::Ancestral => "ancestral",
        YState::NoCall => "no-call",
    };
    // Bound the table to a fixed-height scroll pane — on a WGS these run thousands of rows and would
    // otherwise force endless page scrolling. Status/text filters narrow the list; a cap bounds a
    // pathological profile. The header/filter row above stays fixed; only the grid scrolls.
    const CAP: usize = 2000;
    let (mut shown, mut total_match) = (0usize, 0usize);
    let pane_h = profile_pane_height(ui, profile.variants.len());
    egui::ScrollArea::vertical()
        .id_salt(format!("{id_salt}_scroll"))
        .max_height(pane_h)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Grid::new(format!("{id_salt}_grid"))
                .striped(true)
                .num_columns(5)
                .show(ui, |ui| {
                    for h in [variant_col, "Pos", "State", "Status", "Sources"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for v in &profile.variants {
                        if filter.is_some_and(|f| v.status != f) {
                            continue;
                        }
                        // A catalogued Y-SNP name at this site (for a position-only / novel row).
                        let cataloged = if v.name.is_empty() {
                            snp_names.get(&v.position).map(String::as_str)
                        } else {
                            None
                        };
                        if !q.is_empty()
                            && !v.name.to_ascii_lowercase().contains(&q)
                            && !v.position.to_string().contains(&q)
                            && !cataloged.is_some_and(|n| n.to_ascii_lowercase().contains(&q))
                        {
                            continue;
                        }
                        total_match += 1;
                        if shown >= CAP {
                            continue;
                        }
                        shown += 1;
                        {
                            let conflict = v.status == YVariantStatus::Conflict;
                            // Prefer the consensus name; else the catalogued Y-SNP at this site (teal,
                            // tooltipped); else a bare position marker.
                            let teal = egui::Color32::from_rgb(90, 190, 190);
                            let (display, is_cat) = match (v.name.is_empty(), cataloged) {
                                (false, _) => (v.name.clone(), false),
                                (true, Some(c)) => (c.to_string(), true),
                                (true, None) => (format!("novel@{}", v.position), false),
                            };
                            let mut name_txt = egui::RichText::new(display).strong();
                            if conflict {
                                name_txt = name_txt.color(amber);
                            } else if is_cat {
                                name_txt = name_txt.color(teal);
                            }
                            let resp = ui.label(name_txt);
                            if is_cat {
                                resp.on_hover_text("catalogued Y-SNP at this site (not on the placed lineage)");
                            }
                            ui.label(egui::RichText::new(v.position.to_string()).weak());
                            ui.label(state_label(v.consensus));
                            let (label, color) = consensus_status_badge(v.status);
                            ui.colored_label(color, format!("{label} ({}/{})", v.support, v.total));
                            ui.horizontal(|ui| {
                                for src in &v.sources {
                                    let short = match src.source_type {
                                        SourceType::Chip => "chip",
                                        SourceType::WgsShortRead | SourceType::WgsLongRead => "WGS",
                                        SourceType::Sanger => "Sanger",
                                        SourceType::Imported => "seq",
                                        _ => "src",
                                    };
                                    let glyph = match src.state {
                                        YState::Derived => "✓",
                                        YState::Ancestral => "·",
                                        YState::NoCall => "?",
                                    };
                                    ui.label(egui::RichText::new(format!("{short}{glyph}")).small().weak())
                                        .on_hover_text(format!("{}: {}", src.label, state_label(src.state)));
                                }
                            });
                            ui.end_row();
                        }
                    }
                });
        });
    ui.label(
        egui::RichText::new(format!(
            "{} of {} matching variants",
            shown.min(total_match),
            total_match
        ))
        .weak()
        .small(),
    );
    if total_match > CAP {
        ui.label(egui::RichText::new(format!("…and {} more — filter to narrow", total_match - CAP)).weak());
    }
}

/// Renderer for the autosomal diploid consensus profile — the 0/1/2 sibling of
/// [`draw_consensus_profile`]. Header (confirmed/conflict/single + confidence), a status filter, and a
/// per-site grid `Site (rsID) | GT (0/0,0/1,1/1) | Status | Sources (per-source dosage)`.
fn draw_diploid_profile(
    ui: &mut egui::Ui,
    profile: &navigator_app::DiploidProfile,
    filter: &mut Option<YVariantStatus>,
    query: &mut String,
) {
    if profile.variants.is_empty() {
        ui.label(egui::RichText::new("No autosomal sites across sources.").weak());
        return;
    }
    let s = &profile.summary;
    ui.label(
        egui::RichText::new(format!(
            "{} sites · {} confirmed · {} conflict · {} single-source · confidence {:.0}%",
            s.total,
            s.confirmed,
            s.conflict,
            s.single_source,
            s.overall_confidence * 100.0
        ))
        .weak(),
    );
    if !profile.sources.is_empty() {
        let prov = profile
            .sources
            .iter()
            .map(|src| format!("{} ({})", src.label, src.variant_count))
            .collect::<Vec<_>>()
            .join(" · ");
        ui.label(egui::RichText::new(format!("sources: {prov}")).weak().small());
    }

    let amber = egui::Color32::from_rgb(220, 150, 60);
    ui.horizontal(|ui| {
        ui.label("Show:");
        ui.selectable_value(filter, None, "All");
        ui.selectable_value(filter, Some(YVariantStatus::Conflict), "Conflicts");
        ui.selectable_value(filter, Some(YVariantStatus::Confirmed), "Confirmed");
        ui.add(
            egui::TextEdit::singleline(query)
                .hint_text("filter rsID / site")
                .desired_width(140.0),
        );
        if !query.is_empty() && ui.small_button("✕").clicked() {
            query.clear();
        }
    });
    let q = query.to_ascii_lowercase();
    let matches = |v: &navigator_app::DiploidVariant| {
        q.is_empty() || v.name.to_ascii_lowercase().contains(&q) || format!("{}:{}", v.contig, v.position).contains(&q)
    };

    // Dosage 0/1/2 → diploid genotype string; -1 = no-call.
    let gt = |d: i8| match d {
        0 => "0/0",
        1 => "0/1",
        2 => "1/1",
        _ => "./.",
    };
    // The panel has ~1.2M sites — a non-virtualized Grid lays out every row per frame and beach-balls.
    // Render fixed-width columns through ScrollArea::show_rows, which only builds the visible slice.
    const W_SITE: f32 = 150.0;
    const W_GT: f32 = 44.0;
    const W_STATUS: f32 = 130.0;
    let row_h = ui.text_style_height(&egui::TextStyle::Body) + 4.0;
    ui.horizontal(|ui| {
        ui.add_sized([W_SITE, row_h], egui::Label::new(egui::RichText::new("Site").strong()));
        ui.add_sized([W_GT, row_h], egui::Label::new(egui::RichText::new("GT").strong()));
        ui.add_sized(
            [W_STATUS, row_h],
            egui::Label::new(egui::RichText::new("Status").strong()),
        );
        ui.label(egui::RichText::new("Sources").strong());
    });
    let render_row = |ui: &mut egui::Ui, v: &navigator_app::DiploidVariant| {
        let conflict = v.status == YVariantStatus::Conflict;
        let name = if v.name.is_empty() {
            format!("{}:{}", v.contig, v.position)
        } else {
            v.name.clone()
        };
        let name_txt = egui::RichText::new(name).strong();
        ui.horizontal(|ui| {
            ui.add_sized(
                [W_SITE, row_h],
                egui::Label::new(if conflict { name_txt.color(amber) } else { name_txt }).truncate(),
            );
            ui.add_sized([W_GT, row_h], egui::Label::new(gt(v.consensus_dosage)));
            let (lbl, color) = consensus_status_badge(v.status);
            ui.add_sized(
                [W_STATUS, row_h],
                egui::Label::new(egui::RichText::new(format!("{lbl} ({}/{})", v.support, v.total)).color(color)),
            );
            for src in &v.sources {
                let short = match src.source_type {
                    SourceType::Chip => "chip",
                    SourceType::WgsShortRead | SourceType::WgsLongRead => "WGS",
                    _ => "src",
                };
                ui.label(egui::RichText::new(format!("{short}{}", gt(src.dosage))).small().weak())
                    .on_hover_text(format!("{}: {}", src.label, gt(src.dosage)));
            }
        });
    };
    // Bound the rows to a fixed-height scroll pane (the panel is ~1.2M sites). A hard cap keeps only
    // `shown` widgets built per frame, never the full panel; the status/text filters narrow further.
    const CAP: usize = 2000;
    let (mut shown, mut total_match) = (0usize, 0usize);
    let pane_h = profile_pane_height(ui, profile.variants.len());
    egui::ScrollArea::vertical()
        .id_salt("diploid_profile_scroll")
        .max_height(pane_h)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            for v in &profile.variants {
                if !(filter.map_or(true, |f| v.status == f) && matches(v)) {
                    continue;
                }
                total_match += 1;
                if shown >= CAP {
                    continue;
                }
                shown += 1;
                render_row(ui, v);
            }
        });
    ui.label(
        egui::RichText::new(format!("{} of {} matching sites", shown.min(total_match), total_match))
            .weak()
            .small(),
    );
    if total_match > CAP {
        ui.label(egui::RichText::new(format!("…and {} more — filter to narrow", total_match - CAP)).weak());
    }
}

/// The rCRS-relative mtDNA mutation list, grouped by region (HVR2 / Coding / HVR1) — the classic
/// mtDNA result. `variants` are derived against the bundled rCRS; notation is standard mtDNA form.
fn mtdna_mutations_view(ui: &mut egui::Ui, mtdna_id: i64, variants: &[MtVariant]) {
    if variants.is_empty() {
        ui.label(egui::RichText::new("Identical to rCRS (no mutations).").weak());
        return;
    }
    let (mut hvr1, mut hvr2, mut coding) = (0usize, 0usize, 0usize);
    for v in variants {
        match v.region() {
            MtRegion::Hvr1 => hvr1 += 1,
            MtRegion::Hvr2 => hvr2 += 1,
            MtRegion::Coding => coding += 1,
        }
    }
    ui.label(
        egui::RichText::new(format!(
            "{} mutations vs rCRS  (HVR1 {hvr1} · HVR2 {hvr2} · Coding {coding})",
            variants.len()
        ))
        .weak(),
    );
    egui::ScrollArea::vertical()
        .max_height(300.0)
        .id_salt(("mt_mut", mtdna_id))
        .show(ui, |ui| {
            for region in [MtRegion::Hvr2, MtRegion::Coding, MtRegion::Hvr1] {
                let group: Vec<&MtVariant> = variants.iter().filter(|v| v.region() == region).collect();
                if group.is_empty() {
                    continue;
                }
                ui.add_space(4.0);
                ui.label(egui::RichText::new(format!("{} ({})", region.label(), group.len())).strong());
                egui::Grid::new(("mt_mut_grid", mtdna_id, region.label()))
                    .striped(true)
                    .num_columns(2)
                    .show(ui, |ui| {
                        ui.strong("Mutation");
                        ui.strong("Position");
                        ui.end_row();
                        for v in group {
                            ui.label(egui::RichText::new(v.notation()).monospace());
                            ui.label(v.position.to_string());
                            ui.end_row();
                        }
                    });
            }
        });
}
