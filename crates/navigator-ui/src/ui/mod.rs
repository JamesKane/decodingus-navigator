//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use crate::charts::{
    asset_status_line, coverage_histogram_chart, draw_ancestry_donut, draw_chromosome_painting, draw_color_donut,
    draw_composition_bar, draw_ibd_segments, draw_pca_scatter, draw_population_components, draw_roh,
    draw_variant_track, parse_hex_color, top_populations_for_side, TrackRegion, VariantMark,
};
use crate::widgets::{
    button_hit, capitalize_first, card, chip, combo, empty_state, fmt_depth, fmt_pct, fmt_reads, natural_cmp, opt,
    provider_abbrev, show_assignment, sortable_header, stat_card, variant_change, TableControls,
};
use eframe::egui;
use navigator_app::{
    AncestryResult, AppSettings, AuditEntry, BatchImportSummary, BuildNeed, CallState, ChatTurn, CompatibilityLevel,
    Consensus, Coverage, DenovoCall, DescentReport, DnaType, FtdnaGenealogy, FtdnaImportPlan, FtdnaResolution,
    HaploAssignment, HeteroplasmySite, IbdComparison, IbdSuggestion, IdentityVerification, LineageBrief, LineageKind,
    MatchKind, MatchStrength, MtRegion, MtVariant, NarratedBrief, PackStatus, PaintingResult, PrivateBucket,
    PrivateClass, ProjectBlockTree, ProjectOverview, ProjectSampleReport, ProjectStrChart, ReadMetrics, RefBuildStatus,
    SexInferenceResult, SignalKind, SnpEvidence, SourceType, StrConcordanceRow, SubjectAnalysisStatus, SubjectBrief,
    SvAnalysisResult, UiMode, VerificationStatus, YMatch, YProfile, YSignal, YState, YVariantStatus, YstrClustering,
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

use crate::worker::{self, Command, Event, NewBiosample, YMask};
use rowcache::{ReportRowCache, SubjectRowCache, VariantRows};

#[derive(Default)]
struct Forms {
    /// True when the inline "Add New Subject" form is open.
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
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Nav {
    Dashboard,
    Subjects,
    Projects,
    /// Federated IBD discovery and consent. It is top-level, and not a subject tab, because a
    /// matching conversation belongs to the *account*. Its key is our DID and the request URI of
    /// the broker, and it belongs to no one biosample. The subject comes into it only when it is
    /// time to exchange dosages.
    Matching,
    Community,
}

impl Nav {
    /// Stable persistence key (independent of display strings / i18n).
    fn as_key(self) -> &'static str {
        match self {
            Nav::Dashboard => "dashboard",
            Nav::Subjects => "subjects",
            Nav::Projects => "projects",
            Nav::Matching => "matching",
            Nav::Community => "community",
        }
    }
    fn from_key(s: &str) -> Option<Self> {
        match s {
            "dashboard" => Some(Nav::Dashboard),
            "subjects" => Some(Nav::Subjects),
            "projects" => Some(Nav::Projects),
            "matching" => Some(Nav::Matching),
            "community" => Some(Nav::Community),
            _ => None,
        }
    }
}

/// Sub-tabs of the Matching panel, following one conversation's life: a ranked candidate becomes a
/// request, a request becomes a result.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum MatchingTab {
    #[default]
    Suggestions,
    Requests,
    Results,
}
impl MatchingTab {
    const ALL: [(MatchingTab, &'static str); 3] = [
        (MatchingTab::Suggestions, "matching.tab.suggestions"),
        (MatchingTab::Requests, "matching.tab.requests"),
        (MatchingTab::Results, "matching.tab.results"),
    ];
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

/// Sub-tabs of the project detail panel: the member list, and the analysis report of each sample.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum ProjectTab {
    #[default]
    Members,
    Report,
    Ystr,
    Tree,
}
impl ProjectTab {
    const ALL: [(ProjectTab, &'static str); 4] = [
        (ProjectTab::Members, "project.tab.members"),
        (ProjectTab::Report, "project.tab.report"),
        (ProjectTab::Ystr, "project.tab.ystr"),
        (ProjectTab::Tree, "project.tab.tree"),
    ];
}

/// Sub-tabs of the Settings dialog. They group the preferences by what each one covers: general
/// appearance, server connection, ancestry-painter calibration, AI, reference genomes, one-off
/// tools, and read-only advanced info.
// `pub(crate)`, and its sibling tab enums are not: `open_settings` deep-links into a specific tab,
// so the type appears in a `pub(crate)` signature.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub(crate) enum SettingsTab {
    #[default]
    General,
    Connection,
    Ancestry,
    Ai,
    References,
    Tools,
    Advanced,
}
impl SettingsTab {
    const ALL: [(SettingsTab, &'static str); 7] = [
        (SettingsTab::General, "settings.tab.general"),
        (SettingsTab::Connection, "settings.tab.connection"),
        (SettingsTab::Ancestry, "settings.tab.ancestry"),
        (SettingsTab::Ai, "settings.tab.ai"),
        (SettingsTab::References, "settings.tab.references"),
        (SettingsTab::Tools, "settings.tab.tools"),
        (SettingsTab::Advanced, "settings.tab.advanced"),
    ];
}

/// Sub-tabs of the subject detail panel.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
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
    /// Stable persistence key (independent of the i18n display keys).
    fn as_key(self) -> &'static str {
        match self {
            DetailTab::Overview => "overview",
            DetailTab::YDna => "ydna",
            DetailTab::MtDna => "mtdna",
            DetailTab::Autosomal => "autosomal",
            DetailTab::Ancestry => "ancestry",
            DetailTab::Sources => "sources",
            DetailTab::IbdMatches => "ibd",
        }
    }
    fn from_key(s: &str) -> Option<Self> {
        match s {
            "overview" => Some(DetailTab::Overview),
            "ydna" => Some(DetailTab::YDna),
            "mtdna" => Some(DetailTab::MtDna),
            "autosomal" => Some(DetailTab::Autosomal),
            "ancestry" => Some(DetailTab::Ancestry),
            "sources" => Some(DetailTab::Sources),
            "ibd" => Some(DetailTab::IbdMatches),
            _ => None,
        }
    }

    /// `(tab, i18n key)` in display order. The DNA-type tabs (Y, mt, Autosomal, Ancestry) show the
    /// *consensus* of the subject over all sources. `Sources` is the hub for each sequencing
    /// result.
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

/// Sections of the Simple-mode subject view, in rail order.
///
/// Simple mode used to be one long vertical scroll, with every section at once. This splits it into
/// panels that a left rail reaches, and [`SimplePanel::Story`] is the first synopsis.
///
/// The order holds the narrative the view tells. First who you are, then the two lineages that
/// reach furthest back, then the autosomal ancestry (deep origins → recent populations). Then come
/// the relatives who are alive, and last the test the whole thing rests on.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
enum SimplePanel {
    #[default]
    Story,
    Paternal,
    Maternal,
    Ancestry,
    Relatives,
    Test,
}

impl SimplePanel {
    /// `(panel, icon, i18n label key)` in rail order.
    ///
    /// The **Proportional** font family of egui must cover every icon here (Ubuntu-Light,
    /// NotoEmoji, emoji-icon-font). That set is much narrower than it looks: `◆` U+25C6, `⚭`
    /// U+26AD, `✓` U+2713 and `🧬` U+1F9EC are all absent, and each one draws as a tofu box.
    /// `icon_glyphs_are_renderable` pins this down, so add an icon there when you add one here.
    const ALL: [(SimplePanel, &'static str, &'static str); 6] = [
        (SimplePanel::Story, "📖", "simple.panel.story"),
        (SimplePanel::Paternal, "♂", "simple.panel.paternal"),
        (SimplePanel::Maternal, "♀", "simple.panel.maternal"),
        (SimplePanel::Ancestry, "🌍", "simple.panel.ancestry"),
        (SimplePanel::Relatives, "👥", "simple.panel.relatives"),
        (SimplePanel::Test, "🔬", "simple.panel.test"),
    ];
}

/// Y-DNA sub-tabs: a compact haplogroup first view, the heavy SNP surface, and STR, which is
/// separate from SNP.
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

/// Which Y-STR report view the Y-DNA tab shows.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum StrReportView {
    /// FTDNA/YSEQ-style tier-grouped marker table.
    #[default]
    ByPanel,
    /// Flat, filterable marker table.
    AllMarkers,
    /// The cross-panel consensus value of each marker.
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

/// The Settings-dialog state the user can edit. It loads from `AppSettings`, and the reference rows
/// arrive through `Event::ReferenceSettings`.
#[derive(Clone)]
struct SettingsForm {
    appview_url: String,
    y_tree_provider: String, // "decodingus" | "ftdna"
    /// Prefer a trusted external caller (GATK4 GVCF / 1240K call set) over Navigator's own genotyping.
    prefer_external_calls: bool,
    tree_ttl_days: String,
    prompt_before_download: bool,
    /// UI scale (egui zoom factor); 1.0 = native.
    ui_scale: f32,
    /// Local-LLM (AI assistant) settings.
    llm_enabled: bool,
    llm_base_url: String,
    llm_model: String,
    llm_max_tokens: String,
    references: Vec<RefRow>,
    /// Chromosome-painter (copying-LAI) calibration knobs. Defaults mirror `CopyingLaiParams::default`.
    lai_recomb_per_cm: f64,
    lai_max_ref_haps: u32,
    lai_min_ancestry: f64,
    lai_switch_per_cm: f64,
    lai_min_segment_cm: f64,
    lai_size_normalize: f64,
    lai_mismatch: f64,
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
        // Unset knobs show the painter's own calibrated defaults, never a copy of them.
        let lai = navigator_app::lai_knob_defaults();
        SettingsForm {
            appview_url: s.appview_url.unwrap_or_default(),
            y_tree_provider: s.y_tree_provider.unwrap_or_else(|| "decodingus".to_string()),
            prefer_external_calls: s.prefer_external_calls.unwrap_or(true),
            tree_ttl_days: s
                .tree_ttl_days
                .map(|d| d.to_string())
                .unwrap_or_else(|| "7".to_string()),
            prompt_before_download: s.prompt_before_download.unwrap_or(true),
            ui_scale: s.ui_scale.unwrap_or(1.0),
            llm_enabled: s.llm_enabled.unwrap_or(false),
            llm_base_url: s
                .llm_base_url
                .unwrap_or_else(|| navigator_app::llm::DEFAULT_LLM_BASE_URL.to_string()),
            llm_model: s.llm_model.unwrap_or_default(),
            llm_max_tokens: s
                .llm_max_tokens
                .unwrap_or(navigator_app::llm::DEFAULT_LLM_MAX_TOKENS)
                .to_string(),
            references: Vec::new(),
            lai_recomb_per_cm: s.lai_recomb_per_cm.unwrap_or(lai.recomb_per_cm),
            lai_max_ref_haps: s.lai_max_ref_haps.unwrap_or(lai.max_ref_haps),
            lai_min_ancestry: s.lai_min_ancestry.unwrap_or(lai.min_ancestry),
            lai_switch_per_cm: s.lai_switch_per_cm.unwrap_or(lai.switch_per_cm),
            lai_min_segment_cm: s.lai_min_segment_cm.unwrap_or(lai.min_segment_cm),
            lai_size_normalize: s.lai_size_normalize.unwrap_or(lai.size_normalize),
            lai_mismatch: s.lai_mismatch.unwrap_or(lai.mismatch),
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

/// A one-shot auto UI-scale probe: the "behave like a native app" default. On the first frame that
/// knows the monitor size, it derives a zoom when the OS reports a ~1.0 scale factor on a clearly
/// high-resolution panel. Native 4K is the example, where macOS itself does not scale up. The
/// native scaling of egui already controls a Retina or scaled display, which has a native ppp above
/// 1, so this leaves that at 1.0. It does nothing at all when a manual scale is on disk, because
/// `probed` then starts `true`. The result fills the Settings slider, and nothing persists it until
/// the user saves, so it probes again on each launch.
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

/// The state of a full analysis in progress. It drives the modal dialog.
#[derive(Clone)]
struct AnalysisModal {
    step: usize,
    total: usize,
    label: String,
    detail: String,
    fraction: f32,
    /// The egui time in seconds when this step began, for the elapsed-time display.
    started: f64,
}

/// A realignment in progress, or the result of the last one.
///
/// This is card state on purpose, and not a modal. A realignment runs for hours. A dialog that owns
/// the screen for that long stops the user from a workspace they can perfectly well keep open. The
/// card sits with the alignment it belongs to, and it updates in place.
#[derive(Clone)]
struct RealignState {
    /// The source alignment this job belongs to. The card of another alignment ignores it.
    alignment_id: i64,
    /// The subject it belongs to. The card of Simple mode is about a *person*, and not an
    /// alignment, so it must match on this. With only the alignment id to go on, a page open on
    /// subject A during a job on subject B told A the wrong thing. It reported a job at work on
    /// their genome.
    biosample_guid: Option<SampleGuid>,
    step: usize,
    total: usize,
    label: String,
    detail: String,
    /// Set once the job ends; `None` while it runs.
    finished: Option<RealignFinished>,
}

/// How a realignment ended, for the card to report.
#[derive(Clone)]
enum RealignFinished {
    /// Registered as `new_alignment_id`, with a one-line summary.
    Done {
        new_alignment_id: i64,
        summary: String,
    },
    Cancelled,
    Failed(String),
}

/// A copy of a project the user can edit. It drives the project Edit modal, and `Some` shows the
/// dialog.
#[derive(Clone)]
struct EditProject {
    id: i64,
    name: String,
    description: String,
    administrator: String,
}

/// A copy of a sequence run the user can edit. It drives the run Edit modal, and `Some` shows the
/// dialog. The read-metric columns take no edit here, so this does not carry them.
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

/// This drives the destructive merge-sequence-runs modal, and `Some` shows it. `secondary` is the
/// run this empties and then deletes. `primary` is the chosen merge target, and the default of its
/// picker.
#[derive(Clone)]
struct MergeRuns {
    guid: SampleGuid,
    secondary: i64,
    primary: Option<i64>,
}

/// A copy of an alignment the user can edit. It drives the alignment Edit modal, and `Some` shows
/// the dialog.
#[derive(Clone)]
struct EditAlignment {
    id: i64,
    run_id: i64,
    reference_build: String,
    aligner: String,
    variant_caller: String,
}

/// A copy of a subject the user can edit. It drives the Edit modal, and `Some` shows the dialog.
#[derive(Clone)]
struct EditSubject {
    guid: SampleGuid,
    donor_identifier: String,
    sample_accession: String,
    description: String,
    center_name: String,
    sex: String,
}

/// Drives the "add vendor id" modal (Some ⇒ shown). A new kit association for `guid`; `source` is
/// a well-known vendor (FTDNA/YSEQ/…) or free text, `external_id` is the kit number.
#[derive(Clone)]
struct EditKit {
    guid: SampleGuid,
    source: String,
    external_id: String,
}

/// A copy of an MDKA (most distant known ancestor) the user can edit. It drives the MDKA edit
/// modal, and `Some` shows it. Every field is a string, so that the user can edit it, and a save
/// parses the years and the coordinates. A blank field clears its column. `lineage` does not change
/// after the modal opens (Y, Mt or Auto). The key of the upsert is `(guid, lineage)`.
#[derive(Clone)]
struct EditMdka {
    guid: SampleGuid,
    lineage: String,
    ancestor_name: String,
    birth_year: String,
    death_year: String,
    origin_place: String,
    origin_country: String,
    latitude: String,
    longitude: String,
    notes: String,
}

/// A data-source row that waits for a delete confirmation. The confirm dialog shows the `label`,
/// and the variant carries the ids the worker command needs, plus the parent id to refresh.
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

/// The first inner size of the window, in egui points, on a first run, and when the app remembers
/// nothing. The new layout needs room, because it holds the subjects table, the detail panel and
/// the action bar. The eframe default is far too small.
pub(crate) const DEFAULT_WINDOW: [f32; 2] = [1360.0, 900.0];
/// The minimum inner size of the window: a sane floor, so that the layout never collapses.
pub(crate) const MIN_WINDOW: [f32; 2] = [1024.0, 680.0];

/// Fit a wanted window size, in egui points, to the monitor. It leaves a margin for the menu bar,
/// the dock and the taskbar, and it never goes below `min`. A remembered size that is too large,
/// for example from a bigger display, shrinks to fit. A size that already fits comes back
/// unchanged. It is pure, so a unit test covers the fit.
pub(crate) fn fit_window_to_monitor(desired: [f32; 2], monitor: [f32; 2], min: [f32; 2]) -> [f32; 2] {
    let max_w = (monitor[0] * 0.98).max(min[0]);
    let max_h = (monitor[1] * 0.94).max(min[1]);
    [desired[0].clamp(min[0], max_w), desired[1].clamp(min[1], max_h)]
}

pub struct NavigatorApp {
    tx: UnboundedSender<Command>,
    rx: Receiver<Event>,
    /// The progress of a full analysis in progress. `Some` shows the modal dialog.
    analysis: Option<AnalysisModal>,
    /// The realignment that runs now, or the last one that ended. See [`RealignState`].
    realign: Option<RealignState>,
    /// The realignment confirmation of Simple mode, while it waits. It is the same [`RealignOffer`]
    /// the brief gave, and not a tuple that writes out its two fields again.
    ///
    /// Simple mode gets a confirmation step, and Advanced does not. That difference is deliberate.
    /// The Advanced card sits among alignment internals, and it states its cost in a paragraph its
    /// reader can weigh. The reader of Simple mode has seen a story about their ancestors. One
    /// button, read wrong, must not commit the machine to four hours and 276 GB.
    simple_realign_confirm: Option<navigator_domain::brief::RealignOffer>,
    /// The code sets this the moment the user clicks Cancel, and clears it when the run ends.
    ///
    /// Cancellation is cooperative. The walkers stop at their next check, so there is always a gap
    /// between the click and the end of the run. Without this the UI gave no acknowledgement at
    /// all: the spinner, the timer and the progress bar all continued, and the button stayed live.
    /// That is why the button read as broken, and not as busy.
    cancelling: bool,
    /// The subject the user edits. `Some` shows the Edit modal.
    edit_subject: Option<EditSubject>,
    /// The vendor-id (kit) association the user adds. `Some` shows the add-kit modal.
    edit_kit: Option<EditKit>,
    /// The MDKA the user edits or adds. `Some` shows the MDKA modal.
    edit_mdka: Option<EditMdka>,
    /// The subject that waits for a delete confirmation. `Some` shows the confirm dialog.
    confirm_delete: Option<SampleGuid>,
    /// The subject that waits for a "clear all data" confirmation. `Some` shows the confirm
    /// dialog.
    confirm_clear: Option<SampleGuid>,
    /// The subject that waits for a "reset haplogroup placement" confirmation. `Some` shows the
    /// confirm dialog.
    confirm_reset_haplo: Option<SampleGuid>,
    /// The last batch-import summary, shown in a modal until dismissed.
    batch_import: Option<BatchImportSummary>,
    /// Y-STR-from-sequence concordance for the selected subject: `(guid, source alignment, rows)`.
    str_concordance: Option<(SampleGuid, i64, Vec<StrConcordanceRow>)>,
    /// True while a Y-STR-from-sequence call is in progress (the heavy first pass).
    str_running: bool,
    /// Cross-subject Y matches for the selected subject: `(guid, ranked matches)`. Gap §2.
    y_matches: Option<(SampleGuid, Vec<YMatch>)>,
    /// True while a Y-match search is in progress.
    y_matches_running: bool,
    /// Project filter for the Y-match search (None ⇒ whole workspace).
    y_match_project: Option<i64>,
    /// Text filter over the Y-match table (by donor / haplogroup).
    y_match_query: String,
    /// The data-source row that waits for a delete confirmation. `Some` shows the confirm dialog.
    confirm_data_delete: Option<DataDelete>,
    /// The subject the user assigns to a project: (subject, selected project or None). `Some`
    /// shows the picker.
    assign_project: Option<(SampleGuid, Option<i64>)>,
    /// The project the user edits. `Some` shows the project Edit modal.
    edit_project: Option<EditProject>,
    /// The project that waits for a delete confirmation: (id, name). `Some` shows the confirm
    /// dialog.
    confirm_delete_project: Option<(i64, String)>,
    /// The sequence run the user edits. `Some` shows the run Edit modal.
    edit_run: Option<EditRun>,
    merge_runs: Option<MergeRuns>,
    /// Whether the read-only Y-profile source-audit modal is open (reads the cached `y_profile`).
    audit_y_profile: bool,
    /// The alignment the user edits. `Some` shows the alignment Edit modal.
    edit_alignment: Option<EditAlignment>,
    /// Current frame's egui time (seconds), captured at the top of `update`.
    frame_time: f64,
    /// Window-size persistence, all in egui points. It holds the current inner size, the last size
    /// that went to [`AppSettings`], and when it last changed (`frame_time` seconds), for a
    /// debounced save.
    /// `window_restored` guards the one-time restore and fit to the screen. `startup_frames` lets
    /// that wait until the UI scale (zoom) settles, so that the sizes are in a stable unit.
    window_size: Option<[f32; 2]>,
    saved_window_size: Option<[f32; 2]>,
    window_size_changed_at: f64,
    window_restored: bool,
    startup_frames: u32,
    /// Focused subject remembered from the last session (GUID string), applied once the subject list
    /// loads (`.take()`n so it restores only once).
    pending_restore_subject: Option<String>,
    /// The signature of the last navigation state that went to [`AppSettings`] (view | subject |
    /// tab), so that a save fires only on a real change.
    saved_ui_sig: Option<String>,
    /// Selected primary navigation tab.
    nav: Nav,
    /// Interface mode: Simple (casual single-person briefs) vs. Advanced (full power-user UI).
    ui_mode: UiMode,
    /// True when something pinned the mode explicitly: the environment, the settings, or a user
    /// toggle. With `false`, the first-run workspace heuristic can still change it as data loads.
    ui_mode_pinned: bool,
    /// Precomputed Simple-mode brief for the selected subject `(guid, brief)`; `None` until built.
    subject_brief: Option<(SampleGuid, SubjectBrief)>,
    /// True while a Subject Brief build is in progress.
    subject_brief_loading: bool,
    /// Free-text filter over the Simple-mode "My DNA" subject selector (matches the donor name).
    simple_subject_filter: String,
    /// Which Simple-mode panel the left rail has open. It goes back to the first synopsis on every
    /// subject switch. The view of each person starts from their story, and not from the panel that
    /// held the last person.
    simple_panel: SimplePanel,
    /// True when the local-LLM "AI assistant" is on. It comes from the settings cache, and it gates
    /// "Polish with AI".
    ai_enabled: bool,
    /// AI-assisted narration of the selected subject's brief `(guid, narration)`; `None` until run.
    brief_narration: Option<(SampleGuid, NarratedBrief)>,
    /// The live narration text, which grows while the stream runs, as `(guid, text)`. It clears
    /// when the final narration arrives.
    narration_stream: Option<(SampleGuid, String)>,
    /// True while a brief narration request is in progress.
    narrating: bool,
    /// "Ask my results" chat history for the selected subject (cleared on subject switch).
    chat_history: Vec<ChatTurn>,
    /// The chat input box, and whether an answer is in progress.
    chat_input: String,
    chat_pending: bool,
    /// The "Explain this" (M5) state of each tab, for the selected subject. A subject switch clears
    /// all of it. It holds the final explanations, with `(guid, signal)` as the key. It also holds
    /// the live stream buffer for the one in progress, and which `(guid, signal)` the model narrates
    /// now. `None` means idle. Only one runs at a time, because they share the one worker.
    signal_narration: Vec<(SampleGuid, SignalKind, NarratedBrief)>,
    signal_stream: Option<(SampleGuid, SignalKind, String)>,
    signal_narrating: Option<(SampleGuid, SignalKind)>,
    /// Selected subject-detail sub-tab.
    detail_tab: DetailTab,
    /// Active UI language.
    lang: crate::i18n::Lang,
    /// Dark (default) vs light theme.
    dark_mode: bool,
    /// True after the one-shot auto-UI-scale probe runs. It does not run when a manual scale is on
    /// disk.
    scale_probed: bool,
    /// Settings dialog open + its editable form.
    show_settings: bool,
    settings_form: SettingsForm,
    /// Active Settings dialog sub-tab.
    settings_tab: SettingsTab,
    /// Local-LLM "Test connection" state (Settings → AI assistant): discovered models, an in-flight
    /// flag, and the last plain-language status/error line.
    llm_models: Vec<String>,
    llm_testing: bool,
    llm_test_msg: Option<String>,
    /// The sort state, and the inline filter state of each column, for the subjects table.
    subjects_table_ctl: TableControls,
    /// This grows by one for each worker [`Event`] the code applies. It is the invalidation signal
    /// for the view caches below. Every field they derive from comes from [`Self::drain_events`],
    /// so a change of epoch is the one thing they all have to watch.
    ///
    /// A bump on *every* event invalidates too much, because a progress tick rebuilds a table it
    /// did not touch. That is better than a risk of stale rows. One extra rebuild costs one frame,
    /// and one that never happens shows wrong data until the user clicks something.
    data_epoch: u64,
    /// Derived display rows for the subjects table (see [`Self::subject_rows`]).
    subject_rows: SubjectRowCache,
    /// Derived cell text + row order for the project Report table (see [`Self::report_rows`]).
    report_rows: ReportRowCache,
    /// Matching-row indices for the Y / mtDNA / autosomal variant tables.
    y_profile_rows: VariantRows,
    mt_profile_rows: VariantRows,
    auto_profile_rows: VariantRows,
    /// Collapse the subjects side panel to a thin strip so the detail panel (charts/tables)
    /// gets the full width.
    subjects_collapsed: bool,
    /// Collapse the projects side panel to a thin strip, and give the detail panel the full width.
    projects_collapsed: bool,
    overview: Vec<ProjectOverview>,
    selected_project: Option<i64>,
    /// When a subject came from the report row of a project, the project id to return to. It drives
    /// the "back to project" button in the detail header. Any other navigation clears it.
    return_to_project: Option<i64>,
    /// The coverage and haplogroup report rows of each sample, for the selected project.
    project_report: Vec<ProjectSampleReport>,
    /// Precomputed Y-STR overview (FTDNA-style chart) for the selected project; `None` until the
    /// background build returns. A boolean tracks the in-flight build so the UI can show a spinner.
    project_str_chart: Option<ProjectStrChart>,
    project_str_loading: bool,
    /// The cohort Y **block tree** for the selected project. It is `None` until the background build
    /// returns. It loads **lazily, on the first view of the Tree tab**, and not on a project select
    /// as the STR chart does. A build reads and parses a multi-MB haplotree, and that is too much to
    /// spend on a tab nobody opened.
    project_blocktree: Option<ProjectBlockTree>,
    project_blocktree_loading: bool,
    /// Blocks (by node id) the user expanded to reveal their equivalent SNPs and full member list.
    /// Zoom factor for the block-tree canvas (1.0 = natural size).
    blocktree_zoom: f32,
    /// The candidate branch open for review, as its synthetic node id, when there is one. A
    /// candidate is an inference, so it gets a surface that shows the evidence, and does not ask for
    /// trust.
    blocktree_review: Option<i64>,
    /// The block whose member roster sits beside the tree. The Big Tree keeps the men in a table,
    /// and not in the diagram. This is that table, held to what the user clicked.
    blocktree_selected: Option<i64>,
    /// Centre the canvas on the root again on the next frame. The code sets it when a tree first
    /// arrives. The view then does not open on the empty left margin of a canvas far wider than any
    /// viewport.
    blocktree_recentre: bool,
    samples: Vec<Biosample>,
    /// Every biosample (the project-independent subjects list).
    all_biosamples: Vec<Biosample>,
    /// The terminal Y and mt haplogroups of each subject, for the list columns (`guid → (Y, mt)`).
    haplo_summary: std::collections::HashMap<SampleGuid, (Option<String>, Option<String>)>,
    /// The analysis status of each subject (Pending or Complete), for the Status column of the
    /// subjects list. A subject the map does not hold has no alignment to analyze, and it appears
    /// with no status.
    subject_status: std::collections::HashMap<SampleGuid, SubjectAnalysisStatus>,
    selected_sample: Option<SampleGuid>,
    runs: Vec<SequenceRun>,
    /// Donor-level haplogroup consensus for the selected subject (Y, mtDNA).
    consensus_y: Option<Consensus>,
    consensus_mt: Option<Consensus>,
    /// Descent reports in the YFull style, for the selected subject. Each `DnaType` loads lazily,
    /// and the cache holds `Some(report)`, or `None` when placement reached it and it is empty. A
    /// result that built one time never loads again. This also holds the (guid, dna) pairs that
    /// load now. A subject switch clears all of it.
    descent_reports: Vec<(SampleGuid, DnaType, Option<DescentReport>)>,
    descent_loading: Vec<(SampleGuid, DnaType)>,
    /// The branch report of each marker, for the selected subject. It holds the node-name text
    /// inputs (Y and mt), and the last report the cache took, as `Some(report)`, or `None` when
    /// there is no alignment. It also holds the (guid, dna) pairs that load now. A subject switch
    /// clears it. A node starts it, from a Load button, and it is not lazy.
    branch_node_y: String,
    branch_node_mt: String,
    branch_reports: Vec<(SampleGuid, DnaType, Option<navigator_app::BranchReport>)>,
    branch_loading: Vec<(SampleGuid, DnaType)>,
    /// Reconciliation audit log for the selected subject (Y, mtDNA).
    audit_y: Vec<AuditEntry>,
    audit_mt: Vec<AuditEntry>,
    /// Last mtDNA heteroplasmy scan: (alignment id, sites).
    heteroplasmy: Option<(i64, Vec<HeteroplasmySite>)>,
    /// STR profiles for the selected subject.
    str_profiles: Vec<StrProfile>,
    /// The view state of the Y-STR report: which view, which provider when there is more than one,
    /// and the marker filter.
    str_report_view: StrReportView,
    str_provider: Option<String>,
    str_marker_filter: String,
    /// SNP variant sets for the selected subject.
    variant_sets: Vec<VariantSet>,
    /// Chip/array profiles for the selected subject.
    chip_profiles: Vec<ChipProfile>,
    /// mtDNA sequences for the selected subject.
    mtdna_sequences: Vec<MtdnaSequence>,
    /// The mutation list against rCRS of each mtDNA sequence id. It loads on demand.
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
    /// True while a build of the haplogroup report runs.
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
    /// Reference PC1/PC2 centroids for the PCA scatter, keyed by alignment_id (lazy-loaded).
    pca_reference: Option<(i64, PcaCentroids)>,
    /// The PCA-reference key we already dispatched a load for. It stops a second dispatch on every
    /// frame.
    pca_reference_attempted: Option<i64>,
    /// Donor-level private-Y union across the subject's sources.
    donor_private_y: Option<PrivateBucket>,
    /// The selected subject's multi-source Y-variant profile.
    y_profile: Option<YProfile>,
    /// Y-variant profile status filter (None = all).
    y_profile_filter: Option<YVariantStatus>,
    /// Text search over the variant and SNP tables, by SNP name or site, one for each table.
    y_profile_query: String,
    mt_profile_query: String,
    auto_profile_query: String,
    private_y_query: String,
    str_seq_query: String,
    /// The catalogued Y-SNP names at variant positions (`position → name`). They annotate the
    /// position-only and novel calls of the two Y-SNP tables. The Y-SNP dictionary resolves them one
    /// time for each subject.
    y_snp_names: std::collections::HashMap<i64, String>,
    /// True after we dispatch the Y-SNP-name resolution for the current subject. It stops a second
    /// dispatch.
    y_snp_names_requested: bool,
    /// True while a build of the Y-variant profile runs. That build has a high cost.
    y_profile_loading: bool,
    /// The selected subject's multi-source mtDNA consensus profile.
    mt_profile: Option<navigator_app::ConsensusProfile>,
    /// mtDNA consensus-profile status filter (None = all).
    mt_profile_filter: Option<YVariantStatus>,
    /// True while a build of the mtDNA consensus profile runs. That build has a high cost.
    mt_profile_loading: bool,
    /// The selected subject's multi-source autosomal (diploid 0/1/2) consensus profile.
    auto_profile: Option<navigator_app::DiploidProfile>,
    /// Autosomal consensus status filter (None = all).
    auto_profile_filter: Option<YVariantStatus>,
    /// True while a build of the autosomal consensus profile runs. That build has a high cost.
    auto_profile_loading: bool,
    /// True while the donor ancestry estimate from the consensus is in progress.
    estimating_donor_ancestry: bool,
    /// True while the heavy deep (ancient) ancestry estimate is in progress.
    estimating_deep_ancestry: bool,
    /// Local-ancestry painting: (alignment id, result). The result holds the segments of each side,
    /// and the side labels. `painting_running` is true while the genotyping runs.
    painting: Option<(i64, PaintingResult)>,
    painting_running: bool,
    /// Runs-of-homozygosity result for the selected subject. `roh_running` while the HMM computes.
    roh: Option<navigator_app::RohResult>,
    roh_running: bool,
    /// Archaic (Neanderthal / Denisovan) Tier-A marker count for the selected subject.
    archaic: Option<navigator_app::ArchaicMarkerResult>,
    archaic_running: bool,
    /// Tier B archaic segments for the selected subject.
    archaic_segments: Option<navigator_app::ArchaicSegmentResult>,
    archaic_segments_running: bool,
    /// Last private Y bucket: (alignment id, bucket).
    private_y: Option<(i64, PrivateBucket)>,
    finding_private_y: bool,
    /// Callable-region BED (external mask), reused across private-Y runs.
    y_mask_path: Option<PathBuf>,
    /// Use the callable-Y BED of the sample itself (self-referential), and not an external mask.
    y_self_mask: bool,
    selected_run: Option<i64>,
    alignments: Vec<Alignment>,
    selected_alignment: Option<i64>,
    /// An alignment to auto-select once its run's alignments load (subject-centric default).
    pending_alignment: Option<i64>,
    coverage: Option<Coverage>,
    /// The cached coverage of each alignment, for the Data Sources rows of the selected run. Each
    /// row then shows coverage and callable before anybody selects that alignment. The key is the
    /// alignment id.
    coverage_by_aln: std::collections::HashMap<i64, Coverage>,
    /// Genome-region metadata (the cytoband ideogram) for the build of the selected alignment, as
    /// `(alignment_id, regions)`. It loads lazily when the user opens the Ideogram tab.
    genome_regions: Option<(i64, std::sync::Arc<navigator_app::GenomeRegions>)>,
    /// True while the cytoBand read is in progress.
    loading_regions: bool,
    /// The alignment we already started a region load for, or completed one for. It stops a second
    /// read on every frame, and after a failure too.
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
    all_alignments: Vec<Alignment>,
    /// `(project_id, eligible alignment ids)` for the project realignment card. The app answers it,
    /// because the rule belongs to the app, and the number must cover the project, and not the whole
    /// workspace.
    project_realignable: Option<(i64, Vec<i64>)>,
    /// The project that already has a count request. The card then does not ask again on every
    /// frame, and it does not ask again after a failure. It has the same shape as
    /// `regions_attempted`.
    project_realignable_asked: Option<i64>,
    /// Chip-compatible IBD compare: the two picked sources (each a WGS alignment or an imported chip).
    ibd_src_a: Option<navigator_app::IbdSource>,
    ibd_src_b: Option<navigator_app::IbdSource>,
    /// Subject-level (consensus) IBD compare: the other subject picked for comparison.
    ibd_other_subject: Option<SampleGuid>,
    /// True when the filter and the list of the consensus-compare subject picker are open.
    ibd_other_picking: bool,
    /// Filter text for that picker.
    ibd_other_filter: String,
    ibd_result: Option<IbdComparison>,
    running_ibd: bool,
    /// Identity-verification result for the current IBD pair.
    identity: Option<IdentityVerification>,
    /// Federated IBD: pseudonymous match suggestions fetched from the AppView.
    ibd_suggestions: Vec<IbdSuggestion>,
    /// True while a read of the suggestions is in progress. It drives the spinner.
    loading_ibd_suggestions: bool,
    /// The introduction status of each candidate, with `suggested_sample_guid` as the key, for
    /// example "PENDING".
    ibd_intros: std::collections::HashMap<String, String>,
    /// The selected subject's persisted IBD exchange results.
    exchange_results: Vec<navigator_app::StoredIbdExchange>,
    /// True while an inbox refresh, a consent, or an exchange run is in progress.
    exchange_busy: bool,
    /// The matching ledger: every conversation, whatever its stage. It replaces the view of the
    /// same data on each card. `ibd_intros` does not survive a restart, and this does, because the
    /// app persists it.
    matching: Vec<navigator_app::MatchingEntry>,
    /// Which stage of the Matching panel is on the screen.
    matching_tab: MatchingTab,
    /// The candidates the user dismissed in this session. They go at once, and the code does not
    /// wait for a second read. The AppView keeps the authoritative dismissal.
    dismissed_candidates: std::collections::HashSet<String>,
    /// The local subject whose dosages an exchange will use. Defaults to the selected subject.
    matching_subject: Option<SampleGuid>,
    /// True when the filter and the list of the subject picker are open. It is a reveal, and not a
    /// dropdown, so a workspace of 10k subjects costs only the rows on the screen.
    matching_subject_picking: bool,
    /// Filter text for that picker.
    matching_subject_filter: String,
    /// The request URI whose consent decision waits for a confirmation, with what we know of the
    /// request.
    consent_prompt: Option<navigator_app::MatchingEntry>,
    /// Signed-in account DID, or `None`. Gates the "Publish" actions.
    account: Option<String>,
    /// Whether the last PDS write reached the server (offline indicator).
    online: bool,
    /// The outbox rows that still wait for a push to succeed (the "N pending" sync indicator).
    sync_pending: i64,
    /// True while a PULL reconcile is in progress.
    pulling: bool,
    logging_in: bool,
    publishing: bool,
    /// True while a batch project-directory import is in progress. It disables the button.
    importing: bool,
    /// The dir to import a second time, after the necessary references arrive.
    pending_import_dir: Option<PathBuf>,
    /// The reference builds an import waits on. Prompt the user to download them.
    reference_needs: Vec<BuildNeed>,
    /// In-flight reference download: (build, received, total).
    reference_progress: Option<(String, u64, Option<u64>)>,
    /// In-flight coordinate-index build: (file label, done bytes, total bytes / `None` = indeterminate).
    index_progress: Option<(String, u64, Option<u64>)>,
    /// A newer installer is available (drives the update-notification modal). Set once by the
    /// startup `CheckForUpdate`; cleared when the user dismisses it.
    update_info: Option<navigator_app::UpdateInfo>,
    /// True while an analyze pass over the whole project runs. It disables the analyze button of
    /// the report.
    analyzing: bool,
    /// Deep-analyze progress from the stream: `(done, total, current_sample, fraction)` while the
    /// pass runs.
    deep_progress: Option<(usize, usize, String, f32)>,
    /// Workspace-chore survey (Dashboard → Maintenance). `None` until the user asks for it: two of
    /// the three chores cost real work to measure, one a multi-MB tree fetch.
    maintenance: Option<Vec<navigator_app::ChoreSurvey>>,
    /// True while the survey is in progress, so that the button can say so.
    maintenance_surveying: bool,
    /// The chore that runs now, with its progress line.
    chore_running: Option<(navigator_app::Chore, usize, usize, String, f32)>,
    /// What the last chore did, kept on screen so a finished job is not just a vanished bar.
    chore_last: Option<(navigator_app::Chore, navigator_app::ChoreOutcome)>,
    /// The dry-run FTDNA import plan under review. It drives the review modal.
    ftdna_plan: Option<FtdnaImportPlan>,
    /// The resolution the admin chose for each kit, for the fuzzy rows in [`Self::ftdna_plan`].
    ftdna_resolutions: std::collections::BTreeMap<String, FtdnaResolution>,
    /// The selected subject's imported genealogy (vendor ids + FTDNA member + MDKA), for the
    /// Overview card. `(guid, data)` so a stale bundle from a prior subject is not shown.
    genealogy: Option<(SampleGuid, FtdnaGenealogy)>,
    /// The current project's Y-STR clustering, keyed by project id (so a stale one is not shown).
    project_clustering: Option<(i64, YstrClustering)>,
    /// True while the project Y-STR clustering runs.
    clustering_running: bool,
    /// Active project detail sub-tab (Members vs Report).
    project_tab: ProjectTab,
    /// Filter for the project Members list (kit / name / branch substring).
    member_filter: String,
    /// The sort state, and the inline filter state of each column, for the project Report table.
    report_table_ctl: TableControls,
    // ---- Community (social) ------------------------------------------------
    /// Active Community sub-tab (Support / Feed / Notifications).
    community_tab: CommunityTab,
    /// The signed-in account's support threads.
    support_threads: Vec<navigator_app::SocialThreadSummary>,
    /// The `(conversation_id, messages)` of the open thread. It is `None` while the list is on the
    /// screen.
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
    /// Opt-in: also publish the next community post to the PDS that signed in, as a federated
    /// `feed.post` record (roadmap 3b). It is off by default, because a write to your own repo is an
    /// explicit, portable public act.
    feed_publish_pds: bool,
    forms: Forms,
    status: String,
    /// The file-level diagnosis behind the last failed alignment command, when there was one.
    /// [`Event::Diagnosed`] sets it, a modal shows it, and it clears when the user dismisses it or
    /// starts something new. `Some` is what makes the "Details" control of the status bar appear.
    /// An error with no diagnosis must not offer one.
    diagnosis: Option<String>,
    /// True when the diagnosis modal is open. It is separate from [`Self::diagnosis`], so that a
    /// dismiss of the modal keeps the report reachable from the status bar, and does not destroy
    /// it.
    show_diagnosis: bool,
}

/// Sentinel option in the chip-provider dropdown that means "let the parser guess".
const AUTO_DETECT: &str = "(auto-detect)";

/// The workbench accent (primary buttons, selection, active tabs).
pub(crate) const ACCENT: egui::Color32 = egui::Color32::from_rgb(45, 125, 246);
/// Destructive-action red (Delete buttons, the confirm modals, unconfirmed rows).
const DANGER: egui::Color32 = egui::Color32::from_rgb(220, 60, 60);

/// How many rows the heavy variant and site tables show before they scroll inside themselves. Those
/// tables are the Y, mt and autosomal consensus profiles, private-Y, and the de-novo SNPs. On a WGS
/// they run thousands of rows. A limit keeps the detail page easy to move around, with no endless
/// scroll. But the pane must be tall enough to be useful: the user wants 20-30 rows, and not a slot
/// of 3.
const PROFILE_TABLE_ROWS: usize = 26;

/// An explicit height for a profile table that scrolls. It shows as many as
/// [`PROFILE_TABLE_ROWS`] rows of the given `count`, then it scrolls. With fewer rows it takes the
/// size of its content, so there is no empty pane. The code calculates it, and does not use a fixed
/// constant. A nested `ScrollArea` with vertical `auto_shrink` collapses to a few rows inside the
/// page scroll. Pair this with `auto_shrink([false, false])`.
fn profile_pane_height(ui: &egui::Ui, count: usize) -> f32 {
    let row_h = ui.text_style_height(&egui::TextStyle::Body) + ui.spacing().item_spacing.y + 3.0;
    let rows = count.clamp(1, PROFILE_TABLE_ROWS) as f32;
    row_h * (rows + 1.0) // +1 for the header row
}

/// Apply the Decoding-Us workbench look: a dark or light palette with the accent blue, rounded
/// widgets, and more room between them. That is the visual base that closes most of the gap to the
/// Scala Workbench. A theme toggle applies it again.
fn apply_theme(ctx: &egui::Context, dark: bool) {
    use egui::{Color32, Rounding, Stroke};
    // Pin the preference, so that egui stops following the theme of the OS. If not,
    // `theme_preference` defaults to `System`, and the host writes over our styled visuals. A light
    // macOS would then show light, even with Dark selected in Settings.
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
        v.window_stroke = Stroke::new(1.0_f32, Color32::from_gray(48));
    }
    v.hyperlink_color = ACCENT;
    v.selection.bg_fill = ACCENT.gamma_multiply(0.55);
    v.selection.stroke = Stroke::new(1.0_f32, ACCENT);

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

mod blocktree;
mod branch;
mod central;
mod chrome;
mod community;
mod descent;
mod detail;
mod events;
mod ibd;
mod matching;
mod modals;
mod rowcache;
mod simple;
mod sources;

impl NavigatorApp {
    pub fn new(cc: &eframe::CreationContext<'_>, db_path: PathBuf) -> Self {
        let ctx = cc.egui_ctx.clone();
        let (tx, rx) = worker::spawn(db_path, move || ctx.request_repaint());
        let _ = tx.send(Command::LoadOverview);
        let _ = tx.send(Command::LoadAllBiosamples);
        let _ = tx.send(Command::LoadAllAlignments);
        let _ = tx.send(Command::AuthStatus);
        let _ = tx.send(Command::SyncStatus);
        let _ = tx.send(Command::BackfillLabs); // resolve labs for runs imported before D8 landed
        let _ = tx.send(Command::VerifySourceFiles); // flag any imported file that moved/disappeared
        let _ = tx.send(Command::LoadAssetStatus); // ancestry/IBD "data sources" line

        // Check for a newer installer at startup, unless the user opted out. It is not fatal: a
        // check that fails only logs to the status line, and the app never updates itself.
        //
        // One read of settings.json for the whole constructor. The code used to load it six
        // separate times.
        let settings = AppSettings::load();
        if settings.check_for_updates != Some(false) {
            let _ = tx.send(Command::CheckForUpdate);
        }
        // Persisted theme wins; default dark. (Must match `dark_mode` below.)
        let dark = !matches!(settings.theme.as_deref(), Some("light"));
        apply_theme(&cc.egui_ctx, dark);
        // The UI scale (the egui zoom) from disk. It fixes tiny text on a native-4K display that
        // the OS reports at scale factor 1.0. The keyboard zoom of egui (Cmd +/-/0) also works, and
        // nothing persists that.
        cc.egui_ctx.set_zoom_factor(resolved_ui_scale());
        // Restore the last navigation position: the view, the focused subject, and the detail tab.
        // The subject applies after the list loads (see the `AllBiosamples` handler). The nav and
        // the tab apply at once, and `normalize_for_mode` then reconciles the nav to the interface
        // mode. Fill `saved_ui_sig` with the restored intent, so that a restore that matches does
        // not start a save nobody needs.
        let restore = &settings;
        let restored_nav = restore
            .last_nav
            .as_deref()
            .and_then(Nav::from_key)
            .unwrap_or(Nav::Subjects);
        let restored_tab = restore
            .last_detail_tab
            .as_deref()
            .and_then(DetailTab::from_key)
            .unwrap_or(DetailTab::Overview);
        let restored_subject = restore.last_subject.clone();
        let ui_sig = format!(
            "{}|{}|{}",
            restored_nav.as_key(),
            restored_subject.clone().unwrap_or_default(),
            restored_tab.as_key()
        );
        NavigatorApp {
            tx,
            rx,
            analysis: None,
            realign: None,
            simple_realign_confirm: None,
            cancelling: false,
            edit_subject: None,
            edit_kit: None,
            edit_mdka: None,
            confirm_delete: None,
            confirm_clear: None,
            confirm_reset_haplo: None,
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
            window_size: None,
            // The remembered size, when there is one. It starts here, so that an unchanged window
            // never writes the settings again.
            saved_window_size: settings.window_size,
            window_size_changed_at: 0.0,
            window_restored: false,
            startup_frames: 0,
            pending_restore_subject: restored_subject,
            saved_ui_sig: Some(ui_sig),
            nav: restored_nav,
            // Pinned mode (env / settings) wins; else default Simple provisionally and let the
            // first-run workspace heuristic adjust once subjects/projects load (see
            // `apply_ui_mode_heuristic`).
            ui_mode: navigator_app::configured_ui_mode().unwrap_or(UiMode::Simple),
            ui_mode_pinned: navigator_app::configured_ui_mode().is_some(),
            subject_brief: None,
            subject_brief_loading: false,
            simple_subject_filter: String::new(),
            simple_panel: SimplePanel::default(),
            ai_enabled: navigator_app::llm::llm_config().enabled,
            brief_narration: None,
            narration_stream: None,
            narrating: false,
            chat_history: Vec::new(),
            chat_input: String::new(),
            chat_pending: false,
            signal_narration: Vec::new(),
            signal_stream: None,
            signal_narrating: None,
            detail_tab: restored_tab,
            // Persisted choice wins; else honor $LANG (e.g. "es_ES.UTF-8") when it names a
            // supported locale; else English.
            lang: crate::i18n::load_lang()
                .or_else(|| std::env::var("LANG").ok().and_then(|l| crate::i18n::Lang::parse(&l)))
                .unwrap_or(crate::i18n::Lang::En),
            // Persisted theme wins; default dark.
            dark_mode: dark,
            // A persisted manual scale takes precedence; otherwise probe the monitor on frame 1.
            scale_probed: settings.ui_scale.is_some(),
            show_settings: false,
            settings_form: SettingsForm::from_settings(),
            settings_tab: SettingsTab::default(),
            llm_models: Vec::new(),
            llm_testing: false,
            llm_test_msg: None,
            subjects_table_ctl: TableControls::sorted_by(0),
            data_epoch: 0,
            subject_rows: SubjectRowCache::default(),
            report_rows: ReportRowCache::default(),
            y_profile_rows: VariantRows::default(),
            mt_profile_rows: VariantRows::default(),
            auto_profile_rows: VariantRows::default(),
            subjects_collapsed: false,
            projects_collapsed: false,
            overview: Vec::new(),
            selected_project: None,
            return_to_project: None,
            project_report: Vec::new(),
            project_str_chart: None,
            project_str_loading: false,
            project_blocktree: None,
            project_blocktree_loading: false,
            blocktree_zoom: 1.0,
            blocktree_review: None,
            blocktree_selected: None,
            blocktree_recentre: false,
            samples: Vec::new(),
            all_biosamples: Vec::new(),
            haplo_summary: std::collections::HashMap::new(),
            subject_status: std::collections::HashMap::new(),
            selected_sample: None,
            runs: Vec::new(),
            consensus_y: None,
            consensus_mt: None,
            descent_reports: Vec::new(),
            descent_loading: Vec::new(),
            branch_node_y: String::new(),
            branch_node_mt: String::new(),
            branch_reports: Vec::new(),
            branch_loading: Vec::new(),
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
            estimating_deep_ancestry: false,
            painting: None,
            painting_running: false,
            roh: None,
            roh_running: false,
            archaic: None,
            archaic_running: false,
            archaic_segments: None,
            archaic_segments_running: false,
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
            all_alignments: Vec::new(),
            project_realignable: None,
            project_realignable_asked: None,
            ibd_src_a: None,
            ibd_src_b: None,
            ibd_other_subject: None,
            ibd_other_picking: false,
            ibd_other_filter: String::new(),
            ibd_result: None,
            running_ibd: false,
            identity: None,
            ibd_suggestions: Vec::new(),
            loading_ibd_suggestions: false,
            ibd_intros: std::collections::HashMap::new(),
            exchange_results: Vec::new(),
            exchange_busy: false,
            matching: Vec::new(),
            matching_tab: MatchingTab::default(),
            dismissed_candidates: std::collections::HashSet::new(),
            matching_subject: None,
            matching_subject_picking: false,
            matching_subject_filter: String::new(),
            consent_prompt: None,
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
            index_progress: None,
            update_info: None,
            analyzing: false,
            deep_progress: None,
            maintenance: None,
            maintenance_surveying: false,
            chore_running: None,
            chore_last: None,
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
                run_test_type: "WGS".into(),
                str_panel: "Y-37".into(),
                str_provider: "FTDNA".into(),
                str_source: "DIRECT_TEST".into(),
                chip_provider: AUTO_DETECT.into(),
                variant_source_type: "IMPORTED".into(),
                ..Forms::default()
            },
            status: "Loading…".into(),
            diagnosis: None,
            show_diagnosis: false,
        }
    }

    /// Remember the window size across launches, and fit it to the screen. `screen_rect`,
    /// `monitor_size` and [`egui::ViewportCommand::InnerSize`] are all in egui points at the current
    /// zoom. So a `screen_rect` on disk, restored at the same zoom from disk, round-trips exactly.
    ///
    /// The one-time restore and fit runs through `InnerSize` at runtime. It does not use the startup
    /// `ViewportBuilder`, which sets the size before the UI scale applies. It waits until the zoom
    /// settles (`startup_frames >= 2` and `scale_probed`), and until the monitor size is there. It
    /// then clamps the remembered size to fit the screen, so a size that is too large, from a bigger
    /// display, shrinks.
    ///
    /// After that, a change of size goes to [`AppSettings`]. The code debounces the write until the
    /// size settles, and it writes one more time on window close.
    fn manage_window_geometry(&mut self, ctx: &egui::Context) {
        self.startup_frames = self.startup_frames.saturating_add(1);
        let size = ctx.screen_rect().size();
        let cur = [size.x, size.y];
        // Skip degenerate sizes (e.g. while minimized) so we never persist a collapsed window.
        if cur[0] < 200.0 || cur[1] < 200.0 {
            return;
        }

        // One-time restore + fit-to-screen, after the zoom has settled so `cur`/target share a unit.
        if !self.window_restored {
            let monitor = ctx.input(|i| i.viewport().monitor_size);
            // Wait until the UI scale (zoom) settles, and until the monitor size is there.
            let Some(mon) = monitor.filter(|_| self.scale_probed && self.startup_frames >= 2) else {
                ctx.request_repaint(); // keep frames coming until we can restore
                return; // don't track/save until the restore has run
            };
            let desired = self.saved_window_size.unwrap_or(cur);
            let [tw, th] = fit_window_to_monitor(desired, [mon.x, mon.y], MIN_WINDOW);
            if (cur[0] - tw).abs() > 1.0 || (cur[1] - th).abs() > 1.0 {
                ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(egui::vec2(tw, th)));
            }
            self.window_restored = true;
        }

        // Keep `self.window_size` current, because the on-exit save reads it. Persist shortly after
        // a resize settles, so that one drag becomes one write. The *definitive* final save is
        // [`Self::on_exit`], which a window close reaches, and so does Cmd+Q on macOS. Cmd+Q
        // terminates through `applicationWillTerminate:`, and it never runs a `close_requested`
        // update frame, so that flag alone lost the size. This save inside the loop means that a
        // hard kill still leaves a recent size on disk.
        if self.window_size != Some(cur) {
            self.window_size = Some(cur);
            self.window_size_changed_at = self.frame_time;
            // Make sure a frame fires after the resize stops, so that the settled save runs even
            // when nothing else happens.
            ctx.request_repaint_after(std::time::Duration::from_millis(500));
        }
        if self.frame_time - self.window_size_changed_at >= 0.5 {
            self.persist_window_size(cur);
        }
    }

    /// Write the window size to [`AppSettings`]. It loads, changes and saves, so that the other
    /// settings stay. It drops the write when the value already matches disk. The debounced save
    /// inside the loop, and the on-exit save, both use it.
    fn persist_window_size(&mut self, size: [f32; 2]) {
        if self.saved_window_size == Some(size) {
            return;
        }
        let mut settings = AppSettings::load();
        settings.window_size = Some(size);
        if settings.save().is_ok() {
            self.saved_window_size = Some(size);
        }
    }

    /// Signature of the current navigation state: `view | subject | detail-tab`.
    fn ui_state_sig(&self) -> String {
        format!(
            "{}|{}|{}",
            self.nav.as_key(),
            self.selected_sample.map(|g| g.0.to_string()).unwrap_or_default(),
            self.detail_tab.as_key()
        )
    }

    /// Persist the navigation position to [`AppSettings`] when it changes: the view, the focused
    /// subject, and the detail tab. The next launch then opens where the user left off. This runs on
    /// each frame, and the signature guard means a write happens only on a real navigation change.
    /// It loads, changes and saves, so that the window size and every other setting stay.
    fn persist_ui_state(&mut self) {
        // Wait until the one-time subject restore applies, which happens when the subject list
        // first loads. If not, the early frames, before anything selects a subject, would write over
        // the remembered subject with `None`.
        if self.pending_restore_subject.is_some() {
            return;
        }
        let sig = self.ui_state_sig();
        if self.saved_ui_sig.as_deref() == Some(sig.as_str()) {
            return;
        }
        let mut settings = AppSettings::load();
        settings.last_nav = Some(self.nav.as_key().to_string());
        settings.last_subject = self.selected_sample.map(|g| g.0.to_string());
        settings.last_detail_tab = Some(self.detail_tab.as_key().to_string());
        if settings.save().is_ok() {
            self.saved_ui_sig = Some(sig);
        }
    }
}

impl eframe::App for NavigatorApp {
    /// The final, reliable window-size save. eframe calls it on shutdown, from `save_and_destroy` on
    /// `LoopExiting`. Cmd+Q on macOS reaches that through `applicationWillTerminate:`, and that path
    /// runs no `close_requested` update frame. `manage_window_geometry` keeps `self.window_size`
    /// current.
    fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {
        if let Some(size) = self.window_size {
            self.persist_window_size(size);
        }
        self.persist_ui_state();
    }

    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.frame_time = ctx.input(|i| i.time);
        run_auto_scale(&mut self.scale_probed, &mut self.settings_form, ctx);
        self.manage_window_geometry(ctx);
        // While an analysis runs, keep the paint going, so that the spinner and the elapsed timer
        // animate. That holds even during a long step that emits no event, for example whole-genome
        // coverage.
        if self.analysis.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(120));
        }
        self.drain_events();
        // Mode upkeep: the first-run heuristic, until something pins the mode, and the auto-select
        // of the one subject in Simple.
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
                // PULL reconcile is an operation on a PDS repo, so it needs a real PDS (OAuth)
                // account. A local did:key identity covers federation and exchange only, and it has
                // no PDS repo. So show the control disabled, with a reason. Do not let it fail with
                // a "not signed in" that puzzles the user.
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
                // This appears only when a diagnosis exists, because an error with none must not
                // suggest there is more to see. It opens the modal again after a dismiss.
                if self.diagnosis.is_some() && ui.button(self.tr("status.details")).clicked() {
                    self.show_diagnosis = true;
                }
            });
        });
        // The batch, compare and add-to-project controls of the action bar are power-user
        // features, and they appear in Advanced only.
        if self.nav == Nav::Subjects && self.ui_mode == UiMode::Advanced {
            self.action_bar(ctx);
        }
        self.left_panel(ctx);
        egui::CentralPanel::default().show(ctx, |ui| match self.nav {
            Nav::Dashboard => self.dashboard_central(ui),
            Nav::Subjects => self.subjects_central(ui),
            Nav::Projects => self.projects_central(ui),
            Nav::Matching => self.matching_central(ui),
            Nav::Community => self.community_central(ui),
        });
        self.analysis_modal(ctx);
        self.blocktree_review_modal(ctx);
        self.diagnosis_modal(ctx);
        self.update_modal(ctx);
        self.edit_subject_modal(ctx);
        self.add_kit_modal(ctx);
        self.edit_mdka_modal(ctx);
        self.delete_subject_modal(ctx);
        self.clear_subject_modal(ctx);
        self.simple_realign_confirm_modal(ctx);
        self.consent_modal(ctx);
        self.reset_haplo_modal(ctx);
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

        // Persist the navigation position after the UI has processed this frame's clicks (view /
        // focused subject / detail tab). Guarded so it writes only on an actual change.
        self.persist_ui_state();
    }
}

/// Render a depth histogram (`bin d` = bases observed at depth `d`, top bin = ≥255) as an
const STR_CONFLICT: egui::Color32 = egui::Color32::from_rgb(220, 150, 60);

/// The By-Panel view in the FTDNA or YSEQ style. It groups markers into tiers (Y-12, Y-25, and so
/// on). It draws each tier as transposed mini-grids: a marker-name row over a value row, and no more
/// than 12 markers wide. A marker in conflict is amber.
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
    // No inner scroll area here. One vertical ScrollArea already wraps the detail panel. A second
    // one of fixed height inside it clips the panel tables and takes the wheel, so the page can not
    // scroll. Let the tiers flow into the page scroll.
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

/// A flat marker table with a filter: Marker | Panel | Value. With more than one provider it adds
/// ⚠ | Other, which is the value of the other provider that does not agree. A value in conflict is
/// amber.
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
        if !filter.is_empty() && ui.button("✖").clicked() {
            filter.clear();
        }
    });
    let f = filter.trim().to_uppercase();
    let cols = if multi { 5 } else { 3 };
    // Flow into the outer ScrollArea of the detail panel. No nested vertical scroll, because that
    // clips the table and takes the wheel.
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

/// The canonical short label and badge color for a consensus status. The Y and mt consensus card,
/// and the autosomal diploid card, both use it. `Novel` can not arise on the diploid path over panel
/// sites (see [`navigator_domain::consensus::reconcile_diploid`]), and the map holds it here to be
/// complete.
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

/// The shared draw for a multi-source consensus profile, Y or mtDNA, over the same generic engine.
/// It has a header with the counts, the lineage label and the provenance. Then comes a status
/// filter, then a grid with one row for each variant.
///
/// `variant_col` names the identity column ("SNP" or "Mutation"). `kind` labels the empty state.
/// `id_salt` keeps the scroll and grid ids of the two cards distinct. `snp_names` adds the
/// catalogued Y-SNP name at a site to a position-only or novel row. It is empty for mtDNA, whose
/// mutations already have names.
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
    epoch: u64,
    rows: &mut VariantRows,
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
        if !query.is_empty() && ui.small_button("✖").clicked() {
            query.clear();
        }
    });
    let q = query.to_ascii_lowercase();

    let state_label = |st: YState| match st {
        YState::Derived => "derived",
        YState::Ancestral => "ancestral",
        YState::NoCall => "no-call",
    };
    // Hold the table inside a scroll pane of fixed height. On a WGS these run thousands of rows,
    // and they would otherwise force an endless page scroll. The status filter and the text filter
    // narrow the list, and a cap limits a pathological profile. The header and filter row above
    // stays fixed, and only the grid scrolls.
    const CAP: usize = 2000;
    // A catalogued Y-SNP name at this site (for a position-only / novel row).
    let cataloged_at = |v: &navigator_domain::consensus::ConsensusVariant| {
        if v.name.is_empty() {
            snp_names.get(&v.position).map(String::as_str)
        } else {
            None
        }
    };
    // A cache holds which variants match. The view draws only `CAP` rows, but the total needs a
    // scan of the whole profile, and this function runs on every frame. `snp_names` feeds the
    // search, and an event loads it, so `epoch` covers it too.
    let status = *filter;
    let matching = rows.get(epoch, status, &q, &profile.variants, |v| {
        if status.is_some_and(|f| v.status != f) {
            return false;
        }
        q.is_empty()
            || v.name.to_ascii_lowercase().contains(&q)
            || v.position.to_string().contains(&q)
            || cataloged_at(v).is_some_and(|n| n.to_ascii_lowercase().contains(&q))
    });
    let total_match = matching.len();
    let shown = total_match.min(CAP);
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
                    for &i in &matching[..shown] {
                        let v = &profile.variants[i as usize];
                        let cataloged = cataloged_at(v);
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
                            // Show the actual consensus base alongside its derived/ancestral reading.
                            match v.consensus_base.as_deref().filter(|b| !b.is_empty()) {
                                Some(base) => ui.label(format!("{base}  {}", state_label(v.consensus))),
                                None => ui.label(state_label(v.consensus)),
                            };
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
                                        YState::Derived => "✔",
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
        egui::RichText::new(format!("{shown} of {total_match} matching variants"))
            .weak()
            .small(),
    );
    if total_match > CAP {
        ui.label(egui::RichText::new(format!("…and {} more — filter to narrow", total_match - CAP)).weak());
    }
}

/// The draw for the autosomal diploid consensus profile: the 0/1/2 sibling of
/// [`draw_consensus_profile`]. It has a header with confirmed, conflict, single and the confidence,
/// then a status filter, then a grid with one row for each site:
/// `Site (rsID) | GT (0/0,0/1,1/1) | Status | Sources (dosage of each source)`.
fn draw_diploid_profile(
    ui: &mut egui::Ui,
    profile: &navigator_app::DiploidProfile,
    filter: &mut Option<YVariantStatus>,
    query: &mut String,
    epoch: u64,
    rows: &mut VariantRows,
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
        if !query.is_empty() && ui.small_button("✖").clicked() {
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
    // The panel has ~1.2M sites. A Grid with no virtualization lays out every row on every frame,
    // and the app then stops to respond. Draw fixed-width columns through ScrollArea::show_rows,
    // which builds only the visible slice.
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
    // Hold the rows inside a scroll pane of fixed height, because the panel is ~1.2M sites. A hard
    // cap builds only `shown` widgets on each frame, and never the full panel. The status filter and
    // the text filter narrow it further. A cache holds which sites match (`rows`), because the count
    // needs a scan of the whole panel, and this function runs on every frame.
    const CAP: usize = 2000;
    let status = *filter;
    let matching = rows.get(epoch, status, &q, &profile.variants, |v| {
        status.map_or(true, |f| v.status == f) && matches(v)
    });
    let total_match = matching.len();
    let shown = total_match.min(CAP);
    let pane_h = profile_pane_height(ui, profile.variants.len());
    egui::ScrollArea::vertical()
        .id_salt("diploid_profile_scroll")
        .max_height(pane_h)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            for &i in &matching[..shown] {
                render_row(ui, &profile.variants[i as usize]);
            }
        });
    ui.label(
        egui::RichText::new(format!("{shown} of {total_match} matching sites"))
            .weak()
            .small(),
    );
    if total_match > CAP {
        ui.label(egui::RichText::new(format!("…and {} more — filter to narrow", total_match - CAP)).weak());
    }
}

/// The mtDNA mutation list against rCRS, in groups by region (HVR2, Coding, HVR1). It is the classic
/// mtDNA result. The `variants` come from a derivation against the bundled rCRS, and the notation is
/// the standard mtDNA form.
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

#[cfg(test)]
mod window_geometry_tests {
    use super::{fit_window_to_monitor, DEFAULT_WINDOW, MIN_WINDOW};

    #[test]
    fn size_that_fits_is_unchanged() {
        // A comfortable window on a 2560×1440 monitor stays as it is.
        let got = fit_window_to_monitor([1600.0, 1000.0], [2560.0, 1440.0], MIN_WINDOW);
        assert_eq!(got, [1600.0, 1000.0]);
    }

    #[test]
    fn oversize_is_shrunk_to_fit_the_screen() {
        // A size remembered from a big display, opened on a 1440×900 laptop, must fit within it.
        let mon = [1440.0, 900.0];
        let got = fit_window_to_monitor([3000.0, 2000.0], mon, MIN_WINDOW);
        assert!(got[0] <= mon[0] && got[1] <= mon[1], "must fit: {got:?} in {mon:?}");
        assert!(
            got[0] <= mon[0] * 0.98 + 0.5 && got[1] <= mon[1] * 0.94 + 0.5,
            "margin respected"
        );
    }

    #[test]
    fn never_below_minimum() {
        let got = fit_window_to_monitor([300.0, 200.0], [2560.0, 1440.0], MIN_WINDOW);
        assert_eq!(got, MIN_WINDOW);
    }

    #[test]
    fn tiny_monitor_clamps_to_minimum_not_below() {
        // A monitor smaller than the minimum: the floor wins (the window can't usefully go smaller).
        let got = fit_window_to_monitor(DEFAULT_WINDOW, [800.0, 600.0], MIN_WINDOW);
        assert_eq!(got, MIN_WINDOW);
    }
}

#[cfg(test)]
mod nav_persistence_tests {
    use super::{DetailTab, Nav};

    #[test]
    fn nav_keys_round_trip() {
        for nav in [
            Nav::Dashboard,
            Nav::Subjects,
            Nav::Projects,
            Nav::Matching,
            Nav::Community,
        ] {
            assert_eq!(Nav::from_key(nav.as_key()), Some(nav));
        }
        assert_eq!(Nav::from_key("bogus"), None);
    }

    #[test]
    fn detail_tab_keys_round_trip() {
        for tab in [
            DetailTab::Overview,
            DetailTab::YDna,
            DetailTab::MtDna,
            DetailTab::Autosomal,
            DetailTab::Ancestry,
            DetailTab::Sources,
            DetailTab::IbdMatches,
        ] {
            assert_eq!(DetailTab::from_key(tab.as_key()), Some(tab));
        }
        assert_eq!(DetailTab::from_key("bogus"), None);
    }
}

/// Every icon the chrome renders must have a real glyph in the font family it renders with.
///
/// The **Proportional** family of egui is Ubuntu-Light, NotoEmoji-Regular and emoji-icon-font. That
/// is a much narrower set than "Unicode".
///
/// A character that looks plausible, with no check, is how the Simple-mode rail went out with `◆`,
/// `⚭` and `✓`. The nav went out with `🧬`. All four drew as empty tofu boxes. A glyph that is
/// missing stays invisible to every other test, and to the compiler. Only a look at the live app
/// catches it. So read the fonts instead.
#[cfg(test)]
mod icon_glyph_tests {
    use super::egui::{FontDefinitions, FontFamily};
    use super::SimplePanel;
    use ab_glyph::{Font, FontRef};

    /// True when at least one font in `Proportional`'s fallback chain has a glyph for `c`.
    ///
    /// It reads the `FontDefinitions::default()` of egui itself, and not a vendored copy of the
    /// `.ttf` files. So the test keeps its grip on the fonts the app ships, as egui moves forward.
    /// `glyph_id` returns 0, which is `.notdef`, the tofu box, when a font has no mapping for the
    /// character.
    fn renderable(c: char) -> bool {
        let defs = FontDefinitions::default();
        let chain = &defs.families[&FontFamily::Proportional];
        chain.iter().any(|name| {
            let data = &defs.font_data[name];
            FontRef::try_from_slice(&data.font).is_ok_and(|f| f.glyph_id(c).0 != 0)
        })
    }

    #[test]
    fn probe_agrees_with_known_bad_glyphs() {
        // This guards the test itself. If these ever come back as drawable, the check broke, and
        // the fonts did not improve. The second row is the near-miss set. Each one is a character
        // somebody reached for because it looked right, and each one went out as a box: `✕` beside
        // `✖`, `✎` beside `✏`, `●` beside `⚫`, and `▲▼▸` beside `⏶⏷▶`.
        for c in ['◆', '⚭', '✓', '🧬', '✕', '✎', '✗', '●', '▲', '▼', '▸', '→'] {
            assert!(
                !renderable(c),
                "{c} (U+{:04X}) should be missing from Proportional",
                c as u32
            );
        }
        // And the replacements that took their place, so that this test catches a font change that
        // drops one, and a screenshot does not.
        for c in ['♂', '✖', '✏', '⚫', '⚪', '⏶', '⏷', '▶', '›'] {
            assert!(renderable(c), "sanity: {c} (U+{:04X}) is present", c as u32);
        }
    }

    /// Prints the coverage of the characters the app already uses, plus a candidate set. Run it
    /// when you choose an icon, instead of a character that only looks right:
    /// `cargo test -p navigator-ui --bin navigator report_glyph_coverage -- --ignored --nocapture`
    #[test]
    #[ignore = "diagnostic, not a gate — see the doc comment for how to run it"]
    fn report_glyph_coverage() {
        use std::collections::BTreeSet;
        let mut chars: BTreeSet<char> = BTreeSet::new();
        for lang in navigator_domain::i18n::Lang::all() {
            for (_, v) in navigator_domain::i18n::entries(*lang) {
                chars.extend(v.chars().filter(|c| !c.is_ascii()));
            }
        }
        for c in [
            '•', '·', '▪', '▫', '■', '◻', '◼', '⚫', '⚪', '🔴', '⏺', '☑', '☐', '✔', '✖', '★', '☆', '±', '‹', '›', '«',
            '»', '▶', '⏶', '⏷', '⬇', '🔽', '↘', '☰', '…', '≥', '×', '➕', '📁', '📥', '🔍', '🗑', '✏', '➡',
        ] {
            chars.insert(c);
        }
        for c in &chars {
            println!(
                "  {} U+{:04X} {c}",
                if renderable(*c) { "ok  " } else { "TOFU" },
                *c as u32
            );
        }
    }

    /// Every character of every translated string must be drawable.
    ///
    /// A check on the rail and nav icons below covers each one on its own. That left the ~1,900
    /// strings in the catalogs with no check at all. Those held eight more characters nothing can
    /// draw (`←→▾◆●✓✗🗂`). Among them were the status dot on the dashboard, and the check and the
    /// cross on the exchange screen. A scan of the whole catalog is the only version of this test
    /// that new copy can not outgrow.
    #[test]
    fn every_translated_string_is_renderable() {
        let mut bad: Vec<String> = Vec::new();
        for lang in navigator_domain::i18n::Lang::all() {
            for (key, value) in navigator_domain::i18n::entries(*lang) {
                for c in value.chars().filter(|c| !c.is_ascii() && !renderable(*c)) {
                    bad.push(format!(
                        "  {} / {key} = {value:?} -> {c:?} (U+{:04X})",
                        lang.code(),
                        c as u32
                    ));
                }
            }
        }
        bad.sort();
        bad.dedup();
        assert!(
            bad.is_empty(),
            "these strings contain characters egui cannot draw; they render as empty boxes:\n{}",
            bad.join("\n")
        );
    }

    /// Every string literal this crate draws must be drawable.
    ///
    /// `every_translated_string_is_renderable` covers the catalogs, and that is where copy for a
    /// person belongs. But icons do not live there. A button label like `ui.small_button("✎")` is a
    /// bare literal in the source, and the catalog scan can not see it. That is where the second
    /// round of tofu boxes turned up. They were the MDKA buttons, the kit-remove button on the
    /// Genealogy card, and every clear-filter `✕`. They were also the Y-STR agreement marks, the
    /// arrows of the sortable table, and the match-strength meter.
    ///
    /// A list by hand is the failure mode this test exists to remove. It reads the sources
    /// themselves, so a new icon has a check the moment somebody types it.
    ///
    /// A parse, and not a grep, is what makes that safe. Comments and doc comments are full of
    /// characters that nothing ever draws: `→`, `⇒`, and a `◆` that names the fault it describes.
    /// So the AST is the right input, because it holds literals and no comments. An override of
    /// `visit_attribute` also drops `#[doc = "…"]`.
    ///
    /// Two things are out of scope on purpose, because egui never draws them. The first is `cli.rs`,
    /// whose output goes to a terminal that draws in the font of the user. The second is a
    /// `#[cfg(test)]` module, whose literals are assertion messages.
    ///
    /// So is the rest of the workspace, and that is a real gap, and not an oversight. Two of the
    /// boxes this round fixed came from `navigator-app`: a consensus warning, and the reference
    /// notes on an import summary. Only this crate drew them.
    ///
    /// The invariant can not move to that crate as it stands. That crate also serves the CLI and
    /// the HTML exporter, where `→` and `✓` are correct. Which of its strings reach a window is a
    /// dataflow question, and not a syntactic one. When a lower crate builds a string for the UI,
    /// the author has to check it, and `report_glyph_coverage` is the tool for that.
    #[test]
    fn every_source_string_literal_is_renderable() {
        use syn::visit::Visit;

        /// Collects `(literal, offending char)` for every undrawable character in a file's literals.
        #[derive(Default)]
        struct Scan(Vec<(String, char)>);
        impl<'ast> Visit<'ast> for Scan {
            fn visit_lit_str(&mut self, lit: &'ast syn::LitStr) {
                let value = lit.value();
                for c in value.chars().filter(|c| !c.is_ascii() && !renderable(*c)) {
                    self.0.push((value.clone(), c));
                }
            }
            /// Doc comments arrive as `#[doc = "…"]`; nothing in any attribute is ever drawn.
            fn visit_attribute(&mut self, _: &'ast syn::Attribute) {}
            /// Assertion messages in `#[cfg(test)] mod tests` reach a terminal, not a window.
            fn visit_item_mod(&mut self, m: &'ast syn::ItemMod) {
                if !m.attrs.iter().any(is_cfg_test) {
                    syn::visit::visit_item_mod(self, m);
                }
            }
            /// A macro body is a `TokenStream` that nothing parsed, and the visitor of syn steps
            /// over it. That would hide most of what this test is for, because `format!` builds
            /// almost every label. A walk of the tokens by hand costs one recursion. It also needs
            /// no guess about the grammar of the macro: a `LitStr` is a `LitStr` wherever it
            /// sits.
            fn visit_macro(&mut self, m: &'ast syn::Macro) {
                self.tokens(m.tokens.clone());
            }
        }

        impl Scan {
            fn tokens(&mut self, stream: proc_macro2::TokenStream) {
                for tree in stream {
                    match tree {
                        proc_macro2::TokenTree::Group(g) => self.tokens(g.stream()),
                        tree @ proc_macro2::TokenTree::Literal(_) => {
                            // Only string literals parse; numbers, chars and byte strings just fail.
                            if let Ok(lit) = syn::parse2::<syn::LitStr>(tree.into()) {
                                self.visit_lit_str(&lit);
                            }
                        }
                        _ => {}
                    }
                }
            }
        }

        /// True for `#[cfg(test)]`. It matches on the token text, which is stable enough for an
        /// attribute this conventional, and it avoids a nested-meta parser.
        fn is_cfg_test(attr: &syn::Attribute) -> bool {
            attr.path().is_ident("cfg") && attr.parse_args::<syn::Path>().is_ok_and(|p| p.is_ident("test"))
        }

        /// Every `.rs` file under `dir`, recursively.
        fn sources(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
            for entry in std::fs::read_dir(dir).expect("read src dir").flatten() {
                let path = entry.path();
                if path.is_dir() {
                    sources(&path, out);
                } else if path.extension().is_some_and(|e| e == "rs") {
                    out.push(path);
                }
            }
        }

        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut files = Vec::new();
        sources(&root, &mut files);
        files.retain(|p| p.file_name().is_none_or(|n| n != "cli.rs"));
        files.sort();
        assert!(!files.is_empty(), "found no sources under {}", root.display());

        let mut bad: Vec<String> = Vec::new();
        for path in &files {
            let text = std::fs::read_to_string(path).expect("read source");
            let file = syn::parse_file(&text).expect("parse source");
            let mut scan = Scan::default();
            scan.visit_file(&file);
            let name = path.strip_prefix(&root).unwrap_or(path).display();
            for (value, c) in scan.0 {
                bad.push(format!("  {name} / {value:?} -> {c:?} (U+{:04X})", c as u32));
            }
        }
        bad.sort();
        bad.dedup();
        assert!(
            bad.is_empty(),
            "these source literals contain characters egui cannot draw; they render as empty boxes:\n{}",
            bad.join("\n")
        );
    }

    /// The marks in the asset-status line, which the code builds inline, and no catalog
    /// translates.
    #[test]
    fn asset_status_marks_are_renderable() {
        use crate::charts::{MARK_ABSENT, MARK_PRESENT, MARK_VERIFIED};
        for mark in [MARK_VERIFIED, MARK_PRESENT, MARK_ABSENT] {
            for c in mark.chars() {
                assert!(
                    renderable(c),
                    "asset mark {c:?} (U+{:04X}) has no glyph — it renders as a tofu box",
                    c as u32
                );
            }
        }
    }

    #[test]
    fn icon_glyphs_are_renderable() {
        for (panel, icon, _) in SimplePanel::ALL {
            for c in icon.chars() {
                assert!(
                    renderable(c),
                    "{panel:?} rail icon {c:?} (U+{:04X}) has no glyph — it renders as a tofu box",
                    c as u32
                );
            }
        }
        // The nav strip's icons, which live inline in `chrome::nav_bar`.
        for icon in ['📊', '👤', '👥', '📁', '🔗', '💬'] {
            assert!(
                renderable(icon),
                "nav icon {icon:?} (U+{:04X}) has no glyph — it renders as a tofu box",
                icon as u32
            );
        }
    }
}
