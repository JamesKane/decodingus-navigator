//! The egui front end. Thin by design: it holds immutable view-state, renders it, and
//! dispatches [`Command`]s; all data/logic lives behind the worker + app. The worker's
//! `wake` callback calls `request_repaint`, so events refresh the view promptly.

use std::path::PathBuf;
use std::sync::mpsc::Receiver;

use eframe::egui;
use navigator_app::{
    AncestryResult, AncestrySegment, AuditEntry, BuildNeed, CallState, CompatibilityLevel, Consensus,
    Coverage, DenovoCall, DnaType, HaploAssignment, HeteroplasmySite, IbdComparison,
    IdentityVerification, PanelGenotype, PrivateBucket, PrivateClass, ProjectOverview,
    ProjectSampleReport, ReadMetrics, ReconciledVariant, SexInferenceResult, SourceType,
    SuperPopulationSummary, SvAnalysisResult, VariantStatus, VerificationStatus,
};
use navigator_domain::ancestry::{population_color, population_lonlat, population_name, population_super};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::chipprofile::{self, ChipProfile};
use navigator_domain::mtdna::MtdnaSequence;
use navigator_domain::strprofile::{self, StrProfile};
use navigator_domain::testtype;
use navigator_domain::variants::{VariantCall, VariantSet};
use navigator_domain::workspace::{Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, SequenceRun};
use tokio::sync::mpsc::UnboundedSender;

use crate::worker::{self, Command, Event, NewBiosample, PanelInfo, YMask};

/// PCA scatter backdrop: the alignment it was loaded for + reference centroids `(code, pc1, pc2)`.
type PcaReferenceState = (i64, Vec<(String, f64, f64)>);

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
}

/// Sub-tabs of the subject detail panel.
#[derive(Clone, Copy, PartialEq, Eq)]
enum DetailTab {
    Overview,
    YDna,
    MtDna,
    Ancestry,
    IbdMatches,
    DataSources,
}

impl DetailTab {
    /// `(tab, i18n key)` in display order.
    const ALL: [(DetailTab, &'static str); 6] = [
        (DetailTab::Overview, "detail.overview"),
        (DetailTab::YDna, "detail.ydna"),
        (DetailTab::MtDna, "detail.mtdna"),
        (DetailTab::Ancestry, "detail.ancestry"),
        (DetailTab::IbdMatches, "detail.ibd"),
        (DetailTab::DataSources, "detail.datasources"),
    ];
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
            DataDelete::Run { id, guid, .. } => Command::DeleteSequenceRun { id, biosample_guid: guid },
            DataDelete::Alignment { id, run_id, .. } => Command::DeleteAlignment { id, sequence_run_id: run_id },
            DataDelete::Str { id, guid, .. } => Command::DeleteStrProfile { id, biosample_guid: guid },
            DataDelete::Variant { id, guid, .. } => Command::DeleteVariantSet { id, biosample_guid: guid },
            DataDelete::Chip { id, guid, .. } => Command::DeleteChipProfile { id, biosample_guid: guid },
            DataDelete::Mtdna { id, guid, .. } => Command::DeleteMtdnaSequence { id, biosample_guid: guid },
        }
    }
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
    /// Data-source row pending delete confirmation (Some ⇒ the confirm dialog is shown).
    confirm_data_delete: Option<DataDelete>,
    /// Subject being assigned to a project: (subject, selected project or None). Some ⇒ picker shown.
    assign_project: Option<(SampleGuid, Option<i64>)>,
    /// Project being edited (Some ⇒ the project Edit modal is shown).
    edit_project: Option<EditProject>,
    /// Project pending delete confirmation: (id, name). Some ⇒ the confirm dialog is shown.
    confirm_delete_project: Option<(i64, String)>,
    /// Current frame's egui time (seconds), captured at the top of `update`.
    frame_time: f64,
    /// Selected primary navigation tab.
    nav: Nav,
    /// Selected subject-detail sub-tab.
    detail_tab: DetailTab,
    /// Active UI language.
    lang: crate::i18n::Lang,
    /// Dark (default) vs light theme.
    dark_mode: bool,
    /// Subjects-list filter text.
    subject_search: String,
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
    /// Last mtDNA-from-alignment haplogroup assignment: (alignment id, assignment).
    mt_haplogroup: Option<(i64, HaploAssignment)>,
    /// Last ancestry estimate: (alignment id, result). `None` result = computed, no estimate.
    ancestry: Option<(i64, Option<AncestryResult>)>,
    /// Donor-level ancestry (best across the subject's sources): (source alignment id, result).
    donor_ancestry: Option<(i64, AncestryResult)>,
    /// Donor-level private-Y union across the subject's sources.
    donor_private_y: Option<PrivateBucket>,
    estimating_ancestry: bool,
    /// Live genotyping progress for the in-flight estimate: (alignment id, done, total) contigs.
    ancestry_progress: Option<(i64, usize, usize)>,
    /// Reference population centroids for the PCA scatter: (alignment id, [(code, pc1, pc2)]).
    pca_reference: Option<PcaReferenceState>,
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

/// The workbench accent (primary buttons, selection, active tabs).
const ACCENT: egui::Color32 = egui::Color32::from_rgb(45, 125, 246);
/// Destructive-action red (Delete) — used by the subject header in Phase 2.
#[allow(dead_code)]
const DANGER: egui::Color32 = egui::Color32::from_rgb(220, 60, 60);

/// Apply the Decoding-Us workbench look: a dark (or light) palette with the accent blue,
/// rounded widgets, and roomier spacing — the visual base that closes most of the gap to the
/// Scala Workbench. Re-applied on theme toggle.
fn apply_theme(ctx: &egui::Context, dark: bool) {
    use egui::{Color32, Rounding, Stroke};
    let mut style = (*ctx.style()).clone();
    let mut v = if dark { egui::Visuals::dark() } else { egui::Visuals::light() };

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
const SUBJECT_COLS: [(&str, f32); 7] =
    [("ID", 150.0), ("Name", 150.0), ("Y-DNA", 80.0), ("mtDNA", 80.0), ("Sex", 70.0), ("Center", 130.0), ("Status", 90.0)];

/// Short, stable subject id for the table's ID column (first 8 chars of the guid + ellipsis).
fn short_guid(b: &Biosample) -> String {
    let s = b.guid.0.to_string();
    if s.len() > 9 {
        format!("{}…", &s[..9])
    } else {
        s
    }
}

/// Paint a table header row at the column offsets used by [`table_row`].
fn table_header(ui: &mut egui::Ui, cols: &[(&str, f32)]) {
    let total_w: f32 = cols.iter().map(|c| c.1).sum();
    let (rect, _) = ui.allocate_exact_size(egui::vec2(total_w, 24.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let mut x = rect.left() + 8.0;
    let color = ui.visuals().weak_text_color();
    for (name, w) in cols {
        painter.text(
            egui::pos2(x, rect.center().y),
            egui::Align2::LEFT_CENTER,
            name,
            egui::FontId::proportional(12.5),
            color,
        );
        x += w;
    }
    painter.hline(rect.x_range(), rect.bottom(), egui::Stroke::new(1.0, ui.visuals().widgets.noninteractive.bg_stroke.color));
}

/// Paint one clickable table row; returns true when clicked. `status_col` (if any) is rendered
/// as a small accent-coloured badge (the "Status" cell).
fn table_row(ui: &mut egui::Ui, cols: &[(&str, f32)], cells: &[String], selected: bool, status_col: Option<usize>) -> bool {
    let total_w: f32 = cols.iter().map(|c| c.1).sum();
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(total_w, 28.0), egui::Sense::click());
    let painter = ui.painter_at(rect);
    if selected {
        painter.rect_filled(rect, 4.0, ACCENT.gamma_multiply(0.6));
    } else if resp.hovered() {
        painter.rect_filled(rect, 4.0, ui.visuals().faint_bg_color);
    }
    let text_color = if selected { egui::Color32::WHITE } else { ui.visuals().text_color() };
    let mut x = rect.left() + 8.0;
    for (i, ((_, w), val)) in cols.iter().zip(cells).enumerate() {
        let cy = rect.center().y;
        if Some(i) == status_col && val != "-" {
            // a muted pill behind the status text
            let galley = painter.layout_no_wrap(val.clone(), egui::FontId::proportional(11.5), egui::Color32::from_rgb(225, 190, 90));
            let pad = egui::vec2(7.0, 3.0);
            let pill = egui::Rect::from_min_size(egui::pos2(x, cy - galley.size().y / 2.0 - pad.y), galley.size() + pad * 2.0);
            painter.rect_filled(pill, 8.0, egui::Color32::from_rgb(70, 58, 28));
            painter.galley(pill.min + pad, galley, egui::Color32::PLACEHOLDER);
        } else {
            painter.text(egui::pos2(x, cy), egui::Align2::LEFT_CENTER, val, egui::FontId::proportional(13.0), text_color);
        }
        x += w;
    }
    resp.clicked()
}

/// A rounded section card with an optional bold title (the Data Sources look).
fn card(ui: &mut egui::Ui, title: &str, body: impl FnOnce(&mut egui::Ui)) {
    egui::Frame::group(ui.style())
        .fill(ui.visuals().faint_bg_color)
        .rounding(egui::Rounding::same(8.0))
        .inner_margin(egui::Margin::same(14.0))
        .show(ui, |ui| {
            ui.set_width(ui.available_width());
            if !title.is_empty() {
                ui.label(egui::RichText::new(title).strong().size(15.0));
                ui.add_space(8.0);
            }
            body(ui);
        });
}

/// A small rounded chip/badge (provider tag, Y/mt badge).
fn chip(ui: &mut egui::Ui, text: &str, bg: egui::Color32, fg: egui::Color32) {
    let font = egui::FontId::proportional(11.5);
    let galley = ui.painter().layout_no_wrap(text.to_string(), font, fg);
    let pad = egui::vec2(7.0, 3.0);
    let (rect, _) = ui.allocate_exact_size(galley.size() + pad * 2.0, egui::Sense::hover());
    ui.painter().rect_filled(rect, 6.0, bg);
    ui.painter().galley(rect.min + pad, galley, egui::Color32::PLACEHOLDER);
}

/// 3-letter provider abbreviation for the run chip (PACBIO → PAC).
fn provider_abbrev(platform: &str) -> String {
    let p = platform.trim();
    if p.is_empty() || p.eq_ignore_ascii_case("unknown") {
        "SEQ".into()
    } else {
        p.chars().take(3).collect::<String>().to_uppercase()
    }
}

/// Compact read count: 9_900 → "9.9K", 1_200_000 → "1.2M".
fn fmt_reads(n: Option<i64>) -> String {
    match n {
        None => "—".into(),
        Some(v) if v >= 1_000_000 => format!("{:.1}M", v as f64 / 1e6),
        Some(v) if v >= 1_000 => format!("{:.1}K", v as f64 / 1e3),
        Some(v) => v.to_string(),
    }
}

/// A centered empty-state placeholder for a work area with no selection.
fn empty_state(ui: &mut egui::Ui, title: &str, hint: &str) {
    ui.add_space(56.0);
    ui.vertical_centered(|ui| {
        ui.heading(title);
        ui.label(egui::RichText::new(hint).weak());
    });
}

/// Hint shown in an analysis tab when the subject has no analyzable alignment (the default is
/// auto-selected when one exists, so this means "no sequencing data yet").
fn pick_alignment_hint(ui: &mut egui::Ui, msg: &str) {
    ui.add_space(8.0);
    ui.label(egui::RichText::new(msg).weak());
}

/// A dashboard stat tile: a big number over a muted label, in a rounded card.
fn stat_card(ui: &mut egui::Ui, label: &str, value: usize) {
    egui::Frame::group(ui.style())
        .fill(ui.visuals().faint_bg_color)
        .inner_margin(egui::Margin::symmetric(18.0, 14.0))
        .show(ui, |ui| {
            ui.set_min_width(120.0);
            ui.vertical(|ui| {
                ui.label(egui::RichText::new(value.to_string()).size(28.0).strong());
                ui.label(egui::RichText::new(label).weak());
            });
        });
}

/// Format an optional mean/median depth (one decimal), "—" when not computed.
fn fmt_depth(o: Option<f64>) -> String {
    o.map(|v| format!("{v:.1}")).unwrap_or_else(|| "—".into())
}

/// Format an optional fraction (0–1) as a percentage, "—" when not computed.
fn fmt_pct(o: Option<f64>) -> String {
    o.map(|v| format!("{:.1}%", v * 100.0)).unwrap_or_else(|| "—".into())
}

/// Draw the per-chromosome local-ancestry painting: one horizontal bar per autosome (each
/// normalized to full width), segments colored by ancestry, plus a legend of the ancestries shown.
fn draw_chromosome_painting(ui: &mut egui::Ui, segments: &[AncestrySegment]) {
    use std::collections::BTreeMap;
    // Group by contig, ordered by chromosome number.
    let mut by_contig: BTreeMap<i64, Vec<&AncestrySegment>> = BTreeMap::new();
    for s in segments {
        let n: i64 = s.contig.trim_start_matches("chr").parse().unwrap_or(99);
        by_contig.entry(n).or_default().push(s);
    }
    let label_w = 42.0;
    let bar_w = 300.0;
    let bar_h = 13.0;
    for (n, mut segs) in by_contig {
        segs.sort_by_key(|s| s.start);
        let (lo, hi) = (segs.first().unwrap().start, segs.last().unwrap().end.max(segs.first().unwrap().start + 1));
        let span = (hi - lo).max(1) as f32;
        ui.horizontal(|ui| {
            ui.allocate_ui(egui::vec2(label_w, bar_h), |ui| ui.label(format!("chr{n}")));
            let (rect, _) = ui.allocate_exact_size(egui::vec2(bar_w, bar_h), egui::Sense::hover());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 2.0, egui::Color32::from_gray(30));
            for s in &segs {
                let x0 = rect.left() + (s.start - lo) as f32 / span * rect.width();
                let x1 = rect.left() + (s.end - lo) as f32 / span * rect.width();
                let seg_rect = egui::Rect::from_min_max(egui::pos2(x0, rect.top()), egui::pos2(x1.max(x0 + 1.0), rect.bottom()));
                painter.rect_filled(seg_rect, 0.0, parse_hex_color(&population_color(&s.population_code)));
            }
        });
    }
    // Legend: distinct ancestries present.
    let mut seen: Vec<&str> = Vec::new();
    for s in segments {
        if !seen.contains(&s.population_code.as_str()) {
            seen.push(&s.population_code);
        }
    }
    ui.add_space(2.0);
    ui.horizontal_wrapped(|ui| {
        for code in seen {
            let (r, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
            ui.painter().circle_filled(r.center(), 4.0, parse_hex_color(&population_color(code)));
            ui.label(egui::RichText::new(population_name(code)).small());
            ui.add_space(6.0);
        }
    });
}

/// Points along a circle arc from angle `a0` to `a1` (radians), `steps`+1 samples.
fn arc_points(c: egui::Pos2, r: f32, a0: f32, a1: f32, steps: usize) -> Vec<egui::Pos2> {
    (0..=steps)
        .map(|i| {
            let t = a0 + (a1 - a0) * (i as f32 / steps as f32);
            egui::pos2(c.x + r * t.cos(), c.y + r * t.sin())
        })
        .collect()
}

/// Draw a donut chart of the super-population proportions (one wedge per super-population,
/// colored by continent), with the dominant share in the centre.
fn draw_ancestry_donut(ui: &mut egui::Ui, summary: &[SuperPopulationSummary]) {
    let size = 120.0;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(size, size), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let c = rect.center();
    let (r_out, r_in) = (size * 0.46, size * 0.28);
    let total: f32 = summary.iter().map(|s| s.percentage as f32).sum::<f32>().max(1.0);
    let mut a0 = -std::f32::consts::FRAC_PI_2; // start at 12 o'clock
    for s in summary {
        if s.percentage < 0.5 {
            continue;
        }
        let a1 = a0 + (s.percentage as f32 / total) * std::f32::consts::TAU;
        let code = s.populations.first().and_then(|c| population_super(c)).unwrap_or("");
        let mut pts = arc_points(c, r_out, a0, a1, 32);
        pts.extend(arc_points(c, r_in, a1, a0, 32)); // inner arc, reversed → closed ring sector
        painter.add(egui::epaint::PathShape {
            points: pts,
            closed: true,
            fill: parse_hex_color(&population_color(code)),
            stroke: egui::epaint::PathStroke::NONE,
        });
        a0 = a1;
    }
    if let Some(top) = summary.first() {
        painter.text(
            c,
            egui::Align2::CENTER_CENTER,
            format!("{:.0}%", top.percentage),
            egui::FontId::proportional(18.0),
            egui::Color32::WHITE,
        );
    }
}

/// Draw a schematic world map (equirectangular) with each contributing population plotted at its
/// homeland, the marker area proportional to its share and colored by continent.
fn draw_ancestry_map(ui: &mut egui::Ui, components: &[navigator_app::PopulationComponent]) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(360.0, 180.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 4.0, egui::Color32::from_gray(18));
    let proj = |lon: f32, lat: f32| -> egui::Pos2 {
        egui::pos2(
            rect.left() + (lon + 180.0) / 360.0 * rect.width(),
            rect.top() + (90.0 - lat) / 180.0 * rect.height(),
        )
    };
    // Faint equator + prime meridian for orientation.
    let grid = egui::Stroke::new(1.0, egui::Color32::from_gray(36));
    painter.line_segment([proj(-180.0, 0.0), proj(180.0, 0.0)], grid);
    painter.line_segment([proj(0.0, -90.0), proj(0.0, 90.0)], grid);
    // Largest shares last, so they render on top.
    let mut comps: Vec<&navigator_app::PopulationComponent> =
        components.iter().filter(|c| c.percentage >= 1.0).collect();
    comps.sort_by(|a, b| a.percentage.partial_cmp(&b.percentage).unwrap_or(std::cmp::Ordering::Equal));
    for c in comps {
        if let Some((lon, lat)) = population_lonlat(&c.population_code) {
            let p = proj(lon, lat);
            let radius = (2.5 + (c.percentage as f32).sqrt() * 2.4).clamp(2.5, 24.0);
            let col = parse_hex_color(&population_color(&c.population_code));
            painter.circle_filled(p, radius, col.gamma_multiply(0.55));
            painter.circle_stroke(p, radius, egui::Stroke::new(1.5, col));
        }
    }
}

/// Draw the super-population composition as a single stacked horizontal bar (segment widths =
/// proportions, colored by continent).
fn draw_composition_bar(ui: &mut egui::Ui, summary: &[SuperPopulationSummary]) {
    let w = ui.available_width().min(360.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(w, 16.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 3.0, egui::Color32::from_gray(30));
    let mut x = rect.left();
    for s in summary {
        let seg_w = rect.width() * (s.percentage as f32 / 100.0).clamp(0.0, 1.0);
        let seg = egui::Rect::from_min_size(egui::pos2(x, rect.top()), egui::vec2(seg_w, rect.height()));
        let code = s.populations.first().and_then(|c| population_super(c)).unwrap_or("");
        painter.rect_filled(seg, 0.0, parse_hex_color(&population_color(code)));
        x += seg_w;
    }
}

/// Draw a small PCA scatter: each reference population's centroid (colored, labeled) plus the
/// sample's projected point (white ○). Axes are PC1 (x) × PC2 (y), auto-scaled to the data.
fn draw_pca_scatter(ui: &mut egui::Ui, sample: (f64, f64), refs: &[(String, f64, f64)]) {
    let mut xs: Vec<f64> = refs.iter().map(|r| r.1).collect();
    let mut ys: Vec<f64> = refs.iter().map(|r| r.2).collect();
    xs.push(sample.0);
    ys.push(sample.1);
    let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for &x in &xs {
        xmin = xmin.min(x);
        xmax = xmax.max(x);
    }
    let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
    for &y in &ys {
        ymin = ymin.min(y);
        ymax = ymax.max(y);
    }
    // 8% padding; guard zero-range axes.
    let xpad = ((xmax - xmin) * 0.08).max(1.0);
    let ypad = ((ymax - ymin) * 0.08).max(1.0);
    xmin -= xpad;
    xmax += xpad;
    ymin -= ypad;
    ymax += ypad;

    let (rect, _) = ui.allocate_exact_size(egui::vec2(300.0, 240.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 4.0, egui::Color32::from_gray(18));
    let inset = 10.0;
    let map = |x: f64, y: f64| -> egui::Pos2 {
        let fx = ((x - xmin) / (xmax - xmin)) as f32;
        let fy = ((y - ymin) / (ymax - ymin)) as f32;
        egui::pos2(
            rect.left() + inset + fx * (rect.width() - 2.0 * inset),
            rect.bottom() - inset - fy * (rect.height() - 2.0 * inset), // invert y
        )
    };
    for (code, x, y) in refs {
        let p = map(*x, *y);
        let col = parse_hex_color(&population_color(code));
        painter.circle_filled(p, 5.0, col);
        painter.text(
            p + egui::vec2(7.0, 0.0),
            egui::Align2::LEFT_CENTER,
            code,
            egui::FontId::proportional(11.0),
            col,
        );
    }
    let sp = map(sample.0, sample.1);
    painter.circle_stroke(sp, 6.0, egui::Stroke::new(2.0, egui::Color32::WHITE));
    painter.text(
        sp + egui::vec2(8.0, 0.0),
        egui::Align2::LEFT_CENTER,
        "sample",
        egui::FontId::proportional(11.0),
        egui::Color32::WHITE,
    );
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
        ui.label("No match."); // free helper (no `self`); i18n when it takes a `lang` param
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
        apply_theme(&cc.egui_ctx, true);
        NavigatorApp {
            tx,
            rx,
            analysis: None,
            edit_subject: None,
            confirm_delete: None,
            confirm_data_delete: None,
            assign_project: None,
            edit_project: None,
            confirm_delete_project: None,
            frame_time: 0.0,
            nav: Nav::Subjects,
            detail_tab: DetailTab::Overview,
            // Persisted choice wins; else honor $LANG (e.g. "es_ES.UTF-8") when it names a
            // supported locale; else English.
            lang: crate::i18n::load_lang()
                .or_else(|| std::env::var("LANG").ok().and_then(|l| crate::i18n::Lang::parse(&l)))
                .unwrap_or(crate::i18n::Lang::En),
            dark_mode: true,
            subject_search: String::new(),
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
            mt_haplogroup: None,
            ancestry: None,
            donor_ancestry: None,
            donor_private_y: None,
            estimating_ancestry: false,
            ancestry_progress: None,
            pca_reference: None,
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
                Event::Noop => {}
                Event::Overview(v) => {
                    self.status = format!("{} project(s)", v.len());
                    self.overview = v;
                }
                Event::ProjectCreated(p) => {
                    self.select_project(p.id);
                    let _ = self.tx.send(Command::LoadOverview);
                }
                Event::ProjectsChanged => {
                    let _ = self.tx.send(Command::LoadOverview);
                    let _ = self.tx.send(Command::LoadAllBiosamples); // a deleted project clears assignments
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
                Event::ProjectAnalyzed { project_id, samples, coverage_done, y_done, sex_done, metrics_done, sv_done, errors } => {
                    self.analyzing = false;
                    self.status = format!(
                        "Analyzed {samples} sample(s): {coverage_done} coverage, {y_done} Y, {sex_done} sex, {metrics_done} metrics, {sv_done} SV{}",
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
                Event::MtHaplogroup { alignment_id, assignment } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("mtDNA haplogroup: {} (score {:.3})", top.name, top.score),
                        None => "No mtDNA haplogroup match".into(),
                    };
                    self.mt_haplogroup = Some((alignment_id, assignment));
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensus(guid)); // the mt call was recorded
                    }
                }
                Event::AncestryProgress { alignment_id, done, total } => {
                    self.ancestry_progress = Some((alignment_id, done, total));
                    self.status = format!("Genotyping ancestry panel: {done}/{total} contigs…");
                }
                Event::Ancestry { alignment_id, result } => {
                    self.estimating_ancestry = false;
                    self.ancestry_progress = None;
                    // Lead with the robust super-population rollup (fine-pop components are
                    // indicative but noisier on a continental-AIMs panel).
                    self.status = match &result {
                        Some(r) => match r.super_population_summary.first() {
                            Some(top) => format!(
                                "Ancestry: {} {:.1}% ({}/{} SNPs)",
                                top.super_population, top.percentage, r.snps_with_genotype, r.snps_analyzed
                            ),
                            None => "Ancestry: no estimate".into(),
                        },
                        None => "Ancestry: not yet computed".into(),
                    };
                    self.ancestry = Some((alignment_id, result));
                    // A fresh estimate may change the donor's best — reload the donor rollup.
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadDonorAncestry { biosample_guid: guid });
                    }
                }
                Event::PcaReference { alignment_id, points } => {
                    self.pca_reference = Some((alignment_id, points));
                }
                Event::AncestryPainting { alignment_id, segments } => {
                    self.painting_running = false;
                    self.ancestry_progress = None;
                    self.status = format!("Painted {} ancestry segments", segments.len());
                    self.painting = Some((alignment_id, segments));
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
                    // A fresh (self-masked) bucket was just cached — refresh the donor union.
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadDonorPrivateY { biosample_guid: guid });
                    }
                }
                Event::DataImported { biosample_guid, label } => {
                    self.status = format!("Imported {label}");
                    if self.selected_sample == Some(biosample_guid) {
                        // Reload every data section — detection picked one of them (a BAM/CRAM
                        // import auto-creates a sequencing run + alignment, so reload runs too).
                        let _ = self.tx.send(Command::LoadRuns(biosample_guid));
                        let _ = self.tx.send(Command::LoadStrProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadVariantSets(biosample_guid));
                        let _ = self.tx.send(Command::LoadChipProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadMtdna(biosample_guid));
                    }
                }
                Event::DefaultAlignment { run_id, alignment_id } => {
                    // Only auto-select if the user hasn't already chosen an alignment.
                    if self.selected_alignment.is_none() {
                        self.pending_alignment = Some(alignment_id);
                        self.select_run(run_id); // loads the run's alignments → applied below
                    }
                }
                Event::DonorAncestry { alignment_id, result } => {
                    self.donor_ancestry = Some((alignment_id, result));
                }
                Event::DonorPrivateY { bucket } => {
                    self.donor_private_y = Some(bucket);
                }
                Event::Alignments { sequence_run_id, alignments } => {
                    if self.selected_run == Some(sequence_run_id) {
                        self.alignments = alignments;
                        // Apply a queued subject-default alignment once its run's list is loaded.
                        if let Some(pid) = self.pending_alignment {
                            if self.alignments.iter().any(|a| a.id == pid) {
                                self.pending_alignment = None;
                                self.select_alignment(pid);
                            }
                        }
                    }
                }
                Event::AlignmentProbe(p) => {
                    // Auto-fill the add-alignment form from the BAM/CRAM header.
                    if let Some(b) = p.reference_build {
                        self.forms.aln_reference_build = b;
                    }
                    if let Some(a) = p.aligner {
                        self.forms.aln_aligner = a;
                    }
                    let bits: Vec<String> = [p.platform, p.instrument_model, p.test_type].into_iter().flatten().collect();
                    if !bits.is_empty() {
                        self.status = format!("Detected from header: {}", bits.join(" · "));
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
                Event::Sex { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.sex = result;
                    }
                    self.running_sex = false;
                    // Sex inference may have written the sex back to the biosample — reload the
                    // subjects list so the table + header reflect it instead of "Unknown".
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::ReadMetrics { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.read_metrics = result;
                    }
                    self.running_metrics = false;
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::Sv { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.sv = result;
                    }
                    self.running_sv = false;
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::Denovo { alignment_id, contig, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        match result {
                            Some(calls) => {
                                self.denovo.insert(contig, calls);
                            }
                            None => {
                                self.denovo.remove(&contig);
                            }
                        }
                    }
                    self.running_denovo = false;
                }
                Event::AnalysisProgress { step, total, label, detail, fraction } => {
                    // Reset the elapsed timer only when the step changes (sub-progress within a
                    // step keeps the same start time).
                    let started = match &self.analysis {
                        Some(a) if a.step == step => a.started,
                        _ => self.frame_time,
                    };
                    self.analysis = Some(AnalysisModal { step, total, label, detail, fraction, started });
                }
                Event::AnalysisDone { cancelled } => {
                    self.analysis = None;
                    self.status =
                        if cancelled { "Full analysis cancelled.".into() } else { "Full analysis complete.".into() };
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
                    self.painting_running = false;
                    self.running_sex = false;
                    self.running_metrics = false;
                    self.running_sv = false;
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
        self.pending_alignment = None;
        self.donor_ancestry = None;
        self.donor_private_y = None;
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
        // Subject-centric: auto-select the subject's default alignment so the analysis tabs work
        // without navigating Data Sources, and load the donor-level aggregates (best ancestry +
        // private-Y union across all sources).
        let _ = self.tx.send(Command::DefaultAlignment { biosample_guid: guid });
        let _ = self.tx.send(Command::LoadDonorAncestry { biosample_guid: guid });
        let _ = self.tx.send(Command::LoadDonorPrivateY { biosample_guid: guid });
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
        self.sex = None;
        self.read_metrics = None;
        self.sv = None;
        self.running_sex = false;
        self.running_metrics = false;
        self.running_sv = false;
        self.denovo.clear();
        self.panel_genotypes = None;
        self.ibd_result = None;
        self.identity = None;
        self.y_haplogroup = None;
        self.mt_haplogroup = None;
        self.private_y = None;
        self.ancestry = None;
        self.estimating_ancestry = false;
        self.ancestry_progress = None;
        self.pca_reference = None;
        self.painting = None;
        self.painting_running = false;
        let _ = self.tx.send(Command::LoadCoverage(id));
        let _ = self.tx.send(Command::LoadSex(id));
        let _ = self.tx.send(Command::LoadReadMetrics(id));
        let _ = self.tx.send(Command::LoadSv(id));
        let _ = self.tx.send(Command::LoadAncestry { alignment_id: id });
        let _ = self.tx.send(Command::LoadPcaReference { alignment_id: id });
        // Load cached chrM de-novo (mtDNA tab). chrY variant discovery is the masked private-Y
        // pass, not a raw whole-chrY de-novo, so it isn't loaded here.
        let _ = self.tx.send(Command::LoadDenovo { alignment_id: id, contig: "chrM".into() });
        let _ = self.tx.send(Command::LoadPrivateY { alignment_id: id }); // reload cached private-Y
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
        self.frame_time = ctx.input(|i| i.time);
        // While an analysis runs, keep repainting so the spinner/elapsed timer animate even
        // during a long step that emits no events (e.g. whole-genome coverage).
        if self.analysis.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(120));
        }
        self.drain_events();
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
                ui.separator();
                ui.label(egui::RichText::new(self.tr("status.label")).weak());
                ui.label(&self.status);
            });
        });
        if self.nav == Nav::Subjects {
            self.action_bar(ctx);
        }
        self.left_panel(ctx);
        egui::CentralPanel::default().show(ctx, |ui| match self.nav {
            Nav::Dashboard => self.dashboard_central(ui),
            Nav::Subjects => self.subjects_central(ui),
            Nav::Projects => self.projects_central(ui),
        });
        self.analysis_modal(ctx);
        self.edit_subject_modal(ctx);
        self.delete_subject_modal(ctx);
        self.data_delete_modal(ctx);
        self.assign_project_modal(ctx);
        self.edit_project_modal(ctx);
        self.delete_project_modal(ctx);
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

    /// The Projects work area: the open project's samples + coverage/haplogroup report.
    fn projects_central(&mut self, ui: &mut egui::Ui) {
        let Some(pid) = self.selected_project else {
            empty_state(ui, self.tr("empty.projects.title"), self.tr("empty.projects.hint"));
            return;
        };
        if let Some(ov) = self.overview.iter().find(|o| o.project.id == pid) {
            let proj = ov.project.clone();
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.vertical(|ui| {
                    ui.heading(&proj.name);
                    ui.label(egui::RichText::new(format!("Administrator: {}", proj.administrator)).weak());
                    if let Some(d) = &proj.description {
                        ui.label(egui::RichText::new(d).weak().small());
                    }
                });
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .add(egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE)).fill(DANGER))
                        .clicked()
                    {
                        self.confirm_delete_project = Some((proj.id, proj.name.clone()));
                    }
                    if ui.button(self.tr("common.edit")).clicked() {
                        self.edit_project = Some(EditProject {
                            id: proj.id,
                            name: proj.name.clone(),
                            description: proj.description.clone().unwrap_or_default(),
                            administrator: proj.administrator.clone(),
                        });
                    }
                });
            });
            ui.separator();
        }
        egui::ScrollArea::vertical().show(ui, |ui| {
            self.samples_section(ui);
            self.project_report_section(ui);
        });
    }

    /// The Subjects work area: the selected subject's detail — header + sub-tabs.
    fn subjects_central(&mut self, ui: &mut egui::Ui) {
        let Some(guid) = self.selected_sample else {
            empty_state(ui, self.tr("empty.subjects.title"), self.tr("empty.subjects.hint"));
            return;
        };
        self.subject_detail_header(ui, guid);
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            for (tab, key) in DetailTab::ALL {
                if ui.selectable_label(self.detail_tab == tab, egui::RichText::new(self.tr(key)).strong()).clicked() {
                    self.detail_tab = tab;
                }
            }
        });
        ui.separator();
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.add_space(4.0);
            match self.detail_tab {
                DetailTab::Overview => {
                    if let Some(id) = self.selected_alignment {
                        if ui
                            .add(egui::Button::new(egui::RichText::new(self.tr("action.runFullAnalysis")).color(egui::Color32::WHITE)).fill(ACCENT))
                            .clicked()
                        {
                            self.start_full_analysis(id);
                        }
                        ui.add_space(10.0);
                    }
                    if self.consensus_y.is_some() || self.consensus_mt.is_some() {
                        card(ui, self.tr("card.haploConsensus"), |ui| self.consensus_section(ui));
                        ui.add_space(10.0);
                    }
                    if let Some(id) = self.selected_alignment {
                        card(ui, self.tr("card.coverage"), |ui| self.coverage_section(ui, id));
                        ui.add_space(10.0);
                        card(ui, self.tr("card.sexMetrics"), |ui| self.sex_metrics_section(ui, id));
                    } else {
                        pick_alignment_hint(ui, self.tr("hint.pickAlignment"));
                    }
                }
                DetailTab::YDna => {
                    if let Some(id) = self.selected_alignment {
                        card(ui, self.tr("card.yHaplogroup"), |ui| self.y_haplogroup_section(ui, id));
                        ui.add_space(10.0);
                    } else if self.consensus_y.is_some() {
                        card(ui, self.tr("card.yHaplogroup"), |ui| self.consensus_block(ui, "Y-DNA", DnaType::Y));
                        ui.add_space(10.0);
                    } else {
                        pick_alignment_hint(ui, self.tr("hint.pickAlignment"));
                        ui.add_space(10.0);
                    }
                    if self.donor_private_y.is_some() {
                        ui.add_space(10.0);
                        card(ui, self.tr("card.privateYUnion"), |ui| self.donor_private_y_section(ui));
                    }
                    ui.add_space(10.0);
                    card(ui, self.tr("card.snpVariants"), |ui| self.variants_section(ui, guid));
                    ui.add_space(10.0);
                    card(ui, self.tr("card.ystrConsensus"), |ui| self.str_consensus_section(ui));
                }
                DetailTab::MtDna => {
                    // mtDNA haplogroup: assign standalone from the selected alignment (like Y-DNA);
                    // with no alignment selected, show the donor consensus if one was recorded.
                    if let Some(id) = self.selected_alignment {
                        card(ui, self.tr("card.mtHaplogroup"), |ui| self.mt_haplogroup_section(ui, id));
                        ui.add_space(10.0);
                    } else if self.consensus_mt.is_some() {
                        card(ui, self.tr("card.mtHaplogroup"), |ui| self.consensus_block(ui, "mtDNA", DnaType::Mt));
                        ui.add_space(10.0);
                    }
                    card(ui, self.tr("card.mtSequences"), |ui| self.mtdna_section(ui, guid));
                    if let Some(id) = self.selected_alignment {
                        ui.add_space(10.0);
                        card(ui, self.tr("card.mtDenovo"), |ui| self.denovo_section(ui, id, "chrM"));
                        ui.add_space(10.0);
                        card(ui, self.tr("card.mtHeteroplasmy"), |ui| self.heteroplasmy_section(ui, id));
                    }
                }
                DetailTab::Ancestry => {
                    if self.donor_ancestry.is_some() {
                        card(ui, self.tr("card.donorAncestry"), |ui| self.donor_ancestry_summary(ui));
                        ui.add_space(10.0);
                    }
                    if let Some(id) = self.selected_alignment {
                        card(ui, self.tr("card.ancestry"), |ui| self.ancestry_section(ui, id));
                    } else {
                        pick_alignment_hint(ui, self.tr("hint.pickAlignment"));
                    }
                }
                DetailTab::IbdMatches => {
                    if let Some(id) = self.selected_alignment {
                        card(ui, self.tr("card.panelGenotypingIbd"), |ui| self.genotyping_section(ui, id));
                    } else {
                        pick_alignment_hint(ui, self.tr("hint.pickAlignment"));
                    }
                }
                DetailTab::DataSources => self.data_sources_tab(ui, guid),
            }
        });
    }

    /// The subject-detail header: big name, ID + sex, and Add Data / Edit / Delete actions.
    fn subject_detail_header(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let bio = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|b| b.guid == guid)
            .cloned();
        let Some(bio) = bio else { return };
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.heading(&bio.donor_identifier);
                ui.label(
                    egui::RichText::new(format!("ID: {} • {}", bio.guid.0, bio.sex.as_deref().unwrap_or("Unknown"))).weak(),
                );
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE)).fill(DANGER))
                    .clicked()
                {
                    self.confirm_delete = Some(guid);
                }
                if ui.button(self.tr("common.edit")).clicked() {
                    self.edit_subject = Some(EditSubject {
                        guid,
                        donor_identifier: bio.donor_identifier.clone(),
                        sample_accession: bio.sample_accession.clone().unwrap_or_default(),
                        description: bio.description.clone().unwrap_or_default(),
                        center_name: bio.center_name.clone().unwrap_or_default(),
                        sex: bio.sex.clone().unwrap_or_default(),
                    });
                }
                if ui.button(self.tr("detail.addData")).clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("data files", &["vcf", "csv", "tsv", "txt", "fa", "fasta", "fna", "fas", "bam", "cram"])
                        .pick_file()
                    {
                        let _ = self.tx.send(Command::AddData { biosample_guid: guid, path });
                    }
                }
            });
        });
    }

    /// A simple at-a-glance dashboard: counts + account state.
    fn dashboard_central(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);
        ui.heading(self.tr("dash.title"));
        ui.add_space(12.0);
        let (projects, subjects, alignments, panels) =
            (self.overview.len(), self.all_biosamples.len(), self.all_alignments.len(), self.panels.len());
        let (lp, ls, la, lpn) =
            (self.tr("dash.projects"), self.tr("dash.subjects"), self.tr("dash.alignments"), self.tr("dash.panels"));
        ui.horizontal_wrapped(|ui| {
            stat_card(ui, lp, projects);
            stat_card(ui, ls, subjects);
            stat_card(ui, la, alignments);
            stat_card(ui, lpn, panels);
        });
        ui.add_space(16.0);
        match &self.account {
            Some(did) => {
                ui.label(format!("Signed in as {did}"));
                ui.label(if self.online { "● online" } else { "○ offline" });
            }
            None => {
                ui.label(egui::RichText::new("Not signed in — connect a PDS from the top bar to publish.").weak());
            }
        }
    }

    /// The bottom action bar for the Subjects view: selection count + batch actions.
    fn action_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("actions").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                let selected = self.selected_sample.is_some();
                ui.label(format!("{} {}", if selected { 1 } else { 0 }, self.tr("action.selected")));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.add_enabled(selected, egui::Button::new(self.tr("action.addToProject"))).clicked() {
                        if let Some(guid) = self.selected_sample {
                            let current = self
                                .all_biosamples
                                .iter()
                                .chain(self.samples.iter())
                                .find(|b| b.guid == guid)
                                .and_then(|b| b.project_id);
                            self.assign_project = Some((guid, current));
                        }
                    }
                    if ui.add_enabled(selected, egui::Button::new(self.tr("action.batchAnalyze"))).clicked() {
                        if let Some(id) = self.selected_alignment {
                            self.start_full_analysis(id);
                        } else {
                            self.status = "Select an alignment (Data Sources) to run analysis.".into();
                        }
                    }
                    // Compare needs a second subject (multi-select) — disabled for now.
                    let _ = ui.add_enabled(false, egui::Button::new(self.tr("action.compare")));
                });
            });
            ui.add_space(2.0);
        });
    }

    /// The "Full Analysis" progress modal: a dimmed backdrop + centered card with the current
    /// step, a progress bar + percent, and a Cancel button. Shown while `self.analysis` is set.
    fn analysis_modal(&mut self, ctx: &egui::Context) {
        let Some(p) = self.analysis.clone() else { return };
        // Dim everything behind the dialog.
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("modal_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        egui::Area::new(egui::Id::new("analysis_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(460.0);
                    ui.label(egui::RichText::new("Full Analysis").strong().size(16.0));
                    ui.separator();
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label(egui::RichText::new("Analysis in progress…").weak());
                        let elapsed = (self.frame_time - p.started).max(0.0) as u64;
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(egui::RichText::new(format!("{:02}:{:02}", elapsed / 60, elapsed % 60)).weak());
                        });
                    });
                    ui.add_space(10.0);
                    ui.label(format!("Step {}/{}: {} — {}", p.step, p.total, p.label, p.detail));
                    ui.add_space(10.0);
                    ui.horizontal(|ui| {
                        // `animate` shimmers the bar so a long step reads as working, not stalled.
                        ui.add(egui::ProgressBar::new(p.fraction).desired_width(360.0).rounding(4.0).animate(true));
                        ui.label(format!("{}%", (p.fraction * 100.0).round() as i32));
                    });
                    ui.add_space(6.0);
                    ui.label(
                        egui::RichText::new("Whole-genome steps (coverage) can take several minutes on a WGS BAM.")
                            .weak()
                            .small(),
                    );
                    ui.add_space(12.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button(self.tr("common.cancel")).clicked() {
                            let _ = self.tx.send(Command::CancelAnalysis);
                            self.status = "Cancelling…".into();
                        }
                    });
                });
            });
    }

    /// The Edit-subject modal: editable fields over a dimmed backdrop. Save sends an
    /// `UpdateBiosample` command; the resulting `BiosamplesChanged` event refreshes the lists.
    fn edit_subject_modal(&mut self, ctx: &egui::Context) {
        let Some(mut edit) = self.edit_subject.clone() else { return };
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("edit_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        let mut close = false;
        egui::Area::new(egui::Id::new("edit_subject_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(420.0);
                    ui.label(egui::RichText::new(self.tr("edit.title")).strong().size(16.0));
                    ui.separator();
                    ui.add_space(6.0);
                    let field = |ui: &mut egui::Ui, label: &str, value: &mut String, hint: &str| {
                        ui.label(label);
                        ui.add(egui::TextEdit::singleline(value).hint_text(hint).desired_width(f32::INFINITY));
                        ui.add_space(4.0);
                    };
                    field(ui, self.tr("edit.identifier"), &mut edit.donor_identifier, "donor identifier");
                    field(ui, self.tr("edit.accession"), &mut edit.sample_accession, "accession (optional)");
                    field(ui, self.tr("edit.description"), &mut edit.description, "description (optional)");
                    field(ui, self.tr("edit.center"), &mut edit.center_name, "center (optional)");
                    field(ui, self.tr("edit.sex"), &mut edit.sex, "sex (optional)");
                    ui.add_space(10.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add_enabled(
                                !edit.donor_identifier.trim().is_empty(),
                                egui::Button::new(self.tr("common.save")).fill(ACCENT),
                            )
                            .clicked()
                        {
                            let _ = self.tx.send(Command::UpdateBiosample {
                                guid: edit.guid,
                                donor_identifier: edit.donor_identifier.trim().to_string(),
                                sample_accession: opt(&edit.sample_accession),
                                description: opt(&edit.description),
                                center_name: opt(&edit.center_name),
                                sex: opt(&edit.sex),
                            });
                            close = true;
                        }
                        if ui.button(self.tr("common.cancel")).clicked() {
                            close = true;
                        }
                    });
                });
            });
        if close {
            self.edit_subject = None;
        } else {
            self.edit_subject = Some(edit);
        }
    }

    /// The Delete-subject confirmation modal. Confirm sends a `DeleteBiosample` command; the app
    /// layer refuses (surfaced via the status bar) when the subject still has dependent data.
    fn delete_subject_modal(&mut self, ctx: &egui::Context) {
        let Some(guid) = self.confirm_delete else { return };
        let name = self
            .all_biosamples
            .iter()
            .chain(self.samples.iter())
            .find(|b| b.guid == guid)
            .map(|b| b.donor_identifier.clone())
            .unwrap_or_else(|| guid.0.to_string());
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("delete_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        let mut close = false;
        egui::Area::new(egui::Id::new("delete_subject_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(400.0);
                    ui.label(egui::RichText::new(self.tr("delete.title")).strong().size(16.0));
                    ui.separator();
                    ui.add_space(8.0);
                    ui.label(format!("{} “{}”?", self.tr("delete.confirm"), name));
                    ui.add_space(6.0);
                    ui.label(egui::RichText::new(self.tr("delete.note")).weak().small());
                    ui.add_space(12.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add(egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE)).fill(DANGER))
                            .clicked()
                        {
                            let _ = self.tx.send(Command::DeleteBiosample(guid));
                            if self.selected_sample == Some(guid) {
                                self.selected_sample = None;
                            }
                            close = true;
                        }
                        if ui.button(self.tr("common.cancel")).clicked() {
                            close = true;
                        }
                    });
                });
            });
        if close {
            self.confirm_delete = None;
        }
    }

    /// Confirmation modal for deleting a data-source row (run/alignment/profile). Confirm sends
    /// the variant's worker command; the resulting change event refreshes the affected list.
    fn data_delete_modal(&mut self, ctx: &egui::Context) {
        let Some(target) = self.confirm_data_delete.clone() else { return };
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("data_delete_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        let mut close = false;
        egui::Area::new(egui::Id::new("data_delete_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(400.0);
                    ui.label(egui::RichText::new(self.tr("delete.dataTitle")).strong().size(16.0));
                    ui.separator();
                    ui.add_space(8.0);
                    ui.label(format!("{} {}?", self.tr("delete.confirm"), target.label()));
                    ui.add_space(6.0);
                    ui.label(egui::RichText::new(self.tr("delete.dataNote")).weak().small());
                    ui.add_space(12.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add(egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE)).fill(DANGER))
                            .clicked()
                        {
                            let _ = self.tx.send(target.command());
                            close = true;
                        }
                        if ui.button(self.tr("common.cancel")).clicked() {
                            close = true;
                        }
                    });
                });
            });
        if close {
            self.confirm_data_delete = None;
        }
    }

    /// The Add-to-Project picker: a dropdown of projects (plus "no project"). Save sends
    /// `AssignBiosampleProject`; the resulting `BiosamplesChanged` event refreshes the lists.
    fn assign_project_modal(&mut self, ctx: &egui::Context) {
        let Some((guid, mut chosen)) = self.assign_project else { return };
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("assign_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        let selected_text = match chosen {
            Some(pid) => self
                .overview
                .iter()
                .find(|o| o.project.id == pid)
                .map(|o| o.project.name.clone())
                .unwrap_or_else(|| format!("project {pid}")),
            None => self.tr("projects.noProject").to_string(),
        };
        let mut close = false;
        let mut commit = false;
        egui::Area::new(egui::Id::new("assign_project_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(360.0);
                    ui.label(egui::RichText::new(self.tr("action.addToProject")).strong().size(16.0));
                    ui.separator();
                    ui.add_space(8.0);
                    egui::ComboBox::from_id_salt("assign_project_combo")
                        .selected_text(selected_text)
                        .width(300.0)
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut chosen, None, self.tr("projects.noProject"));
                            for o in &self.overview {
                                ui.selectable_value(&mut chosen, Some(o.project.id), &o.project.name);
                            }
                        });
                    ui.add_space(12.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.add(egui::Button::new(self.tr("common.save")).fill(ACCENT)).clicked() {
                            commit = true;
                            close = true;
                        }
                        if ui.button(self.tr("common.cancel")).clicked() {
                            close = true;
                        }
                    });
                });
            });
        if commit {
            let _ = self.tx.send(Command::AssignBiosampleProject { guid, project_id: chosen });
        }
        if close {
            self.assign_project = None;
        } else {
            self.assign_project = Some((guid, chosen));
        }
    }

    /// The Edit-project modal: name / administrator / description. Save sends `UpdateProject`.
    fn edit_project_modal(&mut self, ctx: &egui::Context) {
        let Some(mut edit) = self.edit_project.clone() else { return };
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("edit_proj_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        let mut close = false;
        egui::Area::new(egui::Id::new("edit_project_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(400.0);
                    ui.label(egui::RichText::new(self.tr("editProject.title")).strong().size(16.0));
                    ui.separator();
                    ui.add_space(6.0);
                    ui.label(self.tr("editProject.name"));
                    ui.add(egui::TextEdit::singleline(&mut edit.name).hint_text("name").desired_width(f32::INFINITY));
                    ui.add_space(4.0);
                    ui.label(self.tr("editProject.admin"));
                    ui.add(egui::TextEdit::singleline(&mut edit.administrator).hint_text("administrator").desired_width(f32::INFINITY));
                    ui.add_space(4.0);
                    ui.label(self.tr("editProject.description"));
                    ui.add(egui::TextEdit::multiline(&mut edit.description).hint_text("description (optional)").desired_width(f32::INFINITY));
                    ui.add_space(10.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add_enabled(!edit.name.trim().is_empty(), egui::Button::new(self.tr("common.save")).fill(ACCENT))
                            .clicked()
                        {
                            let _ = self.tx.send(Command::UpdateProject {
                                id: edit.id,
                                name: edit.name.trim().to_string(),
                                description: opt(&edit.description),
                                administrator: opt(&edit.administrator).unwrap_or_else(|| "unknown".into()),
                            });
                            close = true;
                        }
                        if ui.button(self.tr("common.cancel")).clicked() {
                            close = true;
                        }
                    });
                });
            });
        if close {
            self.edit_project = None;
        } else {
            self.edit_project = Some(edit);
        }
    }

    /// The Delete-project confirmation modal. Confirm sends `DeleteProject`; the app layer
    /// refuses (surfaced via the status bar) while subjects still belong to the project.
    fn delete_project_modal(&mut self, ctx: &egui::Context) {
        let Some((id, name)) = self.confirm_delete_project.clone() else { return };
        let painter = ctx.layer_painter(egui::LayerId::new(egui::Order::Foreground, egui::Id::new("del_proj_dim")));
        painter.rect_filled(ctx.screen_rect(), 0.0, egui::Color32::from_black_alpha(150));

        let mut close = false;
        egui::Area::new(egui::Id::new("delete_project_modal"))
            .order(egui::Order::Foreground)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                egui::Frame::window(ui.style()).inner_margin(egui::Margin::same(18.0)).show(ui, |ui| {
                    ui.set_width(400.0);
                    ui.label(egui::RichText::new(self.tr("editProject.deleteTitle")).strong().size(16.0));
                    ui.separator();
                    ui.add_space(8.0);
                    ui.label(format!("{} “{}”?", self.tr("delete.confirm"), name));
                    ui.add_space(6.0);
                    ui.label(egui::RichText::new(self.tr("editProject.deleteNote")).weak().small());
                    ui.add_space(12.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add(egui::Button::new(egui::RichText::new(self.tr("common.delete")).color(egui::Color32::WHITE)).fill(DANGER))
                            .clicked()
                        {
                            let _ = self.tx.send(Command::DeleteProject(id));
                            if self.selected_project == Some(id) {
                                self.selected_project = None;
                            }
                            close = true;
                        }
                        if ui.button(self.tr("common.cancel")).clicked() {
                            close = true;
                        }
                    });
                });
            });
        if close {
            self.confirm_delete_project = None;
        }
    }

    /// Kick off the full-analysis pipeline for an alignment and show the modal immediately.
    fn start_full_analysis(&mut self, alignment_id: i64) {
        self.analysis = Some(AnalysisModal {
            step: 1,
            total: 8,
            label: "Starting".into(),
            detail: "preparing pipeline".into(),
            fraction: 0.0,
            started: self.frame_time,
        });
        self.status = format!("Running full analysis on alignment #{alignment_id}…");
        let _ = self.tx.send(Command::RunFullAnalysis { alignment_id });
    }

    /// Translate a catalog key for the active language. Returns `&'static str` (catalogs are
    /// embedded), so it never borrows `self` — convenient inside egui closures.
    fn tr(&self, key: &'static str) -> &'static str {
        crate::i18n::tr(self.lang, key)
    }

    /// The top app bar: product title (left), theme + language + account controls (right).
    fn app_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("appbar").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                ui.heading(self.tr("app.name"));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let icon = if self.dark_mode { "☀" } else { "🌙" };
                    if ui.button(icon).on_hover_text(self.tr("theme.toggle")).clicked() {
                        self.dark_mode = !self.dark_mode;
                        apply_theme(ctx, self.dark_mode);
                    }
                    ui.separator();
                    let prev_lang = self.lang;
                    egui::ComboBox::from_id_salt("lang")
                        .selected_text(self.lang.label())
                        .show_ui(ui, |ui| {
                            for &l in crate::i18n::Lang::all() {
                                ui.selectable_value(&mut self.lang, l, l.label());
                            }
                        });
                    if self.lang != prev_lang {
                        crate::i18n::save_lang(self.lang); // persist across restarts
                    }
                    ui.separator();
                    self.account_controls(ui);
                });
            });
            ui.add_space(2.0);
        });
    }

    /// The primary navigation strip (Dashboard / Subjects / Projects).
    fn nav_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("nav").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                for (nav, icon, key) in [
                    (Nav::Dashboard, "📊", "nav.dashboard"),
                    (Nav::Subjects, "👥", "nav.subjects"),
                    (Nav::Projects, "📁", "nav.projects"),
                ] {
                    let label = format!("{icon}  {}", self.tr(key));
                    if ui.selectable_label(self.nav == nav, egui::RichText::new(label).strong()).clicked() {
                        self.nav = nav;
                    }
                }
            });
            ui.add_space(2.0);
        });
    }

    /// Account / sign-in controls (reused in the app bar).
    fn account_controls(&mut self, ui: &mut egui::Ui) {
        match &self.account {
            Some(did) => {
                let did = did.clone();
                if ui.button(self.tr("account.signOut")).clicked() {
                    let _ = self.tx.send(Command::Logout);
                }
                if self.online {
                    ui.colored_label(egui::Color32::from_rgb(80, 190, 120), self.tr("account.online"));
                } else {
                    ui.colored_label(egui::Color32::from_rgb(220, 150, 60), self.tr("account.offline"));
                }
                let short: String = did.chars().take(22).collect();
                ui.label(egui::RichText::new(short).weak());
            }
            None => {
                if self.logging_in {
                    ui.spinner();
                }
                let ready = !self.forms.login_handle.trim().is_empty() && !self.logging_in;
                if ui.add_enabled(ready, egui::Button::new(self.tr("account.signIn"))).clicked() {
                    self.logging_in = true;
                    self.status = "Opening browser to authorize…".into();
                    let _ = self.tx.send(Command::Login { handle: self.forms.login_handle.trim().to_string() });
                }
                let hint = self.tr("account.handleHint");
                ui.add(
                    egui::TextEdit::singleline(&mut self.forms.login_handle)
                        .hint_text(hint)
                        .desired_width(180.0),
                );
                ui.label(self.tr("account.pds"));
            }
        }
    }

    /// The left panel, routed by the active nav tab. Hidden on the Dashboard.
    fn left_panel(&mut self, ctx: &egui::Context) {
        match self.nav {
            Nav::Dashboard => {}
            Nav::Projects => {
                egui::SidePanel::left("left").min_width(240.0).show(ctx, |ui| self.projects_side(ui));
            }
            Nav::Subjects => {
                egui::SidePanel::left("left")
                    .resizable(true)
                    .default_width(680.0)
                    .min_width(420.0)
                    .show(ctx, |ui| self.subjects_side(ui));
            }
        }
    }

    /// Project management: the projects list, new-project form, batch import, and panels.
    fn projects_side(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.heading(self.tr("projects.heading"));
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
        ui.label(self.tr("projects.new"));
        ui.add(egui::TextEdit::singleline(&mut self.forms.project_name).hint_text("name"));
        ui.add(egui::TextEdit::singleline(&mut self.forms.project_admin).hint_text("administrator"));
        if ui
            .add_enabled(!self.forms.project_name.trim().is_empty(), egui::Button::new(self.tr("projects.create")))
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
        ui.label(self.tr("projects.batchImportHint"));
        if ui.add_enabled(!self.importing, egui::Button::new(self.tr("projects.batchImport"))).clicked() {
            if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                self.importing = true;
                self.reference_needs.clear();
                self.status = format!("Importing {}…", dir.display());
                let _ = self.tx.send(Command::ImportProjectDir { dir, reference: None });
            }
        }
        self.reference_prompt(ui);

        ui.add_space(12.0);
        ui.separator();
        self.panels_section(ui);
    }

    /// Subjects browser: a search box + "Add New Subject" on one row, then the subjects table.
    fn subjects_side(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            let btn_w = 160.0;
            let hint = self.tr("subjects.search");
            ui.add(
                egui::TextEdit::singleline(&mut self.subject_search)
                    .hint_text(hint)
                    .desired_width((ui.available_width() - btn_w).max(120.0)),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(egui::Button::new(egui::RichText::new(self.tr("subjects.addNew")).color(egui::Color32::WHITE)).fill(ACCENT))
                    .clicked()
                {
                    self.forms.show_add_subject = !self.forms.show_add_subject;
                }
            });
        });
        if self.forms.show_add_subject {
            self.add_subject_form(ui);
        }
        ui.add_space(6.0);
        self.subjects_table(ui);
    }

    /// All biosamples, independent of any project (the project link is optional). Selecting
    /// one drives the runs → alignments → analysis flow in the central panel; adding one
    /// tags it to the open project if there is one, else leaves it project-less.
    /// The subjects table: columns ID / Name / Y-DNA / mtDNA / Sex / Center / Status, with the
    /// selected row highlighted. Clicking a row selects the subject.
    fn subjects_table(&mut self, ui: &mut egui::Ui) {
        if self.all_biosamples.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.none")).weak());
            return;
        }
        let total_w: f32 = SUBJECT_COLS.iter().map(|c| c.1).sum();
        let needle = self.subject_search.trim().to_lowercase();
        let mut pick = None;
        let mut shown = 0;
        egui::ScrollArea::both().show(ui, |ui| {
            ui.set_min_width(total_w);
            table_header(ui, &SUBJECT_COLS);
            for s in &self.all_biosamples {
                if !needle.is_empty() {
                    let hay = format!(
                        "{} {}",
                        s.donor_identifier.to_lowercase(),
                        s.sample_accession.as_deref().unwrap_or("").to_lowercase()
                    );
                    if !hay.contains(&needle) {
                        continue;
                    }
                }
                shown += 1;
                // Y/mt haplogroups are only loaded for the selected subject (consensus_* state),
                // so fill them in for that row; others stay "-" until selected/analyzed.
                let (y, mt) = if self.selected_sample == Some(s.guid) {
                    (
                        self.consensus_y.as_ref().map(|c| c.haplogroup.clone()).unwrap_or_else(|| "-".into()),
                        self.consensus_mt.as_ref().map(|c| c.haplogroup.clone()).unwrap_or_else(|| "-".into()),
                    )
                } else {
                    ("-".into(), "-".into())
                };
                let cells = [
                    short_guid(s),
                    s.donor_identifier.clone(),
                    y,
                    mt,
                    s.sex.clone().unwrap_or_else(|| "-".into()),
                    s.center_name.clone().unwrap_or_else(|| "-".into()),
                    "Pending".to_string(),
                ];
                if table_row(ui, &SUBJECT_COLS, &cells, self.selected_sample == Some(s.guid), Some(6)) {
                    pick = Some(s.guid);
                }
            }
            if shown == 0 {
                ui.label(egui::RichText::new(self.tr("subjects.noMatch")).weak());
            }
        });
        if let Some(guid) = pick {
            self.select_sample(guid);
        }
    }

    /// The inline new-subject form (toggled by the "Add New Subject" button).
    fn add_subject_form(&mut self, ui: &mut egui::Ui) {
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_donor).hint_text("donor identifier"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_accession).hint_text("accession (optional)"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.sample_sex).hint_text("sex (optional)"));
            if let Some(pid) = self.selected_project {
                let name = self.overview.iter().find(|o| o.project.id == pid).map(|o| o.project.name.as_str());
                ui.label(format!("→ project: {}", name.unwrap_or("(open)")));
            } else {
                ui.label(self.tr("projects.noProject"));
            }
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(!self.forms.sample_donor.trim().is_empty(), egui::Button::new(self.tr("projects.addSubject")).fill(ACCENT))
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
                    self.forms.show_add_subject = false;
                }
                if ui.button(self.tr("common.cancel")).clicked() {
                    self.forms.show_add_subject = false;
                }
            });
        });
    }

    /// Y-STR profiles for the selected subject + an import form (CSV/TSV marker table).
    /// Donor-level Y-STR consensus across all of the subject's panels (Phase 2 rollup): the modal
    /// value per marker, with cross-panel disagreements flagged.
    fn str_consensus_section(&mut self, ui: &mut egui::Ui) {
        if self.str_profiles.is_empty() {
            ui.label(egui::RichText::new("No STR profiles yet — import one under Data Sources.").weak());
            return;
        }
        let consensus = strprofile::consensus_markers(&self.str_profiles);
        ui.label(
            egui::RichText::new(format!("{} markers from {} panel(s)", consensus.len(), self.str_profiles.len()))
                .weak(),
        );
        let conflicts = consensus.iter().filter(|m| m.conflict).count();
        if conflicts > 0 {
            ui.colored_label(
                egui::Color32::from_rgb(220, 150, 60),
                format!("⚠ {conflicts} marker(s) disagree across panels"),
            );
        }
        egui::Grid::new("str_consensus").striped(true).num_columns(3).show(ui, |ui| {
            ui.strong(self.tr("table.marker"));
            ui.strong(self.tr("table.value"));
            ui.strong(self.tr("table.panels"));
            ui.end_row();
            for m in &consensus {
                ui.label(&m.marker);
                if m.conflict {
                    ui.colored_label(egui::Color32::from_rgb(220, 150, 60), &m.value);
                } else {
                    ui.label(&m.value);
                }
                ui.label(m.panels.to_string());
                ui.end_row();
            }
        });
    }

    /// Donor-level ancestry headline (Phase 3): the best estimate across the subject's sources,
    /// with which source + method it came from.
    fn donor_ancestry_summary(&self, ui: &mut egui::Ui) {
        let Some((aln, r)) = &self.donor_ancestry else {
            ui.label(egui::RichText::new("No ancestry estimate for any source yet.").weak());
            return;
        };
        ui.horizontal(|ui| {
            draw_ancestry_donut(ui, &r.super_population_summary);
            ui.add_space(8.0);
            ui.vertical(|ui| {
                if let Some(top) = r.super_population_summary.first() {
                    ui.heading(format!("{} {:.1}%", top.super_population, top.percentage));
                }
                ui.label(format!(
                    "{}/{} SNPs · confidence {:.0}%",
                    r.snps_with_genotype,
                    r.snps_analyzed,
                    r.confidence_level * 100.0
                ));
                ui.label(
                    egui::RichText::new(format!("best source: alignment #{aln} · {} · {}", r.method, r.reference_version))
                        .small()
                        .weak(),
                );
                ui.add_space(4.0);
                draw_composition_bar(ui, &r.super_population_summary);
            });
        });
    }

    /// Donor-level private-Y union (Phase 3): off-backbone calls pooled + deduped across the
    /// subject's Y-bearing sources.
    fn donor_private_y_section(&self, ui: &mut egui::Ui) {
        let Some(bucket) = &self.donor_private_y else {
            ui.label(egui::RichText::new("No private-Y calls across sources yet — run \"Find private Y variants\".").weak());
            return;
        };
        ui.label(format!(
            "{} novel + {} off-path  (union across sources, terminal {})",
            bucket.novel(),
            bucket.off_path(),
            bucket.terminal
        ));
        egui::Grid::new("donor_privy").striped(true).num_columns(4).show(ui, |ui| {
            for h in ["table.position", "table.change", "table.depth", "table.class"] {
                ui.strong(self.tr(h));
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
    }

    fn str_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let mut want_delete: Option<DataDelete> = None;
        for p in &self.str_profiles {
            let provider = p.provider.as_deref().unwrap_or("—");
            let header = format!("{} — {} markers  ({provider})", p.panel_name, p.markers.len());
            egui::CollapsingHeader::new(header).id_salt(("str", p.id)).show(ui, |ui| {
                if ui.small_button(self.tr("delete.thisProfile")).clicked() {
                    want_delete = Some(DataDelete::Str { id: p.id, guid, label: format!("STR profile “{}”", p.panel_name) });
                }
                egui::Grid::new(("str_markers", p.id)).striped(true).num_columns(2).show(ui, |ui| {
                    ui.strong(self.tr("table.marker"));
                    ui.strong(self.tr("table.value"));
                    ui.end_row();
                    for m in &p.markers {
                        ui.label(&m.marker);
                        ui.label(&m.value);
                        ui.end_row();
                    }
                });
            });
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("str.import"), |ui| {
            // Bind labels first — `self.tr()` (immutable) can't share the statement with the
            // `&mut self.forms.*` below (the i18n borrow gotcha).
            let (panel_lbl, provider_lbl, source_lbl) =
                (self.tr("form.panel"), self.tr("form.provider"), self.tr("form.source"));
            combo(ui, panel_lbl, "str_panel", &mut self.forms.str_panel, strprofile::KNOWN_PANELS);
            combo(ui, provider_lbl, "str_provider", &mut self.forms.str_provider, strprofile::KNOWN_PROVIDERS);
            combo(ui, source_lbl, "str_source", &mut self.forms.str_source, strprofile::KNOWN_SOURCES);
            if ui.button(self.tr("str.chooseCsv")).clicked() {
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
            ui.label(self.tr("str.expectsRows"));
        });
    }

    /// SNP variant sets for the selected subject + an import form (VCF or CSV/TSV).
    fn variants_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.variant_sets.is_empty() {
            ui.label(egui::RichText::new("No variants imported yet.").weak());
        }
        const MAX_ROWS: usize = 500;
        let mut want_delete: Option<DataDelete> = None;
        for s in &self.variant_sets {
            let header = format!("{} — {} call(s)", s.source_label, s.calls.len());
            egui::CollapsingHeader::new(header).id_salt(("vset", s.id)).show(ui, |ui| {
                if ui.small_button(self.tr("delete.thisProfile")).clicked() {
                    want_delete = Some(DataDelete::Variant { id: s.id, guid, label: format!("variant set “{}”", s.source_label) });
                }
                egui::Grid::new(("vcalls", s.id)).striped(true).num_columns(4).show(ui, |ui| {
                    for h in ["table.position", "table.change", "table.rsid", "table.genotype"] {
                        ui.strong(self.tr(h));
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
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
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
                        for h in ["table.position", "table.allele", "table.status", "table.sources"] {
                            ui.strong(self.tr(h));
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
        ui.collapsing(self.tr("variants.import"), |ui| {
            let labels: Vec<&str> = SourceType::ALL.iter().map(|t| t.as_str()).collect();
            let source_lbl = self.tr("form.source");
            combo(ui, source_lbl, "variant_source", &mut self.forms.variant_source_type, &labels);
            let source_type = SourceType::from_code(&self.forms.variant_source_type);

            if ui.button(self.tr("chip.import")).clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("variants", &["vcf", "csv", "tsv", "txt"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportVariants { biosample_guid: guid, path, source_type });
                }
            }
            ui.label(self.tr("chip.formatHint"));

            ui.separator();
            ui.label(self.tr("str.pasteCalls"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.variant_manual_label).hint_text("source label (e.g. YSEQ panel)"));
            ui.add(
                egui::TextEdit::multiline(&mut self.forms.variant_manual_text)
                    .hint_text("contig,position,ref,alt per line")
                    .desired_rows(3),
            );
            let ready = !self.forms.variant_manual_text.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new(self.tr("str.addPasted"))).clicked() {
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
        let mut want_delete: Option<DataDelete> = None;
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
                if ui.small_button(self.tr("delete.thisProfile")).clicked() {
                    want_delete = Some(DataDelete::Chip { id: p.id, guid, label: format!("chip profile ({})", p.provider) });
                }
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
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("chip.section"), |ui| {
            ui.horizontal(|ui| {
                ui.label(self.tr("form.provider"));
                egui::ComboBox::from_id_salt("chip_provider")
                    .selected_text(self.forms.chip_provider.clone())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.forms.chip_provider, AUTO_DETECT.to_string(), AUTO_DETECT);
                        for p in chipprofile::KNOWN_PROVIDERS {
                            ui.selectable_value(&mut self.forms.chip_provider, p.to_string(), *p);
                        }
                    });
            });
            if ui.button(self.tr("chip.chooseCsv")).clicked() {
                if let Some(path) = rfd::FileDialog::new().add_filter("array data", &["csv", "txt", "tsv"]).pick_file() {
                    let provider = (self.forms.chip_provider != AUTO_DETECT).then(|| self.forms.chip_provider.clone());
                    let _ = self.tx.send(Command::ImportChipProfile { biosample_guid: guid, provider, path });
                }
            }
            ui.label(self.tr("chip.rawHint"));
        });
    }

    /// mtDNA FASTA sequences for the selected subject + an import form, and a
    /// derive-variants-vs-rCRS action per sequence.
    fn mtdna_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.mtdna_sequences.is_empty() {
            ui.label(egui::RichText::new("No mtDNA sequences yet.").weak());
        }

        // rCRS reference picker (reused for every derivation this session).
        ui.horizontal(|ui| {
            ui.label(self.tr("mt.rcrsRef"));
            let label = self
                .rcrs_path
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "not set".into());
            ui.label(label);
            if ui.button(self.tr("mt.chooseRcrs")).clicked() {
                if let Some(p) = rfd::FileDialog::new().add_filter("FASTA", &["fa", "fasta", "fna", "fas"]).pick_file() {
                    self.rcrs_path = Some(p);
                }
            }
        });

        let rcrs = self.rcrs_path.clone();
        // Bind before the &self loop borrow — used inside the per-row closure.
        let assign_lbl = self.tr("common.assignHaplogroup");
        let derive_lbl = self.tr("mt.deriveVariants");
        let delete_lbl = self.tr("common.delete");
        let mut want_delete: Option<DataDelete> = None;
        for m in &self.mtdna_sequences {
            let name = m.source_file_name.as_deref().or(m.defline.as_deref()).unwrap_or("mtDNA");
            ui.horizontal(|ui| {
                ui.label(format!("{name} — {} bp, {} N", m.length(), m.n_count));
                if ui
                    .add_enabled(rcrs.is_some(), egui::Button::new(derive_lbl))
                    .clicked()
                {
                    if let Some(path) = rcrs.clone() {
                        let _ = self.tx.send(Command::DeriveMtdnaVariants { mtdna_id: m.id, rcrs_path: path });
                    }
                }
                if ui.button(assign_lbl).clicked() {
                    self.status = "Assigning haplogroup (fetching FTDNA tree)…".into();
                    let _ = self.tx.send(Command::AssignMtdnaHaplogroup { mtdna_id: m.id });
                }
                if ui.button(delete_lbl).clicked() {
                    want_delete = Some(DataDelete::Mtdna { id: m.id, guid, label: format!("mtDNA sequence “{name}”") });
                }
            });
            // Show the haplogroup result for this sequence, if any.
            if let Some((id, assignment)) = &self.mtdna_haplogroup {
                if *id == m.id {
                    show_assignment(ui, assignment);
                }
            }
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }

        ui.add_space(6.0);
        ui.collapsing(self.tr("mt.importFasta"), |ui| {
            if ui.button(self.tr("mt.chooseFasta")).clicked() {
                if let Some(path) =
                    rfd::FileDialog::new().add_filter("FASTA", &["fa", "fasta", "fna", "fas"]).pick_file()
                {
                    let _ = self.tx.send(Command::ImportMtdna { biosample_guid: guid, path });
                }
            }
            ui.label(self.tr("mt.fullSeq"));
        });
    }

    fn panels_section(&mut self, ui: &mut egui::Ui) {
        ui.label(self.tr("table.panels"));
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
            .add_enabled(!self.forms.panel_import_name.trim().is_empty(), egui::Button::new(self.tr("mt.importSitesVcf")))
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
                ui.label(self.tr("refdl.required"));
                for b in &self.reference_needs {
                    ui.label(format!("  • {} (~{} MB)", b.build, b.est_bytes / 1_000_000));
                }
                if ui
                    .add_enabled(self.reference_progress.is_none(), egui::Button::new(self.tr("common.downloadContinue")))
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
            ui.heading(self.tr("projects.report"));
            let busy = self.analyzing || self.running;
            if ui.add_enabled(!busy, egui::Button::new(self.tr("projects.analyzeAll"))).clicked() {
                if let Some(pid) = self.selected_project {
                    self.analyzing = true;
                    self.status = "Analyzing project (coverage + Y per sample)…".into();
                    let _ = self.tx.send(Command::AnalyzeProject(pid));
                }
            }
            if ui.button(self.tr("projects.exportCsv")).clicked() {
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
        egui::Grid::new("project_report_grid").striped(true).num_columns(15).show(ui, |ui| {
            for h in ["report.sample", "report.alns", "report.meanCov", "report.median", "report.cov10x", "report.cov20x", "report.callable", "report.y", "report.mtdna", "report.sex", "report.readLen", "report.pctAln", "report.insert", "report.sv", "report.actions"] {
                ui.strong(self.tr(h));
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
                ui.label(r.sex.clone().unwrap_or_else(|| "—".into()));
                ui.label(fmt_depth(r.mean_read_length));
                ui.label(fmt_pct(r.pct_aligned));
                ui.label(fmt_depth(r.median_insert_size));
                ui.label(r.sv_count.map(|v| v.to_string()).unwrap_or_else(|| "—".into()));
                if let Some(aln) = r.primary_alignment_id {
                    ui.horizontal(|ui| {
                        if ui.add_enabled(!running, egui::Button::new(self.tr("btn.cov"))).clicked() {
                            recompute = Some(aln);
                        }
                        if ui.add_enabled(!running, egui::Button::new(self.tr("report.y"))).clicked() {
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
            ui.label(self.tr("projects.noSamples"));
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
        ui.label(self.tr("projects.addSubjectsHint"));
    }

    /// The Data Sources tab: sequencing runs (cards with expandable alignments), chip/array,
    /// and STR profiles — each in a rounded card.
    fn data_sources_tab(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.add_space(4.0);
        card(ui, self.tr("card.sequencingRuns"), |ui| self.runs_card(ui, guid));
        ui.add_space(10.0);
        card(ui, self.tr("card.chipProfiles"), |ui| {
            if self.chip_profiles.is_empty() {
                ui.label(egui::RichText::new("No chip/array data").weak());
            }
            self.chip_section(ui, guid);
        });
        ui.add_space(10.0);
        card(ui, self.tr("card.strProfiles"), |ui| {
            if self.str_profiles.is_empty() {
                ui.label(egui::RichText::new("No STR profiles").weak());
            }
            self.str_section(ui, guid);
        });
    }

    /// The sequencing-runs body: one card per run (provider chip, title, read meta, Y/mt
    /// badges); the selected run expands to its alignment rows + the add-alignment form.
    fn runs_card(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.runs.is_empty() {
            ui.label(egui::RichText::new("No sequencing runs yet.").weak());
        }
        // Clone the small lists so we can call &mut self methods (add forms) inside the loop.
        let runs = self.runs.clone();
        let alignments = self.alignments.clone();
        let coverage = self.coverage.clone();
        let mut pick_run = None;
        let mut pick_aln = None;
        let mut want_delete: Option<DataDelete> = None;

        for r in &runs {
            let selected = self.selected_run == Some(r.id);
            let frame = egui::Frame::group(ui.style())
                .fill(if selected { ACCENT.gamma_multiply(0.18) } else { ui.visuals().extreme_bg_color })
                .stroke(if selected { egui::Stroke::new(1.0, ACCENT) } else { egui::Stroke::NONE })
                .rounding(egui::Rounding::same(6.0))
                .inner_margin(egui::Margin::same(10.0));
            let resp = frame
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        chip(ui, &provider_abbrev(&r.platform_name), ACCENT.gamma_multiply(0.3), ACCENT);
                        ui.add_space(4.0);
                        let tt = testtype::by_code(&r.test_type);
                        ui.vertical(|ui| {
                            let plat = if r.platform_name.is_empty() { "—" } else { r.platform_name.as_str() };
                            let title = format!(
                                "{}  ·  {}  ·  {}",
                                testtype::display_name(&r.test_type),
                                plat,
                                r.instrument_model.as_deref().unwrap_or("—")
                            );
                            ui.label(egui::RichText::new(title).strong());
                            ui.label(
                                egui::RichText::new(format!(
                                    "Reads: {}   Aligned: {}   {}",
                                    fmt_reads(r.total_reads),
                                    fmt_reads(r.pf_reads_aligned),
                                    r.library_layout.as_deref().unwrap_or("SINGLE")
                                ))
                                .weak()
                                .small(),
                            );
                        });
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.small_button("🗑").on_hover_text("Delete run + its alignments").clicked() {
                                want_delete = Some(DataDelete::Run {
                                    id: r.id,
                                    guid,
                                    label: format!("run “{}”", testtype::display_name(&r.test_type)),
                                });
                            }
                            if let Some(t) = tt {
                                let mt = matches!(t.target, testtype::TargetType::WholeGenome | testtype::TargetType::MtDna);
                                let y = matches!(t.target, testtype::TargetType::WholeGenome | testtype::TargetType::YChromosome);
                                if mt {
                                    chip(ui, "mt", egui::Color32::from_rgb(70, 60, 90), egui::Color32::from_rgb(200, 180, 230));
                                }
                                if y {
                                    chip(ui, "Y", egui::Color32::from_rgb(40, 70, 55), egui::Color32::from_rgb(150, 220, 180));
                                }
                            }
                        });
                    });
                })
                .response;
            if resp.interact(egui::Sense::click()).clicked() {
                pick_run = Some(r.id);
            }

            // Selected run → its alignment rows + the add-alignment form.
            if selected {
                ui.indent(("alns", r.id), |ui| {
                    for a in &alignments {
                        let asel = self.selected_alignment == Some(a.id);
                        let cov = if asel { coverage.as_ref() } else { None };
                        let (cov_s, call_s) = match cov {
                            Some(c) => (format!("{:.1}", c.mean_coverage), c.callable_bases.to_string()),
                            None => ("–".to_string(), "–".to_string()),
                        };
                        let row = egui::Frame::group(ui.style())
                            .fill(if asel { ACCENT.gamma_multiply(0.14) } else { ui.visuals().widgets.noninteractive.bg_fill })
                            .rounding(egui::Rounding::same(6.0))
                            .inner_margin(egui::Margin::symmetric(10.0, 8.0))
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.label(egui::RichText::new(&a.reference_build).color(ACCENT).strong());
                                    ui.label(egui::RichText::new(if a.bam_path.is_some() { a.aligner.as_str() } else { "Unknown" }).weak());
                                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                        if ui.small_button("🗑").on_hover_text("Delete alignment").clicked() {
                                            want_delete = Some(DataDelete::Alignment {
                                                id: a.id,
                                                run_id: r.id,
                                                label: format!("alignment {} ({})", a.id, a.reference_build),
                                            });
                                        }
                                        ui.add_space(10.0);
                                        ui.label(egui::RichText::new(format!("Callable: {call_s}")).weak().small());
                                        ui.add_space(10.0);
                                        ui.label(egui::RichText::new(format!("Coverage: {cov_s}")).weak().small());
                                    });
                                });
                            })
                            .response;
                        if row.interact(egui::Sense::click()).clicked() {
                            pick_aln = Some(a.id);
                        }
                    }
                    ui.add_space(4.0);
                    self.add_alignment_form(ui, r.id);
                });
            }
            ui.add_space(6.0);
        }

        if let Some(id) = pick_run {
            self.select_run(id);
        }
        if let Some(id) = pick_aln {
            self.select_alignment(id);
        }
        if want_delete.is_some() {
            self.confirm_data_delete = want_delete;
        }
        self.add_test_form(ui, guid);
    }

    /// The "Add test" (sequencing run) form.
    fn add_test_form(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.collapsing(self.tr("run.addTest"), |ui| {
            ui.horizontal(|ui| {
                ui.label(self.tr("form.testType"));
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
            if ui.add_enabled(ready, egui::Button::new(self.tr("run.addTest"))).clicked() {
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
            }
        });
    }

    /// The "Add alignment" form for a run. Picking a BAM/CRAM probes its header to auto-fill the
    /// reference build + aligner; the reference FASTA is never asked for (resolved from the build).
    fn add_alignment_form(&mut self, ui: &mut egui::Ui, run_id: i64) {
        ui.collapsing(self.tr("aln.add"), |ui| {
            ui.horizontal(|ui| {
                if ui.button(self.tr("common.pickBamCram")).clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("alignment", &["bam", "cram"]).pick_file() {
                        self.forms.aln_bam = p.to_string_lossy().into_owned();
                        // Probe the header to auto-fill build + aligner.
                        let _ = self.tx.send(Command::ProbeAlignment { path: p });
                        self.status = "Reading header…".into();
                    }
                }
                ui.label(if self.forms.aln_bam.is_empty() { "—" } else { self.forms.aln_bam.as_str() });
            });
            ui.add(
                egui::TextEdit::singleline(&mut self.forms.aln_reference_build)
                    .hint_text("reference build (auto-detected; editable)"),
            );
            ui.add(
                egui::TextEdit::singleline(&mut self.forms.aln_aligner).hint_text("aligner (auto-detected; editable)"),
            );
            ui.label(egui::RichText::new("Reference FASTA is resolved from the build automatically.").weak().small());
            let ready = !self.forms.aln_reference_build.trim().is_empty()
                && !self.forms.aln_aligner.trim().is_empty()
                && !self.forms.aln_bam.is_empty();
            if ui.add_enabled(ready, egui::Button::new(self.tr("aln.add"))).clicked() {
                let _ = self.tx.send(Command::AddAlignment(NewAlignment {
                    sequence_run_id: run_id,
                    reference_build: self.forms.aln_reference_build.trim().to_string(),
                    aligner: self.forms.aln_aligner.trim().to_string(),
                    variant_caller: None,
                    bam_path: opt(&self.forms.aln_bam),
                    reference_path: None, // resolved on demand from the build
                }));
                self.forms.aln_reference_build.clear();
                self.forms.aln_aligner.clear();
                self.forms.aln_bam.clear();
            }
        });
    }

    fn coverage_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_paths = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some() && a.reference_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_paths && !self.running, egui::Button::new(self.tr("btn.runCoverage"))).clicked() {
                self.running = true;
                self.status = format!("Running coverage on alignment #{alignment_id}…");
                let _ = self.tx.send(Command::RunCoverage(alignment_id));
            }
            if self.running {
                ui.spinner();
            }
            if !has_paths {
                ui.label(self.tr("hint.noBamRefPath"));
            }
        });

        match &self.coverage {
            None if !self.running => {
                ui.label(self.tr("coverage.none"));
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

    /// Inferred sex + read-level QC metrics for a single alignment.
    fn sex_metrics_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam && !self.running_sex, egui::Button::new(self.tr("btn.inferSex"))).clicked() {
                self.running_sex = true;
                self.status = "Inferring sex…".into();
                let _ = self.tx.send(Command::RunSex(alignment_id));
            }
            if ui.add_enabled(has_bam && !self.running_metrics, egui::Button::new(self.tr("btn.readMetrics"))).clicked() {
                self.running_metrics = true;
                self.status = "Collecting read metrics…".into();
                let _ = self.tx.send(Command::RunReadMetrics(alignment_id));
            }
            if ui.add_enabled(has_bam && !self.running_sv, egui::Button::new(self.tr("btn.callSv"))).clicked() {
                self.running_sv = true;
                self.status = "Calling structural variants (needs ≥10× coverage)…".into();
                let _ = self.tx.send(Command::RunSv(alignment_id));
            }
            if self.running_sex || self.running_metrics || self.running_sv {
                ui.spinner();
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        if let Some(s) = &self.sex {
            let sex = match s.inferred_sex {
                navigator_app::InferredSex::Male => "Male",
                navigator_app::InferredSex::Female => "Female",
                navigator_app::InferredSex::Unknown => "Unknown",
            };
            ui.label(format!(
                "Sex: {sex}  ·  chrX:autosome ratio {:.2}  ·  {:?} confidence",
                s.x_autosome_ratio, s.confidence
            ));
        }
        if let Some(m) = &self.read_metrics {
            egui::Grid::new("read_metrics_grid").striped(true).num_columns(2).show(ui, |ui| {
                let row = |ui: &mut egui::Ui, k: &str, v: String| {
                    ui.label(k);
                    ui.label(v);
                    ui.end_row();
                };
                row(ui, "Total reads", m.total_reads.to_string());
                row(ui, "% PF aligned", format!("{:.1}%", m.pct_pf_reads_aligned * 100.0));
                row(ui, "% proper pairs", format!("{:.1}%", m.pct_proper_pairs * 100.0));
                row(ui, "Mean read length", format!("{:.0}", m.mean_read_length));
                row(ui, "Median insert size", format!("{:.0}", m.median_insert_size));
                row(ui, "Pair orientation", m.pair_orientation.as_str().to_string());
                row(ui, "Mean MAPQ", format!("{:.1}", m.mean_mapping_quality));
            });
        }
        if let Some(sv) = &self.sv {
            ui.label(format!(
                "Structural variants: {} calls ({} CNV segments, {} discordant pairs)",
                sv.sv_calls.len(),
                sv.cnv_segments,
                sv.total_discordant_pairs
            ));
            for c in sv.sv_calls.iter().take(8) {
                ui.label(
                    egui::RichText::new(format!(
                        "  {} {}:{}-{} {}bp q{:.0}",
                        c.sv_type.as_str(), c.chrom, c.start, c.end, c.sv_len, c.quality
                    ))
                    .small()
                    .weak(),
                );
            }
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
                ui.label(self.tr("hint.signInToPublish"));
            }
            if self.publishing {
                ui.spinner();
            }
        });
    }

    /// Donor-level haplogroup consensus across all recorded sources (runs, Sanger, …).
    fn consensus_section(&mut self, ui: &mut egui::Ui) {
        if self.consensus_y.is_none() && self.consensus_mt.is_none() {
            ui.label(egui::RichText::new("No haplogroup consensus yet.").weak());
            return;
        }
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
        // Bind before the &mut self.forms borrow below (used inside the closure).
        let override_lbl = self.tr("form.override");
        let set_lbl = self.tr("common.set");
        let clear_lbl = self.tr("common.clear");
        let (hg_field, reason_field) = match dna_type {
            DnaType::Y => (&mut self.forms.override_y_haplogroup, &mut self.forms.override_y_reason),
            DnaType::Mt => (&mut self.forms.override_mt_haplogroup, &mut self.forms.override_mt_reason),
        };
        ui.horizontal(|ui| {
            ui.label(override_lbl);
            ui.add(egui::TextEdit::singleline(hg_field).hint_text("haplogroup").desired_width(140.0));
            ui.add(egui::TextEdit::singleline(reason_field).hint_text("reason").desired_width(180.0));
            let hg = hg_field.trim().to_string();
            let reason = reason_field.trim().to_string();
            if ui.add_enabled(!hg.is_empty(), egui::Button::new(set_lbl)).clicked() {
                self.status = format!("Overriding {label} consensus → {hg}");
                let _ = self.tx.send(Command::SetHaploOverride {
                    biosample_guid: guid,
                    dna_type,
                    haplogroup: hg,
                    reason: (!reason.is_empty()).then_some(reason),
                });
            }
            if ui.add_enabled(c.overridden, egui::Button::new(clear_lbl)).clicked() {
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
                ui.label(self.tr("hint.signInToPublish"));
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

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new(self.tr("btn.scanHeteroplasmy"))).clicked() {
                self.status = "Scanning chrM pileup for heteroplasmy…".into();
                let _ = self.tx.send(Command::LoadHeteroplasmy { alignment_id });
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        if let Some((id, sites)) = &self.heteroplasmy {
            if *id == alignment_id {
                if sites.is_empty() {
                    ui.label(self.tr("mt.noHeteroplasmy"));
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
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        let busy = self.estimating_ancestry || self.painting_running;
        ui.horizontal(|ui| {
            if ui
                .add_enabled(has_bam && !busy, egui::Button::new(self.tr("btn.estimateAncestry")))
                .clicked()
            {
                self.estimating_ancestry = true;
                self.status = "Estimating ancestry (genotyping AIMs panel)…".into();
                let _ = self.tx.send(Command::EstimateAncestry { alignment_id });
            }
            if ui
                .add_enabled(has_bam && !busy, egui::Button::new(self.tr("ancestry.paint")))
                .clicked()
            {
                self.painting_running = true;
                self.status = "Painting local ancestry (genotyping AIMs panel)…".into();
                let _ = self.tx.send(Command::PaintAncestry { alignment_id });
            }
            if busy {
                ui.spinner();
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        // Live genotyping progress (the slow per-contig pass over the BAM).
        if busy {
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
                // Donut of super-population proportions, beside the headline + composition bar.
                ui.horizontal(|ui| {
                    draw_ancestry_donut(ui, &result.super_population_summary);
                    ui.add_space(8.0);
                    ui.vertical(|ui| {
                        if let Some(top) = result.super_population_summary.first() {
                            ui.heading(format!("{} {:.1}%", top.super_population, top.percentage));
                        }
                        ui.label(format!(
                            "{}/{} panel SNPs · confidence {:.0}%",
                            result.snps_with_genotype,
                            result.snps_analyzed,
                            result.confidence_level * 100.0
                        ));
                        // Donor provenance: which method + reference build produced this estimate.
                        ui.label(
                            egui::RichText::new(format!("{} · {}", result.method, result.reference_version))
                                .small()
                                .weak(),
                        );
                        ui.add_space(4.0);
                        draw_composition_bar(ui, &result.super_population_summary);
                    });
                });
                ui.add_space(8.0);

                // Indented hierarchy: each super-population, then its fine populations (if the
                // panel is fine-grained). Proportions sum to 100% (supervised admixture).
                egui::Grid::new("ancestry_hierarchy").num_columns(2).spacing([24.0, 4.0]).show(ui, |ui| {
                    for s in &result.super_population_summary {
                        if s.percentage < 0.5 {
                            continue;
                        }
                        let super_code = s.populations.first().and_then(|c| population_super(c)).unwrap_or("");
                        ui.horizontal(|ui| {
                            let (r, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
                            ui.painter().circle_filled(r.center(), 5.0, parse_hex_color(&population_color(super_code)));
                            ui.strong(&s.super_population);
                        });
                        ui.strong(format!("{:.1}%", s.percentage));
                        ui.end_row();

                        // Fine sub-populations under this super (skip the self-row of a super panel).
                        let mut fine: Vec<&navigator_app::PopulationComponent> = result
                            .components
                            .iter()
                            .filter(|c| population_super(&c.population_code) == Some(super_code))
                            .filter(|c| c.population_code != super_code && c.percentage >= 0.5)
                            .collect();
                        fine.sort_by(|a, b| b.percentage.partial_cmp(&a.percentage).unwrap_or(std::cmp::Ordering::Equal));
                        for c in fine {
                            ui.horizontal(|ui| {
                                ui.add_space(20.0);
                                let (r, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                                ui.painter().circle_filled(r.center(), 3.5, parse_hex_color(&population_color(&c.population_code)));
                                ui.label(&c.population_name);
                            });
                            ui.label(format!("{:.1}%", c.percentage));
                            ui.end_row();
                        }
                    }
                });

                // Geographic distribution: each contributing population at its homeland,
                // sized by proportion and colored by continent.
                ui.add_space(8.0);
                ui.label(self.tr("ancestry.geo"));
                draw_ancestry_map(ui, &result.components);

                // PCA scatter: reference population centroids (PC1×PC2) + this sample's point.
                if let Some(coords) = &result.pca_coordinates {
                    let refs: &[(String, f64, f64)] = match &self.pca_reference {
                        Some((rid, pts)) if *rid == alignment_id => pts,
                        _ => &[],
                    };
                    if coords.len() >= 2 && !refs.is_empty() {
                        ui.add_space(8.0);
                        ui.label(self.tr("ancestry.pca"));
                        draw_pca_scatter(ui, (coords[0], coords[1]), refs);
                    }
                }
            }
            Some((id, None)) if *id == alignment_id && !busy => {
                ui.label(self.tr("ancestry.none"));
            }
            _ if !busy => {
                ui.label(self.tr("ancestry.none"));
            }
            _ => {}
        }

        // Local-ancestry painting (independent of the global estimate).
        if let Some((id, segments)) = &self.painting {
            if *id == alignment_id && !segments.is_empty() {
                ui.add_space(8.0);
                ui.label(self.tr("ancestry.local"));
                draw_chromosome_painting(ui, segments);
            }
        }
    }

    /// Y-haplogroup assignment for an alignment (calls chrY tree positions; FTDNA tree).
    fn y_haplogroup_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new(self.tr("btn.assignY"))).clicked() {
                self.status = "Assigning Y haplogroup (fetching FTDNA tree)…".into();
                let _ = self.tx.send(Command::AssignYHaplogroup { alignment_id });
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        // The persisted donor consensus (reloaded on select) — so the haplogroup stays visible
        // without re-running. The fresh per-run assignment, when present, adds the SNP detail.
        if let Some(c) = &self.consensus_y {
            ui.label(
                egui::RichText::new(format!(
                    "Haplogroup: {}  ({} run(s), confidence {:.2})",
                    c.haplogroup, c.run_count, c.confidence
                ))
                .strong(),
            );
            if !c.lineage.is_empty() {
                ui.label(egui::RichText::new(c.lineage.join(" › ")).small().weak());
            }
        }
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
            ui.label(self.tr("form.callableMask"));
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
                if ui.button(self.tr("form.chooseBed")).clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("BED", &["bed"]).pick_file() {
                        self.y_mask_path = Some(p);
                    }
                }
            });
        }
        ui.horizontal(|ui| {
            if ui.add_enabled(has_ref && !self.finding_private_y, egui::Button::new(self.tr("btn.findPrivateY"))).clicked() {
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
                ui.label(self.tr("hint.needsBamRef"));
            }
        });
        if let Some((id, bucket)) = &self.private_y {
            if *id == alignment_id {
                ui.label(format!("{} novel + {} off-path, below {}", bucket.novel(), bucket.off_path(), bucket.terminal));
                egui::CollapsingHeader::new("Private variants").id_salt(("privy", alignment_id)).show(ui, |ui| {
                    egui::Grid::new(("privy_grid", alignment_id)).striped(true).num_columns(4).show(ui, |ui| {
                        for h in ["table.position", "table.change", "table.depth", "table.class"] {
                            ui.strong(self.tr(h));
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
        let Some(panel_id) = self.selected_panel else {
            ui.label(self.tr("panel.selectInSidebar"));
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
            ui.label(self.tr("form.ploidy"));
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
                ui.label(self.tr("panel.notGenotyped"));
            }
            None => {}
        }

        // IBD compare against another genotyped alignment.
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            ui.label(self.tr("ibd.vs"));
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
            if ui.add_enabled(ready, egui::Button::new(self.tr("action.compare"))).clicked() {
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
            if ui.add_enabled(self.ibd_other.is_some(), egui::Button::new(self.tr("ibd.verify"))).clicked() {
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
                ui.label(self.tr("ibd.identity"));
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
                    ui.strong(self.tr("table.chr"));
                    ui.strong(self.tr("table.start"));
                    ui.strong(self.tr("table.end"));
                    ui.strong(self.tr("table.cm"));
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

    /// mtDNA haplogroup assigned directly from the alignment's chrM — the standalone counterpart
    /// to the Y-DNA section's "Assign Y haplogroup".
    fn mt_haplogroup_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);
        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new(self.tr("btn.assignMt"))).clicked() {
                self.status = "Assigning mtDNA haplogroup (fetching FTDNA mt tree)…".into();
                let _ = self.tx.send(Command::AssignMtdnaHaplogroupFromAlignment { alignment_id });
            }
            if !has_bam {
                ui.label(egui::RichText::new("(no BAM/CRAM path recorded)").weak());
            }
        });
        if let Some((id, assignment)) = &self.mt_haplogroup {
            if *id == alignment_id {
                show_assignment(ui, assignment);
            }
        }
    }

    /// De-novo haploid SNP calls for a specific `contig` (chrY on the Y-DNA tab, chrM on mtDNA).
    fn denovo_section(&mut self, ui: &mut egui::Ui, alignment_id: i64, contig: &str) {
        // Reference is resolved from the build on demand, so only the BAM is required.
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            let ready = has_bam && !self.running_denovo;
            let label = format!("{} ({contig})", self.tr("btn.runDenovo"));
            if ui.add_enabled(ready, egui::Button::new(label)).clicked() {
                self.running_denovo = true;
                self.denovo.remove(contig);
                self.status = format!("Calling {contig} on alignment #{alignment_id}…");
                let _ = self.tx.send(Command::RunDenovo { alignment_id, contig: contig.to_string() });
            }
            if self.running_denovo {
                ui.spinner();
            }
            if !has_bam {
                ui.label(egui::RichText::new("(no BAM/CRAM recorded)").weak());
            }
        });

        match self.denovo.get(contig) {
            None if !self.running_denovo => {
                ui.label(egui::RichText::new("No calls yet — run for this contig.").weak());
            }
            None => {}
            Some(calls) if calls.is_empty() => {
                ui.label(self.tr("denovo.noCalls"));
            }
            Some(calls) => {
                ui.label(format!("{} SNP call(s)", calls.len()));
                egui::Grid::new(("denovo_calls", contig)).striped(true).num_columns(4).show(ui, |ui| {
                    ui.strong(self.tr("table.position"));
                    ui.strong(self.tr("table.change"));
                    ui.strong(self.tr("table.depth"));
                    ui.strong(self.tr("table.af"));
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

        if self.denovo.get(contig).map(|c| !c.is_empty()).unwrap_or(false) {
            self.publish_row(ui, "Publish variants to PDS", Command::PublishVariants { alignment_id, contig: contig.to_string() });
        }
    }
}
