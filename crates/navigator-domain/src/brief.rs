//! The plain-language **subject brief** model, and the reference-content pack that gives it its
//! narrative.
//!
//! This module is pure, with no I/O. It owns three things: the [`SubjectBrief`] tree, which is
//! ready to draw, the [`BriefPack`] reference-content schema, and the deterministic templates.
//! Those templates turn structured analysis signals (ages, depths, confidences) into sentences for
//! a casual reader. Composition lives in `navigator-app::brief`: it collects the signals, and it
//! loads the pack and adds to it. The UI draws the result, in `navigator-ui`.
//!
//! The narrative content (haplogroup origins, ages, stories, test descriptions) does *not* come
//! from the analysis. It comes from the [`BriefPack`], which ships as a bundled seed, and which a
//! CDN asset refreshes. A lookup falls back up the lineage path, so a compact pack still tells a
//! useful story for a rare terminal haplogroup (see [`BriefPack::lineage_lookup`]).

use crate::ancestry::SuperPopulationSummary;
use crate::i18n::{tr, tr_fmt, Lang};
use crate::roh::RohPattern;
use crate::testtype::TargetType;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ---------------------------------------------------------------------------------------------
// Reference pack (narrative content)
// ---------------------------------------------------------------------------------------------

/// The narrative content of one haplogroup: when it formed, which places it goes with, and a short
/// curated story. Every field is optional, so a sparse pack still gives what it has.
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct HaploEntry {
    /// Years before present when the haplogroup formed, as an estimate.
    #[serde(default)]
    pub formed_ybp: Option<i32>,
    /// Broad geographic / cultural association ("the Pontic-Caspian steppe and early Europe").
    #[serde(default)]
    pub origin: Option<String>,
    /// A 1–4 sentence plain-language narrative.
    #[serde(default)]
    pub story: Option<String>,
    /// Attribution for the narrative content.
    #[serde(default)]
    pub sources: Vec<String>,
}

/// One test type's plain-language description.
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct TestEntry {
    /// What the test tells you ("reads your whole genome, so it covers every lineage and ancestry").
    pub what: String,
    /// Honest limitation, when there is one ("covers only the Y chromosome: no ancestry, and no
    /// maternal line").
    #[serde(default)]
    pub limits: Option<String>,
}

/// One population's plain-language content, keyed by the (super-)population code or name the
/// ancestry estimate reports (e.g. "EUR", "European").
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct PopEntry {
    /// Friendly display name, when the estimate reports a bare code.
    #[serde(default)]
    pub name: Option<String>,
    /// A short plain-language note about the population.
    #[serde(default)]
    pub blurb: Option<String>,
}

/// The reference pack, from the bundle or from a download. The key of a map is a haplogroup name,
/// a test-type code, or a population code.
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct BriefPack {
    pub version: String,
    #[serde(default)]
    pub y_haplogroups: HashMap<String, HaploEntry>,
    #[serde(default)]
    pub mt_haplogroups: HashMap<String, HaploEntry>,
    #[serde(default)]
    pub test_types: HashMap<String, TestEntry>,
    #[serde(default)]
    pub populations: HashMap<String, PopEntry>,
}

impl BriefPack {
    /// Overlay `other` onto `self`, so a downloaded/cached pack augments (and overrides) the bundled
    /// seed entry-by-entry. `other`'s version wins when non-empty.
    pub fn merge(&mut self, other: BriefPack) {
        if !other.version.trim().is_empty() {
            self.version = other.version;
        }
        self.y_haplogroups.extend(other.y_haplogroups);
        self.mt_haplogroups.extend(other.mt_haplogroups);
        self.test_types.extend(other.test_types);
        self.populations.extend(other.populations);
    }

    /// Population content by code/name.
    pub fn population(&self, code: &str) -> Option<&PopEntry> {
        self.populations.get(code)
    }

    /// Y lookup with ancestor fallback (see [`Self::lineage_lookup`]).
    pub fn y_lookup(&self, terminal: &str, lineage: &[String]) -> Option<(String, &HaploEntry)> {
        Self::lineage_lookup(&self.y_haplogroups, terminal, lineage)
    }

    /// mtDNA lookup with ancestor fallback.
    pub fn mt_lookup(&self, terminal: &str, lineage: &[String]) -> Option<(String, &HaploEntry)> {
        Self::lineage_lookup(&self.mt_haplogroups, terminal, lineage)
    }

    /// Test-type description by code.
    pub fn test(&self, code: &str) -> Option<&TestEntry> {
        self.test_types.get(code)
    }

    /// Look up `terminal` in `map`. If it is not there, walk the **root→tip** `lineage` and return
    /// the entry for the haplogroup *closest to the tip* that the pack covers. Returns the name it
    /// matched, which can be an ancestor of `terminal`, and its entry. `None` if the pack covers
    /// nothing on the lineage.
    fn lineage_lookup<'a>(
        map: &'a HashMap<String, HaploEntry>,
        terminal: &str,
        lineage: &[String],
    ) -> Option<(String, &'a HaploEntry)> {
        if let Some(e) = map.get(terminal) {
            return Some((terminal.to_string(), e));
        }
        // lineage is root→tip, so the last match is the deepest covered ancestor.
        let mut found: Option<(String, &HaploEntry)> = None;
        for name in lineage {
            if let Some(e) = map.get(name) {
                found = Some((name.clone(), e));
            }
        }
        found
    }
}

// ---------------------------------------------------------------------------------------------
// Brief model (render-ready)
// ---------------------------------------------------------------------------------------------

/// Which parental line a [`LineageBrief`] describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LineageKind {
    Paternal,
    Maternal,
}

/// Provenance of the loaded reference pack, surfaced so the UI can show how fresh the narrative is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PackStatus {
    /// Refreshed from the CDN this session.
    Downloaded,
    /// Served from the on-disk cache (an earlier download).
    Cached,
    /// The bundled seed only (offline / CDN unavailable).
    Bundled,
    /// No pack at all (even the seed failed to parse). A brief then degrades to structured facts.
    Unavailable,
}

impl PackStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            PackStatus::Downloaded => "downloaded",
            PackStatus::Cached => "cached",
            PackStatus::Bundled => "bundled",
            PackStatus::Unavailable => "unavailable",
        }
    }
}

/// The top-of-brief summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Headline {
    pub name: String,
    /// Friendly test label, e.g. "Whole-genome sequence".
    pub test_chip: String,
    /// One-sentence "who you are" line.
    pub summary: String,
}

/// A paternal or maternal lineage section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LineageBrief {
    pub kind: LineageKind,
    /// Terminal haplogroup, e.g. "R-FGC29071".
    pub haplogroup: String,
    /// Root→tip lineage path (for an optional expandable trail).
    pub lineage_path: Vec<String>,
    /// When the narrative is for an *ancestor* of the terminal (compact-pack fallback), the matched
    /// ancestor's name; `None` when the story is for the terminal itself.
    pub matched_ancestor: Option<String>,
    pub age_phrase: Option<String>,
    pub origin_phrase: Option<String>,
    pub story: Option<String>,
    pub confidence_phrase: String,
    pub sources: Vec<String>,
}

/// One ancient-ancestry component (a prehistoric source population, e.g. "Steppe pastoralist"),
/// with its share, display color, and an optional plain-language explanation from the pack.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AncientComponent {
    pub code: String,
    pub name: String,
    /// 0.0–100.0
    pub percentage: f64,
    /// Display color (hex), so the UI pie matches the rest of the ancestry palette.
    pub color: String,
    /// Plain-language explanation of this ancient source (from the reference pack).
    pub blurb: Option<String>,
}

/// The ancestry-composition section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AncestryBrief {
    /// A one-line summary, for example "Predominantly European".
    pub summary_phrase: String,
    /// Continental breakdown (carried whole so the UI can reuse the existing donut).
    pub super_populations: Vec<SuperPopulationSummary>,
    /// Fine-grained populations `(name, percentage)`, when a detailed estimate exists.
    pub fine_pops: Vec<(String, f64)>,
    /// Ancient-ancestry components (prehistoric source populations), when that estimate exists.
    pub ancient_pops: Vec<AncientComponent>,
    /// Optional plain-language note about the mix (from the reference pack).
    pub interpretation: Option<String>,
    /// What made the estimate, for example "estimated from 412,000 genome-wide markers".
    pub method_note: String,
}

/// The runs-of-homozygosity section (relatedness and endogamy). It is here only after ROH runs for
/// the subject. F_ROH is the share of the genome in long homozygous runs. That share shows how much
/// recent shared ancestry the two parental lines of a person have.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RohBrief {
    /// Inbreeding coefficient F_ROH (0.0–1.0), the physical share of autosomes in runs of homozygosity.
    pub f_roh: f64,
    /// Plain-language pattern label, e.g. "Outbred", "Endogamous background", "Recent shared ancestry".
    pub pattern: String,
    /// One-sentence casual explanation.
    pub summary_phrase: String,
    /// Number of homozygous runs reported.
    pub n_segments: usize,
    /// Total length of runs (Mb) and the longest single run (Mb).
    pub total_mb: f64,
    pub longest_mb: f64,
}

/// The "your test & quality" section. It is always there.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TestBrief {
    pub test_name: String,
    pub what_it_tells: String,
    pub limitations: Option<String>,
    pub quality_phrase: String,
    /// Drives a ✓ / ⚠ chip.
    pub quality_ok: bool,
}

/// An offer to rebuild one alignment against CHM13, for a reader who should not have to know what
/// a reference build is.
///
/// It is narrow on purpose. Realignment costs hours and hundreds of GB, and it gives exactly one
/// thing: Y-chromosome discovery on a reference whose Y is complete.
///
/// So the offer appears only when two things hold. There must be a paternal line to improve, and an
/// alignment on an older reference that no realignment has touched. Everybody else never sees it: a
/// subject with no Y, a chip-only subject, and a subject already on CHM13. For them the answer
/// would not change.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RealignOffer {
    /// The alignment that would be rebuilt.
    pub alignment_id: i64,
    /// The reference it is on today, named for the prompt ("GRCh38").
    pub current_build: String,
}

/// A brief for one subject, for a casual reader. Each section is an `Option`, and it becomes
/// absent when its data is missing. A Y-only test gives no maternal line, and a subject with no
/// haplogroup placed yet gives no lineage section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubjectBrief {
    pub headline: Headline,
    pub paternal: Option<LineageBrief>,
    pub maternal: Option<LineageBrief>,
    pub ancestry: Option<AncestryBrief>,
    /// Relatedness and endogamy, read from runs of homozygosity. Absent until ROH runs.
    #[serde(default)]
    pub roh: Option<RohBrief>,
    /// Archaic (Neanderthal) marker count. Absent until that count runs.
    #[serde(default)]
    pub archaic: Option<ArchaicBrief>,
    pub test: TestBrief,
    /// True when the subject has a sequencing alignment that no analysis has touched: the data is
    /// there, and no coverage exists. This is the signal for the one-click "Analyze" prompt in
    /// Simple mode. False for a subject that an analysis already covered, and for a subject with no
    /// alignment (chip or VCF only, with nothing to analyze).
    #[serde(default)]
    pub needs_analysis: bool,
    /// An alignment that is worth a realignment to CHM13, when that would tell the reader
    /// something new. Absent whenever it would not. See [`RealignOffer`].
    #[serde(default)]
    pub realign_offer: Option<RealignOffer>,
    /// Global uncertainty notes.
    pub caveats: Vec<String>,
    /// Loaded pack version (for display), if any.
    pub pack_version: Option<String>,
    pub pack_status: PackStatus,
    /// True when live AppView or DecodingUs content (haplogroup ages, provenance) went into this.
    pub enriched: bool,
}

// ---------------------------------------------------------------------------------------------
// Templates (deterministic, with unit tests)
// ---------------------------------------------------------------------------------------------

/// Group an integer with thousands separators ("4000" → "4,000"). Small helper to keep the phrase
/// builders dependency-free.
fn group_thousands(n: i64) -> String {
    let s = n.abs().to_string();
    let mut out = String::new();
    let bytes = s.as_bytes();
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(*b as char);
    }
    if n < 0 {
        format!("-{out}")
    } else {
        out
    }
}

/// Round an age to a friendly magnitude so a precise estimate does not read as false precision
/// (4,237 → "about 4,200"; 63,500 → "about 64,000"; 850 → "about 850").
fn round_age(ybp: i32) -> i64 {
    let y = ybp.max(0) as i64;
    let step = if y >= 50_000 {
        1_000
    } else if y >= 10_000 {
        500
    } else if y >= 1_000 {
        100
    } else {
        50
    };
    ((y + step / 2) / step) * step
}

/// "formed roughly 4,200 years ago". `None` when the age is unknown.
pub fn age_phrase(lang: Lang, formed_ybp: Option<i32>) -> Option<String> {
    let ybp = formed_ybp?;
    if ybp <= 0 {
        return None;
    }
    Some(tr_fmt(lang, "brief.agePhrase", &[&group_thousands(round_age(ybp))]))
}

/// "associated with the Pontic-Caspian steppe and early Europe". `None` when unknown.
pub fn origin_phrase(lang: Lang, origin: Option<&str>) -> Option<String> {
    let o = origin?.trim();
    if o.is_empty() {
        return None;
    }
    Some(tr_fmt(lang, "brief.originPhrase", &[o]))
}

/// Plain-language confidence for a haplogroup placement, from the consensus confidence, the number
/// of sources that agree, and whether the sources conflict. Deliberately blunt about weak placements.
pub fn confidence_phrase(lang: Lang, confidence: f64, run_count: usize, conflict: bool) -> String {
    if conflict {
        return tr(lang, "brief.confidenceConflict").to_string();
    }
    let sources = match run_count {
        0 | 1 => tr(lang, "brief.confidenceSingle"),
        _ => tr(lang, "brief.confidenceMultiple"),
    };
    let key = if confidence >= 0.9 {
        "brief.confidenceStrong"
    } else if confidence >= 0.6 {
        "brief.confidenceGood"
    } else {
        "brief.confidenceTentative"
    };
    tr_fmt(lang, key, &[sources])
}

/// Sequencing-depth quality, gated by what the test targets. Returns the phrase and an ok flag,
/// which drives a ✓/⚠ chip. A targeted test (Y or mt) uses its own target depth. That depth is much
/// higher than a WGS average, so the WGS thresholds do not apply.
pub fn quality_phrase(lang: Lang, mean_coverage: f64, target: TargetType) -> (String, bool) {
    let (label_key, ok) = match target {
        // Whole-genome / autosomal / exome: judged on genome-wide average depth.
        TargetType::WholeGenome | TargetType::Autosomal | TargetType::XChromosome | TargetType::Mixed => {
            if mean_coverage >= 25.0 {
                ("brief.qualityHigh", true)
            } else if mean_coverage >= 10.0 {
                ("brief.qualityGood", true)
            } else if mean_coverage >= 4.0 {
                ("brief.qualityShallow", false)
            } else {
                ("brief.qualityVeryShallow", false)
            }
        }
        // Targeted Y / mt: high on-target depth is the norm; be lenient.
        TargetType::YChromosome | TargetType::MtDna => {
            if mean_coverage >= 10.0 {
                ("brief.qualityHigh", true)
            } else if mean_coverage >= 3.0 {
                ("brief.qualityGood", true)
            } else {
                ("brief.qualityTargetedShallow", false)
            }
        }
    };
    let depth = tr_fmt(lang, "brief.depth", &[&format!("{mean_coverage:.0}")]);
    (
        tr_fmt(lang, "brief.qualityWithDepth", &[tr(lang, label_key), &depth]),
        ok,
    )
}

/// A one-line summary of an ancestry mix, from the continental breakdown. It sorts a copy by
/// share, so a caller does not have to sort first. Empty input → a neutral phrase.
pub fn ancestry_summary(lang: Lang, super_pops: &[SuperPopulationSummary]) -> String {
    let mut sorted: Vec<&SuperPopulationSummary> = super_pops.iter().collect();
    sorted.sort_by(|a, b| b.percentage.total_cmp(&a.percentage));
    match sorted.as_slice() {
        [] => tr(lang, "brief.ancestryNone").to_string(),
        [top, rest @ ..] => {
            let name = top.super_population.as_str();
            if top.percentage >= 85.0 {
                tr_fmt(lang, "brief.ancestryPredominantly", &[name])
            } else if top.percentage >= 55.0 {
                match rest.first().filter(|s| s.percentage >= 10.0) {
                    Some(second) => tr_fmt(lang, "brief.ancestryMostlyWith", &[name, &second.super_population]),
                    None => tr_fmt(lang, "brief.ancestryMostly", &[name]),
                }
            } else {
                match rest.first() {
                    Some(second) => tr_fmt(lang, "brief.ancestryMix", &[name, &second.super_population]),
                    None => tr_fmt(lang, "brief.ancestryMostly", &[name]),
                }
            }
        }
    }
}

/// The casual-mode read of the archaic (Neanderthal) marker count.
///
/// This says it the same way as the Advanced card, on purpose. It gives a **count over what the
/// test of the subject covered**, never a "percent Neanderthal", and no Denisovan figure. Outside
/// Oceania that signal sits at the noise floor, and to report it would invent a result (design
/// §7).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicBrief {
    /// Archaic-derived allele copies carried.
    pub total_copies: u32,
    /// Copies assayed: two times the number of marker sites the data of the subject covered.
    pub possible_copies: u32,
    /// Marker sites covered, and how many exist in the panel.
    pub called_sites: usize,
    pub panel_sites: usize,
    /// Percentile within `cohort`, when the comparison is sound.
    pub percentile: Option<f32>,
    pub cohort: Option<String>,
    /// Plain-language pattern label, e.g. "Typical amount".
    pub pattern: String,
    /// One-sentence casual explanation.
    pub summary_phrase: String,
}

/// Build the casual archaic read.
///
/// The pattern comes from the **percentile** when one is available. A raw count says nothing on
/// its own, because it scales with how many sites the test covered. With no percentile, the label
/// stays neutral ("Neanderthal markers found"), and does not suggest a comparison that nobody
/// made.
pub fn archaic_brief(
    lang: Lang,
    total_copies: u32,
    possible_copies: u32,
    called_sites: usize,
    panel_sites: usize,
    percentile: Option<f32>,
    cohort: Option<String>,
) -> ArchaicBrief {
    let (pattern_key, summary_key) = match percentile {
        Some(p) if p < 25.0 => ("brief.archaicFewer", "brief.archaicFewerSummary"),
        Some(p) if p > 75.0 => ("brief.archaicMore", "brief.archaicMoreSummary"),
        Some(_) => ("brief.archaicTypical", "brief.archaicTypicalSummary"),
        None => ("brief.archaicFound", "brief.archaicFoundSummary"),
    };
    ArchaicBrief {
        total_copies,
        possible_copies,
        called_sites,
        panel_sites,
        percentile,
        cohort,
        pattern: tr(lang, pattern_key).to_string(),
        summary_phrase: tr(lang, summary_key).to_string(),
    }
}

/// Put the plain-language words on the runs-of-homozygosity verdict that the analysis engine
/// already reached. `pattern` is the [`RohPattern`] of `navigator_analysis::roh` itself. This does
/// *not* derive the classification again. So the Simple brief and the Advanced ROH chart can never
/// disagree about whether a subject reads as outbred, endogamous, or recently consanguineous. The
/// words are strictly about *shared ancestry between the lines of the parents*. That is a
/// genealogical read, and never a clinical one.
pub fn roh_brief(
    lang: Lang,
    pattern: RohPattern,
    f_roh: f64,
    n_segments: usize,
    total_mb: f64,
    longest_mb: f64,
) -> RohBrief {
    // No runs at all reads as outbred regardless of the classifier's view of an empty distribution.
    let effective = if n_segments == 0 { RohPattern::Outbred } else { pattern };
    let (total, longest, runs) = (
        format!("{total_mb:.0}"),
        format!("{longest_mb:.0}"),
        n_segments.to_string(),
    );
    let (pattern_key, summary) = match effective {
        RohPattern::Outbred => ("brief.rohOutbred", tr(lang, "brief.rohOutbredSummary").to_string()),
        RohPattern::RecentConsanguinity => (
            "brief.rohRecent",
            tr_fmt(lang, "brief.rohRecentSummary", &[&total, &longest]),
        ),
        RohPattern::Endogamy => (
            "brief.rohEndogamy",
            tr_fmt(lang, "brief.rohEndogamySummary", &[&total, &runs]),
        ),
        RohPattern::Mixed => (
            "brief.rohMixed",
            tr_fmt(lang, "brief.rohMixedSummary", &[&total, &runs, &longest]),
        ),
    };
    RohBrief {
        f_roh,
        pattern: tr(lang, pattern_key).to_string(),
        summary_phrase: summary,
        n_segments,
        total_mb,
        longest_mb,
    }
}

/// "estimated from 412,000 genome-wide markers" / "estimated from 220 ancestry-informative markers".
pub fn ancestry_method_note(lang: Lang, snps_with_genotype: usize, panel_type: &str) -> String {
    let kind = if panel_type.eq_ignore_ascii_case("aims") {
        tr(lang, "brief.markersAims")
    } else {
        tr(lang, "brief.markersGenomeWide")
    };
    tr_fmt(
        lang,
        "brief.ancestryMethodNote",
        &[&group_thousands(snps_with_genotype as i64), kind],
    )
}

/// Quality words for a genotyping array (chip) test, which has no sequencing depth. The judgement
/// is on how many markers the test genotyped.
pub fn chip_quality_phrase(lang: Lang, markers: usize) -> (String, bool) {
    let count = group_thousands(markers as i64);
    if markers >= 100_000 {
        (tr_fmt(lang, "brief.chipArrayMarkers", &[&count]), true)
    } else if markers > 0 {
        (tr_fmt(lang, "brief.chipArraySparse", &[&count]), false)
    } else {
        (tr(lang, "brief.chipArray").to_string(), true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thousands_grouping() {
        assert_eq!(group_thousands(0), "0");
        assert_eq!(group_thousands(850), "850");
        assert_eq!(group_thousands(4200), "4,200");
        assert_eq!(group_thousands(64000), "64,000");
        assert_eq!(group_thousands(1234567), "1,234,567");
    }

    #[test]
    fn age_rounding_is_friendly() {
        assert_eq!(
            age_phrase(Lang::En, Some(4237)).unwrap(),
            "formed roughly 4,200 years ago"
        );
        assert_eq!(
            age_phrase(Lang::En, Some(63500)).unwrap(),
            "formed roughly 64,000 years ago"
        );
        assert_eq!(age_phrase(Lang::En, Some(842)).unwrap(), "formed roughly 850 years ago");
        assert_eq!(age_phrase(Lang::En, None), None);
        assert_eq!(age_phrase(Lang::En, Some(0)), None);
    }

    #[test]
    fn origin_phrasing() {
        assert_eq!(
            origin_phrase(Lang::En, Some("the steppe")).unwrap(),
            "associated with the steppe"
        );
        assert_eq!(origin_phrase(Lang::En, None), None);
        assert_eq!(origin_phrase(Lang::En, Some("  ")), None);
    }

    #[test]
    fn roh_brief_reads_the_pattern() {
        // The words follow the verdict of the analysis engine, and nothing derives them again
        // from the numbers. Trace-endogamy or outbred (James: F_ROH ~0.008) gives outbred words,
        // and cites no number that would alarm the reader.
        let outbred = roh_brief(Lang::En, RohPattern::Outbred, 0.008, 6, 22.7, 7.1);
        assert_eq!(outbred.pattern, "Outbred");
        assert!(outbred.summary_phrase.contains("outbred"));
        // Elevated F_ROH from many short runs → endogamous background.
        let endog = roh_brief(Lang::En, RohPattern::Endogamy, 0.035, 40, 90.0, 8.0);
        assert_eq!(endog.pattern, "Endogamous background");
        assert!(endog.summary_phrase.contains("shared ancestry"));
        // Long-run-dominated → recent shared ancestry (consanguinity).
        let recent = roh_brief(Lang::En, RohPattern::RecentConsanguinity, 0.08, 12, 220.0, 40.0);
        assert_eq!(recent.pattern, "Recent shared ancestry");
        assert!(recent.summary_phrase.contains("few generations"));
        // Both classes are here. This used to read as "recent" only because one run was ≥15 Mb.
        let mixed = roh_brief(Lang::En, RohPattern::Mixed, 0.05, 30, 150.0, 18.0);
        assert_eq!(mixed.pattern, "Mixed shared ancestry");
        // No runs at all always reads as outbred, whatever the classifier says of an empty set.
        assert_eq!(
            roh_brief(Lang::En, RohPattern::Mixed, 0.0, 0, 0.0, 0.0).pattern,
            "Outbred"
        );
    }

    #[test]
    fn confidence_phrasing_is_honest() {
        assert!(confidence_phrase(Lang::En, 0.95, 3, false).starts_with("strong placement"));
        assert!(confidence_phrase(Lang::En, 0.95, 3, false).contains("multiple tests"));
        assert!(confidence_phrase(Lang::En, 0.95, 1, false).contains("single test"));
        assert!(confidence_phrase(Lang::En, 0.7, 2, false).starts_with("good placement"));
        assert!(confidence_phrase(Lang::En, 0.3, 2, false).starts_with("tentative placement"));
        assert!(confidence_phrase(Lang::En, 0.99, 5, true).starts_with("tentative"));
    }

    #[test]
    fn quality_thresholds_depend_on_target() {
        let (p, ok) = quality_phrase(Lang::En, 31.0, TargetType::WholeGenome);
        assert!(ok && p.starts_with("high-quality") && p.contains("31×"));
        let (_, ok) = quality_phrase(Lang::En, 6.0, TargetType::WholeGenome);
        assert!(!ok, "6× WGS is shallow");
        // The same 6× on a targeted Y test is fine.
        let (_, ok) = quality_phrase(Lang::En, 6.0, TargetType::YChromosome);
        assert!(ok);
    }

    #[test]
    fn pack_lineage_fallback() {
        let mut y = HashMap::new();
        y.insert(
            "R-M269".to_string(),
            HaploEntry {
                formed_ybp: Some(6400),
                origin: Some("the steppe".into()),
                story: Some("…".into()),
                sources: vec!["YFull".into()],
            },
        );
        let pack = BriefPack {
            version: "test".into(),
            y_haplogroups: y,
            ..Default::default()
        };
        // Direct hit.
        assert_eq!(pack.y_lookup("R-M269", &[]).unwrap().0, "R-M269");
        // Terminal absent → fall back to the deepest covered ancestor on the root→tip lineage.
        let lineage = vec![
            "R".to_string(),
            "R-M269".to_string(),
            "R-P312".to_string(),
            "R-FGC29071".to_string(),
        ];
        let (matched, _) = pack.y_lookup("R-FGC29071", &lineage).unwrap();
        assert_eq!(matched, "R-M269");
        // The pack covers nothing on the lineage.
        assert!(pack.y_lookup("Q-M3", &["Q".into(), "Q-M242".into()]).is_none());
    }

    fn sp(name: &str, pct: f64) -> SuperPopulationSummary {
        SuperPopulationSummary {
            super_population: name.to_string(),
            percentage: pct,
            populations: vec![],
        }
    }

    #[test]
    fn ancestry_summary_framing() {
        assert_eq!(
            ancestry_summary(Lang::En, &[]),
            "Ancestry composition not yet estimated"
        );
        assert_eq!(
            ancestry_summary(Lang::En, &[sp("European", 92.0)]),
            "Predominantly European"
        );
        // The function sorts input by share.
        assert_eq!(
            ancestry_summary(Lang::En, &[sp("African", 30.0), sp("European", 70.0)]),
            "Mostly European, with some African"
        );
        assert_eq!(
            ancestry_summary(Lang::En, &[sp("European", 45.0), sp("East Asian", 40.0)]),
            "A mix of European and East Asian"
        );
        // Dominant but lone.
        assert_eq!(
            ancestry_summary(Lang::En, &[sp("European", 60.0), sp("African", 3.0)]),
            "Mostly European"
        );
    }

    #[test]
    fn ancestry_method_note_kind() {
        assert_eq!(
            ancestry_method_note(Lang::En, 412000, "genome-wide"),
            "estimated from 412,000 genome-wide markers"
        );
        assert_eq!(
            ancestry_method_note(Lang::En, 220, "aims"),
            "estimated from 220 ancestry-informative markers"
        );
    }

    /// Why this prose goes through the catalog: a reader of another language gets it in their own.
    /// Before, every one of these sentences was a hardcoded English literal below the UI layer. The
    /// Simple-mode brief was then permanently English, whatever locale the user chose.
    #[test]
    fn brief_prose_is_localized_not_hardcoded() {
        // Sentences differ by language, and the Spanish is real text, and not a key fallback.
        for (en, es) in [
            (
                age_phrase(Lang::En, Some(4237)).unwrap(),
                age_phrase(Lang::Es, Some(4237)).unwrap(),
            ),
            (
                confidence_phrase(Lang::En, 0.95, 3, false),
                confidence_phrase(Lang::Es, 0.95, 3, false),
            ),
            (
                ancestry_summary(Lang::En, &[sp("European", 92.0)]),
                ancestry_summary(Lang::Es, &[sp("European", 92.0)]),
            ),
            (
                quality_phrase(Lang::En, 31.0, TargetType::WholeGenome).0,
                quality_phrase(Lang::Es, 31.0, TargetType::WholeGenome).0,
            ),
            (
                roh_brief(Lang::En, RohPattern::Endogamy, 0.035, 40, 90.0, 8.0).summary_phrase,
                roh_brief(Lang::Es, RohPattern::Endogamy, 0.035, 40, 90.0, 8.0).summary_phrase,
            ),
        ] {
            assert_ne!(en, es, "not translated: {en}");
            assert!(!es.starts_with("brief."), "rendered a raw key instead of text: {es}");
        }

        // The numbers survive interpolation in both languages. A template that dropped `{0}`
        // would lose the figure the sentence is about, and give no message.
        assert!(age_phrase(Lang::Es, Some(4237)).unwrap().contains("4,200"));
        let es_roh = roh_brief(Lang::Es, RohPattern::Endogamy, 0.035, 40, 90.0, 8.0);
        assert!(es_roh.summary_phrase.contains("90") && es_roh.summary_phrase.contains("40"));
        // And no `{n}` placeholder stays without a substitution.
        for text in [
            age_phrase(Lang::Es, Some(4237)).unwrap(),
            es_roh.summary_phrase,
            quality_phrase(Lang::Es, 31.0, TargetType::WholeGenome).0,
        ] {
            assert!(!text.contains('{'), "unsubstituted placeholder in: {text}");
        }
    }

    #[test]
    fn pack_merge_overlays() {
        let mut base = BriefPack {
            version: "seed".into(),
            ..Default::default()
        };
        base.y_haplogroups.insert("A".into(), HaploEntry::default());
        let mut over = BriefPack {
            version: "2026.07".into(),
            ..Default::default()
        };
        over.y_haplogroups.insert(
            "B".into(),
            HaploEntry {
                formed_ybp: Some(100),
                ..Default::default()
            },
        );
        base.merge(over);
        assert_eq!(base.version, "2026.07");
        assert!(base.y_haplogroups.contains_key("A") && base.y_haplogroups.contains_key("B"));
    }
}
