//! Plain-language **subject brief** model + the reference-content pack that supplies its narrative.
//!
//! This module is pure (no I/O): it owns the render-ready [`SubjectBrief`] tree, the [`BriefPack`]
//! reference-content schema, and the deterministic templating that turns structured analysis signals
//! (ages, depths, confidences) into casual-reader sentences. Composition — pulling the signals and
//! loading/enriching the pack — lives in `navigator-app::brief`; rendering lives in `navigator-ui`.
//!
//! The narrative content (haplogroup origins, ages, stories, test descriptions) is *not* derivable
//! from the analysis; it comes from the [`BriefPack`], shipped as a bundled seed and refreshed from a
//! CDN asset. Lookups fall back up the lineage path so a compact pack still tells a useful story for
//! a rare terminal haplogroup (see [`BriefPack::lineage_lookup`]).

use crate::testtype::TargetType;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ---------------------------------------------------------------------------------------------
// Reference pack (narrative content)
// ---------------------------------------------------------------------------------------------

/// One haplogroup's narrative content: when it formed, where it's associated with, and a short
/// curated story. Every field is optional so a sparse pack still contributes what it has.
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct HaploEntry {
    /// Years before present the haplogroup is estimated to have formed.
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
    /// Honest limitation, when there is one ("covers only the Y chromosome — no ancestry or
    /// maternal line").
    #[serde(default)]
    pub limits: Option<String>,
}

/// The bundled/downloaded reference pack. Maps are keyed by haplogroup name / test-type code.
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct BriefPack {
    pub version: String,
    #[serde(default)]
    pub y_haplogroups: HashMap<String, HaploEntry>,
    #[serde(default)]
    pub mt_haplogroups: HashMap<String, HaploEntry>,
    #[serde(default)]
    pub test_types: HashMap<String, TestEntry>,
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

    /// Look up `terminal` in `map`; if absent, walk the **root→tip** `lineage` and return the entry
    /// for the haplogroup *closest to the tip* that the pack covers. Returns the matched name (which
    /// may be an ancestor of `terminal`) and its entry, or `None` if nothing on the lineage is known.
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
    /// Served from the on-disk cache (a prior download).
    Cached,
    /// The bundled seed only (offline / CDN unavailable).
    Bundled,
    /// No pack at all (even the seed failed to parse) — briefs degrade to structured facts.
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

/// The "your test & quality" section — always present.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TestBrief {
    pub test_name: String,
    pub what_it_tells: String,
    pub limitations: Option<String>,
    pub quality_phrase: String,
    /// Drives a ✓ / ⚠ chip.
    pub quality_ok: bool,
}

/// A casual-reader brief for one subject. Sections are `Option` — each degrades to absent when its
/// data is missing (Y-only test → no maternal line; no haplogroup placed yet → no lineage section).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubjectBrief {
    pub headline: Headline,
    pub paternal: Option<LineageBrief>,
    pub maternal: Option<LineageBrief>,
    pub test: TestBrief,
    /// Global uncertainty notes.
    pub caveats: Vec<String>,
    /// Loaded pack version (for display), if any.
    pub pack_version: Option<String>,
    pub pack_status: PackStatus,
}

// ---------------------------------------------------------------------------------------------
// Templating (deterministic, unit-tested)
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

/// Round an age to a friendly magnitude so a precise estimate doesn't read as false precision
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

/// "formed roughly 4,200 years ago" — `None` when the age is unknown.
pub fn age_phrase(formed_ybp: Option<i32>) -> Option<String> {
    let ybp = formed_ybp?;
    if ybp <= 0 {
        return None;
    }
    Some(format!("formed roughly {} years ago", group_thousands(round_age(ybp))))
}

/// "associated with the Pontic-Caspian steppe and early Europe" — `None` when unknown.
pub fn origin_phrase(origin: Option<&str>) -> Option<String> {
    let o = origin?.trim();
    if o.is_empty() {
        return None;
    }
    Some(format!("associated with {o}"))
}

/// Plain-language confidence for a haplogroup placement, from the consensus confidence, the number
/// of sources that agree, and whether the sources conflict. Deliberately blunt about weak placements.
pub fn confidence_phrase(confidence: f64, run_count: usize, conflict: bool) -> String {
    if conflict {
        return "tentative — your sources don't fully agree on this branch".to_string();
    }
    let sources = match run_count {
        0 | 1 => "from a single test",
        _ => "confirmed across multiple tests",
    };
    if confidence >= 0.9 {
        format!("strong placement, {sources}")
    } else if confidence >= 0.6 {
        format!("good placement, {sources}")
    } else {
        format!("tentative placement, {sources}")
    }
}

/// Sequencing-depth quality, gated by what the test targets. Returns the phrase and an ok flag
/// (drives a ✓/⚠ chip). A targeted test (Y/mt) is judged on its own target depth, which is much
/// higher than a WGS average, so the WGS thresholds don't apply.
pub fn quality_phrase(mean_coverage: f64, target: TargetType) -> (String, bool) {
    let depth = format!("{mean_coverage:.0}× average depth");
    let (label, ok) = match target {
        // Whole-genome / autosomal / exome: judged on genome-wide average depth.
        TargetType::WholeGenome | TargetType::Autosomal | TargetType::XChromosome | TargetType::Mixed => {
            if mean_coverage >= 25.0 {
                ("high-quality", true)
            } else if mean_coverage >= 10.0 {
                ("good-quality", true)
            } else if mean_coverage >= 4.0 {
                ("usable but shallow", false)
            } else {
                ("very shallow — results are preliminary", false)
            }
        }
        // Targeted Y / mt: high on-target depth is the norm; be lenient.
        TargetType::YChromosome | TargetType::MtDna => {
            if mean_coverage >= 10.0 {
                ("high-quality", true)
            } else if mean_coverage >= 3.0 {
                ("good-quality", true)
            } else {
                ("shallow — results are preliminary", false)
            }
        }
    };
    (format!("{label} ({depth})"), ok)
}

/// Quality phrasing for a genotyping array (chip) test, which has no sequencing depth — judged on
/// the number of markers genotyped.
pub fn chip_quality_phrase(markers: usize) -> (String, bool) {
    if markers >= 100_000 {
        (
            format!("genotyping array ({} markers)", group_thousands(markers as i64)),
            true,
        )
    } else if markers > 0 {
        (
            format!("sparse genotyping array ({} markers)", group_thousands(markers as i64)),
            false,
        )
    } else {
        ("genotyping array".to_string(), true)
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
        assert_eq!(age_phrase(Some(4237)).unwrap(), "formed roughly 4,200 years ago");
        assert_eq!(age_phrase(Some(63500)).unwrap(), "formed roughly 64,000 years ago");
        assert_eq!(age_phrase(Some(842)).unwrap(), "formed roughly 850 years ago");
        assert_eq!(age_phrase(None), None);
        assert_eq!(age_phrase(Some(0)), None);
    }

    #[test]
    fn origin_phrasing() {
        assert_eq!(origin_phrase(Some("the steppe")).unwrap(), "associated with the steppe");
        assert_eq!(origin_phrase(None), None);
        assert_eq!(origin_phrase(Some("  ")), None);
    }

    #[test]
    fn confidence_phrasing_is_honest() {
        assert!(confidence_phrase(0.95, 3, false).starts_with("strong placement"));
        assert!(confidence_phrase(0.95, 3, false).contains("multiple tests"));
        assert!(confidence_phrase(0.95, 1, false).contains("single test"));
        assert!(confidence_phrase(0.7, 2, false).starts_with("good placement"));
        assert!(confidence_phrase(0.3, 2, false).starts_with("tentative placement"));
        assert!(confidence_phrase(0.99, 5, true).starts_with("tentative"));
    }

    #[test]
    fn quality_thresholds_depend_on_target() {
        let (p, ok) = quality_phrase(31.0, TargetType::WholeGenome);
        assert!(ok && p.starts_with("high-quality") && p.contains("31×"));
        let (_, ok) = quality_phrase(6.0, TargetType::WholeGenome);
        assert!(!ok, "6× WGS is shallow");
        // The same 6× on a targeted Y test is fine.
        let (_, ok) = quality_phrase(6.0, TargetType::YChromosome);
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
        // Nothing on the lineage is covered.
        assert!(pack.y_lookup("Q-M3", &["Q".into(), "Q-M242".into()]).is_none());
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
