//! Y-STR panel taxonomy and classification: a port of the Scala `str-panels.conf` and
//! `StrPanelService`.
//!
//! Static data. It holds the FTDNA tiers (Y-12 ⊂ Y-25 ⊂ Y-37 ⊂ Y-67 ⊂ Y-111) and the YSEQ panels
//! (Alpha, Beta, Delta, Gamma, and the exclusive and Kittler detection sets). Three jobs use it.
//! They classify a profile into its highest reached tier, group markers into tiers for the
//! FTDNA-style report, and detect the provider from exclusive markers.
//!
//! This does **not** list the FTDNA Big-Y bonus tiers (Y-500/Y-700, ~700 markers), and that is
//! deliberate. A marker outside the tiers here lands in the [`EXTENDED`] group in
//! [`assign_markers_to_panels`], and still counts toward the total. A full list of Y-500 and Y-700
//! is a follow-up.

use std::collections::{HashMap, HashSet};

use crate::strprofile::StrMarker;

/// One panel tier: the markers **new** to this tier (FTDNA panels are cumulative, so each lists only
/// its additions). `marketing_count` is the vendor-advertised size; `actual_count` is the distinct
/// marker-key count (smaller, because a multi-value marker like DYS464 counts as more than one
/// value).
pub struct StrPanelDef {
    pub id: &'static str,
    pub name: &'static str,
    pub provider: &'static str,
    pub marketing_count: u32,
    pub actual_count: u32,
    pub order: u32,
    pub markers: &'static [&'static str],
}

impl StrPanelDef {
    /// Marker count required to classify into this tier (the distinct-key count).
    pub fn threshold(&self) -> u32 {
        self.actual_count
    }
}

pub struct ProviderDef {
    pub key: &'static str,
    pub display_name: &'static str,
    pub cumulative: bool,
    pub exclusive_markers: &'static [&'static str],
}

/// Label for markers that fall outside the enumerated tiers (FTDNA Big-Y bonus STRs, etc.).
pub const EXTENDED: &str = "Extended (Y-500/700)";

// ---- FTDNA cumulative tiers (each lists only markers NEW to that tier) --------------------------

static FTDNA_Y12: &[&str] = &[
    "DYS393", "DYS390", "DYS19", "DYS391", "DYS385", "DYS426", "DYS388", "DYS439", "DYS389I", "DYS392", "DYS389II",
];
static FTDNA_Y25: &[&str] = &[
    "DYS458", "DYS459", "DYS455", "DYS454", "DYS447", "DYS437", "DYS448", "DYS449", "DYS464",
];
static FTDNA_Y37: &[&str] = &[
    "DYS460",
    "Y-GATA-H4",
    "YCAII",
    "DYS456",
    "DYS607",
    "DYS576",
    "DYS570",
    "CDY",
    "DYS442",
    "DYS438",
];
static FTDNA_Y67: &[&str] = &[
    "DYS531", "DYS578", "DYF395S1", "DYS590", "DYS537", "DYS641", "DYS472", "DYF406S1", "DYS511", "DYS425", "DYS413",
    "DYS557", "DYS594", "DYS436", "DYS490", "DYS534", "DYS450", "DYS444", "DYS481", "DYS520", "DYS446", "DYS617",
    "DYS568", "DYS487", "DYS572", "DYS640", "DYS492", "DYS565",
];
static FTDNA_Y111: &[&str] = &[
    "DYS710",
    "DYS485",
    "DYS632",
    "DYS495",
    "DYS540",
    "DYS714",
    "DYS716",
    "DYS717",
    "DYS505",
    "DYS556",
    "DYS549",
    "DYS589",
    "DYS522",
    "DYS494",
    "DYS533",
    "DYS636",
    "DYS575",
    "DYS638",
    "DYS462",
    "DYS452",
    "DYS445",
    "Y-GATA-A10",
    "DYS463",
    "DYS441",
    "Y-GGAAT-1B07",
    "DYS525",
    "DYS712",
    "DYS593",
    "DYS650",
    "DYS532",
    "DYS715",
    "DYS504",
    "DYS513",
    "DYS561",
    "DYS552",
    "DYS726",
    "DYS635",
    "DYS587",
    "DYS643",
    "DYS497",
    "DYS510",
    "DYS434",
    "DYS461",
    "DYS435",
];

// ---- YSEQ panels --------------------------------------------------------------------------------

static YSEQ_ALPHA: &[&str] = &[
    "DYS391",
    "DYS389I",
    "DYS437",
    "DYS439",
    "DYS389II",
    "DYS438",
    "DYS426",
    "DYS393",
    "YCAII",
    "DYS390",
    "DYS385",
    "Y-GATA-H4",
    "DYS388",
    "DYS447",
    "DYS19",
    "DYS392",
];
static YSEQ_BETA: &[&str] = &[
    "DYS458", "DYS455", "DYS454", "DYS464", "DYS448", "DYS449", "DYS456", "DYS576", "CDY", "DYS460", "DYS459",
    "DYS570", "DYS607", "DYS442",
];
static YSEQ_DELTA: &[&str] = &[
    "DYR112",
    "DYS518",
    "DYS614",
    "DYS626",
    "DYS644",
    "DYS684",
    "DYS710",
    "DYS485",
    "DYS632",
    "DYS495",
    "DYS540",
    "DYS714",
    "DYS716",
    "DYS717",
    "DYS505",
    "DYS556",
    "DYS549",
    "DYS589",
    "DYS522",
    "DYS494",
    "DYS533",
    "DYS636",
    "DYS575",
    "DYS638",
    "DYS462",
    "DYS452",
    "DYS445",
    "Y-GATA-A10",
    "DYS463",
    "DYS441",
    "Y-GGAAT-1B07",
    "DYS525",
    "DYS712",
    "DYS593",
    "DYS650",
    "DYS532",
    "DYS715",
    "DYS504",
    "DYS513",
    "DYS561",
    "DYS552",
    "DYS726",
    "DYS635",
    "DYS587",
    "DYS643",
    "DYS497",
    "DYS510",
    "DYS434",
    "DYS461",
    "DYS435",
];
static YSEQ_GAMMA: &[&str] = &[
    "DYS728", "DYS723", "DYS711", "DYR76", "DYR33", "DYS727", "DYR157", "DYS713", "DYS531", "DYS578", "DYF395",
    "DYS590", "DYS537", "DYS641", "DYS472", "DYF406", "DYS511", "DYS557", "DYS490", "DYS446", "DYS481", "DYS413",
    "DYS534", "DYS450", "DYS425", "DYS594", "DYS444", "DYS520", "DYS436", "DYS565", "DYS572", "DYS617", "DYS568",
    "DYS487", "DYS640", "DYS492",
];
static YSEQ_EXCLUSIVE: &[&str] = &[
    "DYS728", "DYS723", "DYR112", "DYS711", "DYR76", "DYR33", "DYS727", "DYR157", "DYS713", "DYS518", "DYS614",
    "DYS626", "DYS644", "DYS684", "DYF397", "DYF399X", "DYS464X", "DYF408",
];
static YSEQ_KITTLER: &[&str] = &["DYS385A(K)", "DYS385B(K)"];

static PANELS: &[StrPanelDef] = &[
    StrPanelDef {
        id: "FTDNA_Y12",
        name: "Y-12",
        provider: "FTDNA",
        marketing_count: 12,
        actual_count: 11,
        order: 1,
        markers: FTDNA_Y12,
    },
    StrPanelDef {
        id: "FTDNA_Y25",
        name: "Y-25",
        provider: "FTDNA",
        marketing_count: 25,
        actual_count: 20,
        order: 2,
        markers: FTDNA_Y25,
    },
    StrPanelDef {
        id: "FTDNA_Y37",
        name: "Y-37",
        provider: "FTDNA",
        marketing_count: 37,
        actual_count: 30,
        order: 3,
        markers: FTDNA_Y37,
    },
    StrPanelDef {
        id: "FTDNA_Y67",
        name: "Y-67",
        provider: "FTDNA",
        marketing_count: 67,
        actual_count: 58,
        order: 4,
        markers: FTDNA_Y67,
    },
    StrPanelDef {
        id: "FTDNA_Y111",
        name: "Y-111",
        provider: "FTDNA",
        marketing_count: 111,
        actual_count: 102,
        order: 5,
        markers: FTDNA_Y111,
    },
    StrPanelDef {
        id: "YSEQ_ALPHA",
        name: "Alpha",
        provider: "YSEQ",
        marketing_count: 16,
        actual_count: 16,
        order: 1,
        markers: YSEQ_ALPHA,
    },
    StrPanelDef {
        id: "YSEQ_BETA",
        name: "Beta",
        provider: "YSEQ",
        marketing_count: 14,
        actual_count: 14,
        order: 2,
        markers: YSEQ_BETA,
    },
    StrPanelDef {
        id: "YSEQ_DELTA",
        name: "Delta",
        provider: "YSEQ",
        marketing_count: 50,
        actual_count: 50,
        order: 3,
        markers: YSEQ_DELTA,
    },
    StrPanelDef {
        id: "YSEQ_GAMMA",
        name: "Gamma",
        provider: "YSEQ",
        marketing_count: 36,
        actual_count: 36,
        order: 4,
        markers: YSEQ_GAMMA,
    },
    // Detection-only panels (order >= 100). They group and assign markers, but no tier badge uses
    // them.
    StrPanelDef {
        id: "YSEQ_EXCLUSIVE",
        name: "YSEQ-Exclusive",
        provider: "YSEQ",
        marketing_count: 18,
        actual_count: 18,
        order: 100,
        markers: YSEQ_EXCLUSIVE,
    },
    StrPanelDef {
        id: "YSEQ_KITTLER",
        name: "YSEQ-Kittler",
        provider: "YSEQ",
        marketing_count: 2,
        actual_count: 2,
        order: 101,
        markers: YSEQ_KITTLER,
    },
];

static PROVIDERS: &[ProviderDef] = &[
    ProviderDef {
        key: "FTDNA",
        display_name: "FamilyTreeDNA",
        cumulative: true,
        exclusive_markers: &[],
    },
    ProviderDef {
        key: "YSEQ",
        display_name: "YSEQ",
        cumulative: true,
        exclusive_markers: YSEQ_EXCLUSIVE,
    },
];

/// All panel definitions (FTDNA + YSEQ), in declaration order.
pub fn panels() -> &'static [StrPanelDef] {
    PANELS
}

/// The canonical FTDNA marker column order (Y-12 → Y-25 → Y-37 → Y-67 → Y-111), as used by the
/// FTDNA "Y-DNA Results Overview" chart. Each tier's markers follow, in tier order. This is the
/// stable left-to-right column order for the project-level STR overview.
pub fn ftdna_marker_order() -> Vec<&'static str> {
    let mut out: Vec<&'static str> = Vec::new();
    for tier in [FTDNA_Y12, FTDNA_Y25, FTDNA_Y37, FTDNA_Y67, FTDNA_Y111] {
        out.extend_from_slice(tier);
    }
    out
}

fn provider_def(key: &str) -> Option<&'static ProviderDef> {
    PROVIDERS.iter().find(|p| p.key == key)
}

/// Normalize a marker name for matching (case/space-insensitive).
pub fn norm(marker: &str) -> String {
    marker.trim().to_uppercase()
}

/// Canonicalize an arbitrary provider string to a known panel provider (defaults to FTDNA).
pub fn canonical_provider(provider: &str) -> &'static str {
    match provider.trim().to_uppercase().as_str() {
        "YSEQ" => "YSEQ",
        _ => "FTDNA",
    }
}

/// The normalized marker-name set for a profile's markers.
pub fn normalized_set(markers: &[StrMarker]) -> HashSet<String> {
    markers.iter().map(|m| norm(&m.marker)).collect()
}

/// Detect the provider from any present exclusive marker (else `None` → caller defaults to FTDNA).
pub fn detect_provider(markers: &HashSet<String>) -> Option<&'static str> {
    PROVIDERS.iter().find_map(|p| {
        (!p.exclusive_markers.is_empty() && p.exclusive_markers.iter().any(|m| markers.contains(&norm(m))))
            .then_some(p.key)
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelClassification {
    /// Highest tier reached, for example "Y-37". `None` if the profile met no tier threshold.
    pub panel_name: Option<String>,
    pub provider: &'static str,
    pub marker_count: usize,
    pub matched: usize,
}

/// Classify a profile into its highest reached tier for `provider` (auto-detected if `None`).
/// Mirrors `StrPanelService`. A cumulative provider needs ≥90% of the `actual_count` markers of a
/// tier. A provider that is not cumulative takes the highest-order tier with ≥80% overlap.
pub fn classify_panel(markers: &HashSet<String>, provider: Option<&str>) -> PanelClassification {
    let prov = provider
        .map(canonical_provider)
        .or_else(|| detect_provider(markers))
        .unwrap_or("FTDNA");
    let mut tiers: Vec<&StrPanelDef> = PANELS.iter().filter(|p| p.provider == prov && p.order < 100).collect();
    tiers.sort_by_key(|p| p.order);
    let cumulative = provider_def(prov).map(|d| d.cumulative).unwrap_or(true);

    let (name, matched) = if cumulative {
        let mut cum: HashSet<String> = HashSet::new();
        let mut best: Option<&StrPanelDef> = None;
        let mut best_matched = 0usize;
        for p in &tiers {
            cum.extend(p.markers.iter().map(|m| norm(m)));
            let match_count = markers.iter().filter(|m| cum.contains(*m)).count();
            let needed = (p.threshold() as f64 * 0.9) as usize;
            if match_count >= needed && match_count > best_matched {
                best = Some(p);
                best_matched = match_count;
            }
        }
        (best.map(|p| p.name.to_string()), best_matched)
    } else {
        let mut best: Option<&StrPanelDef> = None;
        let mut best_matched = 0usize;
        for p in &tiers {
            let pset: HashSet<String> = p.markers.iter().map(|m| norm(m)).collect();
            let overlap = markers.iter().filter(|m| pset.contains(*m)).count();
            if overlap >= (p.markers.len() as f64 * 0.8) as usize
                && p.order as usize >= best.map_or(0, |b| b.order as usize)
            {
                best = Some(p);
                best_matched = overlap;
            }
        }
        (best.map(|p| p.name.to_string()), best_matched)
    };

    PanelClassification {
        panel_name: name,
        provider: prov,
        marker_count: markers.len(),
        matched,
    }
}

/// A "reached" flag for each tier, for the summary badges: `(tier name, filled)`. A tier is full
/// when the distinct-marker count of the profile meets its threshold. This drops detection-only
/// panels (order ≥100).
pub fn tier_badges(provider: &str, marker_count: usize) -> Vec<(String, bool)> {
    let prov = canonical_provider(provider);
    let mut tiers: Vec<&StrPanelDef> = PANELS.iter().filter(|p| p.provider == prov && p.order < 100).collect();
    tiers.sort_by_key(|p| p.order);
    tiers
        .iter()
        .map(|p| (p.name.to_string(), marker_count as u32 >= p.threshold()))
        .collect()
}

/// Group the markers of a profile into tiers (FTDNA-style report layout). Each marker goes to the
/// first tier (by order) that lists it, and the rest go to [`EXTENDED`]. Returns
/// `(tier name, markers)` in tier order, with the Extended group last, and each group not empty.
/// Marker order inside a tier follows the order of the profile.
pub fn assign_markers_to_panels<'a>(markers: &'a [StrMarker], provider: &str) -> Vec<(String, Vec<&'a StrMarker>)> {
    let prov = canonical_provider(provider);
    let mut tiers: Vec<&StrPanelDef> = PANELS.iter().filter(|p| p.provider == prov).collect();
    tiers.sort_by_key(|p| p.order);

    // marker -> the name of the tier that owns it (the first tier by order that lists it).
    let mut owner: HashMap<String, &'static str> = HashMap::new();
    for p in &tiers {
        for m in p.markers {
            owner.entry(norm(m)).or_insert(p.name);
        }
    }

    // Ordered tier labels (the badge tiers in order, then Extended).
    let mut order_index: Vec<&'static str> = tiers.iter().filter(|p| p.order < 100).map(|p| p.name).collect();
    // Include detection panels' names if any markers land there.
    for p in tiers.iter().filter(|p| p.order >= 100) {
        order_index.push(p.name);
    }
    order_index.push(EXTENDED);

    let mut groups: HashMap<&str, Vec<&'a StrMarker>> = HashMap::new();
    for m in markers {
        let tier = owner.get(&norm(&m.marker)).copied().unwrap_or(EXTENDED);
        groups.entry(tier).or_default().push(m);
    }

    order_index
        .into_iter()
        .filter_map(|name| groups.remove(name).map(|ms| (name.to_string(), ms)))
        .filter(|(_, ms)| !ms.is_empty())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strprofile::StrMarker;

    fn marks(names: &[&str]) -> Vec<StrMarker> {
        names
            .iter()
            .map(|n| StrMarker {
                marker: n.to_string(),
                value: "13".to_string(),
            })
            .collect()
    }

    #[test]
    fn classify_ftdna_tiers() {
        // Y-12 + Y-25 + Y-37 additions = 30 distinct keys → Y-37.
        let all: Vec<&str> = FTDNA_Y12.iter().chain(FTDNA_Y25).chain(FTDNA_Y37).copied().collect();
        let set = normalized_set(&marks(&all));
        let c = classify_panel(&set, Some("FTDNA"));
        assert_eq!(c.provider, "FTDNA");
        assert_eq!(c.panel_name.as_deref(), Some("Y-37"));

        // Full 102 → Y-111.
        let full: Vec<&str> = FTDNA_Y12
            .iter()
            .chain(FTDNA_Y25)
            .chain(FTDNA_Y37)
            .chain(FTDNA_Y67)
            .chain(FTDNA_Y111)
            .copied()
            .collect();
        let c2 = classify_panel(&normalized_set(&marks(&full)), Some("FTDNA"));
        assert_eq!(c2.panel_name.as_deref(), Some("Y-111"));

        // Only Y-12 markers → at most Y-12 (not higher).
        let c3 = classify_panel(&normalized_set(&marks(FTDNA_Y12)), Some("FTDNA"));
        assert_eq!(c3.panel_name.as_deref(), Some("Y-12"));
    }

    #[test]
    fn detects_yseq_from_exclusive_marker() {
        // A set that holds a YSEQ-exclusive marker auto-detects provider YSEQ.
        let mut names = FTDNA_Y12.to_vec();
        names.push("DYS728"); // YSEQ exclusive
        let set = normalized_set(&marks(&names));
        assert_eq!(detect_provider(&set), Some("YSEQ"));
        let c = classify_panel(&set, None);
        assert_eq!(c.provider, "YSEQ");
    }

    #[test]
    fn tier_badges_fill_by_count() {
        let badges = tier_badges("FTDNA", 30); // exactly Y-37 threshold
        let map: HashMap<_, _> = badges.into_iter().collect();
        assert!(map["Y-12"]);
        assert!(map["Y-37"]);
        assert!(!map["Y-67"]);
        assert!(!map["Y-111"]);
    }

    #[test]
    fn assign_groups_by_tier_and_extended() {
        let mut names = vec!["DYS393", "DYS390"]; // Y-12
        names.push("DYS458"); // Y-25
        names.push("FTYXYZ"); // unknown → Extended
        let m = marks(&names);
        let groups = assign_markers_to_panels(&m, "FTDNA");
        let by: HashMap<String, Vec<String>> = groups
            .iter()
            .map(|(t, ms)| (t.clone(), ms.iter().map(|x| x.marker.clone()).collect()))
            .collect();
        assert_eq!(by["Y-12"], vec!["DYS393", "DYS390"]);
        assert_eq!(by["Y-25"], vec!["DYS458"]);
        assert_eq!(by[EXTENDED], vec!["FTYXYZ"]);
        // Tier order preserved: Y-12 before Y-25 before Extended.
        let order: Vec<&str> = groups.iter().map(|(t, _)| t.as_str()).collect();
        assert_eq!(order, vec!["Y-12", "Y-25", EXTENDED]);
    }
}
