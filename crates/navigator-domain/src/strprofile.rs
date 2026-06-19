//! Y-STR profiles — a subject's short-tandem-repeat marker calls (e.g. DYS393=13),
//! grouped by the panel/test that produced them (Y-37, Big Y-700 STRs, …). A pragmatic
//! port of the Scala `StrProfile`: marker values + panel provenance, without the AT-URI/
//! sync/derivation metadata (added later if needed). Types are pure; [`parse_csv`] turns
//! exported marker tables (FTDNA/YSEQ-style) into markers without touching the filesystem.

use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// One STR marker call. `value` is kept as text because multi-copy markers report several
/// alleles (e.g. "16-17" for DYS385) and palindromic markers can carry "-"/null.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StrMarker {
    pub marker: String,
    pub value: String,
}

/// One marker's donor consensus across all of a subject's panels: the modal value, how many
/// panels reported it, and whether the panels disagreed.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ConsensusStrMarker {
    pub marker: String,
    pub value: String,
    /// Panels that reported a (non-null) value for this marker.
    pub panels: usize,
    /// True when panels reported differing values (a conflict to surface).
    pub conflict: bool,
}

/// Merge a subject's STR panels into one **donor consensus** profile: per marker, the modal
/// (most common) value across panels, flagged when panels disagree. Null/palindromic-null
/// values are skipped. Markers keep their first-seen order.
pub fn consensus_markers(profiles: &[StrProfile]) -> Vec<ConsensusStrMarker> {
    use std::collections::HashMap;
    let mut order: Vec<String> = Vec::new();
    let mut by_marker: HashMap<String, Vec<String>> = HashMap::new();
    for p in profiles {
        for m in &p.markers {
            let v = m.value.trim();
            if v.is_empty() || v == "-" {
                continue;
            }
            if !by_marker.contains_key(&m.marker) {
                order.push(m.marker.clone());
            }
            by_marker.entry(m.marker.clone()).or_default().push(v.to_string());
        }
    }
    order
        .into_iter()
        .map(|marker| {
            let vals = &by_marker[&marker];
            let mut counts: HashMap<&str, usize> = HashMap::new();
            for v in vals {
                *counts.entry(v.as_str()).or_default() += 1;
            }
            let conflict = counts.len() > 1;
            // Modal value; ties broken by the lexicographically smallest for determinism.
            let value = counts
                .iter()
                .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
                .map(|(v, _)| (*v).to_string())
                .unwrap_or_default();
            ConsensusStrMarker {
                marker,
                value,
                panels: vals.len(),
                conflict,
            }
        })
        .collect()
}

/// A subject's Y-STR profile from one panel/source.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StrProfile {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    /// Panel name (one of [`KNOWN_PANELS`] or custom), e.g. "Y-37".
    pub panel_name: String,
    /// Testing company / source (one of [`KNOWN_PROVIDERS`]).
    pub provider: Option<String>,
    /// How the STRs were obtained (one of [`KNOWN_SOURCES`]).
    pub source: Option<String>,
    pub markers: Vec<StrMarker>,
}

/// Fields for creating an STR profile (the store assigns the id).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NewStrProfile {
    pub biosample_guid: SampleGuid,
    pub panel_name: String,
    pub provider: Option<String>,
    pub source: Option<String>,
    pub markers: Vec<StrMarker>,
}

/// Common Y-STR panel names (for the import form's dropdown). "CUSTOM" for anything else.
pub const KNOWN_PANELS: &[&str] = &[
    "Y-12",
    "Y-25",
    "Y-37",
    "Y-67",
    "Y-111",
    "Y-500",
    "Y-700",
    "YSEQ_ALPHA",
    "CUSTOM",
];

/// Known testing companies / sources.
pub const KNOWN_PROVIDERS: &[&str] = &["FTDNA", "YSEQ", "NEBULA", "DANTE", "WGS_DERIVED", "OTHER"];

/// How a profile's STRs were obtained.
pub const KNOWN_SOURCES: &[&str] = &[
    "DIRECT_TEST",
    "WGS_DERIVED",
    "BIG_Y_DERIVED",
    "IMPORTED",
    "MANUAL_ENTRY",
];

/// Trim whitespace and one layer of surrounding double-quotes from a cell (FTDNA/YSEQ pad
/// values like `" 13"`), then trim again.
fn clean_cell(s: &str) -> &str {
    s.trim().trim_matches('"').trim()
}

/// A value is "missing" when it's blank or a placeholder dash.
fn is_blank_value(v: &str) -> bool {
    v.is_empty() || v == "-"
}

/// A marker-name row reads as names: most non-empty cells contain a letter (DYS393, FTY10,
/// Y-GATA-H4, YCAII, CDY…).
fn looks_like_names(cells: &[&str]) -> bool {
    let non_empty: Vec<&&str> = cells.iter().filter(|c| !c.is_empty()).collect();
    if non_empty.is_empty() {
        return false;
    }
    let with_letter = non_empty
        .iter()
        .filter(|c| c.bytes().any(|b| b.is_ascii_alphabetic()))
        .count();
    with_letter * 10 >= non_empty.len() * 8 // ≥80%
}

/// A value row reads as STR allele values: most non-empty cells are digits, with optional
/// `-`/`.`/`/` (multi-copy like `11-15`, microvariants like `10.2`).
fn looks_like_values(cells: &[&str]) -> bool {
    let non_empty: Vec<&&str> = cells.iter().filter(|c| !is_blank_value(c)).collect();
    if non_empty.is_empty() {
        return false;
    }
    let numeric = non_empty
        .iter()
        .filter(|c| {
            c.bytes().any(|b| b.is_ascii_digit())
                && c.bytes().all(|b| b.is_ascii_digit() || matches!(b, b'-' | b'.' | b'/'))
        })
        .count();
    numeric * 10 >= non_empty.len() * 8 // ≥80%
}

/// Parse an exported STR marker table into markers. Two layouts are accepted:
///
/// * **Tall**: one `marker,value` (a.k.a. `locus`/`allele`/`result`) row per line, with or
///   without a header row.
/// * **Wide**: the FTDNA / YSEQ export shape — a single header row of marker names and a
///   single parallel row of values (often quoted and space-padded, e.g. `" 13"`).
///
/// Comma- or tab-separated; blank lines and leading `#` comments are ignored; markers with an
/// empty/`-` value are skipped. Errors only if no usable markers are found.
pub fn parse_csv(text: &str) -> Result<Vec<StrMarker>, String> {
    let content: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();

    // Wide layout: exactly two rows — a header of marker names and a parallel row of values.
    // Guard against a 2-row tall file by requiring several columns and confirming the first
    // row reads as names (mostly letters) and the second as values (mostly digits/dashes).
    if content.len() == 2 {
        let sep = if content[0].contains('\t') { '\t' } else { ',' };
        let names: Vec<&str> = content[0].split(sep).map(clean_cell).collect();
        let values: Vec<&str> = content[1].split(sep).map(clean_cell).collect();
        if names.len() >= 5 && values.len() >= 5 && looks_like_names(&names) && looks_like_values(&values) {
            let markers: Vec<StrMarker> = names
                .iter()
                .zip(values.iter())
                .filter(|(name, value)| !name.is_empty() && !is_blank_value(value))
                .map(|(name, value)| StrMarker {
                    marker: (*name).to_string(),
                    value: (*value).to_string(),
                })
                .collect();
            if !markers.is_empty() {
                return Ok(markers);
            }
        }
    }

    // Tall layout: one marker per row.
    let mut markers = Vec::new();
    let mut header_checked = false;
    for line in content {
        let sep = if line.contains('\t') { '\t' } else { ',' };
        let mut cols = line.splitn(3, sep).map(clean_cell);
        let (Some(marker), Some(value)) = (cols.next(), cols.next()) else {
            continue;
        };
        // Skip a header row if the first non-comment line looks like column titles.
        if !header_checked {
            header_checked = true;
            let m = marker.to_ascii_lowercase();
            let v = value.to_ascii_lowercase();
            let is_header = matches!(m.as_str(), "marker" | "locus" | "dys" | "name")
                || matches!(v.as_str(), "value" | "allele" | "result" | "alleles");
            if is_header {
                continue;
            }
        }
        if marker.is_empty() || is_blank_value(value) {
            continue;
        }
        markers.push(StrMarker {
            marker: marker.to_string(),
            value: value.to_string(),
        });
    }
    if markers.is_empty() {
        return Err("no STR markers found (expected `marker,value` rows or a wide FTDNA/YSEQ table)".into());
    }
    Ok(markers)
}

/// Y-STR distance between two profiles: differing marker values over markers present in
/// both. Returns (differing, compared). Distance 0 over many shared markers is consistent
/// with a shared paternal line (identity corroboration).
pub fn str_distance(a: &[StrMarker], b: &[StrMarker]) -> (i64, i64) {
    let mut differing = 0;
    let mut compared = 0;
    for ma in a {
        if let Some(mb) = b.iter().find(|m| m.marker.eq_ignore_ascii_case(&ma.marker)) {
            compared += 1;
            if !ma.value.eq_ignore_ascii_case(&mb.value) {
                differing += 1;
            }
        }
    }
    (differing, compared)
}

/// One marker where providers disagree, with each provider's reported value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MarkerConflict {
    pub marker: String,
    /// `(provider, value)` per provider, provider order as first seen.
    pub by_provider: Vec<(String, String)>,
}

/// Cross-provider comparison of a subject's STR profiles (e.g. FTDNA vs YSEQ).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StrComparison {
    /// Markers reported by ≥2 providers with disagreeing values.
    pub conflicts: Vec<MarkerConflict>,
    /// Markers reported by ≥2 providers that agree.
    pub agreement_count: usize,
    /// Distinct provider labels seen.
    pub providers: Vec<String>,
}

/// Whether two STR values represent the same allele(s). Multi-copy values (`"16-15"`) are compared
/// order-independently (`"16-15"` ≡ `"15-16"`); everything else is a trimmed string match.
pub fn values_match(a: &str, b: &str) -> bool {
    let norm = |s: &str| {
        let mut parts: Vec<&str> = s.split('-').map(|p| p.trim()).filter(|p| !p.is_empty()).collect();
        parts.sort_unstable();
        parts.join("-")
    };
    norm(a) == norm(b)
}

/// Compare a subject's STR profiles across providers: flag markers where providers disagree, count
/// agreements, and list the providers. Mirrors the Scala `StrMarkerComparator.compare`. A profile's
/// provider defaults to `"UNKNOWN"` when unset; multiple profiles from the same provider collapse to
/// that provider's first-seen value for a marker.
pub fn compare_profiles(profiles: &[StrProfile]) -> StrComparison {
    // Normalized marker -> ordered list of (provider, value), one entry per provider.
    let mut by_marker: Vec<(String, Vec<(String, String)>)> = Vec::new();
    let mut index: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut providers: Vec<String> = Vec::new();

    for p in profiles {
        let provider = p.provider.clone().unwrap_or_else(|| "UNKNOWN".to_string());
        if !providers.contains(&provider) {
            providers.push(provider.clone());
        }
        for mk in &p.markers {
            let key = mk.marker.trim().to_uppercase();
            let slot = *index.entry(key.clone()).or_insert_with(|| {
                by_marker.push((mk.marker.clone(), Vec::new()));
                by_marker.len() - 1
            });
            let entries = &mut by_marker[slot].1;
            if !entries.iter().any(|(prov, _)| prov == &provider) {
                entries.push((provider.clone(), mk.value.clone()));
            }
        }
    }

    let mut conflicts = Vec::new();
    let mut agreement_count = 0;
    for (marker, entries) in by_marker {
        if entries.len() < 2 {
            continue; // need ≥2 providers to agree or conflict
        }
        let first = &entries[0].1;
        if entries.iter().all(|(_, v)| values_match(first, v)) {
            agreement_count += 1;
        } else {
            conflicts.push(MarkerConflict {
                marker,
                by_provider: entries,
            });
        }
    }
    conflicts.sort_by(|a, b| a.marker.cmp(&b.marker));

    StrComparison {
        conflicts,
        agreement_count,
        providers,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(marker: &str, value: &str) -> StrMarker {
        StrMarker {
            marker: marker.into(),
            value: value.into(),
        }
    }

    #[test]
    fn str_distance_over_shared_markers() {
        let a = [m("DYS393", "13"), m("DYS390", "24"), m("DYS19", "14")];
        let b = [m("DYS393", "13"), m("DYS390", "25"), m("DYS388", "12")]; // DYS19 absent, DYS388 extra
        assert_eq!(str_distance(&a, &b), (1, 2)); // DYS393 match, DYS390 differ; only 2 shared
    }

    #[test]
    fn parses_csv_with_header() {
        let csv = "Marker,Value\nDYS393,13\nDYS390,24\nDYS385,11-14\n";
        let m = parse_csv(csv).unwrap();
        assert_eq!(m.len(), 3);
        assert_eq!(
            m[0],
            StrMarker {
                marker: "DYS393".into(),
                value: "13".into()
            }
        );
        assert_eq!(m[2].value, "11-14"); // multi-copy preserved
    }

    #[test]
    fn parses_tsv_without_header_and_skips_blanks_and_nulls() {
        let tsv = "# YSEQ export\nDYS393\t13\n\nDYS438\t-\nDYS439\t11\n";
        let m = parse_csv(tsv).unwrap();
        assert_eq!(
            m.iter().map(|x| x.marker.as_str()).collect::<Vec<_>>(),
            ["DYS393", "DYS439"]
        );
    }

    #[test]
    fn empty_input_errors() {
        assert!(parse_csv("\n\n# nothing\n").is_err());
    }

    #[test]
    fn parses_wide_ftdna_layout_with_quotes_and_padding() {
        // FTDNA/YSEQ shape: a row of marker names + a parallel row of quoted, space-padded
        // values; multi-copy markers stay dash-joined; empty cells are skipped.
        let csv = "DYS393,DYS390,DYS385,DYS459,DYS464\n\" 13\",\" 24\",\" 11-15\",\" \",\" 14-15-17-17\"\n";
        let m = parse_csv(csv).unwrap();
        assert_eq!(m.len(), 4); // DYS459 (blank) skipped
        assert_eq!(
            m[0],
            StrMarker {
                marker: "DYS393".into(),
                value: "13".into()
            }
        );
        assert_eq!(
            m[2],
            StrMarker {
                marker: "DYS385".into(),
                value: "11-15".into()
            }
        );
        assert_eq!(
            m[3],
            StrMarker {
                marker: "DYS464".into(),
                value: "14-15-17-17".into()
            }
        );
    }

    #[test]
    fn two_row_tall_file_is_not_mistaken_for_wide() {
        // Two data rows, two columns each — still tall (the wide path needs ≥5 name columns).
        let csv = "DYS393,13\nDYS390,24\n";
        let m = parse_csv(csv).unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(
            m[0],
            StrMarker {
                marker: "DYS393".into(),
                value: "13".into()
            }
        );
        assert_eq!(
            m[1],
            StrMarker {
                marker: "DYS390".into(),
                value: "24".into()
            }
        );
    }

    #[test]
    fn consensus_merges_and_flags_conflicts() {
        let guid = SampleGuid(uuid::Uuid::nil());
        let mk = |panel: &str, pairs: &[(&str, &str)]| StrProfile {
            id: 0,
            biosample_guid: guid,
            panel_name: panel.into(),
            provider: None,
            source: None,
            markers: pairs
                .iter()
                .map(|(m, v)| StrMarker {
                    marker: (*m).into(),
                    value: (*v).into(),
                })
                .collect(),
        };
        let profiles = vec![
            mk("Y-37", &[("DYS393", "13"), ("DYS390", "24"), ("DYS385", "-")]),
            mk("Y-111", &[("DYS393", "13"), ("DYS390", "25"), ("DYS19", "14")]),
        ];
        let c = consensus_markers(&profiles);
        let by = |name: &str| c.iter().find(|m| m.marker == name).cloned();

        // Agreement: one value, 2 panels, no conflict.
        let d393 = by("DYS393").unwrap();
        assert_eq!(d393.value, "13");
        assert_eq!(d393.panels, 2);
        assert!(!d393.conflict);

        // Disagreement: flagged conflict.
        assert!(by("DYS390").unwrap().conflict);

        // Null ("-") is skipped; single-panel marker still appears.
        assert!(by("DYS385").is_none());
        assert_eq!(by("DYS19").unwrap().panels, 1);
    }

    #[test]
    fn values_match_handles_multicopy_order() {
        assert!(values_match("16-15", "15-16")); // order-independent
        assert!(values_match("14-15-17-17", "17-14-17-15"));
        assert!(values_match("13", "13"));
        assert!(!values_match("13", "14"));
        assert!(!values_match("16-15", "16-16"));
    }

    #[test]
    fn compare_flags_conflicts_and_agreements() {
        let guid = SampleGuid(uuid::Uuid::nil());
        let mk = |provider: &str, pairs: &[(&str, &str)]| StrProfile {
            id: 0,
            biosample_guid: guid,
            panel_name: "X".into(),
            provider: Some(provider.into()),
            source: None,
            markers: pairs
                .iter()
                .map(|(m, v)| StrMarker {
                    marker: (*m).into(),
                    value: (*v).into(),
                })
                .collect(),
        };
        let profiles = vec![
            mk(
                "FTDNA",
                &[("DYS393", "13"), ("DYS390", "24"), ("DYS385", "11-15"), ("DYS19", "14")],
            ),
            mk("YSEQ", &[("DYS393", "13"), ("DYS390", "25"), ("DYS385", "15-11")]), // DYS390 differs; DYS385 multi-copy reorder
        ];
        let c = compare_profiles(&profiles);

        assert_eq!(c.providers, vec!["FTDNA".to_string(), "YSEQ".to_string()]);
        // DYS390 is the only conflict; DYS393 + DYS385 agree (DYS385 order-independent).
        assert_eq!(c.conflicts.len(), 1);
        assert_eq!(c.conflicts[0].marker, "DYS390");
        assert_eq!(
            c.conflicts[0].by_provider,
            vec![("FTDNA".into(), "24".into()), ("YSEQ".into(), "25".into())]
        );
        assert_eq!(c.agreement_count, 2); // DYS393, DYS385
                                          // DYS19 is single-provider → neither conflict nor agreement.
    }
}
