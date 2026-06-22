//! Aggregation for the FTDNA-style project "Y-DNA Results Overview" chart: per-subgroup, per-marker
//! MIN / MAX / MODE statistics and per-cell deviation from the modal value (the colour coding).
//!
//! Marker values are kept as text (a multi-copy marker like DYS385 reports "11-15", DYS464 reports
//! "14-15-16-17", CDY "37-37"). For ordering we parse a value into its sorted allele tuple and
//! compare tuples; for the modal value we count the canonical (sorted) string. Non-numeric or null
//! values ("-", "") are ignored.

/// Parse an STR marker value into its sorted allele tuple, e.g. "11-15" → [11, 15], "13" → [13].
/// Returns `None` for null/non-numeric values so callers can skip them.
pub fn parse_allele(value: &str) -> Option<Vec<i32>> {
    let v = value.trim();
    if v.is_empty() || v == "-" {
        return None;
    }
    let mut parts: Vec<i32> = Vec::new();
    for p in v.split('-') {
        parts.push(p.trim().parse::<i32>().ok()?);
    }
    if parts.is_empty() {
        return None;
    }
    parts.sort_unstable();
    Some(parts)
}

/// The canonical (sorted) form of a value used as the modal key, so order-independent multi-copy
/// values ("16-15" and "15-16") collapse to one. Falls back to the trimmed original when unparseable.
pub fn canonical(value: &str) -> String {
    match parse_allele(value) {
        Some(parts) => parts.iter().map(|n| n.to_string()).collect::<Vec<_>>().join("-"),
        None => value.trim().to_string(),
    }
}

/// MIN / MAX / MODE for one marker column across a subgroup. Each is the canonical value string;
/// `None` when no member reported a (numeric) value.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MarkerStats {
    pub min: Option<String>,
    pub max: Option<String>,
    pub mode: Option<String>,
}

/// Summarise one marker column over a set of member values: numeric MIN/MAX (by sorted-tuple order)
/// and the MODE (most frequent canonical value, ties → smallest tuple).
pub fn marker_stats<'a, I>(values: I) -> MarkerStats
where
    I: IntoIterator<Item = &'a str>,
{
    use std::collections::HashMap;
    let parsed: Vec<(Vec<i32>, String)> = values
        .into_iter()
        .filter_map(|v| parse_allele(v).map(|t| (t.clone(), tuple_str(&t))))
        .collect();
    if parsed.is_empty() {
        return MarkerStats::default();
    }
    let min = parsed.iter().min_by(|a, b| a.0.cmp(&b.0)).map(|(_, s)| s.clone());
    let max = parsed.iter().max_by(|a, b| a.0.cmp(&b.0)).map(|(_, s)| s.clone());
    // Mode: most frequent canonical string; ties broken by the smaller tuple for determinism.
    let mut counts: HashMap<&str, (usize, &Vec<i32>)> = HashMap::new();
    for (tuple, s) in &parsed {
        let e = counts.entry(s.as_str()).or_insert((0, tuple));
        e.0 += 1;
    }
    let mode = counts
        .iter()
        .max_by(|a, b| a.1 .0.cmp(&b.1 .0).then_with(|| b.1 .1.cmp(a.1 .1)))
        .map(|(s, _)| (*s).to_string());
    MarkerStats { min, max, mode }
}

fn tuple_str(t: &[i32]) -> String {
    t.iter().map(|n| n.to_string()).collect::<Vec<_>>().join("-")
}

/// How a cell's value relates to its subgroup's modal value — drives the colour coding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Deviation {
    /// Equal to (or order-equivalent to) the mode, or no mode / unparseable.
    None,
    /// Strictly below the modal value (fewer repeats).
    Below,
    /// Strictly above the modal value (more repeats).
    Above,
    /// Differs from the mode but isn't strictly orderable (multi-copy with mixed direction).
    Differs,
}

/// Classify a single cell value against the column's modal value.
pub fn deviation(value: &str, mode: &Option<String>) -> Deviation {
    let (Some(val), Some(mode_s)) = (parse_allele(value), mode.as_deref()) else {
        return Deviation::None;
    };
    let Some(modev) = parse_allele(mode_s) else {
        return Deviation::None;
    };
    if val == modev {
        return Deviation::None;
    }
    if val.len() == modev.len() {
        // Element-wise: all below → Below, all above → Above, else mixed → Differs.
        let all_le = val.iter().zip(&modev).all(|(a, b)| a <= b);
        let all_ge = val.iter().zip(&modev).all(|(a, b)| a >= b);
        match (all_le, all_ge) {
            (true, false) => Deviation::Below,
            (false, true) => Deviation::Above,
            _ => Deviation::Differs,
        }
    } else {
        Deviation::Differs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_single_and_multi() {
        assert_eq!(parse_allele("13"), Some(vec![13]));
        assert_eq!(parse_allele("16-15"), Some(vec![15, 16]));
        assert_eq!(parse_allele("-"), None);
        assert_eq!(parse_allele(""), None);
        assert_eq!(parse_allele("n/a"), None);
    }

    #[test]
    fn multi_copy_is_order_independent() {
        assert_eq!(canonical("16-15"), "15-16");
        assert_eq!(canonical("15-16"), "15-16");
    }

    #[test]
    fn stats_min_max_mode() {
        let s = marker_stats(["13", "13", "14", "12"]);
        assert_eq!(s.min.as_deref(), Some("12"));
        assert_eq!(s.max.as_deref(), Some("14"));
        assert_eq!(s.mode.as_deref(), Some("13"));
    }

    #[test]
    fn deviation_direction() {
        let mode = Some("13".to_string());
        assert_eq!(deviation("13", &mode), Deviation::None);
        assert_eq!(deviation("12", &mode), Deviation::Below);
        assert_eq!(deviation("14", &mode), Deviation::Above);
        let multimode = Some("11-15".to_string());
        assert_eq!(deviation("11-14", &multimode), Deviation::Below);
        assert_eq!(deviation("12-16", &multimode), Deviation::Above);
        assert_eq!(deviation("10-16", &multimode), Deviation::Differs);
    }
}
