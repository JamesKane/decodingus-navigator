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
pub const KNOWN_PANELS: &[&str] =
    &["Y-12", "Y-25", "Y-37", "Y-67", "Y-111", "Y-500", "Y-700", "YSEQ_ALPHA", "CUSTOM"];

/// Known testing companies / sources.
pub const KNOWN_PROVIDERS: &[&str] = &["FTDNA", "YSEQ", "NEBULA", "DANTE", "WGS_DERIVED", "OTHER"];

/// How a profile's STRs were obtained.
pub const KNOWN_SOURCES: &[&str] =
    &["DIRECT_TEST", "WGS_DERIVED", "BIG_Y_DERIVED", "IMPORTED", "MANUAL_ENTRY"];

/// Parse an exported STR marker table into markers. Accepts comma- or tab-separated text
/// with a `marker,value` (a.k.a. `locus`/`allele`/`result`) layout, with or without a
/// header row; blank lines and a leading `#` comment are ignored. Markers with an empty
/// value are skipped. Errors only if no usable marker rows are found.
pub fn parse_csv(text: &str) -> Result<Vec<StrMarker>, String> {
    let mut markers = Vec::new();
    let mut header_checked = false;
    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let sep = if line.contains('\t') { '\t' } else { ',' };
        let mut cols = line.splitn(3, sep).map(str::trim);
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
        if marker.is_empty() || value.is_empty() || value == "-" {
            continue;
        }
        markers.push(StrMarker { marker: marker.to_string(), value: value.to_string() });
    }
    if markers.is_empty() {
        return Err("no STR markers found (expected `marker,value` rows)".into());
    }
    Ok(markers)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_csv_with_header() {
        let csv = "Marker,Value\nDYS393,13\nDYS390,24\nDYS385,11-14\n";
        let m = parse_csv(csv).unwrap();
        assert_eq!(m.len(), 3);
        assert_eq!(m[0], StrMarker { marker: "DYS393".into(), value: "13".into() });
        assert_eq!(m[2].value, "11-14"); // multi-copy preserved
    }

    #[test]
    fn parses_tsv_without_header_and_skips_blanks_and_nulls() {
        let tsv = "# YSEQ export\nDYS393\t13\n\nDYS438\t-\nDYS439\t11\n";
        let m = parse_csv(tsv).unwrap();
        assert_eq!(m.iter().map(|x| x.marker.as_str()).collect::<Vec<_>>(), ["DYS393", "DYS439"]);
    }

    #[test]
    fn empty_input_errors() {
        assert!(parse_csv("\n\n# nothing\n").is_err());
    }
}
