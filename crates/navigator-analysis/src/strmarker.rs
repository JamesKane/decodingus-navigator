//! FTDNA Y-STR convention layer — map the caller's HipSTR-locus repeat counts to FTDNA marker values.
//!
//! The HipSTR reference already names the loci (DYS393, DYS19/DYS394, …; see [`crate::strref`]), so
//! the caller emits DYS names directly. What differs is the **counting convention**: FTDNA reports a
//! per-marker value that is the caller's repeat count plus a fixed offset (0 for most, ±1–3 for some),
//! and a set of markers whose HipSTR tract doesn't correspond 1:1 to the FTDNA marker (large tract
//! mismatches, plus multi-copy/nested markers like DYS385/DYS464/DYS389II) which can't be mapped by a
//! single offset.
//!
//! The offset table below was **calibrated against a 14-kit Big Y corpus** (R-CTS4466 project): each
//! offset is the value the caller's count differed from the kit's FTDNA value, *constant across all
//! kits where the marker varied* (100% agreement — the corpus is what distinguishes a real convention
//! from one sample's variation). The harness that derived it is `examples/str_calibrate.rs`; re-run it
//! over more/diverse kits to extend the table (markers not yet calibrated report `Uncalibrated`).

use crate::strcaller::{StrConfidence, StrGenotype};

/// FTDNA value = caller repeat count + offset. Covers the calibrated single-copy markers (offset 0 =
/// "reliable", ±1–3 = a real convention). Markers absent here are either in [`EXCLUDE`] or uncalibrated.
static OFFSETS: &[(&str, i32)] = &[
    ("DYF406S1", 0), ("DYS19", -1), ("DYS388", 0), ("DYS389I", 0), ("DYS390", 0), ("DYS391", -1),
    ("DYS392", 0), ("DYS393", 0), ("DYS426", 0), ("DYS434", 0), ("DYS435", 2), ("DYS436", 0),
    ("DYS438", 2), ("DYS442", -3), ("DYS445", 0), ("DYS446", 0), ("DYS453", 0), ("DYS454", 0),
    ("DYS455", 0), ("DYS456", 0), ("DYS458", 0), ("DYS461", 1), ("DYS462", 0), ("DYS472", 0),
    ("DYS474", -3), ("DYS476", 0), ("DYS477", 0), ("DYS480", 0), ("DYS481", 0), ("DYS485", -1),
    ("DYS487", 0), ("DYS488", 0), ("DYS490", 0), ("DYS492", 0), ("DYS494", 0), ("DYS495", 1),
    ("DYS497", 0), ("DYS499", 0), ("DYS505", 0), ("DYS508", 0), ("DYS511", 0), ("DYS512", -3),
    ("DYS520", -2), ("DYS522", 1), ("DYS525", 0), ("DYS530", 0), ("DYS531", 0), ("DYS533", 0),
    ("DYS537", 0), ("DYS538", 1), ("DYS539", 0), ("DYS549", 0), ("DYS556", 0), ("DYS559", 1),
    ("DYS561", 0), ("DYS565", 0), ("DYS567", 0), ("DYS568", 0), ("DYS569", 0), ("DYS570", 0),
    ("DYS573", 0), ("DYS574", 0), ("DYS575", 0), ("DYS576", 0), ("DYS577", 0), ("DYS578", 0),
    ("DYS579", 0), ("DYS580", 0), ("DYS581", 0), ("DYS583", 0), ("DYS584", 0), ("DYS585", -3),
    ("DYS587", 1), ("DYS588", -3), ("DYS589", -1), ("DYS590", 0), ("DYS593", 0), ("DYS594", 0),
    ("DYS595", 0), ("DYS615", -2), ("DYS618", 0), ("DYS620", 0), ("DYS621", 0), ("DYS629", -3),
    ("DYS632", -2), ("DYS635", 0), ("DYS638", 0), ("DYS640", 0), ("DYS641", 0), ("DYS645", 0),
    ("Y-GATA-A10", 0),
];

/// Markers whose HipSTR tract can't be mapped to the FTDNA value by a single offset: large tract
/// mismatches and multi-copy/nested markers (DYS385/DYS464 split sub-loci, DYS389II nesting). Reported
/// as `Excluded` — the enclosing-read caller doesn't yield a vendor-comparable value here (yet).
static EXCLUDE: &[&str] = &[
    "DYS389II", "DYS448", "DYS449", "DYS450", "DYS460", "DYS470", "DYS475", "DYS484", "DYS502",
    "DYS504", "DYS510", "DYS513", "DYS516", "DYS532", "DYS534", "DYS540", "DYS541", "DYS543",
    "DYS544", "DYS551", "DYS552", "DYS557", "DYS607", "DYS616", "DYS624", "DYS631", "DYS637",
    "DYS717", "Y-GATA-H4",
];

/// How confidently a caller locus maps to an FTDNA marker value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MarkerStatus {
    /// Calibrated, no convention offset (caller count == FTDNA value).
    Reliable,
    /// Calibrated with a fixed FTDNA convention offset (±1–3).
    ConventionOffset,
    /// Tract mismatch / multi-copy / nested — no 1:1 mapping; value is the raw caller count.
    Excluded,
    /// Not in the calibration table — value is the raw caller count, pending calibration.
    Uncalibrated,
}

/// One marker called from sequence, expressed in the FTDNA convention.
#[derive(Debug, Clone, PartialEq)]
pub struct CalledMarker {
    /// FTDNA marker name (normalized — e.g. `DYS19`, not `DYS19/DYS394`).
    pub marker: String,
    /// FTDNA-convention value (caller count + calibrated offset), or the raw count when not calibrated.
    pub value: i32,
    pub status: MarkerStatus,
    /// Enclosing-read depth behind the call.
    pub depth: u32,
}

/// Normalize a caller locus name (from the HipSTR BED) to its base FTDNA marker: take the first of a
/// `/`-alias (`DYS19/DYS394` → `DYS19`), drop a `_N` copy suffix (`DYS385_1` → `DYS385`) and a `.N`
/// partial suffix (`DYS389II.1` → `DYS389II`).
pub fn normalize_marker(caller_name: &str) -> String {
    let n = caller_name.split('/').next().unwrap_or(caller_name);
    let n = n.split('_').next().unwrap_or(n);
    n.split('.').next().unwrap_or(n).to_string()
}

fn offset(marker: &str) -> Option<i32> {
    OFFSETS.iter().find(|(m, _)| *m == marker).map(|(_, o)| *o)
}

/// Map one caller locus + its repeat count to the FTDNA convention: the marker name, the
/// convention-adjusted value, and how trustworthy that mapping is.
pub fn to_ftdna(caller_name: &str, caller_copies: i32) -> CalledMarker {
    let marker = normalize_marker(caller_name);
    let (value, status) = if EXCLUDE.contains(&marker.as_str()) {
        (caller_copies, MarkerStatus::Excluded)
    } else if let Some(off) = offset(&marker) {
        (caller_copies + off, if off == 0 { MarkerStatus::Reliable } else { MarkerStatus::ConventionOffset })
    } else {
        (caller_copies, MarkerStatus::Uncalibrated)
    };
    CalledMarker { marker, value, status, depth: 0 }
}

/// Convert the caller's genotypes to FTDNA-convention marker calls: single-copy (one allele),
/// non-low-confidence loci, deduped per marker keeping the deepest. Multi-copy markers (two alleles)
/// are skipped — they need the (excluded) aggregation/nesting conventions.
pub fn called_markers(genotypes: &[StrGenotype]) -> Vec<CalledMarker> {
    use std::collections::HashMap;
    let mut best: HashMap<String, CalledMarker> = HashMap::new();
    for g in genotypes.iter().filter(|g| g.confidence != StrConfidence::Low && g.alleles.len() == 1) {
        let mut cm = to_ftdna(&g.name, g.alleles[0]);
        cm.depth = g.depth;
        best.entry(cm.marker.clone())
            .and_modify(|cur| {
                if cm.depth > cur.depth {
                    *cur = cm.clone();
                }
            })
            .or_insert(cm);
    }
    let mut out: Vec<CalledMarker> = best.into_values().collect();
    out.sort_by(|a, b| a.marker.cmp(&b.marker));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_strips_aliases_copies_and_partials() {
        assert_eq!(normalize_marker("DYS19/DYS394"), "DYS19");
        assert_eq!(normalize_marker("DYS385_1"), "DYS385");
        assert_eq!(normalize_marker("DYS389II.1"), "DYS389II");
        assert_eq!(normalize_marker("DYS393"), "DYS393");
    }

    #[test]
    fn applies_calibrated_conventions() {
        // Reliable (offset 0): caller == FTDNA.
        let r = to_ftdna("DYS390", 24);
        assert_eq!((r.value, r.status), (24, MarkerStatus::Reliable));
        // Convention offset: DYS438 caller 10 → FTDNA 12 (+2); DYS19 caller 15 → 14 (-1).
        assert_eq!(to_ftdna("DYS438", 10).value, 12);
        assert_eq!(to_ftdna("DYS438", 10).status, MarkerStatus::ConventionOffset);
        assert_eq!(to_ftdna("DYS19/DYS394", 15).value, 14);
        // Excluded (tract mismatch): raw count, flagged.
        assert_eq!(to_ftdna("Y-GATA-H4", 31).status, MarkerStatus::Excluded);
        // Uncalibrated marker → raw count, flagged.
        assert_eq!(to_ftdna("DYS999", 7).status, MarkerStatus::Uncalibrated);
    }
}
