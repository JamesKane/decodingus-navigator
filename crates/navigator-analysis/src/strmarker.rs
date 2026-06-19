//! FTDNA Y-STR convention layer — map the caller's HipSTR-locus repeat counts to FTDNA marker values.
//!
//! The HipSTR reference already names the loci (DYS393, DYS19/DYS394, …; see [`crate::strref`]), so
//! the caller emits DYS names directly. What differs is the **counting convention**: FTDNA reports a
//! per-marker value that is the caller's repeat count plus a fixed offset (0 for most, ±1–3 for some),
//! and a set of markers whose HipSTR tract doesn't correspond 1:1 to the FTDNA marker (large tract
//! mismatches, plus multi-copy/nested markers like DYS385/DYS464/DYS389II) which can't be mapped by a
//! single offset.
//!
//! The offset table below was **calibrated against a 216-kit Big Y corpus** (FTDNA R1b project,
//! CHM13-realigned `chrYM.cram` + each kit's FTDNA DYS CSV): per marker, the offset is the modal
//! `ftdna − caller` difference across kits, kept only where it agrees across the corpus (≥70%, and
//! ≥20 kits to be authoritative). The harness is `examples/str_calibrate.rs`. The calibration cross-
//! validated build-independence: where the CHM13 corpus and the earlier 14-kit GRCh38 corpus overlap,
//! the offsets match (DYS438 +2, DYS435 +2, DYS474 −3, DYS442 −3, DYS520 −2, DYS585 −3, DYS615 −2,
//! DYS629 −3). Because offsets are build-independent, markers the CHM13 HipSTR liftover dropped
//! (DYS19, DYS391, DYS426, DYS445, DYS461, DYS512, DYS549, DYS565, DYS567, DYS578, DYS589, DYS632 …)
//! retain their GRCh38-corpus values — they still serve the GRCh38/BAM calling path, and will serve
//! the CHM13 path once those loci are recovered in the lifted reference. Markers absent here report
//! `Uncalibrated`. Re-run the harness over more kits to extend the table.

use crate::strcaller::{StrConfidence, StrGenotype};

/// FTDNA value = caller repeat count + offset. Covers the calibrated single-copy markers (offset 0 =
/// "reliable", ±1–3 = a real convention). Markers absent here are either in [`EXCLUDE`] or uncalibrated.
static OFFSETS: &[(&str, i32)] = &[
    // Reliable (offset 0). Most are CHM13 216-kit ≥100% agreement; a few (DYS388, DYS426, DYS445,
    // DYS487, DYS494, DYS505, DYS549, DYS556, DYS565, DYS567, DYS577, DYS578) are GRCh38-corpus
    // retentions for loci the CHM13 lift dropped or under-sampled.
    ("DYS388", 0), ("DYS390", 0), ("DYS392", 0), ("DYS426", 0), ("DYS434", 0), ("DYS436", 0),
    ("DYS445", 0), ("DYS446", 0), ("DYS453", 0), ("DYS454", 0), ("DYS455", 0), ("DYS458", 0),
    ("DYS462", 0), ("DYS472", 0), ("DYS476", 0), ("DYS477", 0), ("DYS480", 0), ("DYS487", 0),
    ("DYS488", 0), ("DYS490", 0), ("DYS492", 0), ("DYS494", 0), ("DYS497", 0), ("DYS499", 0),
    ("DYS505", 0), ("DYS508", 0), ("DYS530", 0), ("DYS531", 0), ("DYS533", 0), ("DYS549", 0),
    ("DYS556", 0), ("DYS561", 0), ("DYS565", 0), ("DYS567", 0), ("DYS568", 0), ("DYS569", 0),
    ("DYS573", 0), ("DYS574", 0), ("DYS575", 0), ("DYS577", 0), ("DYS578", 0), ("DYS579", 0),
    ("DYS580", 0), ("DYS581", 0), ("DYS583", 0), ("DYS584", 0), ("DYS590", 0), ("DYS593", 0),
    ("DYS594", 0), ("DYS595", 0), ("DYS618", 0), ("DYS620", 0), ("DYS621", 0), ("DYS635", 0),
    ("DYS638", 0), ("DYS640", 0), ("DYS641", 0), ("DYS645", 0), ("DYS714", 0), ("Y-GATA-A10", 0),
    // Convention offsets (±1–3). DYS19/DYS391/DYS461/DYS512/DYS589/DYS632 are GRCh38-corpus
    // retentions (dropped by the CHM13 lift); the rest are CHM13 216-kit. DYS460 was previously
    // excluded but the larger corpus resolves it to a clean +1 (n=180, 98%).
    ("DYS19", -1), ("DYS389I", 1), ("DYS391", -1), ("DYS425", 2), ("DYS435", 2), ("DYS438", 2),
    ("DYS442", -3), ("DYS456", 1), ("DYS460", 1), ("DYS461", 1), ("DYS463", 2), ("DYS474", -3),
    ("DYS485", -1), ("DYS512", -3), ("DYS520", -2), ("DYS522", 1), ("DYS525", 1), ("DYS537", 1),
    ("DYS538", 2), ("DYS539", 1), ("DYS559", 1), ("DYS572", 1), ("DYS585", -3), ("DYS587", 1),
    ("DYS589", -1), ("DYS615", -2), ("DYS629", -3), ("DYS632", -2), ("DYS642", 1),
];

/// Markers whose HipSTR tract can't be mapped to the FTDNA value by a single offset: large tract
/// mismatches and multi-copy/nested markers (DYS385/DYS464 split sub-loci, DYS389II nesting). Reported
/// as `Excluded` — the enclosing-read caller doesn't yield a vendor-comparable value here (yet).
static EXCLUDE: &[&str] = &[
    // Multi-copy / nested (split sub-loci, never a single vendor-comparable value).
    "DYS385", "DYS389II", "DYS459", "DYS464", "YCAII", "CDY",
    // Large tract mismatch or variable across the 216-kit corpus (<70% offset agreement).
    "DYF406S1", "DYS393", "DYS448", "DYS449", "DYS450", "DYS470", "DYS475", "DYS481", "DYS484",
    "DYS495", "DYS502", "DYS504", "DYS510", "DYS511", "DYS513", "DYS516", "DYS532", "DYS534",
    "DYS540", "DYS541", "DYS543", "DYS544", "DYS551", "DYS552", "DYS557", "DYS570", "DYS576",
    "DYS588", "DYS607", "DYS616", "DYS624", "DYS631", "DYS637", "DYS717", "Y-GATA-H4",
];

/// A handful of offsets are **build-dependent**: the CHM13 HipSTR liftover shifted these tract
/// boundaries by one repeat unit, so the enclosing-read count differs by 1 between CHM13 and GRCh38.
/// [`OFFSETS`] holds the CHM13 value (the primary corpus); this delta is *added* on the GRCh38 path.
/// Verified by two independent GRCh38 corpora (a 14-kit set + kit 27520) agreeing against the 216-kit
/// CHM13 corpus: DYS389I/DYS456/DYS525/DYS537/DYS539 are +1 on CHM13 but 0 on GRCh38; DYS714 is the
/// reverse (0 on CHM13, +1 on GRCh38).
static GRCH38_DELTA: &[(&str, i32)] = &[
    ("DYS389I", -1), ("DYS456", -1), ("DYS525", -1), ("DYS537", -1), ("DYS539", -1), ("DYS714", 1),
];

/// Which reference build the caller's repeat counts came from — selects the convention offset for the
/// build-dependent markers (see [`GRCH38_DELTA`]). Default ([`StrBuild::Chm13`]) is the primary corpus.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StrBuild {
    #[default]
    Chm13,
    Grch38,
}

impl StrBuild {
    /// Classify a stored `reference_build` string. Anything that looks like GRCh38/hg38/b38 is the
    /// GRCh38 path; everything else (CHM13/T2T, and the default) uses the CHM13-calibrated offsets.
    pub fn from_build_str(build: &str) -> Self {
        let b = build.to_ascii_lowercase();
        if b.contains("38") || b.contains("hg38") {
            StrBuild::Grch38
        } else {
            StrBuild::Chm13
        }
    }
}

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

fn grch38_delta(marker: &str) -> i32 {
    GRCH38_DELTA.iter().find(|(m, _)| *m == marker).map_or(0, |(_, d)| *d)
}

/// Map one caller locus + its repeat count to the FTDNA convention against the CHM13 corpus (the
/// primary path). Use [`to_ftdna_build`] when the caller ran on GRCh38.
pub fn to_ftdna(caller_name: &str, caller_copies: i32) -> CalledMarker {
    to_ftdna_build(caller_name, caller_copies, StrBuild::Chm13)
}

/// Map one caller locus + its repeat count to the FTDNA convention for a specific build: the marker
/// name, the convention-adjusted value, and how trustworthy that mapping is. The build selects the
/// offset for the build-dependent markers (see [`GRCH38_DELTA`]).
pub fn to_ftdna_build(caller_name: &str, caller_copies: i32, build: StrBuild) -> CalledMarker {
    let marker = normalize_marker(caller_name);
    let (value, status) = if EXCLUDE.contains(&marker.as_str()) {
        (caller_copies, MarkerStatus::Excluded)
    } else if let Some(base) = offset(&marker) {
        let off = base + if build == StrBuild::Grch38 { grch38_delta(&marker) } else { 0 };
        (caller_copies + off, if off == 0 { MarkerStatus::Reliable } else { MarkerStatus::ConventionOffset })
    } else {
        (caller_copies, MarkerStatus::Uncalibrated)
    };
    CalledMarker { marker, value, status, depth: 0 }
}

/// Convert the caller's genotypes to FTDNA-convention marker calls against the CHM13 corpus.
/// Use [`called_markers_build`] when the caller ran on GRCh38.
pub fn called_markers(genotypes: &[StrGenotype]) -> Vec<CalledMarker> {
    called_markers_build(genotypes, StrBuild::Chm13)
}

/// Convert the caller's genotypes to FTDNA-convention marker calls for a specific build: single-copy
/// (one allele), non-low-confidence loci, deduped per marker keeping the deepest. Multi-copy markers
/// (two alleles) are skipped — they need the (excluded) aggregation/nesting conventions.
pub fn called_markers_build(genotypes: &[StrGenotype], build: StrBuild) -> Vec<CalledMarker> {
    use std::collections::HashMap;
    let mut best: HashMap<String, CalledMarker> = HashMap::new();
    for g in genotypes.iter().filter(|g| g.confidence != StrConfidence::Low && g.alleles.len() == 1) {
        let mut cm = to_ftdna_build(&g.name, g.alleles[0], build);
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

    #[test]
    fn build_dependent_offsets_differ() {
        // DYS389I: +1 on CHM13 (the default), 0 on GRCh38.
        assert_eq!(to_ftdna("DYS389I", 13).value, 14);
        assert_eq!(to_ftdna_build("DYS389I", 13, StrBuild::Chm13).value, 14);
        let g = to_ftdna_build("DYS389I", 13, StrBuild::Grch38);
        assert_eq!((g.value, g.status), (13, MarkerStatus::Reliable));
        // DYS714: 0 on CHM13, +1 on GRCh38.
        assert_eq!(to_ftdna("DYS714", 24).value, 24);
        let h = to_ftdna_build("DYS714", 24, StrBuild::Grch38);
        assert_eq!((h.value, h.status), (25, MarkerStatus::ConventionOffset));
        // Build-independent marker is unaffected by build.
        assert_eq!(to_ftdna_build("DYS438", 10, StrBuild::Grch38).value, 12);
    }

    #[test]
    fn build_classifier() {
        assert_eq!(StrBuild::from_build_str("GRCh38"), StrBuild::Grch38);
        assert_eq!(StrBuild::from_build_str("hg38"), StrBuild::Grch38);
        assert_eq!(StrBuild::from_build_str("chm13v2.0"), StrBuild::Chm13);
        assert_eq!(StrBuild::from_build_str("CHM13v2MaskedRcrs"), StrBuild::Chm13);
    }
}
