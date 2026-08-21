//! The layer that holds the FTDNA Y-STR convention. It maps the repeat count that the caller gives
//! at a HipSTR locus to the value that an FTDNA marker holds.
//!
//! The HipSTR reference already names the loci, as DYS393, DYS19/DYS394 and so on. See
//! [`crate::strref`]. So the caller emits a DYS name directly. What differs is the **convention
//! for the count**.
//!
//! FTDNA reports a value at each marker. That value is the repeat count of the caller, plus a
//! fixed offset. The offset is 0 at most markers, and ±1 to 3 at some. There is also a set of
//! markers whose HipSTR tract does not match the FTDNA marker 1:1. Those are the large tract
//! mismatches, and the markers with more than one copy or with a nest inside, such as DYS385,
//! DYS464 and DYS389II. No single offset can map one of those.
//!
//! The offset table below comes from a **calibration against a corpus of 216 Big Y kits**. That
//! corpus is the FTDNA R1b project, with a `chrYM.cram` realigned to CHM13, and the FTDNA DYS CSV
//! of each kit. At each marker the offset is the modal `ftdna − caller` difference across the
//! kits. The table keeps it only where the corpus agrees, at 70% or more, and over 20 kits or
//! more. The harness is `examples/str_calibrate.rs`.
//!
//! The calibration also checked that the offsets do not depend on the build. Where the CHM13
//! corpus and the earlier corpus of 14 GRCh38 kits overlap, the offsets match: DYS438 +2, DYS435
//! +2, DYS474 −3, DYS442 −3, DYS520 −2, DYS585 −3, DYS615 −2, DYS629 −3.
//!
//! The CHM13 HipSTR liftover dropped some markers: DYS19, DYS391, DYS426, DYS445, DYS461, DYS512,
//! DYS549, DYS565, DYS567, DYS578, DYS589, DYS632 and more. Because an offset does not depend on
//! the build, those markers keep the values of the GRCh38 corpus. They still serve the path that
//! calls over a GRCh38 BAM. They will serve the CHM13 path too, once somebody recovers those loci
//! in the lifted reference.
//!
//! A marker that this table does not hold reports `Uncalibrated`. Run the harness over more kits
//! to extend the table.

use crate::strcaller::{StrConfidence, StrGenotype};

/// FTDNA value = caller repeat count + offset. Covers the calibrated single-copy markers (offset 0 =
/// "reliable", ±1–3 = a real convention). Markers absent here are either in [`EXCLUDE`] or uncalibrated.
static OFFSETS: &[(&str, i32)] = &[
    // These are reliable, at an offset of 0. Most come from the 216-kit CHM13 corpus, at 100%
    // agreement. A few come from the GRCh38 corpus, for loci that the CHM13 lift dropped, or
    // sampled too little: DYS388, DYS426, DYS445, DYS487, DYS494, DYS505, DYS549, DYS556, DYS565,
    // DYS567, DYS577 and DYS578.
    ("DYS388", 0),
    ("DYS390", 0),
    ("DYS392", 0),
    ("DYS426", 0),
    ("DYS434", 0),
    ("DYS436", 0),
    ("DYS445", 0),
    ("DYS446", 0),
    ("DYS453", 0),
    ("DYS454", 0),
    ("DYS455", 0),
    ("DYS458", 0),
    ("DYS462", 0),
    ("DYS472", 0),
    ("DYS476", 0),
    ("DYS477", 0),
    ("DYS480", 0),
    ("DYS487", 0),
    ("DYS488", 0),
    ("DYS490", 0),
    ("DYS492", 0),
    ("DYS494", 0),
    ("DYS497", 0),
    ("DYS499", 0),
    ("DYS505", 0),
    ("DYS508", 0),
    ("DYS530", 0),
    ("DYS531", 0),
    ("DYS533", 0),
    ("DYS549", 0),
    ("DYS556", 0),
    ("DYS561", 0),
    ("DYS565", 0),
    ("DYS567", 0),
    ("DYS568", 0),
    ("DYS569", 0),
    ("DYS573", 0),
    ("DYS574", 0),
    ("DYS575", 0),
    ("DYS577", 0),
    ("DYS578", 0),
    ("DYS579", 0),
    ("DYS580", 0),
    ("DYS581", 0),
    ("DYS583", 0),
    ("DYS584", 0),
    ("DYS590", 0),
    ("DYS593", 0),
    ("DYS594", 0),
    ("DYS595", 0),
    ("DYS618", 0),
    ("DYS620", 0),
    ("DYS621", 0),
    ("DYS635", 0),
    ("DYS638", 0),
    ("DYS640", 0),
    ("DYS641", 0),
    ("DYS645", 0),
    ("DYS714", 0),
    ("Y-GATA-A10", 0),
    // The offsets of the convention, at ±1 to 3. DYS19, DYS391, DYS461, DYS512, DYS589 and DYS632
    // come from the GRCh38 corpus, because the CHM13 lift dropped them. The rest come from the
    // 216-kit CHM13 corpus. The table once left DYS460 out. The larger corpus resolves it to a
    // clean +1, at n=180 and 98%.
    ("DYS19", -1),
    ("DYS389I", 1),
    ("DYS391", -1),
    ("DYS425", 2),
    ("DYS435", 2),
    ("DYS438", 2),
    ("DYS442", -3),
    ("DYS456", 1),
    ("DYS460", 1),
    ("DYS461", 1),
    ("DYS463", 2),
    ("DYS474", -3),
    ("DYS485", -1),
    ("DYS512", -3),
    ("DYS520", -2),
    ("DYS522", 1),
    ("DYS525", 1),
    ("DYS537", 1),
    ("DYS538", 2),
    ("DYS539", 1),
    ("DYS559", 1),
    ("DYS572", 1),
    ("DYS585", -3),
    ("DYS587", 1),
    ("DYS589", -1),
    ("DYS615", -2),
    ("DYS629", -3),
    ("DYS632", -2),
    ("DYS642", 1),
];

/// The markers whose HipSTR tract no single offset can map to the FTDNA value. Those are the large
/// tract mismatches, and the markers with more than one copy or with a nest inside. DYS385 and
/// DYS464 split into sub-loci, and DYS389II holds a nest.
///
/// They report as `Excluded`. The enclosing-read caller does not yet give a value here that you
/// can compare against a vendor.
static EXCLUDE: &[&str] = &[
    // Multi-copy / nested (split sub-loci, never a single vendor-comparable value).
    "DYS385",
    "DYS389II",
    "DYS459",
    "DYS464",
    "YCAII",
    "CDY",
    // Large tract mismatch or variable across the 216-kit corpus (<70% offset agreement).
    "DYF406S1",
    "DYS393",
    "DYS448",
    "DYS449",
    "DYS450",
    "DYS470",
    "DYS475",
    "DYS481",
    "DYS484",
    "DYS495",
    "DYS502",
    "DYS504",
    "DYS510",
    "DYS511",
    "DYS513",
    "DYS516",
    "DYS532",
    "DYS534",
    "DYS540",
    "DYS541",
    "DYS543",
    "DYS544",
    "DYS551",
    "DYS552",
    "DYS557",
    "DYS570",
    "DYS576",
    "DYS588",
    "DYS607",
    "DYS616",
    "DYS624",
    "DYS631",
    "DYS637",
    "DYS717",
    "Y-GATA-H4",
];

/// A few offsets **do depend on the build**. The CHM13 HipSTR liftover moved these tract
/// boundaries by one repeat unit, so the enclosing-read count differs by 1 between CHM13 and
/// GRCh38. [`OFFSETS`] holds the CHM13 value, which is the primary corpus, and the code *adds*
/// this delta on the GRCh38 path.
///
/// Two independent GRCh38 corpora confirmed this against the 216-kit CHM13 corpus: a 14-kit set,
/// and kit 27520. DYS389I, DYS456, DYS525, DYS537 and DYS539 are +1 on CHM13 and 0 on GRCh38.
/// DYS714 is the reverse: 0 on CHM13 and +1 on GRCh38.
static GRCH38_DELTA: &[(&str, i32)] = &[
    ("DYS389I", -1),
    ("DYS456", -1),
    ("DYS525", -1),
    ("DYS537", -1),
    ("DYS539", -1),
    ("DYS714", 1),
];

/// The reference build that the repeat counts of the caller came from. It selects the convention
/// offset of the markers that depend on the build. See [`GRCH38_DELTA`]. The default,
/// [`StrBuild::Chm13`], is the primary corpus.
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
    /// The tract does not match, or the marker holds more than one copy, or it holds a nest.
    /// There is no 1:1 mapping, and the value is the raw count of the caller.
    Excluded,
    /// The calibration table does not hold this marker. The value is the raw count of the caller,
    /// and it waits for a calibration.
    Uncalibrated,
}

/// One marker called from sequence, expressed in the FTDNA convention.
#[derive(Debug, Clone, PartialEq)]
pub struct CalledMarker {
    /// The name of the FTDNA marker, in its normal form. That is `DYS19`, and not
    /// `DYS19/DYS394`.
    pub marker: String,
    /// FTDNA-convention value (caller count + calibrated offset), or the raw count when not calibrated.
    pub value: i32,
    pub status: MarkerStatus,
    /// Enclosing-read depth behind the call.
    pub depth: u32,
}

/// Bring the locus name of the caller, which comes from the HipSTR BED, to its base FTDNA marker.
/// It takes the first name of a `/` alias, so `DYS19/DYS394` gives `DYS19`. It drops a `_N` copy
/// suffix, so `DYS385_1` gives `DYS385`. And it drops a `.N` partial suffix, so `DYS389II.1` gives
/// `DYS389II`.
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

/// Map one locus of the caller, and its repeat count, to the FTDNA convention of one build. It
/// gives the marker name, the value after the convention offset, and how much you can trust that
/// mapping. The build selects the offset of the markers that depend on it. See
/// [`GRCH38_DELTA`].
pub fn to_ftdna_build(caller_name: &str, caller_copies: i32, build: StrBuild) -> CalledMarker {
    let marker = normalize_marker(caller_name);
    let (value, status) = if EXCLUDE.contains(&marker.as_str()) {
        (caller_copies, MarkerStatus::Excluded)
    } else if let Some(base) = offset(&marker) {
        let off = base
            + if build == StrBuild::Grch38 {
                grch38_delta(&marker)
            } else {
                0
            };
        (
            caller_copies + off,
            if off == 0 {
                MarkerStatus::Reliable
            } else {
                MarkerStatus::ConventionOffset
            },
        )
    } else {
        (caller_copies, MarkerStatus::Uncalibrated)
    };
    CalledMarker {
        marker,
        value,
        status,
        depth: 0,
    }
}

/// Convert the caller's genotypes to FTDNA-convention marker calls against the CHM13 corpus.
/// Use [`called_markers_build`] when the caller ran on GRCh38.
pub fn called_markers(genotypes: &[StrGenotype]) -> Vec<CalledMarker> {
    called_markers_build(genotypes, StrBuild::Chm13)
}

/// Turn the genotypes of the caller into marker calls in the FTDNA convention, for one build.
///
/// It takes a locus with one copy, which holds one allele, and whose confidence is not low. Where
/// two loci give the same marker, it keeps the one with the deepest coverage.
///
/// It skips a marker with more than one copy, which holds two alleles. Such a marker needs the
/// conventions for aggregation and for a nest, and this module excludes those.
pub fn called_markers_build(genotypes: &[StrGenotype], build: StrBuild) -> Vec<CalledMarker> {
    use std::collections::HashMap;
    let mut best: HashMap<String, CalledMarker> = HashMap::new();
    for g in genotypes
        .iter()
        .filter(|g| g.confidence != StrConfidence::Low && g.alleles.len() == 1)
    {
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
        // A marker that does not depend on the build does not change with the build.
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
