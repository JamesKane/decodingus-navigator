//! Find out which test made an alignment. This is the Rust port of the Scala `TestType` catalog
//! and its `TestTypeInference`.
//!
//! The header probe ([`crate::probe`]) knows the *platform* alone: PacBio gives HiFi, Illumina
//! gives WGS, and so on. It can not separate a **targeted** test from a whole-genome one, because
//! the two look the same in a SAM header. An FTDNA Big Y, a Full Genomes Y Elite, a YSEQ test and
//! an mtFull run are all targeted.
//!
//! The Scala app separated them by the **shape of the coverage**. A Big Y BAM piles its reads on
//! chrY, and its autosomes are almost empty. An mtFull run piles them on chrM.
//!
//! This module gets that shape at low cost, from the **BAI index**. That index holds the count of
//! mapped records of each reference, at O(contigs), and it is the same fast path that
//! [`crate::sex`] uses. The code normalizes those counts to a coverage proxy. It then puts that
//! together with the platform, and with a vendor hint when there is one, and takes a test-type
//! code.

use std::path::Path;

use noodles::bam;
use noodles::csi::binning_index::ReferenceSequence as _;

use crate::contig;

// A test-type code is one of the canonical strings of the `navigator_domain::testtype` catalog.
// The display name, the target region and the UI picker all live there. This module decides
// *which* code the coverage shape of a BAM implies, and nothing more. It writes those code
// literals out, and a test checks them against the catalog.

/// The coverage proxy of each chromosome group, as reads × read length ÷ group length. That is the
/// same estimate that the Scala `ChromosomeCoverageStats` used. A group of `None` means that the
/// reference holds no such contig.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageProfile {
    pub autosome_depth: f64,
    pub y_depth: f64,
    pub mt_depth: f64,
    /// Were any autosomal contigs present at all? (A Y-only reference has none → targeted by build.)
    pub has_autosomes: bool,
}

/// The `ASSUMED_READ_LENGTH` of the Scala code. The estimate uses it when nobody knows the true
/// mean read length.
const ASSUMED_READ_LENGTH: u64 = 150;

// The coverage thresholds. The Scala code used absolute cutoffs, at `yCov>1 && autoCov<1`. But a
// real Big Y BAM that somebody aligned to the whole genome carries 1 to 2x of autosomal reads that
// are off target. A real FTDNA Big Y here measured Y at 51x and the autosomes at 1.8x. The
// absolute test reads that as a low-pass WGS run.
//
// This code instead keys a targeted-Y test off the **enrichment ratio of Y against the
// autosomes**. That ratio does not depend on the read length, so it also survives the low coverage
// estimate that a long read gives.
//
// For a targeted-MT test it needs the autosomes to be almost absent. mtDNA is naturally
// high-copy. So a WGS sample shows a very large mt depth, and it is not an mtFull test.
const Y_PRESENT: f64 = 1.0; // Y depth floor below which we don't call targeted-Y at all
const Y_ENRICH: f64 = 5.0; // Y:autosome ratio that marks a Y-targeted capture
const MT_PRESENT: f64 = 10.0;
const AUTOSOME_PRESENT: f64 = 1.0;
const LONG_READ_LEN: u64 = 1000;
const WES_AUTOSOME_DEPTH: f64 = 50.0;
const LOW_PASS_AUTOSOME_DEPTH: f64 = 5.0;

/// Build a [`CoverageProfile`] from the BAI index of a BAM. It scans no read.
///
/// `mean_read_length` makes the estimate better when somebody knows it, for example from
/// `library_stats`. Without it the code uses [`ASSUMED_READ_LENGTH`].
///
/// It returns `None` when the index is absent, or when the code can not read it. A CRAM is one
/// such case, because a `.crai` holds no count of each reference. The caller then keeps the result
/// from the header and the platform.
pub fn coverage_profile_from_bai(bam_path: &Path, mean_read_length: Option<u64>) -> Option<CoverageProfile> {
    let header = crate::reader::read_header(bam_path, None).ok()?;
    let bai_path = bam_path.with_extension("bam.bai");
    let index = bam::bai::fs::read(&bai_path).ok()?;
    let counts: Vec<u64> = index
        .reference_sequences()
        .iter()
        .map(|rs| rs.metadata().map_or(0, |m| m.mapped_record_count()))
        .collect();

    let read_len = mean_read_length.filter(|&l| l > 0).unwrap_or(ASSUMED_READ_LENGTH);
    let (mut a_reads, mut a_len) = (0u64, 0u64);
    let (mut y_reads, mut y_len) = (0u64, 0u64);
    let (mut m_reads, mut m_len) = (0u64, 0u64);
    let mut has_autosomes = false;
    for (i, (name_bytes, map)) in header.reference_sequences().iter().enumerate() {
        let name = String::from_utf8_lossy(name_bytes.as_ref());
        let length = map.length().get() as u64;
        let count = counts.get(i).copied().unwrap_or(0);
        if contig::is_autosome(&name) {
            has_autosomes = true;
            a_reads += count;
            a_len += length;
        } else if contig::is_chr_y(&name) {
            y_reads += count;
            y_len += length;
        } else if contig::is_chr_m(&name) {
            m_reads += count;
            m_len += length;
        }
    }
    let depth = |reads: u64, len: u64| {
        if len > 0 {
            (reads * read_len) as f64 / len as f64
        } else {
            0.0
        }
    };
    Some(CoverageProfile {
        autosome_depth: depth(a_reads, a_len),
        y_depth: depth(y_reads, y_len),
        mt_depth: depth(m_reads, m_len),
        has_autosomes,
    })
}

/// Map a free-text vendor hint to a specific targeted-Y test code (else the honest generic).
fn targeted_y_for_vendor(vendor_hint: Option<&str>) -> &'static str {
    match vendor_hint.map(|v| v.to_lowercase()) {
        // FTDNA sells Big Y alone. But the vendor token does not carry the *generation*, which is
        // 500 or 700. That comes from the `@RG LB` label, which [`crate::probe`] reads and passes
        // in as `big_y_label`. On an older header that leaves that label out, it comes instead
        // from the callable-chrY footprint, which the code resolves after the analysis. So stay
        // generic here, and do not guess a generation from the vendor name alone.
        Some(v) if v.contains("ftdna") || v.contains("familytreedna") => "TARGETED_Y",
        Some(v) if v.contains("full genomes") || v.contains("fullgenomes") => "Y_ELITE",
        Some(v) if v.contains("yseq") => "Y_PRIME",
        // An unknown vendor is not mislabeled to a specific product.
        _ => "TARGETED_Y",
    }
}

/// Pick a WGS subtype code from the platform.
fn wgs_for_platform(platform: Option<&str>, mean_read_length: Option<u64>) -> &'static str {
    if platform.is_some_and(|p| p.to_uppercase().contains("PACBIO"))
        || mean_read_length.is_some_and(|l| l > LONG_READ_LEN)
    {
        "WGS_HIFI"
    } else if platform.is_some_and(|p| {
        let u = p.to_uppercase();
        u.contains("NANOPORE") || u == "ONT"
    }) {
        "WGS_NANOPORE"
    } else {
        "WGS"
    }
}

/// Infer the test type from the coverage shape, the platform and the vendor hint. This is the
/// Scala `inferFromCoverage`.
///
/// With no coverage profile it falls back to the WGS guess from the platform alone. That happens
/// for a CRAM, and for a file with no index. The probe behaved that way before this code. It
/// returns `None` only when the code knows nothing at all.
pub fn infer_test_type(
    profile: Option<&CoverageProfile>,
    platform: Option<&str>,
    vendor_hint: Option<&str>,
    mean_read_length: Option<u64>,
    big_y_label: Option<&str>,
) -> Option<String> {
    // An FTDNA Big Y generation that the header states, in `@RG LB`, is authoritative. It is the
    // own product label of FTDNA, so it wins over the guess from the coverage shape.
    if let Some(code) = big_y_label {
        return Some(code.to_string());
    }
    let Some(p) = profile else {
        // No coverage shape available: platform-only, matching the old probe.
        return platform.map(|_| wgs_for_platform(platform, mean_read_length).to_string());
    };

    // "Autosomal coverage present" = depth above the floor AND autosomal contigs exist at all.
    let has_autosome = p.has_autosomes && p.autosome_depth > AUTOSOME_PRESENT;
    // The enrichment of Y against the autosomes. That is the signal of a targeted-Y test. With no
    // autosome at all, on a Y-only reference, it is infinite.
    let y_ratio = if p.autosome_depth > 0.0 {
        p.y_depth / p.autosome_depth
    } else {
        f64::INFINITY
    };

    // Targeted-Y: Y meaningfully covered AND strongly enriched over the autosomes (off-target
    // autosomal reads are normal), or a Y-only reference (no autosomal contigs).
    let targeted_y = p.y_depth > Y_PRESENT && (!has_autosome || y_ratio > Y_ENRICH);
    // A targeted-MT test covers mtDNA alone. The autosomes are almost absent, and there is no Y.
    // mtDNA is naturally high-copy. So a WGS sample shows a very large mt depth, and it is not an
    // mtFull test.
    let targeted_mt = !targeted_y && p.mt_depth > MT_PRESENT && !has_autosome && p.y_depth <= Y_PRESENT;

    let code = if targeted_y {
        targeted_y_for_vendor(vendor_hint)
    } else if targeted_mt {
        "MT_FULL_SEQUENCE"
    } else if has_autosome && p.y_depth <= Y_PRESENT && p.autosome_depth > WES_AUTOSOME_DEPTH {
        // A very high autosomal depth, with no Y signal. That is an exome capture.
        "WES"
    } else if has_autosome && p.autosome_depth < LOW_PASS_AUTOSOME_DEPTH {
        "WGS_LOW_PASS"
    } else {
        wgs_for_platform(platform, mean_read_length)
    };
    Some(code.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prof(a: f64, y: f64, m: f64, has_auto: bool) -> CoverageProfile {
        CoverageProfile {
            autosome_depth: a,
            y_depth: y,
            mt_depth: m,
            has_autosomes: has_auto,
        }
    }

    #[test]
    fn targeted_y_maps_vendor_or_generic() {
        // A Y-only reference, with no autosome. That is a clean targeted-Y test.
        let p = prof(0.0, 35.0, 0.0, false);
        // An FTDNA test with no generation label stays generic. The vendor token does not say 500
        // or 700. The header `@RG LB`, or the callable-chrY footprint, decides that.
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), Some("FamilyTreeDNA"), None, None).as_deref(),
            Some("TARGETED_Y")
        );
        // An explicit header generation label is authoritative.
        assert_eq!(
            infer_test_type(
                Some(&p),
                Some("ILLUMINA"),
                Some("FamilyTreeDNA"),
                None,
                Some("BIG_Y_500")
            )
            .as_deref(),
            Some("BIG_Y_500")
        );
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), Some("Full Genomes"), None, None).as_deref(),
            Some("Y_ELITE")
        );
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), Some("YSEQ"), None, None).as_deref(),
            Some("Y_PRIME")
        );
        // Unknown vendor → TARGETED_Y, not a mislabeled Big Y.
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), None, None, None).as_deref(),
            Some("TARGETED_Y")
        );
    }

    #[test]
    fn explicit_big_y_label_overrides_coverage_shape() {
        // The header wins. Even a coverage shape that looks like WGS gives the Big Y generation
        // of the header, when `@RG LB` states one.
        let wgs_shape = prof(30.0, 15.0, 1000.0, true);
        assert_eq!(
            infer_test_type(Some(&wgs_shape), Some("ILLUMINA"), None, None, Some("BIG_Y_700")).as_deref(),
            Some("BIG_Y_700")
        );
    }

    #[test]
    fn targeted_y_by_enrichment_with_offtarget_autosomes() {
        // The shape of a real Full Genomes Y Elite run, from B6564_Kane.bam. Y is at 51x, over
        // autosomes at 1.8x that are off target. An absolute test of "autosome<1" would read this
        // as WGS_LOW_PASS. The enrichment of 28x marks it correctly.
        let p = prof(1.84, 51.0, 7.0, true);
        assert_eq!(
            infer_test_type(Some(&p), None, Some("Full Genomes"), None, None).as_deref(),
            Some("Y_ELITE")
        );
        assert_eq!(
            infer_test_type(Some(&p), None, None, None, None).as_deref(),
            Some("TARGETED_Y")
        );
    }

    #[test]
    fn targeted_mt_only_when_autosomes_absent() {
        let p = prof(0.0, 0.0, 800.0, false);
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), None, None, None).as_deref(),
            Some("MT_FULL_SEQUENCE")
        );
    }

    #[test]
    fn wgs_not_mislabeled_by_high_copy_mt_or_modest_y() {
        // Real 30× male WGS (60820188481374.bam): mt 1155× (high-copy), Y 8.5× < autosome 29× →
        // not targeted-MT (Y present + autosomes present), not targeted-Y (ratio 0.29) → WGS.
        let male = prof(29.4, 8.5, 1155.0, true);
        assert_eq!(
            infer_test_type(Some(&male), Some("ILLUMINA"), None, None, None).as_deref(),
            Some("WGS")
        );
        // A female WGS run. Y is near 0, and mt is high. The guard that needs the autosomes to be
        // present keeps this as WGS, and not as targeted-MT.
        let female = prof(30.0, 0.02, 1200.0, true);
        assert_eq!(
            infer_test_type(Some(&female), Some("ILLUMINA"), None, None, None).as_deref(),
            Some("WGS")
        );
    }

    #[test]
    fn wgs_when_autosomes_covered() {
        // Balanced coverage (GFX-like): autosome ≈ Y depth → WGS by platform, not targeted.
        let hifi = prof(8.0, 6.0, 40.0, true);
        assert_eq!(
            infer_test_type(Some(&hifi), Some("PACBIO"), None, None, None).as_deref(),
            Some("WGS_HIFI")
        );
        // Long read by length, no platform string.
        assert_eq!(
            infer_test_type(Some(&hifi), None, None, Some(15000), None).as_deref(),
            Some("WGS_HIFI")
        );
    }

    #[test]
    fn low_pass_and_exome() {
        // Low autosomal depth, no enriched Y → low-pass WGS.
        assert_eq!(
            infer_test_type(Some(&prof(2.0, 0.5, 4.0, true)), Some("ILLUMINA"), None, None, None).as_deref(),
            Some("WGS_LOW_PASS")
        );
        // High autosomal, no Y/MT contigs → exome.
        assert_eq!(
            infer_test_type(Some(&prof(80.0, 0.0, 0.0, true)), Some("ILLUMINA"), None, None, None).as_deref(),
            Some("WES")
        );
    }

    #[test]
    fn every_emitted_code_is_in_the_domain_catalog() {
        // The canonical catalog must know every code that this module writes. Else the UI picker,
        // and display_name, would show a raw code. This test covers the output of every branch.
        let shapes = [
            (prof(0.0, 35.0, 0.0, false), Some("FamilyTreeDNA")),
            (prof(0.0, 35.0, 0.0, false), Some("Full Genomes")),
            (prof(0.0, 35.0, 0.0, false), Some("YSEQ")),
            (prof(0.0, 35.0, 0.0, false), None),
            (prof(0.0, 0.0, 800.0, false), None),
            (prof(29.4, 8.5, 1155.0, true), None),
            (prof(2.0, 0.5, 4.0, true), None),
            (prof(80.0, 0.0, 0.0, true), None),
        ];
        for (p, vendor) in shapes {
            let code = infer_test_type(Some(&p), Some("ILLUMINA"), vendor, None, None).unwrap();
            assert!(
                navigator_domain::testtype::by_code(&code).is_some(),
                "code {code} not in catalog"
            );
        }
        for plat in ["PACBIO", "ILLUMINA", "NANOPORE"] {
            let code = infer_test_type(None, Some(plat), None, None, None).unwrap();
            assert!(
                navigator_domain::testtype::by_code(&code).is_some(),
                "code {code} not in catalog"
            );
        }
    }

    #[test]
    fn no_profile_falls_back_to_platform() {
        assert_eq!(
            infer_test_type(None, Some("PACBIO"), None, None, None).as_deref(),
            Some("WGS_HIFI")
        );
        assert_eq!(
            infer_test_type(None, Some("ILLUMINA"), None, None, None).as_deref(),
            Some("WGS")
        );
        assert_eq!(infer_test_type(None, None, None, None, None), None);
    }
}
