//! Test-type identification — Rust port of the Scala `TestType` catalog + `TestTypeInference`.
//!
//! The header probe ([`crate::probe`]) only knows *platform* (PacBio→HiFi, Illumina→WGS, …). It
//! can't tell a **targeted** test (FTDNA Big Y, Full Genomes Y Elite, YSEQ, an mtFull run) from a
//! whole-genome one — those look the same in the SAM header. The Scala app distinguished them by
//! **coverage shape**: a Big Y BAM has reads piled on chrY with the autosomes near-empty; an mtFull
//! run piles on chrM. We reproduce that cheaply from the **BAI index** (per-reference mapped-record
//! counts, O(contigs) — the same fast path [`crate::sex`] uses), normalized to a coverage proxy, and
//! combine it with the platform + an optional vendor hint to pick a test-type code.

use std::path::Path;

use noodles::bam;
use noodles::csi::binning_index::ReferenceSequence as _;

use crate::contig;

// Test-type codes are the canonical `navigator_domain::testtype` catalog strings — display names,
// target region, and the UI picker live there. This module only decides *which* code a BAM's
// coverage shape implies, emitting those code literals (validated against the catalog by a test).

/// Per-chromosome-group coverage proxies (reads × read-length ÷ group length), the same estimate
/// the Scala `ChromosomeCoverageStats` used. `None` group ⇒ no such contig in the reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageProfile {
    pub autosome_depth: f64,
    pub y_depth: f64,
    pub mt_depth: f64,
    /// Were any autosomal contigs present at all? (A Y-only reference has none → targeted by build.)
    pub has_autosomes: bool,
}

/// Scala `ASSUMED_READ_LENGTH` — coverage estimate when the true mean read length is unknown.
const ASSUMED_READ_LENGTH: u64 = 150;

// Coverage thresholds. The Scala used absolute cutoffs (`yCov>1 && autoCov<1`), but real Big Y
// BAMs aligned to the whole genome carry ~1-2× off-target autosomal reads (a real FTDNA Big Y here
// measured Y 51× / autosome 1.8×), which the absolute test mislabels as low-pass WGS. We instead key
// targeted-Y off the **Y:autosome enrichment ratio** — read-length-independent, so it also survives
// the long-read coverage underestimate — and require autosomes essentially absent for targeted-MT
// (mtDNA is naturally high-copy, so a WGS sample shows huge mt depth without being an mtFull test).
const Y_PRESENT: f64 = 1.0; // Y depth floor below which we don't call targeted-Y at all
const Y_ENRICH: f64 = 5.0; // Y:autosome ratio that marks a Y-targeted capture
const MT_PRESENT: f64 = 10.0;
const AUTOSOME_PRESENT: f64 = 1.0;
const LONG_READ_LEN: u64 = 1000;
const WES_AUTOSOME_DEPTH: f64 = 50.0;
const LOW_PASS_AUTOSOME_DEPTH: f64 = 5.0;

/// Build a [`CoverageProfile`] from a BAM's BAI index (no read scan). `mean_read_length` refines
/// the estimate when known (e.g. from `library_stats`); else [`ASSUMED_READ_LENGTH`]. Returns `None`
/// when the index is absent/unreadable (e.g. CRAM — `.crai` carries no per-reference counts), so the
/// caller keeps the header/platform result.
pub fn coverage_profile_from_bai(bam_path: &Path, mean_read_length: Option<u64>) -> Option<CoverageProfile> {
    let header = crate::reader::read_header(bam_path, None).ok()?;
    let bai_path = bam_path.with_extension("bam.bai");
    let index = bam::bai::read(&bai_path).ok()?;
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
        Some(v) if v.contains("ftdna") || v.contains("familytreedna") => "BIG_Y_700",
        Some(v) if v.contains("full genomes") || v.contains("fullgenomes") => "Y_ELITE",
        Some(v) if v.contains("yseq") => "Y_PRIME",
        // Scala defaulted to BIG_Y_700; we return TARGETED_Y so an unknown vendor isn't mislabeled.
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

/// Infer the test type from coverage shape + platform + vendor hint (Scala `inferFromCoverage`).
///
/// With no coverage profile (CRAM / unindexed), falls back to the platform-only WGS guess — the
/// pre-existing probe behavior. Returns `None` only when nothing at all is known.
pub fn infer_test_type(
    profile: Option<&CoverageProfile>,
    platform: Option<&str>,
    vendor_hint: Option<&str>,
    mean_read_length: Option<u64>,
) -> Option<String> {
    let Some(p) = profile else {
        // No coverage shape available: platform-only, matching the old probe.
        return platform.map(|_| wgs_for_platform(platform, mean_read_length).to_string());
    };

    // "Autosomal coverage present" = depth above the floor AND autosomal contigs exist at all.
    let has_autosome = p.has_autosomes && p.autosome_depth > AUTOSOME_PRESENT;
    // Y:autosome enrichment — the targeted-Y signal. No autosomes (Y-only reference) ⇒ infinite.
    let y_ratio = if p.autosome_depth > 0.0 {
        p.y_depth / p.autosome_depth
    } else {
        f64::INFINITY
    };

    // Targeted-Y: Y meaningfully covered AND strongly enriched over the autosomes (off-target
    // autosomal reads are normal), or a Y-only reference (no autosomal contigs).
    let targeted_y = p.y_depth > Y_PRESENT && (!has_autosome || y_ratio > Y_ENRICH);
    // Targeted-MT: only mtDNA covered — autosomes essentially absent and no Y (mtDNA is naturally
    // high-copy, so a WGS sample shows huge mt depth without being an mtFull test).
    let targeted_mt = !targeted_y && p.mt_depth > MT_PRESENT && !has_autosome && p.y_depth <= Y_PRESENT;

    let code = if targeted_y {
        targeted_y_for_vendor(vendor_hint)
    } else if targeted_mt {
        "MT_FULL_SEQUENCE"
    } else if has_autosome && p.y_depth <= Y_PRESENT && p.autosome_depth > WES_AUTOSOME_DEPTH {
        // Very high autosomal depth with no Y signal — exome capture.
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
        // Y-only reference (no autosomes) — clean targeted-Y.
        let p = prof(0.0, 35.0, 0.0, false);
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), Some("FamilyTreeDNA"), None).as_deref(),
            Some("BIG_Y_700")
        );
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), Some("Full Genomes"), None).as_deref(),
            Some("Y_ELITE")
        );
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), Some("YSEQ"), None).as_deref(),
            Some("Y_PRIME")
        );
        // Unknown vendor → TARGETED_Y, not a mislabeled Big Y.
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), None, None).as_deref(),
            Some("TARGETED_Y")
        );
    }

    #[test]
    fn targeted_y_by_enrichment_with_offtarget_autosomes() {
        // Real Full Genomes Y Elite shape (B6564_Kane.bam): Y 51× over autosome 1.8× off-target —
        // absolute "autosome<1" would mislabel this WGS_LOW_PASS; the 28× enrichment marks it.
        let p = prof(1.84, 51.0, 7.0, true);
        assert_eq!(
            infer_test_type(Some(&p), None, Some("Full Genomes"), None).as_deref(),
            Some("Y_ELITE")
        );
        assert_eq!(
            infer_test_type(Some(&p), None, None, None).as_deref(),
            Some("TARGETED_Y")
        );
    }

    #[test]
    fn targeted_mt_only_when_autosomes_absent() {
        let p = prof(0.0, 0.0, 800.0, false);
        assert_eq!(
            infer_test_type(Some(&p), Some("ILLUMINA"), None, None).as_deref(),
            Some("MT_FULL_SEQUENCE")
        );
    }

    #[test]
    fn wgs_not_mislabeled_by_high_copy_mt_or_modest_y() {
        // Real 30× male WGS (60820188481374.bam): mt 1155× (high-copy), Y 8.5× < autosome 29× →
        // not targeted-MT (Y present + autosomes present), not targeted-Y (ratio 0.29) → WGS.
        let male = prof(29.4, 8.5, 1155.0, true);
        assert_eq!(
            infer_test_type(Some(&male), Some("ILLUMINA"), None, None).as_deref(),
            Some("WGS")
        );
        // Female WGS: y≈0, mt high — the autosome-present guard keeps it WGS, not targeted-MT.
        let female = prof(30.0, 0.02, 1200.0, true);
        assert_eq!(
            infer_test_type(Some(&female), Some("ILLUMINA"), None, None).as_deref(),
            Some("WGS")
        );
    }

    #[test]
    fn wgs_when_autosomes_covered() {
        // Balanced coverage (GFX-like): autosome ≈ Y depth → WGS by platform, not targeted.
        let hifi = prof(8.0, 6.0, 40.0, true);
        assert_eq!(
            infer_test_type(Some(&hifi), Some("PACBIO"), None, None).as_deref(),
            Some("WGS_HIFI")
        );
        // Long read by length, no platform string.
        assert_eq!(
            infer_test_type(Some(&hifi), None, None, Some(15000)).as_deref(),
            Some("WGS_HIFI")
        );
    }

    #[test]
    fn low_pass_and_exome() {
        // Low autosomal depth, no enriched Y → low-pass WGS.
        assert_eq!(
            infer_test_type(Some(&prof(2.0, 0.5, 4.0, true)), Some("ILLUMINA"), None, None).as_deref(),
            Some("WGS_LOW_PASS")
        );
        // High autosomal, no Y/MT contigs → exome.
        assert_eq!(
            infer_test_type(Some(&prof(80.0, 0.0, 0.0, true)), Some("ILLUMINA"), None, None).as_deref(),
            Some("WES")
        );
    }

    #[test]
    fn every_emitted_code_is_in_the_domain_catalog() {
        // The codes we emit must be recognized by the canonical catalog (else the UI picker /
        // display_name would show a raw code). Exercise every branch's output.
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
            let code = infer_test_type(Some(&p), Some("ILLUMINA"), vendor, None).unwrap();
            assert!(
                navigator_domain::testtype::by_code(&code).is_some(),
                "code {code} not in catalog"
            );
        }
        for plat in ["PACBIO", "ILLUMINA", "NANOPORE"] {
            let code = infer_test_type(None, Some(plat), None, None).unwrap();
            assert!(
                navigator_domain::testtype::by_code(&code).is_some(),
                "code {code} not in catalog"
            );
        }
    }

    #[test]
    fn no_profile_falls_back_to_platform() {
        assert_eq!(
            infer_test_type(None, Some("PACBIO"), None, None).as_deref(),
            Some("WGS_HIFI")
        );
        assert_eq!(
            infer_test_type(None, Some("ILLUMINA"), None, None).as_deref(),
            Some("WGS")
        );
        assert_eq!(infer_test_type(None, None, None, None), None);
    }
}
