//! Heteroplasmy-detection tests against the `diploid.bam` fixture (tests/fixtures).
//!
//! That fixture is two haplotypes on chr1 at depth 20 (10 reads each):
//!   H1 = ACGTACGAAC, H2 = AGGTTCGAAC
//! so the per-position pileup carries two alleles only at pos2 (C/G) and pos5 (A/T);
//! every other position is homozygous. With the default screening params (min_depth 20,
//! minor fraction ≥ 0.03, ≥3 minor reads) detection must flag exactly those two sites,
//! each at a 50% minor fraction.

use std::path::PathBuf;

use navigator_analysis::heteroplasmy::{detect_heteroplasmy, HeteroplasmyParams};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

#[test]
fn detects_the_two_mixed_sites_on_the_diploid_fixture() {
    let sites = detect_heteroplasmy(
        &fixtures().join("diploid.bam"),
        "chr1",
        &HeteroplasmyParams::default(),
        None,
    )
    .unwrap();

    let positions: Vec<i64> = sites.iter().map(|s| s.position).collect();
    assert_eq!(positions, vec![2, 5], "only pos2 (C/G) and pos5 (A/T) are mixed");

    for s in &sites {
        assert_eq!(s.depth, 20);
        assert_eq!(s.minor_count, 10);
        assert!(
            (s.minor_fraction - 0.5).abs() < 1e-9,
            "50% minor fraction, got {}",
            s.minor_fraction
        );
    }
    // pos2: C major (tie broken to the earlier base), G minor.
    assert_eq!(sites[0].major_base, 'C');
    assert_eq!(sites[0].minor_base, 'G');
    // pos5: A major, T minor.
    assert_eq!(sites[1].major_base, 'A');
    assert_eq!(sites[1].minor_base, 'T');
}

#[test]
fn min_minor_count_suppresses_low_support() {
    // Demanding more minor reads than the fixture supplies (10) yields nothing.
    let strict = HeteroplasmyParams {
        min_minor_count: 11,
        ..HeteroplasmyParams::default()
    };
    let sites = detect_heteroplasmy(&fixtures().join("diploid.bam"), "chr1", &strict, None).unwrap();
    assert!(sites.is_empty());
}
