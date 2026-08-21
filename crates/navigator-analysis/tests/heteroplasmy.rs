//! Tests of the heteroplasmy detection, against the `diploid.bam` fixture in tests/fixtures.
//!
//! That fixture holds two haplotypes on chr1, at a depth of 20, with 10 reads for each:
//!
//! ```text
//! H1 = ACGTACGAAC
//! H2 = AGGTTCGAAC
//! ```
//!
//! The pileup thereby carries two alleles at two positions alone: pos2, at C and G, and pos5, at A
//! and T. Every other position is homozygous.
//!
//! The default parameters are a min_depth of 20, a minor fraction of 0.03 or more, and 3 minor
//! reads or more. With those, the detection must flag exactly those two sites, and each one must
//! show a minor fraction of 50%.

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
    // A request for more minor reads than the fixture holds, which is 10, gives nothing.
    let strict = HeteroplasmyParams {
        min_minor_count: 11,
        ..HeteroplasmyParams::default()
    };
    let sites = detect_heteroplasmy(&fixtures().join("diploid.bam"), "chr1", &strict, None).unwrap();
    assert!(sites.is_empty());
}
