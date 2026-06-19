//! Coverage-walker tests against the synthetic fixture (see tests/fixtures/make_fixture.sh).
//!
//! Reference chrM (50 bp, N at position 25). Reads:
//!   pos 1-10   depth 4  MAPQ 60  -> CALLABLE              (10 bp)
//!   pos 11-20  depth 2  MAPQ 60  -> LOW_COVERAGE          (10 bp)
//!   pos 21-24  depth 0           -> NO_COVERAGE            (4 bp)
//!   pos 25     depth 0  ref N    -> REF_N                  (1 bp)
//!   pos 26-30  depth 5  MAPQ 0   -> POOR_MAPPING_QUALITY   (5 bp)
//!   pos 31-50  depth 0           -> NO_COVERAGE           (20 bp)
//! Base quality is Phred 40 throughout.

use std::path::PathBuf;

use navigator_analysis::coverage::{collect_coverage_callable, CallableLociParams};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn approx(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "expected {b}, got {a}");
}

#[test]
fn coverage_matches_hand_computed_values() {
    let dir = fixtures();
    let result = collect_coverage_callable(
        &dir.join("coverage.bam"),
        &dir.join("ref.fa"),
        &CallableLociParams::default(),
        None,
    )
    .expect("walker should succeed");

    // --- global coverage ---
    // depth sum = 10*4 + 10*2 + 5*5 = 85 over 50 positions
    assert_eq!(result.genome_territory, 50);
    approx(result.mean_coverage, 85.0 / 50.0); // 1.7
    // histogram: 25 positions @0, 10 @2, 10 @4, 5 @5
    assert_eq!(result.coverage_histogram[0], 25);
    assert_eq!(result.coverage_histogram[2], 10);
    assert_eq!(result.coverage_histogram[4], 10);
    assert_eq!(result.coverage_histogram[5], 5);
    // median: cumulative reaches half (25) at depth 0
    approx(result.median_coverage, 0.0);
    // population sd: sqrt(325/50 - 1.7^2) = sqrt(3.61) = 1.9
    approx(result.sd_coverage, 1.9);
    // pct depth>=1 = 25/50; >=5 = 5/50; >=10 = 0
    approx(result.pct_1x, 0.5);
    approx(result.pct_5x, 0.1);
    approx(result.pct_10x, 0.0);

    // --- callable states ---
    assert_eq!(result.contig_callable.len(), 1);
    let cm = &result.contig_callable[0];
    assert_eq!(cm.contig, "chrM");
    assert_eq!(cm.callable, 10);
    assert_eq!(cm.low_coverage, 10);
    assert_eq!(cm.no_coverage, 24); // 4 + 20
    assert_eq!(cm.poor_mapping_quality, 5);
    assert_eq!(cm.ref_n, 1);
    assert_eq!(
        cm.callable + cm.low_coverage + cm.no_coverage + cm.poor_mapping_quality + cm.ref_n + cm.excessive_coverage,
        50
    );
    assert_eq!(result.callable_bases, 10);

    // --- samtools-style per-contig stats (per-base-observation averaging) ---
    let cs = &result.contig_coverage_stats[0];
    assert_eq!(cs.num_reads, 11); // 4 + 2 + 5
    assert_eq!(cs.cov_bases, 25);
    approx(cs.coverage, 50.0);
    approx(cs.mean_depth, 1.7);
    approx(cs.mean_base_q, 40.0); // all bases Phred 40
    // map quality per base obs: (60 obs * 60 + 25 obs * 0) / 85
    approx(cs.mean_map_q, 3600.0 / 85.0);
}

#[test]
fn cram_coverage_matches_bam_field_for_field() {
    // The CRAM is the same reads as coverage.bam, ref-compressed; the reader unification
    // (both decoded to RecordBuf) must yield an identical CoverageResult.
    let dir = fixtures();
    let params = CallableLociParams::default();
    let bam = collect_coverage_callable(&dir.join("coverage.bam"), &dir.join("ref.fa"), &params, None).unwrap();
    let cram = collect_coverage_callable(&dir.join("coverage.cram"), &dir.join("ref.fa"), &params, None).unwrap();
    assert_eq!(cram, bam, "CRAM coverage must equal BAM coverage");
}

#[test]
fn contig_allowlist_excludes_unlisted_contigs() {
    let dir = fixtures();
    let allow = std::collections::HashSet::new(); // empty: chrM not listed
    let result = collect_coverage_callable(
        &dir.join("coverage.bam"),
        &dir.join("ref.fa"),
        &CallableLociParams::default(),
        Some(&allow),
    )
    .expect("walker should succeed");

    assert_eq!(result.genome_territory, 0);
    assert!(result.contig_callable.is_empty());
    assert!(result.contig_coverage_stats.is_empty());
}
