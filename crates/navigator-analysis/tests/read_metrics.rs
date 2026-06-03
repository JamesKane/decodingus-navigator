//! read_metrics tests against paired.bam: two FR proper pairs on chrM, read length 10,
//! insert sizes 40 and 30 (first-of-pair only), MAPQ 60.

use std::path::PathBuf;

use navigator_analysis::read_metrics::{collect_read_metrics, PairOrientation};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn approx(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "expected {b}, got {a}");
}

#[test]
fn read_metrics_match_hand_computed_values() {
    let m = collect_read_metrics(&fixtures().join("paired.bam"), None).expect("should succeed");

    // counts
    assert_eq!(m.total_reads, 4);
    assert_eq!(m.pf_reads, 4);
    assert_eq!(m.pf_reads_aligned, 4);
    assert_eq!(m.reads_aligned_in_pairs, 4);
    assert_eq!(m.proper_pairs, 4);
    approx(m.pct_pf_reads_aligned, 1.0);
    approx(m.pct_reads_aligned_in_pairs, 1.0);
    approx(m.pct_proper_pairs, 1.0);
    approx(m.pct_chimeras, 0.0);
    approx(m.mean_mapping_quality, 60.0);

    // read length: all 10
    assert_eq!(m.min_read_length, 10);
    assert_eq!(m.max_read_length, 10);
    approx(m.mean_read_length, 10.0);
    approx(m.median_read_length, 10.0);
    approx(m.std_read_length, 0.0);
    assert_eq!(m.read_length_histogram.get(&10), Some(&4));

    // insert size: 40 and 30 (first-of-pair only)
    assert_eq!(m.min_insert_size, 30);
    assert_eq!(m.max_insert_size, 40);
    approx(m.mean_insert_size, 35.0);
    approx(m.median_insert_size, 30.0);
    approx(m.std_insert_size, 5.0);
    assert_eq!(m.insert_size_histogram.get(&30), Some(&1));
    assert_eq!(m.insert_size_histogram.get(&40), Some(&1));

    assert_eq!(m.pair_orientation, PairOrientation::Fr);
}

#[test]
fn cram_read_metrics_match_bam() {
    // Same pairs as paired.bam, ref-compressed — the unified reader must give identical metrics.
    let dir = fixtures();
    let bam = collect_read_metrics(&dir.join("paired.bam"), None).unwrap();
    let cram = collect_read_metrics(&dir.join("paired.cram"), Some(&dir.join("ref.fa"))).unwrap();
    assert_eq!(cram, bam, "CRAM read-metrics must equal BAM");
}
