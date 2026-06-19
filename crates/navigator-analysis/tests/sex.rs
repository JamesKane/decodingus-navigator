//! sex inference test against sex.bam: chr1 (100 bp, 10 reads) + chrX (100 bp, 2 reads)
//! -> autosome 10x, chrX 2x, ratio 0.2 -> Male, high confidence.

use std::path::PathBuf;

use navigator_analysis::sex::{infer_from_bam, Confidence, InferredSex};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

#[test]
fn infers_male_from_low_x_coverage() {
    let r = infer_from_bam(&fixtures().join("sex.bam"), None).expect("should succeed");

    assert_eq!(r.inferred_sex, InferredSex::Male);
    assert_eq!(r.confidence, Confidence::High);
    assert!((r.x_autosome_ratio - 0.2).abs() < 1e-9, "ratio {}", r.x_autosome_ratio);
    assert!((r.autosome_mean_coverage - 10.0).abs() < 1e-9);
    assert!((r.x_coverage - 2.0).abs() < 1e-9);
}

#[test]
fn cram_sex_inference_matches_bam() {
    // CRAM has no per-reference counts in the index, so this exercises the record-scan
    // fallback; same reads as sex.bam, so the result must match.
    let dir = fixtures();
    let bam = infer_from_bam(&dir.join("sex.bam"), None).unwrap();
    let cram = infer_from_bam(&dir.join("sex.cram"), Some(&dir.join("sexref.fa"))).unwrap();
    assert_eq!(cram, bam, "CRAM sex inference must equal BAM");
}
