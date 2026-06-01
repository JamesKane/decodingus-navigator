//! Real-data smoke test (seed of the phase-3 §4c parity harness). Ignored by default;
//! runs only when pointed at a local BAM + reference via env vars:
//!
//!   HG002_CHRM_BAM=/tmp/hg002.chrM.bam CHM13_REF=/Users/.../chm13v2.0.fa \
//!     cargo test -p navigator-analysis --test parity_real -- --ignored --nocapture
//!
//! This is a sanity check that noodles handles a real BAM (varied CIGARs/MAPQ) and the
//! chrM numbers are plausible — NOT strict parity, which is measured against the Scala
//! walker / GATK in phase 3.

use std::collections::HashSet;
use std::path::PathBuf;

use navigator_analysis::coverage::{collect_coverage_callable, CallableLociParams};

#[test]
#[ignore = "requires local HG002_CHRM_BAM + CHM13_REF env vars"]
fn hg002_chrm_smoke() {
    let (Ok(bam), Ok(reference)) = (std::env::var("HG002_CHRM_BAM"), std::env::var("CHM13_REF"))
    else {
        eprintln!("set HG002_CHRM_BAM and CHM13_REF to run this test");
        return;
    };

    let allow: HashSet<String> = ["chrM".to_string()].into_iter().collect();
    let result = collect_coverage_callable(
        &PathBuf::from(bam),
        &PathBuf::from(reference),
        &CallableLociParams::default(),
        Some(&allow),
    )
    .expect("walker should succeed on real data");

    eprintln!("genome_territory = {}", result.genome_territory);
    eprintln!("mean_coverage    = {:.3}", result.mean_coverage);
    eprintln!("median_coverage  = {}", result.median_coverage);
    eprintln!("sd_coverage      = {:.3}", result.sd_coverage);
    eprintln!("pct_10x/20x/30x  = {:.4} / {:.4} / {:.4}", result.pct_10x, result.pct_20x, result.pct_30x);
    eprintln!("callable_bases   = {}", result.callable_bases);
    eprintln!("callable metrics = {:?}", result.contig_callable);
    eprintln!("coverage stats   = {:?}", result.contig_coverage_stats);

    // chrM should be fully covered at high depth.
    assert_eq!(result.genome_territory, 16569);
    let cs = &result.contig_coverage_stats[0];
    assert_eq!(cs.contig, "chrM");
    assert_eq!(cs.cov_bases, 16569);
    assert!((cs.coverage - 100.0).abs() < 1e-9);
    assert!(cs.mean_depth > 50.0, "mean depth {} unexpectedly low", cs.mean_depth);

    let cm = &result.contig_callable[0];
    let total = cm.callable + cm.low_coverage + cm.no_coverage + cm.poor_mapping_quality
        + cm.ref_n + cm.excessive_coverage;
    assert_eq!(total, 16569);
}
