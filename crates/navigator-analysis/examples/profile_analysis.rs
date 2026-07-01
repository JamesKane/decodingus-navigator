//! Profile the per-step analysis cost on a real BAM/CRAM to find batch-analysis hotspots.
//! Times the walkers the deep-analyze pipeline drives, whole-genome vs targeted-Y scoped, and the
//! sequential vs the parallel coverage path — plus a chrY region-query genotyping pass (the
//! haplogroup step). Read-only; nothing persisted.
//!
//!   cargo run --release --example profile_analysis -p navigator-analysis -- <bam|cram> <ref.fa>
use std::collections::HashSet;
use std::path::Path;
use std::time::Instant;

use navigator_analysis::caller::{self, HaploidCallerParams};
use navigator_analysis::coverage::{self, CallableLociParams};
use navigator_analysis::unified;

fn timed<T>(label: &str, f: impl FnOnce() -> T) -> T {
    let t = Instant::now();
    let r = f();
    eprintln!("  {:>10.2?}   {label}", t.elapsed());
    r
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: profile_analysis <bam|cram> <ref.fa>");
        std::process::exit(2);
    }
    let bam = Path::new(&args[1]);
    let reference = Path::new(&args[2]);
    let params = CallableLociParams::default();
    let ym: HashSet<String> = ["chrY", "chrM"].iter().map(|s| s.to_string()).collect();
    let mb = std::fs::metadata(bam).map(|m| m.len() / 1_000_000).unwrap_or(0);
    eprintln!("\n=== {} ({mb} MB) ===", bam.display());

    timed("estimate_molecule_lengths (prefix read)", || {
        coverage::estimate_molecule_lengths(bam, Some(reference)).ok()
    });
    timed("coverage SEQUENTIAL  whole-genome", || {
        coverage::collect_coverage_callable(bam, reference, &params, None).map(|_| ()).err()
    });
    timed("coverage SEQUENTIAL  scoped chrY+chrM", || {
        coverage::collect_coverage_callable(bam, reference, &params, Some(&ym)).map(|_| ()).err()
    });
    timed("coverage PARALLEL    whole-genome", || {
        unified::collect_unified_metrics_parallel(bam, reference, &params, None).map(|_| ()).err()
    });
    timed("coverage PARALLEL    scoped chrY+chrM", || {
        unified::collect_unified_metrics_parallel(bam, reference, &params, Some(&ym)).map(|_| ()).err()
    });

    // chrY haplogroup genotyping pass: a region query over chrY tallying ~200k target sites
    // (representative of the Y tree's chrY loci) — the deep-analyze Y step's read pattern.
    let hp = HaploidCallerParams::default();
    let targets: HashSet<i64> = (1..=200_000u32).map(|i| i as i64 * 300).collect();
    timed("chrY genotyping  call_bases_at (200k sites)", || {
        caller::call_bases_at(bam, "chrY", &targets, &hp, Some(reference)).map(|_| ()).err()
    });
}
