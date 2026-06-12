//! End-to-end validation of the **standalone** coverage path (`collect_coverage_callable`, the
//! sequential walker that now consumes lazy `bam::Record` via `records_lazy`) against the trusted
//! per-contig **parallel** walker, on a real BAM. Both produce a `CoverageResult`; the invariant is
//! exact equality (the same one the `unified_matches_standalone_walkers` unit test asserts on a
//! fixture — this runs it on a whole WGS).
//!
//!   cargo build --release --example validate_coverage -p navigator-analysis
//!   ./target/release/examples/validate_coverage <bam> <reference.fa>

use std::path::Path;
use std::time::Instant;

use navigator_analysis::coverage::{self, CallableLociParams, CoverageResult};
use navigator_analysis::unified;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: validate_coverage <bam> <reference.fa>");
        std::process::exit(2);
    }
    let bam = Path::new(&args[1]);
    let reference = Path::new(&args[2]);
    let params = CallableLociParams::default();

    // 1. The path under test: standalone sequential coverage (lazy records).
    let t0 = Instant::now();
    let mut last = 0usize;
    let mut progress = |done: usize, total: usize| {
        if done != last {
            eprintln!("  [{:>7.1?}] standalone coverage contig {done}/{total}", t0.elapsed());
            last = done;
        }
    };
    let standalone = match coverage::collect_coverage_callable_with_progress(
        bam, reference, &params, None, &mut progress,
    ) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("standalone coverage error: {e}");
            std::process::exit(1);
        }
    };
    let standalone_dur = t0.elapsed();
    eprintln!("\nstandalone coverage done in {standalone_dur:.1?}");
    summarize("standalone", &standalone);

    // 2. The oracle: trusted per-contig parallel walker on the same file.
    let t1 = Instant::now();
    let progress2 = |_done: usize, _total: usize| {};
    let unified = match unified::collect_unified_metrics_parallel_with_progress(
        bam, reference, &params, None, &progress2,
    ) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("parallel walker error: {e}");
            std::process::exit(1);
        }
    };
    let parallel_dur = t1.elapsed();
    eprintln!("\nparallel walker done in {parallel_dur:.1?}");
    summarize("parallel ", &unified.coverage);

    // 3. The invariant: byte-for-byte identical coverage.
    eprintln!("\n=== comparison ===");
    if standalone == unified.coverage {
        eprintln!("PASS — standalone == parallel (field-for-field identical)");
    } else {
        eprintln!("FAIL — coverage results differ:");
        diff(&standalone, &unified.coverage);
        std::process::exit(1);
    }
}

fn summarize(tag: &str, c: &CoverageResult) {
    eprintln!(
        "  [{tag}] territory={} mean={:.4} median={:.1} sd={:.4} callable={} pct1x={:.4} pct10x={:.4} pct30x={:.4} contigs={}",
        c.genome_territory, c.mean_coverage, c.median_coverage, c.sd_coverage, c.callable_bases,
        c.pct_1x, c.pct_10x, c.pct_30x, c.contig_coverage_stats.len()
    );
}

fn diff(a: &CoverageResult, b: &CoverageResult) {
    macro_rules! d {
        ($f:ident) => {
            if a.$f != b.$f {
                eprintln!("  {}: standalone={:?} parallel={:?}", stringify!($f), a.$f, b.$f);
            }
        };
    }
    d!(genome_territory);
    d!(mean_coverage);
    d!(median_coverage);
    d!(sd_coverage);
    d!(callable_bases);
    d!(pct_1x);
    d!(pct_5x);
    d!(pct_10x);
    d!(pct_20x);
    d!(pct_30x);
    if a.coverage_histogram != b.coverage_histogram {
        let n = a.coverage_histogram.len().max(b.coverage_histogram.len());
        let mut shown = 0;
        for i in 0..n {
            let (x, y) = (a.coverage_histogram.get(i), b.coverage_histogram.get(i));
            if x != y && shown < 10 {
                eprintln!("  hist[{i}]: standalone={x:?} parallel={y:?}");
                shown += 1;
            }
        }
    }
    if a.contig_coverage_stats != b.contig_coverage_stats {
        for (sa, pa) in a.contig_coverage_stats.iter().zip(&b.contig_coverage_stats) {
            if sa != pa {
                eprintln!("  contig stats differ: standalone={sa:?} parallel={pa:?}");
            }
        }
    }
}
