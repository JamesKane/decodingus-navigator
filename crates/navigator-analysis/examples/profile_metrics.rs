//! Profiling harness for the unified quality-metrics walker's hot per-read loop.
//!
//!   cargo build --release --example profile_metrics -p navigator-analysis
//!   ./target/release/examples/profile_metrics <bam> <reference.fa> <contig>
//!
//! Times one contig in three passes (raw decode / +RecordBuf copy / +metrics) so the per-read
//! cost splits out, and completes in ~a minute instead of walking the whole genome for hours.

use std::path::Path;

fn ns_per_read(d: std::time::Duration, n: u64) -> f64 {
    if n > 0 { d.as_nanos() as f64 / n as f64 } else { 0.0 }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: profile_metrics <bam> <reference.fa> <contig>");
        std::process::exit(2);
    }
    let params = navigator_analysis::coverage::CallableLociParams::default();

    // "FULL" → run the real production walker on the whole BAM with per-contig timestamps.
    if args[3].eq_ignore_ascii_case("full") {
        let start = std::time::Instant::now();
        let progress = move |done: usize, total: usize| {
            eprintln!("  [{:>7.1?}] coverage contig {done}/{total}", start.elapsed());
        };
        match navigator_analysis::unified::collect_unified_metrics_parallel_with_progress(
            Path::new(&args[1]),
            Path::new(&args[2]),
            &params,
            None,
            &progress,
        ) {
            Ok(r) => eprintln!(
                "\nFULL walker done in {:.1?}: mean_cov={:.1} callable={} reads={:?}",
                start.elapsed(),
                r.coverage.mean_coverage,
                r.coverage.callable_bases,
                r.read_metrics.total_reads
            ),
            Err(e) => {
                eprintln!("error: {e}");
                std::process::exit(1);
            }
        }
        return;
    }

    // A comma-separated contig list → parallel mode: run them concurrently (like the real walker)
    // and report aggregate throughput, to expose contention vs the single-threaded rate.
    if args[3].contains(',') {
        let contigs: Vec<String> = args[3].split(',').map(|s| s.to_string()).collect();
        let nthreads = rayon::current_num_threads();
        match navigator_analysis::unified::profile_contigs_parallel(Path::new(&args[1]), Path::new(&args[2]), &contigs, &params) {
            Ok((n, dur)) => {
                let mreads_s = if dur.as_secs_f64() > 0.0 { n as f64 / dur.as_secs_f64() / 1e6 } else { 0.0 };
                eprintln!(
                    "parallel ({} contigs, {} rayon threads): {} reads in {:.2?}  →  {:.2} M reads/s aggregate",
                    contigs.len(), nthreads, n, dur, mreads_s
                );
            }
            Err(e) => {
                eprintln!("error: {e}");
                std::process::exit(1);
            }
        }
        return;
    }

    match navigator_analysis::unified::profile_contig(Path::new(&args[1]), Path::new(&args[2]), &args[3], &params) {
        Ok(p) => {
            let n = p.reads;
            eprintln!("contig {} — {} reads\n", args[3], n);
            eprintln!("  raw decode      {:>8.2?}   {:>6.0} ns/read", p.raw, ns_per_read(p.raw, n));
            eprintln!("  + RecordBuf     {:>8.2?}   {:>6.0} ns/read   (+{:.0} ns/read for the owned copy)",
                p.recordbuf, ns_per_read(p.recordbuf, n), ns_per_read(p.recordbuf, n) - ns_per_read(p.raw, n));
            eprintln!("  + metrics       {:>8.2?}   {:>6.0} ns/read   (+{:.0} ns/read for read-metrics + pileup)",
                p.full, ns_per_read(p.full, n), ns_per_read(p.full, n) - ns_per_read(p.recordbuf, n));
            let mreads_s = if p.full.as_secs_f64() > 0.0 { n as f64 / p.full.as_secs_f64() / 1e6 } else { 0.0 };
            eprintln!("\n  full pass: {:.2} M reads/s single-threaded", mreads_s);
        }
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    }
}
