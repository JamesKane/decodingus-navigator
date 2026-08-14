//! Profiling harness for the mapping stage alone — stage B with nothing else in the sample.
//!
//! The mapping stage runs read → map → write per batch, and the CPU-load graph shows a valley
//! between the peaks: the rayon pool idles while one thread inflates gzip and parses the next
//! batch of reads. This runs `map_pairs` and nothing else, so a profile attributes that valley to
//! a function rather than to a stage.
//!
//! ```sh
//! cargo build --profile profiling -p navigator-align --example map_profile
//! R1=$HOME/.decodingus/profiling/reverted_1.fastq.gz \
//! R2=$HOME/.decodingus/profiling/reverted_2.fastq.gz \
//! INDEX=$HOME/.decodingus/minimap2_index/chm13v2.0/sr.mmi \
//! OUT=$HOME/.decodingus/profiling/mapped.bam \
//!   samply record target/profiling/examples/map_profile
//! ```
//!
//! `LIMIT_SECONDS` stops after a fixed wall-clock budget, so a profile can be taken over a couple
//! of minutes of steady state instead of the three quarters of an hour a whole sample takes.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use navigator_align::{MapParams, OutputFormat, Preset};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let r1 = PathBuf::from(std::env::var("R1")?);
    let r2 = PathBuf::from(std::env::var("R2")?);
    let index = PathBuf::from(std::env::var("INDEX")?);
    let out = PathBuf::from(std::env::var("OUT").unwrap_or_else(|_| "/tmp/map_profile.bam".into()));
    let scratch = out.with_extension("scratch");
    let limit: u64 = std::env::var("LIMIT_SECONDS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(180);

    let params = MapParams {
        preset: Preset::ShortRead,
        threads: 0,
        read_group: None,
        format: OutputFormat::Bam,
        reference: None,
    };

    let started = Instant::now();
    let stop = AtomicBool::new(false);
    let cancelled = || {
        if started.elapsed().as_secs() >= limit {
            stop.store(true, Ordering::Relaxed);
        }
        stop.load(Ordering::Relaxed)
    };

    let mut last = Instant::now();
    let mut reported = 0u64;
    let mut progress = |queries: u64, part: usize, parts: usize| {
        if last.elapsed().as_secs() >= 10 {
            let rate = (queries - reported) as f64 / last.elapsed().as_secs_f64();
            eprintln!(
                "[{:>6.1}s] part {part}/{parts} · {queries} reads · {rate:.0} reads/s",
                started.elapsed().as_secs_f64(),
            );
            last = Instant::now();
            reported = queries;
        }
    };

    let stats = match navigator_align::map_pairs(&index, &r1, &r2, &out, &scratch, &params, &cancelled, &mut progress) {
        Ok(stats) => stats,
        // Hitting the time budget is the normal way this ends.
        Err(navigator_align::AlignError::Cancelled) => {
            eprintln!("stopped at the {limit}s budget");
            return Ok(());
        }
        Err(e) => return Err(e.into()),
    };

    eprintln!(
        "{} reads in {:.1}s ({:.0} reads/s) — {} mapped, {} unmapped",
        stats.queries,
        started.elapsed().as_secs_f64(),
        stats.queries as f64 / started.elapsed().as_secs_f64(),
        stats.mapped,
        stats.unmapped,
    );
    Ok(())
}
