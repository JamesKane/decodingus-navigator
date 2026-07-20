//! Cancellation against a real alignment — the claim that matters is wall-clock, not a flag.
//!
//! `#[ignore]` (live file, like the other `*_real` harnesses). Point it at a BAM/CRAM:
//!
//!   NAV_CANCEL_BAM=/path/sample.cram NAV_CANCEL_REF=/path/GRCh38.fa \
//!     cargo test -p navigator-analysis --test cancel_real -- --ignored --nocapture

use std::path::PathBuf;
use std::time::{Duration, Instant};

use navigator_analysis::{coverage::CallableLociParams, unified, CancelToken};

fn env_path(key: &str) -> Option<PathBuf> {
    std::env::var(key).ok().map(PathBuf::from)
}

/// A whole-genome walk over a real WGS file takes minutes. Cancel it a second in and assert it
/// returns in well under that — the entire point of threading the token into the walkers, and the
/// thing that a unit test on the token alone cannot demonstrate.
#[test]
#[ignore]
fn cancelling_a_whole_genome_walk_returns_promptly() {
    let (Some(bam), Some(reference)) = (env_path("NAV_CANCEL_BAM"), env_path("NAV_CANCEL_REF")) else {
        eprintln!("set NAV_CANCEL_BAM and NAV_CANCEL_REF");
        return;
    };

    let token = CancelToken::new();
    let canceller = token.clone();
    std::thread::spawn(move || {
        std::thread::sleep(Duration::from_secs(1));
        canceller.cancel();
    });

    let started = Instant::now();
    let result = unified::collect_unified_metrics_parallel_with_progress(
        &bam,
        &reference,
        &CallableLociParams::default(),
        None,
        &|_, _| {},
        &token,
    );
    let elapsed = started.elapsed();
    eprintln!("returned after {elapsed:.1?}: {result:?}");

    assert!(result.is_err(), "a cancelled walk must not return a partial result as success");
    assert!(
        matches!(result, Err(navigator_analysis::AnalysisError::Cancelled)),
        "must report cancellation, not a generic failure"
    );
    // Generous bound: the contigs already in flight finish their current record batch, and rayon
    // has to unwind the fan-out. Anything near the full walk time means the token is not reaching
    // the record loops.
    assert!(
        elapsed < Duration::from_secs(30),
        "cancel took {elapsed:.1?} — the walk is not polling the token"
    );
}
