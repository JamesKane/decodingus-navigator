//! A cancel against a real alignment. The claim that matters is the wall time, and not the state
//! of a flag.
//!
//! This test carries `#[ignore]`, because it needs a live file, as the other `*_real` harnesses do.
//! Point it at a BAM or a CRAM:
//!
//! ```text
//! NAV_CANCEL_BAM=/path/sample.cram NAV_CANCEL_REF=/path/GRCh38.fa \
//!   cargo test -p navigator-analysis --test cancel_real -- --ignored --nocapture
//! ```

use std::path::PathBuf;
use std::time::{Duration, Instant};

use navigator_analysis::{coverage::CallableLociParams, unified, CancelToken};

fn env_path(key: &str) -> Option<PathBuf> {
    std::env::var(key).ok().map(PathBuf::from)
}

/// A walk over the whole genome, on a real WGS file, takes minutes. This test cancels it one second
/// in, and it asserts that the walk returns well below that time. That is the whole point of the
/// token inside the walkers, and a unit test on the token alone can not show it.
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

    assert!(
        result.is_err(),
        "a cancelled walk must not return a partial result as success"
    );
    assert!(
        matches!(result, Err(navigator_analysis::AnalysisError::Cancelled)),
        "must report cancellation, not a generic failure"
    );
    // The bound is generous. The contigs that already run finish their current batch of records,
    // and rayon must then unwind the fan-out. A time near that of the full walk means that the
    // token does not get to the record loops.
    assert!(
        elapsed < Duration::from_secs(30),
        "cancel took {elapsed:.1?} — the walk is not polling the token"
    );
}
