//! Drive a whole-genome realignment headlessly, for the phase 5 WGS-scale validation.
//!
//! The GUI can start this job. But a run of many hours must not depend on an open window. The
//! validation also needs a log with a timestamp for each stage, and the progress cards do not keep
//! one.
//!
//! This example calls the same `App::realign_alignment` function that the UI calls. It prints each
//! stage report, and the UI draws it.
//!
//! ```bash
//! cargo run --release -p navigator-app --example realign_wgs -- <alignment_id>
//! # optional: REF=… BUILD=… SCRATCH=… (defaults: cached chm13v2.0.fa, "chm13v2.0", beside the output)
//! # RESUME=1 picks up from a previous attempt's intermediates
//! # NAVIGATOR_REALIGN_KEEP_SCRATCH=1 leaves them behind when a run fails, so RESUME has something
//! #   to find; a run that is killed outright leaves them regardless
//! ```
//!
//! Ctrl-C stops the job through the cancel token of that job. It does not stop the process. So the
//! example still removes the scratch directory, which holds hundreds of GB for a WGS sample.

use std::path::PathBuf;
use std::sync::Mutex;
use std::time::Instant;

use navigator_align::Preset;
use navigator_app::realign_job::{RealignParams, RealignStage};
use navigator_app::{App, CancelToken};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let alignment_id: i64 = std::env::args().nth(1).unwrap_or_else(|| "5".into()).parse()?;
    let home = PathBuf::from(std::env::var("HOME")?);
    let db = home.join(".decodingus/navigator-rs.db");
    let build = std::env::var("BUILD").unwrap_or_else(|_| "chm13v2.0".to_string());
    let reference = std::env::var("REF")
        .map(PathBuf::from)
        .unwrap_or_else(|_| home.join(".decodingus/references/chm13v2.0.fa"));
    let scratch_root = std::env::var("SCRATCH").ok().map(PathBuf::from);
    // `PRESET` replaces the value that the code deduces from the technology. That code refuses a
    // test type that it does not know, and it never makes an estimate. This behaviour is correct
    // for the app. But it also puts a real vendor product, such as `Y_ELITE`, out of the reach of
    // a quick test.
    let preset = match std::env::var("PRESET") {
        Ok(p) => Some(Preset::parse(&p).map_err(|e| format!("PRESET={p}: {e}"))?),
        Err(_) => None,
    };

    if !reference.exists() {
        return Err(format!("reference not found: {}", reference.display()).into());
    }

    let app = App::open(&db).await?;
    let source = app.alignment(alignment_id).await?.ok_or("no such alignment")?;
    let source_path = source.bam_path.clone().ok_or("source alignment has no file")?;
    if !std::path::Path::new(&source_path).exists() {
        return Err(format!("source alignment #{alignment_id} points at a missing file: {source_path}").into());
    }

    println!(
        "realigning alignment #{alignment_id} ({} {}) → {build}\n  source: {source_path}\n  target: {}",
        source.reference_build,
        source.aligner,
        reference.display(),
    );

    let cancel = CancelToken::new();
    {
        let cancel = cancel.clone();
        tokio::spawn(async move {
            if tokio::signal::ctrl_c().await.is_ok() {
                eprintln!("\ncancelling — the scratch directory is removed before exit");
                cancel.cancel();
            }
        });
    }

    let started = Instant::now();
    // Stage boundaries are the measurement: which stage owns the wall clock is the thing a
    // multi-hour job's log has to answer.
    let stage_clock = Mutex::new((Instant::now(), None::<RealignStage>));
    let progress = move |p: navigator_app::realign_job::RealignProgress| {
        {
            let mut clock = stage_clock.lock().expect("stage clock");
            let (last, current) = &mut *clock;
            if *current != Some(p.stage) {
                if let Some(previous) = *current {
                    println!("  [{:>8.1}s] {} done", last.elapsed().as_secs_f64(), previous.label());
                }
                *last = Instant::now();
                *current = Some(p.stage);
            }
        }
        let detail = if p.detail.is_empty() {
            String::new()
        } else {
            format!(" — {}", p.detail)
        };
        println!(
            "[{:>8.1}s] stage {}/{}: {}{detail}",
            started.elapsed().as_secs_f64(),
            p.stage.step(),
            p.total_stages,
            p.stage.label(),
        );
    };

    let params = RealignParams {
        target_build: build,
        target_reference: reference,
        preset,
        scratch_root,
        resume: std::env::var("RESUME")
            .map(|v| v != "0" && !v.is_empty())
            .unwrap_or(false),
    };

    match app.realign_alignment(alignment_id, params, cancel, progress).await {
        Ok(outcome) => {
            // A run that continues an earlier run does not always do the stage that counts a
            // figure. The earlier run can also stop before it writes that figure. So the report
            // shows "not measured". A zero value here looks like a result.
            let count = |n: Option<u64>| n.map(|n| n.to_string()).unwrap_or_else(|| "not measured".into());
            println!(
                "\ndone in {:.1} min\n  alignment #{} at {}\n  reads written: {}\n  duplicates marked: {}\n  source unmapped reads (had a chance to place): {}",
                started.elapsed().as_secs_f64() / 60.0,
                outcome.alignment.id,
                outcome.alignment.bam_path.as_deref().unwrap_or("?"),
                count(outcome.reads_written),
                count(outcome.duplicates_marked),
                count(outcome.source_unmapped_reads),
            );
            Ok(())
        }
        Err(e) => {
            eprintln!("\nfailed after {:.1} min: {e}", started.elapsed().as_secs_f64() / 60.0);
            Err(e.into())
        }
    }
}
