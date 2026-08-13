//! Drive a whole-genome realignment headlessly, for the phase 5 WGS-scale validation.
//!
//! The GUI can start this job, but a run measured in hours should not depend on a window staying
//! open — and the validation wants a timestamped log of where the time went, which the progress
//! cards do not keep. This is the same `App::realign_alignment` the UI calls, with the stage
//! reports printed instead of drawn.
//!
//! ```bash
//! cargo run --release -p navigator-app --example realign_wgs -- <alignment_id>
//! # optional: REF=… BUILD=… SCRATCH=… (defaults: cached chm13v2.0.fa, "chm13v2.0", beside the output)
//! ```
//!
//! Ctrl-C cancels through the job's own token rather than killing the process, so the scratch
//! directory — hundreds of GB at WGS scale — is still cleaned up on the way out.

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
    // `PRESET` overrides the technology inference, which refuses any test type it does not know
    // rather than guessing — correct for the app, but it puts real vendor products (`Y_ELITE`)
    // out of reach of a smoke test.
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
    };

    match app.realign_alignment(alignment_id, params, cancel, progress).await {
        Ok(outcome) => {
            println!(
                "\ndone in {:.1} min\n  alignment #{} at {}\n  reads written: {}\n  duplicates marked: {}\n  sequenceless non-primary records dropped: {}\n  source unmapped reads (had a chance to place): {}",
                started.elapsed().as_secs_f64() / 60.0,
                outcome.alignment.id,
                outcome.alignment.bam_path.as_deref().unwrap_or("?"),
                outcome.reads_written,
                outcome.duplicates_marked,
                outcome.sequenceless_dropped,
                outcome.source_unmapped_reads,
            );
            Ok(())
        }
        Err(e) => {
            eprintln!("\nfailed after {:.1} min: {e}", started.elapsed().as_secs_f64() / 60.0);
            Err(e.into())
        }
    }
}
