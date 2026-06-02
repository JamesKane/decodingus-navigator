//! Sync↔async bridge. egui runs the immediate-mode loop on the main thread; the
//! `App` (tokio + sqlx) runs on a dedicated worker thread with its own runtime. The UI
//! sends [`Command`]s and drains [`Event`]s each frame — no DB calls or domain
//! decisions on the UI thread (plan §5).
//!
//! Each command is handled on its own task so a long analysis run never blocks quick
//! queries. The command→event mapping ([`handle`]) is pure and unit-tested; [`spawn`]
//! is the thread/runtime/channel glue.

use std::path::PathBuf;
use std::sync::mpsc::{Receiver, Sender};
use std::sync::Arc;

use navigator_app::{App, Coverage, ProjectOverview};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Alignment, Biosample, NewProject, Project};
use tokio::sync::mpsc::{unbounded_channel, UnboundedSender};

/// A request from the UI to the worker.
#[derive(Debug, Clone)]
pub enum Command {
    LoadOverview,
    CreateProject(NewProject),
    LoadSamples(i64),
    LoadAlignments(SampleGuid),
    LoadCoverage(i64),
    RunCoverage(i64),
}

/// A result/notification from the worker to the UI.
#[derive(Debug, Clone)]
pub enum Event {
    Overview(Vec<ProjectOverview>),
    ProjectCreated(Project),
    Samples { project_id: i64, samples: Vec<Biosample> },
    Alignments { biosample_guid: SampleGuid, alignments: Vec<Alignment> },
    Coverage { alignment_id: i64, result: Option<Coverage> },
    Error(String),
}

/// Execute one command against the app, mapping success/failure to an [`Event`].
pub async fn handle(app: &App, cmd: Command) -> Event {
    match cmd {
        Command::LoadOverview => match app.project_overview().await {
            Ok(v) => Event::Overview(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::CreateProject(new) => match app.create_project(new).await {
            Ok(p) => Event::ProjectCreated(p),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSamples(project_id) => match app.list_biosamples(project_id).await {
            Ok(samples) => Event::Samples { project_id, samples },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAlignments(biosample_guid) => match app.list_alignments(biosample_guid).await {
            Ok(alignments) => Event::Alignments { biosample_guid, alignments },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadCoverage(alignment_id) => match app.cached_coverage(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunCoverage(alignment_id) => match app.run_coverage_for_alignment(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        },
    }
}

/// Spawn the worker thread: open the workspace at `db_path` inside the worker's runtime
/// (so the connection pool lives there), then serve commands. `wake` is called after
/// each event so the UI can `request_repaint`. Returns the command sender and event
/// receiver the UI holds.
pub fn spawn(
    db_path: PathBuf,
    wake: impl Fn() + Send + Sync + 'static,
) -> (UnboundedSender<Command>, Receiver<Event>) {
    let (cmd_tx, mut cmd_rx) = unbounded_channel::<Command>();
    let (evt_tx, evt_rx) = std::sync::mpsc::channel::<Event>();
    let wake = Arc::new(wake);

    std::thread::Builder::new()
        .name("navigator-worker".into())
        .spawn(move || {
            let rt = match tokio::runtime::Runtime::new() {
                Ok(rt) => rt,
                Err(e) => {
                    let _ = evt_tx.send(Event::Error(format!("runtime: {e}")));
                    wake();
                    return;
                }
            };
            rt.block_on(async move {
                let app = match App::open(&db_path).await {
                    Ok(app) => app,
                    Err(e) => {
                        let _ = evt_tx.send(Event::Error(format!("open workspace: {e}")));
                        wake();
                        return;
                    }
                };
                while let Some(cmd) = cmd_rx.recv().await {
                    let app = app.clone();
                    let evt_tx: Sender<Event> = evt_tx.clone();
                    let wake = wake.clone();
                    tokio::spawn(async move {
                        let event = handle(&app, cmd).await;
                        let _ = evt_tx.send(event);
                        wake();
                    });
                }
            });
        })
        .expect("spawn worker thread");

    (cmd_tx, evt_rx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_store::Store;

    async fn app() -> App {
        App::new(Store::open_in_memory().await.unwrap())
    }

    #[tokio::test]
    async fn handle_maps_commands_to_events() {
        let app = app().await;

        // create a project
        let created = handle(&app, Command::CreateProject(NewProject {
            name: "Trio".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await;
        let pid = match created {
            Event::ProjectCreated(p) => p.id,
            other => panic!("expected ProjectCreated, got {other:?}"),
        };

        // overview reflects it
        match handle(&app, Command::LoadOverview).await {
            Event::Overview(v) => {
                assert_eq!(v.len(), 1);
                assert_eq!(v[0].sample_count, 0);
            }
            other => panic!("expected Overview, got {other:?}"),
        }

        // samples for the project (empty)
        match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { project_id, samples } => {
                assert_eq!(project_id, pid);
                assert!(samples.is_empty());
            }
            other => panic!("expected Samples, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn handle_runs_coverage_for_a_stored_alignment() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};
        use std::path::PathBuf;

        let app = app().await;
        let b = app.add_biosample(None, "HG002", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                biosample_guid: b.guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await
            .unwrap();
        let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../navigator-analysis/tests/fixtures");
        let aln = app
            .record_alignment(NewAlignment {
                sequence_run_id: run.id,
                reference_build: "chrM".into(),
                aligner: "synthetic".into(),
                variant_caller: None,
                bam_path: Some(fixtures.join("coverage.bam").to_string_lossy().into_owned()),
                reference_path: Some(fixtures.join("ref.fa").to_string_lossy().into_owned()),
            })
            .await
            .unwrap();

        // alignments query
        match handle(&app, Command::LoadAlignments(b.guid)).await {
            Event::Alignments { alignments, .. } => assert_eq!(alignments, vec![aln.clone()]),
            other => panic!("expected Alignments, got {other:?}"),
        }

        // cold cache
        match handle(&app, Command::LoadCoverage(aln.id)).await {
            Event::Coverage { result, .. } => assert!(result.is_none()),
            other => panic!("expected Coverage(None), got {other:?}"),
        }

        // run + persist (uses the alignment's stored paths, via spawn_blocking)
        match handle(&app, Command::RunCoverage(aln.id)).await {
            Event::Coverage { alignment_id, result } => {
                assert_eq!(alignment_id, aln.id);
                assert_eq!(result.unwrap().genome_territory, 50);
            }
            other => panic!("expected Coverage(Some), got {other:?}"),
        }

        // now cached
        match handle(&app, Command::LoadCoverage(aln.id)).await {
            Event::Coverage { result, .. } => assert_eq!(result.unwrap().callable_bases, 10),
            other => panic!("expected cached Coverage, got {other:?}"),
        }
    }
}
