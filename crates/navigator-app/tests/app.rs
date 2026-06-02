//! App command/query layer tests against an in-memory store.

use navigator_app::{App, AppError};
use navigator_domain::workspace::{NewAlignment, NewProject, NewSequenceRun};
use navigator_store::Store;
use serde::{Deserialize, Serialize};

async fn app() -> App {
    App::new(Store::open_in_memory().await.unwrap())
}

#[tokio::test]
async fn command_flow_and_overview() {
    let app = app().await;
    let p = app
        .create_project(NewProject { name: "Trio".into(), description: None, administrator: "jk".into() })
        .await
        .unwrap();

    let b1 = app.add_biosample(Some(p.id), "HG002", Some("SAMEA1".into()), Some("male".into())).await.unwrap();
    app.add_biosample(Some(p.id), "HG003", None, Some("female".into())).await.unwrap();

    let overview = app.project_overview().await.unwrap();
    assert_eq!(overview.len(), 1);
    assert_eq!(overview[0].project, p);
    assert_eq!(overview[0].sample_count, 2);

    // chain a run + alignment off the first sample
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b1.guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: Some("PAIRED".into()),
            total_reads: Some(8_000_000),
            pf_reads_aligned: Some(7_956_881),
            mean_read_length: Some(148.0),
            mean_insert_size: Some(580.7),
        })
        .await
        .unwrap();
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
        })
        .await
        .unwrap();
    assert_eq!(aln.sequence_run_id, run.id);
}

#[tokio::test]
async fn add_biosample_to_missing_project_is_not_found() {
    let app = app().await;
    let err = app.add_biosample(Some(123), "HG002", None, None).await;
    assert!(matches!(err, Err(AppError::Store(navigator_store::StoreError::NotFound(_)))), "got {err:?}");
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CoverageSummary {
    mean_coverage: f64,
    callable_bases: u64,
}

#[tokio::test]
async fn typed_analysis_artifact_round_trips_and_versions() {
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
        })
        .await
        .unwrap();

    let summary = CoverageSummary { mean_coverage: 178.81, callable_bases: 16_292 };
    app.save_analysis(aln.id, "coverage", "walker-v1", &summary).await.unwrap();

    let loaded: CoverageSummary = app.load_analysis(aln.id, "coverage", "walker-v1").await.unwrap().unwrap();
    assert_eq!(loaded, summary);

    // a different version is absent until written
    let v2: Option<CoverageSummary> = app.load_analysis(aln.id, "coverage", "walker-v2").await.unwrap();
    assert!(v2.is_none());
}
