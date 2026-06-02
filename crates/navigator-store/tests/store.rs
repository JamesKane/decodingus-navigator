//! Store integration tests against an in-memory SQLite database.

use chrono::Utc;
use du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Biosample, NewAlignment, NewProject, NewSequenceRun};
use navigator_store::{alignment, artifact, biosample, project, sequence_run, Store};
use uuid::Uuid;

async fn store() -> Store {
    Store::open_in_memory().await.expect("open in-memory db")
}

fn sample(project_id: Option<i64>) -> Biosample {
    Biosample {
        guid: SampleGuid(Uuid::new_v4()),
        sample_accession: Some("SAMEA1".into()),
        donor_identifier: "HG002".into(),
        description: None,
        center_name: Some("NIST".into()),
        sex: Some("male".into()),
        project_id,
    }
}

#[tokio::test]
async fn project_round_trips() {
    let s = store().await;
    let p = project::create(
        s.pool(),
        &NewProject { name: "Trio".into(), description: Some("HG002 trio".into()), administrator: "jk".into() },
    )
    .await
    .unwrap();
    assert_eq!(p.id, 1);

    let got = project::get(s.pool(), p.id).await.unwrap().unwrap();
    assert_eq!(got, p);
    assert_eq!(project::list(s.pool()).await.unwrap().len(), 1);
    assert!(project::get(s.pool(), 999).await.unwrap().is_none());
}

#[tokio::test]
async fn biosample_links_to_project_and_round_trips() {
    let s = store().await;
    let p = project::create(s.pool(), &NewProject { name: "P".into(), description: None, administrator: "jk".into() })
        .await
        .unwrap();
    let b = sample(Some(p.id));
    biosample::create(s.pool(), &b).await.unwrap();

    let got = biosample::get(s.pool(), b.guid).await.unwrap().unwrap();
    assert_eq!(got, b);
    assert_eq!(biosample::count_for_project(s.pool(), p.id).await.unwrap(), 1);
    assert_eq!(biosample::list_for_project(s.pool(), p.id).await.unwrap(), vec![b]);
}

#[tokio::test]
async fn foreign_keys_are_enforced() {
    let s = store().await;
    // biosample referencing a non-existent project must fail.
    let err = biosample::create(s.pool(), &sample(Some(42))).await;
    assert!(err.is_err(), "expected FK violation, got {err:?}");

    // sequence_run referencing a non-existent biosample must fail.
    let run = NewSequenceRun {
        biosample_guid: SampleGuid(Uuid::new_v4()),
        platform_name: "ILLUMINA".into(),
        instrument_model: None,
        test_type: "WGS".into(),
        library_layout: Some("PAIRED".into()),
        total_reads: Some(8_000_000),
        pf_reads_aligned: Some(7_956_881),
        mean_read_length: Some(148.0),
        mean_insert_size: Some(580.7),
    };
    assert!(sequence_run::create(s.pool(), &run).await.is_err());
}

#[tokio::test]
async fn run_alignment_chain_persists() {
    let s = store().await;
    let b = sample(None);
    biosample::create(s.pool(), &b).await.unwrap();

    let run = sequence_run::create(
        s.pool(),
        &NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: Some("HiSeq 2500".into()),
            test_type: "WGS".into(),
            library_layout: Some("PAIRED".into()),
            total_reads: Some(8_000_000),
            pf_reads_aligned: Some(7_956_881),
            mean_read_length: Some(148.0),
            mean_insert_size: Some(580.7),
        },
    )
    .await
    .unwrap();
    assert_eq!(sequence_run::list_for_biosample(s.pool(), b.guid).await.unwrap()[0], run);
    assert_eq!(run.mean_insert_size, Some(580.7)); // flat metric column round-trips

    let aln = alignment::create(
        s.pool(),
        &NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa-mem 0.7.19".into(),
            variant_caller: Some("navigator-haploid".into()),
        },
    )
    .await
    .unwrap();
    assert_eq!(alignment::list_for_run(s.pool(), run.id).await.unwrap(), vec![aln]);
}

#[tokio::test]
async fn artifact_upsert_replaces_same_version_and_keeps_distinct_versions() {
    let s = store().await;
    let b = sample(None);
    biosample::create(s.pool(), &b).await.unwrap();
    let run = sequence_run::create(
        s.pool(),
        &NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        },
    )
    .await
    .unwrap();
    let aln = alignment::create(
        s.pool(),
        &NewAlignment { sequence_run_id: run.id, reference_build: "chm13v2.0".into(), aligner: "bwa".into(), variant_caller: None },
    )
    .await
    .unwrap();

    // Same (kind, version) upserts in place.
    artifact::upsert(s.pool(), aln.id, "coverage", "v1", Utc::now(), r#"{"mean":1.0}"#).await.unwrap();
    let updated = artifact::upsert(s.pool(), aln.id, "coverage", "v1", Utc::now(), r#"{"mean":2.0}"#).await.unwrap();
    let got = artifact::get(s.pool(), aln.id, "coverage", "v1").await.unwrap().unwrap();
    assert_eq!(got.id, updated.id);
    assert_eq!(got.payload, r#"{"mean":2.0}"#);

    // A new algorithm version is a distinct entry.
    artifact::upsert(s.pool(), aln.id, "coverage", "v2", Utc::now(), r#"{"mean":3.0}"#).await.unwrap();
    assert_eq!(artifact::list_for_alignment(s.pool(), aln.id).await.unwrap().len(), 2);
}
