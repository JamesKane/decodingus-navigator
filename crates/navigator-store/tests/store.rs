//! Store integration tests against an in-memory SQLite database.

use chrono::Utc;
use du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Biosample, NewAlignment, NewProject, NewSequenceRun};
use navigator_store::{alignment, artifact, biosample, project, sequence_run, source_file, Store};
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
        &NewProject {
            name: "Trio".into(),
            description: Some("HG002 trio".into()),
            administrator: "jk".into(),
        },
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
    let p = project::create(
        s.pool(),
        &NewProject {
            name: "P".into(),
            description: None,
            administrator: "jk".into(),
        },
    )
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
    assert_eq!(
        sequence_run::list_for_biosample(s.pool(), b.guid).await.unwrap()[0],
        run
    );
    assert_eq!(run.mean_insert_size, Some(580.7)); // flat metric column round-trips
                                                   // The lab/instrument identity block is None at create, then filled by set_library_stats.
    assert_eq!(run.instrument_id, None);
    sequence_run::set_library_stats(
        s.pool(),
        run.id,
        Some("A00182"),
        Some("WGS229"),
        Some("WGS229_Lib1"),
        Some("CeGaT_NovaSeq"),
        Some("H5WLTDMXX"),
    )
    .await
    .unwrap();
    let reloaded = sequence_run::get(s.pool(), run.id).await.unwrap().unwrap();
    assert_eq!(reloaded.instrument_id.as_deref(), Some("A00182"));
    assert_eq!(reloaded.sample_name.as_deref(), Some("WGS229"));
    assert_eq!(reloaded.platform_unit.as_deref(), Some("CeGaT_NovaSeq"));
    assert_eq!(reloaded.flowcell_id.as_deref(), Some("H5WLTDMXX"));

    // Library-level read stats can be (re)written from an analysis pass / backfill.
    sequence_run::set_read_stats(
        s.pool(),
        run.id,
        Some(9_000_000),
        Some(150.0),
        Some(602.5),
        Some("PAIRED"),
    )
    .await
    .unwrap();
    let reloaded = sequence_run::get(s.pool(), run.id).await.unwrap().unwrap();
    assert_eq!(reloaded.total_reads, Some(9_000_000));
    assert_eq!(reloaded.mean_read_length, Some(150.0));
    assert_eq!(reloaded.mean_insert_size, Some(602.5));
    assert_eq!(reloaded.library_layout.as_deref(), Some("PAIRED"));
    // A `None` layout preserves the stored value (COALESCE), other columns still update.
    sequence_run::set_read_stats(s.pool(), run.id, Some(9_100_000), Some(150.0), Some(602.5), None)
        .await
        .unwrap();
    let reloaded = sequence_run::get(s.pool(), run.id).await.unwrap().unwrap();
    assert_eq!(reloaded.total_reads, Some(9_100_000));
    assert_eq!(reloaded.library_layout.as_deref(), Some("PAIRED"));
    // The descriptive + identity columns are untouched by the read-stats write.
    assert_eq!(reloaded.instrument_id.as_deref(), Some("A00182"));

    let aln = alignment::create(
        s.pool(),
        &NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa-mem 0.7.19".into(),
            variant_caller: Some("navigator-haploid".into()),
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        },
    )
    .await
    .unwrap();
    assert_eq!(alignment::list_for_run(s.pool(), run.id).await.unwrap(), vec![aln]);
}

#[tokio::test]
async fn str_find_by_panel_and_replace_markers() {
    use navigator_domain::strprofile::{NewStrProfile, StrMarker};
    use navigator_store::str_profile;

    let s = store().await;
    let b = sample(None);
    biosample::create(s.pool(), &b).await.unwrap();
    let mk = |m: &str, v: &str| StrMarker {
        marker: m.into(),
        value: v.into(),
    };

    str_profile::create(
        s.pool(),
        &NewStrProfile {
            biosample_guid: b.guid,
            panel_name: "CUSTOM".into(),
            provider: None,
            source: Some("IMPORTED".into()),
            markers: vec![mk("DYS393", "13"), mk("DYS390", "24")],
        },
    )
    .await
    .unwrap();

    // Find by panel resolves the profile with its markers; a different panel is absent.
    let found = str_profile::find_by_panel(s.pool(), b.guid, "CUSTOM")
        .await
        .unwrap()
        .unwrap();
    assert_eq!(found.markers.len(), 2);
    assert!(str_profile::find_by_panel(s.pool(), b.guid, "Y-111")
        .await
        .unwrap()
        .is_none());

    // Replace markers (merge result: one updated value + one new marker).
    str_profile::replace_markers(s.pool(), found.id, &[mk("DYS393", "14"), mk("DYS19", "15")])
        .await
        .unwrap();
    let after = str_profile::find_by_panel(s.pool(), b.guid, "CUSTOM")
        .await
        .unwrap()
        .unwrap();
    assert_eq!(after.id, found.id, "same profile (no duplicate)");
    assert_eq!(after.markers.len(), 2);
    assert_eq!(after.markers.iter().find(|m| m.marker == "DYS393").unwrap().value, "14");
}

#[tokio::test]
async fn clear_data_resets_subject_but_keeps_the_biosample() {
    use navigator_domain::reconciliation::{AuditEntry, DnaType};
    use navigator_store::reconciliation;

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
        &NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        },
    )
    .await
    .unwrap();
    artifact::upsert(
        s.pool(),
        aln.id,
        "coverage",
        "v1",
        Utc::now(),
        "{}",
        "walk",
        "full",
        None,
    )
    .await
    .unwrap();
    // A biosample-keyed derived row (the kind that orphaned subject 103589).
    reconciliation::append_audit(
        s.pool(),
        b.guid,
        DnaType::Y,
        &AuditEntry {
            timestamp: "2026-06-20T00:00:00Z".into(),
            action: "RUN_RECORDED".into(),
            note: "aln:1".into(),
        },
    )
    .await
    .unwrap();

    biosample::clear_data(s.pool(), b.guid).await.unwrap();

    // The subject survives; everything hanging off it is gone.
    assert!(
        biosample::get(s.pool(), b.guid).await.unwrap().is_some(),
        "biosample kept"
    );
    assert!(
        sequence_run::list_for_biosample(s.pool(), b.guid)
            .await
            .unwrap()
            .is_empty(),
        "runs cleared"
    );
    assert!(
        alignment::list_for_run(s.pool(), run.id).await.unwrap().is_empty(),
        "alignments cleared"
    );
    assert!(
        reconciliation::list_audit(s.pool(), b.guid, DnaType::Y)
            .await
            .unwrap()
            .is_empty(),
        "audit cleared"
    );

    // Idempotent: a second clear is a no-op (no error).
    biosample::clear_data(s.pool(), b.guid).await.unwrap();
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
        &NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        },
    )
    .await
    .unwrap();

    // Same (kind, version) upserts in place.
    artifact::upsert(
        s.pool(),
        aln.id,
        "coverage",
        "v1",
        Utc::now(),
        r#"{"mean":1.0}"#,
        "navigator-walk",
        "full",
        None,
    )
    .await
    .unwrap();
    let updated = artifact::upsert(
        s.pool(),
        aln.id,
        "coverage",
        "v1",
        Utc::now(),
        r#"{"mean":2.0}"#,
        "navigator-walk",
        "full",
        None,
    )
    .await
    .unwrap();
    let got = artifact::get(s.pool(), aln.id, "coverage", "v1")
        .await
        .unwrap()
        .unwrap();
    assert_eq!(got.id, updated.id);
    assert_eq!(got.payload, r#"{"mean":2.0}"#);

    // A new algorithm version is a distinct entry.
    artifact::upsert(
        s.pool(),
        aln.id,
        "coverage",
        "v2",
        Utc::now(),
        r#"{"mean":3.0}"#,
        "navigator-walk",
        "full",
        None,
    )
    .await
    .unwrap();
    assert_eq!(artifact::list_for_alignment(s.pool(), aln.id).await.unwrap().len(), 2);
}

#[tokio::test]
async fn delete_cascades_run_to_alignments_and_artifacts() {
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
        &NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        },
    )
    .await
    .unwrap();
    artifact::upsert(
        s.pool(),
        aln.id,
        "coverage",
        "v1",
        Utc::now(),
        r#"{"mean":1.0}"#,
        "navigator-walk",
        "full",
        None,
    )
    .await
    .unwrap();

    // Deleting a single alignment removes its artifacts but leaves the run.
    let aln2 = alignment::create(
        s.pool(),
        &NewAlignment {
            sequence_run_id: run.id,
            reference_build: "grch38".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        },
    )
    .await
    .unwrap();
    artifact::upsert(
        s.pool(),
        aln2.id,
        "coverage",
        "v1",
        Utc::now(),
        r#"{"mean":2.0}"#,
        "navigator-walk",
        "full",
        None,
    )
    .await
    .unwrap();
    // A content-hash source_file linked to aln2 must not block its delete (regression: FK 787).
    source_file::upsert_by_checksum(s.pool(), "hash-aln2", Some("/data/x.bam"), Some(1), Some("BAM"), "t")
        .await
        .unwrap();
    source_file::link_to_alignment(s.pool(), "hash-aln2", aln2.id, "t")
        .await
        .unwrap();
    assert!(alignment::delete(s.pool(), aln2.id).await.unwrap());
    assert!(artifact::get(s.pool(), aln2.id, "coverage", "v1")
        .await
        .unwrap()
        .is_none());
    // The file record survives but is unlinked from the deleted alignment.
    let sf = source_file::find_by_checksum(s.pool(), "hash-aln2")
        .await
        .unwrap()
        .unwrap();
    assert_eq!(sf.alignment_id, None, "source_file unlinked, not deleted");
    assert_eq!(
        alignment::list_for_run(s.pool(), run.id).await.unwrap(),
        vec![aln.clone()]
    );
    assert!(sequence_run::get(s.pool(), run.id).await.unwrap().is_some());

    // Deleting the run removes the remaining alignment + artifact (FK-enforced cascade), and a
    // source_file linked to that alignment must not block it either.
    source_file::upsert_by_checksum(s.pool(), "hash-aln1", Some("/data/y.bam"), Some(1), Some("BAM"), "t")
        .await
        .unwrap();
    source_file::link_to_alignment(s.pool(), "hash-aln1", aln.id, "t")
        .await
        .unwrap();
    assert!(sequence_run::delete(s.pool(), run.id).await.unwrap());
    assert_eq!(
        source_file::find_by_checksum(s.pool(), "hash-aln1")
            .await
            .unwrap()
            .unwrap()
            .alignment_id,
        None,
        "source_file unlinked when its alignment's run is deleted"
    );
    assert!(sequence_run::get(s.pool(), run.id).await.unwrap().is_none());
    assert!(alignment::list_for_run(s.pool(), run.id).await.unwrap().is_empty());
    assert!(artifact::get(s.pool(), aln.id, "coverage", "v1")
        .await
        .unwrap()
        .is_none());

    // Deleting a non-existent row reports false rather than erroring.
    assert!(!sequence_run::delete(s.pool(), 9999).await.unwrap());
    assert!(!alignment::delete(s.pool(), 9999).await.unwrap());
}

#[tokio::test]
async fn set_sequence_run_reparents_an_alignment() {
    // The merge primitive: an alignment's owning run can be changed (then the empty run deleted),
    // and its artifacts travel with it (they're alignment-keyed).
    let s = store().await;
    let b = sample(None);
    biosample::create(s.pool(), &b).await.unwrap();
    let mk_run = |layout: &str| NewSequenceRun {
        biosample_guid: b.guid,
        platform_name: "ILLUMINA".into(),
        instrument_model: None,
        test_type: "WGS".into(),
        library_layout: Some(layout.to_string()),
        total_reads: None,
        pf_reads_aligned: None,
        mean_read_length: None,
        mean_insert_size: None,
    };
    let primary = sequence_run::create(s.pool(), &mk_run("A")).await.unwrap();
    let secondary = sequence_run::create(s.pool(), &mk_run("B")).await.unwrap();
    let aln = alignment::create(
        s.pool(),
        &NewAlignment {
            sequence_run_id: secondary.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        },
    )
    .await
    .unwrap();
    artifact::upsert(
        s.pool(),
        aln.id,
        "coverage",
        "v1",
        Utc::now(),
        r#"{"mean":3.0}"#,
        "navigator-walk",
        "full",
        None,
    )
    .await
    .unwrap();

    // Reparent the alignment onto the primary run.
    assert!(alignment::set_sequence_run(s.pool(), aln.id, primary.id).await.unwrap());
    assert!(alignment::list_for_run(s.pool(), secondary.id)
        .await
        .unwrap()
        .is_empty());
    assert_eq!(alignment::list_for_run(s.pool(), primary.id).await.unwrap().len(), 1);

    // Deleting the now-empty secondary leaves the moved alignment + its artifact intact under primary.
    assert!(sequence_run::delete(s.pool(), secondary.id).await.unwrap());
    assert!(alignment::get(s.pool(), aln.id).await.unwrap().is_some());
    assert!(artifact::get(s.pool(), aln.id, "coverage", "v1")
        .await
        .unwrap()
        .is_some());

    // Reparenting a non-existent alignment reports false.
    assert!(!alignment::set_sequence_run(s.pool(), 9999, primary.id).await.unwrap());
}

#[tokio::test]
async fn haplogroup_call_fingerprint_round_trips() {
    use navigator_domain::reconciliation::{DnaType, RunHaplogroupCall};
    use navigator_store::haplogroup_call;

    let s = store().await;
    let b = sample(None);
    biosample::create(s.pool(), &b).await.unwrap();

    let call = RunHaplogroupCall {
        source_label: "aln #1 Y".into(),
        haplogroup: "R-FGC29071".into(),
        lineage: vec!["Y".into(), "R".into(), "R-FGC29071".into()],
        score: 0.9,
        matched: 80,
        expected: 100,
    };
    // Upsert stamps the fingerprint; stored_fingerprint reads it back.
    haplogroup_call::upsert(s.pool(), b.guid, DnaType::Y, "aln:1", &call, Some("f:abc|yt:def"))
        .await
        .unwrap();
    assert_eq!(
        haplogroup_call::stored_fingerprint(s.pool(), b.guid, DnaType::Y, "aln:1")
            .await
            .unwrap()
            .as_deref(),
        Some("f:abc|yt:def")
    );
    let got = haplogroup_call::get_one(s.pool(), b.guid, DnaType::Y, "aln:1")
        .await
        .unwrap()
        .unwrap();
    assert_eq!(got.haplogroup, "R-FGC29071");
    assert_eq!(
        got.lineage,
        vec!["Y".to_string(), "R".to_string(), "R-FGC29071".to_string()]
    );

    // Re-upsert with a new fingerprint (e.g. the tree changed) replaces it.
    haplogroup_call::upsert(s.pool(), b.guid, DnaType::Y, "aln:1", &call, Some("f:abc|yt:NEW"))
        .await
        .unwrap();
    assert_eq!(
        haplogroup_call::stored_fingerprint(s.pool(), b.guid, DnaType::Y, "aln:1")
            .await
            .unwrap()
            .as_deref(),
        Some("f:abc|yt:NEW")
    );

    // Unknown source → no fingerprint / no call.
    assert!(
        haplogroup_call::stored_fingerprint(s.pool(), b.guid, DnaType::Y, "nope")
            .await
            .unwrap()
            .is_none()
    );
    assert!(haplogroup_call::get_one(s.pool(), b.guid, DnaType::Y, "nope")
        .await
        .unwrap()
        .is_none());
}
