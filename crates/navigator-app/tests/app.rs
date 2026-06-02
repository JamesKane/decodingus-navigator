//! App command/query layer tests against an in-memory store.

use std::path::PathBuf;

use navigator_analysis::caller::HaploidCallerParams;
use navigator_analysis::coverage::CallableLociParams;
use navigator_app::{App, AppError};
use navigator_domain::workspace::{NewAlignment, NewProject, NewSequenceRun};
use navigator_store::Store;
use serde::{Deserialize, Serialize};

async fn app() -> App {
    App::new(Store::open_in_memory().await.unwrap())
}

/// Reuse the analysis crate's committed fixtures (workspace-relative).
fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../navigator-analysis/tests/fixtures")
}

/// Create a sample → run → alignment chain and return the alignment id.
async fn alignment_id(app: &App) -> i64 {
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
    app.record_alignment(NewAlignment {
        sequence_run_id: run.id,
        reference_build: "chrM-fixture".into(),
        aligner: "synthetic".into(),
        variant_caller: None, bam_path: None, reference_path: None,
    })
    .await
    .unwrap()
    .id
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
            variant_caller: None, bam_path: None, reference_path: None,
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
            variant_caller: None, bam_path: None, reference_path: None,
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

#[tokio::test]
async fn run_coverage_persists_and_reads_back_from_cache() {
    let app = app().await;
    let aln = alignment_id(&app).await;
    let dir = fixtures();

    assert!(app.cached_coverage(aln).await.unwrap().is_none()); // cold

    let result = app
        .run_coverage(aln, dir.join("coverage.bam"), dir.join("ref.fa"), None, CallableLociParams::default())
        .await
        .unwrap();
    assert_eq!(result.genome_territory, 50); // chrM fixture
    assert_eq!(result.callable_bases, 10);

    // now cached for this version (integer fields exact; floats survive round-trip to
    // ~1 ULP, so compare those approximately rather than with fragile float ==)
    let cached = app.cached_coverage(aln).await.unwrap().unwrap();
    assert_eq!(cached.genome_territory, result.genome_territory);
    assert_eq!(cached.callable_bases, result.callable_bases);
    assert_eq!(cached.contig_callable, result.contig_callable);
    assert_eq!(cached.coverage_histogram, result.coverage_histogram);
    assert!((cached.mean_coverage - result.mean_coverage).abs() < 1e-9);

    // re-running is idempotent (upsert in place; store-layer test covers row count)
    let rerun = app
        .run_coverage(aln, dir.join("coverage.bam"), dir.join("ref.fa"), None, CallableLociParams::default())
        .await
        .unwrap();
    assert_eq!(rerun, result);
}

#[tokio::test]
async fn run_denovo_caller_persists_snp_calls() {
    let app = app().await;
    let aln = alignment_id(&app).await;
    let dir = fixtures();

    let calls = app
        .run_denovo_caller(aln, dir.join("coverage.bam"), dir.join("ref.fa"), "chrM".into(), HaploidCallerParams::default())
        .await
        .unwrap();
    // fixture: ref ACGT.. with all-A reads -> SNPs where ref != A at depth>=4
    assert_eq!(calls.iter().map(|c| c.position).collect::<Vec<_>>(), vec![2, 3, 4, 6, 7, 8, 10]);

    let cached = app.cached_denovo(aln, "chrM").await.unwrap().unwrap();
    assert_eq!(cached, calls);
    // a different contig has no cached result
    assert!(app.cached_denovo(aln, "chrY").await.unwrap().is_none());
}

/// Build a sample → run → alignment chain whose BAM is the diploid fixture, return id.
async fn diploid_alignment(app: &App) -> i64 {
    let b = app.add_biosample(None, "diploid", None, None).await.unwrap();
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
    let bam = fixtures().join("diploid.bam").to_string_lossy().into_owned();
    app.record_alignment(NewAlignment {
        sequence_run_id: run.id,
        reference_build: "chr1".into(),
        aligner: "synthetic".into(),
        variant_caller: None,
        bam_path: Some(bam),
        reference_path: None,
    })
    .await
    .unwrap()
    .id
}

#[tokio::test]
async fn panel_genotyping_then_ibd_compare() {
    use navigator_analysis::ibd::IbdDetectorConfig;
    use navigator_domain::workspace::PanelSite;

    let app = app().await;
    let site = |pos, r: &str, a: &str| PanelSite {
        chrom: "chr1".into(),
        position: pos,
        reference_allele: r.into(),
        alternate_allele: a.into(),
        name: format!("s{pos}"),
    };
    // The four informative sites in the diploid fixture: hom-ref, het, het, hom-alt.
    let sites = vec![site(1, "A", "G"), site(2, "C", "G"), site(5, "A", "T"), site(8, "T", "A")];

    let panel = app.import_panel("test-panel", &sites).await.unwrap();
    assert_eq!(app.panel_site_count(panel.id).await.unwrap(), 4);

    let aln = diploid_alignment(&app).await;
    let genos = app.genotype_panel(aln, panel.id, 2).await.unwrap();
    let dosages: Vec<i32> = genos.iter().map(|g| g.dosage).collect();
    assert_eq!(dosages, vec![0, 1, 1, 2]);

    // cached read-back
    assert_eq!(app.cached_panel_genotypes(aln, panel.id, 2).await.unwrap().unwrap(), genos);

    // IBD self-compare with relaxed thresholds (only 4 sites): one fully-shared segment.
    let cfg = IbdDetectorConfig { min_snp_count: 3, window_size: 3, min_segment_cm: 0.0, ..IbdDetectorConfig::default() };
    let cmp = app.compare_ibd(aln, aln, panel.id, 2, cfg).await.unwrap();
    assert_eq!(cmp.segments.len(), 1);
    assert!(cmp.summary.total_shared_cm >= 0.0);

    // comparing against an un-genotyped alignment errors clearly.
    let other = diploid_alignment(&app).await;
    assert!(matches!(app.compare_ibd(aln, other, panel.id, 2, cfg).await, Err(AppError::NotGenotyped(_))));
}
