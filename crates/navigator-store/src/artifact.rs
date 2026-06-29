//! Analysis-artifact queries — a versioned result cache keyed by
//! `(alignment_id, kind, algorithm_version)`. `upsert` replaces a stale entry so a
//! changed algorithm version supersedes the old payload (plan §6 cache versioning).

use chrono::{DateTime, Utc};
use du_domain::ids::SampleGuid;
use navigator_domain::workspace::AnalysisArtifact;
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    alignment_id: i64,
    kind: String,
    algorithm_version: String,
    created_at: String,
    payload: String,
    source: Option<String>,
    completeness: Option<String>,
    source_sig: Option<String>,
}

impl Row {
    fn into_domain(self) -> Result<AnalysisArtifact, StoreError> {
        let created_at = DateTime::parse_from_rfc3339(&self.created_at)
            .map_err(|e| StoreError::Decode(format!("artifact created_at: {e}")))?
            .with_timezone(&Utc);
        Ok(AnalysisArtifact {
            id: self.id,
            alignment_id: self.alignment_id,
            kind: self.kind,
            algorithm_version: self.algorithm_version,
            created_at,
            payload: self.payload,
            source: self.source,
            completeness: self.completeness,
            source_sig: self.source_sig,
        })
    }
}

const COLS: &str = "id, alignment_id, kind, algorithm_version, created_at, payload, source, completeness, source_sig";

/// Insert or replace the artifact for `(alignment_id, kind, algorithm_version)`, recording its
/// provenance (`source` = how produced, `completeness` = full/partial).
#[allow(clippy::too_many_arguments)] // one parameter per artifact column — a DB row, not a refactor target
pub async fn upsert(
    pool: &SqlitePool,
    alignment_id: i64,
    kind: &str,
    algorithm_version: &str,
    created_at: DateTime<Utc>,
    payload: &str,
    source: &str,
    completeness: &str,
    source_sig: Option<&str>,
) -> Result<AnalysisArtifact, StoreError> {
    let created = created_at.to_rfc3339();
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO analysis_artifact (alignment_id, kind, algorithm_version, created_at, payload, source, completeness, source_sig) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT (alignment_id, kind, algorithm_version) \
         DO UPDATE SET created_at = excluded.created_at, payload = excluded.payload, \
                       source = excluded.source, completeness = excluded.completeness, \
                       source_sig = excluded.source_sig \
         RETURNING id",
    )
    .bind(alignment_id)
    .bind(kind)
    .bind(algorithm_version)
    .bind(&created)
    .bind(payload)
    .bind(source)
    .bind(completeness)
    .bind(source_sig)
    .fetch_one(pool)
    .await?;
    Ok(AnalysisArtifact {
        id,
        alignment_id,
        kind: kind.to_string(),
        algorithm_version: algorithm_version.to_string(),
        created_at,
        payload: payload.to_string(),
        source: Some(source.to_string()),
        completeness: Some(completeness.to_string()),
        source_sig: source_sig.map(str::to_string),
    })
}

pub async fn get(
    pool: &SqlitePool,
    alignment_id: i64,
    kind: &str,
    algorithm_version: &str,
) -> Result<Option<AnalysisArtifact>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM analysis_artifact \
         WHERE alignment_id = ? AND kind = ? AND algorithm_version = ?"
    ))
    .bind(alignment_id)
    .bind(kind)
    .bind(algorithm_version)
    .fetch_optional(pool)
    .await?;
    row.map(Row::into_domain).transpose()
}

pub async fn list_for_alignment(pool: &SqlitePool, alignment_id: i64) -> Result<Vec<AnalysisArtifact>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM analysis_artifact WHERE alignment_id = ? ORDER BY id"
    ))
    .bind(alignment_id)
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

/// Per-subject analysis coverage census, in one pass over the whole workspace: for each biosample
/// that owns ≥1 alignment, `(total alignments, alignments with a present `(kind, version)` artifact)`.
/// A NULL `completeness` counts as complete (legacy rows predate the column; the app treats absent
/// provenance as a full walk). Drives the Subjects-list Pending/Complete column — subjects with no
/// alignments are simply absent from the result. `kind`/`version` are passed in so the store stays
/// independent of the analysis crate (e.g. `"coverage"` / `coverage::COVERAGE_VERSION`).
pub async fn analyzed_census(
    pool: &SqlitePool,
    kind: &str,
    version: &str,
) -> Result<Vec<(SampleGuid, i64, i64)>, StoreError> {
    let rows: Vec<(String, i64, i64)> = sqlx::query_as(
        "SELECT sr.biosample_guid, \
                COUNT(DISTINCT a.id) AS total, \
                COUNT(DISTINCT CASE WHEN aa.alignment_id IS NOT NULL THEN a.id END) AS analyzed \
         FROM sequence_run sr \
         JOIN alignment a ON a.sequence_run_id = sr.id \
         LEFT JOIN analysis_artifact aa \
              ON aa.alignment_id = a.id \
             AND aa.kind = ? \
             AND aa.algorithm_version = ? \
             AND (aa.completeness = 'full' OR aa.completeness IS NULL) \
         GROUP BY sr.biosample_guid",
    )
    .bind(kind)
    .bind(version)
    .fetch_all(pool)
    .await?;
    rows.into_iter()
        .map(|(g, total, analyzed)| Ok((parse_sample_guid(&g, "analysis_artifact")?, total, analyzed)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::Utc;
    use navigator_domain::workspace::{Biosample, NewAlignment, NewSequenceRun};

    async fn subject(pool: &SqlitePool, donor: &str) -> SampleGuid {
        let guid = SampleGuid(uuid::Uuid::new_v4());
        crate::biosample::create(
            pool,
            &Biosample {
                guid,
                sample_accession: None,
                donor_identifier: donor.into(),
                description: None,
                center_name: None,
                sex: None,
                project_id: None,
            },
        )
        .await
        .unwrap();
        guid
    }

    async fn alignment(pool: &SqlitePool, guid: SampleGuid) -> i64 {
        let run = crate::sequence_run::create(
            pool,
            &NewSequenceRun {
                biosample_guid: guid,
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
        crate::alignment::create(
            pool,
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
        .unwrap()
        .id
    }

    async fn full_coverage(pool: &SqlitePool, aln: i64) {
        upsert(pool, aln, "coverage", "coverage-1", Utc::now(), "{}", "navigator-walk", "full", None)
            .await
            .unwrap();
    }

    #[tokio::test]
    async fn census_counts_full_coverage_per_subject() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();

        // A: one alignment, fully analyzed → Complete.
        let a = subject(pool, "A").await;
        let a_aln = alignment(pool, a).await;
        full_coverage(pool, a_aln).await;

        // B: two alignments, only one analyzed → Pending.
        let b = subject(pool, "B").await;
        let b1 = alignment(pool, b).await;
        let _b2 = alignment(pool, b).await;
        full_coverage(pool, b1).await;

        // C: one alignment with only a *partial* (sidecar) coverage → does not count → Pending.
        let c = subject(pool, "C").await;
        let c_aln = alignment(pool, c).await;
        upsert(pool, c_aln, "coverage", "coverage-1", Utc::now(), "{}", "pipeline-sidecar", "partial", None)
            .await
            .unwrap();

        // D: a subject with no alignments → absent from the census.
        let _d = subject(pool, "D").await;

        let census: std::collections::HashMap<_, _> = analyzed_census(pool, "coverage", "coverage-1")
            .await
            .unwrap()
            .into_iter()
            .map(|(g, total, analyzed)| (g, (total, analyzed)))
            .collect();

        assert_eq!(census.get(&a), Some(&(1, 1)), "A complete");
        assert_eq!(census.get(&b), Some(&(2, 1)), "B partially analyzed");
        assert_eq!(census.get(&c), Some(&(1, 0)), "partial coverage does not count");
        assert!(!census.contains_key(&_d), "no-alignment subject is absent");
    }
}
