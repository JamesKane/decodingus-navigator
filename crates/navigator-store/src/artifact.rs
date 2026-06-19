//! Analysis-artifact queries — a versioned result cache keyed by
//! `(alignment_id, kind, algorithm_version)`. `upsert` replaces a stale entry so a
//! changed algorithm version supersedes the old payload (plan §6 cache versioning).

use chrono::{DateTime, Utc};
use navigator_domain::workspace::AnalysisArtifact;
use sqlx::SqlitePool;

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
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM analysis_artifact WHERE alignment_id = ? ORDER BY id"))
        .bind(alignment_id)
        .fetch_all(pool)
        .await?;
    rows.into_iter().map(Row::into_domain).collect()
}
