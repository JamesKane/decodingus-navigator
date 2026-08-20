//! Analysis-artifact queries: a result cache with a version, keyed by `(alignment_id, kind,
//! algorithm_version)`. `upsert` replaces a stale entry, so a new algorithm version takes the place
//! of the old payload (plan §6, the cache version).

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

/// Insert or replace the artifact for `(alignment_id, kind, algorithm_version)`, and record its
/// provenance. `source` says how the app made it, and `completeness` says full or partial.
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

/// Remove one `(alignment, kind, version)` artifact when it is there, and do nothing when it is
/// not. The app uses it to clear a transient marker, such as the `error` artifact, once a later run
/// succeeds.
pub async fn delete(
    pool: &SqlitePool,
    alignment_id: i64,
    kind: &str,
    algorithm_version: &str,
) -> Result<(), StoreError> {
    sqlx::query("DELETE FROM analysis_artifact WHERE alignment_id = ? AND kind = ? AND algorithm_version = ?")
        .bind(alignment_id)
        .bind(kind)
        .bind(algorithm_version)
        .execute(pool)
        .await?;
    Ok(())
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

/// The `kind` of every artifact on `alignment_id`, without their payloads.
///
/// Exists so a caller can key a cache on *which* analyses are present without reading them. The
/// genome-wide de-novo calls run to ~1 GB of JSON across 22 contigs, so [`list_for_alignment`] is
/// the wrong tool for a question about coverage.
pub async fn list_kinds(pool: &SqlitePool, alignment_id: i64) -> Result<Vec<String>, StoreError> {
    let rows: Vec<(String,)> =
        sqlx::query_as("SELECT kind FROM analysis_artifact WHERE alignment_id = ? ORDER BY kind")
            .bind(alignment_id)
            .fetch_all(pool)
            .await?;
    Ok(rows.into_iter().map(|(k,)| k).collect())
}

/// Every artifact of every id in `alignment_ids`, in one query. The caller indexes the result by
/// `(alignment_id, kind)` itself. This replaces a `get` for each (alignment, kind) pair, which gave
/// a project report one round trip for each cell. An empty `alignment_ids` runs no query.
pub async fn list_for_alignments(
    pool: &SqlitePool,
    alignment_ids: &[i64],
) -> Result<Vec<AnalysisArtifact>, StoreError> {
    if alignment_ids.is_empty() {
        return Ok(Vec::new());
    }
    // SQLite can not bind an array. The id count gives the length of the placeholder list, and no
    // user text goes into it. The code still binds every id, so this SQL holds no interpolated
    // string.
    let placeholders = vec!["?"; alignment_ids.len()].join(",");
    let sql = format!("SELECT {COLS} FROM analysis_artifact WHERE alignment_id IN ({placeholders}) ORDER BY id");
    let mut q = sqlx::query_as(&sql);
    for id in alignment_ids {
        q = q.bind(id);
    }
    let rows: Vec<Row> = q.fetch_all(pool).await?;
    rows.into_iter().map(Row::into_domain).collect()
}

/// [`list_for_alignments`] narrowed to a single `(kind, version)`.
///
/// Use this whenever the caller wants one artifact kind. The query with no filter selects `payload`
/// for *every* artifact of every listed alignment. Some kinds are very large: a `tree-genotype` row
/// runs to megabytes. To fetch a whole cohort and then pick out one small kind reads gigabytes of
/// JSON for nothing.
pub async fn list_for_alignments_of_kind(
    pool: &SqlitePool,
    alignment_ids: &[i64],
    kind: &str,
    version: &str,
) -> Result<Vec<AnalysisArtifact>, StoreError> {
    if alignment_ids.is_empty() {
        return Ok(Vec::new());
    }
    // As in `list_for_alignments`: the id count gives the placeholders, no interpolated text goes
    // into them, and the code binds every value.
    let placeholders = vec!["?"; alignment_ids.len()].join(",");
    let sql = format!(
        "SELECT {COLS} FROM analysis_artifact \
         WHERE kind = ? AND algorithm_version = ? AND alignment_id IN ({placeholders}) ORDER BY id"
    );
    let mut q = sqlx::query_as(&sql).bind(kind).bind(version);
    for id in alignment_ids {
        q = q.bind(id);
    }
    let rows: Vec<Row> = q.fetch_all(pool).await?;
    rows.into_iter().map(Row::into_domain).collect()
}

/// A census of the analysis coverage of each subject, in one pass over the whole workspace. For
/// each biosample that owns one alignment or more it gives `(total alignments, alignments that have
/// a `(kind, version)` artifact)`.
///
/// A NULL `completeness` counts as complete. A legacy row comes from before that column existed,
/// and the app reads absent provenance as a full walk.
///
/// It feeds the Pending and Complete column of the Subjects list. A subject with no alignments is
/// absent from the result. The caller passes `kind` and `version` in, so this crate needs nothing
/// from the analysis crate. An example pair is `"coverage"` with `coverage::COVERAGE_VERSION`.
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
        crate::biosample::create(pool, &Biosample::new(guid, donor))
            .await
            .unwrap();
        guid
    }

    async fn alignment(pool: &SqlitePool, guid: SampleGuid) -> i64 {
        let run = crate::sequence_run::create(pool, &NewSequenceRun::new(guid, "ILLUMINA", "WGS"))
            .await
            .unwrap();
        crate::alignment::create(pool, &NewAlignment::new(run.id, "chm13v2.0", "bwa"))
            .await
            .unwrap()
            .id
    }

    async fn full_coverage(pool: &SqlitePool, aln: i64) {
        upsert(
            pool,
            aln,
            "coverage",
            "coverage-1",
            Utc::now(),
            "{}",
            "navigator-walk",
            "full",
            None,
        )
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
        upsert(
            pool,
            c_aln,
            "coverage",
            "coverage-1",
            Utc::now(),
            "{}",
            "pipeline-sidecar",
            "partial",
            None,
        )
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
