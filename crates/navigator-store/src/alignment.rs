//! Alignment queries.

use du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Alignment, NewAlignment};
use sqlx::SqlitePool;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    sequence_run_id: i64,
    reference_build: String,
    aligner: String,
    variant_caller: Option<String>,
    bam_path: Option<String>,
    reference_path: Option<String>,
    content_sha256: Option<String>,
}

impl Row {
    fn into_domain(self) -> Alignment {
        Alignment {
            id: self.id,
            sequence_run_id: self.sequence_run_id,
            reference_build: self.reference_build,
            aligner: self.aligner,
            variant_caller: self.variant_caller,
            bam_path: self.bam_path,
            reference_path: self.reference_path,
            content_sha256: self.content_sha256,
        }
    }
}

const COLS: &str =
    "id, sequence_run_id, reference_build, aligner, variant_caller, bam_path, reference_path, content_sha256";

pub async fn create(pool: &SqlitePool, a: &NewAlignment) -> Result<Alignment, StoreError> {
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO alignment (sequence_run_id, reference_build, aligner, variant_caller, bam_path, reference_path, content_sha256) \
         VALUES (?, ?, ?, ?, ?, ?, ?) RETURNING id",
    )
    .bind(a.sequence_run_id)
    .bind(&a.reference_build)
    .bind(&a.aligner)
    .bind(&a.variant_caller)
    .bind(&a.bam_path)
    .bind(&a.reference_path)
    .bind(&a.content_sha256)
    .fetch_one(pool)
    .await?;
    Ok(Alignment {
        id,
        sequence_run_id: a.sequence_run_id,
        reference_build: a.reference_build.clone(),
        aligner: a.aligner.clone(),
        variant_caller: a.variant_caller.clone(),
        bam_path: a.bam_path.clone(),
        reference_path: a.reference_path.clone(),
        content_sha256: a.content_sha256.clone(),
    })
}

/// Set (or update) an alignment's content SHA-256, computed lazily after import.
pub async fn set_content_hash(pool: &SqlitePool, id: i64, sha256: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE alignment SET content_sha256 = ? WHERE id = ?")
        .bind(sha256)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Reparent an alignment to a different sequence run (used by run merge). Just an `UPDATE` of the
/// FK column — no schema change. Returns whether a row was affected.
pub async fn set_sequence_run(pool: &SqlitePool, id: i64, sequence_run_id: i64) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE alignment SET sequence_run_id = ? WHERE id = ?")
        .bind(sequence_run_id)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<Alignment>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM alignment WHERE id = ?"))
        .bind(id)
        .fetch_optional(pool)
        .await?;
    Ok(row.map(Row::into_domain))
}

pub async fn list_for_run(pool: &SqlitePool, sequence_run_id: i64) -> Result<Vec<Alignment>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM alignment WHERE sequence_run_id = ? ORDER BY id"
    ))
    .bind(sequence_run_id)
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(Row::into_domain).collect())
}

/// Every alignment in the workspace (for cross-sample selection, e.g. IBD compare).
pub async fn list_all(pool: &SqlitePool) -> Result<Vec<Alignment>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM alignment ORDER BY id"))
        .fetch_all(pool)
        .await?;
    Ok(rows.into_iter().map(Row::into_domain).collect())
}

/// Update an alignment's descriptive fields (reference build, aligner, variant caller). The
/// file paths are managed by import/probe, not here. Returns whether a row was affected.
pub async fn update(
    pool: &SqlitePool,
    id: i64,
    reference_build: &str,
    aligner: &str,
    variant_caller: Option<&str>,
) -> Result<bool, StoreError> {
    let affected =
        sqlx::query("UPDATE alignment SET reference_build = ?, aligner = ?, variant_caller = ? WHERE id = ?")
            .bind(reference_build)
            .bind(aligner)
            .bind(variant_caller)
            .bind(id)
            .execute(pool)
            .await?
            .rows_affected();
    Ok(affected > 0)
}

/// Delete an alignment and its cached analysis artifacts (FKs are enforced, so the
/// `analysis_artifact` children go first). Returns whether the alignment row was removed.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let mut tx = pool.begin().await?;
    sqlx::query("DELETE FROM analysis_artifact WHERE alignment_id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?;
    let affected = sqlx::query("DELETE FROM alignment WHERE id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?
        .rows_affected();
    tx.commit().await?;
    Ok(affected > 0)
}

/// All alignments for a biosample (joined through its sequence runs).
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<Alignment>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {} FROM alignment a JOIN sequence_run r ON a.sequence_run_id = r.id \
         WHERE r.biosample_guid = ? ORDER BY a.id",
        COLS.split(", ")
            .map(|c| format!("a.{c}"))
            .collect::<Vec<_>>()
            .join(", ")
    ))
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(Row::into_domain).collect())
}
