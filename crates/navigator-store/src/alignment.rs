//! Alignment queries.

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
}

impl Row {
    fn into_domain(self) -> Alignment {
        Alignment {
            id: self.id,
            sequence_run_id: self.sequence_run_id,
            reference_build: self.reference_build,
            aligner: self.aligner,
            variant_caller: self.variant_caller,
        }
    }
}

const COLS: &str = "id, sequence_run_id, reference_build, aligner, variant_caller";

pub async fn create(pool: &SqlitePool, a: &NewAlignment) -> Result<Alignment, StoreError> {
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO alignment (sequence_run_id, reference_build, aligner, variant_caller) \
         VALUES (?, ?, ?, ?) RETURNING id",
    )
    .bind(a.sequence_run_id)
    .bind(&a.reference_build)
    .bind(&a.aligner)
    .bind(&a.variant_caller)
    .fetch_one(pool)
    .await?;
    Ok(Alignment {
        id,
        sequence_run_id: a.sequence_run_id,
        reference_build: a.reference_build.clone(),
        aligner: a.aligner.clone(),
        variant_caller: a.variant_caller.clone(),
    })
}

pub async fn list_for_run(pool: &SqlitePool, sequence_run_id: i64) -> Result<Vec<Alignment>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM alignment WHERE sequence_run_id = ? ORDER BY id"))
        .bind(sequence_run_id)
        .fetch_all(pool)
        .await?;
    Ok(rows.into_iter().map(Row::into_domain).collect())
}
