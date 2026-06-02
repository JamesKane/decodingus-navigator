//! mtDNA-sequence queries: a subject's imported mitochondrial FASTA sequences.

use du_domain::ids::SampleGuid;
use navigator_domain::mtdna::{MtdnaSequence, NewMtdnaSequence};
use sqlx::SqlitePool;
use uuid::Uuid;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    biosample_guid: String,
    defline: Option<String>,
    sequence: String,
    n_count: i64,
    source_file_name: Option<String>,
}

impl Row {
    fn into_domain(self) -> Result<MtdnaSequence, StoreError> {
        let uuid = Uuid::parse_str(&self.biosample_guid)
            .map_err(|e| StoreError::Decode(format!("mtdna_sequence guid {:?}: {e}", self.biosample_guid)))?;
        Ok(MtdnaSequence {
            id: self.id,
            biosample_guid: SampleGuid(uuid),
            defline: self.defline,
            sequence: self.sequence,
            n_count: self.n_count,
            source_file_name: self.source_file_name,
        })
    }
}

const COLS: &str = "id, biosample_guid, defline, sequence, n_count, source_file_name";

/// Insert an mtDNA sequence (the store assigns the id).
pub async fn create(pool: &SqlitePool, new: &NewMtdnaSequence) -> Result<MtdnaSequence, StoreError> {
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO mtdna_sequence (biosample_guid, defline, sequence, n_count, source_file_name) \
         VALUES (?, ?, ?, ?, ?) RETURNING id",
    )
    .bind(new.biosample_guid.0.to_string())
    .bind(&new.defline)
    .bind(&new.sequence)
    .bind(new.n_count)
    .bind(&new.source_file_name)
    .fetch_one(pool)
    .await?;
    Ok(MtdnaSequence {
        id,
        biosample_guid: new.biosample_guid,
        defline: new.defline.clone(),
        sequence: new.sequence.clone(),
        n_count: new.n_count,
        source_file_name: new.source_file_name.clone(),
    })
}

/// Fetch one mtDNA sequence by id.
pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<MtdnaSequence>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM mtdna_sequence WHERE id = ?"))
        .bind(id)
        .fetch_optional(pool)
        .await?;
    row.map(Row::into_domain).transpose()
}

/// All mtDNA sequences for a biosample.
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<MtdnaSequence>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM mtdna_sequence WHERE biosample_guid = ? ORDER BY id"))
        .bind(guid.0.to_string())
        .fetch_all(pool)
        .await?;
    rows.into_iter().map(Row::into_domain).collect()
}
