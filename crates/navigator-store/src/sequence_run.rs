//! Sequence-run queries. Read metrics are flat columns (not a JSON blob).

use du_domain::ids::SampleGuid;
use navigator_domain::workspace::{NewSequenceRun, SequenceRun};
use sqlx::SqlitePool;
use uuid::Uuid;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    biosample_guid: String,
    platform_name: String,
    instrument_model: Option<String>,
    test_type: String,
    library_layout: Option<String>,
    total_reads: Option<i64>,
    pf_reads_aligned: Option<i64>,
    mean_read_length: Option<f64>,
    mean_insert_size: Option<f64>,
}

impl Row {
    fn into_domain(self) -> Result<SequenceRun, StoreError> {
        let uuid = Uuid::parse_str(&self.biosample_guid)
            .map_err(|e| StoreError::Decode(format!("sequence_run biosample_guid: {e}")))?;
        Ok(SequenceRun {
            id: self.id,
            biosample_guid: SampleGuid(uuid),
            platform_name: self.platform_name,
            instrument_model: self.instrument_model,
            test_type: self.test_type,
            library_layout: self.library_layout,
            total_reads: self.total_reads,
            pf_reads_aligned: self.pf_reads_aligned,
            mean_read_length: self.mean_read_length,
            mean_insert_size: self.mean_insert_size,
        })
    }
}

const COLS: &str = "id, biosample_guid, platform_name, instrument_model, test_type, \
    library_layout, total_reads, pf_reads_aligned, mean_read_length, mean_insert_size";

pub async fn create(pool: &SqlitePool, r: &NewSequenceRun) -> Result<SequenceRun, StoreError> {
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO sequence_run (biosample_guid, platform_name, instrument_model, test_type, \
         library_layout, total_reads, pf_reads_aligned, mean_read_length, mean_insert_size) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING id",
    )
    .bind(r.biosample_guid.0.to_string())
    .bind(&r.platform_name)
    .bind(&r.instrument_model)
    .bind(&r.test_type)
    .bind(&r.library_layout)
    .bind(r.total_reads)
    .bind(r.pf_reads_aligned)
    .bind(r.mean_read_length)
    .bind(r.mean_insert_size)
    .fetch_one(pool)
    .await?;
    Ok(SequenceRun {
        id,
        biosample_guid: r.biosample_guid,
        platform_name: r.platform_name.clone(),
        instrument_model: r.instrument_model.clone(),
        test_type: r.test_type.clone(),
        library_layout: r.library_layout.clone(),
        total_reads: r.total_reads,
        pf_reads_aligned: r.pf_reads_aligned,
        mean_read_length: r.mean_read_length,
        mean_insert_size: r.mean_insert_size,
    })
}

pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<SequenceRun>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM sequence_run WHERE biosample_guid = ? ORDER BY id"))
        .bind(guid.0.to_string())
        .fetch_all(pool)
        .await?;
    rows.into_iter().map(Row::into_domain).collect()
}
