//! Chip-profile queries: a subject's genotyping-array QC summaries (one row per import).

use du_domain::ids::SampleGuid;
use navigator_domain::chipprofile::{ChipProfile, ChipSummary, NewChipProfile};
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    biosample_guid: String,
    provider: String,
    chip_version: Option<String>,
    total_markers_possible: i64,
    total_markers_called: i64,
    no_call_rate: f64,
    het_rate: Option<f64>,
    y_markers_called: i64,
    mt_markers_called: i64,
    autosomal_markers_called: i64,
    source_file_name: Option<String>,
    source_path: Option<String>,
}

impl Row {
    fn into_domain(self) -> Result<ChipProfile, StoreError> {
        let biosample_guid = parse_sample_guid(&self.biosample_guid, "chip_profile")?;
        Ok(ChipProfile {
            id: self.id,
            biosample_guid,
            provider: self.provider,
            chip_version: self.chip_version,
            summary: ChipSummary {
                total_markers_possible: self.total_markers_possible,
                total_markers_called: self.total_markers_called,
                no_call_rate: self.no_call_rate,
                het_rate: self.het_rate,
                y_markers_called: self.y_markers_called,
                mt_markers_called: self.mt_markers_called,
                autosomal_markers_called: self.autosomal_markers_called,
            },
            source_file_name: self.source_file_name,
            source_path: self.source_path,
        })
    }
}

const COLS: &str = "id, biosample_guid, provider, chip_version, total_markers_possible, \
    total_markers_called, no_call_rate, het_rate, y_markers_called, mt_markers_called, \
    autosomal_markers_called, source_file_name, source_path";

/// Insert a chip profile (the store assigns the id).
pub async fn create(pool: &SqlitePool, new: &NewChipProfile) -> Result<ChipProfile, StoreError> {
    let s = &new.summary;
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO chip_profile (biosample_guid, provider, chip_version, total_markers_possible, \
         total_markers_called, no_call_rate, het_rate, y_markers_called, mt_markers_called, \
         autosomal_markers_called, source_file_name, source_path) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING id",
    )
    .bind(new.biosample_guid.0.to_string())
    .bind(&new.provider)
    .bind(&new.chip_version)
    .bind(s.total_markers_possible)
    .bind(s.total_markers_called)
    .bind(s.no_call_rate)
    .bind(s.het_rate)
    .bind(s.y_markers_called)
    .bind(s.mt_markers_called)
    .bind(s.autosomal_markers_called)
    .bind(&new.source_file_name)
    .bind(&new.source_path)
    .fetch_one(pool)
    .await?;
    Ok(ChipProfile {
        id,
        biosample_guid: new.biosample_guid,
        provider: new.provider.clone(),
        chip_version: new.chip_version.clone(),
        summary: new.summary,
        source_file_name: new.source_file_name.clone(),
        source_path: new.source_path.clone(),
    })
}

/// Delete a chip profile (no child rows). Returns whether the row was removed.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM chip_profile WHERE id = ?")
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// A single chip profile by id.
pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<ChipProfile>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM chip_profile WHERE id = ?"))
        .bind(id)
        .fetch_optional(pool)
        .await?;
    row.map(Row::into_domain).transpose()
}

/// All chip profiles for a biosample.
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<ChipProfile>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM chip_profile WHERE biosample_guid = ? ORDER BY id"))
        .bind(guid.0.to_string())
        .fetch_all(pool)
        .await?;
    rows.into_iter().map(Row::into_domain).collect()
}
