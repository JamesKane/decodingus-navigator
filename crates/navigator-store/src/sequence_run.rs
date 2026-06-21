//! Sequence-run queries. Read metrics are flat columns (not a JSON blob).

use du_domain::ids::SampleGuid;
use navigator_domain::workspace::{NewSequenceRun, SequenceRun};
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
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
    sequencing_facility: Option<String>,
    instrument_id: Option<String>,
    sample_name: Option<String>,
    library_id: Option<String>,
    platform_unit: Option<String>,
    flowcell_id: Option<String>,
}

impl Row {
    fn into_domain(self) -> Result<SequenceRun, StoreError> {
        let biosample_guid = parse_sample_guid(&self.biosample_guid, "sequence_run")?;
        Ok(SequenceRun {
            id: self.id,
            biosample_guid,
            platform_name: self.platform_name,
            instrument_model: self.instrument_model,
            test_type: self.test_type,
            library_layout: self.library_layout,
            total_reads: self.total_reads,
            pf_reads_aligned: self.pf_reads_aligned,
            mean_read_length: self.mean_read_length,
            mean_insert_size: self.mean_insert_size,
            sequencing_facility: self.sequencing_facility,
            instrument_id: self.instrument_id,
            sample_name: self.sample_name,
            library_id: self.library_id,
            platform_unit: self.platform_unit,
            flowcell_id: self.flowcell_id,
        })
    }
}

const COLS: &str = "id, biosample_guid, platform_name, instrument_model, test_type, \
    library_layout, total_reads, pf_reads_aligned, mean_read_length, mean_insert_size, \
    sequencing_facility, instrument_id, sample_name, library_id, platform_unit, flowcell_id";

/// Fetch one sequence run by id.
pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<SequenceRun>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM sequence_run WHERE id = ?"))
        .bind(id)
        .fetch_optional(pool)
        .await?;
    row.map(Row::into_domain).transpose()
}

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
        // The lab/instrument identity block is filled in post-create by `set_library_stats`.
        sequencing_facility: None,
        instrument_id: None,
        sample_name: None,
        library_id: None,
        platform_unit: None,
        flowcell_id: None,
    })
}

/// Persist the lab/instrument identity block inferred from the alignment at import (read-name
/// scan + `@RG` tags). Does not touch `sequencing_facility` (set separately via [`update`], or by
/// a later instrument→lab resolution). Returns whether a row was affected.
#[allow(clippy::too_many_arguments)]
pub async fn set_library_stats(
    pool: &SqlitePool,
    id: i64,
    instrument_id: Option<&str>,
    sample_name: Option<&str>,
    library_id: Option<&str>,
    platform_unit: Option<&str>,
    flowcell_id: Option<&str>,
) -> Result<bool, StoreError> {
    let affected = sqlx::query(
        "UPDATE sequence_run SET instrument_id = ?, sample_name = ?, library_id = ?, \
         platform_unit = ?, flowcell_id = ? WHERE id = ?",
    )
    .bind(instrument_id)
    .bind(sample_name)
    .bind(library_id)
    .bind(platform_unit)
    .bind(flowcell_id)
    .bind(id)
    .execute(pool)
    .await?
    .rows_affected();
    Ok(affected > 0)
}

/// Set the library-level read stats (`total_reads`, `mean_read_length`, `mean_insert_size`,
/// `library_layout`) — populated after a read-metrics / unified-walker pass (or backfilled from a
/// cached artifact). These describe the run's library; per-alignment counts (e.g. reads aligned)
/// live on the alignment. A `None` `library_layout` leaves the existing value (set at import from
/// the BAM flags). Leaves the descriptive + lab columns untouched. Returns whether a row was
/// affected.
pub async fn set_read_stats(
    pool: &SqlitePool,
    id: i64,
    total_reads: Option<i64>,
    mean_read_length: Option<f64>,
    mean_insert_size: Option<f64>,
    library_layout: Option<&str>,
) -> Result<bool, StoreError> {
    let affected = sqlx::query(
        "UPDATE sequence_run SET total_reads = ?, mean_read_length = ?, mean_insert_size = ?, \
         library_layout = COALESCE(?, library_layout) WHERE id = ?",
    )
    .bind(total_reads)
    .bind(mean_read_length)
    .bind(mean_insert_size)
    .bind(library_layout)
    .bind(id)
    .execute(pool)
    .await?
    .rows_affected();
    Ok(affected > 0)
}

/// Set only the sequencing facility (the lab) — used by the AppView instrument→lab resolution,
/// which leaves the analysis-derived columns untouched. Returns whether a row was affected.
pub async fn set_facility(pool: &SqlitePool, id: i64, facility: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE sequence_run SET sequencing_facility = ? WHERE id = ?")
        .bind(facility)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Set a run's test-type code (e.g. normalizing a generic `TARGETED_Y` to `BIG_Y_700` once the
/// run's vendor is known to be FTDNA, which only sells Big Y). Returns whether a row was affected.
pub async fn set_test_type(pool: &SqlitePool, id: i64, test_type: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE sequence_run SET test_type = ? WHERE id = ?")
        .bind(test_type)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Update a run's descriptive fields. The analysis-derived read-metric columns (total_reads,
/// pf_reads_aligned, mean_read_length, mean_insert_size) are left untouched. Returns whether a
/// row was affected.
pub async fn update(
    pool: &SqlitePool,
    id: i64,
    platform_name: &str,
    instrument_model: Option<&str>,
    test_type: &str,
    library_layout: Option<&str>,
    sequencing_facility: Option<&str>,
) -> Result<bool, StoreError> {
    let affected = sqlx::query(
        "UPDATE sequence_run SET platform_name = ?, instrument_model = ?, test_type = ?, \
         library_layout = ?, sequencing_facility = ? WHERE id = ?",
    )
    .bind(platform_name)
    .bind(instrument_model)
    .bind(test_type)
    .bind(library_layout)
    .bind(sequencing_facility)
    .bind(id)
    .execute(pool)
    .await?
    .rows_affected();
    Ok(affected > 0)
}

/// Delete a sequence run and everything beneath it (its alignments and their cached analysis
/// artifacts), children-first since FKs are enforced. Returns whether the run row was removed.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let mut tx = pool.begin().await?;
    sqlx::query(
        "DELETE FROM analysis_artifact WHERE alignment_id IN \
         (SELECT id FROM alignment WHERE sequence_run_id = ?)",
    )
    .bind(id)
    .execute(&mut *tx)
    .await?;
    // Unlink content-hash file records pointing at this run's alignments (keep the file identity).
    sqlx::query(
        "UPDATE source_file SET alignment_id = NULL WHERE alignment_id IN \
         (SELECT id FROM alignment WHERE sequence_run_id = ?)",
    )
    .bind(id)
    .execute(&mut *tx)
    .await?;
    sqlx::query("DELETE FROM alignment WHERE sequence_run_id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?;
    let affected = sqlx::query("DELETE FROM sequence_run WHERE id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?
        .rows_affected();
    tx.commit().await?;
    Ok(affected > 0)
}

pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<SequenceRun>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM sequence_run WHERE biosample_guid = ? ORDER BY id"
    ))
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}
