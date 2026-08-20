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
    total_bases: Option<i64>,
    read_type: Option<String>,
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
            total_bases: self.total_bases,
            read_type: self.read_type,
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
    total_bases, read_type, \
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
        // `set_library_stats` fills in `read_type` after the create, during the import scan.
        // `set_read_stats` fills in `total_bases`, during the read-metrics pass. Both of them
        // support the standardized test label.
        total_bases: None,
        read_type: None,
        // `set_library_stats` fills in the lab and instrument identity block after the create.
        sequencing_facility: None,
        instrument_id: None,
        sample_name: None,
        library_id: None,
        platform_unit: None,
        flowcell_id: None,
    })
}

/// Store the lab and instrument identity block that the import inferred from the alignment, with a
/// scan of the read names and the `@RG` tags. It does not touch `sequencing_facility`, which
/// [`update`] sets, or which a later instrument→lab resolve sets. It returns whether it changed a
/// row.
#[allow(clippy::too_many_arguments)]
pub async fn set_library_stats(
    pool: &SqlitePool,
    id: i64,
    instrument_id: Option<&str>,
    sample_name: Option<&str>,
    library_id: Option<&str>,
    platform_unit: Option<&str>,
    flowcell_id: Option<&str>,
    read_type: Option<&str>,
) -> Result<bool, StoreError> {
    let affected = sqlx::query(
        "UPDATE sequence_run SET instrument_id = ?, sample_name = ?, library_id = ?, \
         platform_unit = ?, flowcell_id = ?, read_type = COALESCE(?, read_type) WHERE id = ?",
    )
    .bind(instrument_id)
    .bind(sample_name)
    .bind(library_id)
    .bind(platform_unit)
    .bind(flowcell_id)
    .bind(read_type)
    .bind(id)
    .execute(pool)
    .await?
    .rows_affected();
    Ok(affected > 0)
}

/// Set the read stats of the library: `total_reads`, `mean_read_length`, `mean_insert_size`,
/// `library_layout`, and `total_bases`. A read-metrics or unified-walker pass fills them, or a
/// backfill takes them from a cached artifact.
///
/// A `total_bases` of `None` keeps the value that is there. A `library_layout` of `None` does the
/// same, and the import sets that field from the BAM flags.
///
/// These describe the run's library. A count that belongs to one alignment, such as the reads
/// aligned, lives on the alignment. This leaves the descriptive columns and the lab columns as they
/// are. It returns whether it changed a row.
pub async fn set_read_stats(
    pool: &SqlitePool,
    id: i64,
    total_reads: Option<i64>,
    mean_read_length: Option<f64>,
    mean_insert_size: Option<f64>,
    library_layout: Option<&str>,
    total_bases: Option<i64>,
) -> Result<bool, StoreError> {
    let affected = sqlx::query(
        "UPDATE sequence_run SET total_reads = ?, mean_read_length = ?, mean_insert_size = ?, \
         library_layout = COALESCE(?, library_layout), total_bases = COALESCE(?, total_bases) \
         WHERE id = ?",
    )
    .bind(total_reads)
    .bind(mean_read_length)
    .bind(mean_insert_size)
    .bind(library_layout)
    .bind(total_bases)
    .bind(id)
    .execute(pool)
    .await?
    .rows_affected();
    Ok(affected > 0)
}

/// Set the sequencing facility, which is the lab, and nothing else. The AppView instrument→lab
/// resolve uses it, and leaves the columns that come from the analysis as they are. It returns
/// whether it changed a row.
pub async fn set_facility(pool: &SqlitePool, id: i64, facility: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE sequence_run SET sequencing_facility = ? WHERE id = ?")
        .bind(facility)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Set a run's read chemistry or mode: `SHORT`, `HIFI`, `CLR`, or an `ONT_*` value. That is the
/// long-read arm of the standardized test label. The backfill uses it to fill the field on a run
/// that came in before the field existed. It leaves the rest of the library-stats block as it is,
/// and returns whether it changed a row.
pub async fn set_read_type(pool: &SqlitePool, id: i64, read_type: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE sequence_run SET read_type = ? WHERE id = ?")
        .bind(read_type)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Set a run's test-type code. One example: change a generic `TARGETED_Y` to `BIG_Y_700`, once the
/// app knows that the run's vendor is FTDNA, which sells Big Y alone. It returns whether it changed
/// a row.
pub async fn set_test_type(pool: &SqlitePool, id: i64, test_type: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE sequence_run SET test_type = ? WHERE id = ?")
        .bind(test_type)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Update a run's descriptive fields. It leaves the read-metric columns that come from the
/// analysis as they are: total_reads, pf_reads_aligned, mean_read_length, and mean_insert_size. It
/// returns whether it changed a row.
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

/// Delete a sequence run and everything below it, which is its alignments and their cached
/// analysis artifacts. The children go first, because the database enforces the FKs. It returns
/// whether it removed the run row.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let mut tx = pool.begin().await?;
    sqlx::query(
        "DELETE FROM analysis_artifact WHERE alignment_id IN \
         (SELECT id FROM alignment WHERE sequence_run_id = ?)",
    )
    .bind(id)
    .execute(&mut *tx)
    .await?;
    // Unlink the content-hash file records that name this run's alignments. Keep the file
    // identity itself.
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
