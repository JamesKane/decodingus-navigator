//! Biosample queries. The `SampleGuid` (UUID) is stored as its hyphenated TEXT form.

use du_domain::ids::SampleGuid;
use navigator_domain::workspace::Biosample;
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    guid: String,
    sample_accession: Option<String>,
    donor_identifier: String,
    description: Option<String>,
    center_name: Option<String>,
    sex: Option<String>,
    project_id: Option<i64>,
}

impl Row {
    fn into_domain(self) -> Result<Biosample, StoreError> {
        let guid = parse_sample_guid(&self.guid, "biosample")?;
        Ok(Biosample {
            guid,
            sample_accession: self.sample_accession,
            donor_identifier: self.donor_identifier,
            description: self.description,
            center_name: self.center_name,
            sex: self.sex,
            project_id: self.project_id,
        })
    }
}

const COLS: &str = "guid, sample_accession, donor_identifier, description, center_name, sex, project_id";

/// Insert a biosample (the caller assigns the `SampleGuid`).
pub async fn create(pool: &SqlitePool, b: &Biosample) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO biosample (guid, sample_accession, donor_identifier, description, center_name, sex, project_id) \
         VALUES (?, ?, ?, ?, ?, ?, ?)",
    )
    .bind(b.guid.0.to_string())
    .bind(&b.sample_accession)
    .bind(&b.donor_identifier)
    .bind(&b.description)
    .bind(&b.center_name)
    .bind(&b.sex)
    .bind(b.project_id)
    .execute(pool)
    .await?;
    Ok(())
}

pub async fn get(pool: &SqlitePool, guid: SampleGuid) -> Result<Option<Biosample>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM biosample WHERE guid = ?"))
        .bind(guid.0.to_string())
        .fetch_optional(pool)
        .await?;
    row.map(Row::into_domain).transpose()
}

/// Set the biosample's recorded sex (e.g. write back an inferred sex when the user left it blank).
pub async fn set_sex(pool: &SqlitePool, guid: SampleGuid, sex: &str) -> Result<(), StoreError> {
    sqlx::query("UPDATE biosample SET sex = ? WHERE guid = ?")
        .bind(sex)
        .bind(guid.0.to_string())
        .execute(pool)
        .await?;
    Ok(())
}

/// Update the user-editable biosample fields. `donor_identifier` is required; the rest are
/// nullable (an empty value clears them). Returns whether a row was affected.
pub async fn update(
    pool: &SqlitePool,
    guid: SampleGuid,
    donor_identifier: &str,
    sample_accession: Option<&str>,
    description: Option<&str>,
    center_name: Option<&str>,
    sex: Option<&str>,
) -> Result<bool, StoreError> {
    let affected = sqlx::query(
        "UPDATE biosample SET donor_identifier = ?, sample_accession = ?, description = ?, \
         center_name = ?, sex = ? WHERE guid = ?",
    )
    .bind(donor_identifier)
    .bind(sample_accession)
    .bind(description)
    .bind(center_name)
    .bind(sex)
    .bind(guid.0.to_string())
    .execute(pool)
    .await?
    .rows_affected();
    Ok(affected > 0)
}

/// Assign (or clear, with `None`) the biosample's project. Returns whether a row was affected.
pub async fn set_project(pool: &SqlitePool, guid: SampleGuid, project_id: Option<i64>) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE biosample SET project_id = ? WHERE guid = ?")
        .bind(project_id)
        .bind(guid.0.to_string())
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Clear **all sequencing + derived/imported analysis data** for a subject in one transaction,
/// leaving the biosample row itself (and its identity: name/sex/center, vendor external IDs,
/// project memberships, and MDKA genealogy) intact. The "reset this subject" maintenance op —
/// used both by the explicit *Clear data* action and as the pre-step of [`delete`] (so a delete
/// can never orphan rows). Removes: sequencing runs → alignments → cached analysis artifacts
/// (and unlinks source files); Y/mt haplogroup calls + genome consensus + chromosome painting;
/// reconciliation overrides + audit log; ancestry results; IBD exchange results; mtDNA sequences;
/// and chip / STR / variant profiles (with their child rows). Idempotent.
pub async fn clear_data(pool: &SqlitePool, guid: SampleGuid) -> Result<(), StoreError> {
    let g = guid.0.to_string();
    let mut tx = pool.begin().await?;
    // The alignments that belong to this subject (via its runs) — drives the alignment-keyed deletes.
    const ALN: &str =
        "SELECT id FROM alignment WHERE sequence_run_id IN (SELECT id FROM sequence_run WHERE biosample_guid = ?)";
    // Alignment-keyed children first (artifacts, then unlink source files so the file identity survives).
    sqlx::query(&format!("DELETE FROM analysis_artifact WHERE alignment_id IN ({ALN})"))
        .bind(&g)
        .execute(&mut *tx)
        .await?;
    sqlx::query(&format!(
        "UPDATE source_file SET alignment_id = NULL WHERE alignment_id IN ({ALN})"
    ))
    .bind(&g)
    .execute(&mut *tx)
    .await?;
    sqlx::query(
        "DELETE FROM alignment WHERE sequence_run_id IN (SELECT id FROM sequence_run WHERE biosample_guid = ?)",
    )
    .bind(&g)
    .execute(&mut *tx)
    .await?;
    sqlx::query("DELETE FROM sequence_run WHERE biosample_guid = ?")
        .bind(&g)
        .execute(&mut *tx)
        .await?;
    // Profile children (markers/calls) before their parents.
    sqlx::query("DELETE FROM str_marker WHERE str_profile_id IN (SELECT id FROM str_profile WHERE biosample_guid = ?)")
        .bind(&g)
        .execute(&mut *tx)
        .await?;
    sqlx::query(
        "DELETE FROM variant_call WHERE variant_set_id IN (SELECT id FROM variant_set WHERE biosample_guid = ?)",
    )
    .bind(&g)
    .execute(&mut *tx)
    .await?;
    // Biosample-keyed derived + imported tables (the biosample row itself is kept).
    for table in [
        "haplogroup_call",
        "consensus_profile",
        "consensus_painting",
        "reconciliation_override",
        "reconciliation_audit",
        "ancestry_result",
        "ibd_exchange_result",
        "mtdna_sequence",
        "str_profile",
        "variant_set",
        "chip_profile",
    ] {
        sqlx::query(&format!("DELETE FROM {table} WHERE biosample_guid = ?"))
            .bind(&g)
            .execute(&mut *tx)
            .await?;
    }
    tx.commit().await?;
    Ok(())
}

/// Delete a biosample row. Returns whether a row was removed. Callers must ensure no
/// dependent rows reference it (sequence runs, profiles, etc.) — the app layer guards this.
pub async fn delete(pool: &SqlitePool, guid: SampleGuid) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM biosample WHERE guid = ?")
        .bind(guid.0.to_string())
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

pub async fn list_for_project(pool: &SqlitePool, project_id: i64) -> Result<Vec<Biosample>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM biosample WHERE project_id = ? ORDER BY guid"
    ))
    .bind(project_id)
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

/// All members of a project — the union of the M:N membership table (the source of truth) and the
/// legacy `biosample.project_id` home column (older imports never wrote a membership row). A subject
/// merged into a project by FTDNA import gets a membership row but keeps its original home column, so
/// the report must read both. Deduped by guid.
pub async fn list_members_for_project(pool: &SqlitePool, project_id: i64) -> Result<Vec<Biosample>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM biosample WHERE guid IN ( \
           SELECT biosample_guid FROM biosample_project WHERE project_id = ? \
           UNION \
           SELECT guid FROM biosample WHERE project_id = ? \
         ) ORDER BY donor_identifier, guid"
    ))
    .bind(project_id)
    .bind(project_id)
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

/// Every biosample, regardless of project (biosamples are first-class — the project link
/// is optional). Ordered by donor identifier for a stable subjects list.
pub async fn list_all(pool: &SqlitePool) -> Result<Vec<Biosample>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM biosample ORDER BY donor_identifier, guid"))
        .fetch_all(pool)
        .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

pub async fn count_for_project(pool: &SqlitePool, project_id: i64) -> Result<i64, StoreError> {
    let n: i64 = sqlx::query_scalar("SELECT COUNT(*) FROM biosample WHERE project_id = ?")
        .bind(project_id)
        .fetch_one(pool)
        .await?;
    Ok(n)
}

/// Count of all project members (M:N membership ∪ legacy home column), deduped by guid — matches
/// [`list_members_for_project`]. Used for the projects-list sample badge.
pub async fn count_members_for_project(pool: &SqlitePool, project_id: i64) -> Result<i64, StoreError> {
    let n: i64 = sqlx::query_scalar(
        "SELECT COUNT(*) FROM biosample WHERE guid IN ( \
           SELECT biosample_guid FROM biosample_project WHERE project_id = ? \
           UNION \
           SELECT guid FROM biosample WHERE project_id = ? \
         )",
    )
    .bind(project_id)
    .bind(project_id)
    .fetch_one(pool)
    .await?;
    Ok(n)
}
