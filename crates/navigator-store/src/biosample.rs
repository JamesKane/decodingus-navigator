//! Biosample queries. The `SampleGuid` (UUID) is stored as its hyphenated TEXT form.

use du_domain::ids::SampleGuid;
use navigator_domain::workspace::Biosample;
use sqlx::SqlitePool;
use uuid::Uuid;

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
        let uuid = Uuid::parse_str(&self.guid)
            .map_err(|e| StoreError::Decode(format!("biosample guid {:?}: {e}", self.guid)))?;
        Ok(Biosample {
            guid: SampleGuid(uuid),
            sample_accession: self.sample_accession,
            donor_identifier: self.donor_identifier,
            description: self.description,
            center_name: self.center_name,
            sex: self.sex,
            project_id: self.project_id,
        })
    }
}

const COLS: &str =
    "guid, sample_accession, donor_identifier, description, center_name, sex, project_id";

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

pub async fn list_for_project(pool: &SqlitePool, project_id: i64) -> Result<Vec<Biosample>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM biosample WHERE project_id = ? ORDER BY guid"))
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
