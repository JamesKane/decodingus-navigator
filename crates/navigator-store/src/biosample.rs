//! Biosample queries. The `SampleGuid` is a UUID, and it goes in as its hyphenated TEXT form.

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

/// Find a subject that exists, by its donor identifier. A person is one subject across every
/// project, so an importer reuses that subject and makes no duplicate. The result is deterministic,
/// in `guid` order, when more than one subject shares an identifier. The match uses case, which is
/// how the store holds an identifier.
pub async fn find_by_donor(pool: &SqlitePool, donor_identifier: &str) -> Result<Option<Biosample>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM biosample WHERE donor_identifier = ? ORDER BY guid LIMIT 1"
    ))
    .bind(donor_identifier)
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

/// Update the biosample fields that a user can edit. `donor_identifier` must have a value. Each
/// other field can be null, and an empty value clears it. It returns whether it changed a row.
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

/// Set the biosample's `sample_accession`, and nothing else. One example: change a friendly-name
/// placeholder to the authoritative catalog accession that the app fetched from the AppView. It
/// leaves every other field as it is, and returns whether it changed a row.
pub async fn set_accession(pool: &SqlitePool, guid: SampleGuid, accession: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE biosample SET sample_accession = ? WHERE guid = ?")
        .bind(accession)
        .bind(guid.0.to_string())
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Assign the biosample's project, or clear it with `None`. It returns whether it changed a row.
pub async fn set_project(pool: &SqlitePool, guid: SampleGuid, project_id: Option<i64>) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE biosample SET project_id = ? WHERE guid = ?")
        .bind(project_id)
        .bind(guid.0.to_string())
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Clear the legacy home-project column of every subject whose home is `project_id`. A delete of a
/// project uses it, and the subjects survive with no home. It returns the count that it cleared.
pub async fn clear_home_project(pool: &SqlitePool, project_id: i64) -> Result<u64, StoreError> {
    let affected = sqlx::query("UPDATE biosample SET project_id = NULL WHERE project_id = ?")
        .bind(project_id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected)
}

/// Clear **every piece of sequencing data, and every derived or imported analysis**, for a
/// subject, in one transaction.
///
/// It keeps the biosample row itself, and its identity. That identity is the name, the sex, the
/// center, the vendor external IDs, the project memberships, and the MDKA genealogy.
///
/// This is the "reset this subject" maintenance operation. The explicit *Clear data* action uses
/// it, and so does [`delete`], as its first step, so a delete can never leave an orphan row.
///
/// It removes the sequencing runs, then their alignments, then the cached analysis artifacts, and
/// it unlinks the source files.
///
/// It also removes the Y and mt haplogroup calls, and the genome consensus. Then every cache that a
/// signature keys, from [`crate::sig_cache::ALL`], which is the painting, the ROH, and both archaic
/// tiers. Then the reconciliation overrides and the audit log, the ancestry results, the IBD
/// exchange results, and the mtDNA sequences. Last, the chip, STR, and variant profiles, with their
/// child rows. It is idempotent.
pub async fn clear_data(pool: &SqlitePool, guid: SampleGuid) -> Result<(), StoreError> {
    let g = guid.0.to_string();
    let mut tx = pool.begin().await?;
    // The alignments that belong to this subject, through its runs. They drive the deletes that
    // key on an alignment.
    const ALN: &str =
        "SELECT id FROM alignment WHERE sequence_run_id IN (SELECT id FROM sequence_run WHERE biosample_guid = ?)";
    // The children that key on an alignment go first: the artifacts, and then the source-file
    // rows. The code DELETES a source file, and does not only unlink it. A clear removes the
    // subject's sequencing data, so the file identity must go with it.
    //
    // A `source_file` that stays behind keeps a `UNIQUE content_sha256`. A second import of the
    // same file would then dedup to a dead row. That orphan blocks the import, and says nothing,
    // after a delete. A delete needs the clear first.
    sqlx::query(&format!("DELETE FROM analysis_artifact WHERE alignment_id IN ({ALN})"))
        .bind(&g)
        .execute(&mut *tx)
        .await?;
    sqlx::query(&format!("DELETE FROM source_file WHERE alignment_id IN ({ALN})"))
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
    // The derived and imported tables that key on a biosample. The biosample row itself stays.
    //
    // The caches that a signature keys come from `sig_cache::ALL`, and this code does not list
    // them. When a hand-written list stood here, the loop named `consensus_painting` alone. The ROH
    // cache, and both archaic caches, survived a "clear this subject's data". Each one still keyed
    // to a consensus signature that no longer existed.
    for table in [
        "haplogroup_call",
        "consensus_profile",
        "reconciliation_override",
        "reconciliation_audit",
        "ancestry_result",
        "ibd_exchange_result",
        "mtdna_sequence",
        "str_profile",
        "variant_set",
        "chip_profile",
    ]
    .into_iter()
    .chain(crate::sig_cache::ALL.iter().map(|c| c.table()))
    {
        sqlx::query(&format!("DELETE FROM {table} WHERE biosample_guid = ?"))
            .bind(&g)
            .execute(&mut *tx)
            .await?;
    }
    tx.commit().await?;
    Ok(())
}

/// Clear a subject's haplogroup placement state, and *nothing else*, for Y and for mtDNA. That
/// state is the call of each source, the pooled consensus snapshot, and the reconciliation override
/// with its audit.
///
/// It leaves the coverage, the ancestry, the imported profiles, the runs, and the alignments as
/// they are.
///
/// The app uses it to drop a stale legacy lineage, such as an old trail with ISOGG names. The next
/// analysis then places the subject again, consistently.
pub async fn clear_haplogroup_data(pool: &SqlitePool, guid: SampleGuid) -> Result<(), StoreError> {
    let g = guid.0.to_string();
    let mut tx = pool.begin().await?;
    for table in [
        "haplogroup_call",
        "consensus_profile",
        "reconciliation_override",
        "reconciliation_audit",
    ] {
        sqlx::query(&format!("DELETE FROM {table} WHERE biosample_guid = ?"))
            .bind(&g)
            .execute(&mut *tx)
            .await?;
    }
    tx.commit().await?;
    Ok(())
}

/// Delete a biosample row. It returns whether it removed a row. A caller must make sure that no
/// dependent row names it, such as a sequence run or a profile. The app layer guards that.
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

/// Every member of a project. It is the union of two things. The first is the M:N membership
/// table, which is the source of truth. The second is the legacy `biosample.project_id` home
/// column, because an older import wrote no membership row.
///
/// An FTDNA import that merges a subject into a project gives it a membership row, and the subject
/// keeps its original home column. So the report must read both. The guid dedupes the result.
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

/// Every biosample, whatever its project. A biosample is an entity in its own right, and the link
/// to a project is optional. The order is by donor identifier, which keeps the subjects list
/// stable.
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

/// The member count of **every** project at once, as `(project_id, count)`.
///
/// The membership rule is the same as in [`count_members_for_project`]. It takes the M:N table
/// together with the legacy home column, dedupes by guid, and counts only a guid that still exists
/// in `biosample`.
///
/// It is one `GROUP BY`, and not a query for each project. A project with no members is absent from
/// the result, and the query does not report it as zero.
pub async fn member_counts(pool: &SqlitePool) -> Result<Vec<(i64, i64)>, StoreError> {
    // The inner UNION, which is not a UNION ALL, dedupes the `(project, guid)` pairs. So a subject
    // that is both an M:N member and has the project as its legacy home counts once. That matches
    // the `IN (...)` behaviour of the single-project query.
    //
    // The join to `biosample` drops a membership row whose subject is gone. The
    // `SELECT COUNT(*) FROM biosample WHERE guid IN (...)` form also left those out.
    let rows: Vec<(i64, i64)> = sqlx::query_as(
        "SELECT m.project_id, COUNT(*) FROM ( \
             SELECT project_id, biosample_guid AS guid FROM biosample_project \
             UNION \
             SELECT project_id, guid FROM biosample WHERE project_id IS NOT NULL \
           ) m \
           JOIN biosample b ON b.guid = m.guid \
         GROUP BY m.project_id",
    )
    .fetch_all(pool)
    .await?;
    Ok(rows)
}

/// The count of every project member, from the M:N membership table together with the legacy home
/// column, deduped by guid. It matches [`list_members_for_project`]. The sample badge of the
/// projects list uses it.
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
