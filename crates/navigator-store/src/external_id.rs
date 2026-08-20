//! Vendor-neutral Subject identifiers (FTDNA project-import design §4.2). `(source, external_id)`
//! is UNIQUE, so one vendor id belongs to one Subject. The match and dedup engine looks that pair
//! up.
//!
//! This is PII, and it never goes to the federation. See the header of the migration
//! `0029_subject_identity`.

use du_domain::ids::SampleGuid;
use navigator_domain::identity::ExternalId;
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    biosample_guid: String,
    source: String,
    external_id: String,
}

impl Row {
    fn into_domain(self) -> Result<ExternalId, StoreError> {
        Ok(ExternalId {
            id: self.id,
            biosample_guid: parse_sample_guid(&self.biosample_guid, "external_id")?,
            source: self.source,
            external_id: self.external_id,
        })
    }
}

const COLS: &str = "id, biosample_guid, source, external_id";

/// Attach a vendor id to a Subject. It is idempotent on `(source, external_id)`.
///
/// A second add for the **same** biosample does nothing. A `(source, external_id)` that already
/// belongs to a **different** biosample stays as it is. The match engine resolves such a case, and
/// this code must never move an id to another subject with no word.
///
/// It returns the row, whether that row was there before or is new.
pub async fn add(
    pool: &SqlitePool,
    guid: SampleGuid,
    source: &str,
    external_id: &str,
) -> Result<ExternalId, StoreError> {
    sqlx::query(
        "INSERT INTO external_id (biosample_guid, source, external_id) VALUES (?, ?, ?) \
         ON CONFLICT(source, external_id) DO NOTHING",
    )
    .bind(guid.0.to_string())
    .bind(source)
    .bind(external_id)
    .execute(pool)
    .await?;
    // Read back the canonical row (may belong to a pre-existing owner on conflict).
    find(pool, source, external_id)
        .await?
        .ok_or_else(|| StoreError::NotFound(format!("external_id {source}:{external_id}")))
}

/// Detach a vendor id by its row id. It returns `true` when it removed a row. The subject editor
/// uses it to drop a kit association. The row is PII that stays in the workspace, so a hard delete
/// is correct.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let res = sqlx::query("DELETE FROM external_id WHERE id = ?")
        .bind(id)
        .execute(pool)
        .await?;
    Ok(res.rows_affected() > 0)
}

/// Fetch one vendor-id row by its id. The caller uses it to find the Subject that owns the row,
/// before a delete. It can then publish the biosample anchor again, with the new set of ids.
pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<ExternalId>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM external_id WHERE id = ?"))
        .bind(id)
        .fetch_optional(pool)
        .await?;
    row.map(Row::into_domain).transpose()
}

/// Look up the Subject that a `(source, external_id)` pair belongs to. This is the exact-match
/// step of the dedup (§5.1).
pub async fn find(pool: &SqlitePool, source: &str, external_id: &str) -> Result<Option<ExternalId>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM external_id WHERE source = ? AND external_id = ?"
    ))
    .bind(source)
    .bind(external_id)
    .fetch_optional(pool)
    .await?;
    row.map(Row::into_domain).transpose()
}

/// All vendor ids for a Subject.
pub async fn list_for(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<ExternalId>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM external_id WHERE biosample_guid = ? ORDER BY source, external_id"
    ))
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_domain::identity::IdSource;
    use navigator_domain::workspace::Biosample;

    async fn seed(pool: &SqlitePool, donor: &str) -> SampleGuid {
        let guid = SampleGuid(uuid::Uuid::new_v4());
        crate::biosample::create(pool, &Biosample::new(guid, donor))
            .await
            .unwrap();
        guid
    }

    #[tokio::test]
    async fn unique_vendor_id_and_lookup() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        let a = seed(pool, "GFX").await;
        let b = seed(pool, "other").await;

        let row = add(pool, a, IdSource::FTDNA, "B5163").await.unwrap();
        assert_eq!(row.biosample_guid, a);
        // Re-add for the same subject is a no-op (same row).
        let again = add(pool, a, IdSource::FTDNA, "B5163").await.unwrap();
        assert_eq!(again.id, row.id);
        // A (source, id) that another subject already holds does NOT move to this subject.
        let conflict = add(pool, b, IdSource::FTDNA, "B5163").await.unwrap();
        assert_eq!(
            conflict.biosample_guid, a,
            "exact-match must resolve to the original owner"
        );

        // Exact-match lookup (the dedup anchor).
        assert_eq!(
            find(pool, IdSource::FTDNA, "B5163")
                .await
                .unwrap()
                .unwrap()
                .biosample_guid,
            a
        );
        assert!(find(pool, IdSource::FTDNA, "NOPE").await.unwrap().is_none());
        // Same id under a different source is distinct.
        add(pool, b, IdSource::YSEQ, "B5163").await.unwrap();
        assert_eq!(
            find(pool, IdSource::YSEQ, "B5163")
                .await
                .unwrap()
                .unwrap()
                .biosample_guid,
            b
        );
        assert_eq!(list_for(pool, a).await.unwrap().len(), 1);

        // A delete by row id detaches the id, and frees the (source, id) pair for another
        // subject.
        let row = find(pool, IdSource::FTDNA, "B5163").await.unwrap().unwrap();
        assert!(delete(pool, row.id).await.unwrap(), "row removed");
        assert!(!delete(pool, row.id).await.unwrap(), "second delete is a no-op");
        assert!(find(pool, IdSource::FTDNA, "B5163").await.unwrap().is_none());
        let rebound = add(pool, b, IdSource::FTDNA, "B5163").await.unwrap();
        assert_eq!(rebound.biosample_guid, b, "freed id rebinds to the new owner");
    }
}
