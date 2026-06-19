//! Vendor-neutral Subject identifiers (FTDNA project-import design §4.2). `(source, external_id)` is
//! UNIQUE — one Subject per vendor id — and is the key the matching/dedup engine looks up.
//!
//! PII / never-federated: see the migration `0029_subject_identity` header.

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

/// Attach a vendor id to a Subject. Idempotent on `(source, external_id)`: a re-add for the **same**
/// biosample is a no-op; a conflicting `(source, external_id)` already bound to a **different**
/// biosample is left untouched (the matching engine resolves that — never silently re-point an id).
/// Returns the resulting row (existing or new).
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

/// Look up the Subject bound to a `(source, external_id)` — the exact-match step of dedup (§5.1).
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
        crate::biosample::create(
            pool,
            &Biosample {
                guid,
                sample_accession: None,
                donor_identifier: donor.into(),
                description: None,
                center_name: None,
                sex: None,
                project_id: None,
            },
        )
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
        // A conflicting (source, id) for a different subject does NOT steal the id.
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
    }
}
