//! The PDS-assigned identity of each published entity (gap §5-p2). On first publish the PDS assigns a
//! TID; we keep it here so the next publish UPDATES that record (putRecord at `rkey`) instead of
//! creating a duplicate. `payload_hash` is the published JSON's sha256 at push time — a PULL compares
//! it to the current local payload to tell whether local changed since the last push (conflict).

use sqlx::SqlitePool;

use crate::StoreError;

/// One published entity's remote identity + last-push fingerprint.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredSyncState {
    pub account_did: String,
    pub entity_ref: String,
    pub kind: String,
    pub collection: String,
    /// The TID the PDS assigned (the record key to putRecord at on update).
    pub rkey: String,
    pub at_uri: String,
    pub at_cid: String,
    pub payload_hash: String,
    pub pushed_at: String,
}

/// Insert or replace the sync-state for a published entity.
pub async fn upsert(pool: &SqlitePool, s: &StoredSyncState) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO sync_state (account_did, entity_ref, kind, collection, rkey, at_uri, at_cid, payload_hash, pushed_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(account_did, entity_ref) DO UPDATE SET \
         kind = excluded.kind, collection = excluded.collection, rkey = excluded.rkey, \
         at_uri = excluded.at_uri, at_cid = excluded.at_cid, payload_hash = excluded.payload_hash, \
         pushed_at = excluded.pushed_at",
    )
    .bind(&s.account_did)
    .bind(&s.entity_ref)
    .bind(&s.kind)
    .bind(&s.collection)
    .bind(&s.rkey)
    .bind(&s.at_uri)
    .bind(&s.at_cid)
    .bind(&s.payload_hash)
    .bind(&s.pushed_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// The sync-state for one entity, if it has been published.
pub async fn get(pool: &SqlitePool, account_did: &str, entity_ref: &str) -> Result<Option<StoredSyncState>, StoreError> {
    let row = sqlx::query_as("SELECT * FROM sync_state WHERE account_did = ? AND entity_ref = ?")
        .bind(account_did)
        .bind(entity_ref)
        .fetch_optional(pool)
        .await?;
    Ok(row)
}

/// All published entities for an account in one collection (for a PULL reconcile pass).
pub async fn list_for_collection(pool: &SqlitePool, account_did: &str, collection: &str) -> Result<Vec<StoredSyncState>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM sync_state WHERE account_did = ? AND collection = ?")
        .bind(account_did)
        .bind(collection)
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

/// All published entities for an account.
pub async fn list_for_account(pool: &SqlitePool, account_did: &str) -> Result<Vec<StoredSyncState>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM sync_state WHERE account_did = ?")
        .bind(account_did)
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(entity: &str, cid: &str) -> StoredSyncState {
        StoredSyncState {
            account_did: "did:plc:abc".into(),
            entity_ref: entity.into(),
            kind: "coverage".into(),
            collection: "com.decodingus.alignment".into(),
            rkey: "3kab2c".into(),
            at_uri: "at://did:plc:abc/com.decodingus.alignment/3kab2c".into(),
            at_cid: cid.into(),
            payload_hash: "hash1".into(),
            pushed_at: "2026-06-17T00:00:00Z".into(),
        }
    }

    #[tokio::test]
    async fn upsert_get_list_round_trip() {
        let store = crate::Store::open_in_memory().await.unwrap();
        assert!(get(store.pool(), "did:plc:abc", "alignment:1").await.unwrap().is_none());
        upsert(store.pool(), &row("alignment:1", "cidA")).await.unwrap();
        let got = get(store.pool(), "did:plc:abc", "alignment:1").await.unwrap().unwrap();
        assert_eq!(got.rkey, "3kab2c");
        assert_eq!(got.at_cid, "cidA");
        // Upsert replaces the cid (a re-push) without duplicating.
        let mut updated = row("alignment:1", "cidB");
        updated.payload_hash = "hash2".into();
        upsert(store.pool(), &updated).await.unwrap();
        let got = get(store.pool(), "did:plc:abc", "alignment:1").await.unwrap().unwrap();
        assert_eq!(got.at_cid, "cidB");
        assert_eq!(got.payload_hash, "hash2");
        assert_eq!(list_for_account(store.pool(), "did:plc:abc").await.unwrap().len(), 1);
        assert_eq!(list_for_collection(store.pool(), "did:plc:abc", "com.decodingus.alignment").await.unwrap().len(), 1);
    }
}
