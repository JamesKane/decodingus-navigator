//! Append-only audit trail of completed PDS-publish attempts (sync durability, gap §5). One row is
//! written per terminal outcome — a successful push (with the resulting `at://` URI + CID) or a
//! non-transient failure. Transient retries don't write history (they stay in [`crate::sync_outbox`]).

use sqlx::SqlitePool;

use crate::StoreError;

/// A recorded push outcome.
#[derive(Debug, Clone, PartialEq, Eq, sqlx::FromRow)]
pub struct HistoryEntry {
    pub id: i64,
    pub account_did: String,
    pub kind: String,
    pub entity_ref: String,
    pub collection: String,
    pub direction: String,
    pub status: String,
    pub at_uri: Option<String>,
    pub at_cid: Option<String>,
    pub attempt_count: i64,
    pub error: Option<String>,
    pub created_at: String,
}

/// The fields to log for a push outcome.
#[derive(Debug, Clone)]
pub struct NewHistoryEntry {
    pub account_did: String,
    pub kind: String,
    pub entity_ref: String,
    pub collection: String,
    pub status: String,
    pub at_uri: Option<String>,
    pub at_cid: Option<String>,
    pub attempt_count: i64,
    pub error: Option<String>,
}

/// Append an outcome row (direction is always `PUSH` for now).
pub async fn record(pool: &SqlitePool, e: &NewHistoryEntry, now: &str) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO sync_history \
         (account_did, kind, entity_ref, collection, direction, status, at_uri, at_cid, \
          attempt_count, error, created_at) \
         VALUES (?, ?, ?, ?, 'PUSH', ?, ?, ?, ?, ?, ?)",
    )
    .bind(&e.account_did)
    .bind(&e.kind)
    .bind(&e.entity_ref)
    .bind(&e.collection)
    .bind(&e.status)
    .bind(&e.at_uri)
    .bind(&e.at_cid)
    .bind(e.attempt_count)
    .bind(&e.error)
    .bind(now)
    .execute(pool)
    .await?;
    Ok(())
}

/// The most recent `limit` outcomes for `account_did`, newest first.
pub async fn recent(pool: &SqlitePool, account_did: &str, limit: i64) -> Result<Vec<HistoryEntry>, StoreError> {
    let rows = sqlx::query_as::<_, HistoryEntry>(
        "SELECT * FROM sync_history WHERE account_did = ? ORDER BY created_at DESC, id DESC LIMIT ?",
    )
    .bind(account_did)
    .bind(limit)
    .fetch_all(pool)
    .await?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn record_and_read_back_recent() {
        let s = crate::Store::open_in_memory().await.unwrap();
        let p = s.pool();
        let ok = NewHistoryEntry {
            account_did: "did:a".into(),
            kind: "coverage".into(),
            entity_ref: "alignment:1".into(),
            collection: "com.decodingus.alignment".into(),
            status: "SUCCESS".into(),
            at_uri: Some("at://did:a/c/1".into()),
            at_cid: Some("bafy...".into()),
            attempt_count: 0,
            error: None,
        };
        record(p, &ok, "2026-06-13T00:00:00Z").await.unwrap();
        let mut fail = ok.clone();
        fail.status = "FAILED".into();
        fail.at_uri = None;
        fail.at_cid = None;
        fail.error = Some("invalid record".into());
        record(p, &fail, "2026-06-13T00:01:00Z").await.unwrap();

        let rows = recent(p, "did:a", 10).await.unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].status, "FAILED"); // newest first
        assert_eq!(rows[1].at_uri.as_deref(), Some("at://did:a/c/1"));
        assert!(recent(p, "did:b", 10).await.unwrap().is_empty());
    }
}
