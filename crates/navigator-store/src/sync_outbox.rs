//! The durable PDS-publish outbox (sync durability, gap §5). A publish puts a complete record in
//! the queue, and a background drain pushes it, with an exponential backoff.
//!
//! A transient or offline failure gives the row a new time, so nothing is lost. A failure that is
//! not transient marks the row `FAILED`. A success removes the row, and
//! [`crate::sync_history`] logs the outcome.

use sqlx::SqlitePool;

use crate::StoreError;

/// A queued publish. `rkey` is `Some` only for idempotent singletons (server-assigned TID
/// otherwise). `next_retry_at` is `None` when the row is ready to send now.
#[derive(Debug, Clone, PartialEq, Eq, sqlx::FromRow)]
pub struct OutboxEntry {
    pub id: i64,
    pub account_did: String,
    pub kind: String,
    pub entity_ref: String,
    pub collection: String,
    pub rkey: Option<String>,
    pub payload: String,
    pub status: String,
    pub attempt_count: i64,
    pub next_retry_at: Option<String>,
    pub last_error: Option<String>,
    pub created_at: String,
    pub updated_at: String,
}

/// The fields needed to enqueue (or coalesce onto an existing) outbox row.
#[derive(Debug, Clone)]
pub struct NewOutboxEntry {
    pub account_did: String,
    pub kind: String,
    pub entity_ref: String,
    pub collection: String,
    pub rkey: Option<String>,
    pub payload: String,
}

/// Put a publish in the queue. A second publish of the same `(account_did, collection,
/// entity_ref)` joins the row that is already there. The newest payload wins, and the row goes back
/// to `PENDING`, with its backoff and its error cleared. So a manual re-publish tries a failed
/// entry again, immediately.
pub async fn enqueue(pool: &SqlitePool, e: &NewOutboxEntry, now: &str) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO sync_outbox \
         (account_did, kind, entity_ref, collection, rkey, payload, status, attempt_count, \
          next_retry_at, last_error, created_at, updated_at) \
         VALUES (?, ?, ?, ?, ?, ?, 'PENDING', 0, NULL, NULL, ?, ?) \
         ON CONFLICT(account_did, collection, entity_ref) DO UPDATE SET \
         kind = excluded.kind, rkey = excluded.rkey, payload = excluded.payload, \
         status = 'PENDING', attempt_count = 0, next_retry_at = NULL, last_error = NULL, \
         updated_at = excluded.updated_at",
    )
    .bind(&e.account_did)
    .bind(&e.kind)
    .bind(&e.entity_ref)
    .bind(&e.collection)
    .bind(&e.rkey)
    .bind(&e.payload)
    .bind(now)
    .bind(now)
    .execute(pool)
    .await?;
    Ok(())
}

/// The next batch of ready rows for `account_did`: `PENDING` and due (no backoff, or `next_retry_at`
/// at/before `now`), oldest first.
pub async fn ready(
    pool: &SqlitePool,
    account_did: &str,
    now: &str,
    limit: i64,
) -> Result<Vec<OutboxEntry>, StoreError> {
    let rows = sqlx::query_as::<_, OutboxEntry>(
        "SELECT * FROM sync_outbox \
         WHERE account_did = ? AND status = 'PENDING' \
         AND (next_retry_at IS NULL OR next_retry_at <= ?) \
         ORDER BY created_at ASC, id ASC LIMIT ?",
    )
    .bind(account_did)
    .bind(now)
    .bind(limit)
    .fetch_all(pool)
    .await?;
    Ok(rows)
}

/// Reschedule a row after a transient failure: bump `attempt_count`, set `next_retry_at`, keep it
/// `PENDING`.
pub async fn reschedule(
    pool: &SqlitePool,
    id: i64,
    attempt_count: i64,
    next_retry_at: &str,
    last_error: &str,
    now: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "UPDATE sync_outbox SET attempt_count = ?, next_retry_at = ?, last_error = ?, \
         status = 'PENDING', updated_at = ? WHERE id = ?",
    )
    .bind(attempt_count)
    .bind(next_retry_at)
    .bind(last_error)
    .bind(now)
    .bind(id)
    .execute(pool)
    .await?;
    Ok(())
}

/// Mark a row `FAILED`, which is final. The cause is an error that is not transient, such as a
/// validation error or an auth error. The drain never tries it again. A manual re-publish puts it
/// in the queue afresh, and resets it.
pub async fn mark_failed(
    pool: &SqlitePool,
    id: i64,
    attempt_count: i64,
    last_error: &str,
    now: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "UPDATE sync_outbox SET status = 'FAILED', attempt_count = ?, last_error = ?, updated_at = ? WHERE id = ?",
    )
    .bind(attempt_count)
    .bind(last_error)
    .bind(now)
    .bind(id)
    .execute(pool)
    .await?;
    Ok(())
}

/// Remove a row after a successful push (the outcome lives in `sync_history`).
pub async fn complete(pool: &SqlitePool, id: i64) -> Result<(), StoreError> {
    sqlx::query("DELETE FROM sync_outbox WHERE id = ?")
        .bind(id)
        .execute(pool)
        .await?;
    Ok(())
}

/// The count of rows that still wait for a successful push, for `account_did`. The UI indicator
/// reads it.
pub async fn pending_count(pool: &SqlitePool, account_did: &str) -> Result<i64, StoreError> {
    let (n,): (i64,) = sqlx::query_as("SELECT COUNT(*) FROM sync_outbox WHERE account_did = ? AND status = 'PENDING'")
        .bind(account_did)
        .fetch_one(pool)
        .await?;
    Ok(n)
}

/// Every row for `account_did` that is not complete, which is PENDING and FAILED, newest first.
/// It feeds a sync detail view.
pub async fn list(pool: &SqlitePool, account_did: &str) -> Result<Vec<OutboxEntry>, StoreError> {
    let rows = sqlx::query_as::<_, OutboxEntry>(
        "SELECT * FROM sync_outbox WHERE account_did = ? ORDER BY updated_at DESC, id DESC",
    )
    .bind(account_did)
    .fetch_all(pool)
    .await?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(did: &str, reff: &str) -> NewOutboxEntry {
        NewOutboxEntry {
            account_did: did.into(),
            kind: "coverage".into(),
            entity_ref: reff.into(),
            collection: "com.decodingus.alignment".into(),
            rkey: None,
            payload: r#"{"a":1}"#.into(),
        }
    }

    #[tokio::test]
    async fn enqueue_coalesces_and_drains_in_order() {
        let s = crate::Store::open_in_memory().await.unwrap();
        let p = s.pool();
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T00:00:00Z")
            .await
            .unwrap();
        enqueue(p, &entry("did:a", "alignment:2"), "2026-06-13T00:00:01Z")
            .await
            .unwrap();
        // A second publish of alignment:1, with a new payload, joins the same row.
        let mut again = entry("did:a", "alignment:1");
        again.payload = r#"{"a":2}"#.into();
        enqueue(p, &again, "2026-06-13T00:00:02Z").await.unwrap();
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 2);

        let batch = ready(p, "did:a", "2026-06-13T01:00:00Z", 10).await.unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0].entity_ref, "alignment:1"); // oldest created_at first
        assert_eq!(batch[0].payload, r#"{"a":2}"#); // newest payload won
                                                    // Another account has its own queue.
        assert_eq!(pending_count(p, "did:b").await.unwrap(), 0);
    }

    #[tokio::test]
    async fn reschedule_hides_until_due_then_complete_removes() {
        let s = crate::Store::open_in_memory().await.unwrap();
        let p = s.pool();
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T00:00:00Z")
            .await
            .unwrap();
        let id = ready(p, "did:a", "2026-06-13T00:00:00Z", 10).await.unwrap()[0].id;

        reschedule(p, id, 1, "2026-06-13T02:00:00Z", "timeout", "2026-06-13T00:01:00Z")
            .await
            .unwrap();
        // Not due yet → not returned, but still counts as pending.
        assert!(ready(p, "did:a", "2026-06-13T01:00:00Z", 10).await.unwrap().is_empty());
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 1);
        // The row is due, so it comes back with a higher try count.
        let due = ready(p, "did:a", "2026-06-13T03:00:00Z", 10).await.unwrap();
        assert_eq!(due[0].attempt_count, 1);

        complete(p, id).await.unwrap();
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 0);
    }

    #[tokio::test]
    async fn mark_failed_drops_from_ready_but_re_enqueue_revives() {
        let s = crate::Store::open_in_memory().await.unwrap();
        let p = s.pool();
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T00:00:00Z")
            .await
            .unwrap();
        let id = ready(p, "did:a", "2026-06-13T00:00:00Z", 10).await.unwrap()[0].id;
        mark_failed(p, id, 1, "invalid record", "2026-06-13T00:01:00Z")
            .await
            .unwrap();
        assert!(ready(p, "did:a", "2026-06-13T02:00:00Z", 10).await.unwrap().is_empty());
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 0); // FAILED not counted as pending
        assert_eq!(list(p, "did:a").await.unwrap().len(), 1); // but still visible

        // A manual re-publish resets it to PENDING.
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T03:00:00Z")
            .await
            .unwrap();
        assert_eq!(ready(p, "did:a", "2026-06-13T03:00:00Z", 10).await.unwrap().len(), 1);
    }
}
