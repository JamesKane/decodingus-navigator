//! Persistent PDS-publish outbox (sync durability, gap §5). A publish enqueues a fully-built
//! record; a background drain pushes it with exponential backoff. A transient/offline failure
//! reschedules the row (so it isn't lost); a non-transient failure marks it `FAILED`; a success
//! removes it (its outcome is logged in [`crate::sync_history`]).

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

/// Enqueue a publish. Re-publishing the same `(account_did, collection, entity_ref)` coalesces
/// onto the existing row: the newest payload wins and the row is reset to `PENDING` with a cleared
/// backoff/error, so a manual re-publish retries a previously-failed entry immediately.
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
pub async fn ready(pool: &SqlitePool, account_did: &str, now: &str, limit: i64) -> Result<Vec<OutboxEntry>, StoreError> {
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

/// Mark a row terminally `FAILED` (a non-transient error — e.g. validation/auth). It is not
/// auto-retried; a manual re-publish (which re-enqueues) resets it.
pub async fn mark_failed(pool: &SqlitePool, id: i64, attempt_count: i64, last_error: &str, now: &str) -> Result<(), StoreError> {
    sqlx::query("UPDATE sync_outbox SET status = 'FAILED', attempt_count = ?, last_error = ?, updated_at = ? WHERE id = ?")
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
    sqlx::query("DELETE FROM sync_outbox WHERE id = ?").bind(id).execute(pool).await?;
    Ok(())
}

/// Count of rows still awaiting a successful push for `account_did` (drives the UI indicator).
pub async fn pending_count(pool: &SqlitePool, account_did: &str) -> Result<i64, StoreError> {
    let (n,): (i64,) =
        sqlx::query_as("SELECT COUNT(*) FROM sync_outbox WHERE account_did = ? AND status = 'PENDING'")
            .bind(account_did)
            .fetch_one(pool)
            .await?;
    Ok(n)
}

/// All non-completed rows for `account_did` (PENDING + FAILED), newest first — for a sync detail view.
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
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T00:00:00Z").await.unwrap();
        enqueue(p, &entry("did:a", "alignment:2"), "2026-06-13T00:00:01Z").await.unwrap();
        // Re-publish alignment:1 with a new payload — coalesces onto the same row.
        let mut again = entry("did:a", "alignment:1");
        again.payload = r#"{"a":2}"#.into();
        enqueue(p, &again, "2026-06-13T00:00:02Z").await.unwrap();
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 2);

        let batch = ready(p, "did:a", "2026-06-13T01:00:00Z", 10).await.unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0].entity_ref, "alignment:1"); // oldest created_at first
        assert_eq!(batch[0].payload, r#"{"a":2}"#); // newest payload won
        // A different account's queue is isolated.
        assert_eq!(pending_count(p, "did:b").await.unwrap(), 0);
    }

    #[tokio::test]
    async fn reschedule_hides_until_due_then_complete_removes() {
        let s = crate::Store::open_in_memory().await.unwrap();
        let p = s.pool();
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T00:00:00Z").await.unwrap();
        let id = ready(p, "did:a", "2026-06-13T00:00:00Z", 10).await.unwrap()[0].id;

        reschedule(p, id, 1, "2026-06-13T02:00:00Z", "timeout", "2026-06-13T00:01:00Z").await.unwrap();
        // Not due yet → not returned, but still counts as pending.
        assert!(ready(p, "did:a", "2026-06-13T01:00:00Z", 10).await.unwrap().is_empty());
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 1);
        // Due → returned again with the bumped attempt count.
        let due = ready(p, "did:a", "2026-06-13T03:00:00Z", 10).await.unwrap();
        assert_eq!(due[0].attempt_count, 1);

        complete(p, id).await.unwrap();
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 0);
    }

    #[tokio::test]
    async fn mark_failed_drops_from_ready_but_re_enqueue_revives() {
        let s = crate::Store::open_in_memory().await.unwrap();
        let p = s.pool();
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T00:00:00Z").await.unwrap();
        let id = ready(p, "did:a", "2026-06-13T00:00:00Z", 10).await.unwrap()[0].id;
        mark_failed(p, id, 1, "invalid record", "2026-06-13T00:01:00Z").await.unwrap();
        assert!(ready(p, "did:a", "2026-06-13T02:00:00Z", 10).await.unwrap().is_empty());
        assert_eq!(pending_count(p, "did:a").await.unwrap(), 0); // FAILED not counted as pending
        assert_eq!(list(p, "did:a").await.unwrap().len(), 1); // but still visible

        // A manual re-publish resets it to PENDING.
        enqueue(p, &entry("did:a", "alignment:1"), "2026-06-13T03:00:00Z").await.unwrap();
        assert_eq!(ready(p, "did:a", "2026-06-13T03:00:00Z", 10).await.unwrap().len(), 1);
    }
}
