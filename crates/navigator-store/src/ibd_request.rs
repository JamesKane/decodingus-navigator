//! The federated-IBD **request ledger** — durable state for a matching conversation from the
//! introduction request through consent to a completed exchange. Its completed counterpart is
//! [`crate::ibd_exchange`], which stores the *result*; this table is what makes the in-flight
//! middle of that story survive a restart.
//!
//! Keyed by the broker's `request_uri`. Rows are plain data: the lifecycle `status` and
//! `direction` are stored as TEXT and given meaning by `navigator-app` (the same convention as
//! [`crate::ibd_exchange::StoredIbdExchange::relationship`]).

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// One matching conversation. See `migrations/0041_ibd_request.up.sql` for the field semantics —
/// in particular that `my_sample_ref` / `partner_sample_ref` are **AppView** sample handles while
/// `biosample_guid` is the local subject, and that `consent_given` records only our own decision.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredIbdRequest {
    pub request_uri: String,
    pub direction: String,
    pub purpose: String,
    pub status: String,
    pub partner_did: Option<String>,
    pub session_id: Option<String>,
    pub biosample_guid: Option<String>,
    pub my_sample_ref: Option<String>,
    pub partner_sample_ref: Option<String>,
    pub consent_given: Option<bool>,
    pub consent_at: Option<String>,
    pub attested_at: Option<String>,
    pub last_error: Option<String>,
    pub created_at: String,
    pub updated_at: String,
}

/// Insert or replace a request. `created_at` is preserved from the existing row (an update never
/// rewrites when the conversation began).
pub async fn upsert(pool: &SqlitePool, r: &StoredIbdRequest) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO ibd_request (request_uri, direction, purpose, status, partner_did, session_id, \
         biosample_guid, my_sample_ref, partner_sample_ref, consent_given, consent_at, attested_at, \
         last_error, created_at, updated_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(request_uri) DO UPDATE SET \
         direction = excluded.direction, purpose = excluded.purpose, status = excluded.status, \
         partner_did = excluded.partner_did, session_id = excluded.session_id, \
         biosample_guid = excluded.biosample_guid, my_sample_ref = excluded.my_sample_ref, \
         partner_sample_ref = excluded.partner_sample_ref, consent_given = excluded.consent_given, \
         consent_at = excluded.consent_at, attested_at = excluded.attested_at, \
         last_error = excluded.last_error, updated_at = excluded.updated_at",
    )
    .bind(&r.request_uri)
    .bind(&r.direction)
    .bind(&r.purpose)
    .bind(&r.status)
    .bind(&r.partner_did)
    .bind(&r.session_id)
    .bind(&r.biosample_guid)
    .bind(&r.my_sample_ref)
    .bind(&r.partner_sample_ref)
    .bind(r.consent_given)
    .bind(&r.consent_at)
    .bind(&r.attested_at)
    .bind(&r.last_error)
    .bind(&r.created_at)
    .bind(&r.updated_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// Insert only if the request is unknown — the reconciler's primitive for adopting a request the
/// broker reports. An existing row keeps every local field (notably our consent decision and the
/// subject we chose), so re-polling can never walk them back.
pub async fn insert_if_absent(pool: &SqlitePool, r: &StoredIbdRequest) -> Result<bool, StoreError> {
    let res = sqlx::query(
        "INSERT INTO ibd_request (request_uri, direction, purpose, status, partner_did, session_id, \
         biosample_guid, my_sample_ref, partner_sample_ref, consent_given, consent_at, attested_at, \
         last_error, created_at, updated_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(request_uri) DO NOTHING",
    )
    .bind(&r.request_uri)
    .bind(&r.direction)
    .bind(&r.purpose)
    .bind(&r.status)
    .bind(&r.partner_did)
    .bind(&r.session_id)
    .bind(&r.biosample_guid)
    .bind(&r.my_sample_ref)
    .bind(&r.partner_sample_ref)
    .bind(r.consent_given)
    .bind(&r.consent_at)
    .bind(&r.attested_at)
    .bind(&r.last_error)
    .bind(&r.created_at)
    .bind(&r.updated_at)
    .execute(pool)
    .await?;
    Ok(res.rows_affected() > 0)
}

/// One request by its broker URI.
pub async fn get(pool: &SqlitePool, request_uri: &str) -> Result<Option<StoredIbdRequest>, StoreError> {
    let row = sqlx::query_as("SELECT * FROM ibd_request WHERE request_uri = ?")
        .bind(request_uri)
        .fetch_optional(pool)
        .await?;
    Ok(row)
}

/// All requests, newest first.
pub async fn list(pool: &SqlitePool) -> Result<Vec<StoredIbdRequest>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM ibd_request ORDER BY created_at DESC")
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

/// Requests bound to one local subject, newest first.
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<StoredIbdRequest>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM ibd_request WHERE biosample_guid = ? ORDER BY created_at DESC")
        .bind(guid.0.to_string())
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

/// Forget a request (the user dismissed it locally). The broker keeps its own record.
pub async fn delete(pool: &SqlitePool, request_uri: &str) -> Result<(), StoreError> {
    sqlx::query("DELETE FROM ibd_request WHERE request_uri = ?")
        .bind(request_uri)
        .execute(pool)
        .await?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    fn row(uri: &str, status: &str) -> StoredIbdRequest {
        StoredIbdRequest {
            request_uri: uri.into(),
            direction: "OUTBOUND".into(),
            purpose: "IBD_AUTOSOMAL".into(),
            status: status.into(),
            partner_did: None,
            session_id: None,
            biosample_guid: None,
            my_sample_ref: Some("sample-mine".into()),
            partner_sample_ref: Some("sample-theirs".into()),
            consent_given: None,
            consent_at: None,
            attested_at: None,
            last_error: None,
            created_at: "2026-08-01T00:00:00Z".into(),
            updated_at: "2026-08-01T00:00:00Z".into(),
        }
    }

    async fn store_with_subject() -> (crate::Store, SampleGuid) {
        let store = crate::Store::open_in_memory().await.unwrap();
        let g = SampleGuid(Uuid::new_v4());
        let bio = navigator_domain::workspace::Biosample {
            guid: g,
            sample_accession: None,
            donor_identifier: "S1".into(),
            description: None,
            center_name: None,
            sex: None,
            project_id: None,
        };
        crate::biosample::create(store.pool(), &bio).await.unwrap();
        (store, g)
    }

    #[tokio::test]
    async fn upsert_advances_status_but_keeps_created_at() {
        let (store, g) = store_with_subject().await;
        let mut r = row("urn:ibd:abc", "REQUESTED");
        r.biosample_guid = Some(g.0.to_string());
        upsert(store.pool(), &r).await.unwrap();

        r.status = "READY".into();
        r.partner_did = Some("did:key:zB".into());
        r.session_id = Some("sess-1".into());
        r.updated_at = "2026-08-02T00:00:00Z".into();
        upsert(store.pool(), &r).await.unwrap();

        let got = get(store.pool(), "urn:ibd:abc").await.unwrap().unwrap();
        assert_eq!(got.status, "READY");
        assert_eq!(got.partner_did.as_deref(), Some("did:key:zB"));
        assert_eq!(got.created_at, "2026-08-01T00:00:00Z", "created_at is never rewritten");
        assert_eq!(got.updated_at, "2026-08-02T00:00:00Z");
        assert_eq!(list_for_biosample(store.pool(), g).await.unwrap().len(), 1);
    }

    /// Re-polling the broker must not walk back a decision we already made locally.
    #[tokio::test]
    async fn insert_if_absent_preserves_local_consent() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let mut mine = row("urn:ibd:xyz", "DECLINED");
        mine.direction = "INBOUND".into();
        mine.consent_given = Some(false);
        mine.consent_at = Some("2026-08-01T12:00:00Z".into());
        assert!(insert_if_absent(store.pool(), &mine).await.unwrap());

        // The broker still lists it as awaiting consent; adopting it again is a no-op.
        let fresh = row("urn:ibd:xyz", "AWAITING_CONSENT");
        assert!(!insert_if_absent(store.pool(), &fresh).await.unwrap());

        let got = get(store.pool(), "urn:ibd:xyz").await.unwrap().unwrap();
        assert_eq!(got.status, "DECLINED");
        assert_eq!(got.consent_given, Some(false));
        assert_eq!(got.direction, "INBOUND");
    }

    #[tokio::test]
    async fn list_and_delete() {
        let store = crate::Store::open_in_memory().await.unwrap();
        upsert(store.pool(), &row("urn:ibd:a", "REQUESTED")).await.unwrap();
        upsert(store.pool(), &row("urn:ibd:b", "EXCHANGED")).await.unwrap();
        assert_eq!(list(store.pool()).await.unwrap().len(), 2);
        delete(store.pool(), "urn:ibd:a").await.unwrap();
        let rest = list(store.pool()).await.unwrap();
        assert_eq!(rest.len(), 1);
        assert_eq!(rest[0].request_uri, "urn:ibd:b");
        assert!(get(store.pool(), "urn:ibd:a").await.unwrap().is_none());
    }
}
