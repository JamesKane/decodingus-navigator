//! `AsyncSync`: the completed sync engine (plan §6 and §7).
//!
//! It gives a PDS write the resilience that the old `AsyncSyncService` only sketched. There are
//! three parts. **Refresh-token rotation**, when the server refuses an access token. **A second
//! try with exponential backoff**, on a transient failure such as offline, a timeout, or a 5xx.
//! And an **offline indicator** that the UI can show. It returns a validation error (4xx) at once,
//! because a second try could never succeed.
//!
//! The conflict policy: a write goes through `createRecord` with a TID rkey that the server
//! generates. So two creates never collide, and each one gets its own key. An idempotent create,
//! update, or delete on an rkey that the caller chooses came later, once a record carried a stable
//! identity.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Duration;

use crate::error::SyncError;
use crate::oauth::refresh;
use crate::publish::{PdsClient, RecordRef, RemoteRecord};
use crate::tokens::{Session, TokenStore};

/// Retry/backoff schedule for transient failures.
#[derive(Debug, Clone, Copy)]
pub struct RetryPolicy {
    /// How many times to try again after the first try fails for a transient reason.
    pub max_retries: u32,
    /// Delay before the first retry; doubles each subsequent retry.
    pub base_delay: Duration,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        RetryPolicy {
            max_retries: 3,
            base_delay: Duration::from_millis(500),
        }
    }
}

impl RetryPolicy {
    /// Exponential backoff for the `attempt`-th retry (0-based): `base * 2^attempt`.
    pub fn backoff(&self, attempt: u32) -> Duration {
        self.base_delay * 2u32.saturating_pow(attempt)
    }
}

/// A resilient PDS writer for one authenticated account. It holds the live [`Session`] in memory,
/// refreshes it when it expires, and stores it again. It also holds an offline flag that it shares
/// with the app, so the indicator survives from one call to the next.
pub struct AsyncSync {
    http: reqwest::Client,
    tokens: TokenStore,
    /// The account DID. It is the keychain key that the rotated session goes under.
    did: String,
    session: Session,
    policy: RetryPolicy,
    online: Arc<AtomicBool>,
}

impl AsyncSync {
    /// Build an engine for `session`. It stores a rotated session under its DID, through `tokens`.
    /// It shares `online` with the app, for the offline indicator.
    pub fn new(
        http: reqwest::Client,
        tokens: TokenStore,
        session: Session,
        policy: RetryPolicy,
        online: Arc<AtomicBool>,
    ) -> Self {
        AsyncSync {
            http,
            tokens,
            did: session.did.clone(),
            session,
            policy,
            online,
        }
    }

    /// Whether the last write reached the server (true until a transport/5xx failure).
    pub fn is_online(&self) -> bool {
        self.online.load(Ordering::Relaxed)
    }

    /// The current session, which a refresh may have replaced. The caller reads it to see a
    /// rotation.
    pub fn session(&self) -> &Session {
        &self.session
    }

    /// Create a record. On a 401 it refreshes once, by itself. On a transient failure it tries
    /// again, with a backoff. On success it clears the offline flag.
    pub async fn push_create(&mut self, collection: &str, record: serde_json::Value) -> Result<RecordRef, SyncError> {
        self.push_create_inner(collection, record, None).await
    }

    /// Like [`push_create`](Self::push_create), but with an explicit record key. It is for an
    /// idempotent record of which there is only one. An example is the signing key of a device,
    /// keyed by its own `did:key`. A second registration then overwrites, and makes no duplicate.
    pub async fn push_create_rkey(
        &mut self,
        collection: &str,
        record: serde_json::Value,
        rkey: &str,
    ) -> Result<RecordRef, SyncError> {
        self.push_create_inner(collection, record, Some(rkey)).await
    }

    /// Upsert a record at a known `rkey`, with `putRecord`. This is the idempotent path for a
    /// second publish. It refreshes on a 401 and backs off on a transient failure, the same as
    /// [`push_create`](Self::push_create).
    pub async fn push_put(
        &mut self,
        collection: &str,
        rkey: &str,
        record: serde_json::Value,
    ) -> Result<RecordRef, SyncError> {
        let record = &record;
        self.with_resilience(|c| async move { c.put_record(collection, rkey, record.clone()).await })
            .await
    }

    /// Delete a record at `rkey`, with `deleteRecord`. This is the path that prunes an orphan. It
    /// keeps the same discipline as [`push_put`](Self::push_put): refresh on a 401, and back off on
    /// a transient failure.
    pub async fn push_delete(&mut self, collection: &str, rkey: &str) -> Result<(), SyncError> {
        self.with_resilience(|c| async move { c.delete_record(collection, rkey).await })
            .await
    }

    /// Fetch one page of the account's own records in `collection` (`listRecords`, for a PULL). Same
    /// refresh/backoff discipline. Returns the records + the next cursor.
    pub async fn pull_list(
        &mut self,
        collection: &str,
        cursor: Option<&str>,
    ) -> Result<(Vec<RemoteRecord>, Option<String>), SyncError> {
        self.with_resilience(|c| async move { c.list_records(collection, cursor).await })
            .await
    }

    async fn push_create_inner(
        &mut self,
        collection: &str,
        record: serde_json::Value,
        rkey: Option<&str>,
    ) -> Result<RecordRef, SyncError> {
        let record = &record;
        self.with_resilience(|c| async move { c.create_record(collection, record.clone(), rkey).await })
            .await
    }

    /// Run one PDS call under the engine's whole resilience discipline. It calls `op` again,
    /// against a client that it builds afresh, until the call succeeds or the engine stops.
    ///
    /// This is the *only* place that holds the policy. Every public method above is one line over
    /// it. So the refresh on a 401, the backoff, and the offline flag can never become different
    /// between the create, put, delete, and list paths.
    ///
    /// The engine calls `op` again for each try, which is why it is an `Fn`, and why each caller
    /// clones its record. A second try needs a client that the code builds from the session, and a
    /// rotation may have replaced that session.
    async fn with_resilience<T, F, Fut>(&mut self, op: F) -> Result<T, SyncError>
    where
        F: Fn(PdsClient) -> Fut,
        Fut: std::future::Future<Output = Result<T, SyncError>>,
    {
        let mut refreshed = false;
        let mut attempt = 0u32;
        loop {
            let client = PdsClient::from_session(self.http.clone(), &self.session)?;
            match op(client).await {
                Ok(v) => {
                    self.online.store(true, Ordering::Relaxed);
                    return Ok(v);
                }
                // Token expired/revoked: refresh once, persist, and retry immediately.
                Err(SyncError::Unauthorized) if !refreshed => {
                    self.session = refresh(&self.http, &self.session).await?;
                    self.tokens.save(&self.did, &self.session)?;
                    refreshed = true;
                }
                // Transient (offline/timeout/5xx): mark offline and back off, up to the cap.
                Err(e) if e.is_transient() && attempt < self.policy.max_retries => {
                    self.online.store(false, Ordering::Relaxed);
                    tokio::time::sleep(self.policy.backoff(attempt)).await;
                    attempt += 1;
                }
                Err(e) => {
                    if e.is_transient() {
                        self.online.store(false, Ordering::Relaxed);
                    }
                    return Err(e);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backoff_is_exponential() {
        let p = RetryPolicy {
            max_retries: 4,
            base_delay: Duration::from_millis(100),
        };
        assert_eq!(p.backoff(0), Duration::from_millis(100));
        assert_eq!(p.backoff(1), Duration::from_millis(200));
        assert_eq!(p.backoff(2), Duration::from_millis(400));
        assert_eq!(p.backoff(3), Duration::from_millis(800));
    }

    #[test]
    fn transient_classification() {
        // 5xx is transient; 401/oauth/validation are not.
        assert!(SyncError::Server(503, "down".into()).is_transient());
        assert!(SyncError::Server(500, "boom".into()).is_transient());
        assert!(!SyncError::Unauthorized.is_transient());
        assert!(!SyncError::Oauth("bad request".into()).is_transient());
    }

    /// A session that points at a dead endpoint. `push_create` tries the transient connect failure
    /// again, up to the cap, sets the offline flag, and then gives the error. That drives the path
    /// of the second try, the backoff, and the offline flag, with no network and no keychain.
    #[tokio::test]
    async fn push_create_retries_then_goes_offline() {
        use du_atproto::oauth::EcKey;

        let session = Session {
            did: "did:plc:test".into(),
            pds: "http://127.0.0.1:1".into(), // nothing listens here → connection refused
            access_token: "at".into(),
            refresh_token: "rt".into(),
            dpop_key_b64: EcKey::generate().to_base64(),
            scope: "atproto".into(),
            client_id: "http://localhost".into(),
        };
        let online = Arc::new(AtomicBool::new(true));
        let policy = RetryPolicy {
            max_retries: 2,
            base_delay: Duration::from_millis(1),
        };
        let mut engine = AsyncSync::new(
            reqwest::Client::new(),
            TokenStore::new("navigator-sync-test-offline"),
            session,
            policy,
            online.clone(),
        );

        let err = engine
            .push_create("com.decodingus.test", serde_json::json!({"x": 1}))
            .await;
        assert!(err.is_err(), "expected the dead endpoint to fail");
        assert!(err.unwrap_err().is_transient(), "connect failure should be transient");
        assert!(
            !engine.is_online(),
            "offline flag should be set after transient failures"
        );
        assert!(
            !online.load(Ordering::Relaxed),
            "shared flag should be visible to the app"
        );
    }
}
