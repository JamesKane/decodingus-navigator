//! `AsyncSync` — the completed sync engine (plan §6/§7). Wraps PDS writes with the
//! resilience the old `AsyncSyncService` only stubbed: **refresh-token rotation** on a
//! rejected access token, **retry with exponential backoff** on transient failures
//! (offline, timeout, 5xx), and an **offline indicator** the UI can surface. Validation
//! errors (4xx) are returned immediately — retrying them would never succeed.
//!
//! Conflict policy: writes go through `createRecord` with a server-generated TID rkey, so
//! two creates never collide (each gets its own key). Idempotent create/update/delete on a
//! caller-chosen rkey is a later addition once records carry stable identities.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Duration;

use crate::error::SyncError;
use crate::oauth::refresh;
use crate::publish::{PdsClient, RecordRef};
use crate::tokens::{Session, TokenStore};

/// Retry/backoff schedule for transient failures.
#[derive(Debug, Clone, Copy)]
pub struct RetryPolicy {
    /// How many times to retry after the first attempt fails transiently.
    pub max_retries: u32,
    /// Delay before the first retry; doubles each subsequent retry.
    pub base_delay: Duration,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        RetryPolicy { max_retries: 3, base_delay: Duration::from_millis(500) }
    }
}

impl RetryPolicy {
    /// Exponential backoff for the `attempt`-th retry (0-based): `base * 2^attempt`.
    pub fn backoff(&self, attempt: u32) -> Duration {
        self.base_delay * 2u32.saturating_pow(attempt)
    }
}

/// A resilient PDS writer for one authenticated account. Holds the live [`Session`] in
/// memory (refreshing and re-persisting it on expiry) and an offline flag shared with the
/// app so the indicator survives across calls.
pub struct AsyncSync {
    http: reqwest::Client,
    tokens: TokenStore,
    /// Account DID — the keychain key the rotated session is saved under.
    did: String,
    session: Session,
    policy: RetryPolicy,
    online: Arc<AtomicBool>,
}

impl AsyncSync {
    /// Build an engine for `session`, persisting any rotated session under its DID via
    /// `tokens`. `online` is shared with the app for the offline indicator.
    pub fn new(
        http: reqwest::Client,
        tokens: TokenStore,
        session: Session,
        policy: RetryPolicy,
        online: Arc<AtomicBool>,
    ) -> Self {
        AsyncSync { http, tokens, did: session.did.clone(), session, policy, online }
    }

    /// Whether the last write reached the server (true until a transport/5xx failure).
    pub fn is_online(&self) -> bool {
        self.online.load(Ordering::Relaxed)
    }

    /// The current (possibly refreshed) session — so the caller can observe rotation.
    pub fn session(&self) -> &Session {
        &self.session
    }

    /// Create a record, transparently refreshing on 401 (once) and retrying transient
    /// failures with backoff. On success the offline flag is cleared.
    pub async fn push_create(
        &mut self,
        collection: &str,
        record: serde_json::Value,
    ) -> Result<RecordRef, SyncError> {
        self.push_create_inner(collection, record, None).await
    }

    /// Like [`push_create`](Self::push_create) but with an explicit record key — for
    /// idempotent singleton-style records (e.g. the per-device signing key, keyed by its
    /// own `did:key` so re-registration overwrites rather than duplicates).
    pub async fn push_create_rkey(
        &mut self,
        collection: &str,
        record: serde_json::Value,
        rkey: &str,
    ) -> Result<RecordRef, SyncError> {
        self.push_create_inner(collection, record, Some(rkey)).await
    }

    async fn push_create_inner(
        &mut self,
        collection: &str,
        record: serde_json::Value,
        rkey: Option<&str>,
    ) -> Result<RecordRef, SyncError> {
        let mut refreshed = false;
        let mut attempt = 0u32;
        loop {
            let client = PdsClient::from_session(self.http.clone(), &self.session)?;
            match client.create_record(collection, record.clone(), rkey).await {
                Ok(r) => {
                    self.online.store(true, Ordering::Relaxed);
                    return Ok(r);
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
        let p = RetryPolicy { max_retries: 4, base_delay: Duration::from_millis(100) };
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

    /// A session pointed at a dead endpoint: push_create retries the transient connect
    /// failure up to the cap, flips the offline flag, and finally surfaces the error —
    /// exercising the retry/backoff/offline path without a network or keychain.
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
        let policy = RetryPolicy { max_retries: 2, base_delay: Duration::from_millis(1) };
        let mut engine = AsyncSync::new(
            reqwest::Client::new(),
            TokenStore::new("navigator-sync-test-offline"),
            session,
            policy,
            online.clone(),
        );

        let err = engine.push_create("com.decodingus.test", serde_json::json!({"x": 1})).await;
        assert!(err.is_err(), "expected the dead endpoint to fail");
        assert!(err.unwrap_err().is_transient(), "connect failure should be transient");
        assert!(!engine.is_online(), "offline flag should be set after transient failures");
        assert!(!online.load(Ordering::Relaxed), "shared flag should be visible to the app");
    }
}
