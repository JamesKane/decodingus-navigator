//! `impl App` methods extracted from `lib.rs` (the `auth` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- authentication ----------------------------------------------------

    /// Run the public-client OAuth login for `handle` (handle or DID): browser authorize →
    /// loopback callback → token exchange. On success the DPoP-bound session is persisted
    /// to the OS keychain and becomes the active account. Returns the authenticated DID.
    pub async fn login(&self, handle: &str) -> Result<String, AppError> {
        let session = login_default(&self.auth.http, &self.auth.config, handle).await?;
        let did = session.did.clone();
        self.auth.tokens.save(&did, &session)?;
        self.auth.tokens.set_active(&did)?;
        *self.auth.active.lock().unwrap() = Some(did.clone());
        Ok(did)
    }

    /// The signed-in account's DID, or `None`.
    pub fn current_account(&self) -> Option<String> {
        self.auth.active.lock().unwrap().clone()
    }

    /// The signed-in account's DID, or [`AppError::NotAuthenticated`] — the cheap auth guard publish
    /// methods run before building a record / touching the DB.
    pub(crate) fn require_account(&self) -> Result<String, AppError> {
        self.current_account().ok_or(AppError::NotAuthenticated)
    }

    /// Adopt a **local `did:key` identity** as the active account: the device key *is* the identity,
    /// so AppView calls self-certify (`verify_signed` accepts `did:key` directly — no PDS record).
    /// This is the desktop bootstrap for the federated edge: device-key-signed calls (IBD suggestions,
    /// the encrypted exchange) work with no OAuth/PDS. Reuses an existing local identity if one is
    /// active; otherwise generates + persists a fresh device key. Returns the `did:key`.
    pub fn use_local_identity(&self) -> Result<String, AppError> {
        if let Some(did) = self.current_account() {
            if did.starts_with("did:key:") && DeviceKey::load(KEYCHAIN_SERVICE, &did)?.is_some() {
                return Ok(did);
            }
        }
        let key = DeviceKey::generate();
        let did = key.did_key();
        key.save(KEYCHAIN_SERVICE, &did)?;
        let _ = self.auth.tokens.set_active(&did);
        *self.auth.active.lock().unwrap() = Some(did.clone());
        Ok(did)
    }

    /// Switch the active account to an already-known DID (in-memory; the keychain marker too). For
    /// multi-identity flows — e.g. driving both sides of an exchange from one process.
    pub fn set_active_account(&self, did: &str) {
        let _ = self.auth.tokens.set_active(did);
        *self.auth.active.lock().unwrap() = Some(did.to_string());
    }

    /// Sign out: drop the active account and delete its stored session.
    pub async fn logout(&self) -> Result<(), AppError> {
        let did = self.auth.active.lock().unwrap().take();
        if let Some(did) = did {
            self.auth.tokens.delete(&did)?;
        }
        self.auth.tokens.clear_active()?;
        Ok(())
    }

    /// Build the resilient sync engine for the active account, loading its session from the
    /// keychain. Errors with [`AppError::NotAuthenticated`] when no one is signed in. The
    /// engine auto-refreshes on 401 and retries transient failures with backoff.
    pub(crate) fn sync_engine(&self) -> Result<AsyncSync, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let session = self.auth.tokens.load(&did)?.ok_or(AppError::NotAuthenticated)?;
        Ok(AsyncSync::new(
            self.auth.http.clone(),
            self.auth.tokens.clone(),
            session,
            RetryPolicy::default(),
            self.auth.online.clone(),
        ))
    }

    /// Whether the last PDS write reached the server. Drives the UI's offline indicator;
    /// optimistic (`true`) until a transient write failure.
    pub fn is_online(&self) -> bool {
        self.auth.online.load(Ordering::Relaxed)
    }
}
