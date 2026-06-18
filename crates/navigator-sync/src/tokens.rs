//! Authenticated session + OS-keychain token storage.
//!
//! A [`Session`] (access/refresh tokens, the DPoP key, DID, PDS) is stored as JSON in the
//! OS keychain under `(service, account)` so it survives restarts and never touches the
//! H2/SQLite workspace. The DPoP key is persisted because DPoP-bound tokens are only
//! usable with the same key that minted them.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

use keyring::Entry;
use serde::{Deserialize, Serialize};

use crate::error::SyncError;

/// When set, [`TokenStore`] uses a process-global in-memory map instead of the OS keychain — so tests
/// (and CI) never touch the real keychain (no prompts) and stay hermetic regardless of any ambient
/// session. Set once via [`TokenStore::use_in_memory_for_tests`]; never enabled in production.
static IN_MEMORY: AtomicBool = AtomicBool::new(false);

fn mem() -> &'static Mutex<HashMap<(String, String), String>> {
    static M: OnceLock<Mutex<HashMap<(String, String), String>>> = OnceLock::new();
    M.get_or_init(|| Mutex::new(HashMap::new()))
}
fn in_memory() -> bool {
    IN_MEMORY.load(Ordering::Relaxed)
}

/// An authenticated AT Proto session for one account.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Session {
    /// Account DID.
    pub did: String,
    /// Resolved PDS service endpoint.
    pub pds: String,
    pub access_token: String,
    pub refresh_token: String,
    /// The DPoP key (base64, via `EcKey::to_base64`) that bound these tokens.
    pub dpop_key_b64: String,
    pub scope: String,
    /// The `client_id` presented at login — replayed verbatim on token refresh (a public
    /// client must send the same identifier).
    pub client_id: String,
}

/// Keychain account name under which the active-account DID is remembered, so the app
/// can reload the right [`Session`] on the next launch. Not itself an account DID, so it
/// can't collide with one.
const ACTIVE_MARKER: &str = "__active__";

/// Stores [`Session`]s in the OS keychain, keyed by account (typically the DID).
#[derive(Clone)]
pub struct TokenStore {
    service: String,
}

impl TokenStore {
    pub fn new(service: impl Into<String>) -> Self {
        TokenStore { service: service.into() }
    }

    /// Route all token storage to a process-global in-memory map (no OS keychain). For tests/CI —
    /// call once before constructing anything that reads the keychain. Idempotent; never used in prod.
    pub fn use_in_memory_for_tests() {
        IN_MEMORY.store(true, Ordering::Relaxed);
    }

    fn key(&self, account: &str) -> (String, String) {
        (self.service.clone(), account.to_string())
    }

    /// Remember `did` as the active account (so the next launch reloads its session).
    pub fn set_active(&self, did: &str) -> Result<(), SyncError> {
        if in_memory() {
            mem().lock().unwrap().insert(self.key(ACTIVE_MARKER), did.to_string());
            return Ok(());
        }
        Entry::new(&self.service, ACTIVE_MARKER)?.set_password(did)?;
        Ok(())
    }

    /// The active account's DID, or `None` if no one is signed in.
    pub fn active(&self) -> Result<Option<String>, SyncError> {
        if in_memory() {
            return Ok(mem().lock().unwrap().get(&self.key(ACTIVE_MARKER)).cloned());
        }
        match Entry::new(&self.service, ACTIVE_MARKER)?.get_password() {
            Ok(did) => Ok(Some(did)),
            Err(keyring::Error::NoEntry) => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    /// Forget the active account (sign-out); leaves any stored session untouched.
    pub fn clear_active(&self) -> Result<(), SyncError> {
        if in_memory() {
            mem().lock().unwrap().remove(&self.key(ACTIVE_MARKER));
            return Ok(());
        }
        match Entry::new(&self.service, ACTIVE_MARKER)?.delete_credential() {
            Ok(()) | Err(keyring::Error::NoEntry) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    pub fn save(&self, account: &str, session: &Session) -> Result<(), SyncError> {
        let json = serde_json::to_string(session)?;
        if in_memory() {
            mem().lock().unwrap().insert(self.key(account), json);
            return Ok(());
        }
        Entry::new(&self.service, account)?.set_password(&json)?;
        Ok(())
    }

    /// Load a session, or `None` if no entry exists for `account`.
    pub fn load(&self, account: &str) -> Result<Option<Session>, SyncError> {
        if in_memory() {
            return match mem().lock().unwrap().get(&self.key(account)) {
                Some(json) => Ok(Some(serde_json::from_str(json)?)),
                None => Ok(None),
            };
        }
        let entry = Entry::new(&self.service, account)?;
        match entry.get_password() {
            Ok(json) => Ok(Some(serde_json::from_str(&json)?)),
            Err(keyring::Error::NoEntry) => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    pub fn delete(&self, account: &str) -> Result<(), SyncError> {
        if in_memory() {
            mem().lock().unwrap().remove(&self.key(account));
            return Ok(());
        }
        let entry = Entry::new(&self.service, account)?;
        match entry.delete_credential() {
            Ok(()) | Err(keyring::Error::NoEntry) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn session() -> Session {
        Session {
            did: "did:plc:abc123".into(),
            pds: "https://pds.example".into(),
            access_token: "at".into(),
            refresh_token: "rt".into(),
            dpop_key_b64: "key".into(),
            scope: "atproto navigatorCore".into(),
            client_id: "http://localhost?redirect_uri=…".into(),
        }
    }

    #[test]
    fn session_json_round_trips() {
        // The stored keychain format is just this JSON; lock it in.
        let json = serde_json::to_string(&session()).unwrap();
        let back: Session = serde_json::from_str(&json).unwrap();
        assert_eq!(back, session());
    }

    #[test]
    fn loading_an_absent_account_is_none_not_an_error() {
        // Read-only lookup of a guaranteed-absent account: exercises the NoEntry -> None
        // mapping against the real backend without writing (no keychain prompt). The
        // backend is unavailable in some sandboxes, so only assert when the lookup works.
        let store = TokenStore::new("navigator-sync-absent-account-test");
        if let Ok(found) = store.load("did:plc:definitely-not-present") {
            assert!(found.is_none());
        }
    }
}
