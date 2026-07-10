//! Authenticated session + keychain token storage.
//!
//! A [`Session`] (access/refresh tokens, the DPoP key, DID, PDS) is stored as JSON under
//! `(service, account)` so it survives restarts and never touches the H2/SQLite workspace.
//! The DPoP key is persisted because DPoP-bound tokens are only usable with the same key
//! that minted them.
//!
//! Storage goes through [`crate::secret_store`], which is in-memory unless the production
//! binary opted into the OS keychain — so tests never reach the real one.

use serde::{Deserialize, Serialize};

use crate::error::SyncError;
use crate::secret_store;

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

/// Stores [`Session`]s in the keychain, keyed by account (typically the DID).
#[derive(Clone)]
pub struct TokenStore {
    service: String,
}

impl TokenStore {
    pub fn new(service: impl Into<String>) -> Self {
        TokenStore {
            service: service.into(),
        }
    }

    /// Remember `did` as the active account (so the next launch reloads its session).
    pub fn set_active(&self, did: &str) -> Result<(), SyncError> {
        secret_store::set(&self.service, ACTIVE_MARKER, did)
    }

    /// The active account's DID, or `None` if no one is signed in.
    pub fn active(&self) -> Result<Option<String>, SyncError> {
        secret_store::get(&self.service, ACTIVE_MARKER)
    }

    /// Forget the active account (sign-out); leaves any stored session untouched.
    pub fn clear_active(&self) -> Result<(), SyncError> {
        secret_store::delete(&self.service, ACTIVE_MARKER)
    }

    pub fn save(&self, account: &str, session: &Session) -> Result<(), SyncError> {
        secret_store::set(&self.service, account, &serde_json::to_string(session)?)
    }

    /// Load a session, or `None` if no entry exists for `account`.
    pub fn load(&self, account: &str) -> Result<Option<Session>, SyncError> {
        match secret_store::get(&self.service, account)? {
            Some(json) => Ok(Some(serde_json::from_str(&json)?)),
            None => Ok(None),
        }
    }

    pub fn delete(&self, account: &str) -> Result<(), SyncError> {
        secret_store::delete(&self.service, account)
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
        let store = TokenStore::new("navigator-sync-absent-account-test");
        assert!(store.load("did:plc:definitely-not-present").unwrap().is_none());
    }

    #[test]
    fn session_survives_a_save_load_round_trip() {
        let store = TokenStore::new("navigator-sync-round-trip-test");
        store.save("did:plc:abc123", &session()).unwrap();
        assert_eq!(store.load("did:plc:abc123").unwrap(), Some(session()));

        store.delete("did:plc:abc123").unwrap();
        assert!(store.load("did:plc:abc123").unwrap().is_none());
    }

    /// The active marker is a separate account name, so it can't be clobbered by — or clobber —
    /// a session stored under a DID.
    #[test]
    fn active_marker_is_independent_of_the_stored_session() {
        let store = TokenStore::new("navigator-sync-active-marker-test");
        assert!(store.active().unwrap().is_none(), "nobody signed in initially");

        store.save("did:plc:abc123", &session()).unwrap();
        store.set_active("did:plc:abc123").unwrap();
        assert_eq!(store.active().unwrap().as_deref(), Some("did:plc:abc123"));

        // Sign-out forgets who was active but leaves the session recoverable.
        store.clear_active().unwrap();
        assert!(store.active().unwrap().is_none());
        assert_eq!(store.load("did:plc:abc123").unwrap(), Some(session()));
    }
}
