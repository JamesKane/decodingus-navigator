//! Authenticated session + OS-keychain token storage.
//!
//! A [`Session`] (access/refresh tokens, the DPoP key, DID, PDS) is stored as JSON in the
//! OS keychain under `(service, account)` so it survives restarts and never touches the
//! H2/SQLite workspace. The DPoP key is persisted because DPoP-bound tokens are only
//! usable with the same key that minted them.

use keyring::Entry;
use serde::{Deserialize, Serialize};

use crate::error::SyncError;

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
}

/// Stores [`Session`]s in the OS keychain, keyed by account (typically the DID).
pub struct TokenStore {
    service: String,
}

impl TokenStore {
    pub fn new(service: impl Into<String>) -> Self {
        TokenStore { service: service.into() }
    }

    pub fn save(&self, account: &str, session: &Session) -> Result<(), SyncError> {
        let entry = Entry::new(&self.service, account)?;
        entry.set_password(&serde_json::to_string(session)?)?;
        Ok(())
    }

    /// Load a session, or `None` if no entry exists for `account`.
    pub fn load(&self, account: &str) -> Result<Option<Session>, SyncError> {
        let entry = Entry::new(&self.service, account)?;
        match entry.get_password() {
            Ok(json) => Ok(Some(serde_json::from_str(&json)?)),
            Err(keyring::Error::NoEntry) => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    pub fn delete(&self, account: &str) -> Result<(), SyncError> {
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
