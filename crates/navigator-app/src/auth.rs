//! `impl App` methods extracted from `lib.rs` (the `auth` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- authentication ----------------------------------------------------

    /// Do the public-client OAuth login for `handle`. The value can be a handle or a DID.
    ///
    /// The sequence is: the browser authorizes, the loopback receives the callback, and the client
    /// exchanges the token. After a good login, the code writes the DPoP-bound session to the OS
    /// keychain, and that session becomes the active account. The function returns the DID.
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

    /// The DID of the active account, or [`AppError::NotAuthenticated`]. This is the small guard
    /// that the publish methods call first, before they make a record or read the database.
    pub(crate) fn require_account(&self) -> Result<String, AppError> {
        self.current_account().ok_or(AppError::NotAuthenticated)
    }

    /// Use a **local `did:key` identity** as the active account. The device key is the identity.
    /// So a call to the AppView certifies itself, because `verify_signed` accepts a `did:key`
    /// directly and needs no PDS record.
    ///
    /// This function is the desktop start point for the federated edge. Calls that the device key
    /// signs, such as IBD suggestions and the encrypted exchange, then work with no OAuth and no
    /// PDS. If a local identity is already active, this function uses it. If not, it makes a new
    /// device key and writes it to the keychain. It returns the `did:key`.
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

    /// Change the active account to a known DID. The change applies to memory and to the keychain
    /// marker. Use this for a flow with more than one identity. One example is a test that operates
    /// both sides of an exchange in one process.
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

    /// Make the sync engine for the active account and read its session from the keychain.
    ///
    /// The function returns [`AppError::NotAuthenticated`] if no account is active. On a 401
    /// response, the engine refreshes the token. It also tries again after a temporary failure.
    /// The delay becomes longer after each try.
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
