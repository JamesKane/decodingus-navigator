//! The one and only gateway to the OS keychain.
//!
//! Three things need durable secret storage: the OAuth [`Session`](crate::Session), the Ed25519
//! [`DeviceKey`](crate::DeviceKey), and the X25519 [`ExchangeKey`](crate::ExchangeKey). Each one
//! used to call `keyring::Entry` directly. So the question "does this code path touch the real
//! keychain?" had three answers, and no one place could control them. It has one place now: every
//! read, write, and delete goes through [`get`], [`set`], and [`delete`] here.
//!
//! **The backend is in memory, until a process asks for the real one.** So a test binary, a CI
//! runner, a doctest, and an `examples/` probe *can not* reach the login keychain. That holds
//! whatever they build, and in whatever order. The capability is off. The production binary turns
//! it on once, at the top of `main`, with [`use_os_keychain`].
//!
//! This is the safe direction for the default to fail. Under the old opt-*out* scheme, a test had
//! to remember to call an escape hatch. If it did not, it read the user's real credentials under
//! the production service name, with no warning. On macOS it also raised an unlock prompt that
//! hangs CI.
//!
//! Under this scheme, the same mistake means a session does not survive a restart. That is loud,
//! it is local to the one binary that owns `main`, and nobody can miss it at the first launch.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

use keyring::Entry;

use crate::error::SyncError;

/// Off by default: everything stays in memory until [`use_os_keychain`] sets this. Only the `main`
/// of the shipped binary may set it. A test must not, and a library initializer must not.
static OS_KEYCHAIN: AtomicBool = AtomicBool::new(false);

/// Process-global stand-in for the keychain: `(service, account) -> secret`.
fn mem() -> &'static Mutex<HashMap<(String, String), String>> {
    static M: OnceLock<Mutex<HashMap<(String, String), String>>> = OnceLock::new();
    M.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Route secret storage to the real OS keychain. Call it once from the `main` of the production
/// binary, before the code builds any `App`. It is idempotent.
///
/// Nothing else may call this. A test that calls it reaches for the developer's login keychain.
pub fn use_os_keychain() {
    OS_KEYCHAIN.store(true, Ordering::Relaxed);
}

/// Whether secrets are going to the real keychain. Exposed so tests can assert they are *not*.
pub fn os_keychain_enabled() -> bool {
    OS_KEYCHAIN.load(Ordering::Relaxed)
}

fn key(service: &str, account: &str) -> (String, String) {
    (service.to_string(), account.to_string())
}

/// The secret stored for `(service, account)`, or `None` if there is no entry.
pub(crate) fn get(service: &str, account: &str) -> Result<Option<String>, SyncError> {
    if !os_keychain_enabled() {
        return Ok(mem().lock().unwrap().get(&key(service, account)).cloned());
    }
    match Entry::new(service, account)?.get_password() {
        Ok(secret) => Ok(Some(secret)),
        Err(keyring::Error::NoEntry) => Ok(None),
        Err(e) => Err(e.into()),
    }
}

/// Store `secret` under `(service, account)`. It replaces an entry that is already there.
pub(crate) fn set(service: &str, account: &str, secret: &str) -> Result<(), SyncError> {
    if !os_keychain_enabled() {
        mem().lock().unwrap().insert(key(service, account), secret.to_string());
        return Ok(());
    }
    Entry::new(service, account)?.set_password(secret)?;
    Ok(())
}

/// Remove the entry for `(service, account)`. An absent entry counts as success, because the
/// delete is idempotent.
pub(crate) fn delete(service: &str, account: &str) -> Result<(), SyncError> {
    if !os_keychain_enabled() {
        mem().lock().unwrap().remove(&key(service, account));
        return Ok(());
    }
    match Entry::new(service, account)?.delete_credential() {
        Ok(()) | Err(keyring::Error::NoEntry) => Ok(()),
        Err(e) => Err(e.into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The one rule that this module must keep. If this test ever fails, every test in the
    /// workspace reads and writes the developer's real login keychain.
    #[test]
    fn os_keychain_is_off_unless_a_binary_opts_in() {
        assert!(
            !os_keychain_enabled(),
            "the OS keychain must stay off in tests — only the shipped binary's main() may enable it"
        );
    }

    #[test]
    fn absent_entry_reads_as_none_not_an_error() {
        assert_eq!(get("svc-absent", "nobody").unwrap(), None);
    }

    #[test]
    fn set_then_get_round_trips_and_delete_is_idempotent() {
        set("svc-rt", "acct", "s3cret").unwrap();
        assert_eq!(get("svc-rt", "acct").unwrap().as_deref(), Some("s3cret"));

        set("svc-rt", "acct", "rotated").unwrap();
        assert_eq!(get("svc-rt", "acct").unwrap().as_deref(), Some("rotated"));

        delete("svc-rt", "acct").unwrap();
        assert_eq!(get("svc-rt", "acct").unwrap(), None);
        delete("svc-rt", "acct").unwrap(); // deleting an absent entry is not an error
    }

    #[test]
    fn entries_are_namespaced_by_service_and_account() {
        set("svc-a", "acct", "a").unwrap();
        set("svc-b", "acct", "b").unwrap();
        set("svc-a", "other", "c").unwrap();
        assert_eq!(get("svc-a", "acct").unwrap().as_deref(), Some("a"));
        assert_eq!(get("svc-b", "acct").unwrap().as_deref(), Some("b"));
        assert_eq!(get("svc-a", "other").unwrap().as_deref(), Some("c"));
    }
}
