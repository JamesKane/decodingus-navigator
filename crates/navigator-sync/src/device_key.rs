//! Per-device Ed25519 signing key for authenticated Edge↔AppView calls.
//!
//! Navigator cannot sign as its `did:plc` account (the PDS custodies that signing key), so
//! each installation generates its own Ed25519 *device key*, persists the 32-byte seed in
//! the OS keychain beside the OAuth session (keyed by account DID), and publishes the
//! public half once as a `com.decodingus.atmosphere.deviceKey` record in the user's PDS
//! repo. The AppView ingests that record via Jetstream and verifies every signed call
//! against the key — the *same* `du_atproto::signature::verify_did_key` code path proven by
//! the round-trip test below. Revocation = deleting the record.
//!
//! Wire contract (settled with the AppView team): the signature is base64-**STANDARD** of
//! the 64-byte Ed25519 signature over a canonical `\n`-joined UTF-8 message; the public key
//! travels as a `did:key:z…` identifier.

use base64::engine::general_purpose::STANDARD;
use base64::Engine;
use ed25519_dalek::{Signer, SigningKey};
use keyring::Entry;

use crate::error::SyncError;

/// Keychain account-name prefix for the device-key seed, namespaced so it can't collide
/// with a session entry (those are keyed by bare DID) or the active-account marker.
const DEVICE_KEY_PREFIX: &str = "__devicekey__";

/// PDS collection NSID for the published device-key record. Locked with the AppView team —
/// must match their Jetstream consumer exactly. The record value is
/// `{ "publicKey": "did:key:z…", "createdAt": "<rfc3339>" }`.
pub const DEVICE_KEY_COLLECTION: &str = "com.decodingus.atmosphere.deviceKey";

/// An installation's Ed25519 signing key, used to authenticate Edge→AppView calls.
#[derive(Clone)]
pub struct DeviceKey {
    signing: SigningKey,
}

impl DeviceKey {
    /// Generate a fresh random device key.
    pub fn generate() -> Self {
        let mut seed = [0u8; 32];
        rand_core::RngCore::fill_bytes(&mut rand_core::OsRng, &mut seed);
        DeviceKey { signing: SigningKey::from_bytes(&seed) }
    }

    /// Reconstruct from a stored 32-byte seed.
    fn from_seed(seed: &[u8]) -> Result<Self, SyncError> {
        let arr: [u8; 32] = seed
            .try_into()
            .map_err(|_| SyncError::Crypto("device key seed must be 32 bytes".into()))?;
        Ok(DeviceKey { signing: SigningKey::from_bytes(&arr) })
    }

    /// The `did:key:z…` identifier carrying this key's public half — what gets published in
    /// the PDS record and what the AppView verifies signed calls against. Reuses the shared
    /// `du-atproto` encoder (the inverse of its verify-side decoder), so the two can't drift.
    pub fn did_key(&self) -> String {
        du_atproto::did::did_key_from_ed25519(&self.signing.verifying_key())
    }

    /// Sign `message`, returning base64-**STANDARD** of the 64-byte signature — the exact
    /// wire format `du_atproto::signature::verify_did_key` (the AppView's verifier) expects.
    pub fn sign(&self, message: &str) -> String {
        STANDARD.encode(self.signing.sign(message.as_bytes()).to_bytes())
    }

    /// Deterministic PDS record key for this key's published record: the `did:key` multibase
    /// body (`z…`, base58btc — a valid record-key alphabet) with the `did:key:` scheme
    /// stripped. Stable per key, so re-publishing overwrites the same record (idempotent) and
    /// distinct devices get distinct records.
    pub fn record_rkey(&self) -> String {
        self.did_key().strip_prefix("did:key:").unwrap_or_default().to_string()
    }

    // --- keychain persistence (seed stored beside the OAuth session) ---

    fn account(did: &str) -> String {
        format!("{DEVICE_KEY_PREFIX}{did}")
    }

    /// Load this account's device key, or `None` if none has been generated yet.
    pub fn load(service: &str, did: &str) -> Result<Option<Self>, SyncError> {
        let entry = Entry::new(service, &Self::account(did))?;
        match entry.get_password() {
            Ok(seed_b64) => {
                let seed = STANDARD
                    .decode(seed_b64.trim())
                    .map_err(|e| SyncError::Crypto(format!("device key seed base64: {e}")))?;
                Ok(Some(Self::from_seed(&seed)?))
            }
            Err(keyring::Error::NoEntry) => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    /// Persist this key's 32-byte seed for `did` (base64 in the keychain).
    pub fn save(&self, service: &str, did: &str) -> Result<(), SyncError> {
        let seed_b64 = STANDARD.encode(self.signing.to_bytes());
        Entry::new(service, &Self::account(did))?.set_password(&seed_b64)?;
        Ok(())
    }

    /// Forget this account's device key (revocation companion: also delete the PDS record).
    pub fn delete(service: &str, did: &str) -> Result<(), SyncError> {
        let entry = Entry::new(service, &Self::account(did))?;
        match entry.delete_credential() {
            Ok(()) | Err(keyring::Error::NoEntry) => Ok(()),
            Err(e) => Err(e.into()),
        }
    }

    /// Load the device key for `did`, generating + persisting one on first use.
    pub fn load_or_generate(service: &str, did: &str) -> Result<Self, SyncError> {
        if let Some(key) = Self::load(service, did)? {
            return Ok(key);
        }
        let key = DeviceKey::generate();
        key.save(service, did)?;
        Ok(key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sign_verifies_against_the_appview_verifier() {
        // The whole point: sign here, verify with the *same code the AppView runs*.
        let key = DeviceKey::generate();
        let did = key.did_key();
        assert!(did.starts_with("did:key:z"), "did:key multibase prefix");

        let msg = "ibd-poll\ndid:plc:abc123\n1718000000";
        let sig = key.sign(msg);
        assert!(
            du_atproto::signature::verify_did_key(&did, msg.as_bytes(), &sig).is_ok(),
            "AppView verifier must accept our signature"
        );

        // Tampered message must be rejected.
        assert!(
            du_atproto::signature::verify_did_key(&did, b"tampered", &sig).is_err(),
            "verifier must reject a tampered message"
        );
    }

    #[test]
    fn seed_round_trips_through_base64() {
        let key = DeviceKey::generate();
        let seed_b64 = STANDARD.encode(key.signing.to_bytes());
        let restored = DeviceKey::from_seed(&STANDARD.decode(&seed_b64).unwrap()).unwrap();
        // Same key reconstructed → same public identity + same signatures.
        assert_eq!(key.did_key(), restored.did_key());
        assert_eq!(key.sign("ibd-poll\nx\n1"), restored.sign("ibd-poll\nx\n1"));
    }

    #[test]
    fn signature_is_standard_base64_not_url_safe() {
        // STANDARD alphabet uses '+' '/' '='; URL-safe uses '-' '_'. The AppView decodes
        // with STANDARD, so we must never emit the URL-safe-only characters.
        let sig = DeviceKey::generate().sign("hello\nworld\n42");
        assert!(!sig.contains('-') && !sig.contains('_'), "must be STANDARD, not URL-safe base64");
    }

    #[test]
    fn from_seed_rejects_wrong_length() {
        assert!(matches!(DeviceKey::from_seed(&[0u8; 16]), Err(SyncError::Crypto(_))));
    }

    #[test]
    fn record_rkey_is_did_key_without_scheme() {
        let key = DeviceKey::generate();
        let rkey = key.record_rkey();
        assert!(rkey.starts_with('z'), "multibase base58btc body");
        assert_eq!(format!("did:key:{rkey}"), key.did_key());
        // Valid record-key alphabet (no ':' from the scheme).
        assert!(!rkey.contains(':'));
    }
}
