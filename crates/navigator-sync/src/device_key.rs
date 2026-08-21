//! An Ed25519 signing key for each device, for authenticated Edge↔AppView calls.
//!
//! Navigator can not sign as its `did:plc` account, because the PDS holds that signing key. So
//! each installation makes its own Ed25519 *device key*. It keeps the 32-byte seed in the OS
//! keychain, beside the OAuth session, keyed by the account DID. It publishes the public half
//! once, as a `com.decodingus.atmosphere.deviceKey` record in the user's PDS repo.
//!
//! The AppView ingests that record through Jetstream, and checks every signed call against the
//! key. That is the *same* `du_atproto::signature::verify_did_key` code path that the round-trip
//! test below proves. To revoke the key, delete the record.
//!
//! The wire contract, which the AppView team agreed. The signature is base64-**STANDARD** of the
//! 64-byte Ed25519 signature. It covers a canonical UTF-8 message that `\n` joins. The public key
//! travels as a `did:key:z…` identifier.

use base64::engine::general_purpose::STANDARD;
use base64::Engine;
use ed25519_dalek::{Signer, SigningKey};

use crate::error::SyncError;
use crate::secret_store;

/// The keychain account-name prefix for the device-key seed. The prefix gives it its own space,
/// so it can never match a session entry, which keys on the bare DID, or the active-account
/// marker.
const DEVICE_KEY_PREFIX: &str = "__devicekey__";

/// PDS collection NSID for the published device-key record. The AppView team agreed it, and it
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
        DeviceKey {
            signing: SigningKey::from_bytes(&seed),
        }
    }

    /// Reconstruct from a stored 32-byte seed.
    fn from_seed(seed: &[u8]) -> Result<Self, SyncError> {
        let arr: [u8; 32] = seed
            .try_into()
            .map_err(|_| SyncError::Crypto("device key seed must be 32 bytes".into()))?;
        Ok(DeviceKey {
            signing: SigningKey::from_bytes(&arr),
        })
    }

    /// The `did:key:z…` identifier that holds this key's public half. The PDS record publishes it,
    /// and the AppView checks a signed call against it. It reuses the shared `du-atproto` encoder,
    /// which is the inverse of that crate's decoder, so the two can never become different.
    pub fn did_key(&self) -> String {
        du_atproto::did::did_key_from_ed25519(&self.signing.verifying_key())
    }

    /// Sign `message`, and return base64-**STANDARD** of the 64-byte signature. That is the exact
    /// wire format that `du_atproto::signature::verify_did_key`, the AppView's checker, expects.
    pub fn sign(&self, message: &str) -> String {
        STANDARD.encode(self.signing.sign(message.as_bytes()).to_bytes())
    }

    /// Sign an Edge request that **changes** something, in the frame with a replay guard that the
    /// AppView needs. That frame puts the caller's timestamp in front of the base message, as
    /// `{ts}\n{base}`.
    ///
    /// It mirrors `du_web::sig::fresh_message` exactly. The AppView refuses a stale `ts`, and it
    /// burns the signature, so nobody can use it twice. The caller sends `ts` in the request, with
    /// `did` and `signature`, and signs again with a new `ts` on a second try.
    pub fn sign_fresh(&self, ts: i64, base_message: &str) -> String {
        self.sign(&format!("{ts}\n{base_message}"))
    }

    /// Deterministic PDS record key for this key's published record. It is the `did:key`
    /// multibase body (`z…`, base58btc, which is a valid record-key alphabet), with the `did:key:`
    /// scheme removed. It is stable for each key. So a second publish overwrites the same record,
    /// which makes it idempotent, and two different devices get two different records.
    pub fn record_rkey(&self) -> String {
        self.did_key().strip_prefix("did:key:").unwrap_or_default().to_string()
    }

    // --- keychain persistence (seed stored beside the OAuth session) ---

    fn account(did: &str) -> String {
        format!("{DEVICE_KEY_PREFIX}{did}")
    }

    /// Load this account's device key. It gives `None` when no key exists yet.
    pub fn load(service: &str, did: &str) -> Result<Option<Self>, SyncError> {
        match secret_store::get(service, &Self::account(did))? {
            Some(seed_b64) => {
                let seed = STANDARD
                    .decode(seed_b64.trim())
                    .map_err(|e| SyncError::Crypto(format!("device key seed base64: {e}")))?;
                Ok(Some(Self::from_seed(&seed)?))
            }
            None => Ok(None),
        }
    }

    /// Persist this key's 32-byte seed for `did` (base64 in the keychain).
    pub fn save(&self, service: &str, did: &str) -> Result<(), SyncError> {
        secret_store::set(service, &Self::account(did), &STANDARD.encode(self.signing.to_bytes()))
    }

    /// Forget this account's device key (revocation companion: also delete the PDS record).
    pub fn delete(service: &str, did: &str) -> Result<(), SyncError> {
        secret_store::delete(service, &Self::account(did))
    }

    /// Load the device key for `did`. On the first use it makes one and stores it.
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
        // The whole point: sign here, and check with the *same code that the AppView runs*.
        let key = DeviceKey::generate();
        let did = key.did_key();
        assert!(did.starts_with("did:key:z"), "did:key multibase prefix");

        let msg = "ibd-poll\ndid:plc:abc123\n1718000000";
        let sig = key.sign(msg);
        assert!(
            du_atproto::signature::verify_did_key(&did, msg.as_bytes(), &sig).is_ok(),
            "AppView verifier must accept our signature"
        );

        // The check must refuse a message that somebody changed.
        assert!(
            du_atproto::signature::verify_did_key(&did, b"tampered", &sig).is_err(),
            "verifier must reject a tampered message"
        );
    }

    #[test]
    fn sign_fresh_frames_ts_and_verifies() {
        // `sign_fresh(ts, base)` signs exactly `{ts}\n{base}`. That is the frame that the
        // AppView's `verify_signed_fresh` builds. Prove it two ways. It must equal a signature over
        // the framed string by hand, and the AppView checker must accept it over that frame.
        let key = DeviceKey::generate();
        let did = key.did_key();
        let ts = 1_718_000_000_i64;
        let base = "exchange-consent\nurn:x\ndid:key:zA\ntrue";
        let framed = format!("{ts}\n{base}");
        let sig = key.sign_fresh(ts, base);
        assert_eq!(sig, key.sign(&framed), "sign_fresh must frame as {{ts}}\\n{{base}}");
        assert!(
            du_atproto::signature::verify_did_key(&did, framed.as_bytes(), &sig).is_ok(),
            "AppView verifier must accept the framed signature"
        );
        // A signature over the bare base, with no frame, must NOT pass against the framed bytes.
        let bare = key.sign(base);
        assert!(du_atproto::signature::verify_did_key(&did, framed.as_bytes(), &bare).is_err());
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
        assert!(
            !sig.contains('-') && !sig.contains('_'),
            "must be STANDARD, not URL-safe base64"
        );
    }

    #[test]
    fn from_seed_rejects_wrong_length() {
        assert!(matches!(DeviceKey::from_seed(&[0u8; 16]), Err(SyncError::Crypto(_))));
    }

    /// First call generates + persists; every later call must return the *same* identity, or the
    /// AppView would reject signatures made by the regenerated key.
    #[test]
    fn load_or_generate_persists_and_is_stable() {
        let (svc, did) = ("device-key-stable-test", "did:plc:abc123");
        assert!(DeviceKey::load(svc, did).unwrap().is_none(), "nothing stored yet");

        let first = DeviceKey::load_or_generate(svc, did).unwrap();
        let second = DeviceKey::load_or_generate(svc, did).unwrap();
        assert_eq!(first.did_key(), second.did_key(), "must reuse the persisted seed");
        assert_eq!(first.sign("m"), second.sign("m"));

        DeviceKey::delete(svc, did).unwrap();
        assert!(DeviceKey::load(svc, did).unwrap().is_none(), "revoked");
        // A key generated after revocation is a genuinely new identity.
        assert_ne!(
            DeviceKey::load_or_generate(svc, did).unwrap().did_key(),
            first.did_key()
        );
    }

    /// Each account gives its keys their own space, so two identities on one installation never
    /// share a name.
    #[test]
    fn device_keys_are_per_account() {
        let svc = "device-key-per-account-test";
        let a = DeviceKey::load_or_generate(svc, "did:plc:aaa").unwrap();
        let b = DeviceKey::load_or_generate(svc, "did:plc:bbb").unwrap();
        assert_ne!(a.did_key(), b.did_key());
        assert_eq!(
            DeviceKey::load(svc, "did:plc:aaa").unwrap().unwrap().did_key(),
            a.did_key()
        );
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
