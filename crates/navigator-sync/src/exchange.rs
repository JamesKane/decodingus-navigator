//! Encrypted edge-to-edge exchange (IBD Phase 2, AppView design **D1**).
//!
//! The AppView brokers the discovery and the consent, and relays *opaque* ciphertext. It never
//! holds a decryption key. Two edges that both consent do four things:
//!   1. publish a static X25519 **identity key** (IK), signed by their Ed25519 device key;
//!   2. make an ephemeral X25519 key (EK) for each session, and exchange the public halves;
//!   3. derive a shared session key with **X3DH-lite**, a triple ECDH into HKDF-SHA-256;
//!   4. seal a payload with **AES-256-GCM**, and relay the ciphertext through the AppView.
//!
//! The AppView fixes the broker contract, which is the endpoints and the signed-message formats.
//! The **envelope and the key derivation** here belong to the edges alone, and the broker never
//! parses them. `EXCHANGE_VERSION` gives them a version, so both edges agree.

use aes_gcm::aead::{Aead, KeyInit, Payload};
use aes_gcm::{Aes256Gcm, Nonce};
use base64::engine::general_purpose::STANDARD;
use base64::Engine;
use hkdf::Hkdf;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use x25519_dalek::{PublicKey, StaticSecret};

use crate::error::SyncError;
use crate::secret_store;

/// Edge-convention version for the envelope + key derivation (bump on any change to either).
pub const EXCHANGE_VERSION: u8 = 1;
/// HKDF salt for the session key (edge convention).
const HKDF_SALT: &[u8] = b"decodingus/exchange/v1/salt";
/// HKDF info for the session key (edge convention).
const HKDF_INFO: &[u8] = b"decodingus/exchange/v1/session";
/// Keychain account-name prefix for the X25519 identity-key secret (namespaced like the device key).
const X25519_PREFIX: &str = "__x25519__";

fn random_32() -> [u8; 32] {
    let mut b = [0u8; 32];
    rand_core::RngCore::fill_bytes(&mut rand_core::OsRng, &mut b);
    b
}

fn decode_pub(b64: &str) -> Result<PublicKey, SyncError> {
    let bytes = STANDARD
        .decode(b64.trim())
        .map_err(|_| SyncError::Crypto("x25519 public-key base64".into()))?;
    let arr: [u8; 32] = bytes
        .as_slice()
        .try_into()
        .map_err(|_| SyncError::Crypto("x25519 public key must be 32 bytes".into()))?;
    Ok(PublicKey::from(arr))
}

/// A stored static X25519 **identity key** (IK) for an account. Each installation makes one, once.
/// The OS keychain holds it, beside the device key, and the app publishes its public half to the
/// AppView.
#[derive(Clone)]
pub struct ExchangeKey {
    secret: StaticSecret,
}

impl ExchangeKey {
    /// Generate a fresh identity key.
    pub fn generate() -> Self {
        ExchangeKey {
            secret: StaticSecret::from(random_32()),
        }
    }

    fn from_bytes(b: &[u8]) -> Result<Self, SyncError> {
        let arr: [u8; 32] = b
            .try_into()
            .map_err(|_| SyncError::Crypto("x25519 secret must be 32 bytes".into()))?;
        Ok(ExchangeKey {
            secret: StaticSecret::from(arr),
        })
    }

    fn account(did: &str) -> String {
        format!("{X25519_PREFIX}{did}")
    }

    /// Load the stored identity key for `did`, if present.
    pub fn load(service: &str, did: &str) -> Result<Option<Self>, SyncError> {
        match secret_store::get(service, &Self::account(did))? {
            Some(b64) => {
                let bytes = STANDARD
                    .decode(b64.trim())
                    .map_err(|e| SyncError::Crypto(e.to_string()))?;
                Ok(Some(Self::from_bytes(&bytes)?))
            }
            None => Ok(None),
        }
    }

    /// Persist this identity key for `did`.
    pub fn save(&self, service: &str, did: &str) -> Result<(), SyncError> {
        secret_store::set(service, &Self::account(did), &STANDARD.encode(self.secret.to_bytes()))
    }

    /// Load the identity key for `did`. On the first use it makes one and stores it.
    pub fn load_or_generate(service: &str, did: &str) -> Result<Self, SyncError> {
        if let Some(k) = Self::load(service, did)? {
            return Ok(k);
        }
        let k = Self::generate();
        k.save(service, did)?;
        Ok(k)
    }

    /// The public half, as STANDARD base64 of 32 raw bytes. That is the wire form that the AppView
    /// stores.
    pub fn public_b64(&self) -> String {
        STANDARD.encode(PublicKey::from(&self.secret).as_bytes())
    }
}

/// An ephemeral X25519 key for one session. It gives the forward secrecy of X3DH-lite, and the
/// code discards it when the session ends.
pub struct EphemeralKey {
    secret: StaticSecret,
}

impl Default for EphemeralKey {
    fn default() -> Self {
        Self::generate()
    }
}

impl EphemeralKey {
    pub fn generate() -> Self {
        EphemeralKey {
            secret: StaticSecret::from(random_32()),
        }
    }

    /// The ephemeral public half as STANDARD base64 (sent in the handshake envelope).
    pub fn public_b64(&self) -> String {
        STANDARD.encode(PublicKey::from(&self.secret).as_bytes())
    }
}

/// Whether this edge is party "A" in the canonical X3DH-lite order. The two DIDs decide it alone,
/// in lexicographic order. So both sides agree, and neither one finds out who started.
pub fn role_is_a(my_did: &str, partner_did: &str) -> bool {
    my_did < partner_did
}

/// Derive the 32-byte AES session key with X3DH-lite. The canonical secret is
/// `DH(IK_A,EK_B) ‖ DH(EK_A,IK_B) ‖ DH(EK_A,EK_B)`, through HKDF-SHA-256. `i_am_a` selects which
/// of my keys take the A role and which take the B role, see [`role_is_a`]. Both edges compute the
/// same key.
pub fn derive_session_key(
    my_ik: &ExchangeKey,
    my_ek: &EphemeralKey,
    their_ik_b64: &str,
    their_ek_b64: &str,
    i_am_a: bool,
) -> Result<[u8; 32], SyncError> {
    let their_ik = decode_pub(their_ik_b64)?;
    let their_ek = decode_pub(their_ek_b64)?;
    // The three canonical DH values, computed from whichever halves this edge holds.
    let (d1, d2, d3) = if i_am_a {
        // A holds IK_A, EK_A (mine) and IK_B, EK_B (theirs).
        (
            my_ik.secret.diffie_hellman(&their_ek),
            my_ek.secret.diffie_hellman(&their_ik),
            my_ek.secret.diffie_hellman(&their_ek),
        )
    } else {
        // B holds IK_B, EK_B (mine) and IK_A, EK_A (theirs): same three values, mirrored.
        (
            my_ek.secret.diffie_hellman(&their_ik),
            my_ik.secret.diffie_hellman(&their_ek),
            my_ek.secret.diffie_hellman(&their_ek),
        )
    };
    let mut ikm = Vec::with_capacity(96);
    ikm.extend_from_slice(d1.as_bytes());
    ikm.extend_from_slice(d2.as_bytes());
    ikm.extend_from_slice(d3.as_bytes());
    let hk = Hkdf::<Sha256>::new(Some(HKDF_SALT), &ikm);
    let mut key = [0u8; 32];
    hk.expand(HKDF_INFO, &mut key)
        .map_err(|_| SyncError::Crypto("hkdf expand".into()))?;
    Ok(key)
}

/// An exchange envelope: the payload, in a cleartext frame, that becomes the relay `blob`. A
/// `Handshake` carries the sender's ephemeral public key. A `Data` carries an AES-256-GCM
/// ciphertext.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "t", rename_all = "snake_case")]
pub enum Envelope {
    /// Session handshake: the sender's ephemeral X25519 public key (STANDARD base64).
    Handshake { v: u8, ek: String },
    /// Encrypted data: GCM nonce + ciphertext (both STANDARD base64).
    Data { v: u8, iv: String, ct: String },
}

impl Envelope {
    /// The handshake envelope that gives out `ek`'s public half.
    pub fn handshake(ek: &EphemeralKey) -> Self {
        Envelope::Handshake {
            v: EXCHANGE_VERSION,
            ek: ek.public_b64(),
        }
    }

    /// The relay `blob`: STANDARD base64 of the JSON-serialized envelope.
    pub fn to_blob(&self) -> Result<String, SyncError> {
        let bytes = serde_json::to_vec(self).map_err(|e| SyncError::Crypto(e.to_string()))?;
        Ok(STANDARD.encode(bytes))
    }

    /// Parse a relay `blob` back into an envelope.
    pub fn from_blob(blob_b64: &str) -> Result<Self, SyncError> {
        let bytes = STANDARD
            .decode(blob_b64.trim())
            .map_err(|e| SyncError::Crypto(e.to_string()))?;
        serde_json::from_slice(&bytes).map_err(|e| SyncError::Crypto(e.to_string()))
    }
}

/// The STANDARD-base64 SHA-256 the AppView signs over for a relay (it hashes the *decoded* blob
/// bytes). Matches the server: `STANDARD(SHA256(STANDARD::decode(blob)))`.
pub fn blob_sha256_b64(blob_b64: &str) -> Result<String, SyncError> {
    let bytes = STANDARD
        .decode(blob_b64.trim())
        .map_err(|e| SyncError::Crypto(e.to_string()))?;
    Ok(STANDARD.encode(Sha256::digest(&bytes)))
}

/// The AES-GCM additional authenticated data, which ties an envelope to its relay metadata.
/// Nobody can then replay a blob into another session, or change where it goes.
pub fn relay_aad(session_id: &str, from_did: &str, to_did: &str, seq: i32) -> String {
    format!("{session_id}\n{from_did}\n{to_did}\n{seq}")
}

/// Seal `plaintext` into a `Data` envelope under `session_key`, and authenticate `aad`.
// `Nonce::from_slice` routes through generic-array's `from_slice`, which the resolved transitive
// generic-array marks deprecated (a 0.14→1.x churn artifact, not a real deprecation for aes-gcm
// 0.10). The call is correct; silence the lint locally.
#[allow(deprecated)]
pub fn seal(session_key: &[u8; 32], aad: &str, plaintext: &[u8]) -> Result<Envelope, SyncError> {
    let cipher = Aes256Gcm::new_from_slice(session_key).map_err(|_| SyncError::Crypto("aes-gcm key".into()))?;
    let mut iv = [0u8; 12];
    rand_core::RngCore::fill_bytes(&mut rand_core::OsRng, &mut iv);
    let ct = cipher
        .encrypt(
            Nonce::from_slice(&iv),
            Payload {
                msg: plaintext,
                aad: aad.as_bytes(),
            },
        )
        .map_err(|_| SyncError::Crypto("aes-gcm seal".into()))?;
    Ok(Envelope::Data {
        v: EXCHANGE_VERSION,
        iv: STANDARD.encode(iv),
        ct: STANDARD.encode(ct),
    })
}

/// Open a `Data` envelope under `session_key`, and check `aad`. It gives an error for an envelope
/// that is not data, and for a failed authentication. That failure means a changed ciphertext, a
/// wrong key, or a wrong aad.
#[allow(deprecated)] // see `seal` — generic-array `from_slice` churn lint
pub fn open(session_key: &[u8; 32], aad: &str, env: &Envelope) -> Result<Vec<u8>, SyncError> {
    let Envelope::Data { iv, ct, .. } = env else {
        return Err(SyncError::Crypto("not a data envelope".into()));
    };
    let iv = STANDARD
        .decode(iv.trim())
        .map_err(|e| SyncError::Crypto(e.to_string()))?;
    let ct = STANDARD
        .decode(ct.trim())
        .map_err(|e| SyncError::Crypto(e.to_string()))?;
    let cipher = Aes256Gcm::new_from_slice(session_key).map_err(|_| SyncError::Crypto("aes-gcm key".into()))?;
    cipher
        .decrypt(
            Nonce::from_slice(&iv),
            Payload {
                msg: &ct,
                aad: aad.as_bytes(),
            },
        )
        .map_err(|_| SyncError::Crypto("aes-gcm open (auth failed)".into()))
}

/// The canonical signed-message builders. They mirror the AppView's `du_db::exchange::messages`
/// **exactly**, and a check against the server source proved that. The device key signs them.
pub mod messages {
    /// `exchange-publickey\n{did}\n{x25519_pub_b64}\n{key_uri_or_empty}`
    pub fn publickey(did: &str, x25519_pub_b64: &str, key_uri: Option<&str>) -> String {
        format!("exchange-publickey\n{did}\n{x25519_pub_b64}\n{}", key_uri.unwrap_or(""))
    }
    /// `exchange-request\n{request_uri}\n{initiator_did}\n{partner_did}\n{purpose}\n{scope_or_empty}`
    pub fn request(
        request_uri: &str,
        initiator_did: &str,
        partner_did: &str,
        purpose: &str,
        scope: Option<&str>,
    ) -> String {
        format!(
            "exchange-request\n{request_uri}\n{initiator_did}\n{partner_did}\n{purpose}\n{}",
            scope.unwrap_or("")
        )
    }
    /// `exchange-consent\n{request_uri}\n{consenting_did}\n{given}`
    pub fn consent(request_uri: &str, consenting_did: &str, given: bool) -> String {
        format!("exchange-consent\n{request_uri}\n{consenting_did}\n{given}")
    }
    /// `exchange-poll\n{did}\n{ts}`: the GET auth for these endpoints: the one that lists what
    /// arrived, the one that lists what waits, and the relay pull.
    pub fn poll(did: &str, ts: i64) -> String {
        format!("exchange-poll\n{did}\n{ts}")
    }
    /// `exchange-relay\n{session_id}\n{from_did}\n{to_did}\n{seq}\n{blob_sha256_b64}`
    pub fn relay(session_id: &str, from_did: &str, to_did: &str, seq: i32, blob_sha256_b64: &str) -> String {
        format!("exchange-relay\n{session_id}\n{from_did}\n{to_did}\n{seq}\n{blob_sha256_b64}")
    }
    /// `exchange-ack\n{did}\n{envelope_id}`
    pub fn ack(did: &str, envelope_id: i64) -> String {
        format!("exchange-ack\n{did}\n{envelope_id}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_edges_derive_the_same_session_key() {
        // Alice (A) and Bob (B). Each one has a static IK, and an EK for the session.
        let (a_did, b_did) = ("did:plc:aaa", "did:plc:bbb");
        let a_ik = ExchangeKey::generate();
        let b_ik = ExchangeKey::generate();
        let a_ek = EphemeralKey::generate();
        let b_ek = EphemeralKey::generate();
        // The DID order decides the roles, the same way on both sides.
        let a_is_a = role_is_a(a_did, b_did);
        let b_is_a = role_is_a(b_did, a_did);
        assert!(a_is_a && !b_is_a);

        let ka = derive_session_key(&a_ik, &a_ek, &b_ik.public_b64(), &b_ek.public_b64(), a_is_a).unwrap();
        let kb = derive_session_key(&b_ik, &b_ek, &a_ik.public_b64(), &a_ek.public_b64(), b_is_a).unwrap();
        assert_eq!(ka, kb, "X3DH-lite must yield the same key on both edges");
    }

    #[test]
    fn seal_open_round_trip_and_tamper_rejection() {
        let key = [7u8; 32];
        let aad = relay_aad("sess-1", "did:plc:aaa", "did:plc:bbb", 1);
        let env = seal(&key, &aad, b"variant positions").unwrap();
        assert_eq!(open(&key, &aad, &env).unwrap(), b"variant positions");

        // Wrong AAD (e.g. replayed under a different seq) fails authentication.
        let other_aad = relay_aad("sess-1", "did:plc:aaa", "did:plc:bbb", 2);
        assert!(open(&key, &other_aad, &env).is_err());
        // Wrong key fails.
        assert!(open(&[9u8; 32], &aad, &env).is_err());
    }

    #[test]
    fn full_handshake_then_encrypted_exchange() {
        // End to end, in one process. Derive on both sides from the handshakes that they
        // exchanged, then send data from A to B.
        let (a_did, b_did) = ("did:plc:aaa", "did:plc:bbb");
        let (a_ik, b_ik) = (ExchangeKey::generate(), ExchangeKey::generate());
        let (a_ek, b_ek) = (EphemeralKey::generate(), EphemeralKey::generate());

        // Each posts a handshake blob; the other parses out the ephemeral public key.
        let a_hs = Envelope::handshake(&a_ek).to_blob().unwrap();
        let b_hs = Envelope::handshake(&b_ek).to_blob().unwrap();
        let a_sees_b_ek = match Envelope::from_blob(&b_hs).unwrap() {
            Envelope::Handshake { ek, .. } => ek,
            _ => panic!("expected handshake"),
        };
        let b_sees_a_ek = match Envelope::from_blob(&a_hs).unwrap() {
            Envelope::Handshake { ek, .. } => ek,
            _ => panic!("expected handshake"),
        };

        let ka = derive_session_key(&a_ik, &a_ek, &b_ik.public_b64(), &a_sees_b_ek, role_is_a(a_did, b_did)).unwrap();
        let kb = derive_session_key(&b_ik, &b_ek, &a_ik.public_b64(), &b_sees_a_ek, role_is_a(b_did, a_did)).unwrap();
        assert_eq!(ka, kb);

        // A seals a data envelope → blob; the relay-signed hash matches the decoded bytes.
        let aad = relay_aad("sess-1", a_did, b_did, 1);
        let blob = seal(&ka, &aad, b"hello bob").unwrap().to_blob().unwrap();
        let expected_hash = STANDARD.encode(Sha256::digest(STANDARD.decode(&blob).unwrap()));
        assert_eq!(blob_sha256_b64(&blob).unwrap(), expected_hash);

        // B pulls the blob, parses, and opens it with the shared key + same AAD.
        let env = Envelope::from_blob(&blob).unwrap();
        assert_eq!(open(&kb, &aad, &env).unwrap(), b"hello bob");
    }

    #[test]
    fn canonical_messages_match_the_appview() {
        // These strings are the cross-repo contract (du_db::exchange::messages).
        assert_eq!(
            messages::publickey("did:plc:a", "KEY==", None),
            "exchange-publickey\ndid:plc:a\nKEY==\n"
        );
        assert_eq!(
            messages::publickey("did:plc:a", "KEY==", Some("at://x")),
            "exchange-publickey\ndid:plc:a\nKEY==\nat://x"
        );
        assert_eq!(
            messages::request("urn:ibd:h", "did:plc:a", "did:plc:b", "IBD_AUTOSOMAL", None),
            "exchange-request\nurn:ibd:h\ndid:plc:a\ndid:plc:b\nIBD_AUTOSOMAL\n"
        );
        assert_eq!(
            messages::consent("urn:ibd:h", "did:plc:b", true),
            "exchange-consent\nurn:ibd:h\ndid:plc:b\ntrue"
        );
        assert_eq!(
            messages::poll("did:plc:a", 1718000000),
            "exchange-poll\ndid:plc:a\n1718000000"
        );
        assert_eq!(
            messages::relay("s1", "did:plc:a", "did:plc:b", 3, "HASH"),
            "exchange-relay\ns1\ndid:plc:a\ndid:plc:b\n3\nHASH"
        );
        assert_eq!(messages::ack("did:plc:a", 42), "exchange-ack\ndid:plc:a\n42");
    }
}
