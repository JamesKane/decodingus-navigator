//! IBD segment-exchange payload + match attestation (gap §4 application layer).
//!
//! IBD detection needs **both** peers' genotypes, so over the encrypted edge channel each peer sends
//! its dosages at the canonical IBD-panel sites ([`IbdSite`], an [`IbdExchangeMsg::Dosages`]); each
//! then runs the symmetric [`crate::ibd::PairwiseIbdDetector`] locally → an identical
//! [`crate::ibd::MatchSummary`]. Each peer signs an [`IbdAttestation`] over its computed summary and
//! exchanges it; agreement = both signed attestations carry the same `summary_hash` (proof both
//! computed the same result). Only panel dosages cross the wire — encrypted, never seen by the broker.
//!
//! This module is pure: it defines the wire types, the canonical signing string, and the summary
//! hash. Signing (the device key) and verification (`du_atproto::verify_did_key`) happen in the app.

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::ibd::MatchSummary;

/// One panel-site dosage on the wire — the minimal input the IBD detector consumes (the heavy
/// [`crate::caller::SiteGenotype`] fields aren't sent). `dosage` is 0/1/2, or -1 for no-call.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct IbdSite {
    pub contig: String,
    pub position: i64,
    pub dosage: i32,
}

/// A signed claim that the attester computed a given IBD match summary for a session/partner. Both
/// peers publish one; equal `summary_hash` across the pair proves they computed identical results.
/// `signature` is STANDARD-base64 Ed25519 over [`IbdAttestation::canonical`], by the device key whose
/// `did:key` is `signing_public_key` (verified with `du_atproto::verify_did_key`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IbdAttestation {
    /// The exchange request URI this match belongs to.
    pub match_request_uri: String,
    /// The opened session id.
    pub session_id: String,
    /// The attester's account DID.
    pub attesting_did: String,
    /// Opaque references to each side's biosample (never an identifier — a guid/at-uri).
    pub attesting_sample_ref: Option<String>,
    pub partner_sample_ref: Option<String>,
    pub total_shared_cm: f64,
    pub segment_count: u32,
    pub longest_segment_cm: f64,
    /// The relationship band (`RelationshipEstimate` debug name).
    pub relationship: String,
    /// SHA-256 (base64) of the attester's match summary — the agreement fingerprint.
    pub summary_hash: String,
    /// RFC3339 timestamp.
    pub attested_at: String,
    /// STANDARD-base64 Ed25519 signature over [`canonical`](Self::canonical) (empty until signed).
    pub signature: String,
    /// The signer's `did:key` (the device key).
    pub signing_public_key: String,
}

impl IbdAttestation {
    /// Build an unsigned attestation from a computed summary. The app fills `signature` +
    /// `signing_public_key` after signing [`canonical`](Self::canonical).
    #[allow(clippy::too_many_arguments)]
    pub fn unsigned(
        match_request_uri: impl Into<String>,
        session_id: impl Into<String>,
        attesting_did: impl Into<String>,
        attesting_sample_ref: Option<String>,
        partner_sample_ref: Option<String>,
        summary: &MatchSummary,
        attested_at: impl Into<String>,
    ) -> Self {
        IbdAttestation {
            match_request_uri: match_request_uri.into(),
            session_id: session_id.into(),
            attesting_did: attesting_did.into(),
            attesting_sample_ref,
            partner_sample_ref,
            total_shared_cm: summary.total_shared_cm,
            segment_count: summary.segment_count as u32,
            longest_segment_cm: summary.longest_segment_cm,
            relationship: format!("{:?}", summary.relationship),
            summary_hash: summary_hash(summary),
            attested_at: attested_at.into(),
            signature: String::new(),
            signing_public_key: String::new(),
        }
    }

    /// The canonical `\n`-joined message that gets signed/verified — every field except the signature
    /// and the signing key (which carries the signature). Both peers build it byte-identically.
    pub fn canonical(&self) -> String {
        [
            self.match_request_uri.as_str(),
            self.session_id.as_str(),
            self.attesting_did.as_str(),
            self.attesting_sample_ref.as_deref().unwrap_or(""),
            self.partner_sample_ref.as_deref().unwrap_or(""),
            &format!("{:.4}", self.total_shared_cm),
            &self.segment_count.to_string(),
            &format!("{:.4}", self.longest_segment_cm),
            self.relationship.as_str(),
            self.summary_hash.as_str(),
            self.attested_at.as_str(),
        ]
        .join("\n")
    }
}

/// SHA-256 (lowercase hex) of a match summary's salient fields — the cross-peer agreement
/// fingerprint. Deterministic so both peers, computing the same segments, get the same hash.
pub fn summary_hash(s: &MatchSummary) -> String {
    let canon = format!(
        "{:.4}|{}|{:.4}|{:?}",
        s.total_shared_cm, s.segment_count, s.longest_segment_cm, s.relationship
    );
    Sha256::digest(canon.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// A message exchanged inside the encrypted channel during an IBD exchange.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum IbdExchangeMsg {
    /// The sender's panel dosages (`build` = the panel's reference build, e.g. `hs1`).
    Dosages { build: String, sites: Vec<IbdSite> },
    /// The sender's signed match attestation (boxed — much larger than the other variant).
    Attest(Box<IbdAttestation>),
}

impl IbdExchangeMsg {
    /// Gzipped JSON bytes for the channel. The dosage payload is large (a panel of sites), and the
    /// relay caps an envelope at 1 MiB, so it's compressed (the dosage vector compresses well).
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        use flate2::{write::GzEncoder, Compression};
        use std::io::Write;
        let json = serde_json::to_vec(self).map_err(|e| e.to_string())?;
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&json).map_err(|e| e.to_string())?;
        enc.finish().map_err(|e| e.to_string())
    }
    pub fn from_bytes(b: &[u8]) -> Result<Self, String> {
        use flate2::read::GzDecoder;
        use std::io::Read;
        let mut dec = GzDecoder::new(b);
        let mut json = Vec::new();
        dec.read_to_end(&mut json).map_err(|e| e.to_string())?;
        serde_json::from_slice(&json).map_err(|e| e.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ibd::{IbdSegment, MatchSummary};

    fn summary(cm: f64) -> MatchSummary {
        MatchSummary::from_segments(&[IbdSegment {
            chromosome: "chr1".into(),
            start_position: 1,
            end_position: 10_000_000,
            length_cm: cm,
            snp_count: Some(500),
            is_half_identical: None,
        }])
    }

    #[test]
    fn summary_hash_is_deterministic_and_distinguishing() {
        assert_eq!(summary_hash(&summary(20.0)), summary_hash(&summary(20.0)));
        assert_ne!(summary_hash(&summary(20.0)), summary_hash(&summary(40.0)));
    }

    #[test]
    fn canonical_is_stable_and_excludes_signature() {
        let mut a = IbdAttestation::unsigned(
            "exchange:r1",
            "s1",
            "did:key:zA",
            None,
            None,
            &summary(20.0),
            "2026-06-17T00:00:00Z",
        );
        let c1 = a.canonical();
        a.signature = "sig".into();
        a.signing_public_key = "did:key:zA".into();
        assert_eq!(a.canonical(), c1, "signing must not change the canonical message");
        assert!(c1.contains("exchange:r1") && c1.contains(&a.summary_hash));
    }

    #[test]
    fn exchange_msg_round_trips() {
        let m = IbdExchangeMsg::Dosages {
            build: "hs1".into(),
            sites: vec![IbdSite {
                contig: "chr1".into(),
                position: 100,
                dosage: 2,
            }],
        };
        assert_eq!(IbdExchangeMsg::from_bytes(&m.to_bytes().unwrap()).unwrap(), m);
        let att = IbdExchangeMsg::Attest(Box::new(IbdAttestation::unsigned(
            "r",
            "s",
            "d",
            None,
            None,
            &summary(7.0),
            "t",
        )));
        assert_eq!(IbdExchangeMsg::from_bytes(&att.to_bytes().unwrap()).unwrap(), att);
    }
}
