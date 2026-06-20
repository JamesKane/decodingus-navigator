//! Canonical signing strings for the AppView's signed recruitment Edge API
//! (`/api/v1/recruitment/*`, social roadmap 3c).
//!
//! These mirror `du_db::recruitment::messages` on the AppView **byte-for-byte** — the device-key
//! signature is verified against the exact string the server reconstructs, so any drift here is an
//! instant 403. Respond-only: poll the caller's open invitations + accept/decline one. Signed with
//! the Ed25519 device key ([`crate::device_key::DeviceKey::sign`]); base64-standard signature.

pub mod messages {
    /// `recruitment-poll\n{did}\n{ts}` — replay-guarded poll for the caller's open invitations.
    pub fn poll(did: &str, ts: i64) -> String {
        format!("recruitment-poll\n{did}\n{ts}")
    }
    /// `recruitment-respond\n{did}\n{campaign_id}\n{accept}` — accept (`true`) / decline (`false`).
    pub fn respond(did: &str, campaign_id: i64, accept: bool) -> String {
        format!("recruitment-respond\n{did}\n{campaign_id}\n{accept}")
    }
}

#[cfg(test)]
mod tests {
    use super::messages;

    /// The strings match the AppView's `du_db::recruitment::messages` literals exactly.
    #[test]
    fn canonical_strings() {
        assert_eq!(messages::poll("did:plc:x", 1700), "recruitment-poll\ndid:plc:x\n1700");
        assert_eq!(
            messages::respond("did:plc:x", 42, true),
            "recruitment-respond\ndid:plc:x\n42\ntrue"
        );
        assert_eq!(
            messages::respond("did:plc:x", 42, false),
            "recruitment-respond\ndid:plc:x\n42\nfalse"
        );
    }
}
