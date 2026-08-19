//! Canonical signing strings for the AppView's signed recruitment Edge API
//! (`/api/v1/recruitment/*`, social roadmap 3c).
//!
//! These mirror `du_db::recruitment::messages` on the AppView, **exactly**. The server
//! checks the device-key signature against the exact string that it builds itself, so any
//! difference here gives an instant 403.
//!
//! The API is respond-only: poll the caller's open invitations, and accept or decline one. The
//! Ed25519 device key signs it ([`crate::device_key::DeviceKey::sign`]), and the signature is
//! base64-standard.

pub mod messages {
    /// `recruitment-poll\n{did}\n{ts}`: a poll for the caller's open invitations, with a replay
    /// guard.
    pub fn poll(did: &str, ts: i64) -> String {
        format!("recruitment-poll\n{did}\n{ts}")
    }
    /// `recruitment-respond\n{did}\n{campaign_id}\n{accept}`: accept (`true`) or decline
    /// (`false`).
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
