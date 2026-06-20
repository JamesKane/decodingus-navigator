//! Canonical signing strings for the AppView's signed social Edge API (`/api/v1/social/*`).
//!
//! These mirror `du_db::social::messages` (+ `du_db::notification::messages::read`) on the AppView
//! **byte-for-byte** — the device-key signature is verified against the exact same string the server
//! reconstructs, so any drift here is an instant 403. Each is `\n`-joined; an absent optional field
//! is the empty string (never omitted). Signed with the Ed25519 device key
//! ([`crate::device_key::DeviceKey::sign`]); base64-standard signature.

pub mod messages {
    /// `social-thread\n{did}\n{conversation_id_or_empty}` — open a new support thread
    /// (`conversation_id` absent) or reply to one the caller owns.
    pub fn thread(did: &str, conversation_id: Option<&str>) -> String {
        format!("social-thread\n{did}\n{}", conversation_id.unwrap_or(""))
    }

    /// `social-poll\n{did}\n{ts}` — replay-guarded read poll (list threads / read feed / read
    /// notifications); proves the caller is `did` at unix-seconds `ts`.
    pub fn poll(did: &str, ts: i64) -> String {
        format!("social-poll\n{did}\n{ts}")
    }

    /// `social-thread-read\n{did}\n{conversation_id}\n{ts}` — replay-guarded read of one thread's
    /// messages (also marks it read on the user side).
    pub fn thread_read(did: &str, conversation_id: &str, ts: i64) -> String {
        format!("social-thread-read\n{did}\n{conversation_id}\n{ts}")
    }

    /// `social-post\n{did}\n{parent_or_empty}` — create a community feed post (`parent` absent) or a
    /// reply (`parent` set).
    pub fn post(did: &str, parent: Option<&str>) -> String {
        format!("social-post\n{did}\n{}", parent.unwrap_or(""))
    }

    /// `social-notif-read\n{did}\n{id_or_empty}\n{ts}` — mark one notification read (`id` set) or all
    /// (`id` absent); replay-guarded by `ts`.
    pub fn notif_read(did: &str, id: Option<&str>, ts: i64) -> String {
        format!("social-notif-read\n{did}\n{}\n{ts}", id.unwrap_or(""))
    }
}

#[cfg(test)]
mod tests {
    use super::messages;

    #[test]
    fn canonical_strings_match_the_appview() {
        assert_eq!(messages::thread("did:k", None), "social-thread\ndid:k\n");
        assert_eq!(messages::thread("did:k", Some("c1")), "social-thread\ndid:k\nc1");
        assert_eq!(messages::poll("did:k", 42), "social-poll\ndid:k\n42");
        assert_eq!(messages::thread_read("did:k", "c1", 42), "social-thread-read\ndid:k\nc1\n42");
        assert_eq!(messages::post("did:k", None), "social-post\ndid:k\n");
        assert_eq!(messages::post("did:k", Some("p1")), "social-post\ndid:k\np1");
        assert_eq!(messages::notif_read("did:k", None, 42), "social-notif-read\ndid:k\n\n42");
        assert_eq!(messages::notif_read("did:k", Some("n1"), 42), "social-notif-read\ndid:k\nn1\n42");
    }
}
