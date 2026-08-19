//! Canonical signing strings for the AppView's signed social Edge API (`/api/v1/social/*`).
//!
//! These mirror `du_db::social::messages`, and `du_db::notification::messages::read`, on the
//! AppView, **exactly**. The server checks the device-key signature against the exact string
//! that it builds itself, so any difference here gives an instant 403.
//!
//! A `\n` joins the parts of each string. An optional field that has no value becomes the empty
//! string, and the code never leaves it out. The Ed25519 device key signs it
//! ([`crate::device_key::DeviceKey::sign`]), and the signature is base64-standard.

pub mod messages {
    /// `social-thread\n{did}\n{conversation_id_or_empty}`: open a new support thread, with no
    /// `conversation_id`, or reply to one that the caller owns.
    pub fn thread(did: &str, conversation_id: Option<&str>) -> String {
        format!("social-thread\n{did}\n{}", conversation_id.unwrap_or(""))
    }

    /// `social-poll\n{did}\n{ts}`: a read poll with a replay guard, which lists threads, reads the
    /// feed, or reads notifications. It proves that the caller is `did` at unix-seconds `ts`.
    pub fn poll(did: &str, ts: i64) -> String {
        format!("social-poll\n{did}\n{ts}")
    }

    /// `social-thread-read\n{did}\n{conversation_id}\n{ts}`: a read of one thread's messages, with
    /// a replay guard. It also marks that thread read on the user side.
    pub fn thread_read(did: &str, conversation_id: &str, ts: i64) -> String {
        format!("social-thread-read\n{did}\n{conversation_id}\n{ts}")
    }

    /// `social-post\n{did}\n{parent_or_empty}`: create a community feed post, with no `parent`, or
    /// a reply, with `parent` set.
    pub fn post(did: &str, parent: Option<&str>) -> String {
        format!("social-post\n{did}\n{}", parent.unwrap_or(""))
    }

    /// `social-notif-read\n{did}\n{id_or_empty}\n{ts}`: mark one notification read, with `id` set,
    /// or mark all of them read, with no `id`. `ts` is the replay guard.
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
        assert_eq!(
            messages::thread_read("did:k", "c1", 42),
            "social-thread-read\ndid:k\nc1\n42"
        );
        assert_eq!(messages::post("did:k", None), "social-post\ndid:k\n");
        assert_eq!(messages::post("did:k", Some("p1")), "social-post\ndid:k\np1");
        assert_eq!(
            messages::notif_read("did:k", None, 42),
            "social-notif-read\ndid:k\n\n42"
        );
        assert_eq!(
            messages::notif_read("did:k", Some("n1"), 42),
            "social-notif-read\ndid:k\nn1\n42"
        );
    }
}
