//! Canonical signing strings for the signed Grid Edge API of the AppView (`/api/v1/grid/*`).
//!
//! These mirror `du_db::grid::messages` on the AppView, **exactly**. The server checks the
//! device-key signature against the string that it builds itself. Any difference here gives an
//! immediate 403, and the message tells the user nothing about the cause.
//!
//! Each string starts with its own operation name. So an attacker can not take a signature from one
//! endpoint and use it on a different endpoint.
//!
//! A change to a string here is a change to a published contract. A desktop version that signs the
//! old bytes stops work at the moment the server changes. The tests below hold each string, so an
//! accidental change fails here and not against a released version.
//!
//! Each call that changes data goes through
//! [`DeviceKey::sign_fresh`](crate::device_key::DeviceKey::sign_fresh), which puts the timestamp in
//! front as `{ts}\n{base}`. One signature then holds both the operation and the time. A read poll
//! signs the string here directly and puts its own `ts` in the string.

pub mod messages {
    /// `grid-register\n{did}\n{software_version}\n{caps_sha256_b64}`: announce a node and what it
    /// can do.
    ///
    /// The hash of the capabilities is in the signed string. So a node can not have capabilities
    /// that it did not send. That is important, because the claim path filters on them. A false
    /// claim of FASTQ ability would give the node work that it can not do.
    pub fn register(did: &str, software_version: &str, caps_sha256_b64: &str) -> String {
        format!("grid-register\n{did}\n{software_version}\n{caps_sha256_b64}")
    }

    /// `grid-poll\n{did}\n{ts}`: a read of what the caller has done, with a replay guard.
    pub fn poll(did: &str, ts: i64) -> String {
        format!("grid-poll\n{did}\n{ts}")
    }

    /// `grid-claim\n{did}\n{kinds}\n{count}\n{lease_secs}`: reserve up to `count` work units.
    ///
    /// `kinds` is the comma-joined list, in upper case and in alphabetical order. The client makes
    /// that form before it signs, and the server makes the same form before it checks. Without one
    /// agreed form, `["CRAM","cram"]` and `["cram","CRAM"]` give two different signed strings for
    /// one request.
    ///
    /// The signed values are the values that the node asks for. The server can reduce `count` or
    /// `lease_secs` to its own limits after the check. A node does not know those limits, and a
    /// signature over a value that the node can not calculate is not possible to make.
    pub fn claim(did: &str, kinds: &str, count: i32, lease_secs: i64) -> String {
        format!("grid-claim\n{did}\n{kinds}\n{count}\n{lease_secs}")
    }

    /// `grid-heartbeat\n{did}\n{lease_id}\n{stage}`: liveness for one lease, with the stage.
    pub fn heartbeat(did: &str, lease_id: i64, stage: &str) -> String {
        format!("grid-heartbeat\n{did}\n{lease_id}\n{stage}")
    }

    /// `grid-release\n{did}\n{lease_id}\n{reason}`: give a lease back with no result.
    ///
    /// The reason is in the signed string. So the server can not record a release reason that the
    /// node did not send.
    pub fn release(did: &str, lease_id: i64, reason: &str) -> String {
        format!("grid-release\n{did}\n{lease_id}\n{reason}")
    }

    /// `grid-submit\n{did}\n{work_unit_id}\n{digest_sha256_b64}`: send a result.
    ///
    /// The signature covers the **hash of the digest** and not the digest itself. The server
    /// calculates that hash again from the body that arrives. A signature over a hash proves only
    /// that the signer knew the hash. Without the second calculation, a node could sign the hash of
    /// a good result and send a different result.
    pub fn submit(did: &str, work_unit_id: i64, digest_sha256_b64: &str) -> String {
        format!("grid-submit\n{did}\n{work_unit_id}\n{digest_sha256_b64}")
    }

    /// Put the data kinds of a node into the one agreed form: upper case, no repeats, alphabetical
    /// order, joined with commas.
    ///
    /// The client and the server must make the same form. This function is the client half. The
    /// AppView handler does the same operation before it checks the signature.
    pub fn normalize_kinds(kinds: &[String]) -> String {
        let mut k: Vec<String> = kinds.iter().map(|s| s.trim().to_ascii_uppercase()).collect();
        k.sort();
        k.dedup();
        k.join(",")
    }
}

#[cfg(test)]
mod tests {
    use super::messages;

    /// The strings match the `du_db::grid::messages` literals of the AppView exactly. A change to
    /// one side only gives a 403 with no explanation, so this test is the guard.
    #[test]
    fn canonical_strings() {
        assert_eq!(
            messages::poll("did:plc:abc", 1_724_500_000),
            "grid-poll\ndid:plc:abc\n1724500000"
        );
        assert_eq!(
            messages::claim("did:plc:abc", "CRAM,FASTQ", 4, 259_200),
            "grid-claim\ndid:plc:abc\nCRAM,FASTQ\n4\n259200"
        );
        assert_eq!(
            messages::heartbeat("did:plc:abc", 7, "align"),
            "grid-heartbeat\ndid:plc:abc\n7\nalign"
        );
        assert_eq!(
            messages::release("did:plc:abc", 7, "cancelled"),
            "grid-release\ndid:plc:abc\n7\ncancelled"
        );
        assert_eq!(
            messages::submit("did:plc:abc", 12, "3q2+7w=="),
            "grid-submit\ndid:plc:abc\n12\n3q2+7w=="
        );
        assert_eq!(
            messages::register("did:plc:abc", "0.1.0-alpha.18", "3q2+7w=="),
            "grid-register\ndid:plc:abc\n0.1.0-alpha.18\n3q2+7w=="
        );
    }

    /// Each string starts with a different operation name. A signature from one endpoint is then of
    /// no use on a different endpoint.
    #[test]
    fn each_message_has_its_own_operation_name() {
        let all = [
            messages::poll("d", 1),
            messages::claim("d", "CRAM", 1, 1),
            messages::heartbeat("d", 1, "s"),
            messages::release("d", 1, "r"),
            messages::submit("d", 1, "h"),
            messages::register("d", "v", "h"),
        ];
        let names: Vec<&str> = all.iter().map(|m| m.split('\n').next().unwrap()).collect();
        let mut sorted = names.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            names.len(),
            "two grid messages share an operation name: {names:?}"
        );
    }

    /// The one agreed form of the data kinds. The AppView handler makes the same form, so these
    /// results must not change without a change there.
    #[test]
    fn data_kinds_have_one_agreed_form() {
        let k = |v: &[&str]| messages::normalize_kinds(&v.iter().map(|s| s.to_string()).collect::<Vec<_>>());
        assert_eq!(k(&["FASTQ", "CRAM"]), "CRAM,FASTQ");
        assert_eq!(k(&["cram", " CRAM ", "FASTQ"]), "CRAM,FASTQ");
        assert_eq!(k(&["CRAM"]), "CRAM");
        assert_eq!(k(&[]), "");
    }
}
