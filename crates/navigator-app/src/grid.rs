//! The Grid client: `impl App` methods for the signed Grid Edge API of the AppView
//! (`/api/v1/grid/*`).
//!
//! A node announces itself, reserves work units, reports progress, gives a lease back, and sends a
//! signed result. The device key signs each call, as it does for the exchange client and the
//! recruitment client. This module uses the shared [`appview_post`](App::appview_post) and
//! [`appview_get_signed`](App::appview_get_signed) transport.
//!
//! The canonical strings are in [`navigator_sync::grid::messages`], which mirrors
//! `du_db::grid::messages` on the AppView. Design
//! `documents/design/distributed-compute-grid.md` §4.4 and §7.1.
//!
//! # What a signature covers, and why the digest is a hash
//!
//! A call that changes data signs `{ts}\n{base}` through
//! [`DeviceKey::sign_fresh`](navigator_sync::device_key::DeviceKey::sign_fresh). One signature then
//! holds the operation and the time. The AppView keeps each accepted signature for a short period.
//! It then refuses the same bytes a second time.
//!
//! The submit call signs the **hash** of the digest and sends the digest with it. The AppView
//! calculates the hash again from the body that arrives, and refuses a difference.
//!
//! The canonical bytes are what `serde_json` writes. This workspace does not use the
//! `preserve_order` feature. So the keys are in alphabetical order, and the output has no spaces.
//! The AppView uses the same crate with the same setting. That gives the smallest possible
//! agreement between the two repositories. There is no field order to agree, and no number format
//! rules.

use super::*;
use crate::ena::ManifestFile;
use navigator_sync::grid::messages;

/// A work unit as the AppView gives it at claim time. It holds everything that a node needs, so the
/// node asks ENA for nothing.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct ClaimedUnit {
    pub lease_id: i64,
    pub work_unit_id: i64,
    pub sample_accession: String,
    #[serde(default)]
    pub study_accession: Option<String>,
    /// `CRAM` for a unit that needs no new alignment, or `FASTQ` for a unit that the node maps.
    pub data_kind: String,
    #[serde(default)]
    pub manifest: Vec<ManifestFile>,
    #[serde(default)]
    pub est_bases: Option<i64>,
    #[serde(default)]
    pub total_bytes: Option<i64>,
    pub expires_at: chrono::DateTime<chrono::Utc>,
}

/// What a contributor has done, and where the contributor is on the public board.
#[derive(Debug, Clone, Default, serde::Deserialize)]
pub struct GridStanding {
    #[serde(default)]
    pub leases: Vec<ClaimedUnit>,
    #[serde(default)]
    pub agreed: i64,
    #[serde(default)]
    pub divergent: i64,
    #[serde(default)]
    pub cobblestones: f64,
    #[serde(default)]
    pub units_credited: i64,
    /// The position on the board. It is `None` for a contributor with no credit. Such a
    /// contributor has no row on the board, so a number here would be an answer to a question that
    /// nobody asked.
    #[serde(default)]
    pub rank: Option<i64>,
}

/// What this node can do. The AppView keeps it and uses it to select work.
#[derive(Debug, Clone, serde::Serialize)]
pub struct NodeCapabilities {
    /// `["CRAM"]`, or `["CRAM","FASTQ"]` for a node that can map reads.
    pub data_kinds: Vec<String>,
    pub threads: u32,
    /// The disk space, in bytes, that the user gives to this work.
    pub disk_budget: u64,
    pub memory_bytes: u64,
}

impl App {
    /// Announce this node, or send its capabilities again.
    ///
    /// This call is also the heartbeat of the node. The AppView row keeps one `last_heartbeat`
    /// value, and this call sets it. A second endpoint that writes the same row would let the two
    /// values disagree.
    pub async fn grid_register(&self, caps: &NodeCapabilities) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = chrono::Utc::now().timestamp();
        let caps_value = serde_json::to_value(caps).map_err(|e| AppError::Import(e.to_string()))?;
        let caps_hash = canonical_sha256_b64(&caps_value);
        let version = env!("CARGO_PKG_VERSION");
        let sig = dev.sign_fresh(ts, &messages::register(&did, version, &caps_hash));
        let body = serde_json::json!({
            "did": did,
            "software_version": version,
            "capabilities": caps_value,
            "os_info": os_info(),
            "ts": ts,
            "signature": sig,
        });
        let v = self.appview_post("grid/node/register", body).await?;
        Ok(v.get("node_id").and_then(|x| x.as_i64()).unwrap_or_default())
    }

    /// Reserve up to `count` work units for `lease_secs` seconds.
    ///
    /// The node sends only the kinds that it can process. The AppView never gives a FASTQ unit to a
    /// node that can not map reads.
    ///
    /// The result can hold fewer units than `count`, or none. An empty result is the usual answer
    /// when the catalogue holds no more work of those kinds. It is not an error.
    pub async fn grid_claim(
        &self,
        data_kinds: &[String],
        count: i32,
        lease_secs: i64,
    ) -> Result<Vec<ClaimedUnit>, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = chrono::Utc::now().timestamp();
        let kinds = messages::normalize_kinds(data_kinds);
        if kinds.is_empty() {
            return Err(AppError::Import("this node advertises no data kinds".into()));
        }
        let sig = dev.sign_fresh(ts, &messages::claim(&did, &kinds, count, lease_secs));
        let body = serde_json::json!({
            "did": did,
            // The list goes on the wire in the same form that the signature covers.
            "data_kinds": kinds.split(',').collect::<Vec<_>>(),
            "count": count,
            "lease_secs": lease_secs,
            "ts": ts,
            "signature": sig,
        });
        let v = self.appview_post("grid/claim", body).await?;
        let units = v.get("units").cloned().unwrap_or_else(|| serde_json::json!([]));
        serde_json::from_value(units).map_err(|e| AppError::Import(e.to_string()))
    }

    /// Report progress on a lease that this node holds.
    ///
    /// Returns `false` when the lease is no longer the lease of this node. The node must then stop
    /// work on that unit. Without this answer, a node can spend hours on a unit that it lost, and
    /// it receives no credit for that work.
    ///
    /// This call does **not** make the lease longer. A node can send a heartbeat and still not
    /// finish. Such a node would hold a unit for ever. A limited lease prevents that fault.
    pub async fn grid_heartbeat(&self, lease_id: i64, stage: &str, fraction: Option<f32>) -> Result<bool, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = chrono::Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &messages::heartbeat(&did, lease_id, stage));
        let body = serde_json::json!({
            "did": did,
            "lease_id": lease_id,
            "stage": stage,
            "progress": { "stage": stage, "fraction": fraction },
            "ts": ts,
            "signature": sig,
        });
        let v = self.appview_post("grid/heartbeat", body).await?;
        // Only an explicit `false` means that this node lost the lease.
        //
        // An absent field, a new name for it, or a value of another type gives `None` here. An
        // earlier version read each of those as `false`, and the node then stopped work of many
        // hours.
        //
        // The lease has its own time limit. So the safe answer to an unclear reply is to continue.
        // At worst, the node finishes a unit that another node also finished.
        Ok(v.get("held").and_then(|x| x.as_bool()).unwrap_or(true))
    }

    /// Give a lease back with no result, so another node can take the unit immediately.
    ///
    /// A second call for the same lease is safe. It returns `false`, and that is not an error: a
    /// node that sends the call again after a lost answer did nothing wrong.
    pub async fn grid_release(&self, lease_id: i64, reason: &str) -> Result<bool, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = chrono::Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &messages::release(&did, lease_id, reason));
        let body = serde_json::json!({
            "did": did,
            "lease_id": lease_id,
            "reason": reason,
            "ts": ts,
            "signature": sig,
        });
        let v = self.appview_post("grid/release", body).await?;
        Ok(v.get("released").and_then(|x| x.as_bool()).unwrap_or(false))
    }

    /// Send the result of a unit, and close the lease that made it.
    ///
    /// `digest` holds **raw** values. The AppView puts the continuous values into groups when it
    /// compares two results. A client that made the groups itself would put the group rule into two
    /// repositories. A difference between them would then give `DIVERGENT` results against nodes
    /// that did nothing wrong, and the message would give no cause.
    ///
    /// A second call for the same unit replaces the first result. It does not add a second vote.
    #[allow(clippy::too_many_arguments)]
    pub async fn grid_submit(
        &self,
        work_unit_id: i64,
        lease_id: Option<i64>,
        digest: &serde_json::Value,
        stack_version: &str,
        reference_build: &str,
        aligner: Option<&str>,
        record_refs: &[String],
    ) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = chrono::Utc::now().timestamp();
        let hash = canonical_sha256_b64(digest);
        // Two signatures, for two different purposes. The request signature proves who sent this
        // call now. The digest signature stays in the row. A later check can then prove which node
        // made this result, after the AppView copies it to other places.
        let sig = dev.sign_fresh(ts, &messages::submit(&did, work_unit_id, &hash));
        let digest_sig = dev.sign(&canonical_bytes_string(digest));
        let body = serde_json::json!({
            "did": did,
            "work_unit_id": work_unit_id,
            "lease_id": lease_id,
            "digest": digest,
            "digest_sig": digest_sig,
            "stack_version": stack_version,
            "reference_build": reference_build,
            "aligner": aligner,
            "record_refs": record_refs,
            "ts": ts,
            "signature": sig,
        });
        let v = self.appview_post("grid/submit", body).await?;
        Ok(v.get("submission_id").and_then(|x| x.as_i64()).unwrap_or_default())
    }

    /// The leases, the history and the board position of this node.
    pub async fn grid_standing(&self) -> Result<GridStanding, AppError> {
        self.appview_get_signed("grid/mine", messages::poll, &[]).await
    }
}

/// What the digest signature covers: the canonical bytes of the digest, as text.
///
/// `serde_json` writes the keys in alphabetical order and adds no spaces, because this workspace
/// does not use the `preserve_order` feature. The AppView uses the same crate with the same
/// setting, so both sides make the same bytes with no rules to agree.
fn canonical_bytes_string(value: &serde_json::Value) -> String {
    serde_json::to_string(value).unwrap_or_default()
}

/// The SHA-256 of the canonical bytes of a JSON value, as standard base64.
///
/// This must give the same answer as `du_db::grid::digest::canonical_sha256_b64` on the AppView.
/// The submit handler calculates it again from the body that arrives, and refuses a difference.
pub(crate) fn canonical_sha256_b64(value: &serde_json::Value) -> String {
    use base64::Engine as _;
    use sha2::{Digest as _, Sha256};
    let bytes = serde_json::to_vec(value).unwrap_or_default();
    base64::engine::general_purpose::STANDARD.encode(Sha256::digest(bytes))
}

/// A short description of this machine, for the fleet view.
fn os_info() -> String {
    format!("{} {}", std::env::consts::OS, std::env::consts::ARCH)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The hash must not change with the order of the keys in the source text. If it did, a node
    /// and the AppView could disagree about what the signature covers.
    #[test]
    fn the_canonical_hash_ignores_key_order() {
        let a = serde_json::json!({"calls": {"sex": "XY", "y_terminal": "R-A"}, "unit": "SAMEA1"});
        let b = serde_json::json!({"unit": "SAMEA1", "calls": {"y_terminal": "R-A", "sex": "XY"}});
        assert_eq!(canonical_sha256_b64(&a), canonical_sha256_b64(&b));
    }

    /// A different result must give a different hash, or the check has no value.
    #[test]
    fn a_different_result_gives_a_different_hash() {
        let a = serde_json::json!({"calls": {"sex": "XY"}});
        let b = serde_json::json!({"calls": {"sex": "XX"}});
        assert_ne!(canonical_sha256_b64(&a), canonical_sha256_b64(&b));
    }

    /// The bytes that the digest signature covers are the bytes that go on the wire.
    #[test]
    fn the_signed_bytes_are_the_bytes_that_are_sent() {
        let v = serde_json::json!({"b": 2, "a": 1});
        assert_eq!(canonical_bytes_string(&v), r#"{"a":1,"b":2}"#);
    }

    /// A unit as the AppView sends it, with the manifest that `grid-curate` made.
    #[test]
    fn a_claimed_unit_decodes_with_its_manifest() {
        let json = serde_json::json!({
            "lease_id": 7,
            "work_unit_id": 12,
            "sample_accession": "SAMEA0000001",
            "study_accession": "PRJEB00000",
            "data_kind": "CRAM",
            "manifest": [{
                "run_accession": "ERR0000001",
                "url": "ftp.sra.ebi.ac.uk/vol1/s1.cram",
                "md5": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                "bytes": 12000000000i64,
                "format": "CRAM"
            }],
            "est_bases": 90000000000i64,
            "total_bytes": 12000000000i64,
            "expires_at": "2026-08-28T00:00:00Z"
        });
        let u: ClaimedUnit = serde_json::from_value(json).expect("decode");
        assert_eq!(u.manifest.len(), 1);
        assert_eq!(u.manifest[0].file_name(), "s1.cram");
        assert_eq!(u.data_kind, "CRAM");
    }

    /// A contributor with no credit has no board position. The field must arrive as `None` and not
    /// as a number.
    #[test]
    fn an_uncredited_contributor_has_no_rank() {
        let json = serde_json::json!({
            "leases": [], "agreed": 0, "divergent": 0,
            "cobblestones": 0.0, "units_credited": 0, "rank": null
        });
        let s: GridStanding = serde_json::from_value(json).expect("decode");
        assert!(s.rank.is_none());
        assert_eq!(s.agreed, 0);
    }
}
