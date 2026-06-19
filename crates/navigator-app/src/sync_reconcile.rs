//! Pure PULL reconcile planner (gap §5-p2). Given what we last published (`sync_state` rows, each with
//! the PDS CID + payload fingerprint at push time) and the records currently on the PDS, decide what to
//! do per record — with **no I/O**, so it's exhaustively unit-tested. The app executes the plan.
//!
//! Policy: **last-write-wins, remote authoritative on divergence** (the confirmed §5-p2 decision). A
//! record that changed on the PDS since our push is applied locally; if our local copy *also* changed
//! since the push (we can detect that via the payload hash), it's still applied but flagged a conflict.

use navigator_store::sync_state::StoredSyncState;
use navigator_sync::RemoteRecord;

/// One reconcile decision for a record.
#[derive(Debug, Clone, PartialEq)]
pub enum ReconcileAction {
    /// Remote matches what we published and local is unchanged — nothing to do.
    InSync { entity_ref: String },
    /// Remote changed since our push (or we have a local-only edit) — apply remote→local. `conflict`
    /// ⇒ local *also* changed since the push (both diverged; remote wins, logged).
    ApplyRemote {
        entity_ref: String,
        collection: String,
        remote: RemoteRecord,
        conflict: bool,
    },
    /// Local was published but the record is gone on the PDS — re-publish our copy.
    RePush { entity_ref: String },
    /// A record exists on the PDS we have no local sync-state for — adopt it locally.
    AdoptRemote { collection: String, remote: RemoteRecord },
}

/// Plan the reconcile for one collection. `local` pairs each published entity with its *current* local
/// payload hash (`None` = not recomputed / assume clean); compared to the stored push-time hash to tell
/// whether local changed. `remote` is the PDS's current records for the same collection.
pub fn plan(local: &[(StoredSyncState, Option<String>)], remote: &[RemoteRecord]) -> Vec<ReconcileAction> {
    use std::collections::HashSet;
    let mut actions = Vec::new();
    let mut matched_rkeys: HashSet<&str> = HashSet::new();

    for (ss, local_hash) in local {
        let local_dirty = local_hash.as_ref().is_some_and(|h| h != &ss.payload_hash);
        match remote.iter().find(|r| r.rkey() == ss.rkey) {
            Some(r) => {
                matched_rkeys.insert(r.rkey());
                if r.cid == ss.at_cid {
                    // Remote unchanged since our push.
                    if local_dirty {
                        actions.push(ReconcileAction::RePush {
                            entity_ref: ss.entity_ref.clone(),
                        });
                    } else {
                        actions.push(ReconcileAction::InSync {
                            entity_ref: ss.entity_ref.clone(),
                        });
                    }
                } else {
                    // Remote moved on (edited elsewhere) → apply; conflict if we also have local edits.
                    actions.push(ReconcileAction::ApplyRemote {
                        entity_ref: ss.entity_ref.clone(),
                        collection: ss.collection.clone(),
                        remote: r.clone(),
                        conflict: local_dirty,
                    });
                }
            }
            None => {
                // We published it but it's gone on the PDS — re-publish.
                actions.push(ReconcileAction::RePush {
                    entity_ref: ss.entity_ref.clone(),
                });
            }
        }
    }

    for r in remote {
        if !matched_rkeys.contains(r.rkey()) {
            // Pick the collection off the at-uri (…/<collection>/<rkey>).
            let collection = r.uri.rsplit('/').nth(1).unwrap_or("").to_string();
            actions.push(ReconcileAction::AdoptRemote {
                collection,
                remote: r.clone(),
            });
        }
    }
    actions
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ss(entity: &str, rkey: &str, cid: &str, hash: &str) -> StoredSyncState {
        StoredSyncState {
            account_did: "did:plc:me".into(),
            entity_ref: entity.into(),
            kind: "biosample".into(),
            collection: "com.decodingus.biosample".into(),
            rkey: rkey.into(),
            at_uri: format!("at://did:plc:me/com.decodingus.biosample/{rkey}"),
            at_cid: cid.into(),
            payload_hash: hash.into(),
            pushed_at: "t".into(),
        }
    }
    fn rr(rkey: &str, cid: &str) -> RemoteRecord {
        RemoteRecord {
            uri: format!("at://did:plc:me/com.decodingus.biosample/{rkey}"),
            cid: cid.into(),
            value: serde_json::json!({}),
        }
    }

    #[test]
    fn in_sync_when_cid_matches_and_local_clean() {
        let local = vec![(ss("biosample:1", "rk1", "cidA", "h1"), Some("h1".to_string()))];
        let remote = vec![rr("rk1", "cidA")];
        assert_eq!(
            plan(&local, &remote),
            vec![ReconcileAction::InSync {
                entity_ref: "biosample:1".into()
            }]
        );
    }

    #[test]
    fn apply_remote_when_remote_changed() {
        let local = vec![(ss("biosample:1", "rk1", "cidA", "h1"), Some("h1".to_string()))];
        let remote = vec![rr("rk1", "cidB")];
        let a = plan(&local, &remote);
        assert!(matches!(&a[0], ReconcileAction::ApplyRemote { conflict: false, .. }));
    }

    #[test]
    fn conflict_when_both_changed() {
        // remote cid differs AND local hash differs from the push-time hash.
        let local = vec![(
            ss("biosample:1", "rk1", "cidA", "h1"),
            Some("h2-local-edit".to_string()),
        )];
        let remote = vec![rr("rk1", "cidB")];
        let a = plan(&local, &remote);
        assert!(
            matches!(&a[0], ReconcileAction::ApplyRemote { conflict: true, .. }),
            "both diverged → conflict (remote wins)"
        );
    }

    #[test]
    fn repush_when_local_dirty_but_remote_unchanged() {
        let local = vec![(ss("biosample:1", "rk1", "cidA", "h1"), Some("h2".to_string()))];
        let remote = vec![rr("rk1", "cidA")];
        assert_eq!(
            plan(&local, &remote),
            vec![ReconcileAction::RePush {
                entity_ref: "biosample:1".into()
            }]
        );
    }

    #[test]
    fn repush_when_remote_missing() {
        let local = vec![(ss("biosample:1", "rk1", "cidA", "h1"), None)];
        assert_eq!(
            plan(&local, &[]),
            vec![ReconcileAction::RePush {
                entity_ref: "biosample:1".into()
            }]
        );
    }

    #[test]
    fn adopt_remote_when_unknown_locally() {
        let a = plan(&[], &[rr("rk9", "cidZ")]);
        assert!(
            matches!(&a[0], ReconcileAction::AdoptRemote { collection, remote } if collection == "com.decodingus.biosample" && remote.rkey() == "rk9")
        );
    }
}
