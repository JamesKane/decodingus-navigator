//! The PULL reconcile planner (gap §5-p2). This module has no I/O, so a unit test can cover every
//! case. The app does the plan that this module makes.
//!
//! The planner reads two inputs. The first input is the record of our last publish. That record is
//! the `sync_state` rows, and each row holds the PDS CID and the payload fingerprint at the time of
//! the push. The second input is the set of records that the PDS holds now. The planner then
//! decides the action for each record.
//!
//! The policy is **last-write-wins**, and the remote copy has authority when the two copies differ.
//! This was the §5-p2 decision. The app applies a record that changed on the PDS after our push. If
//! our local copy also changed after the push, the app still applies the remote record, but it
//! marks a conflict. The payload hash shows a local change.

use navigator_store::sync_state::StoredSyncState;
use navigator_sync::RemoteRecord;

/// One reconcile decision for a record.
#[derive(Debug, Clone, PartialEq)]
pub enum ReconcileAction {
    /// The remote record is the same as our published record, and the local record did not change.
    /// There is no action.
    InSync { entity_ref: String },
    /// The remote record changed after our push, or there is a local edit. Apply the remote record
    /// to the local record. `conflict` shows that the local record also changed after the push.
    /// Both copies changed, the remote copy wins, and the app writes a log entry.
    ApplyRemote {
        entity_ref: String,
        collection: String,
        remote: RemoteRecord,
        conflict: bool,
    },
    /// The app published the local record, but the PDS no longer holds it. Publish our copy
    /// again.
    RePush { entity_ref: String },
    /// The PDS holds a record, and we have no local sync-state row for it. Add the record to the
    /// local store.
    AdoptRemote { collection: String, remote: RemoteRecord },
}

/// Make the reconcile plan for one collection.
///
/// `local` gives the current local payload hash of each published entity. A value of `None` means
/// that the app did not calculate the hash again, and the planner treats the record as clean. The
/// planner compares this hash with the hash from the time of the push. A difference shows a local
/// change. `remote` is the set of records that the PDS holds now for the same collection.
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
                // The app published the record, but the PDS no longer holds it. Publish again.
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
