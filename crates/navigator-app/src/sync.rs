//! `impl App` methods extracted from `lib.rs` (the `sync` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- sync durability: outbox enqueue + drain (gap §5) -------------------

    /// Put a record in the queue for a publish to the PDS of the active account.
    ///
    /// The publish is then durable. It continues after a restart, and it tries again with a longer
    /// delay after a temporary failure or an offline failure. The app does not lose it.
    ///
    /// A second entry for the same `entity_ref` replaces the first entry, and the newest record
    /// wins. The method returns [`AppError::NotAuthenticated`] when no account is active, because
    /// the queue needs the DID of the destination.
    pub(crate) async fn enqueue_publish(
        &self,
        kind: &str,
        entity_ref: &str,
        collection: &str,
        rkey: Option<&str>,
        value: serde_json::Value,
    ) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let entry = sync_outbox::NewOutboxEntry {
            account_did: did,
            kind: kind.to_string(),
            entity_ref: entity_ref.to_string(),
            collection: collection.to_string(),
            rkey: rkey.map(str::to_string),
            payload: serde_json::to_string(&value)?,
        };
        sync_outbox::enqueue(self.store.pool(), &entry, &Utc::now().to_rfc3339()).await?;
        Ok(())
    }

    /// The count of outbox rows for the active account that the app did not publish. The UI shows
    /// this count in its "N pending" indicator. The value is `0` when no account is active.
    pub async fn outbox_pending_count(&self) -> Result<i64, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(0);
        };
        Ok(sync_outbox::pending_count(self.store.pool(), &did).await?)
    }

    /// Each outbox row for the active account that is not complete. The rows have the state
    /// PENDING or FAILED. The sync detail view shows them.
    pub async fn outbox_entries(&self) -> Result<Vec<sync_outbox::OutboxEntry>, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(Vec::new());
        };
        Ok(sync_outbox::list(self.store.pool(), &did).await?)
    }

    /// The recent results of a publish for the active account. A result is a success or a failure.
    /// Together the results are the audit trail.
    pub async fn sync_history(&self, limit: i64) -> Result<Vec<sync_history::HistoryEntry>, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(Vec::new());
        };
        Ok(sync_history::recent(self.store.pool(), &did, limit).await?)
    }

    /// Remove orphan alignment records, which hold a coverage summary, from the PDS repository of
    /// the active account.
    ///
    /// These orphans are duplicates. An earlier `create` call chose the record key itself, and two
    /// calls could race. The repository then held two records for one alignment, and `sync_state`
    /// held only one of them.
    ///
    /// The method lists each alignment record in the repository. It removes a record when the rkey
    /// of that record has no source. A rkey has a source when `sync_state` holds it, or when it is
    /// the `aln-{id}` key of a local alignment that still exists.
    ///
    /// With `apply == false`, the method only reports the records that it would delete and changes
    /// nothing. The method returns the result.
    pub async fn prune_orphan_alignments(&self, apply: bool) -> Result<PruneReport, AppError> {
        let did = self.require_account()?;
        // The rkeys with a source. These are each rkey in sync_state for the alignment
        // collection, and the fixed key of each local alignment that still exists. The second group
        // keeps a row that the app did not yet publish.
        let mut keep: std::collections::HashSet<String> =
            sync_state::list_for_collection(self.store.pool(), &did, NS_ALIGNMENT)
                .await?
                .into_iter()
                .map(|s| s.rkey)
                .collect();
        for a in alignment::list_all(self.store.pool()).await? {
            keep.insert(alignment_rkey(a.id));
        }

        let mut engine = self.sync_engine()?;
        let mut report = PruneReport {
            applied: apply,
            ..PruneReport::default()
        };
        let mut cursor: Option<String> = None;
        loop {
            let (records, next) = engine.pull_list(NS_ALIGNMENT, cursor.as_deref()).await?;
            for r in &records {
                report.examined += 1;
                // rkey is the last path segment of at://did/collection/rkey.
                let rkey = r.uri.rsplit('/').next().unwrap_or_default().to_string();
                if rkey.is_empty() || keep.contains(&rkey) {
                    continue;
                }
                report.orphans.push(rkey.clone());
                if apply {
                    engine.push_delete(NS_ALIGNMENT, &rkey).await?;
                    report.deleted += 1;
                }
            }
            match next {
                Some(c) => cursor = Some(c),
                None => break,
            }
        }
        Ok(report)
    }

    /// Try to publish the ready outbox rows for the active account.
    ///
    /// After a success, the method writes a history row and removes the outbox row. After a
    /// temporary failure, it sets a later time for the row and stops the batch, because the machine
    /// is probably offline. The delay becomes longer after each try. After a permanent failure, it
    /// sets the row to `FAILED`.
    ///
    /// The method does nothing and returns `Ok` when no account is active. The caller can call it
    /// many times, at an interval and after a publish.
    pub async fn drain_outbox(&self) -> Result<DrainOutcome, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(DrainOutcome::default());
        };
        let mut outcome = DrainOutcome::default();
        // Build the resilient engine once (loads the session). Signed-out / no session → nothing to do.
        let mut engine = match self.sync_engine() {
            Ok(e) => e,
            Err(_) => return Ok(outcome),
        };
        let now = Utc::now();
        let batch = sync_outbox::ready(self.store.pool(), &did, &now.to_rfc3339(), OUTBOX_BATCH).await?;
        for entry in batch {
            let value: serde_json::Value = serde_json::from_str(&entry.payload)?;
            // If the app published this entity before, change the record that the PDS holds. Use
            // putRecord at the rkey that the app kept. Do not make a duplicate record.
            let known = sync_state::get(self.store.pool(), &did, &entry.entity_ref).await?;
            let result = match (&known, &entry.rkey) {
                (Some(ss), _) => engine.push_put(&entry.collection, &ss.rkey, value).await,
                (None, Some(rk)) => engine.push_put(&entry.collection, rk, value).await,
                (None, None) => engine.push_create(&entry.collection, value).await,
            };
            let attempt = entry.attempt_count + 1;
            match result {
                Ok(rref) => {
                    // Record the PDS-assigned identity + payload fingerprint, so the next publish
                    // updates this record and a PULL can detect divergence.
                    let state = sync_state::StoredSyncState {
                        account_did: did.clone(),
                        entity_ref: entry.entity_ref.clone(),
                        kind: entry.kind.clone(),
                        collection: entry.collection.clone(),
                        rkey: rref.rkey().to_string(),
                        at_uri: rref.uri.clone(),
                        at_cid: rref.cid.clone(),
                        payload_hash: sha256_str(&entry.payload),
                        pushed_at: now.to_rfc3339(),
                    };
                    sync_state::upsert(self.store.pool(), &state).await?;
                    self.log_history(&entry, "SUCCESS", Some(&rref), attempt, None).await?;
                    sync_outbox::complete(self.store.pool(), entry.id).await?;
                    outcome.published.push((entry.kind.clone(), rref.uri));
                }
                Err(e) if e.is_transient() => {
                    // The machine is offline, or the server gave a 5xx or a timeout. Wait, and
                    // stop. The other rows of the batch also wait.
                    let next = now + chrono::Duration::seconds(backoff_secs(attempt));
                    sync_outbox::reschedule(
                        self.store.pool(),
                        entry.id,
                        attempt,
                        &next.to_rfc3339(),
                        &e.to_string(),
                        &now.to_rfc3339(),
                    )
                    .await?;
                    outcome.retry_scheduled += 1;
                    break;
                }
                Err(e) => {
                    // Validation / auth / other terminal error: give up on this row (visible as FAILED).
                    self.log_history(&entry, "FAILED", None, attempt, Some(&e.to_string()))
                        .await?;
                    sync_outbox::mark_failed(self.store.pool(), entry.id, attempt, &e.to_string(), &now.to_rfc3339())
                        .await?;
                    outcome.failed += 1;
                }
            }
        }
        outcome.pending = sync_outbox::pending_count(self.store.pool(), &did).await?;
        Ok(outcome)
    }

    /// Add a sync-history row for a push that is complete.
    async fn log_history(
        &self,
        entry: &sync_outbox::OutboxEntry,
        status: &str,
        rref: Option<&RecordRef>,
        attempt_count: i64,
        error: Option<&str>,
    ) -> Result<(), AppError> {
        let h = sync_history::NewHistoryEntry {
            account_did: entry.account_did.clone(),
            kind: entry.kind.clone(),
            entity_ref: entry.entity_ref.clone(),
            collection: entry.collection.clone(),
            status: status.to_string(),
            at_uri: rref.map(|r| r.uri.clone()),
            at_cid: rref.map(|r| r.cid.clone()),
            attempt_count,
            error: error.map(str::to_string),
        };
        sync_history::record(self.store.pool(), &h, &Utc::now().to_rfc3339()).await?;
        Ok(())
    }

    /// Do a **PULL reconcile** (gap §5-p2). The method reads the records of the account from the
    /// PDS and compares them with the records that the app published, which `sync_state` holds. The
    /// policy is last-write-wins, and the remote copy has authority.
    ///
    /// The app knows a record by the rkey that it kept. When such a record changed on the PDS, the
    /// method applies the remote values to the local record where the data model permits it. Today
    /// this is the sex and the center of a biosample. The method then stores the new CID.
    ///
    /// The method marks a record for a new publish when the PDS no longer holds it. It counts a
    /// remote record that has no local record. A federated record is a summary with no personal
    /// data. It holds no local guid, so the app can not make a local entity from it.
    pub async fn pull_sync(&self) -> Result<PullOutcome, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        if did.starts_with("did:key:") {
            // A local did:key identity has no PDS repository. A PULL and a publish need an OAuth
            // account with a did:plc identity.
            return Err(AppError::Import(
                "PDS sync needs a signed-in PDS account — the local did:key identity has no PDS repo".into(),
            ));
        }
        let mut engine = self.sync_engine()?;
        let mut out = PullOutcome::default();
        for &collection in PUBLISHED_COLLECTIONS {
            // Page through the account's records in this collection.
            let mut remote = Vec::new();
            let mut cursor: Option<String> = None;
            loop {
                let (recs, next) = engine
                    .pull_list(collection, cursor.as_deref())
                    .await
                    .map_err(AppError::Sync)?;
                remote.extend(recs);
                match next {
                    Some(c) => cursor = Some(c),
                    None => break,
                }
            }
            let local: Vec<_> = sync_state::list_for_collection(self.store.pool(), &did, collection)
                .await?
                .into_iter()
                .map(|s| (s, None)) // local-hash recompute is future work — treat local as clean for now
                .collect();
            for action in sync_reconcile::plan(&local, &remote) {
                use sync_reconcile::ReconcileAction::*;
                match action {
                    InSync { .. } => out.in_sync += 1,
                    RePush { .. } => out.repushed += 1,
                    AdoptRemote { .. } => out.adopted += 1,
                    ApplyRemote {
                        entity_ref,
                        collection,
                        remote,
                        conflict,
                    } => {
                        self.apply_remote(&collection, &entity_ref, &remote.value).await?;
                        self.track_remote(&did, &entity_ref, &remote).await?;
                        out.applied += 1;
                        if conflict {
                            self.log_conflict(&did, &entity_ref, &collection).await?;
                            out.conflicts += 1;
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    /// Apply a remote record to the local state.
    ///
    /// The method applies only the values that the user can edit and that the local store owns. A
    /// federated record has no personal data and carries only those values. The app calculates a
    /// derived summary again on this machine. So the app tracks such a collection but does not
    /// write to it.
    pub(crate) async fn apply_remote(
        &self,
        collection: &str,
        entity_ref: &str,
        value: &serde_json::Value,
    ) -> Result<(), AppError> {
        if collection == NS_BIOSAMPLE {
            if let Some(guid) = entity_ref
                .strip_prefix("biosample:")
                .and_then(|s| Uuid::parse_str(s).ok())
                .map(SampleGuid)
            {
                if let Some(bio) = biosample::get(self.store.pool(), guid).await? {
                    let sex = value.get("sex").and_then(|v| v.as_str()).map(String::from).or(bio.sex);
                    let center = value
                        .get("center_name")
                        .or_else(|| value.get("centerName"))
                        .and_then(|v| v.as_str())
                        .map(String::from)
                        .or(bio.center_name);
                    self.update_biosample(
                        guid,
                        bio.donor_identifier,
                        bio.sample_accession,
                        bio.description,
                        center,
                        sex,
                    )
                    .await?;
                }
            }
        }
        Ok(())
    }

    /// Re-track a reconciled record's PDS identity so the next PULL sees it in sync.
    async fn track_remote(
        &self,
        did: &str,
        entity_ref: &str,
        remote: &navigator_sync::RemoteRecord,
    ) -> Result<(), AppError> {
        if let Some(mut ss) = sync_state::get(self.store.pool(), did, entity_ref).await? {
            ss.at_cid = remote.cid.clone();
            ss.at_uri = remote.uri.clone();
            ss.payload_hash = sha256_str(&remote.value.to_string());
            ss.pushed_at = Utc::now().to_rfc3339();
            sync_state::upsert(self.store.pool(), &ss).await?;
        }
        Ok(())
    }

    /// Log a both-sides-diverged conflict (remote won) to the sync history.
    async fn log_conflict(&self, did: &str, entity_ref: &str, collection: &str) -> Result<(), AppError> {
        let h = sync_history::NewHistoryEntry {
            account_did: did.to_string(),
            kind: "pull".into(),
            entity_ref: entity_ref.to_string(),
            collection: collection.to_string(),
            status: "RESOLVED_REMOTE".into(),
            at_uri: None,
            at_cid: None,
            attempt_count: 0,
            error: Some("local and remote both changed since last push; remote applied".into()),
        };
        sync_history::record_dir(self.store.pool(), &h, "CONFLICT", &Utc::now().to_rfc3339()).await?;
        Ok(())
    }

    /// Read the Ed25519 **device key** of this installation. At the first use, the method makes
    /// the key and publishes it.
    ///
    /// This key signs a call from the edge to the AppView. Federated IBD uses it today, and later
    /// the full signed surface will use it. The OS keychain holds the seed of the key under the
    /// active DID.
    ///
    /// The method publishes the public half of the key one time to the PDS of the user, as a
    /// [`DEVICE_KEY_COLLECTION`] record. The AppView reads that record through Jetstream, and it
    /// can then check our signatures.
    ///
    /// A second call is safe. The record key is the `did:key` value itself. So a second publish
    /// replaces the record and does not add a duplicate, and the method does not change a record
    /// that already exists. The method returns [`AppError::NotAuthenticated`] when no account is
    /// active.
    ///
    /// The method does *not* wait for the AppView to read the record. A signed call to the AppView
    /// absorbs the delay from 403 to 200 with a limited count of tries. See the IBD client.
    pub async fn ensure_device_key(&self) -> Result<DeviceKey, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = DeviceKey::load_or_generate(KEYCHAIN_SERVICE, &did)?;

        // A local did:key identity certifies itself, because the AppView checks the signature
        // against the DID. So there is no PDS record to publish. There is also no OAuth session
        // for such a publish.
        if did.starts_with("did:key:") {
            return Ok(key);
        }

        // Publish the public key once. A public getRecord on the deterministic rkey tells us
        // whether it already exists; only create it when absent (keeps re-launches quiet).
        let rkey = key.record_rkey();
        let session = self.auth.tokens.load(&did)?.ok_or(AppError::NotAuthenticated)?;
        let client = PdsClient::from_session(self.auth.http.clone(), &session)?;
        let already_published = client.get_record(DEVICE_KEY_COLLECTION, &rkey).await.is_ok();
        if !already_published {
            let record = serde_json::json!({
                "publicKey": key.did_key(),
                "createdAt": Utc::now().to_rfc3339(),
            });
            let mut engine = self.sync_engine()?;
            engine.push_create_rkey(DEVICE_KEY_COLLECTION, record, &rkey).await?;
        }
        Ok(key)
    }

    /// Federated IBD, **Step 1**. Fetch this account's pseudonymous match suggestions from
    /// the AppView (`GET /api/v1/ibd/suggestions`).
    ///
    /// The AppView reads the `fed.*` records that we published and makes a list of the best
    /// candidates. No genotype leaves the device in this step.
    ///
    /// To authenticate the call, the device key signs `"ibd-poll\n<DID>\n<ts>"`. The app registers
    /// that key at its first use. A 403 response directly after the first registration shows that
    /// the AppView did not yet read the device-key record. The client then tries again, and the
    /// delay becomes longer after each try.
    pub async fn ibd_suggestions(&self) -> Result<Vec<IbdSuggestion>, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let url = self.appview_url("ibd/suggestions");

        let mut attempt = 0u32;
        loop {
            let ts = Utc::now().timestamp().to_string();
            let sig = key.sign(&format!("ibd-poll\n{did}\n{ts}"));
            // reqwest applies URL encoding to a query value. So it escapes the `+`, `/`, and `=`
            // characters of a STANDARD base64 signature.
            let resp = self
                .auth
                .http
                .get(&url)
                .query(&[("did", did.as_str()), ("ts", ts.as_str()), ("sig", sig.as_str())])
                .send()
                .await
                .map_err(appview::transport)?;
            let status = resp.status();
            if status.is_success() {
                let body: serde_json::Value = resp.json().await.map_err(appview::transport)?;
                return Ok(parse_ibd_suggestions(&body));
            }
            if status.as_u16() == 403 && attempt < DEVICE_KEY_INGEST_RETRIES {
                tokio::time::sleep(std::time::Duration::from_secs(1u64 << attempt)).await;
                attempt += 1;
                continue;
            }
            return Err(appview::status_error("ibd/suggestions", resp).await);
        }
    }

    /// Federated IBD, **Step 2**. Request an introduction to a suggested candidate
    /// (`POST /api/v1/ibd/introduce`).
    ///
    /// The method signs `"ibd-introduce\n<DID>\n<suggested_sample_guid>"` and sends
    /// `{ did, suggestedSampleGuid, signature }`. It returns the `request_uri` of the AppView and
    /// the status, which is `PENDING`.
    ///
    /// This endpoint only opens the request. It exchanges no genetic data. After both parties
    /// agree, the consent messages and the encrypted segment exchange use the separate edge channel
    /// in `ibd_exchange`.
    pub async fn ibd_introduce(&self, suggested_sample_guid: &str) -> Result<IbdIntroResult, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = key.sign_fresh(ts, &format!("ibd-introduce\n{did}\n{suggested_sample_guid}"));
        // The IntroduceBody type of the AppView reads plain snake_case names and has no serde
        // rename. It parses the guid as a UUID. Send the guid exactly as the suggestion gives
        // it.
        let body = serde_json::json!({
            "did": did,
            "suggested_sample_guid": suggested_sample_guid,
            "ts": ts,
            "signature": sig,
        });
        let v = self.appview_post("ibd/introduce", body).await?;
        let request_uri = v
            .get("requestUri")
            .or_else(|| v.get("request_uri"))
            .and_then(|x| x.as_str())
            .unwrap_or_default()
            .to_string();
        let status = v
            .get("status")
            .and_then(|x| x.as_str())
            .unwrap_or("PENDING")
            .to_string();
        let purpose = v
            .get("purpose")
            .and_then(|x| x.as_str())
            .unwrap_or_default()
            .to_string();
        Ok(IbdIntroResult {
            request_uri,
            status,
            purpose,
        })
    }
}
