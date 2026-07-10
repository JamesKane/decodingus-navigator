//! `impl App` methods extracted from `lib.rs` (the `sync` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- sync durability: outbox enqueue + drain (gap §5) -------------------

    /// Enqueue a built record for publishing to the signed-in account's PDS. The publish becomes
    /// durable: it survives restart and retries automatically (with backoff) on a transient/offline
    /// failure instead of being lost. Re-enqueuing the same `entity_ref` coalesces (newest wins).
    /// Errors [`AppError::NotAuthenticated`] when signed out (we need the destination DID).
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

    /// Pending (not-yet-published) outbox rows for the signed-in account — drives the UI's
    /// "N pending" indicator. `0` when signed out.
    pub async fn outbox_pending_count(&self) -> Result<i64, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(0);
        };
        Ok(sync_outbox::pending_count(self.store.pool(), &did).await?)
    }

    /// All non-completed outbox rows (PENDING + FAILED) for the signed-in account — a sync detail view.
    pub async fn outbox_entries(&self) -> Result<Vec<sync_outbox::OutboxEntry>, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(Vec::new());
        };
        Ok(sync_outbox::list(self.store.pool(), &did).await?)
    }

    /// Recent publish outcomes (success/failure) for the signed-in account — the audit trail.
    pub async fn sync_history(&self, limit: i64) -> Result<Vec<sync_history::HistoryEntry>, AppError> {
        let Some(did) = self.current_account() else {
            return Ok(Vec::new());
        };
        Ok(sync_history::recent(self.store.pool(), &did, limit).await?)
    }

    /// Prune orphaned alignment (coverage-summary) records from the signed-in account's PDS repo —
    /// the duplicates left by the pre-deterministic-rkey `create` race (two records for one
    /// alignment, only one tracked in `sync_state`). Lists every alignment record in the repo and
    /// removes any whose rkey is **not accounted for** — i.e. neither a `sync_state`-tracked rkey nor
    /// the deterministic `aln-{id}` of a live local alignment. With `apply == false` it's a dry run
    /// (reports what it would delete, touches nothing). Returns the outcome.
    pub async fn prune_orphan_alignments(&self, apply: bool) -> Result<PruneReport, AppError> {
        let did = self.require_account()?;
        // Accounted-for rkeys: everything tracked in sync_state for the alignment collection, plus
        // the deterministic key for every live local alignment (so a not-yet-drained one isn't culled).
        let mut keep: std::collections::HashSet<String> = sync_state::list_for_collection(self.store.pool(), &did, NS_ALIGNMENT)
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

    /// Attempt to publish the ready outbox rows for the signed-in account. Each success is logged to
    /// history and its row removed; a transient failure reschedules the row with exponential backoff
    /// and stops the batch (we're likely offline); a non-transient failure marks the row `FAILED`.
    /// A no-op (and `Ok`) when signed out. Safe to call repeatedly (periodically + after a publish).
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
            // Idempotency: if we've published this entity before, update the PDS-assigned record in
            // place (putRecord at the kept rkey) instead of creating a duplicate.
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
                    // Offline / 5xx / timeout: back off and stop — the rest of the batch will wait too.
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

    /// Append a sync-history row for a finished push attempt.
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

    /// **PULL reconcile** (gap §5-p2): fetch the account's own records from the PDS and reconcile
    /// against what we published (`sync_state`), last-write-wins / remote-authoritative. For records
    /// we recognise (via the kept rkey) that changed on the PDS, apply remote→local where the data
    /// model allows (today: a biosample's sex / center) and re-track the CID. Records missing remotely
    /// are flagged for re-publish; remote records with no local mapping are counted (the fed records
    /// are PII-free *summaries* and carry no local guid, so they can't reconstruct a local entity).
    pub async fn pull_sync(&self) -> Result<PullOutcome, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        if did.starts_with("did:key:") {
            // A local did:key identity has no PDS repo — PULL/publish need a real OAuth (did:plc) account.
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

    /// Apply a remote record onto local state. Only the editable, locally-authoritative bits the
    /// PII-free fed record carries can be applied; derived-summary collections are recomputed locally,
    /// so they're tracked but not overwritten.
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

    /// Load (or, on first use, generate + publish) this installation's Ed25519 **device key**
    /// — the signing key that authenticates Edge→AppView calls (federated IBD and, later, the
    /// whole signed surface). The key seed lives in the OS keychain scoped to the signed-in
    /// DID; its public half is published once to the user's PDS as a
    /// [`DEVICE_KEY_COLLECTION`] record so the AppView (which ingests it via Jetstream) can
    /// verify our signatures. Idempotent: the record is keyed by its own `did:key`, so a
    /// re-publish overwrites rather than duplicates, and an already-present record is left
    /// alone. Errors [`AppError::NotAuthenticated`] when signed out.
    ///
    /// This does *not* wait for ingest — the signed AppView calls absorb the 403→200 lag with
    /// bounded retries (see the IBD client).
    pub async fn ensure_device_key(&self) -> Result<DeviceKey, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = DeviceKey::load_or_generate(KEYCHAIN_SERVICE, &did)?;

        // A local did:key identity self-certifies (the AppView verifies the signature against the DID
        // itself), so there is no PDS record to publish — and no OAuth session to do it with.
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

    /// Federated IBD — **Step 1**: fetch this account's pseudonymous match suggestions from
    /// the AppView (`GET /api/v1/ibd/suggestions`).
    ///
    /// The AppView mines our already-published `fed.*` records into a top-K candidate list;
    /// no genotypes leave the device here. The call is authenticated by signing
    /// `"ibd-poll\n<DID>\n<ts>"` with the device key (registered on first use). A 403 right
    /// after first-time registration means the AppView hasn't ingested the device-key record
    /// yet, so it's retried with exponential backoff.
    pub async fn ibd_suggestions(&self) -> Result<Vec<IbdSuggestion>, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/ibd/suggestions", decodingus_appview_url());

        let mut attempt = 0u32;
        loop {
            let ts = Utc::now().timestamp().to_string();
            let sig = key.sign(&format!("ibd-poll\n{did}\n{ts}"));
            // reqwest URL-encodes query values, so the STANDARD-base64 sig (`+` `/` `=`) is
            // safely escaped.
            let resp = self
                .auth
                .http
                .get(&url)
                .query(&[("did", did.as_str()), ("ts", ts.as_str()), ("sig", sig.as_str())])
                .send()
                .await
                .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
            let status = resp.status();
            if status.is_success() {
                let body: serde_json::Value = resp
                    .json()
                    .await
                    .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
                return Ok(parse_ibd_suggestions(&body));
            }
            if status.as_u16() == 403 && attempt < DEVICE_KEY_INGEST_RETRIES {
                tokio::time::sleep(std::time::Duration::from_secs(1u64 << attempt)).await;
                attempt += 1;
                continue;
            }
            return Err(appview_status_error("ibd/suggestions", resp).await);
        }
    }

    /// Federated IBD — **Step 2**: request an introduction to a suggested candidate
    /// (`POST /api/v1/ibd/introduce`).
    ///
    /// Signs `"ibd-introduce\n<DID>\n<suggested_sample_guid>"` and posts
    /// `{ did, suggestedSampleGuid, signature }`. Returns the AppView's `request_uri` and
    /// status (`PENDING`). The downstream consent round-trip + key exchange are deferred
    /// (gated on the AppView's symmetric-blind counterpart discovery), so this only opens the
    /// request — it does not exchange any genetic data.
    pub async fn ibd_introduce(&self, suggested_sample_guid: &str) -> Result<IbdIntroResult, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/ibd/introduce", decodingus_appview_url());
        let ts = Utc::now().timestamp();
        let sig = key.sign_fresh(ts, &format!("ibd-introduce\n{did}\n{suggested_sample_guid}"));
        // The AppView's IntroduceBody deserializes plain snake_case (no serde rename), and
        // parses the guid as a UUID — send it verbatim from the suggestion.
        let body = serde_json::json!({
            "did": did,
            "suggested_sample_guid": suggested_sample_guid,
            "ts": ts,
            "signature": sig,
        });
        let resp = self
            .auth
            .http
            .post(&url)
            .json(&body)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error("ibd/introduce", resp).await);
        }
        let v: serde_json::Value = resp
            .json()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
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
        Ok(IbdIntroResult { request_uri, status })
    }
}
