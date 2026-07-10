//! `impl App` methods extracted from `lib.rs` (the `ibd_exchange` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- IBD Phase 2: encrypted edge-to-edge exchange (D1 substrate) -------
    //
    // The AppView brokers discovery/consent + relays opaque ciphertext (never decrypts). These
    // wrap the `/api/v1/exchange/*` endpoints; the crypto (X25519/X3DH-lite/AES-GCM) lives in
    // `navigator_sync::exchange`. All calls are device-key-signed (no per-call OAuth).

    /// The signed-in account's X25519 identity key (load-or-generate), with its public half
    /// published to the AppView (`POST /exchange/key`, idempotent upsert) so partners can fetch it.
    pub async fn ensure_exchange_key(&self) -> Result<ExchangeKey, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ik = ExchangeKey::load_or_generate(KEYCHAIN_SERVICE, &did)?;
        let pub_b64 = ik.public_b64();
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &exchange::messages::publickey(&did, &pub_b64, None));
        let body = serde_json::json!({ "did": did, "x25519_pub": pub_b64, "ts": ts, "signature": sig });
        let v = self.exchange_post("exchange/key", body).await?;
        let _ = v; // { did, status: "published" }
        Ok(ik)
    }

    /// Fetch a peer's published X25519 public key (STANDARD base64), or `None` if they haven't
    /// published one. Public read — no signature.
    pub async fn fetch_exchange_key(&self, did: &str) -> Result<Option<String>, AppError> {
        let url = format!("{}/api/v1/exchange/key", decodingus_appview_url());
        let resp = self
            .auth
            .http
            .get(&url)
            .query(&[("did", did)])
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if resp.status().as_u16() == 404 {
            return Ok(None);
        }
        if !resp.status().is_success() {
            return Err(appview_status_error("exchange/key", resp).await);
        }
        let v: serde_json::Value = resp
            .json()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        Ok(v.get("x25519_pub").and_then(|x| x.as_str()).map(str::to_string))
    }

    /// Open an exchange request to a specific partner DID — the direct counterpart to the
    /// suggestion-mediated [`ibd_introduce`] (`POST /api/v1/exchange/request`). Generates an opaque
    /// request URI, signs the canonical request message, and posts it. The partner discovers it via
    /// [`exchange_incoming`] (symmetric-blind) and consents; on mutual consent a session opens. Returns
    /// the request URI to track. `scope` carries an optional project scope (team-ACL-gated server-side).
    pub async fn exchange_request(
        &self,
        partner_did: &str,
        purpose: &str,
        scope: Option<&str>,
    ) -> Result<String, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let request_uri = format!("exchange:{}", Uuid::new_v4());
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &exchange::messages::request(
            &request_uri,
            &did,
            partner_did,
            purpose,
            scope,
        ));
        let body = serde_json::json!({
            "request_uri": request_uri,
            "initiator_did": did,
            "partner_did": partner_did,
            "purpose": purpose,
            "scope": scope,
            "ts": ts,
            "signature": sig,
        });
        self.exchange_post("exchange/request", body).await?;
        Ok(request_uri)
    }

    /// Consent to (or decline) an exchange request. On mutual consent the AppView opens a session
    /// and returns its id.
    pub async fn exchange_consent(&self, request_uri: &str, given: bool) -> Result<ConsentOutcome, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &exchange::messages::consent(request_uri, &did, given));
        let body = serde_json::json!({
            "request_uri": request_uri,
            "consenting_did": did,
            "consent_given": given,
            "ts": ts,
            "signature": sig,
        });
        let v = self.exchange_post("exchange/consent", body).await?;
        Ok(ConsentOutcome {
            status: v
                .get("status")
                .and_then(|x| x.as_str())
                .unwrap_or("PENDING")
                .to_string(),
            session_id: v.get("session_id").and_then(|x| x.as_str()).map(str::to_string),
        })
    }

    /// Poll for inbound (symmetric-blind) exchange requests awaiting this account's consent.
    pub async fn exchange_incoming(&self) -> Result<Vec<IncomingRequest>, AppError> {
        let v = self.exchange_get_poll("exchange/incoming", &[]).await?;
        Ok(v.get("items")
            .and_then(|x| x.as_array())
            .map(|items| {
                items
                    .iter()
                    .map(|i| IncomingRequest {
                        request_uri: i
                            .get("request_uri")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        purpose: i
                            .get("purpose")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        created_at: i
                            .get("created_at")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                    })
                    .collect()
            })
            .unwrap_or_default())
    }

    /// Poll for consent-ready sessions (both parties consented; partner identity now revealed).
    pub async fn exchange_pending(&self) -> Result<Vec<ExchangeSessionInfo>, AppError> {
        let v = self.exchange_get_poll("exchange/pending", &[]).await?;
        Ok(v.get("items")
            .and_then(|x| x.as_array())
            .map(|items| {
                items
                    .iter()
                    .map(|i| ExchangeSessionInfo {
                        session_id: i
                            .get("session_id")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        request_uri: i
                            .get("request_uri")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        purpose: i
                            .get("purpose")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        partner_did: i
                            .get("partner_did")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        partner_key_uri: i.get("partner_key_uri").and_then(|x| x.as_str()).map(str::to_string),
                    })
                    .collect()
            })
            .unwrap_or_default())
    }

    /// Relay an opaque ciphertext `blob` to `to_did` in a session. The signed hash binds the blob to
    /// its routing (the broker stores ciphertext only). Returns the broker envelope id.
    pub async fn exchange_relay(&self, session_id: &str, to_did: &str, seq: i32, blob: &str) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let hash = exchange::blob_sha256_b64(blob).map_err(AppError::Sync)?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &exchange::messages::relay(session_id, &did, to_did, seq, &hash));
        let body = serde_json::json!({
            "session_id": session_id,
            "from_did": did,
            "to_did": to_did,
            "seq": seq,
            "blob": blob,
            "ts": ts,
            "signature": sig,
        });
        let v = self.exchange_post("exchange/relay", body).await?;
        Ok(v.get("id").and_then(|x| x.as_i64()).unwrap_or_default())
    }

    /// Pull undelivered relay envelopes for a session (ordered by seq).
    pub async fn exchange_relay_pull(&self, session_id: &str) -> Result<Vec<RelayEnvelope>, AppError> {
        let v = self
            .exchange_get_poll("exchange/relay/pull", &[("session_id", session_id)])
            .await?;
        Ok(v.get("items")
            .and_then(|x| x.as_array())
            .map(|items| {
                items
                    .iter()
                    .map(|i| RelayEnvelope {
                        id: i.get("id").and_then(|x| x.as_i64()).unwrap_or_default(),
                        from_did: i
                            .get("from_did")
                            .and_then(|x| x.as_str())
                            .unwrap_or_default()
                            .to_string(),
                        seq: i.get("seq").and_then(|x| x.as_i64()).unwrap_or_default() as i32,
                        blob: i.get("blob").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                    })
                    .collect()
            })
            .unwrap_or_default())
    }

    /// Acknowledge a delivered relay envelope (the broker drops it).
    pub async fn exchange_relay_ack(&self, envelope_id: i64) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &exchange::messages::ack(&did, envelope_id));
        let body = serde_json::json!({ "envelope_id": envelope_id, "did": did, "ts": ts, "signature": sig });
        self.exchange_post("exchange/ack", body).await.map(|_| ())
    }

    /// Establish a shared session key for a consent-ready session: publish/load our identity key,
    /// fetch the partner's, exchange ephemeral keys via the relay (handshake, seq 0), and derive the
    /// X3DH-lite session key. Polls the relay up to ~15s for the partner's handshake. The returned
    /// [`EstablishedSession`] then seals/opens payloads. (Live-only — needs a running AppView + the
    /// partner edge online to complete the handshake.)
    pub async fn open_exchange_session(&self, info: &ExchangeSessionInfo) -> Result<EstablishedSession, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let ik = self.ensure_exchange_key().await?;
        let partner_ik = self.fetch_exchange_key(&info.partner_did).await?.ok_or_else(|| {
            AppError::AppView(format!("partner {} has not published an X25519 key", info.partner_did))
        })?;

        let ek = exchange::EphemeralKey::generate();
        let hs = exchange::Envelope::handshake(&ek).to_blob().map_err(AppError::Sync)?;
        self.exchange_relay(&info.session_id, &info.partner_did, 0, &hs).await?;

        // Wait for the partner's handshake (seq 0 / a Handshake envelope), acking just it.
        let mut their_ek: Option<String> = None;
        for _ in 0..15 {
            for env in self.exchange_relay_pull(&info.session_id).await? {
                if let Ok(exchange::Envelope::Handshake { ek, .. }) = exchange::Envelope::from_blob(&env.blob) {
                    their_ek = Some(ek);
                    let _ = self.exchange_relay_ack(env.id).await;
                    break;
                }
            }
            if their_ek.is_some() {
                break;
            }
            tokio::time::sleep(std::time::Duration::from_secs(1)).await;
        }
        let their_ek =
            their_ek.ok_or_else(|| AppError::AppView("partner handshake not received (peer offline?)".into()))?;

        let key = exchange::derive_session_key(
            &ik,
            &ek,
            &partner_ik,
            &their_ek,
            exchange::role_is_a(&did, &info.partner_did),
        )
        .map_err(AppError::Sync)?;
        Ok(EstablishedSession {
            session_id: info.session_id.clone(),
            partner_did: info.partner_did.clone(),
            key,
        })
    }

    /// Seal `plaintext` and relay it on an established session (data starts at seq 1).
    pub async fn exchange_send(
        &self,
        session: &EstablishedSession,
        seq: i32,
        plaintext: &[u8],
    ) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let aad = exchange::relay_aad(&session.session_id, &did, &session.partner_did, seq);
        let blob = exchange::seal(&session.key, &aad, plaintext)
            .and_then(|e| e.to_blob())
            .map_err(AppError::Sync)?;
        self.exchange_relay(&session.session_id, &session.partner_did, seq, &blob)
            .await
    }

    /// Pull + decrypt + ack the data payloads waiting on an established session (returns plaintexts
    /// in pull order). Non-data / undecryptable envelopes are left un-acked.
    pub async fn exchange_receive(&self, session: &EstablishedSession) -> Result<Vec<Vec<u8>>, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let mut out = Vec::new();
        for env in self.exchange_relay_pull(&session.session_id).await? {
            let Ok(parsed) = exchange::Envelope::from_blob(&env.blob) else {
                continue;
            };
            // AAD binds the sender's routing: from = the partner (sender), to = us.
            let aad = exchange::relay_aad(&session.session_id, &env.from_did, &did, env.seq);
            if let Ok(pt) = exchange::open(&session.key, &aad, &parsed) {
                out.push(pt);
                let _ = self.exchange_relay_ack(env.id).await;
            }
        }
        Ok(out)
    }

    /// Run a **federated IBD exchange** over an established session (gap §4): send our IBD-panel
    /// dosages, receive the partner's, detect IBD locally (both peers run the symmetric detector →
    /// identical summary), then exchange + verify signed [`IbdAttestation`]s. `agreed` ⇒ the partner's
    /// signature verified and both summary hashes match. Only panel dosages cross the wire (encrypted;
    /// the broker never sees them). `my_source` supplies our dosages; the refs are opaque biosample
    /// pointers carried in the attestation. Live-only — needs the partner edge online.
    pub async fn exchange_ibd(
        &self,
        session: &EstablishedSession,
        my_source: IbdSource,
        request_uri: &str,
        my_sample_ref: Option<String>,
        partner_sample_ref: Option<String>,
        config: IbdDetectorConfig,
    ) -> Result<IbdExchangeResult, AppError> {
        let dosages = self.ibd_panel_dosages(my_source).await?;
        let sites = dosages
            .into_iter()
            .map(|g| IbdSite {
                contig: g.contig,
                position: g.position,
                dosage: g.dosage,
            })
            .collect();
        self.exchange_ibd_with_dosages(session, sites, request_uri, my_sample_ref, partner_sample_ref, config)
            .await
    }

    /// The dosage-level core of [`exchange_ibd`] — takes the panel dosages directly (e.g. from a
    /// consensus profile, or synthetic vectors in tests) rather than resolving an [`IbdSource`].
    pub async fn exchange_ibd_with_dosages(
        &self,
        session: &EstablishedSession,
        my_sites: Vec<IbdSite>,
        request_uri: &str,
        my_sample_ref: Option<String>,
        partner_sample_ref: Option<String>,
        config: IbdDetectorConfig,
    ) -> Result<IbdExchangeResult, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;

        // Fit the relay's 1 MiB envelope: decimate a large panel (both peers apply the same
        // position-based rule, so the intersection is preserved). Detect on the decimated set we send.
        let my_sites = decimate_for_exchange(my_sites);

        // 1. Send our dosages (the IBD panel is on CHM13 / hs1).
        let dos = IbdExchangeMsg::Dosages {
            build: "hs1".into(),
            sites: my_sites.clone(),
        };
        self.exchange_send(session, 1, &dos.to_bytes().map_err(AppError::Import)?)
            .await?;

        // 2. Receive the partner's dosages (buffering any attestation that arrives early).
        let mut partner_sites: Option<Vec<IbdSite>> = None;
        let mut partner_att: Option<IbdAttestation> = None;
        for _ in 0..EXCHANGE_POLL_ROUNDS {
            for pt in self.exchange_receive(session).await? {
                match IbdExchangeMsg::from_bytes(&pt) {
                    Ok(IbdExchangeMsg::Dosages { sites, .. }) => partner_sites = Some(sites),
                    Ok(IbdExchangeMsg::Attest(a)) => partner_att = Some(*a),
                    Err(_) => {}
                }
            }
            if partner_sites.is_some() {
                break;
            }
            tokio::time::sleep(std::time::Duration::from_secs(1)).await;
        }
        let partner_sites = partner_sites
            .ok_or_else(|| AppError::AppView("partner IBD dosages not received (peer offline?)".into()))?;

        // 3. Detect IBD locally (symmetric — the partner computes the same summary).
        let comparison = detect_ibd_sites(&my_sites, &partner_sites, ReferenceBuild::Chm13v2, config);

        // 4. Sign our attestation over the computed summary.
        let mut att = IbdAttestation::unsigned(
            request_uri,
            &session.session_id,
            &did,
            my_sample_ref,
            partner_sample_ref,
            &comparison.summary,
            Utc::now().to_rfc3339(),
        );
        att.signature = dev.sign(&att.canonical());
        att.signing_public_key = dev.did_key();

        // 5. Send our attestation.
        let att_msg = IbdExchangeMsg::Attest(Box::new(att.clone()));
        self.exchange_send(session, 2, &att_msg.to_bytes().map_err(AppError::Import)?)
            .await?;

        // 6. Receive the partner's attestation.
        for _ in 0..EXCHANGE_POLL_ROUNDS {
            if partner_att.is_some() {
                break;
            }
            for pt in self.exchange_receive(session).await? {
                if let Ok(IbdExchangeMsg::Attest(a)) = IbdExchangeMsg::from_bytes(&pt) {
                    partner_att = Some(*a);
                }
            }
            if partner_att.is_some() {
                break;
            }
            tokio::time::sleep(std::time::Duration::from_secs(1)).await;
        }
        let partner_att =
            partner_att.ok_or_else(|| AppError::AppView("partner attestation not received (peer offline?)".into()))?;

        // 7. Verify the partner's signature + summary-hash agreement.
        let sig_ok = du_atproto::verify_did_key(
            &partner_att.signing_public_key,
            partner_att.canonical().as_bytes(),
            &partner_att.signature,
        )
        .is_ok();
        let agreed = sig_ok && partner_att.summary_hash == att.summary_hash;

        Ok(IbdExchangeResult {
            summary: comparison.summary,
            segments: comparison.segments,
            overlapping_sites: comparison.overlapping_sites,
            my_attestation: att,
            partner_attestation: partner_att,
            agreed,
        })
    }

    /// The subject's best IBD data source: the highest-coverage alignment (the
    /// [`default_alignment_for_subject`](Self::default_alignment_for_subject) pattern), else its first
    /// chip profile. `Ok(None)` when the subject has neither.
    pub async fn best_ibd_source_for_subject(&self, guid: SampleGuid) -> Result<Option<IbdSource>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), guid).await?;
        let mut best: Option<(f64, i64)> = None;
        for a in &alignments {
            let cov = self
                .cached_coverage(a.id)
                .await
                .ok()
                .flatten()
                .map(|c| c.mean_coverage)
                .unwrap_or(0.0);
            if best.as_ref().map_or(true, |(bc, _)| cov > *bc) {
                best = Some((cov, a.id));
            }
        }
        if let Some((_, id)) = best {
            return Ok(Some(IbdSource::Alignment(id)));
        }
        Ok(self
            .list_chip_profiles(guid)
            .await?
            .first()
            .map(|c| IbdSource::Chip(c.id)))
    }

    /// The subject's IBD-panel dosages from its best source (panel-restricted — only the canonical IBD
    /// sites, not the whole genome, so that's all that can leave the device).
    pub async fn ibd_dosages_for_subject(&self, guid: SampleGuid) -> Result<Vec<SiteGenotype>, AppError> {
        let source = self.best_ibd_source_for_subject(guid).await?.ok_or_else(|| {
            AppError::Import("no IBD-capable data for this subject (need an alignment or a chip profile)".into())
        })?;
        self.ibd_panel_dosages(source).await
    }

    /// Run a federated IBD exchange for a **subject** (resolves its real IBD-panel dosages), persists
    /// the result, and best-effort publishes our attestation. The UI entry point.
    pub async fn exchange_ibd_for_subject(
        &self,
        session: &EstablishedSession,
        guid: SampleGuid,
        request_uri: &str,
        partner_sample_ref: Option<String>,
        config: IbdDetectorConfig,
    ) -> Result<IbdExchangeResult, AppError> {
        let dosages = self.ibd_dosages_for_subject(guid).await?;
        let sites = dosages
            .into_iter()
            .map(|g| IbdSite {
                contig: g.contig,
                position: g.position,
                dosage: g.dosage,
            })
            .collect();
        let result = self
            .exchange_ibd_with_dosages(
                session,
                sites,
                request_uri,
                Some(guid.to_string()),
                partner_sample_ref,
                config,
            )
            .await?;
        self.record_ibd_exchange(guid, session, request_uri, &result).await?;
        // Best-effort: publish our attestation to the PDS (skipped for did:key; never fails the exchange).
        let _ = self.publish_ibd_attestation(&result.my_attestation).await;
        Ok(result)
    }

    /// Persist an exchange result (upsert by session id).
    pub(crate) async fn record_ibd_exchange(
        &self,
        guid: SampleGuid,
        session: &EstablishedSession,
        request_uri: &str,
        r: &IbdExchangeResult,
    ) -> Result<(), AppError> {
        let row = navigator_store::ibd_exchange::StoredIbdExchange {
            session_id: session.session_id.clone(),
            request_uri: request_uri.to_string(),
            my_did: self.current_account().unwrap_or_default(),
            partner_did: session.partner_did.clone(),
            biosample_guid: guid.to_string(),
            partner_sample_ref: r.partner_attestation.attesting_sample_ref.clone(),
            total_shared_cm: r.summary.total_shared_cm,
            segment_count: r.summary.segment_count as i64,
            longest_segment_cm: r.summary.longest_segment_cm,
            relationship: format!("{:?}", r.summary.relationship),
            agreed: r.agreed,
            segments: serde_json::to_string(&r.segments).unwrap_or_else(|_| "[]".into()),
            my_attestation: serde_json::to_string(&r.my_attestation).unwrap_or_else(|_| "{}".into()),
            partner_attestation: serde_json::to_string(&r.partner_attestation).unwrap_or_else(|_| "{}".into()),
            created_at: Utc::now().to_rfc3339(),
        };
        navigator_store::ibd_exchange::upsert(self.store.pool(), &row).await?;
        Ok(())
    }

    /// All persisted IBD exchange results (newest first).
    pub async fn list_ibd_exchanges(&self) -> Result<Vec<StoredIbdExchange>, AppError> {
        Ok(navigator_store::ibd_exchange::list(self.store.pool()).await?)
    }

    /// Persisted IBD exchange results for one subject (newest first).
    pub async fn list_ibd_exchanges_for_subject(&self, guid: SampleGuid) -> Result<Vec<StoredIbdExchange>, AppError> {
        Ok(navigator_store::ibd_exchange::list_for_biosample(self.store.pool(), guid).await?)
    }

    /// Publish a signed attestation to the PDS (the AppView indexes it via Jetstream). No-op for a
    /// did:key local identity (self-certifying, no repo to write). Idempotent via a session-derived rkey.
    pub async fn publish_ibd_attestation(&self, att: &IbdAttestation) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        if did.starts_with("did:key:") {
            return Ok(()); // did:key has no PDS repo
        }
        let rkey: String = att
            .session_id
            .chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .take(64)
            .collect();
        let record = serde_json::to_value(att).map_err(|e| AppError::Import(e.to_string()))?;
        let mut engine = self.sync_engine()?;
        engine
            .push_create_rkey(IBD_ATTESTATION_COLLECTION, record, &rkey)
            .await?;
        Ok(())
    }

    /// POST a JSON body to an `/api/v1/<path>` exchange endpoint, mapping non-2xx to an AppView error.
    async fn exchange_post(&self, path: &str, body: serde_json::Value) -> Result<serde_json::Value, AppError> {
        let url = format!("{}/api/v1/{path}", decodingus_appview_url());
        let resp = self
            .auth
            .http
            .post(&url)
            .json(&body)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error(path, resp).await);
        }
        resp.json()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))
    }

    /// Issue a device-key-signed `exchange-poll` GET to an `/api/v1/<path>` endpoint, with `extra`
    /// query params appended. Shared by incoming / pending / relay-pull.
    async fn exchange_get_poll(&self, path: &str, extra: &[(&str, &str)]) -> Result<serde_json::Value, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/{path}", decodingus_appview_url());
        let ts = Utc::now().timestamp();
        let sig = dev.sign(&exchange::messages::poll(&did, ts));
        let ts_s = ts.to_string();
        let mut query: Vec<(&str, &str)> = vec![("did", did.as_str()), ("ts", ts_s.as_str()), ("sig", sig.as_str())];
        query.extend_from_slice(extra);
        let resp = self
            .auth
            .http
            .get(&url)
            .query(&query)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error(path, resp).await);
        }
        resp.json()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))
    }

    /// Enqueue the anchor records every child record references: the subject's biosample summary
    /// and each of its sequence runs, at their **deterministic** rkeys so the at:// URIs resolve.
    /// Idempotent — the outbox coalesces per `entity_ref` and re-publishing overwrites in place — so
    /// child publishes call this freely to guarantee there's always a biosample to tie back to.
    async fn ensure_subject_anchor(&self, did: &str, biosample_guid: SampleGuid) -> Result<(), AppError> {
        // Sequence runs first (the biosample record links to them).
        for run in self.list_sequence_runs(biosample_guid).await? {
            let value = self.sequence_run_record(did, &run).await?;
            self.enqueue_publish(
                "seqrun",
                &format!("seqrun:{}", run.id),
                NS_SEQUENCERUN,
                Some(&seqrun_rkey(run.id)),
                value,
            )
            .await?;
        }
        let value = self.biosample_record(did, biosample_guid).await?;
        self.enqueue_publish(
            "biosample",
            &format!("biosample:{biosample_guid}"),
            NS_BIOSAMPLE,
            Some(&biosample_rkey(biosample_guid)),
            value,
        )
        .await
    }

    /// Publish the alignment's coverage summary to the signed-in account's PDS (with
    /// refresh-on-expiry and retry/backoff via [`AsyncSync`]). Anchors the subject first so the
    /// record's biosample/sequence-run refs resolve.
    pub async fn publish_coverage(&self, alignment_id: i64) -> Result<(), AppError> {
        let did = self.require_account()?; // auth check before touching the DB
        let guid = self.biosample_of_alignment(alignment_id).await?;
        self.ensure_subject_anchor(&did, guid).await?;
        let value = self.coverage_record(&did, alignment_id).await?;
        self.enqueue_publish(
            "coverage",
            &format!("alignment:{alignment_id}"),
            NS_ALIGNMENT,
            // Deterministic rkey → the idempotent put path (never a fresh create), so re-publishing
            // or two concurrent drains converge on one record instead of duplicating.
            Some(&alignment_rkey(alignment_id)),
            value,
        )
        .await
    }

    /// Publish a subject's **consensus** ancestry estimate to the signed-in account's PDS — one
    /// populationBreakdown record per method (ADMIXTURE / PCA_PROJECTION_GMM / FINE_ADMIXTURE /
    /// G25_NMONTE), each linked to the biosample. Subject-level (the breakdown is computed from the
    /// pooled autosomal consensus, not per alignment), so one authoritative record set per subject
    /// rather than a conflicting set per sequencing run. The researcher opt-in act for the ancestry
    /// section — anonymized population proportions only.
    pub async fn publish_ancestry(&self, biosample_guid: SampleGuid) -> Result<(), AppError> {
        let did = self.require_account()?; // auth check before touching the DB
        self.ensure_subject_anchor(&did, biosample_guid).await?; // the breakdown links back to it
        let biosample_ref = biosample_at_uri(&did, biosample_guid);
        // One outbox row per method, keyed by subject+method so re-publishing coalesces per estimate.
        for r in &self.consensus_ancestry_results(biosample_guid).await? {
            let value =
                serde_json::to_value(population_breakdown_record(r).with_biosample_ref(Some(biosample_ref.clone())))?;
            let entity_ref = format!("ancestry:{biosample_guid}:{}", r.method);
            self.enqueue_publish("ancestry", &entity_ref, NS_POPULATION_BREAKDOWN, None, value)
                .await?;
        }
        Ok(())
    }

    /// Publish the anonymized biosample summary (sex, haplogroups) **and its sequence runs** to the
    /// signed-in account's PDS — the subject anchor every derived record ties back to. Deterministic
    /// rkeys make it idempotent (a re-publish overwrites rather than duplicating).
    pub async fn publish_biosample(&self, biosample_guid: SampleGuid) -> Result<(), AppError> {
        let did = self.require_account()?; // auth check before touching the DB
        self.ensure_subject_anchor(&did, biosample_guid).await
    }

    /// Publish a single sequence-run characterization to the signed-in account's PDS.
    pub async fn publish_sequence_run(&self, run: &SequenceRun) -> Result<(), AppError> {
        let did = self.require_account()?; // auth check before touching the DB
        let value = self.sequence_run_record(&did, run).await?;
        self.enqueue_publish(
            "seqrun",
            &format!("seqrun:{}", run.id),
            NS_SEQUENCERUN,
            Some(&seqrun_rkey(run.id)),
            value,
        )
        .await
    }

    /// Publish the alignment's de-novo calls for `contig` to the signed-in account's PDS.
    pub async fn publish_variants(&self, alignment_id: i64, contig: &str) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self.variants_record(alignment_id, contig).await?;
        let entity_ref = format!("variants:{alignment_id}:{contig}");
        self.enqueue_publish("variants", &entity_ref, PRIVATE_VARIANTS_COLLECTION, None, value)
            .await
    }
}
