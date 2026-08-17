//! `impl App` methods for the **matching ledger** — the durable state behind federated-IBD
//! discovery and consent.
//!
//! The pieces this coordinates already existed: the AppView's candidate engine
//! ([`App::ibd_suggestions`]), the blind introduction broker ([`App::ibd_introduce`]), the
//! consent round-trip and encrypted channel (`ibd_exchange.rs`), and the stored result
//! (`navigator_store::ibd_exchange`). What was missing is the thread between them — every
//! in-flight request lived in UI memory, so a restart forgot that we had asked anyone anything.
//!
//! [`App::refresh_matching`] is the single reconcile: it adopts what the broker reports, advances
//! rows the broker has moved on, and never overwrites a decision we made locally.

use super::*;

/// Region tag an attestation is filed under, derived from the exchange purpose
/// (`ibd.ibd_discovery_index.match_region_type` is `AUTOSOMAL`/`X`/`Y`/`MT`).
fn region_type_for(purpose: &str) -> &str {
    match purpose {
        "IBD_Y" => "Y",
        "IBD_MT" => "MT",
        "IBD_X" => "X",
        _ => "AUTOSOMAL",
    }
}

/// The ledger row a fresh conversation starts as.
fn new_row(request_uri: &str, direction: MatchingDirection, purpose: &str, status: MatchingStatus) -> StoredIbdRequest {
    let now = Utc::now().to_rfc3339();
    StoredIbdRequest {
        request_uri: request_uri.to_string(),
        direction: direction.as_str().to_string(),
        purpose: purpose.to_string(),
        status: status.as_str().to_string(),
        partner_did: None,
        session_id: None,
        biosample_guid: None,
        my_sample_ref: None,
        partner_sample_ref: None,
        consent_given: None,
        consent_at: None,
        attested_at: None,
        last_error: None,
        created_at: now.clone(),
        updated_at: now,
    }
}

impl App {
    /// Every matching conversation, newest first, with its exchange result attached when it has one.
    pub async fn matching_entries(&self) -> Result<Vec<MatchingEntry>, AppError> {
        let rows = navigator_store::ibd_request::list(self.store.pool()).await?;
        let results = navigator_store::ibd_exchange::list(self.store.pool()).await?;
        Ok(rows
            .into_iter()
            .map(|r| {
                let result = results.iter().find(|x| x.request_uri == r.request_uri).cloned();
                MatchingEntry {
                    direction: MatchingDirection::parse(&r.direction),
                    status: MatchingStatus::parse(&r.status),
                    biosample_guid: r.biosample_guid.as_deref().and_then(parse_sample_guid),
                    request_uri: r.request_uri,
                    purpose: r.purpose,
                    partner_did: r.partner_did,
                    session_id: r.session_id,
                    my_sample_ref: r.my_sample_ref,
                    partner_sample_ref: r.partner_sample_ref,
                    consent_given: r.consent_given,
                    attested: r.attested_at.is_some(),
                    last_error: r.last_error,
                    created_at: r.created_at,
                    updated_at: r.updated_at,
                    result,
                }
            })
            .collect())
    }

    /// Reconcile the ledger with the broker, then return the full list.
    ///
    /// Three passes, in order of increasing knowledge:
    /// 1. `/exchange/incoming` — inbound requests we have not seen. Adopted with
    ///    `insert_if_absent`, so a request we already declined stays declined even though the
    ///    broker keeps listing it.
    /// 2. `/exchange/pending` — mutual consent happened: the partner DID and session id are now
    ///    known. Advances anything not already terminal (a completed exchange stays completed).
    /// 3. Stored results — a completed exchange marks its request `EXCHANGED`.
    ///
    /// A pass that fails does not abort the others: a broker hiccup should degrade the view, not
    /// empty it. The first error is returned once the local state is consistent.
    pub async fn refresh_matching(&self) -> Result<Vec<MatchingEntry>, AppError> {
        let mut first_err: Option<AppError> = None;

        match self.exchange_incoming().await {
            Ok(incoming) => {
                for r in incoming {
                    let mut row = new_row(
                        &r.request_uri,
                        MatchingDirection::Inbound,
                        &r.purpose,
                        MatchingStatus::AwaitingConsent,
                    );
                    if !r.created_at.is_empty() {
                        row.created_at = r.created_at.clone();
                    }
                    navigator_store::ibd_request::insert_if_absent(self.store.pool(), &row).await?;
                }
            }
            Err(e) => first_err = Some(e),
        }

        match self.exchange_pending().await {
            Ok(pending) => {
                for info in pending {
                    // An unknown session means the request was opened on another device (or the
                    // ledger predates it) — adopt it so the session is still runnable here.
                    let existing = navigator_store::ibd_request::get(self.store.pool(), &info.request_uri).await?;
                    let mut row = existing.unwrap_or_else(|| {
                        new_row(
                            &info.request_uri,
                            MatchingDirection::Inbound,
                            &info.purpose,
                            MatchingStatus::Ready,
                        )
                    });
                    if !MatchingStatus::parse(&row.status).is_terminal() {
                        row.status = MatchingStatus::Ready.as_str().to_string();
                    }
                    row.partner_did = Some(info.partner_did.clone());
                    row.session_id = Some(info.session_id.clone());
                    if row.purpose.is_empty() {
                        row.purpose = info.purpose.clone();
                    }
                    row.updated_at = Utc::now().to_rfc3339();
                    navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
                }
            }
            Err(e) => first_err = first_err.or(Some(e)),
        }

        for done in navigator_store::ibd_exchange::list(self.store.pool()).await? {
            let Some(mut row) = navigator_store::ibd_request::get(self.store.pool(), &done.request_uri).await? else {
                continue;
            };
            if MatchingStatus::parse(&row.status) == MatchingStatus::Exchanged {
                continue;
            }
            row.status = MatchingStatus::Exchanged.as_str().to_string();
            row.session_id = Some(done.session_id.clone());
            row.partner_did = Some(done.partner_did.clone());
            row.updated_at = Utc::now().to_rfc3339();
            navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        }

        let entries = self.matching_entries().await?;
        match first_err {
            Some(e) => Err(e),
            None => Ok(entries),
        }
    }

    /// Ask to be introduced to a candidate and record the conversation.
    ///
    /// Carries both AppView sample handles into the ledger: `target_sample_guid` (ours) and
    /// `suggested_sample_guid` (theirs) are the only two identifiers an attestation can be filed
    /// under, and the suggestion is the one place we ever see them.
    pub async fn request_introduction(
        &self,
        suggestion: &IbdSuggestion,
        biosample_guid: Option<SampleGuid>,
    ) -> Result<MatchingEntry, AppError> {
        let intro = self.ibd_introduce(&suggestion.suggested_sample_guid).await?;
        let mut row = new_row(
            &intro.request_uri,
            MatchingDirection::Outbound,
            &intro.purpose,
            MatchingStatus::Requested,
        );
        row.my_sample_ref = suggestion.target_sample_guid.clone();
        row.partner_sample_ref = Some(suggestion.suggested_sample_guid.clone());
        row.biosample_guid = biosample_guid.map(|g| g.to_string());
        navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        self.matching_entry(&intro.request_uri).await
    }

    /// Consent to (or decline) an inbound request, recording our decision durably. The decision is
    /// written whatever the broker says next — re-polling must never resurrect a request we
    /// turned down.
    pub async fn matching_consent(
        &self,
        request_uri: &str,
        given: bool,
        biosample_guid: Option<SampleGuid>,
    ) -> Result<MatchingEntry, AppError> {
        let outcome = self.exchange_consent(request_uri, given).await?;
        let mut row = navigator_store::ibd_request::get(self.store.pool(), request_uri)
            .await?
            .unwrap_or_else(|| {
                new_row(
                    request_uri,
                    MatchingDirection::Inbound,
                    "",
                    MatchingStatus::AwaitingConsent,
                )
            });
        row.consent_given = Some(given);
        row.consent_at = Some(Utc::now().to_rfc3339());
        row.status = if !given {
            MatchingStatus::Declined
        } else if outcome.session_id.is_some() {
            MatchingStatus::Ready
        } else {
            // Recorded, but the counterpart has not consented yet.
            MatchingStatus::Requested
        }
        .as_str()
        .to_string();
        if let Some(sid) = outcome.session_id {
            row.session_id = Some(sid);
        }
        if biosample_guid.is_some() {
            row.biosample_guid = biosample_guid.map(|g| g.to_string());
        }
        row.updated_at = Utc::now().to_rfc3339();
        navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        self.matching_entry(request_uri).await
    }

    /// Bind a conversation to the local subject whose dosages it will exchange.
    pub async fn set_matching_subject(&self, request_uri: &str, guid: SampleGuid) -> Result<(), AppError> {
        let Some(mut row) = navigator_store::ibd_request::get(self.store.pool(), request_uri).await? else {
            return Ok(());
        };
        row.biosample_guid = Some(guid.to_string());
        row.updated_at = Utc::now().to_rfc3339();
        navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        Ok(())
    }

    /// Record the AppView sample handles for a conversation that did not come from a suggestion
    /// (or whose suggestion predated the AppView returning our own handle). Without both, a
    /// completed comparison cannot be attested.
    pub async fn set_matching_sample_refs(
        &self,
        request_uri: &str,
        mine: Option<&str>,
        theirs: Option<&str>,
    ) -> Result<(), AppError> {
        let Some(mut row) = navigator_store::ibd_request::get(self.store.pool(), request_uri).await? else {
            return Ok(());
        };
        if mine.is_some() {
            row.my_sample_ref = mine.map(str::to_string);
        }
        if theirs.is_some() {
            row.partner_sample_ref = theirs.map(str::to_string);
        }
        row.updated_at = Utc::now().to_rfc3339();
        navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        Ok(())
    }

    /// Record that an exchange attempt failed, so the row reads as `FAILED` with the reason rather
    /// than sitting at `READY` forever.
    pub async fn record_matching_failure(&self, request_uri: &str, err: &str) -> Result<(), AppError> {
        let Some(mut row) = navigator_store::ibd_request::get(self.store.pool(), request_uri).await? else {
            return Ok(());
        };
        row.status = MatchingStatus::Failed.as_str().to_string();
        row.last_error = Some(err.to_string());
        row.updated_at = Utc::now().to_rfc3339();
        navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        Ok(())
    }

    /// Drop a conversation from the local ledger. The broker keeps its own record, so a still-live
    /// request can reappear on the next refresh — this forgets, it does not cancel.
    pub async fn forget_matching_request(&self, request_uri: &str) -> Result<(), AppError> {
        navigator_store::ibd_request::delete(self.store.pool(), request_uri).await?;
        Ok(())
    }

    /// Mark a conversation complete once its exchange result is stored, adopting the request if the
    /// ledger has never seen it (a session opened on another device, or one predating the ledger).
    pub(crate) async fn mark_matching_exchanged(
        &self,
        guid: SampleGuid,
        session: &EstablishedSession,
        request_uri: &str,
    ) -> Result<(), AppError> {
        let mut row = navigator_store::ibd_request::get(self.store.pool(), request_uri)
            .await?
            .unwrap_or_else(|| new_row(request_uri, MatchingDirection::Inbound, "", MatchingStatus::Ready));
        row.status = MatchingStatus::Exchanged.as_str().to_string();
        row.session_id = Some(session.session_id.clone());
        row.partner_did = Some(session.partner_did.clone());
        row.biosample_guid = Some(guid.to_string());
        row.last_error = None;
        row.updated_at = Utc::now().to_rfc3339();
        navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        Ok(())
    }

    /// One conversation by request URI.
    pub async fn matching_entry(&self, request_uri: &str) -> Result<MatchingEntry, AppError> {
        self.matching_entries()
            .await?
            .into_iter()
            .find(|e| e.request_uri == request_uri)
            .ok_or_else(|| AppError::AppView(format!("no matching request {request_uri}")))
    }

    /// Tell the AppView to stop suggesting a candidate (`POST /api/v1/ibd/dismiss`). The dismissal
    /// is kept server-side across recomputes, so it survives without any local mirror.
    pub async fn ibd_dismiss(&self, suggested_sample_guid: &str) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = key.sign_fresh(ts, &format!("ibd-dismiss\n{did}\n{suggested_sample_guid}"));
        let body = serde_json::json!({
            "did": did,
            "suggested_sample_guid": suggested_sample_guid,
            "ts": ts,
            "signature": sig,
        });
        self.appview_post("ibd/dismiss", body).await?;
        Ok(())
    }

    /// Report a completed comparison to the AppView (`POST /api/v1/ibd/attest`) — the step that
    /// turns a private, edge-computed match into a discovery signal. Only coarse totals travel:
    /// two opaque sample handles, a region tag, cM and segment count. Never coordinates, never
    /// genotypes.
    ///
    /// The signed `cm` is formatted `{:.1}` to match the AppView's own canonical string byte for
    /// byte; a mismatch there fails signature verification, not parsing.
    pub async fn ibd_attest(
        &self,
        request_uri: &str,
        claimed_sample: &str,
        counterpart_sample: &str,
        region_type: &str,
        total_shared_cm: f64,
        num_segments: i32,
    ) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let cm = format!("{total_shared_cm:.1}");
        let ts = Utc::now().timestamp();
        let sig = key.sign_fresh(
            ts,
            &format!("ibd-attest\n{did}\n{request_uri}\n{claimed_sample}\n{counterpart_sample}\n{region_type}\n{cm}"),
        );
        let body = serde_json::json!({
            "did": did,
            "request_uri": request_uri,
            "claimed_sample": claimed_sample,
            "counterpart_sample": counterpart_sample,
            "region_type": region_type,
            "total_shared_cm": total_shared_cm,
            "num_segments": num_segments,
            "ts": ts,
            "signature": sig,
        });
        self.appview_post("ibd/attest", body).await?;
        Ok(())
    }

    /// Attest a completed exchange if — and only if — it is attestable: both parties agreed on the
    /// summary, and we know both AppView sample handles.
    ///
    /// Agreement is the gate because the AppView confirms an edge only when both parties report a
    /// compatible total; filing a one-sided figure from a comparison our partner disputes would put
    /// a claim on the discovery graph that our own run says is wrong. Handles are missing whenever
    /// the conversation did not come from a suggestion (a direct request never names them), so this
    /// is a no-op rather than an error.
    pub(crate) async fn attest_exchange_if_possible(&self, request_uri: &str) -> Result<bool, AppError> {
        let entry = self.matching_entry(request_uri).await?;
        let (Some(mine), Some(theirs)) = (entry.my_sample_ref.clone(), entry.partner_sample_ref.clone()) else {
            return Ok(false);
        };
        let Some(result) = entry.result.as_ref().filter(|r| r.agreed) else {
            return Ok(false);
        };
        self.ibd_attest(
            request_uri,
            &mine,
            &theirs,
            region_type_for(&entry.purpose),
            result.total_shared_cm,
            result.segment_count as i32,
        )
        .await?;
        if let Some(mut row) = navigator_store::ibd_request::get(self.store.pool(), request_uri).await? {
            row.attested_at = Some(Utc::now().to_rfc3339());
            row.updated_at = Utc::now().to_rfc3339();
            navigator_store::ibd_request::upsert(self.store.pool(), &row).await?;
        }
        Ok(true)
    }
}

/// A stored guid string back into a `SampleGuid` (a malformed one just reads as unbound).
fn parse_sample_guid(s: &str) -> Option<SampleGuid> {
    Uuid::parse_str(s).ok().map(SampleGuid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_type_follows_the_exchange_purpose() {
        assert_eq!(region_type_for("IBD_AUTOSOMAL"), "AUTOSOMAL");
        assert_eq!(region_type_for("IBD_Y"), "Y");
        assert_eq!(region_type_for("IBD_MT"), "MT");
        assert_eq!(region_type_for("IBD_X"), "X");
        // An unknown/empty purpose must not invent a uniparental claim.
        assert_eq!(region_type_for(""), "AUTOSOMAL");
        assert_eq!(region_type_for("GENEALOGY_PII"), "AUTOSOMAL");
    }

    /// Cross-repo contract: these strings are signed, and the AppView rebuilds them byte for byte
    /// (`du_db::ibd::messages`). A drift here fails as a signature rejection, not a parse error.
    #[test]
    fn canonical_dismiss_and_attest_messages() {
        let did = "did:plc:abc123";
        assert_eq!(
            format!("ibd-dismiss\n{did}\n{}", "sample-xyz"),
            "ibd-dismiss\ndid:plc:abc123\nsample-xyz"
        );
        let cm = format!("{:.1}", 75.04_f64);
        assert_eq!(cm, "75.0", "cM is signed at one decimal place");
        assert_eq!(
            format!(
                "ibd-attest\n{did}\n{}\n{}\n{}\n{}\n{}",
                "urn:ibd:r", "s-mine", "s-theirs", "AUTOSOMAL", cm
            ),
            "ibd-attest\ndid:plc:abc123\nurn:ibd:r\ns-mine\ns-theirs\nAUTOSOMAL\n75.0"
        );
    }

    #[test]
    fn status_round_trips_and_marks_terminal() {
        for s in [
            MatchingStatus::Requested,
            MatchingStatus::AwaitingConsent,
            MatchingStatus::Declined,
            MatchingStatus::Ready,
            MatchingStatus::Exchanged,
            MatchingStatus::Failed,
        ] {
            assert_eq!(MatchingStatus::parse(s.as_str()), s);
        }
        // An unknown token degrades to the least-committal state, never to a terminal one.
        assert_eq!(MatchingStatus::parse("WAT"), MatchingStatus::Requested);
        assert!(MatchingStatus::Exchanged.is_terminal());
        assert!(MatchingStatus::Declined.is_terminal());
        assert!(!MatchingStatus::Ready.is_terminal());
        assert!(
            !MatchingStatus::Failed.is_terminal(),
            "a failure is retryable, not done"
        );
        assert_eq!(MatchingDirection::parse("INBOUND"), MatchingDirection::Inbound);
        assert_eq!(MatchingDirection::parse("OUTBOUND"), MatchingDirection::Outbound);
    }
}
