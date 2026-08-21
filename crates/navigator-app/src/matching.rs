//! `impl App` methods for the **matching ledger**. The ledger is the durable state behind
//! federated-IBD discovery and consent.
//!
//! The parts that this module coordinates already existed. They are the candidate engine of the
//! AppView ([`App::ibd_suggestions`]) and the blind introduction broker
//! ([`App::ibd_introduce`]). They are also the consent messages with the encrypted channel
//! (`ibd_exchange.rs`), and the stored result (`navigator_store::ibd_exchange`).
//!
//! The connection between these parts was absent. Each open request stayed in the memory of the
//! UI. So after a restart, the app did not know that it had sent a request to anybody.
//!
//! [`App::refresh_matching`] is the one reconcile. It adopts the state that the broker reports. It
//! advances a row when the broker moves that row forward. It never writes over a decision that the
//! user made on this machine.

use super::*;

/// The region tag of an attestation. The code takes the tag from the purpose of the exchange. The
/// field `ibd.ibd_discovery_index.match_region_type` holds `AUTOSOMAL`, `X`, `Y`, or `MT`.
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
    /// There are three passes. Each pass knows more than the pass before it.
    /// 1. `/exchange/incoming` gives the requests that arrived and that the app did not see. The
    ///    app adopts each one with `insert_if_absent`. So a request that the user declined stays
    ///    declined, and the broker can continue to list it.
    /// 2. `/exchange/pending` shows that both parties agreed. The partner DID and the session id
    ///    are now known. This pass advances each row that is not yet in a final state. A complete
    ///    exchange keeps its state.
    /// 3. The stored results set the request of a complete exchange to `EXCHANGED`.
    ///
    /// A pass that fails does not stop the other passes. A short broker fault must make the view
    /// less complete, but it must not empty the view. The method returns the first error after the
    /// local state is consistent.
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
                    // An unknown session shows that another device opened the request. The
                    // ledger can also be older than the request. Adopt the session, so the user
                    // can still run it on this machine.
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

    /// Ask the broker for an introduction to a candidate, and record the conversation.
    ///
    /// The method puts both AppView sample handles in the ledger. `target_sample_guid` is our
    /// handle and `suggested_sample_guid` is the handle of the partner. An attestation can use only
    /// these two identifiers. The suggestion is the one place where the app sees them.
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

    /// Agree to a request that arrived, or decline it, and write the decision to the store. The
    /// app writes the decision, and a later report from the broker does not change it. A new poll
    /// must never return a request that the user declined.
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

    /// Record the AppView sample handles of a conversation that has no suggestion. A suggestion
    /// can also be older than the AppView change that added our own handle. The app needs both
    /// handles. Without them, it can not attest a complete comparison.
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

    /// Record a failed try of an exchange. The row then shows `FAILED` with the reason. Without
    /// this record, the row stays at `READY` for all time.
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

    /// Remove a conversation from the local ledger. The broker keeps its own record. So an open
    /// request can come back at the next refresh. This method removes a local row. It does not
    /// cancel the request.
    pub async fn forget_matching_request(&self, request_uri: &str) -> Result<(), AppError> {
        navigator_store::ibd_request::delete(self.store.pool(), request_uri).await?;
        Ok(())
    }

    /// Mark a conversation as complete after the app stores the result of its exchange. If the
    /// ledger has no row for the request, the method adds one. Another device can open a session,
    /// and a session can also be older than the ledger.
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

    /// Tell the AppView to remove a candidate from its suggestions (`POST /api/v1/ibd/dismiss`).
    /// The server keeps this decision when it calculates the candidates again. So the app needs no
    /// local copy of the decision.
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

    /// Report a complete comparison to the AppView (`POST /api/v1/ibd/attest`). This step changes
    /// a private match, which the device calculated, into a discovery signal.
    ///
    /// Only approximate totals cross the network. They are two opaque sample handles, a region tag,
    /// the cM value, and the count of segments. A coordinate never crosses. A genotype never
    /// crosses.
    ///
    /// The code formats the signed `cm` value as `{:.1}`, which gives the same canonical string as
    /// the AppView. A difference here fails the signature check. It does not fail the parser.
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

    /// Attest a complete exchange, but only when the app can attest it. Two conditions apply. Both
    /// parties must agree on the summary, and the app must know both AppView sample handles.
    ///
    /// Agreement is a condition because the AppView confirms an edge only when both parties report
    /// a compatible total. The partner can dispute a comparison. A total from such a comparison
    /// puts a claim on the discovery graph that our own run says is wrong.
    ///
    /// The handles are absent when the conversation has no suggestion, because a direct request
    /// never names them. In that case the method does nothing and gives no error.
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

    /// A contract between two repositories. The device key signs these strings, and the AppView
    /// makes the same strings in `du_db::ibd::messages`. A change on one side only fails as a
    /// rejected signature. It does not fail as a parse error.
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
