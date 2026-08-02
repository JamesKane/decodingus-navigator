//! `impl NavigatorApp` — the **Matching** tab: federated-IBD discovery and consent.
//!
//! This is the front door for machinery that was already complete but had no coherent surface. It
//! is account-scoped, not subject-scoped: a conversation is keyed by our DID and the broker's
//! request URI, and a local subject is chosen only when it is time to exchange dosages. The
//! subject's own IBD tab keeps the *results* for that person; this tab owns the conversation.
//!
//! Three sub-tabs follow one conversation's life — a ranked candidate (Suggestions) becomes a
//! request awaiting consent (Requests) and then a result (Results).
use super::*;

impl NavigatorApp {
    /// The Matching work area. Gated on sign-in: every call here is device-key-signed.
    pub(crate) fn matching_central(&mut self, ui: &mut egui::Ui) {
        if self.account.is_none() {
            empty_state(
                ui,
                self.tr("matching.signedout.title"),
                self.tr("matching.signedout.hint"),
            );
            return;
        }
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.heading(self.tr("nav.matching"));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add_enabled(!self.exchange_busy, egui::Button::new(self.tr("common.refresh")))
                    .clicked()
                {
                    self.exchange_busy = true;
                    let _ = self.tx.send(Command::RefreshMatching);
                }
                if self.exchange_busy {
                    ui.spinner();
                }
            });
        });
        ui.label(egui::RichText::new(self.tr("matching.intro")).weak().small());
        ui.separator();
        ui.add_space(4.0);
        self.matching_subject_picker(ui);
        ui.add_space(4.0);
        self.matching_tab = self.sub_bar(ui, self.matching_tab, &MatchingTab::ALL);
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.add_space(4.0);
            match self.matching_tab {
                MatchingTab::Suggestions => self.matching_suggestions(ui),
                MatchingTab::Requests => self.matching_requests(ui),
                MatchingTab::Results => self.matching_results(ui),
            }
        });
    }

    /// Which local subject an exchange speaks for. Explicit here rather than implied by whichever
    /// subject tab happened to be open — the same account can hold several people's data, and
    /// sending the wrong one's genotypes is not a recoverable mistake.
    ///
    /// Shown as the current choice plus a *Change* toggle rather than a dropdown: a workspace can
    /// hold tens of thousands of subjects, and a `ComboBox` builds a widget per entry every frame
    /// its popup is open. The reveal is the same filter-then-virtualized-list the subjects rail
    /// uses, so the cost is the number of rows on screen, not the number in the workspace.
    fn matching_subject_picker(&mut self, ui: &mut egui::Ui) {
        if self.matching_subject.is_none() {
            self.matching_subject = self.selected_sample.or(self.all_biosamples.first().map(|b| b.guid));
        }
        let label = self
            .matching_subject
            .and_then(|g| self.find_subject(g).map(|b| b.donor_identifier.clone()))
            .unwrap_or_else(|| self.tr("matching.noSubject").to_string());
        ui.horizontal(|ui| {
            ui.label(self.tr("matching.exchangeAs"));
            ui.label(egui::RichText::new(label).strong());
            let toggle = if self.matching_subject_picking {
                self.tr("common.cancel")
            } else {
                self.tr("common.change")
            };
            if ui.button(toggle).clicked() {
                self.matching_subject_picking = !self.matching_subject_picking;
                self.matching_subject_filter.clear();
            }
        });
        if !self.matching_subject_picking {
            return;
        }

        ui.add_space(4.0);
        let hint = self.tr("subjects.filter");
        ui.add(
            egui::TextEdit::singleline(&mut self.matching_subject_filter)
                .hint_text(hint)
                .desired_width(280.0),
        );
        // Build the filtered view from immutable reads first, so the scroll closure borrows locals.
        let needle = self.matching_subject_filter.trim().to_lowercase();
        let rows: Vec<(SampleGuid, String)> = self
            .all_biosamples
            .iter()
            .filter(|b| needle.is_empty() || b.donor_identifier.to_lowercase().contains(&needle))
            .map(|b| (b.guid, b.donor_identifier.clone()))
            .collect();
        ui.label(egui::RichText::new(format!("{}", rows.len())).weak().small());
        if rows.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.noMatch")).weak());
            return;
        }

        let selected = self.matching_subject;
        let mut pick = None;
        let row_h = ui.spacing().interact_size.y;
        egui::ScrollArea::vertical()
            .id_salt("matching_subject_list")
            .max_height(180.0)
            .auto_shrink([false, false])
            .show_rows(ui, row_h, rows.len(), |ui, range| {
                for i in range {
                    let (guid, name) = &rows[i];
                    if ui.selectable_label(selected == Some(*guid), name).clicked() {
                        pick = Some(*guid);
                    }
                }
            });
        if let Some(guid) = pick {
            self.matching_subject = Some(guid);
            self.matching_subject_picking = false;
            self.matching_subject_filter.clear();
        }
    }

    /// Ranked candidates from the AppView's engine. Pseudonymous: a candidate is an opaque sample
    /// handle plus the signals behind its score — never a DID, never a name.
    fn matching_suggestions(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !self.loading_ibd_suggestions,
                    egui::Button::new(self.tr("network.find")),
                )
                .clicked()
            {
                self.loading_ibd_suggestions = true;
                self.status = self.tr("network.finding").to_string();
                let _ = self.tx.send(Command::LoadIbdSuggestions);
            }
            if self.loading_ibd_suggestions {
                ui.spinner();
            }
        });
        ui.label(egui::RichText::new(self.tr("network.note")).weak().small());

        // Requesting an introduction and dismissing both remove a row, so filter against what the
        // ledger already knows rather than trusting the fetched list to be current.
        let requested: std::collections::HashSet<String> = self
            .matching
            .iter()
            .filter_map(|e| e.partner_sample_ref.clone())
            .collect();
        let rows: Vec<navigator_app::IbdSuggestion> = self
            .ibd_suggestions
            .iter()
            .filter(|s| {
                !requested.contains(&s.suggested_sample_guid)
                    && !self.dismissed_candidates.contains(&s.suggested_sample_guid)
            })
            .cloned()
            .collect();
        if rows.is_empty() {
            if !self.loading_ibd_suggestions {
                ui.add_space(6.0);
                ui.weak(self.tr("network.empty"));
            }
            return;
        }

        ui.add_space(6.0);
        let mut introduce: Option<navigator_app::IbdSuggestion> = None;
        let mut dismiss: Option<String> = None;
        egui::Grid::new("matching_suggestions")
            .striped(true)
            .num_columns(5)
            .spacing([14.0, 4.0])
            .show(ui, |ui| {
                ui.strong(self.tr("network.col.candidate"));
                ui.strong(self.tr("network.col.type"));
                ui.strong(self.tr("network.col.score"));
                ui.strong(self.tr("network.col.signals"));
                ui.strong("");
                ui.end_row();
                for s in &rows {
                    let short: String = s.suggested_sample_guid.chars().take(12).collect();
                    ui.label(short).on_hover_text(&s.suggested_sample_guid);
                    ui.label(&s.suggestion_type);
                    ui.label(format!("{:.2}", s.score));
                    ui.label(s.signals.join(", "));
                    ui.horizontal(|ui| {
                        if ui.button(self.tr("network.introduce")).clicked() {
                            introduce = Some(s.clone());
                        }
                        if ui
                            .button(self.tr("matching.dismiss"))
                            .on_hover_text(self.tr("matching.dismissHint"))
                            .clicked()
                        {
                            dismiss = Some(s.suggested_sample_guid.clone());
                        }
                    });
                    ui.end_row();
                }
            });
        if let Some(suggestion) = introduce {
            self.status = self.tr("network.introducing").to_string();
            let _ = self.tx.send(Command::RequestIntroduction {
                suggestion,
                biosample_guid: self.matching_subject,
            });
        }
        if let Some(suggested_sample_guid) = dismiss {
            self.exchange_busy = true;
            let _ = self.tx.send(Command::DismissCandidate { suggested_sample_guid });
        }
    }

    /// Every conversation that has not produced a result yet, with the action it is waiting on.
    fn matching_requests(&mut self, ui: &mut egui::Ui) {
        let rows: Vec<navigator_app::MatchingEntry> = self
            .matching
            .iter()
            .filter(|e| e.status != navigator_app::MatchingStatus::Exchanged)
            .cloned()
            .collect();
        if rows.is_empty() {
            ui.add_space(6.0);
            ui.weak(self.tr("matching.noRequests"));
            return;
        }
        let mut consent_for: Option<navigator_app::MatchingEntry> = None;
        let mut run: Option<navigator_app::MatchingEntry> = None;
        let mut forget: Option<String> = None;
        egui::Grid::new("matching_requests")
            .striped(true)
            .num_columns(5)
            .spacing([14.0, 4.0])
            .show(ui, |ui| {
                ui.strong(self.tr("matching.col.who"));
                ui.strong(self.tr("matching.col.direction"));
                ui.strong(self.tr("matching.col.purpose"));
                ui.strong(self.tr("matching.col.status"));
                ui.strong("");
                ui.end_row();
                for e in &rows {
                    // Before mutual consent there is no partner identity to show — the broker is
                    // symmetric-blind by design, so the request URI is all either side has.
                    match &e.partner_did {
                        Some(did) => {
                            let short: String = did.chars().take(20).collect();
                            ui.label(short).on_hover_text(did);
                        }
                        None => {
                            let short: String = e.request_uri.chars().take(16).collect();
                            ui.label(egui::RichText::new(short).italics())
                                .on_hover_text(self.tr("matching.blindHint"));
                        }
                    }
                    ui.label(self.tr(direction_key(e.direction)));
                    ui.label(if e.purpose.is_empty() { "—" } else { &e.purpose });
                    let (text, color) = self.status_chip(e.status);
                    ui.colored_label(color, text).on_hover_text(match &e.last_error {
                        Some(err) => err.clone(),
                        None => self.tr(status_hint_key(e.status)).to_string(),
                    });
                    ui.horizontal(|ui| {
                        if e.status == navigator_app::MatchingStatus::AwaitingConsent {
                            if ui.button(self.tr("matching.review")).clicked() {
                                consent_for = Some(e.clone());
                            }
                        } else if matches!(
                            e.status,
                            navigator_app::MatchingStatus::Ready | navigator_app::MatchingStatus::Failed
                        ) && e.session_id.is_some()
                        {
                            let label = if e.status == navigator_app::MatchingStatus::Failed {
                                self.tr("matching.retry")
                            } else {
                                self.tr("exchange.run")
                            };
                            if ui.add_enabled(!self.exchange_busy, egui::Button::new(label)).clicked() {
                                run = Some(e.clone());
                            }
                        }
                        if ui
                            .button(self.tr("matching.forget"))
                            .on_hover_text(self.tr("matching.forgetHint"))
                            .clicked()
                        {
                            forget = Some(e.request_uri.clone());
                        }
                    });
                    ui.end_row();
                }
            });
        if let Some(e) = consent_for {
            self.consent_prompt = Some(e);
        }
        if let Some(e) = run {
            self.run_matching_exchange(&e);
        }
        if let Some(request_uri) = forget {
            let _ = self.tx.send(Command::ForgetMatchingRequest { request_uri });
        }
    }

    /// Start the encrypted exchange for a consent-ready conversation, using the chosen subject.
    fn run_matching_exchange(&mut self, e: &navigator_app::MatchingEntry) {
        let (Some(session_id), Some(partner_did)) = (e.session_id.clone(), e.partner_did.clone()) else {
            self.status = self.tr("matching.notReady").to_string();
            return;
        };
        // Prefer the subject already bound to this conversation; the picker only supplies one when
        // the conversation has never chosen.
        let Some(guid) = e.biosample_guid.or(self.matching_subject) else {
            self.status = self.tr("matching.noSubject").to_string();
            return;
        };
        self.exchange_busy = true;
        self.status = self.tr("exchange.running").to_string();
        let _ = self.tx.send(Command::RunIbdExchange {
            info: navigator_app::ExchangeSessionInfo {
                session_id,
                request_uri: e.request_uri.clone(),
                purpose: e.purpose.clone(),
                partner_did,
                partner_key_uri: None,
            },
            biosample_guid: guid,
        });
    }

    /// Completed comparisons, across every subject in the workspace.
    fn matching_results(&mut self, ui: &mut egui::Ui) {
        let rows: Vec<navigator_app::MatchingEntry> =
            self.matching.iter().filter(|e| e.result.is_some()).cloned().collect();
        if rows.is_empty() {
            ui.add_space(6.0);
            ui.weak(self.tr("exchange.noResults"));
            return;
        }
        let mut message_partner: Option<String> = None;
        egui::Grid::new("matching_results")
            .striped(true)
            .num_columns(6)
            .spacing([14.0, 4.0])
            .show(ui, |ui| {
                ui.strong(self.tr("exchange.col.partner"));
                ui.strong(self.tr("exchange.col.shared"));
                ui.strong(self.tr("exchange.col.relationship"));
                ui.strong(self.tr("exchange.col.agreed"));
                ui.strong(self.tr("matching.col.reported"));
                ui.strong("");
                ui.end_row();
                for e in &rows {
                    let Some(r) = e.result.as_ref() else { continue };
                    let short: String = r.partner_did.chars().take(20).collect();
                    ui.label(short).on_hover_text(&r.partner_did);
                    ui.label(format!("{:.1} cM · {} seg", r.total_shared_cm, r.segment_count));
                    ui.label(&r.relationship);
                    if r.agreed {
                        ui.colored_label(OK_GREEN, self.tr("exchange.agreedYes"));
                    } else {
                        ui.colored_label(WARN_RED, self.tr("exchange.agreedNo"));
                    }
                    // Whether the AppView has this match on the discovery graph. Not every result
                    // can be reported: a disputed summary, or a conversation with no AppView sample
                    // handles, is deliberately kept private.
                    if e.attested {
                        ui.colored_label(OK_GREEN, self.tr("matching.reportedYes"))
                            .on_hover_text(self.tr("matching.reportedHint"));
                    } else {
                        ui.weak(self.tr("matching.reportedNo"))
                            .on_hover_text(self.tr("matching.notReportedHint"));
                    }
                    if ui.button(self.tr("dm.message")).clicked() {
                        message_partner = Some(r.partner_did.clone());
                    }
                    ui.end_row();
                }
            });
        if let Some(partner_did) = message_partner {
            let _ = self.tx.send(Command::DmInitiate { partner_did });
            self.nav = Nav::Community;
            self.community_tab = CommunityTab::Messages;
            self.dm_loaded = false;
        }
    }

    /// Colour + label for a lifecycle status.
    fn status_chip(&self, s: navigator_app::MatchingStatus) -> (&'static str, egui::Color32) {
        use navigator_app::MatchingStatus as S;
        let neutral = egui::Color32::from_rgb(170, 150, 40);
        match s {
            S::Requested => (self.tr("matching.status.requested"), neutral),
            S::AwaitingConsent => (
                self.tr("matching.status.awaiting"),
                egui::Color32::from_rgb(90, 140, 200),
            ),
            S::Declined => (self.tr("matching.status.declined"), egui::Color32::GRAY),
            S::Ready => (self.tr("matching.status.ready"), OK_GREEN),
            S::Exchanged => (self.tr("matching.status.exchanged"), OK_GREEN),
            S::Failed => (self.tr("matching.status.failed"), WARN_RED),
        }
    }
}

/// Agreement / success green, matching the exchange card's existing verdict colour.
const OK_GREEN: egui::Color32 = egui::Color32::from_rgb(60, 160, 60);
/// Disagreement / failure red (softer than [`DANGER`], which is reserved for destructive buttons).
const WARN_RED: egui::Color32 = egui::Color32::from_rgb(200, 90, 90);

/// i18n key for a direction.
fn direction_key(d: navigator_app::MatchingDirection) -> &'static str {
    match d {
        navigator_app::MatchingDirection::Outbound => "matching.dir.outbound",
        navigator_app::MatchingDirection::Inbound => "matching.dir.inbound",
    }
}

/// i18n key for the tooltip explaining what a status is waiting on.
fn status_hint_key(s: navigator_app::MatchingStatus) -> &'static str {
    use navigator_app::MatchingStatus as S;
    match s {
        S::Requested => "matching.hint.requested",
        S::AwaitingConsent => "matching.hint.awaiting",
        S::Declined => "matching.hint.declined",
        S::Ready => "matching.hint.ready",
        S::Exchanged => "matching.hint.exchanged",
        S::Failed => "matching.hint.failed",
    }
}
