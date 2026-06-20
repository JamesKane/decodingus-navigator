//! `impl NavigatorApp` — the Community tab: the signed-in tester's social surface over the AppView's
//! signed Edge API (Support threads to the team / community Feed / Notifications). Account-global
//! (not per-subject); requires a signed-in identity. Mirrors the `central` rendering idioms.
use super::*;

impl NavigatorApp {
    /// The Community work area: a sub-tab bar over Support / Feed / Notifications. Gated on sign-in
    /// (the API is device-key-signed). Lazily loads all three sections on first entry this session.
    pub(crate) fn community_central(&mut self, ui: &mut egui::Ui) {
        if self.account.is_none() {
            empty_state(
                ui,
                self.tr("community.signedout.title"),
                self.tr("community.signedout.hint"),
            );
            return;
        }
        if !self.community_loaded {
            self.community_loaded = true;
            self.refresh_community();
        }
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.heading(self.tr("nav.community"));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button(self.tr("common.refresh")).clicked() {
                    self.refresh_community();
                }
            });
        });
        ui.separator();
        ui.add_space(4.0);
        self.community_tab = self.sub_bar(ui, self.community_tab, &CommunityTab::ALL);
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.add_space(4.0);
            match self.community_tab {
                CommunityTab::Support => self.community_support(ui),
                CommunityTab::Feed => self.community_feed_view(ui),
                CommunityTab::Messages => self.community_messages(ui),
                CommunityTab::Notifications => self.community_notifications(ui),
            }
        });
    }

    /// Re-poll all sections (also drives the app-bar bell via the Notifications event).
    fn refresh_community(&self) {
        let _ = self.tx.send(Command::LoadSupportThreads);
        let _ = self.tx.send(Command::LoadCommunityFeed);
        let _ = self.tx.send(Command::LoadNotifications);
        let _ = self.tx.send(Command::LoadRecruitmentInvitations);
    }

    // ---- Support: team↔tester threads --------------------------------------
    fn community_support(&mut self, ui: &mut egui::Ui) {
        // Reading one thread: back button + transcript + reply box.
        if let Some((cid, messages)) = self.open_thread.clone() {
            if ui.button(self.tr("community.back")).clicked() {
                self.open_thread = None;
                self.thread_reply.clear();
                return;
            }
            ui.add_space(6.0);
            for m in &messages {
                let team = self.tr("community.team");
                let who = if m.from_team {
                    team
                } else {
                    m.author.as_deref().unwrap_or("you")
                };
                let color = if m.from_team { ACCENT } else { egui::Color32::GRAY };
                ui.horizontal(|ui| {
                    ui.colored_label(color, egui::RichText::new(who).strong());
                    if let Some(at) = &m.at {
                        ui.label(egui::RichText::new(at).weak().small());
                    }
                });
                ui.label(&m.body);
                ui.add_space(6.0);
            }
            ui.separator();
            let reply_label = self.tr("community.reply");
            ui.add(
                egui::TextEdit::multiline(&mut self.thread_reply)
                    .hint_text("reply…")
                    .desired_rows(2)
                    .desired_width(f32::INFINITY),
            );
            let ready = !self.thread_reply.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new(reply_label)).clicked() {
                let _ = self.tx.send(Command::ReplySupportThread {
                    conversation_id: cid.clone(),
                    body: self.thread_reply.trim().to_string(),
                });
                self.thread_reply.clear();
            }
            return;
        }

        // New-thread composer.
        ui.group(|ui| {
            ui.label(egui::RichText::new(self.tr("community.newThread")).strong());
            ui.add(egui::TextEdit::singleline(&mut self.new_thread_subject).hint_text("subject (optional)"));
            ui.add(
                egui::TextEdit::multiline(&mut self.new_thread_body)
                    .hint_text("message to the team")
                    .desired_rows(3)
                    .desired_width(f32::INFINITY),
            );
            let send_label = self.tr("community.send");
            let ready = !self.new_thread_body.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new(send_label)).clicked() {
                let _ = self.tx.send(Command::OpenSupportThread {
                    subject: self.new_thread_subject.trim().to_string(),
                    body: self.new_thread_body.trim().to_string(),
                });
                self.new_thread_subject.clear();
                self.new_thread_body.clear();
            }
        });
        ui.add_space(8.0);

        if self.support_threads.is_empty() {
            ui.label(egui::RichText::new(self.tr("community.noThreads")).weak());
            return;
        }
        for t in self.support_threads.clone() {
            let dot = if t.unread { "● " } else { "" };
            let subject = t.subject.as_deref().unwrap_or("(no subject)");
            let label = format!("{dot}{subject}   [{}]", t.status);
            if ui.selectable_label(false, label).clicked() {
                let _ = self.tx.send(Command::LoadSupportThread {
                    conversation_id: t.conversation_id.clone(),
                });
            }
        }
    }

    // ---- Feed: announcements + community + federated ------------------------
    fn community_feed_view(&mut self, ui: &mut egui::Ui) {
        // Composer.
        ui.group(|ui| {
            ui.label(egui::RichText::new(self.tr("community.postToFeed")).strong());
            ui.add(
                egui::TextEdit::multiline(&mut self.feed_content)
                    .hint_text("share something with the community")
                    .desired_rows(2)
                    .desired_width(f32::INFINITY),
            );
            ui.horizontal(|ui| {
                ui.add(
                    egui::TextEdit::singleline(&mut self.feed_topic)
                        .hint_text("topic (optional, e.g. haplogroup:R-M269)"),
                );
                let post_label = self.tr("community.post");
                let ready = !self.feed_content.trim().is_empty();
                if ui.add_enabled(ready, egui::Button::new(post_label)).clicked() {
                    let topic = self.feed_topic.trim();
                    let _ = self.tx.send(Command::PostCommunity {
                        content: self.feed_content.trim().to_string(),
                        topic: (!topic.is_empty()).then(|| topic.to_string()),
                        publish_pds: self.feed_publish_pds,
                    });
                    self.feed_content.clear();
                    self.feed_topic.clear();
                }
            });
            // Opt-in federation: publishing to your own PDS makes the post a portable, public
            // `feed.post` record (mirrored back as a "via Atmosphere" entry). Only a real PDS
            // account has a repo to write to — a local did:key identity can't federate.
            let can_federate = self.account.as_deref().is_some_and(|d| !d.starts_with("did:key:"));
            if can_federate {
                // Bind the labels first: `tr()` borrows `&self`, which would clash with the
                // `&mut self.feed_publish_pds` checkbox binding (the i18n borrow gotcha).
                let label = self.tr("community.publishPds");
                let hint = self.tr("community.publishPds.hint");
                ui.checkbox(&mut self.feed_publish_pds, label).on_hover_text(hint);
            }
        });
        ui.add_space(8.0);

        let Some(feed) = self.feed.clone() else {
            ui.label(egui::RichText::new(self.tr("community.loading")).weak());
            return;
        };
        if feed.announcements.is_empty() && feed.community.is_empty() && feed.federated.is_empty() {
            ui.label(egui::RichText::new(self.tr("community.emptyFeed")).weak());
            return;
        }
        for a in &feed.announcements {
            feed_card(
                ui,
                FeedCard {
                    author: a.author.as_deref(),
                    topic: a.topic.as_deref(),
                    body: &a.content,
                    at: a.at.as_deref(),
                    pinned: a.pinned,
                    badge: None,
                },
            );
        }
        for c in &feed.community {
            feed_card(
                ui,
                FeedCard {
                    author: c.author.as_deref(),
                    topic: c.topic.as_deref(),
                    body: &c.content,
                    at: c.at.as_deref(),
                    pinned: false,
                    badge: None,
                },
            );
        }
        let via = self.tr("community.viaAtmosphere");
        for f in &feed.federated {
            feed_card(
                ui,
                FeedCard {
                    author: f.author.as_deref(),
                    topic: f.topic.as_deref(),
                    body: &f.text,
                    at: f.at.as_deref(),
                    pinned: false,
                    badge: Some(via),
                },
            );
        }
    }

    // ---- Messages: peer DMs over the encrypted D1 relay (social 3a) ---------
    fn community_messages(&mut self, ui: &mut egui::Ui) {
        if !self.dm_loaded {
            self.dm_loaded = true;
            let _ = self.tx.send(Command::LoadDmInbox);
            let _ = self.tx.send(Command::LoadDmConversations);
        }

        // Open conversation: transcript + composer.
        if let Some((session_id, messages)) = self.open_dm.clone() {
            ui.horizontal(|ui| {
                if ui.button(self.tr("community.back")).clicked() {
                    self.open_dm = None;
                    self.dm_compose.clear();
                }
                if ui.button(self.tr("common.refresh")).clicked() {
                    let _ = self.tx.send(Command::DmSync {
                        session_id: session_id.clone(),
                    });
                }
            });
            if self.open_dm.is_none() {
                return; // backed out this frame
            }
            ui.add_space(6.0);
            let me = self.account.clone().unwrap_or_default();
            for m in &messages {
                let mine = m.from_did == me;
                let who = if mine {
                    self.tr("dm.you").to_string()
                } else {
                    short_did(&m.from_did)
                };
                let color = if mine { egui::Color32::GRAY } else { ACCENT };
                ui.horizontal(|ui| {
                    ui.colored_label(color, egui::RichText::new(who).strong());
                    ui.label(egui::RichText::new(&m.created_at).weak().small());
                });
                ui.label(&m.body);
                ui.add_space(6.0);
            }
            ui.separator();
            let send_label = self.tr("dm.send");
            ui.add(
                egui::TextEdit::multiline(&mut self.dm_compose)
                    .hint_text("message…")
                    .desired_rows(2)
                    .desired_width(f32::INFINITY),
            );
            let ready = !self.dm_compose.trim().is_empty();
            if ui.add_enabled(ready, egui::Button::new(send_label)).clicked() {
                let _ = self.tx.send(Command::DmSend {
                    session_id,
                    text: self.dm_compose.trim().to_string(),
                });
            }
            return;
        }

        // Start a DM by partner DID.
        ui.group(|ui| {
            ui.label(egui::RichText::new(self.tr("dm.startTitle")).strong());
            ui.add(egui::TextEdit::singleline(&mut self.dm_partner_did).hint_text("did:plc:…"));
            let start_label = self.tr("dm.start");
            let ready = self.dm_partner_did.trim().starts_with("did:");
            if ui.add_enabled(ready, egui::Button::new(start_label)).clicked() {
                let _ = self.tx.send(Command::DmInitiate {
                    partner_did: self.dm_partner_did.trim().to_string(),
                });
            }
        });
        ui.add_space(8.0);

        // Inbound requests awaiting our consent (symmetric-blind).
        if !self.dm_incoming.is_empty() {
            ui.label(egui::RichText::new(self.tr("dm.incoming")).strong());
            for r in self.dm_incoming.clone() {
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new(&r.created_at).weak().small());
                    if ui.button(self.tr("dm.accept")).clicked() {
                        let _ = self.tx.send(Command::DmConsent {
                            request_uri: r.request_uri.clone(),
                            given: true,
                        });
                    }
                    if ui.button(self.tr("dm.decline")).clicked() {
                        let _ = self.tx.send(Command::DmConsent {
                            request_uri: r.request_uri.clone(),
                            given: false,
                        });
                    }
                });
            }
            ui.add_space(8.0);
        }

        // Consent-ready sessions to connect (one-time handshake; both peers must be online).
        if !self.dm_ready.is_empty() {
            ui.label(egui::RichText::new(self.tr("dm.ready")).strong());
            for info in self.dm_ready.clone() {
                ui.horizontal(|ui| {
                    ui.label(short_did(&info.partner_did));
                    if ui.button(self.tr("dm.connect")).clicked() {
                        let _ = self.tx.send(Command::DmConnect { info: info.clone() });
                    }
                });
            }
            ui.add_space(8.0);
        }

        // Conversation list.
        ui.label(egui::RichText::new(self.tr("dm.conversations")).strong());
        if self.dm_conversations.is_empty() {
            ui.label(egui::RichText::new(self.tr("dm.empty")).weak());
            return;
        }
        for c in self.dm_conversations.clone() {
            let dot = if c.unread > 0 { "● " } else { "" };
            let preview = c.last_body.as_deref().unwrap_or("");
            let label = format!("{dot}{}  —  {}", short_did(&c.partner_did), preview);
            if ui.selectable_label(false, label).clicked() {
                let _ = self.tx.send(Command::LoadDmMessages {
                    session_id: c.session_id.clone(),
                });
            }
        }
    }

    // ---- Recruitment invitations (3c, respond-only) ------------------------
    /// Open recruitment invitations the user can accept/decline (they also arrive as SYSTEM
    /// notifications below; this is the actionable view). Hidden when there are none.
    fn recruitment_invitations_section(&mut self, ui: &mut egui::Ui) {
        if self.recruitment_invitations.is_empty() {
            return;
        }
        ui.label(egui::RichText::new(self.tr("recruit.title")).strong());
        let accept_label = self.tr("recruit.accept");
        let decline_label = self.tr("recruit.decline");
        for inv in self.recruitment_invitations.clone() {
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.colored_label(ACCENT, egui::RichText::new(&inv.title).strong());
                    ui.label(egui::RichText::new(&inv.project_name).weak().small());
                });
                ui.label(&inv.message);
                ui.horizontal(|ui| {
                    if ui.button(accept_label).clicked() {
                        let _ = self.tx.send(Command::RespondRecruitment {
                            campaign_id: inv.campaign_id,
                            accept: true,
                        });
                    }
                    if ui.button(decline_label).clicked() {
                        let _ = self.tx.send(Command::RespondRecruitment {
                            campaign_id: inv.campaign_id,
                            accept: false,
                        });
                    }
                });
            });
            ui.add_space(4.0);
        }
        ui.separator();
        ui.add_space(4.0);
    }

    // ---- Notifications -----------------------------------------------------
    fn community_notifications(&mut self, ui: &mut egui::Ui) {
        self.recruitment_invitations_section(ui);
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(format!("{}: {}", self.tr("community.unread"), self.notif_unread)).weak());
            if ui.button(self.tr("community.markAllRead")).clicked() {
                let _ = self.tx.send(Command::MarkNotificationRead { id: None });
            }
        });
        ui.separator();
        if self.notifications.is_empty() {
            ui.label(egui::RichText::new(self.tr("community.noNotifications")).weak());
            return;
        }
        for n in self.notifications.clone() {
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    if n.unread {
                        ui.colored_label(ACCENT, "●");
                    }
                    ui.label(egui::RichText::new(&n.title).strong());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if n.unread && ui.button(self.tr("community.markRead")).clicked() {
                            let _ = self.tx.send(Command::MarkNotificationRead { id: Some(n.id.clone()) });
                        }
                    });
                });
                if let Some(b) = &n.body {
                    ui.label(egui::RichText::new(b).weak().small());
                }
                ui.horizontal(|ui| {
                    if let Some(actor) = &n.actor {
                        ui.label(egui::RichText::new(actor).weak().small());
                    }
                    if let Some(at) = &n.at {
                        ui.label(egui::RichText::new(at).weak().small());
                    }
                });
            });
            ui.add_space(4.0);
        }
    }
}

/// A truncated DID for display (pseudonymous handle, not PII) — full value goes in a hover.
fn short_did(did: &str) -> String {
    did.chars().take(20).collect()
}

/// A feed entry's display fields (announcement / community / federated).
struct FeedCard<'a> {
    author: Option<&'a str>,
    topic: Option<&'a str>,
    body: &'a str,
    at: Option<&'a str>,
    pinned: bool,
    /// Provenance badge (e.g. "via Atmosphere" for a federated post).
    badge: Option<&'a str>,
}

/// Render one feed entry card.
fn feed_card(ui: &mut egui::Ui, card: FeedCard) {
    ui.group(|ui| {
        ui.horizontal(|ui| {
            if card.pinned {
                ui.label("📌");
            }
            ui.colored_label(ACCENT, egui::RichText::new(card.author.unwrap_or("?")).strong());
            if let Some(t) = card.topic {
                ui.label(egui::RichText::new(format!("#{t}")).weak().small());
            }
            if let Some(b) = card.badge {
                ui.label(egui::RichText::new(b).weak().small().italics());
            }
            if let Some(at) = card.at {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(egui::RichText::new(at).weak().small());
                });
            }
        });
        ui.label(card.body);
    });
    ui.add_space(4.0);
}
