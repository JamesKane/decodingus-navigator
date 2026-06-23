//! `impl NavigatorApp` methods extracted from `ui.rs` (the `ibd` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    pub(crate) fn genotyping_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let Some(panel_id) = self.selected_panel else {
            ui.label(self.tr("panel.selectInSidebar"));
            return;
        };
        let panel_name = self
            .panels
            .iter()
            .find(|p| p.panel.id == panel_id)
            .map(|p| p.panel.name.clone())
            .unwrap_or_default();
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            ui.label(self.tr("form.ploidy"));
            ui.add(egui::TextEdit::singleline(&mut self.forms.ploidy).desired_width(32.0));
            if ui
                .add_enabled(
                    has_bam && !self.running_genotype,
                    egui::Button::new(format!("Genotype vs {panel_name}")),
                )
                .clicked()
            {
                self.running_genotype = true;
                let _ = self.tx.send(Command::GenotypePanel {
                    alignment_id,
                    panel_id,
                    ploidy: self.ploidy(),
                });
            }
            if self.running_genotype {
                ui.spinner();
            }
        });

        match &self.panel_genotypes {
            Some(genos) => {
                let (mut hr, mut het, mut ha, mut nc) = (0, 0, 0, 0);
                for g in genos {
                    match g.dosage {
                        0 => hr += 1,
                        1 => het += 1,
                        2 => ha += 1,
                        _ => nc += 1,
                    }
                }
                ui.label(format!(
                    "{} sites — {hr} hom-ref, {het} het, {ha} hom-alt, {nc} no-call",
                    genos.len()
                ));
            }
            None if !self.running_genotype => {
                ui.label(self.tr("panel.notGenotyped"));
            }
            None => {}
        }

        // IBD compare against another genotyped alignment.
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            ui.label(self.tr("ibd.vs"));
            let current = self
                .all_alignments
                .iter()
                .find(|a| Some(a.id) == self.ibd_other)
                .map(|a| format!("#{} {}", a.id, a.reference_build))
                .unwrap_or_else(|| "(pick alignment)".into());
            egui::ComboBox::from_id_salt("ibd_other")
                .selected_text(current)
                .show_ui(ui, |ui| {
                    for a in &self.all_alignments {
                        if a.id != alignment_id {
                            ui.selectable_value(
                                &mut self.ibd_other,
                                Some(a.id),
                                format!("#{} {}", a.id, a.reference_build),
                            );
                        }
                    }
                });
            let ready = self.ibd_other.is_some() && !self.running_ibd;
            if ui
                .add_enabled(ready, egui::Button::new(self.tr("action.compare")))
                .clicked()
            {
                self.running_ibd = true;
                self.ibd_result = None;
                self.identity = None;
                let _ = self.tx.send(Command::CompareIbd {
                    a: alignment_id,
                    b: self.ibd_other.unwrap(),
                    panel_id,
                    ploidy: self.ploidy(),
                });
            }
            if ui
                .add_enabled(self.ibd_other.is_some(), egui::Button::new(self.tr("ibd.verify")))
                .clicked()
            {
                self.identity = None;
                let _ = self.tx.send(Command::VerifyIdentity {
                    a: alignment_id,
                    b: self.ibd_other.unwrap(),
                    panel_id,
                    ploidy: self.ploidy(),
                });
            }
            if self.running_ibd {
                ui.spinner();
            }
        });

        // Chip-compatible IBD: pick two sources (each a WGS alignment or an imported chip) and
        // compare over the multi-build IBD panel — the chip↔WGS / chip↔chip volume path. Needs the
        // `ibd_panel` asset built; a chip source needs its raw file (source_path).
        let mut sources: Vec<(navigator_app::IbdSource, String)> = Vec::new();
        for a in &self.all_alignments {
            sources.push((
                navigator_app::IbdSource::Alignment(a.id),
                format!("WGS #{} {}", a.id, a.reference_build),
            ));
        }
        for c in self.chip_profiles.iter().filter(|c| c.source_path.is_some()) {
            sources.push((
                navigator_app::IbdSource::Chip(c.id),
                format!("{} chip #{}", c.provider, c.id),
            ));
        }
        let label_of = |src: Option<navigator_app::IbdSource>| -> String {
            src.and_then(|s| sources.iter().find(|(x, _)| *x == s).map(|(_, l)| l.clone()))
                .unwrap_or_else(|| "(pick source)".into())
        };
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            ui.label(self.tr("ibd.chipCompatible"));
            egui::ComboBox::from_id_salt("ibd_src_a")
                .selected_text(label_of(self.ibd_src_a))
                .show_ui(ui, |ui| {
                    for (s, l) in &sources {
                        ui.selectable_value(&mut self.ibd_src_a, Some(*s), l);
                    }
                });
            ui.label(self.tr("ibd.vs"));
            egui::ComboBox::from_id_salt("ibd_src_b")
                .selected_text(label_of(self.ibd_src_b))
                .show_ui(ui, |ui| {
                    for (s, l) in &sources {
                        ui.selectable_value(&mut self.ibd_src_b, Some(*s), l);
                    }
                });
            let ready = self.ibd_src_a.is_some()
                && self.ibd_src_b.is_some()
                && self.ibd_src_a != self.ibd_src_b
                && !self.running_ibd;
            if ui
                .add_enabled(ready, egui::Button::new(self.tr("action.compare")))
                .clicked()
            {
                self.running_ibd = true;
                self.ibd_result = None;
                self.identity = None;
                let _ = self.tx.send(Command::CompareIbdSources {
                    a: self.ibd_src_a.unwrap(),
                    b: self.ibd_src_b.unwrap(),
                });
            }
        });

        self.render_identity(ui);
        self.render_ibd_result(ui);
    }

    /// The identity-verification verdict (shared by the per-source + subject-level compare paths).
    fn render_identity(&self, ui: &mut egui::Ui) {
        let Some(v) = &self.identity else { return };
        let (txt, col) = match v.status {
            VerificationStatus::VerifiedSame => ("same individual", egui::Color32::from_rgb(60, 160, 60)),
            VerificationStatus::LikelySame => ("likely same", egui::Color32::from_rgb(120, 160, 60)),
            VerificationStatus::Uncertain => ("uncertain", egui::Color32::from_rgb(170, 150, 40)),
            VerificationStatus::LikelyDifferent => ("likely different", egui::Color32::from_rgb(200, 120, 40)),
            VerificationStatus::VerifiedDifferent => ("different individuals", egui::Color32::from_rgb(200, 60, 60)),
        };
        ui.horizontal(|ui| {
            ui.label(self.tr("ibd.identity"));
            ui.colored_label(col, txt);
            if let Some(c) = v.snp_concordance {
                ui.label(format!("SNP concordance {:.3} over {} sites", c, v.sites_compared));
            }
            if v.y_str_markers > 0 {
                ui.label(format!(
                    "· Y-STR {}/{} differ",
                    v.y_str_distance.unwrap_or(0),
                    v.y_str_markers
                ));
            }
        });
    }

    /// Render the current IBD comparison result (summary line + segment table), if any. Shared by the
    /// per-source picker and the subject-level consensus comparison.
    fn render_ibd_result(&mut self, ui: &mut egui::Ui) {
        // Clone out of the borrow so the export button can touch `self.status` / `self.tx` below.
        let Some(cmp) = self.ibd_result.clone() else { return };
        ui.label(format!(
            "{:?} — total {:.1} cM, {} segment(s), longest {:.1} cM  ·  {} overlapping sites",
            cmp.summary.relationship,
            cmp.summary.total_shared_cm,
            cmp.summary.segment_count,
            cmp.summary.longest_segment_cm,
            cmp.overlapping_sites,
        ));
        if cmp.segments.is_empty() {
            return;
        }
        // Per-chromosome segment ideogram (true chr lengths when genome regions are loaded).
        ui.add_space(6.0);
        ui.label(egui::RichText::new(self.tr("ibd.segmentMap")).strong().small());
        let regions = self.genome_regions.as_ref().map(|(_, r)| r.as_ref());
        draw_ibd_segments(ui, &cmp.segments, regions);

        ui.add_space(4.0);
        if ui.button(self.tr("ibd.exportSegments")).clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name("ibd_segments.tsv")
                .add_filter("TSV", &["tsv"])
                .save_file()
            {
                self.status = format!("Exporting {}…", self.tr("ibd.exportSegments"));
                let _ = self.tx.send(Command::ExportIbdSegments {
                    segments: cmp.segments.clone(),
                    path,
                });
            }
        }

        ui.add_space(4.0);
        ui.collapsing(self.tr("ibd.segmentTable"), |ui| {
            egui::Grid::new("ibd_segments")
                .striped(true)
                .num_columns(4)
                .show(ui, |ui| {
                    ui.strong(self.tr("table.chr"));
                    ui.strong(self.tr("table.start"));
                    ui.strong(self.tr("table.end"));
                    ui.strong(self.tr("table.cm"));
                    ui.end_row();
                    for s in &cmp.segments {
                        ui.label(&s.chromosome);
                        ui.label(s.start_position.to_string());
                        ui.label(s.end_position.to_string());
                        ui.label(format!("{:.1}", s.length_cm));
                        ui.end_row();
                    }
                });
        });
    }

    /// Subject-level IBD: compare this subject's autosomal consensus against another subject's — the
    /// pooled-genotype path (no per-source genotyping). A near-complete match is the dedup/identity
    /// signal (read off the relationship).
    pub(crate) fn consensus_ibd_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let others: Vec<(SampleGuid, String)> = self
            .all_biosamples
            .iter()
            .filter(|b| b.guid != guid)
            .map(|b| (b.guid, b.donor_identifier.clone()))
            .collect();
        if others.is_empty() {
            ui.label(egui::RichText::new(self.tr("hint.ibdNoOtherSubjects")).weak());
            return;
        }
        let sel = self
            .ibd_other_subject
            .and_then(|g| others.iter().find(|(og, _)| *og == g).map(|(_, l)| l.clone()))
            .unwrap_or_else(|| "—".to_string());
        ui.horizontal(|ui| {
            ui.label(self.tr("ibd.otherSubject"));
            egui::ComboBox::from_id_salt("ibd_subject")
                .selected_text(sel)
                .show_ui(ui, |ui| {
                    for (og, l) in &others {
                        ui.selectable_value(&mut self.ibd_other_subject, Some(*og), l);
                    }
                });
            let ready = self.ibd_other_subject.is_some() && !self.running_ibd;
            if ui
                .add_enabled(ready, egui::Button::new(self.tr("ibd.compare")))
                .clicked()
            {
                self.running_ibd = true;
                self.ibd_result = None;
                self.status = "Comparing consensuses…".into();
                let _ = self.tx.send(Command::CompareIbdConsensus {
                    a: guid,
                    b: self.ibd_other_subject.unwrap(),
                });
            }
            // Same-individual check (duplicate detection) over the same pooled consensus — no panel.
            if ui
                .add_enabled(ready, egui::Button::new(self.tr("ibd.verifyIdentity")))
                .clicked()
            {
                self.identity = None;
                self.status = "Verifying identity…".into();
                let _ = self.tx.send(Command::VerifyIdentityConsensus {
                    a: guid,
                    b: self.ibd_other_subject.unwrap(),
                });
            }
            if self.running_ibd {
                ui.spinner();
            }
        });
        ui.label(egui::RichText::new(self.tr("hint.ibdConsensus")).weak().small());
        self.render_identity(ui);
        self.render_ibd_result(ui);
    }

    /// Federated IBD: the AppView's pseudonymous "people who may share DNA with you" list,
    /// mined from the records we've published. Distinct from the local 1:1 compare above —
    /// these are network candidates we haven't exchanged any genotypes with. Requesting an
    /// introduction opens a PENDING request; the consent round-trip + actual segment detection
    /// are a later phase (the AppView's symmetric counterpart-discovery is still being speced).
    /// Simple-mode "Genetic relatives" card: the live network match suggestions, framed in plain
    /// language (strength + shared-signal summary, pseudonymous handles), with a privacy note and a
    /// per-match Connect action. Sign-in gated. Kept separate from the precomputed brief because
    /// matches are live/online; the brief itself stays local.
    pub(crate) fn simple_relatives_section(&mut self, ui: &mut egui::Ui) {
        let title = self.tr("brief.relatives");
        if self.account.is_none() {
            card(ui, title, |ui| {
                ui.label(self.tr("brief.relativesSignIn"));
            });
            return;
        }

        // Pre-build rows (owned) so the card body can mutate self (button) without a borrow clash.
        let rows: Vec<(String, String, String, Option<String>)> = self
            .ibd_suggestions
            .iter()
            .map(|s| {
                let handle: String = s.suggested_sample_guid.chars().take(12).collect();
                let strength = self.tr(match s.score {
                    x if x >= 0.8 => "brief.matchStrong",
                    x if x >= 0.5 => "brief.matchLikely",
                    _ => "brief.matchPossible",
                });
                let why = if s.signals.is_empty() {
                    String::new()
                } else {
                    format!("{} {}", self.tr("brief.relativesWhy"), s.signals.join(", "))
                };
                (
                    handle,
                    strength.to_string(),
                    why,
                    self.ibd_intros.get(&s.suggested_sample_guid).cloned(),
                )
            })
            .collect();
        let loading = self.loading_ibd_suggestions;
        let guids: Vec<String> = self
            .ibd_suggestions
            .iter()
            .map(|s| s.suggested_sample_guid.clone())
            .collect();

        let mut do_find = false;
        let mut introduce: Option<String> = None;
        card(ui, title, |ui| {
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(!loading, egui::Button::new(self.tr("brief.relativesFind")))
                    .clicked()
                {
                    do_find = true;
                }
                if loading {
                    ui.spinner();
                }
            });
            ui.label(egui::RichText::new(self.tr("brief.relativesNote")).weak().small());
            ui.add_space(6.0);

            if rows.is_empty() {
                if !loading {
                    ui.label(egui::RichText::new(self.tr("brief.relativesEmpty")).weak());
                }
                return;
            }
            for (i, (handle, strength, why, intro)) in rows.iter().enumerate() {
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.horizontal(|ui| {
                        chip(
                            ui,
                            strength,
                            egui::Color32::from_rgb(40, 52, 70),
                            egui::Color32::from_rgb(150, 190, 240),
                        );
                        ui.label(egui::RichText::new(handle).monospace())
                            .on_hover_text(&guids[i]);
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| match intro {
                            Some(status) => {
                                ui.label(egui::RichText::new(status).weak().small());
                            }
                            None => {
                                if ui.button(self.tr("brief.relativesConnect")).clicked() {
                                    introduce = Some(guids[i].clone());
                                }
                            }
                        });
                    });
                    if !why.is_empty() {
                        ui.label(egui::RichText::new(why).weak().small());
                    }
                });
            }
        });

        if do_find {
            self.loading_ibd_suggestions = true;
            self.status = self.tr("network.finding").to_string();
            let _ = self.tx.send(Command::LoadIbdSuggestions);
        }
        if let Some(guid) = introduce {
            self.status = self.tr("network.introducing").to_string();
            let _ = self.tx.send(Command::IbdIntroduce {
                suggested_sample_guid: guid,
            });
        }
    }

    pub(crate) fn network_suggestions_section(&mut self, ui: &mut egui::Ui) {
        if self.account.is_none() {
            ui.label(self.tr("network.signInRequired"));
            return;
        }
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
        ui.label(self.tr("network.note"));

        if self.ibd_suggestions.is_empty() {
            if !self.loading_ibd_suggestions {
                ui.add_space(4.0);
                ui.weak(self.tr("network.empty"));
            }
            return;
        }

        ui.add_space(6.0);
        // Collect the rows first so the table closure doesn't borrow `self` immutably while we
        // also need `self.tx` / `self.ibd_intros` (and to send commands without a borrow clash).
        let rows: Vec<(String, String, f64, String, Option<String>)> = self
            .ibd_suggestions
            .iter()
            .map(|s| {
                let signals = s.signals.join(", ");
                (
                    s.suggested_sample_guid.clone(),
                    s.suggestion_type.clone(),
                    s.score,
                    signals,
                    self.ibd_intros.get(&s.suggested_sample_guid).cloned(),
                )
            })
            .collect();

        let mut introduce: Option<String> = None;
        egui::Grid::new("ibd_suggestions")
            .striped(true)
            .num_columns(5)
            .show(ui, |ui| {
                ui.strong(self.tr("network.col.candidate"));
                ui.strong(self.tr("network.col.type"));
                ui.strong(self.tr("network.col.score"));
                ui.strong(self.tr("network.col.signals"));
                ui.strong("");
                ui.end_row();
                for (guid, ty, score, signals, intro) in &rows {
                    // Pseudonymous guid, shown truncated (it's an opaque AppView handle, not PII).
                    let short: String = guid.chars().take(12).collect();
                    ui.label(short).on_hover_text(guid);
                    ui.label(ty);
                    ui.label(format!("{score:.2}"));
                    ui.label(signals);
                    if let Some(status) = intro {
                        ui.label(status);
                    } else if ui.button(self.tr("network.introduce")).clicked() {
                        introduce = Some(guid.clone());
                    }
                    ui.end_row();
                }
            });
        if let Some(guid) = introduce {
            self.status = self.tr("network.introducing").to_string();
            let _ = self.tx.send(Command::IbdIntroduce {
                suggested_sample_guid: guid,
            });
        }
    }

    /// The encrypted edge-to-edge exchange (gap §4): inbound requests awaiting consent, consent-ready
    /// sessions to run an IBD exchange over, and this subject's saved results. Requires an active
    /// account (real PDS or did:key). Flows into the page scroll (no nested ScrollArea).
    pub(crate) fn exchange_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.account.is_none() {
            ui.label(self.tr("network.signInRequired"));
            return;
        }
        ui.horizontal(|ui| {
            if ui
                .add_enabled(!self.exchange_busy, egui::Button::new(self.tr("exchange.refresh")))
                .clicked()
            {
                self.exchange_busy = true;
                let _ = self.tx.send(Command::ExchangeInbox);
            }
            if self.exchange_busy {
                ui.spinner();
            }
            ui.label(egui::RichText::new(self.tr("hint.encryptedExchange")).weak().small());
        });

        // Inbound requests awaiting our consent (symmetric-blind: no initiator DID until consent).
        let incoming: Vec<(String, String, String)> = self
            .exchange_incoming
            .iter()
            .map(|r| (r.request_uri.clone(), r.purpose.clone(), r.created_at.clone()))
            .collect();
        if !incoming.is_empty() {
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("exchange.incoming")).strong());
            let mut consent: Option<(String, bool)> = None;
            egui::Grid::new(("exchange_incoming", guid))
                .striped(true)
                .num_columns(4)
                .spacing([12.0, 2.0])
                .show(ui, |ui| {
                    for (req, purpose, created) in &incoming {
                        let short: String = req.chars().take(16).collect();
                        ui.label(short).on_hover_text(req);
                        ui.label(purpose);
                        ui.label(egui::RichText::new(created).weak().small());
                        ui.horizontal(|ui| {
                            if ui.button(self.tr("exchange.accept")).clicked() {
                                consent = Some((req.clone(), true));
                            }
                            if ui.button(self.tr("exchange.decline")).clicked() {
                                consent = Some((req.clone(), false));
                            }
                        });
                        ui.end_row();
                    }
                });
            if let Some((request_uri, given)) = consent {
                self.exchange_busy = true;
                let _ = self.tx.send(Command::ExchangeConsent { request_uri, given });
            }
        }

        // Consent-ready sessions → run the IBD exchange (handshake + dosage exchange + attestation).
        let ready: Vec<navigator_app::ExchangeSessionInfo> = self.exchange_ready.clone();
        if !ready.is_empty() {
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("exchange.ready")).strong());
            let mut run: Option<navigator_app::ExchangeSessionInfo> = None;
            egui::Grid::new(("exchange_ready", guid))
                .striped(true)
                .num_columns(3)
                .spacing([12.0, 2.0])
                .show(ui, |ui| {
                    for info in &ready {
                        let short: String = info.partner_did.chars().take(20).collect();
                        ui.label(short).on_hover_text(&info.partner_did);
                        ui.label(&info.purpose);
                        if ui
                            .add_enabled(!self.exchange_busy, egui::Button::new(self.tr("exchange.run")))
                            .clicked()
                        {
                            run = Some(info.clone());
                        }
                        ui.end_row();
                    }
                });
            if let Some(info) = run {
                self.exchange_busy = true;
                self.status = self.tr("exchange.running").to_string();
                let _ = self.tx.send(Command::RunIbdExchange {
                    info,
                    biosample_guid: guid,
                });
            }
        }

        // This subject's saved results.
        ui.add_space(6.0);
        if self.exchange_results.is_empty() {
            ui.weak(self.tr("exchange.noResults"));
        } else {
            ui.label(egui::RichText::new(self.tr("exchange.results")).strong());
            let mut message_partner: Option<String> = None;
            egui::Grid::new(("exchange_results", guid))
                .striped(true)
                .num_columns(5)
                .spacing([14.0, 2.0])
                .show(ui, |ui| {
                    ui.strong(self.tr("exchange.col.partner"));
                    ui.strong(self.tr("exchange.col.shared"));
                    ui.strong(self.tr("exchange.col.relationship"));
                    ui.strong(self.tr("exchange.col.agreed"));
                    ui.end_row();
                    for r in &self.exchange_results {
                        let short: String = r.partner_did.chars().take(20).collect();
                        ui.label(short).on_hover_text(&r.partner_did);
                        ui.label(format!("{:.1} cM · {} seg", r.total_shared_cm, r.segment_count));
                        ui.label(&r.relationship);
                        if r.agreed {
                            ui.colored_label(egui::Color32::from_rgb(60, 160, 60), self.tr("exchange.agreedYes"));
                        } else {
                            ui.colored_label(egui::Color32::from_rgb(200, 90, 90), self.tr("exchange.agreedNo"));
                        }
                        // Open an encrypted DM with this match (social 3a) — sends a DM request and
                        // jumps to Community → Messages, where the conversation appears once accepted.
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
                self.dm_loaded = false; // force a fresh inbox/conversation load on entry
            }
            // Per-tab AI explanation of these matches (M5) — additive, below the structured table.
            ui.add_space(6.0);
            self.ai_explain(ui, guid, SignalKind::Ibd);
        }
    }

    /// mtDNA haplogroup assigned directly from the alignment's chrM — the standalone counterpart
    /// to the Y-DNA section's "Assign Y haplogroup".
    pub(crate) fn mt_haplogroup_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);
        ui.horizontal(|ui| {
            if ui
                .add_enabled(has_bam, egui::Button::new(self.tr("btn.assignMt")))
                .clicked()
            {
                self.status = "Assigning mtDNA haplogroup (fetching FTDNA mt tree)…".into();
                let _ = self
                    .tx
                    .send(Command::AssignMtdnaHaplogroupFromAlignment { alignment_id });
            }
            if !has_bam {
                ui.label(egui::RichText::new("(no BAM/CRAM path recorded)").weak());
            }
        });
        if let Some((id, assignment)) = &self.mt_haplogroup {
            if *id == alignment_id {
                show_assignment(ui, assignment);
            }
        }
    }

    /// De-novo haploid SNP calls for a specific `contig` (chrY on the Y-DNA tab, chrM on mtDNA).
    pub(crate) fn denovo_section(&mut self, ui: &mut egui::Ui, alignment_id: i64, contig: &str) {
        // Reference is resolved from the build on demand, so only the BAM is required.
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            let ready = has_bam && !self.running_denovo;
            let label = format!("{} ({contig})", self.tr("btn.runDenovo"));
            if ui.add_enabled(ready, egui::Button::new(label)).clicked() {
                self.running_denovo = true;
                self.denovo.remove(contig);
                self.status = format!("Calling {contig} on alignment #{alignment_id}…");
                let _ = self.tx.send(Command::RunDenovo {
                    alignment_id,
                    contig: contig.to_string(),
                });
            }
            if self.running_denovo {
                ui.spinner();
            }
            if !has_bam {
                ui.label(egui::RichText::new("(no BAM/CRAM recorded)").weak());
            }
        });

        match self.denovo.get(contig) {
            None if !self.running_denovo => {
                ui.label(egui::RichText::new("No calls yet — run for this contig.").weak());
            }
            None => {}
            Some(calls) if calls.is_empty() => {
                ui.label(self.tr("denovo.noCalls"));
            }
            Some(calls) => {
                ui.label(format!("{} SNP call(s)", calls.len()));
                egui::Grid::new(("denovo_calls", contig))
                    .striped(true)
                    .num_columns(4)
                    .show(ui, |ui| {
                        ui.strong(self.tr("table.position"));
                        ui.strong(self.tr("table.change"));
                        ui.strong(self.tr("table.depth"));
                        ui.strong(self.tr("table.af"));
                        ui.end_row();
                        for c in calls {
                            ui.label(c.position.to_string());
                            ui.label(format!("{}>{}", c.reference_allele, c.alternate_allele));
                            ui.label(c.depth.to_string());
                            ui.label(format!("{:.2}", c.allele_fraction));
                            ui.end_row();
                        }
                    });
            }
        }

        if self.denovo.get(contig).map(|c| !c.is_empty()).unwrap_or(false) {
            self.publish_row(
                ui,
                "Publish variants to PDS",
                Command::PublishVariants {
                    alignment_id,
                    contig: contig.to_string(),
                },
            );
        }
    }
}
