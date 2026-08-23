//! `impl NavigatorApp` methods extracted from `ui.rs` (the `ibd` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    /// IBD that works with a chip. Pick two sources, each a WGS alignment or an imported chip, and
    /// compare them over the multi-build IBD panel. This is the volume path for chip↔WGS and
    /// chip↔chip, it knows the build, and an asset backs it. It needs the `ibd_panel` asset, and a
    /// chip source needs its raw file (source_path).
    pub(crate) fn genotyping_section(&mut self, ui: &mut egui::Ui) {
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

    /// The identity-verification verdict. The compare path for each source, and the one at subject
    /// level, share it.
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

    /// Draw the current IBD comparison result, if there is one: a summary line and a segment table.
    /// The picker of each source, and the consensus comparison at subject level, share it.
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
        // A segment ideogram for each chromosome. It uses true chromosome lengths when the genome
        // regions are in memory.
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

    /// Subject-level IBD: compare the autosomal consensus of this subject against that of another
    /// subject. This is the pooled-genotype path, and it genotypes no single source. A match that
    /// is almost complete is the dedup and identity signal, which you read off the relationship.
    ///
    /// A *Change* reveal picks the comparison target: it shows the current choice, then a filter
    /// over a virtualized list. It is not a dropdown. A `ComboBox` builds one widget for each entry
    /// on every frame its popup is open, and this list is every other subject in the workspace. At
    /// 10k that stops each frame. The subject picker of the Matching tab has this shape for the
    /// same reason.
    pub(crate) fn consensus_ibd_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if !self.all_biosamples.iter().any(|b| b.guid != guid) {
            ui.label(egui::RichText::new(self.tr("hint.ibdNoOtherSubjects")).weak());
            return;
        }
        // A lookup, and not a copy of every subject. The old build allocated the whole roster on
        // each frame.
        let sel = self
            .ibd_other_subject
            .and_then(|g| self.find_subject(g).map(|b| b.donor_identifier.clone()))
            .unwrap_or_else(|| "—".to_string());
        let mut toggle_picker = false;
        ui.horizontal(|ui| {
            ui.label(self.tr("ibd.otherSubject"));
            ui.label(egui::RichText::new(sel).strong());
            let label = if self.ibd_other_picking {
                self.tr("common.cancel")
            } else {
                self.tr("common.change")
            };
            if ui.button(label).clicked() {
                toggle_picker = true;
            }
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
            // The same-individual check (duplicate detection) over the same pooled consensus, with
            // no panel.
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
        if toggle_picker {
            self.ibd_other_picking = !self.ibd_other_picking;
            self.ibd_other_filter.clear();
        }
        self.ibd_other_picker(ui, guid);
        ui.label(egui::RichText::new(self.tr("hint.ibdConsensus")).weak().small());
        self.render_identity(ui);
        self.render_ibd_result(ui);
    }

    /// The filter and the virtualized subject list that the *Change* button of
    /// [`Self::consensus_ibd_section`] reveals. It builds only the visible rows, so the cost does
    /// not depend on the size of the workspace. It assembles the filtered `Vec` from immutable
    /// reads first, so the scroll closure borrows only locals.
    fn ibd_other_picker(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if !self.ibd_other_picking {
            return;
        }
        ui.add_space(4.0);
        let hint = self.tr("subjects.filter");
        ui.add(
            egui::TextEdit::singleline(&mut self.ibd_other_filter)
                .hint_text(hint)
                .desired_width(280.0),
        );
        let needle = self.ibd_other_filter.trim().to_lowercase();
        let rows: Vec<(SampleGuid, String)> = self
            .all_biosamples
            .iter()
            .filter(|b| b.guid != guid)
            .filter(|b| needle.is_empty() || b.donor_identifier.to_lowercase().contains(&needle))
            .map(|b| (b.guid, b.donor_identifier.clone()))
            .collect();
        ui.label(egui::RichText::new(format!("{}", rows.len())).weak().small());
        if rows.is_empty() {
            ui.label(egui::RichText::new(self.tr("subjects.noMatch")).weak());
            return;
        }

        let selected = self.ibd_other_subject;
        let mut pick = None;
        let row_h = ui.spacing().interact_size.y;
        egui::ScrollArea::vertical()
            .id_salt("ibd_other_list")
            .max_height(180.0)
            .auto_shrink([false, false])
            .show_rows(ui, row_h, rows.len(), |ui, range| {
                for i in range {
                    let (g, name) = &rows[i];
                    if ui.selectable_label(selected == Some(*g), name).clicked() {
                        pick = Some(*g);
                    }
                }
            });
        if let Some(g) = pick {
            self.ibd_other_subject = Some(g);
            self.ibd_other_picking = false;
            self.ibd_other_filter.clear();
        }
    }

    /// The federated exchanges of this subject that are complete. Discovery and consent are **not**
    /// here. They are account-scoped, and they live in the top-level Matching tab. This card
    /// answers one question: "what did the network find for *this person*". It goes into the page
    /// scroll, with no nested ScrollArea.
    pub(crate) fn exchange_section(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        if self.account.is_none() {
            ui.label(self.tr("network.signInRequired"));
            return;
        }
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(self.tr("hint.encryptedExchange")).weak().small());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button(self.tr("matching.openTab")).clicked() {
                    self.nav = Nav::Matching;
                    self.matching_subject = Some(guid);
                    let _ = self.tx.send(Command::RefreshMatching);
                }
            });
        });

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
                        // Open an encrypted DM with this match (social 3a). It sends a DM request
                        // and goes to Community → Messages, where the conversation appears after
                        // the other side accepts.
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
            // An AI explanation of these matches on this tab (M5). It is additive, and it sits
            // below the structured table.
            ui.add_space(6.0);
            self.ai_explain(ui, guid, SignalKind::Ibd);
        }
    }

    /// An mtDNA haplogroup that comes directly from the chrM of the alignment. It is the standalone
    /// equivalent of "Assign Y haplogroup" in the Y-DNA section.
    pub(crate) fn mt_haplogroup_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self.alignment_has_bam(alignment_id);
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
        // The code resolves the reference from the build on demand, so it needs only the BAM.
        let has_bam = self.alignment_has_bam(alignment_id);

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
                self.cancel_button(ui);
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
