//! `impl NavigatorApp` methods extracted from `ui.rs` (the `sources` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    pub(crate) fn coverage_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_paths = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some() && a.reference_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_paths && !self.running, egui::Button::new(self.tr("btn.runCoverage"))).clicked() {
                self.running = true;
                self.status = format!("Running coverage on alignment #{alignment_id}…");
                let _ = self.tx.send(Command::RunCoverage(alignment_id));
            }
            if self.running {
                ui.spinner();
            }
            if !has_paths {
                ui.label(self.tr("hint.noBamRefPath"));
            }
        });

        // Local copy so the table's row clicks can update the histogram selection without
        // borrowing `self` mutably while `&self.coverage` is held; written back after the match.
        let mut sel = self.coverage_hist_contig;
        match &self.coverage {
            None if !self.running => {
                ui.label(self.tr("coverage.none"));
            }
            None => {}
            Some(c) => {
                egui::Grid::new("coverage_metrics").striped(true).num_columns(2).show(ui, |ui| {
                    let row = |ui: &mut egui::Ui, k: &str, v: String| {
                        ui.label(k);
                        ui.label(v);
                        ui.end_row();
                    };
                    row(ui, "Genome territory", c.genome_territory.to_string());
                    row(ui, "Mean coverage", format!("{:.2}", c.mean_coverage));
                    row(ui, "Median coverage", format!("{:.0}", c.median_coverage));
                    row(ui, "MAD coverage", format!("{:.0}", c.mad_coverage));
                    row(ui, "Callable bases", c.callable_bases.to_string());
                    row(ui, "% ≥10x", format!("{:.1}%", c.pct_10x * 100.0));
                    row(ui, "% ≥20x", format!("{:.1}%", c.pct_20x * 100.0));
                    row(ui, "% ≥30x", format!("{:.1}%", c.pct_30x * 100.0));
                    row(ui, "% excl. low MAPQ", format!("{:.1}%", c.pct_exc_mapq * 100.0));
                    row(ui, "% excl. low base-Q", format!("{:.1}%", c.pct_exc_baseq * 100.0));
                });

                // Drop a stale selection (e.g. fewer contigs than a prior result).
                if let Some(i) = sel {
                    if i >= c.contig_coverage_stats.len() {
                        sel = None;
                    }
                }

                ui.separator();
                ui.strong("Per-contig coverage");
                ui.horizontal(|ui| {
                    ui.label("Histogram:");
                    if ui.selectable_label(sel.is_none(), "Whole genome").clicked() {
                        sel = None;
                    }
                    ui.weak("or click a contig row");
                });

                // Per-contig table: stats joined with the GATK/callable breakdown by header order.
                egui::ScrollArea::vertical().max_height(240.0).id_salt("cov_contig_table").show(ui, |ui| {
                    egui::Grid::new("coverage_contig_grid").striped(true).num_columns(12).show(ui, |ui| {
                        for h in [
                            "Contig", "Length", "Reads", "Mean depth", "Cov %", "Callable",
                            "NoCov", "LowCov", "ExcessCov", "PoorMQ", "Mean BQ", "Mean MQ",
                        ] {
                            ui.strong(h);
                        }
                        ui.end_row();

                        for (i, s) in c.contig_coverage_stats.iter().enumerate() {
                            if ui.selectable_label(sel == Some(i), &s.contig).clicked() {
                                sel = Some(i);
                            }
                            ui.label(s.end_pos.to_string());
                            ui.label(s.num_reads.to_string());
                            ui.label(format!("{:.2}", s.mean_depth));
                            ui.label(format!("{:.1}%", s.coverage));
                            match c.contig_callable.get(i) {
                                Some(cm) => {
                                    ui.label(cm.callable.to_string());
                                    ui.label(cm.no_coverage.to_string());
                                    ui.label(cm.low_coverage.to_string());
                                    ui.label(cm.excessive_coverage.to_string());
                                    ui.label(cm.poor_mapping_quality.to_string());
                                }
                                None => {
                                    for _ in 0..5 {
                                        ui.label("–");
                                    }
                                }
                            }
                            ui.label(format!("{:.1}", s.mean_base_q));
                            ui.label(format!("{:.1}", s.mean_map_q));
                            ui.end_row();
                        }
                    });
                });

                // Histogram chart for the current selection (whole-genome or a contig).
                ui.separator();
                let (title, hist): (String, &[u64]) = match sel {
                    None => ("whole genome".to_string(), c.coverage_histogram.as_slice()),
                    Some(i) => (
                        c.contig_coverage_stats[i].contig.clone(),
                        c.contig_coverage_stats[i].histogram.as_slice(),
                    ),
                };
                if hist.is_empty() {
                    // No histogram persisted (fast-path / pipeline-sidecar import).
                    ui.label(format!(
                        "Per-contig depth histogram unavailable for {title} (pipeline-sidecar import) — \
                         the GATK CallableLoci breakdown is shown in the table above."
                    ));
                } else if hist.iter().skip(1).any(|&v| v > 0) {
                    coverage_histogram_chart(ui, hist, &title);
                } else {
                    ui.label(format!("{title}: no covered positions (every base at depth 0)."));
                }
            }
        }
        self.coverage_hist_contig = sel;

        if self.coverage.is_some() {
            self.export_row(
                ui,
                &[
                    navigator_app::ExportRequest::CoverageTsv(alignment_id),
                    navigator_app::ExportRequest::CoverageHtml(alignment_id),
                    navigator_app::ExportRequest::CallableBed(alignment_id),
                    navigator_app::ExportRequest::DiploidVcf(alignment_id),
                ],
            );
            self.publish_row(ui, "Publish summary to PDS", Command::PublishCoverage(alignment_id));
        }
    }

    /// Inferred sex + read-level QC metrics for a single alignment.
    pub(crate) fn sex_metrics_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam && !self.running_sex, egui::Button::new(self.tr("btn.inferSex"))).clicked() {
                self.running_sex = true;
                self.status = "Inferring sex…".into();
                let _ = self.tx.send(Command::RunSex(alignment_id));
            }
            if ui.add_enabled(has_bam && !self.running_metrics, egui::Button::new(self.tr("btn.readMetrics"))).clicked() {
                self.running_metrics = true;
                self.status = "Collecting read metrics…".into();
                let _ = self.tx.send(Command::RunReadMetrics(alignment_id));
            }
            if ui.add_enabled(has_bam && !self.running_sv, egui::Button::new(self.tr("btn.callSv"))).clicked() {
                self.running_sv = true;
                self.status = "Calling structural variants (needs ≥10× coverage)…".into();
                let _ = self.tx.send(Command::RunSv(alignment_id));
            }
            if self.running_sex || self.running_metrics || self.running_sv {
                ui.spinner();
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        if let Some(s) = &self.sex {
            let sex = match s.inferred_sex {
                navigator_app::InferredSex::Male => "Male",
                navigator_app::InferredSex::Female => "Female",
                navigator_app::InferredSex::Unknown => "Unknown",
            };
            ui.label(format!(
                "Sex: {sex}  ·  chrX:autosome ratio {:.2}  ·  {:?} confidence",
                s.x_autosome_ratio, s.confidence
            ));
        }
        if let Some(m) = &self.read_metrics {
            egui::Grid::new("read_metrics_grid").striped(true).num_columns(2).show(ui, |ui| {
                let row = |ui: &mut egui::Ui, k: &str, v: String| {
                    ui.label(k);
                    ui.label(v);
                    ui.end_row();
                };
                row(ui, "Total reads", m.total_reads.to_string());
                row(ui, "% PF aligned", format!("{:.1}%", m.pct_pf_reads_aligned * 100.0));
                row(ui, "% proper pairs", format!("{:.1}%", m.pct_proper_pairs * 100.0));
                row(ui, "Mean read length", format!("{:.0}", m.mean_read_length));
                row(ui, "Median insert size", format!("{:.0}", m.median_insert_size));
                row(ui, "Pair orientation", m.pair_orientation.as_str().to_string());
                row(ui, "Mean MAPQ", format!("{:.1}", m.mean_mapping_quality));
            });
        }
        if let Some(sv) = &self.sv {
            ui.label(format!(
                "Structural variants: {} calls ({} CNV segments, {} discordant pairs)",
                sv.sv_calls.len(),
                sv.cnv_segments,
                sv.total_discordant_pairs
            ));
            for c in sv.sv_calls.iter().take(8) {
                ui.label(
                    egui::RichText::new(format!(
                        "  {} {}:{}-{} {}bp q{:.0}",
                        c.sv_type.as_str(), c.chrom, c.start, c.end, c.sv_len, c.quality
                    ))
                    .small()
                    .weak(),
                );
            }
        }
        if self.read_metrics.is_some() {
            self.export_row(ui, &[navigator_app::ExportRequest::ReadMetricsTsv(alignment_id)]);
        }
    }

    /// A "Publish to PDS" button + sign-in hint/spinner, shared by the result sections.
    pub(crate) fn publish_row(&mut self, ui: &mut egui::Ui, label: &str, cmd: Command) {
        ui.horizontal(|ui| {
            let ready = self.account.is_some() && !self.publishing;
            if ui.add_enabled(ready, egui::Button::new(label)).clicked() {
                self.publishing = true;
                self.status = "Queueing for publish…".into();
                let _ = self.tx.send(cmd);
            }
            if self.account.is_none() {
                ui.label(self.tr("hint.signInToPublish"));
            }
            if self.publishing {
                ui.spinner();
            }
        });
    }

    /// One or more "Export" buttons in a row. Each opens a native save dialog (suggested filename +
    /// extension filter); on confirm it dispatches the export to the worker, which writes the file.
    pub(crate) fn export_row(&mut self, ui: &mut egui::Ui, requests: &[navigator_app::ExportRequest]) {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(self.tr("export.label")).weak());
            for req in requests {
                if ui.button(req.label()).clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .set_file_name(req.default_filename())
                        .add_filter(req.extension().to_uppercase(), &[req.extension()])
                        .save_file()
                    {
                        self.status = format!("Exporting {}…", req.label());
                        let _ = self.tx.send(Command::Export { request: *req, path });
                    }
                }
            }
        });
    }

    /// One DNA type's consensus row plus its override controls, audit log, and publish
    /// button. The consensus/audit are cloned up front so the form fields can be borrowed
    /// mutably for the override inputs.
    pub(crate) fn consensus_block(&mut self, ui: &mut egui::Ui, label: &str, dna_type: DnaType) {
        let cons = match dna_type {
            DnaType::Y => self.consensus_y.clone(),
            DnaType::Mt => self.consensus_mt.clone(),
        };
        let Some(c) = cons else { return };
        let Some(guid) = self.selected_sample else { return };

        let (compat, col) = match c.compatibility {
            CompatibilityLevel::Compatible => ("compatible", egui::Color32::from_rgb(60, 160, 60)),
            CompatibilityLevel::MinorDivergence => ("minor divergence", egui::Color32::from_rgb(170, 150, 40)),
            CompatibilityLevel::MajorDivergence => ("major divergence", egui::Color32::from_rgb(200, 120, 40)),
            CompatibilityLevel::Incompatible => ("incompatible", egui::Color32::from_rgb(200, 60, 60)),
        };
        ui.horizontal(|ui| {
            ui.strong(format!("{label}: {}", c.haplogroup));
            ui.label(format!("({} source(s), conf {:.3})", c.run_count, c.confidence));
            ui.colored_label(col, compat);
            if c.overridden {
                ui.colored_label(egui::Color32::from_rgb(120, 120, 220), "manual override");
            }
        });
        for w in &c.warnings {
            ui.label(format!("  ⚠ {w}"));
        }

        // Manual override: set a curator-corrected terminal (e.g. Sanger-confirmed), or clear.
        // Bind before the &mut self.forms borrow below (used inside the closure).
        let override_lbl = self.tr("form.override");
        let set_lbl = self.tr("common.set");
        let clear_lbl = self.tr("common.clear");
        let (hg_field, reason_field) = match dna_type {
            DnaType::Y => (&mut self.forms.override_y_haplogroup, &mut self.forms.override_y_reason),
            DnaType::Mt => (&mut self.forms.override_mt_haplogroup, &mut self.forms.override_mt_reason),
        };
        ui.horizontal(|ui| {
            ui.label(override_lbl);
            ui.add(egui::TextEdit::singleline(hg_field).hint_text("haplogroup").desired_width(140.0));
            ui.add(egui::TextEdit::singleline(reason_field).hint_text("reason").desired_width(180.0));
            let hg = hg_field.trim().to_string();
            let reason = reason_field.trim().to_string();
            if ui.add_enabled(!hg.is_empty(), egui::Button::new(set_lbl)).clicked() {
                self.status = format!("Overriding {label} consensus → {hg}");
                let _ = self.tx.send(Command::SetHaploOverride {
                    biosample_guid: guid,
                    dna_type,
                    haplogroup: hg,
                    reason: (!reason.is_empty()).then_some(reason),
                });
            }
            if ui.add_enabled(c.overridden, egui::Button::new(clear_lbl)).clicked() {
                self.status = format!("Clearing {label} override");
                let _ = self.tx.send(Command::ClearHaploOverride { biosample_guid: guid, dna_type });
            }
        });

        // mtDNA heteroplasmy observations from the last scan (folded into the published record).
        let het: Vec<HeteroplasmySite> = if dna_type == DnaType::Mt {
            self.heteroplasmy.as_ref().map(|(_, s)| s.clone()).unwrap_or_default()
        } else {
            Vec::new()
        };
        if dna_type == DnaType::Mt && !het.is_empty() {
            egui::CollapsingHeader::new(format!("heteroplasmy — {} site(s)", het.len()))
                .id_salt(("het", label))
                .show(ui, |ui| {
                    for h in &het {
                        ui.label(format!(
                            "  pos {}: {}/{} minor {:.1}% (depth {})",
                            h.position, h.major_base, h.minor_base, h.minor_fraction * 100.0, h.depth
                        ));
                    }
                });
        }

        // Audit log.
        let audit = match dna_type {
            DnaType::Y => &self.audit_y,
            DnaType::Mt => &self.audit_mt,
        };
        if !audit.is_empty() {
            egui::CollapsingHeader::new(format!("audit log — {} entr{}", audit.len(), if audit.len() == 1 { "y" } else { "ies" }))
                .id_salt(("audit", label))
                .show(ui, |ui| {
                    for e in audit {
                        ui.label(format!("  {} · {} — {}", e.timestamp, e.action, e.note));
                    }
                });
        }

        // Publish the donor-level reconciliation record (gated on sign-in).
        ui.horizontal(|ui| {
            let signed_in = self.account.is_some();
            if ui
                .add_enabled(signed_in && !self.publishing, egui::Button::new(format!("Publish {label} reconciliation")))
                .clicked()
            {
                self.publishing = true;
                self.status = format!("Publishing {label} reconciliation…");
                let _ = self.tx.send(Command::PublishReconciliation {
                    biosample_guid: guid,
                    dna_type,
                    heteroplasmy: het,
                    identity: self.identity.clone(),
                });
            }
            if !signed_in {
                ui.label(self.tr("hint.signInToPublish"));
            }
        });
        ui.add_space(6.0);
    }

    /// mtDNA heteroplasmy scan for an alignment (chrM pileup → mixed positions). Results
    /// feed the mtDNA reconciliation record's heteroplasmy observations.
    pub(crate) fn heteroplasmy_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new(self.tr("btn.scanHeteroplasmy"))).clicked() {
                self.status = "Scanning chrM pileup for heteroplasmy…".into();
                let _ = self.tx.send(Command::LoadHeteroplasmy { alignment_id });
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        if let Some((id, sites)) = &self.heteroplasmy {
            if *id == alignment_id {
                if sites.is_empty() {
                    ui.label(self.tr("mt.noHeteroplasmy"));
                } else {
                    ui.label(format!("{} heteroplasmic position(s):", sites.len()));
                    for h in sites {
                        ui.label(format!(
                            "  pos {}: {} (major) / {} (minor) — minor {:.1}%, depth {}",
                            h.position, h.major_base, h.minor_base, h.minor_fraction * 100.0, h.depth
                        ));
                    }
                }
            }
        }
    }

    /// Y-haplogroup assignment for an alignment (calls chrY tree positions; FTDNA tree).
    pub(crate) fn y_haplogroup_section(&mut self, ui: &mut egui::Ui, alignment_id: i64) {
        let has_bam = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some())
            .unwrap_or(false);

        ui.horizontal(|ui| {
            if ui.add_enabled(has_bam, egui::Button::new(self.tr("btn.assignY"))).clicked() {
                self.status = "Assigning Y haplogroup (fetching FTDNA tree)…".into();
                let _ = self.tx.send(Command::AssignYHaplogroup { alignment_id });
            }
            if ui.add_enabled(has_bam && !self.y_report_running, egui::Button::new(self.tr("haplo.fullReport"))).clicked() {
                self.y_report_running = true;
                self.status = "Building haplogroup report…".into();
                let _ = self.tx.send(Command::YHaploReport { alignment_id });
            }
            if self.y_report_running {
                ui.spinner();
            }
            if !has_bam {
                ui.label(self.tr("hint.noBamPath"));
            }
        });

        // The persisted donor consensus (reloaded on select) — so the haplogroup stays visible
        // without re-running. The fresh per-run assignment, when present, adds the SNP detail.
        if let Some(c) = &self.consensus_y {
            ui.label(
                egui::RichText::new(format!(
                    "Haplogroup: {}  ({} run(s), confidence {:.2})",
                    c.haplogroup, c.run_count, c.confidence
                ))
                .strong(),
            );
            if !c.lineage.is_empty() {
                ui.label(egui::RichText::new(c.lineage.join(" › ")).small().weak());
            }
        }
        if let Some((id, assignment)) = &self.y_haplogroup {
            if *id == alignment_id {
                show_assignment(ui, assignment);
            }
        }
        self.haplo_report_section(ui, alignment_id);

        // Private bucket: de-novo chrY calls off the assigned backbone (branch candidates).
        let has_ref = self
            .alignments
            .iter()
            .find(|a| a.id == alignment_id)
            .map(|a| a.bam_path.is_some() && a.reference_path.is_some())
            .unwrap_or(false);
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.label(self.tr("form.callableMask"));
            ui.checkbox(&mut self.y_self_mask, "self-referential (this sample)");
        });
        if !self.y_self_mask {
            ui.horizontal(|ui| {
                let label = self
                    .y_mask_path
                    .as_ref()
                    .and_then(|p| p.file_name())
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_else(|| "none (noisy)".into());
                ui.label(format!("External BED: {label}"));
                if ui.button(self.tr("form.chooseBed")).clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("BED", &["bed"]).pick_file() {
                        self.y_mask_path = Some(p);
                    }
                }
            });
        }
        ui.horizontal(|ui| {
            if ui.add_enabled(has_ref && !self.finding_private_y, egui::Button::new(self.tr("btn.findPrivateY"))).clicked() {
                self.finding_private_y = true;
                self.status = "Finding private Y variants (de-novo chrY)…".into();
                let mask = if self.y_self_mask {
                    YMask::SelfReferential
                } else {
                    self.y_mask_path.clone().map(YMask::Bed).unwrap_or(YMask::None)
                };
                let _ = self.tx.send(Command::FindPrivateY { alignment_id, mask });
            }
            if self.finding_private_y {
                ui.spinner();
            }
            if !has_ref {
                ui.label(self.tr("hint.needsBamRef"));
            }
        });
        if let Some((id, bucket)) = &self.private_y {
            if *id == alignment_id {
                ui.label(format!("{} novel + {} off-path, below {}", bucket.novel(), bucket.off_path(), bucket.terminal));
                egui::CollapsingHeader::new("Private variants").id_salt(("privy", alignment_id)).show(ui, |ui| {
                    egui::Grid::new(("privy_grid", alignment_id)).striped(true).num_columns(4).show(ui, |ui| {
                        for h in ["table.position", "table.change", "table.depth", "table.class"] {
                            ui.strong(self.tr(h));
                        }
                        ui.end_row();
                        for v in bucket.variants.iter().take(500) {
                            ui.label(v.position.to_string());
                            ui.label(format!("{}>{}", v.reference, v.alternate));
                            ui.label(v.depth.to_string());
                            match &v.class {
                                PrivateClass::Novel => ui.colored_label(egui::Color32::from_rgb(60, 160, 60), "novel"),
                                PrivateClass::OffPathKnown(name) => ui.label(format!("off-path: {name}")),
                            };
                            ui.end_row();
                        }
                    });
                    if bucket.variants.len() > 500 {
                        ui.label(format!("…and {} more", bucket.variants.len() - 500));
                    }
                });
            }
        }
    }

    /// The full Y-haplogroup placement report (gap §8): the ranked candidate haplogroups + the
    /// defining-SNP evidence along the reported lineage. Shown once "Full report" is run, flowing into
    /// the page scroll (no nested ScrollArea).
    fn haplo_report_section(&self, ui: &mut egui::Ui, alignment_id: i64) {
        let Some(r) = &self.y_report else { return };
        if r.alignment_id != alignment_id {
            return;
        }
        ui.add_space(4.0);
        egui::CollapsingHeader::new(self.tr("haplo.report")).id_salt(("haplo_report", alignment_id)).default_open(true).show(ui, |ui| {
            ui.label(egui::RichText::new(self.tr("haplo.candidates")).strong().small());
            egui::Grid::new(("haplo_ranked", alignment_id)).striped(true).num_columns(4).spacing([14.0, 2.0]).show(ui, |ui| {
                for h in ["Haplogroup", "Score", "Depth", "Matched/Expected"] {
                    ui.label(egui::RichText::new(h).strong().small());
                }
                ui.end_row();
                for c in r.assignment.ranked.iter().take(12) {
                    ui.label(&c.name);
                    ui.label(format!("{:.3}", c.score));
                    ui.label(c.depth.to_string());
                    ui.label(format!("{}/{}", c.matched, c.expected));
                    ui.end_row();
                }
            });
            ui.add_space(6.0);
            ui.label(egui::RichText::new(format!("{} ({})", self.tr("haplo.lineageSnps"), r.lineage.len())).strong().small());
            egui::Grid::new(("haplo_lineage", alignment_id)).striped(true).num_columns(3).spacing([14.0, 2.0]).show(ui, |ui| {
                for s in &r.lineage {
                    ui.label(&s.name);
                    ui.label(format!("{}{}>{}", s.position, s.ancestral, s.derived));
                    let (txt, col) = match s.state {
                        CallState::Derived => ("derived", egui::Color32::from_rgb(60, 160, 60)),
                        CallState::Ancestral => ("ancestral", egui::Color32::from_rgb(170, 120, 40)),
                        CallState::NoCall => ("no-call", egui::Color32::GRAY),
                    };
                    ui.colored_label(col, txt);
                    ui.end_row();
                }
            });
        });
    }

}
