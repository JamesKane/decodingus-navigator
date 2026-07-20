//! `impl NavigatorApp` methods extracted from `ui.rs` (the `events` group). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + helpers.
use super::*;

impl NavigatorApp {
    pub(crate) fn drain_events(&mut self) {
        while let Ok(event) = self.rx.try_recv() {
            match event {
                Event::Noop => {}
                Event::Overview(v) => {
                    self.status = format!("{} project(s)", v.len());
                    self.overview = v;
                }
                Event::AssetStatus(v) => self.asset_status = v,
                Event::ProjectCreated(p) => {
                    self.select_project(p.id);
                    let _ = self.tx.send(Command::LoadOverview);
                }
                Event::ProjectsChanged => {
                    let _ = self.tx.send(Command::LoadOverview);
                    let _ = self.tx.send(Command::LoadAllBiosamples); // a deleted project clears assignments
                }
                Event::ProjectImported(summary) => {
                    let mut msg = format!(
                        "Imported {}: {} sample(s), {} alignment(s)",
                        summary.project.name, summary.samples_total, summary.alignments_created
                    );
                    if summary.alignments_skipped > 0 {
                        msg.push_str(&format!(" ({} already present)", summary.alignments_skipped));
                    }
                    if !summary.missing_index.is_empty() {
                        msg.push_str(&format!("; {} sample(s) missing an index", summary.missing_index.len()));
                    }
                    // Per-sample failures that were skipped so the rest could import (recoverable).
                    if !summary.sample_errors.is_empty() {
                        msg.push_str(&format!(
                            "; {} sample(s) skipped on error: {}",
                            summary.sample_errors.len(),
                            summary.sample_errors.join(" | ")
                        ));
                    }
                    // Which reference each build resolved to (so a wrong/defaulted build is visible).
                    if !summary.reference_notes.is_empty() {
                        msg.push_str(&format!(". References: {}", summary.reference_notes.join("; ")));
                    }
                    // Fast path: what the pipeline sidecars filled without walking the CRAM.
                    let fp = &summary.fast_path;
                    if fp.samples_with_sidecars > 0 {
                        msg.push_str(&format!(
                            ". Fast path on {} sample(s): {} Y, {} mt, {} sex, {} metrics, {} coverage",
                            fp.samples_with_sidecars,
                            fp.y_placed,
                            fp.mt_placed,
                            fp.sex_filled,
                            fp.metrics_filled,
                            fp.coverage_filled,
                        ));
                        if !fp.errors.is_empty() {
                            msg.push_str(&format!(" ({} fast-path error(s))", fp.errors.len()));
                        }
                    }
                    self.status = msg;
                    self.importing = false;
                    self.pending_import_dir = None;
                    self.reference_needs.clear();
                    self.select_project(summary.project.id);
                    let _ = self.tx.send(Command::LoadOverview);
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                }
                Event::FtdnaPlan(plan) => {
                    let (new, merge, confirm) = plan.counts();
                    self.status = format!("FTDNA plan: {new} new, {merge} auto-merge, {confirm} to confirm");
                    self.importing = false;
                    self.ftdna_resolutions.clear();
                    self.ftdna_plan = Some(plan);
                }
                Event::FtdnaImported(summary) => {
                    self.status = format!(
                        "FTDNA import: {} merged, {} created, {} Y-STR, {} MDKA{}{}",
                        summary.merged,
                        summary.created,
                        summary.str_profiles,
                        summary.mdka_written,
                        if summary.skipped > 0 {
                            format!(", {} skipped", summary.skipped)
                        } else {
                            String::new()
                        },
                        if summary.errors.is_empty() {
                            String::new()
                        } else {
                            format!(" ({} error(s))", summary.errors.len())
                        },
                    );
                    self.ftdna_plan = None;
                    self.ftdna_resolutions.clear();
                    let _ = self.tx.send(Command::LoadOverview);
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                    // Surface the project the kits landed in (created or targeted).
                    if summary.project_id > 0 {
                        self.select_project(summary.project_id);
                    }
                    // Refresh the open subject's genealogy card if a merge/create touched it.
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadGenealogy(guid));
                    }
                }
                Event::Genealogy { guid, data } => {
                    if self.selected_sample == Some(guid) {
                        self.genealogy = Some((guid, data));
                    }
                }
                Event::ProjectClustering { project_id, clustering } => {
                    self.clustering_running = false;
                    if self.selected_project == Some(project_id) {
                        let (clusters, suggested) = (
                            clustering.clusters.len(),
                            clustering.clusters.iter().map(|c| c.suggested_count()).sum::<usize>(),
                        );
                        self.status =
                            format!("Y-STR clustering: {clusters} cluster(s), {suggested} branch suggestion(s)");
                        self.project_clustering = Some((project_id, clustering));
                    }
                }
                Event::ReferenceNeeded { dir, builds } => {
                    self.importing = false;
                    self.pending_import_dir = Some(dir);
                    self.status = format!("{} reference build(s) need downloading", builds.len());
                    self.reference_needs = builds;
                }
                Event::ReferenceProgress { build, received, total } => {
                    // Mirror the download into the always-visible status bar. The progress bar is only
                    // drawn in a couple of views (and none in Simple mode), so without this a slow
                    // multi-GB reference pull — kicked off in the background after import — looks like
                    // the app is stuck. The status line is the one surface visible in every view.
                    let recv_mb = received / 1_000_000;
                    self.status = match total {
                        Some(t) if t > 0 => format!(
                            "Downloading {build} reference in the background — {recv_mb} / {} MB…",
                            t / 1_000_000
                        ),
                        _ => format!("Downloading {build} reference in the background — {recv_mb} MB…"),
                    };
                    self.reference_progress = Some((build, received, total));
                }
                Event::ReferenceReady { build, path } => {
                    self.status = format!("Reference {build} ready ({})", path.display());
                    self.reference_progress = None;
                    self.reference_needs.retain(|b| b.build != build);
                    // When every needed build is in, retry the import automatically.
                    if self.reference_needs.is_empty() {
                        if let Some(dir) = self.pending_import_dir.take() {
                            self.importing = true;
                            self.status = format!("Importing {}…", dir.display());
                            let _ = self.tx.send(Command::ImportProjectDir { dir, reference: None });
                        }
                    }
                }
                Event::IndexProgress { file, done, total } => {
                    self.status = format!("Building index for {file}…");
                    self.index_progress = Some((file, done, total));
                }
                Event::IndexReady { built } => {
                    self.index_progress = None;
                    if built.is_some() {
                        self.status = "Index built.".into();
                    }
                }
                Event::UpdateAvailable(info) => {
                    self.status = format!("Update available: {} → {}", info.current_version, info.latest_version);
                    self.update_info = Some(*info);
                }
                Event::UpToDate => {
                    // Quietly current — no nagging. (A failed check surfaces via Event::Error.)
                }
                Event::Samples { project_id, samples } => {
                    if self.selected_project == Some(project_id) {
                        self.samples = samples;
                    }
                }
                Event::ProjectReport { project_id, rows } => {
                    if self.selected_project == Some(project_id) {
                        self.project_report = rows;
                    }
                }
                Event::ProjectStrChart { project_id, chart } => {
                    if self.selected_project == Some(project_id) {
                        self.project_str_chart = Some(chart);
                        self.project_str_loading = false;
                    }
                }
                Event::SubjectBrief { guid, brief } => {
                    if self.selected_sample == Some(guid) {
                        self.subject_brief = Some((guid, *brief));
                        self.subject_brief_loading = false;
                    }
                }
                Event::DescentReportLoaded { guid, dna, result } => {
                    self.descent_loading.retain(|(g, d)| !(*g == guid && *d == dna));
                    if self.selected_sample == Some(guid) {
                        match result {
                            Ok(report) => self.descent_reports.push((guid, dna, report)),
                            Err(msg) => self.status = format!("{} {msg}", self.tr("descent.failed")),
                        }
                    }
                }
                Event::BranchReportLoaded { guid, dna, result } => {
                    self.branch_loading.retain(|(g, d)| !(*g == guid && *d == dna));
                    if self.selected_sample == Some(guid) {
                        // Replace any prior report for this (guid, dna) — a new node was queried.
                        self.branch_reports.retain(|(g, d, _)| !(*g == guid && *d == dna));
                        match result {
                            Ok(report) => self.branch_reports.push((guid, dna, report)),
                            Err(msg) => self.status = format!("{} {msg}", self.tr("branch.failed")),
                        }
                    }
                }
                Event::BriefNarrationChunk { guid, text } => {
                    if self.selected_sample == Some(guid) {
                        match &mut self.narration_stream {
                            Some((g, buf)) if *g == guid => buf.push_str(&text),
                            _ => self.narration_stream = Some((guid, text)),
                        }
                    }
                }
                Event::BriefNarration { guid, result } => {
                    if self.selected_sample == Some(guid) {
                        self.narrating = false;
                        self.narration_stream = None; // the final result is authoritative
                        match result {
                            Ok(narration) => self.brief_narration = Some((guid, narration)),
                            // Fallback: keep the deterministic brief; surface why in the status line.
                            Err(msg) => {
                                self.brief_narration = None;
                                self.status = format!("{} {msg}", self.tr("brief.aiUnavailable"));
                            }
                        }
                    }
                }
                Event::ChatAnswerChunk { guid, text } => {
                    if self.selected_sample == Some(guid) {
                        // Append to the pending (last) assistant turn.
                        if let Some(turn) = self.chat_history.last_mut().filter(|t| !t.from_user) {
                            turn.text.push_str(&text);
                        }
                    }
                }
                Event::ChatAnswer { guid, result } => {
                    if self.selected_sample == Some(guid) {
                        self.chat_pending = false;
                        let text = match result {
                            Ok(answer) => answer,
                            Err(msg) => format!("{} {msg}", self.tr("brief.aiUnavailable")),
                        };
                        // Set the authoritative answer on the pending assistant turn (pre-pushed on
                        // send); fall back to appending one if it's missing.
                        match self.chat_history.last_mut().filter(|t| !t.from_user) {
                            Some(turn) => turn.text = text,
                            None => self.chat_history.push(ChatTurn { from_user: false, text }),
                        }
                    }
                }
                Event::SignalNarrationChunk { guid, kind, text } => {
                    if self.selected_sample == Some(guid) {
                        match &mut self.signal_stream {
                            Some((g, k, buf)) if *g == guid && *k == kind => buf.push_str(&text),
                            _ => self.signal_stream = Some((guid, kind, text)),
                        }
                    }
                }
                Event::SignalNarration { guid, kind, result } => {
                    if self.selected_sample == Some(guid) {
                        if self.signal_narrating == Some((guid, kind)) {
                            self.signal_narrating = None;
                        }
                        self.signal_stream = None; // the final result is authoritative
                        match result {
                            Ok(narration) => {
                                self.signal_narration.retain(|(g, k, _)| !(*g == guid && *k == kind));
                                self.signal_narration.push((guid, kind, narration));
                            }
                            // Fallback: keep the structured facts; surface why in the status line.
                            Err(msg) => self.status = format!("{} {msg}", self.tr("brief.aiUnavailable")),
                        }
                    }
                }
                Event::ProjectAnalyzed {
                    project_id,
                    samples,
                    coverage_done,
                    y_done,
                    sex_done,
                    metrics_done,
                    sv_done,
                    errors,
                    cancelled,
                } => {
                    self.analyzing = false;
                    self.cancelling = false;
                    self.deep_progress = None;
                    self.status = format!(
                        "{} {samples} sample(s): {coverage_done} coverage, {y_done} Y, {sex_done} sex, {metrics_done} metrics, {sv_done} SV{}",
                        if cancelled { "Deep analysis cancelled after" } else { "Analyzed" },
                        if errors > 0 { format!(", {errors} error(s)") } else { String::new() }
                    );
                    if self.selected_project == Some(project_id) {
                        let _ = self.tx.send(Command::LoadProjectReport(project_id));
                        // Haplogroups may have been assigned — regroup the STR chart.
                        self.reload_project_str();
                    }
                    // Coverage was (re)computed — refresh the subjects-list Status column.
                    let _ = self.tx.send(Command::LoadSubjectStatus);
                }
                Event::DeepAnalyzeProgress {
                    project_id,
                    done,
                    total,
                    sample,
                    fraction,
                } => {
                    if self.selected_project == Some(project_id) {
                        self.deep_progress = Some((done, total, sample, fraction));
                    }
                }
                Event::ImportProgress {
                    done,
                    total,
                    sample,
                    fraction,
                } => {
                    self.importing = true;
                    let pct = (fraction * 100.0).round() as u32;
                    self.status = format!("Importing: {done}/{total} ({pct}%) — {sample}…");
                }
                Event::AllBiosamples(v) => {
                    // Drop a dangling selection: after deleting the last subject the async list reload
                    // lands here empty, but `selected_sample` may still point at the deleted (or any
                    // now-removed) subject. Left set, the per-frame auto-select and brief-load keep
                    // re-fetching a brief that errors — never clearing `subject_brief_loading` — so the
                    // Simple view spins on "Building your brief…" forever. Clear it so the empty-state
                    // (or a valid re-selection) renders instead.
                    if let Some(sel) = self.selected_sample {
                        if !v.iter().any(|b| b.guid == sel) {
                            self.selected_sample = None;
                            self.subject_brief = None;
                            self.subject_brief_loading = false;
                        }
                    }
                    self.all_biosamples = v;
                    let _ = self.tx.send(Command::LoadHaploSummary); // fill the Y/mt columns
                    let _ = self.tx.send(Command::LoadSubjectStatus); // fill the Status column
                }
                Event::HaploSummary(map) => self.haplo_summary = map,
                Event::SubjectStatus(map) => self.subject_status = map,
                Event::BiosamplesChanged => {
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadSamples(pid));
                    }
                    let _ = self.tx.send(Command::LoadOverview); // project sample counts changed
                }
                Event::Runs { biosample_guid, runs } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.runs = runs;
                    }
                }
                Event::RunsChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadRuns(guid));
                    }
                }
                Event::BiosampleDataCleared(guid) => {
                    self.status = self.tr("clear.done").to_string();
                    // Fully reload the subject view from the now-empty DB + refresh the list columns.
                    if self.selected_sample == Some(guid) {
                        self.select_sample(guid);
                    }
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                }
                Event::HaplogroupDataReset(guid) => {
                    self.status = self.tr("resetHaplo.done").to_string();
                    // Reload the subject (re-reads the now-empty placement → brief refreshes).
                    if self.selected_sample == Some(guid) {
                        self.select_sample(guid);
                    }
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                }
                Event::StrProfiles {
                    biosample_guid,
                    profiles,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.str_profiles = profiles;
                        // Reset the report view-state for the new subject's data.
                        self.str_provider = None;
                        self.str_marker_filter.clear();
                    }
                }
                Event::StrProfilesChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadStrProfiles(guid));
                    }
                    // A member's STR data changed — refresh the project chart (best-effort; the
                    // builder only includes members of the open project).
                    self.reload_project_str();
                    self.status = "STR profile imported".into();
                }
                Event::VariantSets { biosample_guid, sets } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.variant_sets = sets;
                    }
                }
                Event::VariantSetsChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadVariantSets(guid));
                    }
                    self.status = "Variants imported".into();
                }
                Event::ChipProfiles {
                    biosample_guid,
                    profiles,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.chip_profiles = profiles;
                    }
                }
                Event::ChipProfilesChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadChipProfiles(guid));
                        // A chip import also places Y (and, for 23andMe, mtDNA) haplogroups —
                        // refresh the consensus so they appear without a manual reload.
                        let _ = self.tx.send(Command::LoadConsensus(guid));
                    }
                    let _ = self.tx.send(Command::LoadHaploSummary); // subjects-list Y/mt columns
                    self.status = "Chip data imported".into();
                }
                Event::MtdnaSequences {
                    biosample_guid,
                    sequences,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.mtdna_sequences = sequences;
                    }
                }
                Event::MtdnaChanged(guid) => {
                    if self.selected_sample == Some(guid) {
                        let _ = self.tx.send(Command::LoadMtdna(guid));
                    }
                    self.status = "mtDNA sequence imported".into();
                }
                Event::MtdnaVariants { mtdna_id, variants } => {
                    self.status = format!("mtDNA: {} mutations vs rCRS", variants.len());
                    self.mtdna_variants.insert(mtdna_id, variants);
                }
                Event::Haplogroup { mtdna_id, assignment } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("mtDNA haplogroup: {} (score {:.3})", top.name, top.score),
                        None => "No haplogroup match".into(),
                    };
                    self.mtdna_haplogroup = Some((mtdna_id, assignment));
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensus(guid)); // a call was recorded
                    }
                }
                Event::YHaplogroup {
                    alignment_id,
                    assignment,
                } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("Y haplogroup: {} (score {:.3})", top.name, top.score),
                        None => "No Y haplogroup match".into(),
                    };
                    self.y_haplogroup = Some((alignment_id, assignment));
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensus(guid));
                    }
                    // A per-row "Assign Y" from the project report just recorded a call —
                    // refresh the report so its Y column fills in.
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::YHaploReport {
                    alignment_id,
                    assignment,
                    lineage,
                } => {
                    self.y_report_running = false;
                    self.status = format!(
                        "Haplogroup report: {} candidate(s), {} lineage SNP(s)",
                        assignment.ranked.len(),
                        lineage.len()
                    );
                    self.y_report = Some(YReport {
                        alignment_id,
                        assignment,
                        lineage,
                    });
                }
                Event::YBisdnaHaplogroup {
                    biosample_guid,
                    assignment,
                } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("Y haplogroup (panel): {} (score {:.3})", top.name, top.score),
                        None => "No Y haplogroup match from the panel".into(),
                    };
                    // The call was recorded — refresh the donor consensus so the Y-DNA card fills in.
                    let _ = self.tx.send(Command::LoadConsensus(biosample_guid));
                }
                Event::MtHaplogroup {
                    alignment_id,
                    assignment,
                } => {
                    self.status = match assignment.ranked.first() {
                        Some(top) => format!("mtDNA haplogroup: {} (score {:.3})", top.name, top.score),
                        None => "No mtDNA haplogroup match".into(),
                    };
                    self.mt_haplogroup = Some((alignment_id, assignment));
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensus(guid)); // the mt call was recorded
                    }
                }
                Event::AncestryPainting { alignment_id, segments } => {
                    self.painting_running = false;
                    self.status = format!("Painted {} ancestry segments", segments.len());
                    self.painting = Some((alignment_id, segments));
                }
                Event::Consensus { biosample_guid, y, mt } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.consensus_y = y;
                        self.consensus_mt = mt;
                        // Consensus drives the Simple-mode brief — (re)build it now (no-op in Advanced).
                        self.reload_subject_brief();
                    }
                }
                Event::Audit {
                    biosample_guid,
                    dna_type,
                    entries,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        match dna_type {
                            DnaType::Y => self.audit_y = entries,
                            DnaType::Mt => self.audit_mt = entries,
                        }
                    }
                }
                Event::Heteroplasmy { alignment_id, sites } => {
                    self.status = format!("mtDNA heteroplasmy: {} site(s)", sites.len());
                    self.heteroplasmy = Some((alignment_id, sites));
                }
                Event::ReconciliationChanged {
                    biosample_guid,
                    dna_type,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        let _ = self.tx.send(Command::LoadConsensus(biosample_guid));
                        let _ = self.tx.send(Command::LoadAudit {
                            biosample_guid,
                            dna_type,
                        });
                    }
                    // An assigned haplogroup changed — regroup the project STR chart.
                    if matches!(dna_type, DnaType::Y) {
                        self.reload_project_str();
                    }
                }
                Event::PrivateY { alignment_id, bucket } => {
                    self.status = format!("Private Y: {} novel, {} off-path", bucket.novel(), bucket.off_path());
                    self.private_y = Some((alignment_id, bucket));
                    self.finding_private_y = false;
                    // A fresh (self-masked) bucket was just cached — refresh the donor union.
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadDonorPrivateY { biosample_guid: guid });
                    }
                }
                Event::LabsResolved(count) => {
                    if count > 0 {
                        self.status = format!("Resolved {count} sequencing lab(s) from instrument ids");
                        // Refresh the open subject's run cards so the new lab chips appear.
                        if let Some(guid) = self.selected_sample {
                            let _ = self.tx.send(Command::LoadRuns(guid));
                        }
                    }
                }
                Event::ReferenceSettings(rows) => {
                    self.settings_form.references = rows
                        .into_iter()
                        .map(|r: RefBuildStatus| RefRow {
                            build: r.build,
                            status: r.status,
                            local_path: r.local_path.unwrap_or_default(),
                            auto_download: r.auto_download,
                            verify: String::new(),
                        })
                        .collect();
                }
                Event::LlmConnection(result) => {
                    self.llm_testing = false;
                    match result {
                        Ok(models) => {
                            let n = models.len();
                            // Prefill the model field with the server's single loaded model.
                            if self.settings_form.llm_model.trim().is_empty() && models.len() == 1 {
                                self.settings_form.llm_model = models[0].clone();
                            }
                            self.llm_models = models;
                            self.llm_test_msg = Some(format!("{} {}", self.tr("settings.ai.ok"), n));
                        }
                        Err(msg) => {
                            self.llm_models.clear();
                            self.llm_test_msg = Some(msg);
                        }
                    }
                }
                Event::ReferenceSettingsChanged => {
                    self.status = "Reference settings saved".into();
                    let _ = self.tx.send(Command::LoadReferenceSettings); // refresh statuses
                }
                Event::ReferenceVerified { build, status } => {
                    if let Some(row) = self.settings_form.references.iter_mut().find(|r| r.build == build) {
                        row.verify = status;
                    }
                }
                Event::VcfLifted { summary } => {
                    self.status = summary;
                }
                Event::DataBatchImported {
                    biosample_guid,
                    summary,
                } => {
                    self.status = format!(
                        "Imported {} file(s){}",
                        summary.imported.len(),
                        if summary.skipped.is_empty() {
                            String::new()
                        } else {
                            format!(", {} skipped", summary.skipped.len())
                        }
                    );
                    self.batch_import = Some(summary);
                    // The import may have added an alignment — refresh the analysis-status map so the
                    // Subjects Status column and the Simple-mode "Analyze" prompt (`Pending` = has data,
                    // not analyzed) pick it up. Without this, adding data to an existing subject leaves
                    // both stale, so the analyze prompt never appears.
                    let _ = self.tx.send(Command::LoadSubjectStatus);
                    if self.selected_sample == Some(biosample_guid) {
                        let _ = self.tx.send(Command::LoadRuns(biosample_guid));
                        let _ = self.tx.send(Command::LoadStrProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadVariantSets(biosample_guid));
                        let _ = self.tx.send(Command::LoadChipProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadMtdna(biosample_guid));
                        // Rebuild the brief so the "Your test" card + not-analyzed state reflect the
                        // new file (Simple mode was showing the stale empty-subject brief).
                        self.reload_subject_brief();
                    }
                }
                Event::SubjectCreatedAndImported {
                    biosample_guid,
                    summary,
                } => {
                    self.status = format!("Created subject and imported {} file(s)", summary.imported.len());
                    self.batch_import = Some(summary);
                    // Refresh the list so the new subject appears, then select it — `select_sample`
                    // loads its runs/profiles and (in Simple mode) triggers the brief build.
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                    let _ = self.tx.send(Command::LoadOverview);
                    self.forms.show_add_subject = false;
                    self.select_sample(biosample_guid);
                }
                Event::StrConcordance {
                    biosample_guid,
                    alignment_id,
                    rows,
                } => {
                    self.str_running = false;
                    let calls = rows.iter().filter(|r| r.called.is_some()).count();
                    let agree = rows.iter().filter(|r| r.agree).count();
                    self.status = format!("Y-STR from sequence: {calls} markers called, {agree} agree with vendor");
                    self.str_concordance = Some((biosample_guid, alignment_id, rows));
                }
                Event::YMatches {
                    biosample_guid,
                    matches,
                } => {
                    self.y_matches_running = false;
                    self.status = format!("Y matches: {} ranked", matches.len());
                    self.y_matches = Some((biosample_guid, matches));
                }
                Event::DefaultAlignment { run_id, alignment_id } => {
                    // Only auto-select if the user hasn't already chosen an alignment.
                    if self.selected_alignment.is_none() {
                        self.pending_alignment = Some(alignment_id);
                        self.select_run(run_id); // loads the run's alignments → applied below
                    }
                }
                Event::DonorAncestry { alignment_id, result } => {
                    self.estimating_donor_ancestry = false;
                    self.donor_ancestry = Some((alignment_id, result));
                    // A fresh consensus estimate persisted the detailed methods too — refresh them.
                    if let Some(g) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadConsensusAncestryDetail { biosample_guid: g });
                    }
                }
                Event::ConsensusAncestryDetail {
                    biosample_guid,
                    fine,
                    ancient,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.fine_ancestry = fine.map(|b| *b);
                        self.ancient_ancestry = ancient.map(|b| *b);
                    }
                }
                Event::DeepAncestryEstimated { biosample_guid, result } => {
                    self.estimating_deep_ancestry = false;
                    if self.selected_sample == Some(biosample_guid) {
                        match result {
                            Some(r) => {
                                self.ancient_ancestry = Some(*r);
                                self.status = "Deep ancestry estimated.".into();
                            }
                            None => {
                                self.ancient_ancestry = None;
                                self.status =
                                    "Deep ancestry: not applicable (non-European, no whole-genome CHM13 alignment, or model rejected).".into();
                            }
                        }
                    }
                }
                Event::DonorPrivateY { bucket } => {
                    self.donor_private_y = Some(bucket);
                    self.y_snp_names_requested = false; // re-resolve names incl. the new positions
                }
                Event::YProfile {
                    biosample_guid,
                    profile,
                } => {
                    self.y_profile_loading = false;
                    if self.selected_sample == Some(biosample_guid) {
                        if let Some(p) = &profile {
                            self.status = format!("Y variant profile: {} variants", p.summary.total);
                        }
                        self.y_profile = profile;
                        self.y_snp_names_requested = false; // re-resolve names incl. the new positions
                                                            // A rebuild re-places the genome consensus (consensus_label); refresh the
                                                            // Overview's cached Y/mt consensus so it doesn't lag until the next reload.
                        let _ = self.tx.send(Command::LoadConsensus(biosample_guid));
                        // The descent report is drawn from this profile — drop its cache so it rebuilds.
                        self.descent_reports.retain(|(g, d, _)| !(*g == biosample_guid && *d == DnaType::Y));
                    }
                }
                Event::YSnpNames { names } => {
                    self.y_snp_names = names;
                }
                Event::MtProfile {
                    biosample_guid,
                    profile,
                } => {
                    self.mt_profile_loading = false;
                    if self.selected_sample == Some(biosample_guid) {
                        if let Some(p) = &profile {
                            self.status = format!("mtDNA consensus profile: {} mutations", p.summary.total);
                        }
                        self.mt_profile = profile;
                        // A rebuild re-places the mt genome consensus; refresh the Overview's cache.
                        let _ = self.tx.send(Command::LoadConsensus(biosample_guid));
                        // The descent report is drawn from this profile — drop its cache so it rebuilds.
                        self.descent_reports.retain(|(g, d, _)| !(*g == biosample_guid && *d == DnaType::Mt));
                    }
                }
                Event::AutosomalProfile {
                    biosample_guid,
                    profile,
                } => {
                    self.auto_profile_loading = false;
                    if self.selected_sample == Some(biosample_guid) {
                        if let Some(p) = &profile {
                            self.status = format!("autosomal consensus profile: {} sites", p.summary.total);
                        }
                        self.auto_profile = profile;
                    }
                }
                Event::Alignments {
                    sequence_run_id,
                    alignments,
                } => {
                    if self.selected_run == Some(sequence_run_id) {
                        self.alignments = alignments;
                        // Load cached coverage for every alignment so each Data Sources row shows
                        // coverage/callable without first being selected.
                        let ids: Vec<i64> = self.alignments.iter().map(|a| a.id).collect();
                        if !ids.is_empty() {
                            let _ = self.tx.send(Command::LoadCoverageBulk(ids));
                        }
                        // Apply a queued subject-default alignment once its run's list is loaded.
                        if let Some(pid) = self.pending_alignment {
                            if self.alignments.iter().any(|a| a.id == pid) {
                                self.pending_alignment = None;
                                self.select_alignment(pid);
                            }
                        }
                    }
                }
                Event::CoverageBulk(results) => {
                    for (id, result) in results {
                        if let Some(c) = result {
                            self.coverage_by_aln.insert(id, c);
                        }
                    }
                }
                Event::GenomeRegions { alignment_id, regions } => {
                    self.loading_regions = false;
                    if self.selected_alignment == Some(alignment_id) {
                        self.genome_regions = regions.map(|r| (alignment_id, r));
                    }
                }
                Event::AlignmentProbe(p) => {
                    // Auto-fill the add-alignment form from the BAM/CRAM header.
                    if let Some(b) = p.reference_build {
                        self.forms.aln_reference_build = b;
                    }
                    if let Some(a) = p.aligner {
                        self.forms.aln_aligner = a;
                    }
                    let bits: Vec<String> = [p.platform, p.instrument_model, p.test_type]
                        .into_iter()
                        .flatten()
                        .collect();
                    if !bits.is_empty() {
                        self.status = format!("Detected from header: {}", bits.join(" · "));
                    }
                }
                Event::AlignmentsChanged(run_id) => {
                    if self.selected_run == Some(run_id) {
                        let _ = self.tx.send(Command::LoadAlignments(run_id));
                    }
                    let _ = self.tx.send(Command::LoadAllAlignments); // keep IBD picker current
                }
                Event::Coverage { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.coverage = result.clone();
                        self.coverage_hist_contig = None; // reset histogram selection to whole-genome
                    }
                    // Keep the per-row map current after a (re)compute.
                    match result {
                        Some(c) => {
                            self.coverage_by_aln.insert(alignment_id, c);
                        }
                        None => {
                            self.coverage_by_aln.remove(&alignment_id);
                        }
                    }
                    self.running = false;
                    // A recompute (possibly from the project report) may have filled a cell.
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::Sex { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.sex = result;
                    }
                    self.running_sex = false;
                    // Sex inference may have written the sex back to the biosample — reload the
                    // subjects list so the table + header reflect it instead of "Unknown".
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::ReadMetrics { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.read_metrics = result;
                    }
                    self.running_metrics = false;
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::Sv { alignment_id, result } => {
                    if self.selected_alignment == Some(alignment_id) {
                        self.sv = result;
                    }
                    self.running_sv = false;
                    if let Some(pid) = self.selected_project {
                        let _ = self.tx.send(Command::LoadProjectReport(pid));
                    }
                }
                Event::Denovo {
                    alignment_id,
                    contig,
                    result,
                } => {
                    if self.selected_alignment == Some(alignment_id) {
                        match result {
                            Some(calls) => {
                                self.denovo.insert(contig, calls);
                            }
                            None => {
                                self.denovo.remove(&contig);
                            }
                        }
                    }
                    self.running_denovo = false;
                }
                Event::AnalysisProgress {
                    step,
                    total,
                    label,
                    detail,
                    fraction,
                } => {
                    // Reset the elapsed timer only when the step changes (sub-progress within a
                    // step keeps the same start time).
                    let started = match &self.analysis {
                        Some(a) if a.step == step => a.started,
                        _ => self.frame_time,
                    };
                    self.analysis = Some(AnalysisModal {
                        step,
                        total,
                        label,
                        detail,
                        fraction,
                        started,
                    });
                }
                Event::AnalysisDone { cancelled } => {
                    self.analysis = None;
                    self.cancelling = false;
                    self.status = if cancelled {
                        "Full analysis cancelled.".into()
                    } else {
                        "Full analysis complete.".into()
                    };
                    // The subject's coverage just changed — refresh the Status column.
                    let _ = self.tx.send(Command::LoadSubjectStatus);
                    // In Simple mode, rebuild the brief so the just-computed lineages/ancestry replace
                    // the "not analyzed yet" prompt (no-op in Advanced / when nothing is selected).
                    self.reload_subject_brief();
                }
                Event::AllAlignments(a) => self.all_alignments = a,
                Event::Ibd(cmp) => {
                    self.ibd_result = Some(cmp);
                    self.running_ibd = false;
                }
                Event::IbdSuggestions(items) => {
                    self.status = format!("{} network match suggestion(s)", items.len());
                    self.ibd_suggestions = items;
                    self.loading_ibd_suggestions = false;
                }
                Event::IbdIntroduced {
                    suggested_sample_guid,
                    request_uri,
                    status,
                } => {
                    self.status = format!("Introduction requested: {status} ({request_uri})");
                    let label = if request_uri.is_empty() {
                        status
                    } else {
                        format!("{status} · {request_uri}")
                    };
                    self.ibd_intros.insert(suggested_sample_guid, label);
                }
                Event::ExchangeInbox { incoming, ready } => {
                    self.exchange_busy = false;
                    self.status = format!(
                        "Exchange inbox: {} request(s), {} ready session(s)",
                        incoming.len(),
                        ready.len()
                    );
                    self.exchange_incoming = incoming;
                    self.exchange_ready = ready;
                }
                Event::ExchangeConsented => {
                    self.exchange_busy = false;
                    self.status = "Consent recorded".into();
                    let _ = self.tx.send(Command::ExchangeInbox); // refresh
                }
                Event::IbdExchangeDone {
                    biosample_guid,
                    total_shared_cm,
                    segment_count,
                    relationship,
                    agreed,
                } => {
                    self.exchange_busy = false;
                    self.status = format!(
                        "IBD exchange: {total_shared_cm:.1} cM, {segment_count} segment(s), {relationship}{}",
                        if agreed { " · agreed" } else { " · NOT agreed" }
                    );
                    let _ = self.tx.send(Command::LoadIbdExchanges { biosample_guid });
                }
                Event::IbdExchanges { biosample_guid, rows } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.exchange_results = rows;
                    }
                }
                Event::Identity(v) => {
                    self.status = format!("Identity: {:?} ({} sites)", v.status, v.sites_compared);
                    self.identity = Some(v);
                }
                Event::Authenticated(account) => {
                    self.status = match &account {
                        Some(did) => format!("Signed in as {did}"),
                        None => "Signed out".into(),
                    };
                    let signed_in = account.is_some();
                    self.account = account;
                    self.logging_in = false;
                    if signed_in {
                        // Flush anything queued while signed out / on a previous session.
                        let _ = self.tx.send(Command::DrainOutbox);
                    } else {
                        self.sync_pending = 0;
                    }
                }
                Event::Published { kind, uri } => {
                    self.status = format!("Published {kind}: {uri}");
                    self.publishing = false;
                }
                Event::Queued { kind } => {
                    // The publish is durably queued; it sends now if online, else on reconnect.
                    self.status = format!("Queued {kind} for publish");
                    self.publishing = false;
                }
                Event::SyncPending(n) => self.sync_pending = n,
                Event::Exported { label, path } => {
                    self.status = format!("Exported {label} → {}", path.display());
                }
                Event::SyncOnline(online) => self.online = online,
                Event::PullDone {
                    in_sync,
                    applied,
                    adopted,
                    repushed,
                    conflicts,
                } => {
                    self.pulling = false;
                    self.status = format!(
                        "Pull: {in_sync} in sync, {applied} applied, {adopted} remote-only, {repushed} to re-publish, {conflicts} conflict(s)"
                    );
                }
                Event::PcaReference { alignment_id, points } => {
                    self.pca_reference = Some((alignment_id, points));
                }
                Event::SourceFilesVerified { missing } => {
                    // Don't clobber a live import's progress status with this workspace-wide sweep
                    // (the sweep and the import are unrelated; overwriting made imports look stalled).
                    if !self.importing {
                        self.status = if missing == 0 {
                            "All source files present".into()
                        } else {
                            format!("{missing} source file(s) moved or missing")
                        };
                    }
                }
                // ---- social (Community tab) --------------------------------
                Event::SupportThreads(items) => self.support_threads = items,
                Event::SupportThread {
                    conversation_id,
                    messages,
                } => {
                    self.open_thread = Some((conversation_id, messages));
                }
                Event::SupportThreadPosted { conversation_id } => {
                    // Reload the list, and refresh the open thread (a reply) or open the new one.
                    let _ = self.tx.send(Command::LoadSupportThreads);
                    let _ = self.tx.send(Command::LoadSupportThread { conversation_id });
                }
                Event::CommunityFeed(feed) => self.feed = Some(feed),
                Event::CommunityPosted => {
                    self.status = self.tr("community.posted").to_string();
                    let _ = self.tx.send(Command::LoadCommunityFeed);
                }
                Event::Notifications { items, unread } => {
                    self.notifications = items;
                    self.notif_unread = unread;
                }
                Event::NotificationsMarked => {
                    let _ = self.tx.send(Command::LoadNotifications);
                }
                Event::DmInitiated => {
                    self.dm_partner_did.clear();
                    self.status = self.tr("dm.requestSent").to_string();
                    let _ = self.tx.send(Command::LoadDmInbox);
                }
                Event::DmInbox { incoming, ready } => {
                    self.dm_incoming = incoming;
                    self.dm_ready = ready;
                }
                Event::DmConsented => {
                    self.status = self.tr("dm.consentRecorded").to_string();
                    let _ = self.tx.send(Command::LoadDmInbox);
                    let _ = self.tx.send(Command::LoadDmConversations);
                }
                Event::DmConnected => {
                    self.status = self.tr("dm.connected").to_string();
                    let _ = self.tx.send(Command::LoadDmInbox);
                    let _ = self.tx.send(Command::LoadDmConversations);
                }
                Event::DmConversations(rows) => self.dm_conversations = rows,
                Event::DmMessages { session_id, rows } => {
                    self.open_dm = Some((session_id, rows));
                    let _ = self.tx.send(Command::LoadDmConversations); // unread cleared on read
                }
                Event::DmSent { session_id } => {
                    self.dm_compose.clear();
                    let _ = self.tx.send(Command::LoadDmMessages { session_id });
                }
                Event::DmSynced { session_id, new_count } => {
                    self.status = if new_count > 0 {
                        format!("{} {}", new_count, self.tr("dm.newMessages"))
                    } else {
                        self.tr("dm.upToDate").to_string()
                    };
                    if self.open_dm.as_ref().is_some_and(|(s, _)| *s == session_id) {
                        let _ = self.tx.send(Command::LoadDmMessages { session_id });
                    }
                    let _ = self.tx.send(Command::LoadDmConversations);
                }
                Event::RecruitmentInvitations(items) => self.recruitment_invitations = items,
                Event::RecruitmentResponded => {
                    self.status = self.tr("recruit.responded").to_string();
                    let _ = self.tx.send(Command::LoadRecruitmentInvitations);
                    let _ = self.tx.send(Command::LoadNotifications);
                }
                Event::TreesRefreshed(n) => {
                    // Drop the interpreted caches so open profiles/descent re-interpret against the
                    // freshly-pulled tree (observation-first: no re-genotyping needed).
                    self.status = format!("Refreshed haplotrees ({n} cached file(s) cleared) — re-interpreting");
                    self.descent_reports.clear();
                    self.descent_loading.clear();
                    self.branch_reports.clear();
                    self.branch_loading.clear();
                    self.y_profile = None;
                    self.mt_profile = None;
                    if let Some(guid) = self.selected_sample {
                        let _ = self.tx.send(Command::LoadYProfile { biosample_guid: guid });
                        let _ = self.tx.send(Command::LoadMtProfile { biosample_guid: guid });
                    }
                }
                Event::Error(e) => {
                    self.status = format!("Error: {e}");
                    // This failure carries no file-level cause; drop any report from a previous
                    // one so the status bar can't offer a "Details" that describes the wrong error.
                    self.diagnosis = None;
                    self.show_diagnosis = false;
                    self.clear_in_flight();
                }
                Event::Cancelled => {
                    self.status = self.tr("analysis.cancelled").to_string();
                    self.clear_in_flight();
                }
                Event::Diagnosed { message, report } => {
                    self.status = format!("Error: {message}");
                    self.diagnosis = Some(report);
                    // Open it unprompted: the whole point is that the one-line message is the part
                    // that isn't actionable, so making the user go find the detail would reproduce
                    // the original problem.
                    self.show_diagnosis = true;
                    self.clear_in_flight();
                }
            }
        }
    }

    pub(crate) fn select_project(&mut self, id: i64) {
        self.selected_project = Some(id);
        self.samples.clear();
        self.project_report.clear();
        self.project_str_chart = None;
        self.project_clustering = None; // stale for the new project until recomputed
        self.clustering_running = false;
        self.clear_sample_selection();
        let _ = self.tx.send(Command::LoadSamples(id));
        let _ = self.tx.send(Command::LoadProjectReport(id));
        self.reload_project_str();
    }

    /// (Re)build the Y-STR overview chart for the open project off the UI thread. Called on project
    /// select and whenever a member's STR data or assigned haplogroup changes (so the grouping stays
    /// in sync). No-op when no project is open.
    pub(crate) fn reload_project_str(&mut self) {
        if let Some(id) = self.selected_project {
            self.project_str_loading = true;
            let _ = self.tx.send(Command::LoadProjectStrChart(id));
        }
    }

    /// (Re)build the Simple-mode Subject Brief off the UI thread. Only meaningful in Simple mode;
    /// called on subject select and whenever the subject's haplogroups/coverage change. No-op when
    /// no subject is selected. Cheap (cache reads + pack lookups).
    pub(crate) fn reload_subject_brief(&mut self) {
        if self.ui_mode != UiMode::Simple {
            return;
        }
        if let Some(guid) = self.selected_sample {
            self.subject_brief_loading = true;
            let _ = self.tx.send(Command::LoadSubjectBrief(guid));
        }
    }

    /// Open a subject from a project's report row: select it, switch to the Subjects view, and
    /// remember the project so the detail header's "back to project" button can return there.
    pub(crate) fn open_sample_from_project(&mut self, guid: SampleGuid) {
        let pid = self.selected_project;
        self.select_sample(guid); // clears return_to_project
        self.return_to_project = pid;
        self.nav = Nav::Subjects;
    }

    pub(crate) fn select_sample(&mut self, guid: SampleGuid) {
        self.selected_sample = Some(guid);
        // A plain selection isn't "from a project" — the project opener re-sets this after.
        self.return_to_project = None;
        self.y_sub = YSub::default();
        self.y_snp_sub = YSnpSub::default();
        self.mt_sub = MtSub::default();
        self.auto_sub = AutoSub::default();
        self.y_snp_names.clear();
        self.y_snp_names_requested = false;
        self.pending_alignment = None;
        self.donor_ancestry = None;
        self.fine_ancestry = None;
        self.ancient_ancestry = None;
        // pca_reference is the global CHM13 centroid cloud (subject-independent) — keep it loaded
        // across subject switches rather than re-fetching the asset each time.
        self.estimating_donor_ancestry = false;
        self.painting = None;
        self.painting_running = false;
        self.donor_private_y = None;
        self.y_profile = None;
        self.y_profile_loading = false;
        self.mt_profile = None;
        self.mt_profile_loading = false;
        self.auto_profile = None;
        self.auto_profile_loading = false;
        self.genealogy = None;
        self.clear_run_selection();
        self.runs.clear();
        self.str_profiles.clear();
        self.variant_sets.clear();
        self.chip_profiles.clear();
        self.mtdna_sequences.clear();
        self.mtdna_variants.clear();
        self.mtdna_haplogroup = None;
        self.consensus_y = None;
        self.consensus_mt = None;
        self.descent_reports.clear();
        self.descent_loading.clear();
        self.branch_reports.clear();
        self.branch_loading.clear();
        self.branch_node_y.clear();
        self.branch_node_mt.clear();
        self.str_concordance = None;
        self.str_running = false;
        self.y_matches = None;
        self.y_matches_running = false;
        self.y_match_query.clear();
        self.y_profile_query.clear();
        self.mt_profile_query.clear();
        self.auto_profile_query.clear();
        self.private_y_query.clear();
        self.str_seq_query.clear();
        self.exchange_results.clear();
        self.coverage_by_aln.clear();
        self.audit_y.clear();
        self.audit_mt.clear();
        self.heteroplasmy = None;
        self.subject_brief = None;
        self.brief_narration = None;
        self.narration_stream = None;
        self.narrating = false;
        self.chat_history.clear();
        self.chat_input.clear();
        self.chat_pending = false;
        self.signal_narration.clear();
        self.signal_stream = None;
        self.signal_narrating = None;
        // The Simple-mode brief is (re)built from the Consensus event below (always fired by
        // LoadConsensus), so a later analysis refreshes it too. Show the spinner meanwhile.
        self.subject_brief_loading = self.ui_mode == UiMode::Simple;
        let _ = self.tx.send(Command::LoadConsensus(guid));
        // Refresh the list's Y/mt columns (picks up an assignment made on another row).
        let _ = self.tx.send(Command::LoadHaploSummary);
        let _ = self.tx.send(Command::LoadAudit {
            biosample_guid: guid,
            dna_type: DnaType::Y,
        });
        let _ = self.tx.send(Command::LoadAudit {
            biosample_guid: guid,
            dna_type: DnaType::Mt,
        });
        let _ = self.tx.send(Command::LoadRuns(guid));
        let _ = self.tx.send(Command::LoadStrProfiles(guid));
        let _ = self.tx.send(Command::LoadVariantSets(guid));
        let _ = self.tx.send(Command::LoadChipProfiles(guid));
        let _ = self.tx.send(Command::LoadMtdna(guid));
        // Subject-centric: auto-select the subject's default alignment so the analysis tabs work
        // without navigating Data Sources, and load the donor-level aggregates (best ancestry +
        // private-Y union across all sources).
        let _ = self.tx.send(Command::DefaultAlignment { biosample_guid: guid });
        let _ = self.tx.send(Command::LoadDonorAncestry { biosample_guid: guid });
        let _ = self
            .tx
            .send(Command::LoadConsensusAncestryDetail { biosample_guid: guid });
        // A cached chromosome painting (current for the consensus signature) shows without a click.
        let _ = self.tx.send(Command::LoadPainting { biosample_guid: guid });
        let _ = self.tx.send(Command::LoadDonorPrivateY { biosample_guid: guid });
        // The Y-variant profile is *built* on explicit request (re-genotypes each alignment), but a
        // previously-built snapshot loads cheaply — fetch it so the Y-DNA tab shows it immediately.
        let _ = self.tx.send(Command::LoadYProfile { biosample_guid: guid });
        // Likewise the mtDNA consensus profile (cheap cached snapshot for the mtDNA tab).
        let _ = self.tx.send(Command::LoadMtProfile { biosample_guid: guid });
        // The subject's persisted federated IBD exchange results (cheap; shown in the IBD tab).
        let _ = self.tx.send(Command::LoadIbdExchanges { biosample_guid: guid });
        // And the autosomal (diploid) consensus snapshot for the Autosomal tab.
        let _ = self.tx.send(Command::LoadAutosomalProfile { biosample_guid: guid });
        // Imported FTDNA genealogy (vendor ids + member labels + MDKA) for the Overview card.
        let _ = self.tx.send(Command::LoadGenealogy(guid));
    }

    pub(crate) fn select_run(&mut self, id: i64) {
        self.selected_run = Some(id);
        self.alignments.clear();
        self.selected_alignment = None;
        self.coverage = None;
        self.coverage_by_aln.clear();
        let _ = self.tx.send(Command::LoadAlignments(id));
    }

    pub(crate) fn select_alignment(&mut self, id: i64) {
        self.selected_alignment = Some(id);
        self.coverage = None;
        // Ideogram regions are fetched lazily when its tab opens; reset for the new alignment.
        self.genome_regions = None;
        self.loading_regions = false;
        self.regions_attempted = None;
        self.sex = None;
        self.read_metrics = None;
        self.sv = None;
        self.running_sex = false;
        self.running_metrics = false;
        self.running_sv = false;
        self.denovo.clear();
        self.ibd_result = None;
        self.identity = None;
        self.y_haplogroup = None;
        self.y_report = None;
        self.y_report_running = false;
        self.mt_haplogroup = None;
        self.private_y = None;
        let _ = self.tx.send(Command::LoadCoverage(id));
        let _ = self.tx.send(Command::LoadSex(id));
        let _ = self.tx.send(Command::LoadReadMetrics(id));
        let _ = self.tx.send(Command::LoadSv(id));
        // Load cached chrM de-novo (mtDNA tab). chrY variant discovery is the masked private-Y
        // pass, not a raw whole-chrY de-novo, so it isn't loaded here.
        let _ = self.tx.send(Command::LoadDenovo {
            alignment_id: id,
            contig: "chrM".into(),
        });
        let _ = self.tx.send(Command::LoadPrivateY { alignment_id: id }); // reload cached private-Y
    }

    fn clear_sample_selection(&mut self) {
        self.selected_sample = None;
        self.runs.clear();
        self.clear_run_selection();
    }

    fn clear_run_selection(&mut self) {
        self.selected_run = None;
        self.alignments.clear();
        self.selected_alignment = None;
        self.coverage = None;
    }

    /// Drop every in-flight spinner after a failure. A command failure is reported by whichever
    /// worker arm was running, but the UI has no way to tell which flag that arm owned, so all of
    /// them clear — a stuck spinner outlives the error message that explains it.
    fn clear_in_flight(&mut self) {
        self.cancelling = false;
        self.running = false;
        self.running_denovo = false;
        self.running_ibd = false;
        self.loading_ibd_suggestions = false;
        self.logging_in = false;
        self.publishing = false;
        self.finding_private_y = false;
        self.estimating_donor_ancestry = false;
        self.estimating_deep_ancestry = false;
        self.y_profile_loading = false;
        self.painting_running = false;
        self.running_sex = false;
        self.running_metrics = false;
        self.running_sv = false;
        self.loading_regions = false;
        let _ = self.tx.send(Command::SyncStatus); // a failed publish may have gone offline
    }
}
