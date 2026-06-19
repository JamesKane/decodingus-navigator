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
                    self.import_progress = None;
                    self.pending_import_dir = None;
                    self.reference_needs.clear();
                    self.select_project(summary.project.id);
                    let _ = self.tx.send(Command::LoadOverview);
                    let _ = self.tx.send(Command::LoadAllBiosamples);
                }
                Event::ReferenceNeeded { dir, builds } => {
                    self.importing = false;
                    self.import_progress = None;
                    self.pending_import_dir = Some(dir);
                    self.status = format!("{} reference build(s) need downloading", builds.len());
                    self.reference_needs = builds;
                }
                Event::ReferenceProgress { build, received, total } => {
                    self.reference_progress = Some((build, received, total));
                }
                Event::ImportProgress { done, total, label } => {
                    self.import_progress = Some((done, total, label));
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
                    self.deep_progress = None;
                    self.status = format!(
                        "{} {samples} sample(s): {coverage_done} coverage, {y_done} Y, {sex_done} sex, {metrics_done} metrics, {sv_done} SV{}",
                        if cancelled { "Deep analysis cancelled after" } else { "Analyzed" },
                        if errors > 0 { format!(", {errors} error(s)") } else { String::new() }
                    );
                    if self.selected_project == Some(project_id) {
                        let _ = self.tx.send(Command::LoadProjectReport(project_id));
                    }
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
                Event::AllBiosamples(v) => {
                    self.all_biosamples = v;
                    let _ = self.tx.send(Command::LoadHaploSummary); // fill the Y/mt columns
                }
                Event::HaploSummary(map) => self.haplo_summary = map,
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
                    if self.selected_sample == Some(biosample_guid) {
                        let _ = self.tx.send(Command::LoadRuns(biosample_guid));
                        let _ = self.tx.send(Command::LoadStrProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadVariantSets(biosample_guid));
                        let _ = self.tx.send(Command::LoadChipProfiles(biosample_guid));
                        let _ = self.tx.send(Command::LoadMtdna(biosample_guid));
                    }
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
                    nmonte,
                } => {
                    if self.selected_sample == Some(biosample_guid) {
                        self.fine_ancestry = fine.map(|b| *b);
                        self.ancient_ancestry = ancient.map(|b| *b);
                        self.nmonte_ancestry = nmonte.map(|b| *b);
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
                    self.status = if cancelled {
                        "Full analysis cancelled.".into()
                    } else {
                        "Full analysis complete.".into()
                    };
                }
                Event::Panels(p) => self.panels = p,
                Event::PanelImported => {
                    self.status = "Panel imported".into();
                    let _ = self.tx.send(Command::LoadPanels);
                }
                Event::AllAlignments(a) => self.all_alignments = a,
                Event::PanelGenotypes {
                    alignment_id,
                    panel_id,
                    ploidy,
                    genotypes,
                } => {
                    if self.selected_alignment == Some(alignment_id)
                        && self.selected_panel == Some(panel_id)
                        && self.ploidy() == ploidy
                    {
                        self.panel_genotypes = (!genotypes.is_empty()).then_some(genotypes);
                    }
                    self.running_genotype = false;
                }
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
                    self.status = if missing == 0 {
                        "All source files present".into()
                    } else {
                        format!("{missing} source file(s) moved or missing")
                    };
                }
                Event::Error(e) => {
                    self.status = format!("Error: {e}");
                    self.importing = false;
                    self.import_progress = None;
                    self.running = false;
                    self.running_denovo = false;
                    self.running_genotype = false;
                    self.running_ibd = false;
                    self.loading_ibd_suggestions = false;
                    self.logging_in = false;
                    self.publishing = false;
                    self.finding_private_y = false;
                    self.estimating_donor_ancestry = false;
                    self.y_profile_loading = false;
                    self.painting_running = false;
                    self.running_sex = false;
                    self.running_metrics = false;
                    self.running_sv = false;
                    self.loading_regions = false;
                    let _ = self.tx.send(Command::SyncStatus); // a failed publish may have gone offline
                }
            }
        }
    }

    pub(crate) fn select_project(&mut self, id: i64) {
        self.selected_project = Some(id);
        self.samples.clear();
        self.project_report.clear();
        self.clear_sample_selection();
        let _ = self.tx.send(Command::LoadSamples(id));
        let _ = self.tx.send(Command::LoadProjectReport(id));
    }

    pub(crate) fn select_sample(&mut self, guid: SampleGuid) {
        self.selected_sample = Some(guid);
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
        self.nmonte_ancestry = None;
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
        self.panel_genotypes = None;
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
        if let Some(panel_id) = self.selected_panel {
            let _ = self.tx.send(Command::LoadPanelGenotypes {
                alignment_id: id,
                panel_id,
                ploidy: self.ploidy(),
            });
        }
    }

    pub(crate) fn select_panel(&mut self, panel_id: i64) {
        self.selected_panel = Some(panel_id);
        self.panel_genotypes = None;
        self.ibd_result = None;
        self.identity = None;
        if let Some(aln) = self.selected_alignment {
            let _ = self.tx.send(Command::LoadPanelGenotypes {
                alignment_id: aln,
                panel_id,
                ploidy: self.ploidy(),
            });
        }
    }

    pub(crate) fn ploidy(&self) -> u8 {
        self.forms.ploidy.trim().parse().unwrap_or(2)
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
}
