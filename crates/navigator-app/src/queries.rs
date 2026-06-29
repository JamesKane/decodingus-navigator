//! `impl App` methods extracted from `lib.rs` (the `queries` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- queries -----------------------------------------------------------

    /// Biosamples belonging to a project (M:N membership ∪ legacy home column).
    pub async fn list_biosamples(&self, project_id: i64) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_members_for_project(self.store.pool(), project_id).await?)
    }

    /// Every biosample (subject), regardless of project association.
    pub async fn list_all_biosamples(&self) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_all(self.store.pool()).await?)
    }

    /// Bulk per-subject analysis status for the Subjects list, in one query (mirrors
    /// [`haplogroup_terminals`](Self::haplogroup_terminals)). A subject is `Complete` once every
    /// alignment it owns has a full `coverage` artifact at the current version; otherwise `Pending`.
    /// Subjects with no alignments are omitted (the list shows no status for them).
    pub async fn subject_analysis_status(
        &self,
    ) -> Result<HashMap<SampleGuid, SubjectAnalysisStatus>, AppError> {
        let census =
            artifact::analyzed_census(self.store.pool(), "coverage", coverage::COVERAGE_VERSION).await?;
        Ok(census
            .into_iter()
            .map(|(guid, total, analyzed)| {
                let status = if total > 0 && analyzed >= total {
                    SubjectAnalysisStatus::Complete
                } else {
                    SubjectAnalysisStatus::Pending
                };
                (guid, status)
            })
            .collect())
    }

    /// Sequence runs for a biosample.
    pub async fn list_sequence_runs(&self, biosample_guid: SampleGuid) -> Result<Vec<SequenceRun>, AppError> {
        let mut runs = sequence_run::list_for_biosample(self.store.pool(), biosample_guid).await?;
        // One-time backfill: runs analyzed before read stats were mirrored onto the run carry no
        // `total_reads` (and older imports no `library_layout`). Recover them from a cached
        // `read_metrics` artifact on any of the run's alignments and persist, so the card shows
        // library stats + PE/SE without a re-walk.
        for run in &mut runs {
            if run.total_reads.is_some() && run.library_layout.is_some() {
                continue;
            }
            let alns = alignment::list_for_run(self.store.pool(), run.id).await?;
            for a in &alns {
                if let Some(m) = self.cached_read_metrics(a.id).await? {
                    self.write_back_read_stats(a.id, &m).await?;
                    run.total_reads = Some(m.total_reads as i64);
                    run.mean_read_length = (m.mean_read_length > 0.0).then_some(m.mean_read_length);
                    run.mean_insert_size = (m.mean_insert_size > 0.0).then_some(m.mean_insert_size);
                    if m.pf_reads_aligned > 0 {
                        run.library_layout = Some(
                            if m.reads_aligned_in_pairs > 0 {
                                "PAIRED"
                            } else {
                                "SINGLE"
                            }
                            .into(),
                        );
                    }
                    break;
                }
            }
        }
        Ok(runs)
    }

    /// Cached coverage for several alignments at once (Data Sources alignment rows). `None` for any
    /// alignment without a persisted coverage artifact. No genotyping/walking — pure cache reads.
    pub async fn cached_coverage_bulk(
        &self,
        alignment_ids: &[i64],
    ) -> Result<Vec<(i64, Option<CoverageResult>)>, AppError> {
        let mut out = Vec::with_capacity(alignment_ids.len());
        for &id in alignment_ids {
            out.push((id, self.cached_coverage(id).await?));
        }
        Ok(out)
    }

    /// Alignments for a sequence run.
    /// The best alignment to drive a subject's analysis tabs (subject-centric default): the
    /// highest mean-coverage alignment with a cached coverage result, else the first with a BAM,
    /// else the first. Returns `(sequence_run_id, alignment_id)` so the UI can select the run then
    /// the alignment without the user navigating Data Sources.
    pub async fn default_alignment_for_subject(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<(i64, i64)>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        if alignments.is_empty() {
            return Ok(None);
        }
        let mut best: Option<(f64, &Alignment)> = None;
        for a in &alignments {
            if let Some(c) = self.cached_coverage(a.id).await? {
                if best.as_ref().map_or(true, |(cov, _)| c.mean_coverage > *cov) {
                    best = Some((c.mean_coverage, a));
                }
            }
        }
        let chosen = best
            .map(|(_, a)| a)
            .or_else(|| alignments.iter().find(|a| a.bam_path.is_some()))
            .or_else(|| alignments.first());
        Ok(chosen.map(|a| (a.sequence_run_id, a.id)))
    }

    /// Donor-level ancestry: the **consensus** estimate ([`CONSENSUS_SOURCE_ID`]) when present —
    /// it pools all sources, so it's authoritative — else the best-quality per-alignment estimate
    /// (most genotyped SNPs) for back-compat with results predating the consensus path.
    pub async fn donor_ancestry(&self, biosample_guid: SampleGuid) -> Result<Option<(i64, AncestryResult)>, AppError> {
        let all = ancestry_result::for_biosample(self.store.pool(), biosample_guid).await?;
        if let Some(c) = all.iter().find(|(id, _)| *id == CONSENSUS_SOURCE_ID) {
            return Ok(Some(c.clone()));
        }
        Ok(all.into_iter().max_by_key(|(_, r)| r.snps_with_genotype))
    }

    /// A specific persisted consensus ancestry estimate (keyed on the consensus pseudo-source +
    /// `method`) — e.g. `"FINE_ADMIXTURE"` (detailed modern populations) or `"PCA_PROJECTION_GMM"`
    /// (ancient components). Filtered per-subject (alignment_id 0 isn't biosample-unique on its own).
    pub async fn consensus_ancestry(
        &self,
        biosample_guid: SampleGuid,
        method: &str,
    ) -> Result<Option<AncestryResult>, AppError> {
        let all = ancestry_result::for_biosample(self.store.pool(), biosample_guid).await?;
        Ok(all
            .into_iter()
            .find(|(id, r)| *id == CONSENSUS_SOURCE_ID && r.method == method)
            .map(|(_, r)| r))
    }

    /// Donor-level private-Y: the **union** of cached (self-masked) private-Y calls across all of
    /// the subject's alignments, deduped by position (keeping the deepest observation). The
    /// terminal is taken from the deepest-covered source bucket.
    pub async fn donor_private_y(&self, biosample_guid: SampleGuid) -> Result<Option<PrivateBucket>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut by_pos: std::collections::HashMap<i64, PrivateVariant> = std::collections::HashMap::new();
        let mut terminal: Option<String> = None;
        let mut any = false;
        for a in &alignments {
            let Some(bucket) = self.cached_private_y(a.id).await? else {
                continue;
            };
            any = true;
            terminal.get_or_insert_with(|| bucket.terminal.clone());
            for v in bucket.variants {
                by_pos
                    .entry(v.position)
                    .and_modify(|cur| {
                        if v.depth > cur.depth {
                            *cur = v.clone();
                        }
                    })
                    .or_insert(v);
            }
        }
        if !any {
            return Ok(None);
        }
        let mut variants: Vec<PrivateVariant> = by_pos.into_values().collect();
        variants.sort_by_key(|v| v.position);
        Ok(Some(PrivateBucket {
            terminal: terminal.unwrap_or_default(),
            variants,
        }))
    }

    pub async fn list_alignments(&self, sequence_run_id: i64) -> Result<Vec<Alignment>, AppError> {
        Ok(alignment::list_for_run(self.store.pool(), sequence_run_id).await?)
    }

    /// Every alignment in the workspace (for cross-sample selection like IBD compare).
    pub async fn list_all_alignments(&self) -> Result<Vec<Alignment>, AppError> {
        Ok(alignment::list_all(self.store.pool()).await?)
    }

    /// Projects with their sample counts, for a dashboard/list view.
    pub async fn project_overview(&self) -> Result<Vec<ProjectOverview>, AppError> {
        let mut out = Vec::new();
        for project in project::list(self.store.pool()).await? {
            let sample_count = biosample::count_members_for_project(self.store.pool(), project.id).await?;
            out.push(ProjectOverview { project, sample_count });
        }
        Ok(out)
    }

    /// Per-sample report for a project: each biosample's alignment count, coverage roll-up
    /// (the first alignment with cached coverage), and Y/mtDNA haplogroup consensus.
    /// Composes existing per-subject queries (no new join) — coverage/haplogroup cells are
    /// `None` until those analyses have run.
    pub async fn project_report(&self, project_id: i64) -> Result<Vec<ProjectSampleReport>, AppError> {
        let mut out = Vec::new();
        for biosample in biosample::list_members_for_project(self.store.pool(), project_id).await? {
            let alignments = alignment::list_for_biosample(self.store.pool(), biosample.guid).await?;
            let mut coverage = None;
            let mut coverage_aln = None;
            for a in &alignments {
                if let Some(c) = self.cached_coverage(a.id).await? {
                    coverage = Some(c);
                    coverage_aln = Some(a.id);
                    break;
                }
            }
            // A lite (sidecar) coverage is flagged so the UI can badge it and offer a deep walk.
            let coverage_partial = match coverage_aln {
                Some(id) => matches!(
                    self.analysis_provenance(id, "coverage", coverage::COVERAGE_VERSION).await?,
                    Some((_, ref c)) if c == "partial"
                ),
                None => false,
            };
            // Prefer the coverage-bearing alignment; else fall back to the first.
            let primary_alignment_id = coverage_aln.or_else(|| alignments.first().map(|a| a.id));
            let y_haplogroup = self
                .haplogroup_consensus(biosample.guid, DnaType::Y)
                .await?
                .map(|c| c.haplogroup);
            let mt_haplogroup = self
                .haplogroup_consensus(biosample.guid, DnaType::Mt)
                .await?
                .map(|c| c.haplogroup);
            // Sex + read-metrics from whichever alignment has them cached.
            let mut sex = None;
            let mut metrics = None;
            let mut sv_count = None;
            for a in &alignments {
                if sex.is_none() {
                    sex = self.cached_sex(a.id).await?;
                }
                if metrics.is_none() {
                    metrics = self.cached_read_metrics(a.id).await?;
                }
                if sv_count.is_none() {
                    sv_count = self.cached_sv(a.id).await?.map(|s| s.sv_calls.len());
                }
            }
            let sex = sex.map(|s| match s.inferred_sex {
                navigator_analysis::sex::InferredSex::Male => "M".to_string(),
                navigator_analysis::sex::InferredSex::Female => "F".to_string(),
                navigator_analysis::sex::InferredSex::Unknown => "U".to_string(),
            });
            out.push(ProjectSampleReport {
                primary_alignment_id,
                alignment_count: alignments.len(),
                mean_coverage: coverage.as_ref().map(|c| c.mean_coverage),
                median_coverage: coverage.as_ref().map(|c| c.median_coverage),
                pct_10x: coverage.as_ref().map(|c| c.pct_10x),
                pct_20x: coverage.as_ref().map(|c| c.pct_20x),
                callable_bases: coverage.as_ref().map(|c| c.callable_bases),
                y_haplogroup,
                mt_haplogroup,
                sex,
                mean_read_length: metrics.as_ref().map(|m| m.mean_read_length),
                pct_aligned: metrics.as_ref().map(|m| m.pct_pf_reads_aligned),
                median_insert_size: metrics.as_ref().map(|m| m.median_insert_size),
                sv_count,
                coverage_partial,
                biosample,
            });
        }
        Ok(out)
    }

    /// Per-member Y-STR overview for a project (the FTDNA-style "Y-DNA Results Overview"): each
    /// member that has at least one STR profile, with identity columns, terminal Y haplogroup, the
    /// reached STR panel/tier, and the consensus marker values (uppercase marker → value). Members
    /// with no STR data are omitted. Composes existing per-subject queries (no new join).
    pub async fn project_str_overview(&self, project_id: i64) -> Result<Vec<ProjectStrMember>, AppError> {
        use navigator_domain::{strpanel, strprofile};
        let mut out = Vec::new();
        for biosample in biosample::list_members_for_project(self.store.pool(), project_id).await? {
            let profiles = self.list_str_profiles(biosample.guid).await?;
            if profiles.is_empty() {
                continue;
            }
            // Consensus marker map, keyed by normalized (uppercase) marker name.
            let mut markers = std::collections::HashMap::new();
            for cm in strprofile::consensus_markers(&profiles) {
                if !cm.value.trim().is_empty() && cm.value.trim() != "-" {
                    markers.insert(strpanel::norm(&cm.marker), cm.value);
                }
            }
            // Reached panel/tier across all of the subject's markers (the "Test" column).
            let all_markers: Vec<navigator_domain::strprofile::StrMarker> =
                profiles.iter().flat_map(|p| p.markers.clone()).collect();
            let normed = strpanel::normalized_set(&all_markers);
            let provider = profiles
                .iter()
                .find_map(|p| p.provider.as_deref().filter(|s| !s.is_empty()))
                .unwrap_or("FTDNA");
            let test = strpanel::classify_panel(&normed, Some(provider)).panel_name;

            let consensus = self.haplogroup_consensus(biosample.guid, DnaType::Y).await?;
            let y_confirmed = consensus.is_some();
            let y_haplogroup = consensus.map(|c| c.haplogroup);

            out.push(ProjectStrMember {
                guid: biosample.guid,
                name: biosample.donor_identifier.clone(),
                kit: biosample.sample_accession.clone(),
                origin: biosample.center_name.clone(),
                ancestor: biosample.description.clone(),
                y_haplogroup,
                y_confirmed,
                test,
                markers,
            });
        }
        Ok(out)
    }

    /// Build the precomputed FTDNA-style Y-STR overview for a project: members grouped by their
    /// **assigned** (consensus) Y haplogroup, ordered by tree topology (basal → derived, children
    /// nested under their ancestor subgroups), with per-subgroup MIN/MAX/MODE and per-cell deviation
    /// from the modal value precomputed. Members without a SNP haplogroup fall into an "Unassigned"
    /// bucket at the base. All heavy work happens here (off the UI thread); the renderer just
    /// iterates [`ProjectStrChart::rows`].
    pub async fn project_str_chart(&self, project_id: i64) -> Result<ProjectStrChart, AppError> {
        use navigator_domain::{strchart, strpanel};
        use std::collections::{BTreeMap, HashMap, HashSet};

        let members = self.project_str_overview(project_id).await?;
        if members.is_empty() {
            return Ok(ProjectStrChart::default());
        }

        // Marker columns: canonical FTDNA order restricted to markers anyone reported, then extras.
        let mut present: HashSet<String> = HashSet::new();
        for m in &members {
            present.extend(m.markers.keys().cloned());
        }
        let mut markers: Vec<String> = Vec::new();
        for name in strpanel::ftdna_marker_order() {
            let n = strpanel::norm(name);
            if present.remove(&n) {
                markers.push(n);
            }
        }
        let mut extras: Vec<String> = present.into_iter().collect();
        extras.sort();
        markers.extend(extras);

        // Group members by assigned Y haplogroup; the unplaced share a bucket keyed by None.
        let mut groups: HashMap<Option<String>, Vec<&ProjectStrMember>> = HashMap::new();
        for m in &members {
            groups.entry(m.y_haplogroup.clone()).or_default().push(m);
        }

        // Tree topology for ordering (best-effort; alphabetical fallback when unavailable).
        let tree = self.chip_y_tree("GRCh38").await.ok();
        let (preorder, name_idx, parents) = match &tree {
            Some(t) => {
                let names = tree_name_index(t);
                (tree_preorder(t), names, tree_parent_map(t))
            }
            None => (HashMap::new(), HashMap::new(), HashMap::new()),
        };
        let group_keys: HashSet<String> = groups.keys().flatten().map(|h| norm_hg(h)).collect();
        let node_of = |hg: &str| -> Option<i64> { name_idx.get(&norm_hg(hg)).copied() };

        // Order the placed groups by tree pre-order (basal → derived); unmatched names sort after,
        // alphabetically. Depth = count of ancestor haplogroups that are themselves groups here.
        let mut placed: Vec<(String, usize, i64)> = Vec::new(); // (haplogroup, depth, sort_rank)
        for hg in groups.keys().flatten() {
            let depth = match node_of(hg) {
                Some(id) => {
                    let mut d = 0usize;
                    let mut cur = parents.get(&id).copied();
                    while let Some(p) = cur {
                        if let Some(node) = tree.as_ref().and_then(|t| t.nodes.get(&p)) {
                            if group_keys.contains(&norm_hg(&node.name)) {
                                d += 1;
                            }
                        }
                        cur = parents.get(&p).copied();
                    }
                    d
                }
                None => 0,
            };
            let rank = node_of(hg)
                .and_then(|id| preorder.get(&id))
                .map(|r| *r as i64)
                .unwrap_or(i64::MAX);
            placed.push((hg.clone(), depth, rank));
        }
        // Sort: matched (rank < MAX) by pre-order; unmatched by name.
        placed.sort_by(|a, b| a.2.cmp(&b.2).then_with(|| a.0.cmp(&b.0)));

        let dev_cells = |m: &ProjectStrMember, stats: &BTreeMap<&str, strchart::MarkerStats>| -> Vec<StrChartCell> {
            markers
                .iter()
                .map(|c| match m.markers.get(c) {
                    Some(v) => StrChartCell {
                        dev: strchart::deviation(v, &stats[c.as_str()].mode),
                        text: v.clone(),
                    },
                    None => StrChartCell {
                        text: String::new(),
                        dev: strchart::Deviation::None,
                    },
                })
                .collect()
        };

        let mut rows: Vec<StrChartRow> = Vec::new();

        // Emit one subgroup: banner + MIN/MAX/MODE + members (members sorted by name).
        let emit_group =
            |rows: &mut Vec<StrChartRow>, label: String, depth: usize, mut members: Vec<&ProjectStrMember>| {
                members.sort_by(|a, b| a.name.cmp(&b.name));
                let mut stats: BTreeMap<&str, strchart::MarkerStats> = BTreeMap::new();
                for c in &markers {
                    let vals = members.iter().filter_map(|m| m.markers.get(c).map(String::as_str));
                    stats.insert(c.as_str(), strchart::marker_stats(vals));
                }
                rows.push(StrChartRow {
                    kind: StrRowKind::Group,
                    depth,
                    label: format!("{label}  ({})", members.len()),
                    kit: String::new(),
                    haplogroup: String::new(),
                    confirmed: false,
                    test: String::new(),
                    cells: Vec::new(),
                });
                for (kind, pick) in [(StrRowKind::Min, 0u8), (StrRowKind::Max, 1u8), (StrRowKind::Mode, 2u8)] {
                    let cells = markers
                        .iter()
                        .map(|c| {
                            let s = &stats[c.as_str()];
                            let t = match pick {
                                0 => s.min.clone(),
                                1 => s.max.clone(),
                                _ => s.mode.clone(),
                            };
                            StrChartCell {
                                text: t.unwrap_or_default(),
                                dev: strchart::Deviation::None,
                            }
                        })
                        .collect();
                    rows.push(StrChartRow {
                        kind,
                        depth,
                        label: String::new(),
                        kit: String::new(),
                        haplogroup: String::new(),
                        confirmed: false,
                        test: String::new(),
                        cells,
                    });
                }
                for m in &members {
                    rows.push(StrChartRow {
                        kind: StrRowKind::Member,
                        depth,
                        label: m.name.clone(),
                        kit: m.kit.clone().unwrap_or_default(),
                        haplogroup: m.y_haplogroup.clone().unwrap_or_default(),
                        confirmed: m.y_confirmed,
                        test: m.test.clone().unwrap_or_default(),
                        cells: dev_cells(m, &stats),
                    });
                }
            };

        // Unassigned bucket first (the base), then the placed clades in tree order.
        if let Some(unplaced) = groups.get(&None) {
            emit_group(&mut rows, "Unassigned".to_string(), 0, unplaced.clone());
        }
        for (hg, depth, _) in &placed {
            if let Some(ms) = groups.get(&Some(hg.clone())) {
                emit_group(&mut rows, hg.clone(), *depth, ms.clone());
            }
        }

        Ok(ProjectStrChart {
            markers,
            rows,
            member_count: members.len(),
            group_count: placed.len() + usize::from(groups.contains_key(&None)),
        })
    }

    /// Analyze every sample in a project: compute coverage and assign the Y haplogroup on each
    /// sample's primary (first BAM-bearing) alignment, so the project report fills in. Coverage
    /// already cached and Y already recorded are skipped (idempotent re-run). Best-effort: one
    /// sample's failure is recorded and the rest continue. mtDNA is intentionally not assigned
    /// here (provisional on CHM13 — see the reconciliation/liftover notes).
    pub async fn analyze_project(&self, project_id: i64) -> Result<AnalyzeSummary, AppError> {
        let mut summary = AnalyzeSummary {
            project_id,
            samples: 0,
            coverage_done: 0,
            y_done: 0,
            sex_done: 0,
            metrics_done: 0,
            sv_done: 0,
            errors: Vec::new(),
        };
        for biosample in biosample::list_members_for_project(self.store.pool(), project_id).await? {
            let o = self.analyze_biosample(&biosample).await?;
            if !o.had_alignment {
                continue;
            }
            summary.samples += 1;
            summary.coverage_done += o.coverage_done as usize;
            summary.y_done += o.y_done as usize;
            summary.sex_done += o.sex_done as usize;
            summary.metrics_done += o.metrics_done as usize;
            summary.sv_done += o.sv_done as usize;
            summary.errors.extend(o.errors);
        }
        Ok(summary)
    }

    /// Deep-analyze one biosample's primary (first BAM-bearing) alignment: coverage, Y
    /// haplogroup, sex, read metrics, and SV (≥10× only). Idempotent — a *full* coverage and a
    /// recorded Y/sex/metrics/SV are skipped; a `partial` (lite sidecar) coverage is upgraded by
    /// the per-base walk, which overwrites it. Best-effort: a per-step failure is recorded in
    /// `errors` (prefixed with the donor id) and the remaining steps still run. This is the
    /// per-sample unit the project pass and the streaming deep-analyze job both drive.
    pub async fn analyze_biosample(&self, biosample: &Biosample) -> Result<SampleAnalyzeOutcome, AppError> {
        let mut o = SampleAnalyzeOutcome::default();
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample.guid).await?;
        let Some(aln) = alignments.iter().find(|a| a.bam_path.is_some()) else {
            return Ok(o); // had_alignment stays false
        };
        o.had_alignment = true;
        let label = &biosample.donor_identifier;
        // Drop any prior subject's local alignment copy so the batch holds at most one file's worth
        // of cache; this subject's passes share the single copy `localize` makes below.
        Self::clear_align_cache();

        // Coverage + read-metrics + sex in ONE pass (the unified walker) instead of three separate
        // reads of the BAM/CRAM — a 3x I/O cut per subject, which dominates the batch on a slow /
        // network volume (the single-subject Full Analysis already does this; the batch path didn't).
        // Walk only when something's missing: a full, correctly-scoped coverage (a stale whole-genome
        // result for a targeted-Y test is recomputed) plus cached read-metrics and sex = all done.
        let coverage_full = matches!(
            self.analysis_provenance(aln.id, "coverage", coverage::COVERAGE_VERSION).await?,
            Some((_, ref c)) if c == "full"
        ) && match self.cached_coverage(aln.id).await? {
            Some(cov) => self.coverage_is_correctly_scoped(aln.id, &cov).await?,
            None => false,
        };
        if coverage_full
            && self.cached_read_metrics(aln.id).await?.is_some()
            && self.cached_sex(aln.id).await?.is_some()
        {
            o.coverage_done = true;
            o.metrics_done = true;
            o.sex_done = true;
        } else {
            match self.run_unified_metrics(aln.id).await {
                Ok(_) => {
                    o.coverage_done = true;
                    o.metrics_done = true;
                    o.sex_done = true;
                }
                Err(e) => o.errors.push(format!("{label} metrics: {e}")),
            }
        }

        if self.haplogroup_consensus(biosample.guid, DnaType::Y).await?.is_some() {
            o.y_done = true;
        } else {
            match self.assign_y_haplogroup(aln.id).await {
                Ok(_) => o.y_done = true,
                Err(e) => o.errors.push(format!("{label} Y: {e}")),
            }
        }

        // SV is a whole-genome analysis that walks every read in the file. Skip it for targeted
        // tests (Big Y / Y Elite / mtFull): SV is meaningless there, and over a targeted CRAM's
        // millions of off-target reads the whole-file walk is pathologically slow. Gate on the
        // run's target being whole-genome, plus the existing ≥10× depth threshold (avoids logging a
        // "coverage too low" error for every low-coverage sample).
        let is_wgs = match sequence_run::get(self.store.pool(), aln.sequence_run_id).await? {
            Some(run) => matches!(
                navigator_domain::testtype::target_of(&run.test_type),
                Some(navigator_domain::testtype::TargetType::WholeGenome)
            ),
            None => false,
        };
        if self.cached_sv(aln.id).await?.is_some() {
            o.sv_done = true;
        } else if is_wgs
            && self
                .cached_coverage(aln.id)
                .await?
                .map(|c| c.mean_coverage >= 10.0)
                .unwrap_or(false)
        {
            match self.run_sv(aln.id).await {
                Ok(_) => o.sv_done = true,
                Err(e) => o.errors.push(format!("{label} SV: {e}")),
            }
        }
        Ok(o)
    }
}

// ---- Y-tree topology helpers for the project STR chart ordering --------------------------------

/// Normalize a haplogroup / node name for matching: uppercase, and (since our consensus labels and
/// the tree nodes both use the "R-CTS4466" convention) keep the full name. Tolerates a bare SNP by
/// also being comparable to the suffix after the last '-' (callers index both forms).
fn norm_hg(name: &str) -> String {
    name.trim().to_ascii_uppercase()
}

/// Map every tree node name (and its bare-SNP suffix) to a node id, for resolving haplogroup labels.
/// Full names win over suffix aliases on collision.
fn tree_name_index(tree: &navigator_analysis::haplo::HaploTree) -> std::collections::HashMap<String, i64> {
    let mut idx = std::collections::HashMap::new();
    // Pass 1: suffix aliases (lower priority).
    for (id, node) in &tree.nodes {
        if let Some(suffix) = node.name.rsplit('-').next() {
            idx.entry(norm_hg(suffix)).or_insert(*id);
        }
    }
    // Pass 2: full names (override).
    for (id, node) in &tree.nodes {
        idx.insert(norm_hg(&node.name), *id);
    }
    idx
}

/// child → parent map over the tree.
fn tree_parent_map(tree: &navigator_analysis::haplo::HaploTree) -> std::collections::HashMap<i64, i64> {
    let mut parent = std::collections::HashMap::new();
    for (id, node) in &tree.nodes {
        for c in &node.children {
            parent.insert(*c, *id);
        }
    }
    parent
}

/// Pre-order DFS rank for every node (basal → derived; children follow their parent, siblings in
/// stored order) so groups can be ordered to mirror the tree.
fn tree_preorder(tree: &navigator_analysis::haplo::HaploTree) -> std::collections::HashMap<i64, usize> {
    let mut rank = std::collections::HashMap::new();
    let mut next = 0usize;
    // Roots: explicit is_root flag, else nodes with no parent.
    let parents = tree_parent_map(tree);
    let mut roots: Vec<i64> = tree
        .nodes
        .values()
        .filter(|n| n.is_root || !parents.contains_key(&n.id))
        .map(|n| n.id)
        .collect();
    roots.sort_unstable();
    let mut stack: Vec<i64> = roots.into_iter().rev().collect();
    while let Some(id) = stack.pop() {
        if rank.contains_key(&id) {
            continue;
        }
        rank.insert(id, next);
        next += 1;
        if let Some(node) = tree.nodes.get(&id) {
            for c in node.children.iter().rev() {
                stack.push(*c);
            }
        }
    }
    rank
}
