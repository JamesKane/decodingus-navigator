//! `impl App` methods extracted from `lib.rs` (the `queries` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// Each analysis artifact of a set of alignments. The code reads them in advance and indexes them
/// for the report builders.
///
/// A read of one cached result through [`App::load_analysis`] costs two queries and one stat call.
/// The queries read the artifact and then the `alignment` row, and the code needs that row to stat
/// the BAM file.
///
/// A project report reads five kinds of artifact for each alignment of each member. So the earlier
/// form, one read for each cell, sent thousands of queries to open one tab. This type reads each
/// artifact in one `IN` query and stats each BAM file one time.
///
/// The rule for an old result is the rule of [`App::load_analysis`]. A cached payload is absent when
/// the `mtime:size` value of the source file changed after the calculation.
struct AlignmentArtifacts {
    /// `(alignment id, kind, algorithm version)` → the stored artifact.
    by_key: HashMap<(i64, String, String), AnalysisArtifact>,
    /// The current source signature of each alignment. A value of `None` means that the file is
    /// absent, or that the operating system can not read its metadata. The code then trusts the
    /// cache, because it has no value to compare.
    sigs: HashMap<i64, Option<String>>,
}

impl AlignmentArtifacts {
    async fn load(store: &Store, alignments: &[&Alignment]) -> Result<Self, AppError> {
        let ids: Vec<i64> = alignments.iter().map(|a| a.id).collect();
        let by_key = artifact::list_for_alignments(store.pool(), &ids)
            .await?
            .into_iter()
            .map(|a| ((a.alignment_id, a.kind.clone(), a.algorithm_version.clone()), a))
            .collect();
        // One stat call for each alignment, from the row that the code already holds. The code
        // does not make one stat call for each artifact.
        let sigs = alignments
            .iter()
            .map(|a| {
                let sig = a.bam_path.as_deref().and_then(|p| file_signature(Path::new(p)));
                (a.id, sig)
            })
            .collect();
        Ok(Self { by_key, sigs })
    }

    /// The stored artifact, with no check for an old result. The method does the same work as a
    /// plain `artifact::get` call.
    fn raw(&self, alignment_id: i64, kind: &str, version: &str) -> Option<&AnalysisArtifact> {
        self.by_key.get(&(alignment_id, kind.to_string(), version.to_string()))
    }

    /// The decoded payload. The method returns `None` when the artifact is absent, when it is out
    /// of date, or when the decoder refuses it. [`App::load_analysis`] behaves in the same way.
    fn fresh<T: DeserializeOwned>(&self, alignment_id: i64, kind: &str, version: &str) -> Option<T> {
        let a = self.raw(alignment_id, kind, version)?;
        let current = self.sigs.get(&alignment_id).and_then(|s| s.as_deref());
        if !artifact_is_fresh(a.source_sig.as_deref(), current) {
            return None;
        }
        serde_json::from_str(&a.payload).ok()
    }

    /// `(source, completeness)` with the same legacy-row defaults as [`App::analysis_provenance`].
    fn provenance(&self, alignment_id: i64, kind: &str, version: &str) -> Option<(String, String)> {
        self.raw(alignment_id, kind, version).map(|a| {
            (
                a.source.clone().unwrap_or_else(|| "navigator-walk".into()),
                a.completeness.clone().unwrap_or_else(|| "full".into()),
            )
        })
    }
}

impl App {
    // ---- queries -----------------------------------------------------------

    /// The biosamples of a project. The set is the union of the M:N memberships and the old home
    /// column.
    pub async fn list_biosamples(&self, project_id: i64) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_members_for_project(self.store.pool(), project_id).await?)
    }

    /// Every biosample (subject), regardless of project association.
    pub async fn list_all_biosamples(&self) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_all(self.store.pool()).await?)
    }

    /// The analysis status of each subject for the Subjects list, in one query. The method has the
    /// same shape as [`haplogroup_terminals`](Self::haplogroup_terminals).
    ///
    /// A subject is `Complete` when each of its alignments has a full `coverage` artifact at the
    /// current version. If not, the subject is `Pending`.
    ///
    /// The result holds no subject with no alignment, and the list then shows no status for such a
    /// subject.
    pub async fn subject_analysis_status(&self) -> Result<HashMap<SampleGuid, SubjectAnalysisStatus>, AppError> {
        let census = artifact::analyzed_census(self.store.pool(), "coverage", coverage::COVERAGE_VERSION).await?;
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
        // A backfill that runs one time. A run that the app analyzed before it copied the read
        // statistics to the run row holds no `total_reads` value. An older import also holds no
        // `library_layout` value.
        //
        // The code reads those values from a cached `read_metrics` artifact on any alignment of the
        // run, and writes them to the run row. The card then shows the library statistics and the
        // PE or SE value with no second walk.
        for run in &mut runs {
            if run.total_reads.is_some() && run.library_layout.is_some() && run.total_bases.is_some() {
                continue;
            }
            let alns = alignment::list_for_run(self.store.pool(), run.id).await?;
            for a in &alns {
                if let Some(m) = self.cached_read_metrics(a.id).await? {
                    self.write_back_read_stats(a.id, &m).await?;
                    run.total_reads = Some(m.total_reads as i64);
                    run.total_bases = m.total_bases().or(run.total_bases);
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

    /// Fill the read-profile fields of many runs at once. Those fields support the standard test
    /// label in [`du_domain::testprofile`]. The CLI command `backfill-profiles` calls this method
    /// for a run that the app imported before the fields existed. A second call is safe, and the
    /// method fills only an empty field.
    ///
    /// - **`total_bases`** comes from a cached `read_metrics` artifact on any alignment of the run.
    ///   The value is the sum of `read_length_histogram`, and the method walks no file.
    /// - **`read_type`** comes from `platform_name` and `test_type` at a low cost. The values are
    ///   `SHORT`, `ONT_SIMPLEX`, and `HIFI` or `CLR` when the code already holds one of them.
    ///
    ///   If those fields give nothing, the method reads the cached mean read length. A short mean
    ///   gives `SHORT`, and that rule resolves a run from a sidecar import with an `UNKNOWN`
    ///   platform.
    ///
    ///   With `rescan`, a run that still has no value gets a limited
    ///   [`library_stats`](Self::library_stats) scan of one alignment file. A long-read run needs
    ///   that scan, because only the read names separate HiFi from CLR.
    ///
    /// `project_id` limits the work to the subjects of one project. The method reads the old home
    /// column, as `rebuild_signatures` does. It returns a count for each field.
    pub async fn backfill_read_profiles(
        &self,
        project_id: Option<i64>,
        rescan: bool,
    ) -> Result<ReadProfileBackfill, AppError> {
        let mut out = ReadProfileBackfill::default();
        for b in self.list_all_biosamples().await? {
            if let Some(pid) = project_id {
                if b.project_id != Some(pid) {
                    continue;
                }
            }
            let runs = sequence_run::list_for_biosample(self.store.pool(), b.guid).await?;
            for run in &runs {
                out.runs_examined += 1;
                let alns = alignment::list_for_run(self.store.pool(), run.id).await?;

                // One cached read-metrics artifact stands for the run. The code uses it for the
                // yield and for the read-length value below, and it reads the artifact one time.
                let mut metrics = None;
                for a in &alns {
                    if let Some(m) = self.cached_read_metrics(a.id).await? {
                        metrics = Some((a.id, m));
                        break;
                    }
                }

                // total_bases: recovered from the cached metrics (write_back persists it).
                if run.total_bases.is_none() {
                    if let Some((aid, m)) = &metrics {
                        if m.total_bases().is_some() {
                            self.write_back_read_stats(*aid, m).await?;
                            out.total_bases_filled += 1;
                        }
                    }
                }

                // The `read_type` value. The code first reads the platform and the test type,
                // which costs little. It then reads the cached mean read length, and that value
                // resolves a run from a sidecar import with an UNKNOWN platform. It can then scan
                // the file. A long-read run needs that scan, because only the read names separate
                // HiFi from CLR.
                if run.read_type.is_none() {
                    let inferred = infer_read_type_cheap(&run.platform_name, &run.test_type).or_else(|| {
                        metrics
                            .as_ref()
                            .and_then(|(_, m)| read_type_from_mean_len(m.mean_read_length))
                    });
                    match inferred {
                        Some(rt) => {
                            sequence_run::set_read_type(self.store.pool(), run.id, rt).await?;
                            out.read_type_filled += 1;
                        }
                        None if rescan => match self.rescan_read_type(&alns).await {
                            Some(rt) => {
                                sequence_run::set_read_type(self.store.pool(), run.id, &rt).await?;
                                out.read_type_rescanned += 1;
                            }
                            None => out.read_type_unresolved += 1,
                        },
                        None => out.read_type_unresolved += 1,
                    }
                }
            }
        }
        Ok(out)
    }

    /// A limited library-stats scan of the first alignment file in `alns` that the code can read.
    /// The method returns the `read_type` value that it deduces, and the read names separate HiFi
    /// from CLR.
    ///
    /// The method returns `None` when it can open no file, and when it finds nothing that it can
    /// decode.
    async fn rescan_read_type(&self, alns: &[navigator_domain::workspace::Alignment]) -> Option<String> {
        for a in alns {
            if a.bam_path.is_none() {
                continue;
            }
            // Find the reference for the decoder. See alignment_reference_for_decode. A CRAM file
            // needs it, and a BAM file uses None. The step is optional, and the code skips an
            // alignment with no reference.
            let Ok((path, reference)) = self.alignment_reference_for_decode(a.id).await else {
                continue;
            };
            if !path.exists() {
                continue;
            }
            if let Ok(stats) = self.library_stats(path, reference).await {
                if stats.read_type.is_some() {
                    return stats.read_type;
                }
            }
        }
        None
    }

    /// The cached coverage of many alignments at once, for the alignment rows of the Data Sources
    /// tab. The value is `None` for an alignment with no stored coverage artifact. The method reads
    /// the cache only. It calls no caller and walks no file.
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

    /// The alignments of a sequence run.
    ///
    /// The method returns the alignment that drives the analysis tabs of a subject. It takes the
    /// alignment with the highest mean coverage that also has a cached coverage result. If there is
    /// none, it takes the first alignment with a BAM file, and then the first alignment.
    ///
    /// The method returns `(sequence_run_id, alignment_id)`. The UI can then select the run and the
    /// alignment, and the user does not open the Data Sources tab.
    ///
    /// When a realignment exists, the method takes its output before the source of that
    /// realignment. This rule applies only when the breadth and the depth are equal. See the
    /// comment on the order below.
    pub async fn default_alignment_for_subject(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<(i64, i64)>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        if alignments.is_empty() {
            return Ok(None);
        }
        // The order is the breadth, then the depth, then the presence of a file.
        //
        // A whole-genome test, such as WGS or HiFi, represents the subject before a targeted Y test
        // or mt test. This rule applies **even when the targeted test has more depth**.
        //
        // The mean depth of a Y-only test covers chrY alone, because the coverage covers the target
        // contigs only. So an order on the depth alone puts a deep Y Elite test above a WGS test.
        // The app then names the Y test as "your test", and that name disagrees with the autosomal
        // ancestry that the brief shows beside it.
        //
        // The depth and the presence of a file apply only inside one breadth class.
        let mut best: Option<(u8, f64, bool, bool, &Alignment)> = None;
        for a in &alignments {
            let target = match navigator_store::sequence_run::get(self.store.pool(), a.sequence_run_id).await? {
                Some(run) => navigator_domain::testtype::target_of(&run.test_type),
                None => None,
            };
            let breadth = test_breadth_rank(target);
            let depth = self.cached_coverage(a.id).await?.map_or(0.0, |c| c.mean_coverage);
            // The last rule, and it applies only when the values above are equal. A realigned
            // alignment goes before the source that made it.
            //
            // Both rows describe the same library, at the same breadth and at almost the same
            // depth. Without this rule the first row in the list wins, and that row can change
            // between two runs. A default that changes is worse than either choice.
            //
            // The realigned row goes first because the user asked for it, and because it is on the
            // newer reference. This rule sets a *default*. It is not a limit, and the user can
            // still select the source.
            let derived = a.is_derived();
            let key = (breadth, depth, a.bam_path.is_some(), derived);
            if best.as_ref().map_or(true, |(b, d, f, r, _)| key > (*b, *d, *f, *r)) {
                best = Some((breadth, depth, a.bam_path.is_some(), derived, a));
            }
        }
        Ok(best.map(|(_, _, _, _, a)| (a.sequence_run_id, a.id)))
    }

    /// The ancestry of a donor. The value is the modern super-population **`ADMIXTURE`** estimate.
    ///
    /// The method takes the consensus estimate, [`CONSENSUS_SOURCE_ID`], which pools each source.
    /// If the store holds none, it takes the estimate of one alignment with the best quality, which
    /// is the estimate with the most genotyped SNPs. That second path supports a result from before
    /// the consensus feature.
    ///
    /// The method must filter on `ADMIXTURE`. The consensus source now also holds a
    /// `FINE_ADMIXTURE` row and an `ANCIENT_ADMIXTURE` row. A read of the *first* consensus row
    /// gives the deep, or ancient, breakdown, which is a separate report. This method must give the
    /// modern super-population breakdown.
    pub async fn donor_ancestry(&self, biosample_guid: SampleGuid) -> Result<Option<(i64, AncestryResult)>, AppError> {
        let all = ancestry_result::for_biosample(self.store.pool(), biosample_guid).await?;
        if let Some(c) = all
            .iter()
            .find(|(id, r)| *id == CONSENSUS_SOURCE_ID && r.method == "ADMIXTURE")
        {
            return Ok(Some(c.clone()));
        }
        Ok(all
            .into_iter()
            .filter(|(_, r)| r.method == "ADMIXTURE")
            .max_by_key(|(_, r)| r.snps_with_genotype))
    }

    /// One stored consensus ancestry estimate. The key is the consensus source with the `method`
    /// value. Two examples are `"FINE_ADMIXTURE"`, which gives the detailed modern populations, and
    /// `"PCA_PROJECTION_GMM"`, which gives the ancient components.
    ///
    /// The query also filters on the subject. An `alignment_id` of 0 does not name one biosample.
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

    /// The private-Y calls of a donor. The value is the **union** of the cached private-Y calls of
    /// each alignment of the subject, and the code applied the self-mask to those calls.
    ///
    /// The method removes a duplicate position and keeps the observation with the most depth. The
    /// terminal comes from the source bucket with the most coverage.
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

    /// The work of [`donor_private_y`](Self::donor_private_y) for **many** subjects at once.
    ///
    /// The method sends two queries in total, and not two for each subject. The first reads the
    /// alignments, and the second reads the artifact rows. The method then joins the rows of each
    /// subject in memory.
    ///
    /// The map holds no entry for a subject with no cached private-Y data. That state means "no
    /// analysis ran". It does not mean "the analysis found nothing".
    ///
    /// This method does **not** use `AlignmentArtifacts`, by design. That type stats each alignment
    /// file first, to check the age of the cache. The app calculates the private-Y data of only a
    /// small part of a cohort.
    ///
    /// A stat call on each other file gives no result and costs much. On a collection that sits on
    /// an external volume, those stat calls were most of the time of the block-tree build. This
    /// method stats only the alignments with a `private_y` row.
    pub(crate) async fn private_y_for_biosamples(
        &self,
        guids: &[SampleGuid],
    ) -> Result<HashMap<SampleGuid, PrivateBucket>, AppError> {
        let mut by_subject: HashMap<SampleGuid, Vec<Alignment>> = HashMap::new();
        for (guid, aln) in alignment::list_for_biosamples(self.store.pool(), guids).await? {
            by_subject.entry(guid).or_default().push(aln);
        }
        let ids: Vec<i64> = by_subject.values().flatten().map(|a| a.id).collect();
        let stored: HashMap<i64, AnalysisArtifact> =
            artifact::list_for_alignments_of_kind(self.store.pool(), &ids, "private_y", "4")
                .await?
                .into_iter()
                .map(|a| (a.alignment_id, a))
                .collect();
        if stored.is_empty() {
            return Ok(HashMap::new());
        }

        // The buckets from a VCF file. They cover a subject whose Y data arrived with no
        // alignment, and such subjects are most of a Y project.
        //
        // The code reads the cache only, as the alignment path does. A user who opens a tab must
        // not start a classification of thousands of call sets.
        let mut out: HashMap<SampleGuid, PrivateBucket> = HashMap::new();
        for guid in guids {
            let rows =
                navigator_store::variant_set_private_y::list_for_biosample(self.store.pool(), &guid.0.to_string())
                    .await
                    .unwrap_or_default();
            for (_, json) in rows {
                let Ok(bucket) = serde_json::from_str::<PrivateBucket>(&json) else {
                    continue;
                };
                let entry = out.entry(*guid).or_insert_with(|| PrivateBucket {
                    terminal: bucket.terminal.clone(),
                    variants: Vec::new(),
                });
                merge_bucket(entry, bucket);
            }
        }
        for (guid, alignments) in &by_subject {
            // Same union-by-position, keep-the-deepest merge as the single-subject path.
            let mut by_pos: HashMap<i64, PrivateVariant> = HashMap::new();
            let mut terminal: Option<String> = None;
            let mut any = false;
            for a in alignments {
                let Some(row) = stored.get(&a.id) else { continue };
                // A stat call is worth its cost only at this point.
                let current = a.bam_path.as_deref().and_then(|p| file_signature(Path::new(p)));
                if !artifact_is_fresh(row.source_sig.as_deref(), current.as_deref()) {
                    continue;
                }
                let Ok(bucket) = serde_json::from_str::<PrivateBucket>(&row.payload) else {
                    continue;
                };
                any = true;
                terminal.get_or_insert(bucket.terminal);
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
                continue;
            }
            let mut variants: Vec<PrivateVariant> = by_pos.into_values().collect();
            variants.sort_by_key(|v| v.position);
            let from_alignments = PrivateBucket {
                terminal: terminal.unwrap_or_default(),
                variants,
            };
            match out.entry(*guid) {
                std::collections::hash_map::Entry::Occupied(mut e) => merge_bucket(e.get_mut(), from_alignments),
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(from_alignments);
                }
            }
        }
        Ok(out)
    }

    /// Every alignment a subject owns, across all of its sequencing runs.
    pub async fn list_alignments_for_biosample(&self, biosample_guid: SampleGuid) -> Result<Vec<Alignment>, AppError> {
        Ok(alignment::list_for_biosample(self.store.pool(), biosample_guid).await?)
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
        // One grouped count covers the full workspace. The code does not send a COUNT query for
        // each project. This code runs at each load of the projects list, and at each `--project`
        // name lookup in the CLI.
        let counts: HashMap<i64, i64> = biosample::member_counts(self.store.pool()).await?.into_iter().collect();
        Ok(project::list(self.store.pool())
            .await?
            .into_iter()
            .map(|project| {
                // Absent from the grouped result means no members.
                let sample_count = counts.get(&project.id).copied().unwrap_or(0);
                ProjectOverview { project, sample_count }
            })
            .collect())
    }

    /// A report of each sample in a project. Each row holds the count of alignments of one
    /// biosample, a coverage summary, and the Y and mtDNA haplogroup consensus. The coverage comes
    /// from the first alignment with a cached coverage result.
    ///
    /// The method calls the queries for one subject that already exist, and it adds no join. A
    /// coverage cell and a haplogroup cell hold `None` until those analyses run.
    pub async fn project_report(&self, project_id: i64) -> Result<Vec<ProjectSampleReport>, AppError> {
        let members = biosample::list_members_for_project(self.store.pool(), project_id).await?;
        // Four queries read each value that this report needs. The code does not send one query
        // for each cell. The queries read the members, their alignments, each artifact of those
        // alignments, and the haplogroup reconciliation. Each cell below is then a lookup in
        // memory.
        let guids: Vec<SampleGuid> = members.iter().map(|b| b.guid).collect();
        let mut by_subject: HashMap<SampleGuid, Vec<Alignment>> = HashMap::new();
        for (guid, aln) in alignment::list_for_biosamples(self.store.pool(), &guids).await? {
            by_subject.entry(guid).or_default().push(aln);
        }
        let all_alignments: Vec<&Alignment> = by_subject.values().flatten().collect();
        let artifacts = AlignmentArtifacts::load(&self.store, &all_alignments).await?;
        // The order is the same as the order in `haplogroup_consensus`. It is the vote of each
        // run, then the placed label, then a value that the user set.
        let terminals = self.haplogroup_terminals().await?;

        let mut out = Vec::new();
        for biosample in members {
            let alignments = by_subject.remove(&biosample.guid).unwrap_or_default();
            let mut coverage: Option<CoverageResult> = None;
            let mut coverage_aln = None;
            for a in &alignments {
                if let Some(c) = artifacts.fresh(a.id, "coverage", coverage::COVERAGE_VERSION) {
                    coverage = Some(c);
                    coverage_aln = Some(a.id);
                    break;
                }
            }
            // The row marks a small coverage result from a sidecar. The UI can then show a badge
            // and offer a full walk.
            let coverage_partial = match coverage_aln {
                Some(id) => matches!(
                    artifacts.provenance(id, "coverage", coverage::COVERAGE_VERSION),
                    Some((_, ref c)) if c == "partial"
                ),
                None => false,
            };
            // Take the alignment with a coverage result. If there is none, take the first
            // alignment.
            let primary_alignment_id = coverage_aln.or_else(|| alignments.first().map(|a| a.id));
            let (y_haplogroup, mt_haplogroup) = terminals.get(&biosample.guid).cloned().unwrap_or_default();
            // Sex + read-metrics from whichever alignment has them cached.
            let mut sex = None;
            let mut metrics: Option<ReadMetrics> = None;
            let mut sv_count = None;
            for a in &alignments {
                if sex.is_none() {
                    sex = artifacts.fresh::<navigator_analysis::sex::SexInferenceResult>(a.id, "sex", "1");
                }
                if metrics.is_none() {
                    metrics = artifacts.fresh(a.id, "read_metrics", "1");
                }
                if sv_count.is_none() {
                    sv_count = artifacts
                        .fresh::<navigator_analysis::sv::types::SvAnalysisResult>(a.id, "sv", "1")
                        .map(|s| s.sv_calls.len());
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
                // Show a stored failure only when the row has no coverage. Such a failure comes
                // from a file that the decoder refuses. A good walk removes the mark.
                //
                // The code reads the mark and does not check its age, as `analysis_error` does. The
                // mark stays until a good walk removes it.
                decode_error: match (coverage.is_none(), primary_alignment_id) {
                    (true, Some(id)) => artifacts
                        .raw(id, ERROR_KIND, ERROR_VERSION)
                        .and_then(|a| serde_json::from_str::<AnalysisError>(&a.payload).ok())
                        .map(|e| e.message),
                    _ => None,
                },
                biosample,
            });
        }
        Ok(out)
    }

    /// A Y-STR overview of each member of a project. The table has the shape of the FTDNA "Y-DNA
    /// Results Overview".
    ///
    /// The result holds each member with one STR profile or more. A row holds the identity columns,
    /// the terminal Y haplogroup, the STR panel that the test reached, and the consensus marker
    /// values. Each marker name is in upper case.
    ///
    /// The result holds no member with no STR data. The method calls the queries for one subject
    /// that already exist, and it adds no join.
    pub async fn project_str_overview(&self, project_id: i64) -> Result<Vec<ProjectStrMember>, AppError> {
        use navigator_domain::{strpanel, strprofile};
        let mut out = Vec::new();
        // One reconciliation covers the full workspace. The code does not send two queries for
        // each member.
        let terminals = self.haplogroup_terminals().await?;
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

            let y_haplogroup = terminals.get(&biosample.guid).and_then(|(y, _)| y.clone());
            let y_confirmed = y_haplogroup.is_some();

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

    /// Build the Y-STR overview of a project, in the FTDNA form, with each value calculated in
    /// advance.
    ///
    /// The method groups the members by their **assigned** Y haplogroup, which is the consensus
    /// value. It then orders the groups by the shape of the tree, from the basal node to the
    /// derived nodes. Each child group goes below its ancestor group.
    ///
    /// For each group, the method calculates the MIN, MAX, and MODE values. For each cell, it
    /// calculates the difference from the modal value.
    ///
    /// A member with no SNP haplogroup goes into an "Unassigned" group at the base.
    ///
    /// Each large calculation happens here, away from the UI thread. The renderer only reads
    /// [`ProjectStrChart::rows`].
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

    /// Analyze each sample in a project. The method calculates the coverage and assigns the Y
    /// haplogroup on the primary alignment of each sample. That alignment is the first one with a
    /// BAM file. The project report then holds a value in each cell.
    ///
    /// The method skips a coverage result that the cache holds, and a Y value that the store
    /// already holds. So a second run is safe.
    ///
    /// A failure on one sample goes into the report, and the method continues with the other
    /// samples.
    ///
    /// The method does not assign mtDNA, by design. That value is not final on CHM13. See the notes
    /// on the reconciliation and the liftover.
    pub async fn analyze_project(
        &self,
        project_id: i64,
        cancel: navigator_analysis::CancelToken,
    ) -> Result<AnalyzeSummary, AppError> {
        let mut summary = AnalyzeSummary {
            project_id,
            samples: 0,
            coverage_done: 0,
            y_done: 0,
            sex_done: 0,
            metrics_done: 0,
            errors: Vec::new(),
        };
        for biosample in biosample::list_members_for_project(self.store.pool(), project_id).await? {
            // The code checks for a stop between two samples, and also inside one sample. A stop
            // that arrives at the end of the walk of one sample must not start the next sample.
            if cancel.is_cancelled() {
                break;
            }
            let o = self.analyze_biosample(&biosample, cancel.clone()).await?;
            if !o.had_alignment {
                continue;
            }
            summary.samples += 1;
            summary.coverage_done += o.coverage_done as usize;
            summary.y_done += o.y_done as usize;
            summary.sex_done += o.sex_done as usize;
            summary.metrics_done += o.metrics_done as usize;
            summary.errors.extend(o.errors);
        }
        Ok(summary)
    }

    /// Do the full analysis of the primary alignment of one biosample. That alignment is its first
    /// alignment with a BAM file. The steps are the coverage, the Y haplogroup, the sex, and the
    /// read metrics.
    ///
    /// The method does **not** call structural variants. See the note at the end of the body.
    ///
    /// A second run is safe. The method skips a *full* coverage result, and a Y value, a sex value,
    /// and a metrics value that the store holds. It replaces a `partial` coverage result from a
    /// sidecar, because the walk across each base gives a better result.
    ///
    /// A failure in one step goes into `errors`, with the donor id at the start of the message. The
    /// other steps still run.
    ///
    /// This method is the unit of work for one sample. The project pass and the deep-analyze job
    /// both call it.
    pub async fn analyze_biosample(
        &self,
        biosample: &Biosample,
        cancel: navigator_analysis::CancelToken,
    ) -> Result<SampleAnalyzeOutcome, AppError> {
        let mut o = SampleAnalyzeOutcome::default();
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample.guid).await?;
        // Take an alignment whose file is still on disk.
        //
        // An earlier rule took any alignment with a recorded path. Take a subject with two
        // alignments, where a user removed the vendor download of one and kept the other. That rule
        // could take the absent one. The preflight then failed, and the code skipped the full
        // subject. The other alignment would have given a good result.
        //
        // When no file is on disk, the code takes any recorded path. The preflight then gives a
        // real diagnosis. Without that step, the report reads as "this subject has no alignment",
        // which is not true.
        let Some(aln) = alignments
            .iter()
            .find(|a| Self::alignment_file(a).is_ok())
            .or_else(|| alignments.iter().find(|a| a.bam_path.is_some()))
        else {
            return Ok(o); // had_alignment stays false
        };
        o.had_alignment = true;
        let label = &biosample.donor_identifier;

        // Run the preflight before any I/O in the steps below. A batch is the worst place to find
        // a file problem the slow way.
        //
        // Without the preflight, an alignment that the code can not read gives one
        // `io error on <the alignment>` message for each step. Each message arrives after a walk
        // that had to fail first, and no message names the file at fault.
        //
        // The preflight does *not* skip the sample on each failure. A broken index stops only the
        // steps that query a region. The unified metrics walk then reads the file from start to end
        // and still gives the coverage, the read metrics, and the sex. A skip would remove results
        // that the user can get. So only a failure that stops a sequential read ends the work.
        //
        // Both destinations receive the *first failure* and not the full report. `o.errors` shows
        // one line for each entry. `record_analysis_error` keeps the first 500 characters. For a
        // full report, those characters are the path and the checks that passed, so the tool would
        // cut the diagnosis itself. The command `navigator doctor` still gives the full report.
        match self.diagnose_alignment(aln.id).await {
            Ok(report) if report.failed() => {
                let cause = report
                    .first_failure()
                    .map(|c| match &c.path {
                        Some(p) => format!("{} ({}): {}", c.name, p.display(), c.detail),
                        None => format!("{}: {}", c.name, c.detail),
                    })
                    .unwrap_or_else(|| "failed".to_string());
                o.errors.push(format!("{label} preflight: {cause}"));
                if report.blocks_sequential_reads() {
                    self.record_analysis_error(aln.id, "preflight", &cause).await;
                    return Ok(o);
                }
            }
            // A preflight that did not run says nothing about the alignment. Continue, and let
            // the real steps report the fault that they find.
            Ok(_) | Err(_) => {}
        }

        // The unified walker gives the coverage, the read metrics, and the sex in ONE pass. It
        // does not read the BAM file or the CRAM file three times. This change divides the I/O of
        // each subject by three. That I/O is most of the time of a batch on a slow volume or a
        // network volume. The Full Analysis of one subject already used this walker, and the batch
        // path did not.
        //
        // The code walks the file only when a value is absent. The work is complete when the store
        // holds a full coverage result with the correct scope, the read metrics, and the sex. A
        // whole-genome coverage result for a targeted Y test has the wrong scope, and the code
        // calculates it again.
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
            match self
                .run_unified_metrics_with_progress(aln.id, |_, _| {}, cancel.clone())
                .await
            {
                Ok(_) => {
                    o.coverage_done = true;
                    o.metrics_done = true;
                    o.sex_done = true;
                    // An earlier run can leave a failure mark, from a bad file that the user then
                    // replaced. Remove that mark.
                    self.clear_analysis_error(aln.id).await;
                }
                // A stop is the decision of the user. It says nothing about the file.
                //
                // A record of it writes a "Failed" mark that stays after the run, and the sample
                // then looks broken for all time. A count of it as an error also makes the batch
                // summary wrong.
                //
                // So the code stops the work on this sample here. The user would stop each
                // remaining step also.
                Err(e) if e.is_cancellation() => return Ok(o),
                Err(e) => {
                    // Persist the failure so the report can show "Failed" instead of a silent blank
                    // (a corrupt/undecodable CRAM otherwise looks identical to an un-analyzed one).
                    self.record_analysis_error(aln.id, "metrics", &e.to_string()).await;
                    o.errors.push(format!("{label} metrics: {e}"));
                }
            }
        }

        if self.haplogroup_consensus(biosample.guid, DnaType::Y).await?.is_some() {
            o.y_done = true;
        } else {
            match self.assign_y_haplogroup(aln.id).await {
                Ok(_) => o.y_done = true,
                Err(e) if e.is_cancellation() => return Ok(o),
                Err(e) => o.errors.push(format!("{label} Y: {e}")),
            }
        }

        // Build the genome-consensus Y signature here. That work is the deep placement and the
        // variant profile, and it gives the descent report.
        //
        // So a batch analysis fills the Y-DNA descent report, and the user does not press "Build
        // descent report". That button stays, and it rebuilds the report at any time.
        //
        // The code builds the report one time, and it skips a profile that already exists. The step
        // is optional. The Y assignment above wrote the chrY genotypes to the cache a moment ago,
        // so this step reads the file no more times.
        if o.y_done && self.cached_y_profile(biosample.guid).await?.is_none() {
            if let Err(e) = self.build_y_profile(biosample.guid).await {
                o.errors.push(format!("{label} Y signature: {e}"));
            }
        }

        // The SV step does NOT run here, by design. It is experimental, and no other step reads
        // its output. It is also the one step that reads each record in the file for its own
        // result.
        //
        // A measurement on the CRAM files of this workspace gave 2 to 5 hours for one whole-genome
        // sample. Each step above needs about 1 hour in total. In a project of 148 samples, that
        // difference is a batch that completes in one night against a batch that needs weeks.
        //
        // The user starts this step: the "Call SV" button, or `analyze --sv`. Both call
        // `plan_full_analysis(.., include_sv = true)`.
        Ok(o)
    }
}

// ---- Y-tree topology helpers for the project STR chart ordering --------------------------------

/// Change a haplogroup name or a node name into the form that the code compares. The function makes
/// each letter upper case and keeps the full name. Our consensus labels and the tree nodes both use
/// the "R-CTS4466" form.
///
/// The result also compares with the part after the last `-`. So the code accepts a plain SNP name,
/// and a caller indexes both forms.
fn norm_hg(name: &str) -> String {
    name.trim().to_ascii_uppercase()
}

/// A map from each tree node name to a node id. The map also holds the plain SNP part of each name.
/// The code uses it to find the node of a haplogroup label. When two entries have the same key, the
/// full name wins.
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

/// The pre-order rank of each node, from a depth-first search. The order goes from the basal node
/// to the derived nodes. Each child follows its parent, and the code keeps the stored order of the
/// children. A caller then orders its groups in the same shape as the tree.
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

/// Find the `read_type` of a run and read no alignment file. The function reads the `test_type`
/// value, which can already name the chemistry, and then the platform.
///
/// The function returns `None` for a PacBio run with a plain WGS test type. Only the read names
/// separate HiFi from CLR, and a scan of the file gives them.
fn infer_read_type_cheap(platform_name: &str, test_type: &str) -> Option<&'static str> {
    let tt = test_type.to_ascii_uppercase();
    if tt.contains("HIFI") {
        return Some("HIFI");
    }
    if tt.contains("CLR") {
        return Some("CLR");
    }
    if tt.contains("NANOPORE") {
        return Some("ONT_SIMPLEX");
    }
    let p = platform_name.to_ascii_uppercase();
    if p.contains("ILLUMINA") || p.contains("MGI") || p.contains("BGI") || p.contains("DNBSEQ") {
        Some("SHORT")
    } else if p.contains("NANOPORE") || p == "ONT" {
        Some("ONT_SIMPLEX")
    } else {
        // A PacBio platform, or a platform that the code does not know. The code can not separate
        // HiFi from CLR here.
        None
    }
}

/// Evidence-based `read_type` from a run's mean read length, for when the platform is uninformative
/// (e.g. sidecar-imported runs with `UNKNOWN` platform): a short mean ⇒ `SHORT`. Long reads
/// (> 1000 bp) can't be split into HiFi vs CLR by length alone, so they stay unresolved until a
/// rescan reads the names. `None` for a non-positive mean (no metrics).
fn read_type_from_mean_len(mean: f64) -> Option<&'static str> {
    (mean > 0.0 && mean <= 1000.0).then_some("SHORT")
}

/// The rank of a test by the part of a person that it covers. The code uses this rank to select the
/// default alignment of a subject. See [`App::default_alignment_for_subject`].
///
/// A higher value covers more. A whole-genome test carries the paternal ancestry, the maternal
/// ancestry, *and* the autosomal ancestry. An autosomal test or a chip test carries the ancestry
/// composition, which is the first section of the brief. A Y test, an mt test, or an X test covers
/// one lineage only.
///
/// A `None` value is a test type that the code does not know. Its rank is above a targeted test,
/// because the test can be a wide test with a label that the code does not hold. Its rank is below
/// each test that the code knows to be genome-wide.
fn test_breadth_rank(target: Option<navigator_domain::testtype::TargetType>) -> u8 {
    use navigator_domain::testtype::TargetType::*;
    match target {
        Some(WholeGenome) => 4,
        Some(Autosomal) | Some(Mixed) => 3,
        None => 2,
        Some(XChromosome) => 1,
        Some(YChromosome) | Some(MtDna) => 0,
    }
}

#[cfg(test)]
mod breadth_tests {
    use super::test_breadth_rank;
    use navigator_domain::testtype::{target_of, TargetType};

    #[test]
    fn whole_genome_outranks_targeted_regardless_of_depth() {
        // The fault that a user reported. A deep Y Elite test must not rank above a genome-wide
        // WGS test or HiFi test. A Y-only test can not give the ancestry that the app shows beside
        // it.
        let wgs = test_breadth_rank(target_of("WGS"));
        let hifi = test_breadth_rank(target_of("WGS_HIFI"));
        let y_elite = test_breadth_rank(target_of("Y_ELITE"));
        assert!(wgs > y_elite, "WGS should outrank Y Elite");
        assert!(hifi > y_elite, "even the shallow HiFi should outrank Y Elite");
        // Y/mt targeted tests sit at the bottom.
        assert_eq!(y_elite, 0);
        assert_eq!(test_breadth_rank(target_of("MT_FULL_SEQUENCE")), 0);
        assert_eq!(test_breadth_rank(target_of("BIG_Y_700")), 0);
    }

    #[test]
    fn ancestry_bearing_and_unknown_tiers() {
        // A chip and an exome carry the autosomal ancestry, which is the first section of the
        // brief. Their rank is above a targeted test and below a WGS test.
        assert!(test_breadth_rank(Some(TargetType::WholeGenome)) > test_breadth_rank(target_of("ARRAY_23ANDME_V5")));
        assert!(test_breadth_rank(target_of("ARRAY_23ANDME_V5")) > test_breadth_rank(target_of("Y_ELITE")));
        // An unrecognized test type ranks above targeted but below a known genome-wide test.
        let unknown = test_breadth_rank(None);
        assert!(unknown > test_breadth_rank(target_of("Y_ELITE")));
        assert!(unknown < test_breadth_rank(target_of("WGS")));
    }
}

#[cfg(test)]
mod read_profile_tests {
    use super::infer_read_type_cheap;

    #[test]
    fn cheap_read_type_inference() {
        // test_type code names the chemistry.
        assert_eq!(infer_read_type_cheap("PACBIO", "WGS_HIFI"), Some("HIFI"));
        assert_eq!(infer_read_type_cheap("PACBIO", "WGS_CLR"), Some("CLR"));
        assert_eq!(infer_read_type_cheap("NANOPORE", "WGS_NANOPORE"), Some("ONT_SIMPLEX"));
        // Short-read platforms.
        assert_eq!(infer_read_type_cheap("ILLUMINA", "WGS"), Some("SHORT"));
        assert_eq!(infer_read_type_cheap("MGI", "WGS"), Some("SHORT"));
        // Nanopore by platform alone.
        assert_eq!(infer_read_type_cheap("NANOPORE", "WGS"), Some("ONT_SIMPLEX"));
        // A PacBio platform with a plain WGS test type. The code needs a scan of the file.
        assert_eq!(infer_read_type_cheap("PACBIO", "WGS"), None);
    }

    #[test]
    fn read_type_from_mean_len_resolves_short() {
        use super::read_type_from_mean_len;
        // Sidecar-imported short-read WGS (mean 62–150 bp) → SHORT.
        assert_eq!(read_type_from_mean_len(150.0), Some("SHORT"));
        assert_eq!(read_type_from_mean_len(62.5), Some("SHORT"));
        // The length of a long read does not separate HiFi from CLR. The code needs a scan of the
        // file.
        assert_eq!(read_type_from_mean_len(21_563.0), None);
        // No metrics.
        assert_eq!(read_type_from_mean_len(0.0), None);
    }
}

/// Add `extra` to `into`. The function removes a duplicate position and keeps the observation with
/// the most depth. The union of the alignments of one subject uses the same rule.
///
/// A donor with a CRAM file and a vendor VCF must have one private set. Two sets that disagree are
/// not correct.
fn merge_bucket(into: &mut PrivateBucket, extra: PrivateBucket) {
    if into.terminal.is_empty() {
        into.terminal = extra.terminal;
    }
    let mut by_pos: HashMap<i64, PrivateVariant> = into.variants.drain(..).map(|v| (v.position, v)).collect();
    for v in extra.variants {
        by_pos
            .entry(v.position)
            .and_modify(|cur| {
                if v.depth > cur.depth {
                    *cur = v.clone();
                }
            })
            .or_insert(v);
    }
    into.variants = by_pos.into_values().collect();
    into.variants.sort_by_key(|v| v.position);
}
