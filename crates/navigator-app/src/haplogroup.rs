//! `impl App` methods extracted from `lib.rs` (the `haplogroup` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// Callable chrY bases from a coverage result — the FTDNA Big Y generation discriminator (see
/// [`App::refine_big_y_generation`]). A Big Y has reads only on chrY, so this is its whole callable
/// footprint. Build-agnostic (`chrY`/`Y`).
fn callable_chr_y_bases(cov: &Coverage) -> u64 {
    cov.contig_callable
        .iter()
        .filter(|c| matches!(c.contig.as_str(), "chrY" | "Y"))
        .map(|c| c.callable)
        .sum()
}

impl App {
    // ---- result exports (gap §6) -------------------------------------------

    /// Format a cached result as a shareable file body (TSV / HTML / BED). The UI writes the
    /// returned string to the user-chosen path. Errors when the source result hasn't been computed
    /// yet (`NotFound`). [`ExportRequest::CallableBed`] re-walks the BAM (no cached intervals).
    pub async fn export_content(&self, req: &ExportRequest) -> Result<String, AppError> {
        match req {
            ExportRequest::CoverageTsv(id) => Ok(export::coverage_tsv(&self.require_coverage(*id).await?)),
            ExportRequest::CoverageHtml(id) => Ok(export::coverage_html(
                &self.require_coverage(*id).await?,
                &format!("alignment {id}"),
            )),
            ExportRequest::ReadMetricsTsv(id) => {
                let m = self
                    .cached_read_metrics(*id)
                    .await?
                    .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("read metrics for alignment {id}"))))?;
                Ok(export::read_metrics_tsv(&m))
            }
            ExportRequest::AncestryTsv(id) => Ok(export::ancestry_tsv(&self.require_ancestry(*id).await?)),
            ExportRequest::AncestryHtml(id) => Ok(export::ancestry_html(&self.require_ancestry(*id).await?)),
            ExportRequest::MtdnaTsv(id) => Ok(export::mtdna_variants_tsv(&self.mtdna_variants(*id).await?)),
            ExportRequest::CallableBed(id) => {
                let per_contig = self.callable_intervals_all(*id).await?;
                Ok(export::callable_bed(&per_contig))
            }
            ExportRequest::DiploidVcf(id) => self.diploid_vcf_genome(*id).await,
            ExportRequest::ConsensusDiploidVcf(guid) => self.consensus_diploid_vcf(*guid).await,
            ExportRequest::SubjectBriefHtml(guid) => {
                let brief = self.subject_brief(*guid).await?;
                // Fold in the AI story if one is already cached (no generation during an export).
                let narration = self.cached_narration(&brief);
                Ok(export::subject_brief_html(&brief, narration.as_ref()))
            }
        }
    }

    async fn require_coverage(&self, alignment_id: i64) -> Result<CoverageResult, AppError> {
        self.cached_coverage(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("coverage for alignment {alignment_id}"))))
    }

    async fn require_ancestry(&self, alignment_id: i64) -> Result<AncestryResult, AppError> {
        self.ancestry_for_alignment(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("ancestry for alignment {alignment_id}"))))
    }

    /// Walk each analyzed contig for its CALLABLE intervals (BED export). Re-reads the BAM — the
    /// coverage artifact stores only per-contig callable *counts*, not the intervals. Uses the
    /// contig list from the cached coverage result, so coverage must have been run first.
    async fn callable_intervals_all(&self, alignment_id: i64) -> Result<Vec<(String, Vec<(i64, i64)>)>, AppError> {
        let cov = self.require_coverage(alignment_id).await?;
        let contigs: Vec<String> = cov.contig_coverage_stats.iter().map(|s| s.contig.clone()).collect();
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let out = tokio::task::spawn_blocking(move || {
            let mut params = CallableLociParams::default();
            if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(&bam, reference.as_deref()) {
                params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            }
            let mut per_contig = Vec::new();
            for contig in contigs {
                // A contig with no aligned reads / bad region just contributes no intervals.
                let intervals =
                    coverage::callable_intervals(&bam, &contig, &params, 1, reference.as_deref()).unwrap_or_default();
                if !intervals.is_empty() {
                    per_contig.push((contig, intervals));
                }
            }
            per_contig
        })
        .await?;
        Ok(out)
    }

    pub async fn derive_mtdna_variants(&self, mtdna_id: i64, rcrs_path: &Path) -> Result<VariantSet, AppError> {
        let seq = mtdna_store::get(self.store.pool(), mtdna_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("mtDNA sequence {mtdna_id}"))))?;
        let rcrs_text = std::fs::read_to_string(rcrs_path)?;
        let rcrs = mtdna::parse_fasta(&rcrs_text).map_err(|e| AppError::Import(format!("rCRS reference: {e}")))?;

        let derived = navigator_analysis::mtvariants::derive(&rcrs.sequence, &seq.sequence);
        let calls = derived
            .iter()
            .map(|v| variants::VariantCall {
                contig: "rCRS".to_string(),
                position: v.position,
                reference: v.reference.to_string(),
                alternate: v.alternate.to_string(),
                rs_id: None,
                genotype: None,
            })
            .collect();
        let label = format!("mtDNA vs rCRS ({} variants)", derived.len());
        let new = NewVariantSet {
            biosample_guid: seq.biosample_guid,
            source_label: label,
            source_type: variants::SourceType::Imported,
            reference_build: None,
            calls,
        };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// Assign an mtDNA haplogroup to a stored sequence: fetch (and cache) the FTDNA mt-DNA
    /// haplotree and rank haplogroups by the Kulczynski measure over the sample's base
    /// calls. RSRS-anchored and reference-free (no rCRS needed). Best first.
    pub async fn assign_mtdna_haplogroup(&self, mtdna_id: i64) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let assignment = self.assign_mtdna_haplogroup_with_tree(mtdna_id, &tree_json).await?;
        if let Some(seq) = mtdna_store::get(self.store.pool(), mtdna_id).await? {
            self.record_call(
                seq.biosample_guid,
                DnaType::Mt,
                &format!("mtseq:{mtdna_id}"),
                format!("mtDNA seq #{mtdna_id}"),
                &assignment,
            )
            .await?;
        }
        Ok(assignment)
    }

    /// The biosample an alignment belongs to (alignment → sequencing run → biosample).
    pub async fn biosample_of_alignment(&self, alignment_id: i64) -> Result<SampleGuid, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let run = sequence_run::get(self.store.pool(), aln.sequence_run_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("sequence run {}", aln.sequence_run_id))))?;
        Ok(run.biosample_guid)
    }

    /// Record (upsert) a source's haplogroup call for donor-level reconciliation.
    pub async fn record_haplogroup_call(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        call: &RunHaplogroupCall,
    ) -> Result<(), AppError> {
        self.record_haplogroup_call_fp(biosample_guid, dna_type, source_key, call, None)
            .await
    }

    /// Like [`record_haplogroup_call`](Self::record_haplogroup_call) but stamps the input
    /// fingerprint (file + tree content hashes) so a later run can skip re-scoring.
    async fn record_haplogroup_call_fp(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        call: &RunHaplogroupCall,
        fingerprint: Option<&str>,
    ) -> Result<(), AppError> {
        haplogroup_call::upsert(
            self.store.pool(),
            biosample_guid,
            dna_type,
            source_key,
            call,
            fingerprint,
        )
        .await?;
        self.audit(
            biosample_guid,
            dna_type,
            "RUN_RECORDED",
            &format!("{source_key}: {}", call.haplogroup),
        )
        .await?;
        Ok(())
    }

    /// Record an assignment's top candidate as a per-source call (no-op if no match).
    async fn record_call(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        source_label: String,
        assignment: &HaploAssignment,
    ) -> Result<(), AppError> {
        self.record_call_fp(biosample_guid, dna_type, source_key, source_label, assignment, None)
            .await
    }

    /// Like [`record_call`](Self::record_call) but stamps the input fingerprint.
    pub(crate) async fn record_call_fp(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        source_label: String,
        assignment: &HaploAssignment,
        fingerprint: Option<&str>,
    ) -> Result<(), AppError> {
        if let Some(top) = assignment.ranked.first() {
            let call = RunHaplogroupCall {
                source_label,
                haplogroup: top.name.clone(),
                lineage: top.lineage.clone(),
                score: top.score,
                matched: top.matched as i64,
                expected: top.expected as i64,
            };
            self.record_haplogroup_call_fp(biosample_guid, dna_type, source_key, &call, fingerprint)
                .await?;
        }
        Ok(())
    }

    /// The reconciled donor-level haplogroup consensus across all recorded sources. A user
    /// manual override, when set, replaces the computed terminal (flagged `overridden`).
    pub async fn haplogroup_consensus(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<Consensus>, AppError> {
        let calls = haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?;
        // Per-run label reconciliation supplies the lineage / compatibility / divergence warnings…
        let mut consensus = reconciliation::reconcile(&calls);

        // …but the authoritative terminal is the genome-level PLACED call persisted by
        // build_{y,mt}_profile (consensus_profile.consensus_label) — cheap to read, no genotyping.
        // When it disagrees with the per-run vote, keep the placed call and record why.
        if matches!(dna_type, DnaType::Y | DnaType::Mt) {
            if let Some(stored) =
                navigator_store::consensus_profile::get(self.store.pool(), biosample_guid, dna_type.as_str()).await?
            {
                if let Some(placed) = stored.consensus_label.filter(|s| !s.is_empty()) {
                    let mut c = consensus.unwrap_or_else(|| Consensus {
                        haplogroup: placed.clone(),
                        lineage: vec![placed.clone()],
                        compatibility: CompatibilityLevel::Compatible,
                        divergence_point: None,
                        confidence: 1.0,
                        run_count: calls.len(),
                        overridden: false,
                        warnings: Vec::new(),
                    });
                    if c.haplogroup != placed {
                        c.warnings
                            .push(format!("per-run calls vary; genome-consensus placement → {placed}"));
                        c.haplogroup = placed;
                    }
                    consensus = Some(c);
                }
            }
        }

        if let Some((hg, reason)) = recon_store::get_override(self.store.pool(), biosample_guid, dna_type).await? {
            let mut c = consensus.unwrap_or(Consensus {
                haplogroup: hg.clone(),
                lineage: vec![hg.clone()],
                compatibility: CompatibilityLevel::Compatible,
                divergence_point: None,
                confidence: 1.0,
                run_count: 0,
                overridden: true,
                warnings: Vec::new(),
            });
            c.haplogroup = hg;
            c.overridden = true;
            c.confidence = 1.0;
            c.warnings.push(match reason {
                Some(r) => format!("manual override: {r}"),
                None => "manual override".to_string(),
            });
            consensus = Some(c);
        }
        Ok(consensus)
    }

    /// Donor-level Y and mtDNA terminal haplogroups for **every** subject, for the subjects
    /// list. Reconciles each subject's recorded calls (and applies any manual override) in
    /// memory from two bulk queries. `(guid → (Y terminal, mt terminal))`; either is `None`
    /// when nothing is recorded.
    pub async fn haplogroup_terminals(
        &self,
    ) -> Result<HashMap<SampleGuid, (Option<String>, Option<String>)>, AppError> {
        let mut groups: HashMap<(SampleGuid, DnaType), Vec<RunHaplogroupCall>> = HashMap::new();
        for (guid, dna_type, call) in haplogroup_call::list_all(self.store.pool()).await? {
            groups.entry((guid, dna_type)).or_default().push(call);
        }
        let mut out: HashMap<SampleGuid, (Option<String>, Option<String>)> = HashMap::new();
        for ((guid, dna_type), calls) in groups {
            if let Some(c) = reconciliation::reconcile(&calls) {
                let entry = out.entry(guid).or_default();
                match dna_type {
                    DnaType::Y => entry.0 = Some(c.haplogroup),
                    DnaType::Mt => entry.1 = Some(c.haplogroup),
                }
            }
        }
        // The genome-level placed terminal (build_{y,mt}_profile) wins over the per-run label vote,
        // so the subjects table matches the detail tab.
        for (guid_s, dna_type_s, label) in navigator_store::consensus_profile::list_labels(self.store.pool()).await? {
            let Ok(uuid) = guid_s.parse::<uuid::Uuid>() else {
                continue;
            };
            let entry = out.entry(SampleGuid(uuid)).or_default();
            match dna_type_s.as_str() {
                "Y" => entry.0 = Some(label),
                "Mt" => entry.1 = Some(label),
                _ => {}
            }
        }
        // Manual overrides win over everything.
        for (guid, dna_type, hg) in recon_store::list_all_overrides(self.store.pool()).await? {
            let entry = out.entry(guid).or_default();
            match dna_type {
                DnaType::Y => entry.0 = Some(hg),
                DnaType::Mt => entry.1 = Some(hg),
            }
        }
        Ok(out)
    }

    /// Manually override the consensus haplogroup for a subject + DNA type.
    pub async fn set_manual_override(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        haplogroup: &str,
        reason: Option<&str>,
    ) -> Result<(), AppError> {
        recon_store::set_override(self.store.pool(), biosample_guid, dna_type, haplogroup, reason).await?;
        self.audit(
            biosample_guid,
            dna_type,
            "MANUAL_OVERRIDE",
            &format!("override to {haplogroup}"),
        )
        .await
    }

    /// Clear a manual override.
    pub async fn clear_manual_override(&self, biosample_guid: SampleGuid, dna_type: DnaType) -> Result<(), AppError> {
        recon_store::clear_override(self.store.pool(), biosample_guid, dna_type).await?;
        self.audit(biosample_guid, dna_type, "OVERRIDE_CLEARED", "cleared manual override")
            .await
    }

    /// The reconciliation audit log for a subject + DNA type.
    pub async fn reconciliation_audit(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Vec<AuditEntry>, AppError> {
        Ok(recon_store::list_audit(self.store.pool(), biosample_guid, dna_type).await?)
    }

    async fn audit(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        action: &str,
        note: &str,
    ) -> Result<(), AppError> {
        let entry = AuditEntry {
            timestamp: Utc::now().to_rfc3339(),
            action: action.to_string(),
            note: note.to_string(),
        };
        recon_store::append_audit(self.store.pool(), biosample_guid, dna_type, &entry).await?;
        Ok(())
    }

    /// The persisted consensus-profile snapshot for a subject + DNA type, if built — cheap (no
    /// genotyping). The shared loader behind [`cached_y_profile`](Self::cached_y_profile); the future
    /// mtDNA / autosomal tabs reuse it with a different [`DnaType`]. `None` until a build runs.
    pub async fn cached_consensus_profile(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<ConsensusProfile>, AppError> {
        match navigator_store::consensus_profile::get(self.store.pool(), biosample_guid, dna_type.as_str()).await? {
            Some(row) => Ok(Some(serde_json::from_str(&row.payload)?)),
            None => Ok(None),
        }
    }

    /// Persist a reconciled consensus snapshot — the low-level row writer shared by every DNA type
    /// (Y / mt key on [`DnaType`], autosomal keys on `"Auto"`; the payload is whatever profile shape
    /// that type uses). The scalar columns mirror the summary header for quick listing.
    #[allow(clippy::too_many_arguments)]
    async fn persist_consensus_row(
        &self,
        biosample_guid: SampleGuid,
        dna_type: &str,
        consensus_label: Option<String>,
        summary: &navigator_domain::consensus::ConsensusSummary,
        source_count: usize,
        tree_provider: Option<String>,
        payload: String,
    ) -> Result<(), AppError> {
        let stored = navigator_store::consensus_profile::StoredConsensusProfile {
            biosample_guid: biosample_guid.0.to_string(),
            dna_type: dna_type.to_string(),
            consensus_label,
            overall_confidence: summary.overall_confidence,
            source_count: source_count as i64,
            total: summary.total as i64,
            confirmed: summary.confirmed as i64,
            novel: summary.novel as i64,
            conflict: summary.conflict as i64,
            single_source: summary.single_source as i64,
            tree_provider,
            payload,
            last_reconciled_at: Utc::now().to_rfc3339(),
        };
        navigator_store::consensus_profile::upsert(self.store.pool(), &stored).await?;
        Ok(())
    }

    /// Persist a Y/mt [`ConsensusProfile`] snapshot via [`persist_consensus_row`].
    async fn persist_consensus_profile(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        profile: &ConsensusProfile,
        tree_provider: Option<String>,
    ) -> Result<(), AppError> {
        self.persist_consensus_row(
            biosample_guid,
            dna_type.as_str(),
            profile.terminal.clone(),
            &profile.summary,
            profile.sources.len(),
            tree_provider,
            serde_json::to_string(profile)?,
        )
        .await
    }

    /// The persisted Y-profile snapshot for a subject, if one has been built — cheap (no
    /// genotyping). `None` until [`build_y_profile`](Self::build_y_profile) runs.
    pub async fn cached_y_profile(&self, biosample_guid: SampleGuid) -> Result<Option<YProfile>, AppError> {
        self.cached_consensus_profile(biosample_guid, DnaType::Y).await
    }

    /// Build (and persist) the multi-source Y-variant profile: reconcile each Y-bearing source's
    /// per-SNP calls — every alignment's haplogroup placement, the combined chip/BISDNA placement,
    /// and the private-Y bucket — into one concordance view (confirmed / novel / conflict /
    /// single-source per SNP, with per-source provenance + per-observation quality weighting).
    /// Expensive (re-genotypes each alignment), so it's an explicit action; the result is persisted
    /// so [`cached_y_profile`](Self::cached_y_profile) reloads it instantly. Sources without Y data
    /// are skipped.
    pub async fn build_y_profile(&self, biosample_guid: SampleGuid) -> Result<YProfile, AppError> {
        let mut sources: Vec<(String, SourceType, Vec<YObsInput>)> = Vec::new();

        // One source per alignment — a *fresh* placement (the cached terminal-only path lacks the
        // per-SNP branch evidence we reconcile here). Expensive; this is why the profile is built
        // on explicit request. Alignments that error / lack chrY are skipped.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            let Ok(assignment) = self.y_assignment_full(a.id).await else {
                continue;
            };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
            }
        }

        // One source *per chip/BISDNA panel* (a distinct VariantSet per import — 23andMe,
        // AncestryDNA, BISDNA chromo2, …), so the profile shows which test confirmed each SNP and a
        // single mistyped panel surfaces as a conflict rather than being averaged into "consumer tests".
        let vsets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let chip_sets: Vec<&VariantSet> = vsets.iter().filter(|s| s.source_type == SourceType::Chip).collect();
        if !chip_sets.is_empty() {
            // Resolve the placement build once (a chip set's stored build, else the alignment's), and
            // fetch the tree once for all panels.
            let build = chip_sets
                .iter()
                .find_map(|s| s.reference_build.clone())
                .unwrap_or(self.bisdna_target_build(biosample_guid).await);
            if let Ok(tree) = self.chip_y_tree(&build).await {
                for set in &chip_sets {
                    let calls: HashMap<i64, char> = set
                        .calls
                        .iter()
                        .filter(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"))
                        .filter_map(|c| c.alternate.chars().next().map(|b| (c.position, b.to_ascii_uppercase())))
                        .collect();
                    if calls.is_empty() {
                        continue;
                    }
                    let assignment = Self::place_chip_panel(&tree, calls);
                    let obs = snp_obs_from_assignment(&assignment, true);
                    if !obs.is_empty() {
                        sources.push((set.source_label.clone(), SourceType::Chip, obs));
                    }
                }
            }
        }

        // One source *per vendor Y-NGS VCF* (FTDNA Big Y / YSEQ / Full Genomes / Nebula / Dante —
        // every non-chip VariantSet with chrY calls). Placed against each set's stored build so a
        // GRCh38 Big Y reconciles alongside any WGS alignment; tagged with the set's real source
        // type (TargetedNgs / WgsShortRead / …) so it carries the right concordance weight.
        let ngs_sets: Vec<&VariantSet> = vsets.iter().filter(|s| s.source_type != SourceType::Chip).collect();
        if !ngs_sets.is_empty() {
            let mut tree_cache: HashMap<String, navigator_analysis::haplo::HaploTree> = HashMap::new();
            for set in &ngs_sets {
                let calls = Self::vset_chr_y_calls(set);
                if calls.is_empty() {
                    continue;
                }
                let build = set.reference_build.clone().unwrap_or_else(|| "GRCh38".to_string());
                if !tree_cache.contains_key(&build) {
                    match self.chip_y_tree(&build).await {
                        Ok(t) => {
                            tree_cache.insert(build.clone(), t);
                        }
                        Err(_) => continue,
                    }
                }
                let assignment = Self::place_chip_panel(&tree_cache[&build], calls);
                let obs = snp_obs_from_assignment(&assignment, true);
                if !obs.is_empty() {
                    sources.push((set.source_label.clone(), set.source_type, obs));
                }
            }
        }

        // Private-Y union: off-path / novel calls (not in the tree).
        if let Some(bucket) = self.donor_private_y(biosample_guid).await? {
            let obs: Vec<YObsInput> = bucket
                .variants
                .iter()
                .map(|v| {
                    let name = match &v.class {
                        PrivateClass::OffPathKnown(n) => n.clone(),
                        PrivateClass::Novel => String::new(), // keyed by position
                    };
                    let mut o = YObsInput::snp(
                        name,
                        v.position,
                        v.reference.to_string(),
                        v.alternate.to_string(),
                        YState::Derived,
                        false,
                    );
                    // De-novo calls carry read depth; a structural-region (palindrome/amplicon) call
                    // is paralog-suspect → down-weight via the region modifier.
                    o.depth = Some(v.depth);
                    // Down-weight by the structural-region quality modifier (palindrome 0.4,
                    // ampliconic 0.3, heterochromatin/centromere 0.1…); unique sequence = 1.0.
                    o.region_modifier = v.region.map_or(1.0, |c| c.modifier());
                    o
                })
                .collect();
            if !obs.is_empty() {
                sources.push(("private".to_string(), SourceType::WgsShortRead, obs));
            }
        }

        // Provenance: one entry per contributing source (label, type, SNP count).
        let source_summaries: Vec<YSourceSummary> = sources
            .iter()
            .map(|(label, st, obs)| YSourceSummary {
                label: label.clone(),
                source_type: *st,
                variant_count: obs.len(),
            })
            .collect();

        let variants = yprofile::reconcile_y(&sources);
        let summary = yprofile::summarize(&variants);
        // Genome-level placement: place the pooled call set on the tree once (not a vote among the
        // per-run terminal labels). Falls back to the label reconciliation when nothing places.
        let terminal = match self.place_y_consensus(biosample_guid).await? {
            Some(a) => a.ranked.first().map(|r| r.name.clone()),
            None => self
                .haplogroup_consensus(biosample_guid, DnaType::Y)
                .await?
                .map(|c| c.haplogroup),
        };
        let profile = ConsensusProfile {
            variants,
            summary,
            terminal,
            sources: source_summaries,
        };

        // Persist the snapshot (keyed dna_type='Y') so the tab reloads it without re-genotyping.
        let provider = Some(match y_tree_provider() {
            YTreeProvider::DecodingUs => "decodingus".to_string(),
            YTreeProvider::Ftdna => "ftdna".to_string(),
        });
        self.persist_consensus_profile(biosample_guid, DnaType::Y, &profile, provider)
            .await?;
        Ok(profile)
    }

    /// Fresh mtDNA placement against the FTDNA mt tree (chrM) with full branch evidence — the mt
    /// counterpart to [`y_assignment_full`](Self::y_assignment_full). Bypasses
    /// [`assign_mtdna_haplogroup_from_alignment`](Self::assign_mtdna_haplogroup_from_alignment)'s
    /// cached terminal-only path (which has `branches: []`) so the consensus has per-mutation evidence.
    async fn mt_assignment_full(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        self.assign_haplogroup_from_alignment(alignment_id, "chrM", &tree_json)
            .await
    }

    /// The persisted mtDNA consensus-profile snapshot for a subject, if one has been built — cheap
    /// (no genotyping). `None` until [`build_mt_profile`](Self::build_mt_profile) runs.
    pub async fn cached_mt_profile(&self, biosample_guid: SampleGuid) -> Result<Option<ConsensusProfile>, AppError> {
        self.cached_consensus_profile(biosample_guid, DnaType::Mt).await
    }

    /// Build (and persist) the multi-source mtDNA consensus profile — the mtDNA adapter over the
    /// generic [`navigator_domain::consensus`] engine. Reconciles each mt-bearing source's
    /// defining-mutation calls (every alignment's chrM placement, each imported mtDNA FASTA
    /// sequence's placement, and the combined chip mtDNA placement) into one concordance view,
    /// keyed by phylotree **mutation name** (rCRS-coordinate, build-independent). Persisted with
    /// `dna_type='Mt'` so [`cached_mt_profile`](Self::cached_mt_profile) reloads it instantly.
    /// Expensive (re-places each alignment's chrM), so it's an explicit action; mt-less sources skip.
    pub async fn build_mt_profile(&self, biosample_guid: SampleGuid) -> Result<ConsensusProfile, AppError> {
        let mut sources: Vec<(String, SourceType, Vec<YObsInput>)> = Vec::new();

        // One source per alignment with chrM — a fresh placement (branch evidence, not the cached
        // terminal). Alignments that error / lack chrM are skipped.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            let Ok(assignment) = self.mt_assignment_full(a.id).await else {
                continue;
            };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
            }
        }

        // One source per imported mtDNA FASTA sequence (FTDNA mtFull / YSEQ) — a finished consensus
        // sequence we ingested, so weight it as `Imported` (provenance/method not ours to vouch for).
        let seqs = self.list_mtdna_sequences(biosample_guid).await?;
        for s in &seqs {
            let Ok(assignment) = self.assign_mtdna_haplogroup(s.id).await else {
                continue;
            };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                let vendor = mt_vendor_label(s.source_file_name.as_deref(), s.defline.as_deref());
                sources.push((format!("{vendor} (mt seq #{})", s.id), SourceType::Imported, obs));
            }
        }

        // The combined chip mtDNA panel (23andMe carries a sparse mt panel). One source — the
        // per-panel split is a deferred follow-on (same as the Y profile).
        if let Ok(assignment) = self.assign_mt_chip(biosample_guid).await {
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push(("Chip mtDNA panel".to_string(), SourceType::Chip, obs));
            }
        }

        // Provenance: one entry per contributing source (label, type, mutation count).
        let source_summaries: Vec<YSourceSummary> = sources
            .iter()
            .map(|(label, st, obs)| YSourceSummary {
                label: label.clone(),
                source_type: *st,
                variant_count: obs.len(),
            })
            .collect();

        let variants = yprofile::reconcile_y(&sources);
        let summary = yprofile::summarize(&variants);
        // Genome-level placement of the pooled chrM call set (see place_mt_consensus); label vote
        // is the fallback when nothing places.
        let terminal = match self.place_mt_consensus(biosample_guid).await? {
            Some(a) => a.ranked.first().map(|r| r.name.clone()),
            None => self
                .haplogroup_consensus(biosample_guid, DnaType::Mt)
                .await?
                .map(|c| c.haplogroup),
        };
        let profile = ConsensusProfile {
            variants,
            summary,
            terminal,
            sources: source_summaries,
        };

        // Persist (keyed dna_type='Mt'); the mt tree is FTDNA-sourced.
        self.persist_consensus_profile(biosample_guid, DnaType::Mt, &profile, Some("ftdna".to_string()))
            .await?;
        Ok(profile)
    }

    /// **Genome-level Y placement**: pool every source's tree-locus genotype (each alignment's
    /// native-build placement calls + each chip/BISDNA panel's chrY calls) into one call set by a
    /// weighted [`pool_bases`] vote — keyed by SNP **name** so sources on different builds merge —
    /// then place that pooled set on one canonical tree **once** via [`assemble_assignment`]. This
    /// replaces voting among the per-run terminal *labels*: a sparse run no longer drags the call
    /// shallow, and a branch confirmed by any source informs the placement. `Ok(None)` when the
    /// subject has no Y-bearing source. Re-genotypes each source (like [`build_y_profile`]), so it's
    /// only run as part of that explicit action.
    pub async fn place_y_consensus(&self, biosample_guid: SampleGuid) -> Result<Option<HaploAssignment>, AppError> {
        // Genotype every WGS alignment against **one** tree in **one** coordinate system — the FTDNA
        // GRCh38 Y tree, with the existing base_calls liftover for CHM13/GRCh37 sources — then pool
        // by **position** and place once. A single tree+coordinate space keeps allele polarity and
        // coverage consistent across sources; pooling across the DecodingUs hs1 and GRCh38 trees (which
        // can disagree on a SNP's polarity / hs1 coverage) corrupts the backbone and the parsimony
        // guard then vetoes the true deep clade. Chips (sparse, various builds) stay in the variant
        // profile + label reconciliation; the genome placement pools the dense WGS evidence.
        let tree_json = self.fetch_ftdna_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;

        let mut sources: Vec<(SourceType, HashMap<i64, char>)> = Vec::new();
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            // assign_haplogroup_detail returns the GRCh38-coordinate calls (lifted from the
            // alignment's build); sources that lack chrY / a reference are skipped.
            let Ok((_, _, calls)) = self.assign_haplogroup_detail(a.id, "chrY", &tree_json).await else {
                continue;
            };
            if !calls.is_empty() {
                sources.push((SourceType::WgsShortRead, calls));
            }
        }

        // Vendor Y-NGS VCFs (FTDNA Big Y / YSEQ / Full Genomes / Nebula) are dense, direct Y-SNP
        // calls on GRCh38 — the FTDNA tree's native coordinate space — so pool them into the genome
        // consensus alongside the WGS alignments, weighted by their source type. Chips stay out
        // (sparse, various builds; they live in the variant profile + label reconciliation). A
        // non-GRCh38 set is skipped here (its positions wouldn't match the GRCh38 tree).
        let vsets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for set in &vsets {
            if set.source_type == SourceType::Chip || !is_grch38_build(&set.reference_build) {
                continue;
            }
            let calls = Self::vset_chr_y_calls(set);
            if !calls.is_empty() {
                sources.push((set.source_type, strand_reconcile_to_tree(&tree, calls)));
            }
        }

        if sources.is_empty() {
            return Ok(None);
        }
        let pooled = pool_votes(&sources);
        Ok(Some(assemble_assignment(&tree, &pooled)))
    }

    /// Assemble a subject's lightweight [`YMatchProfile`] from **cached** data only (no re-genotyping):
    /// the persisted consensus Y profile (derived/novel SNP-name sets + terminal), the terminal's
    /// root→tip lineage from `tree`, and the first imported Y-STR panel's markers. `Ok(None)` when the
    /// subject has neither a placed Y profile nor an STR panel — nothing to match on.
    async fn y_match_profile(
        &self,
        b: &Biosample,
        tree: Option<&navigator_analysis::haplo::HaploTree>,
    ) -> Result<Option<navigator_domain::ymatch::YMatchProfile>, AppError> {
        use navigator_domain::consensus::ConsensusState;

        let str_markers = self
            .list_str_profiles(b.guid)
            .await?
            .first()
            .map(|p| p.markers.clone())
            .unwrap_or_default();

        let (terminal, lineage, derived, novel) = match self.cached_y_profile(b.guid).await? {
            Some(p) => {
                // Lineage (for divergence/LCA) needs the tree; without it, SNP-name sets + STR still match.
                let lineage = match (tree, p.terminal.as_deref()) {
                    (Some(t), Some(term)) => lineage_names(t, term),
                    _ => Vec::new(),
                };
                let mut derived = std::collections::HashSet::new();
                let mut novel = std::collections::HashSet::new();
                for v in &p.variants {
                    if v.consensus == ConsensusState::Derived {
                        if v.in_tree {
                            derived.insert(v.name.clone());
                        } else {
                            novel.insert(v.name.clone());
                        }
                    }
                }
                (p.terminal, lineage, derived, novel)
            }
            None => (
                None,
                Vec::new(),
                std::collections::HashSet::new(),
                std::collections::HashSet::new(),
            ),
        };

        // Nothing matchable: no Y-SNP calls and no STR markers (lineage alone never matches).
        if derived.is_empty() && novel.is_empty() && str_markers.is_empty() {
            return Ok(None);
        }
        Ok(Some(navigator_domain::ymatch::YMatchProfile {
            guid: b.guid,
            donor: b.donor_identifier.clone(),
            terminal,
            lineage,
            derived,
            novel,
            str_markers,
        }))
    }

    /// Rank every other workspace subject against `query_guid` by Y relatedness (gap §2) — shared
    /// derived/novel SNPs, divergence haplogroup, Y-STR genetic distance, and rough SNP/STR TMRCA.
    /// One-vs-all over the workspace (or one project when `project_id` is set); local-only. Consumes
    /// **cached** profiles so it's cheap over hundreds of subjects (no re-genotyping). `Ok(vec![])`
    /// when the query subject has no matchable Y data.
    pub async fn y_matches(&self, query_guid: SampleGuid, project_id: Option<i64>) -> Result<Vec<YMatch>, AppError> {
        // The tree only supplies the divergence haplogroup; shared-SNP and STR matching work without
        // it, so a fetch failure degrades gracefully rather than failing the whole search.
        let tree = match self.fetch_ftdna_y_tree().await {
            Ok(json) => navigator_analysis::haplo::parse_ftdna_json(&json).ok(),
            Err(_) => None,
        };

        let candidates = match project_id {
            Some(pid) => self.list_biosamples(pid).await?,
            None => self.list_all_biosamples().await?,
        };

        // The query subject may sit outside the chosen project — load it directly if so.
        let query_bio = match candidates.iter().find(|b| b.guid == query_guid).cloned() {
            Some(b) => b,
            None => biosample::get(self.store.pool(), query_guid)
                .await?
                .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("biosample {query_guid}"))))?,
        };
        let Some(query) = self.y_match_profile(&query_bio, tree.as_ref()).await? else {
            return Ok(Vec::new());
        };

        let mut profiles = Vec::new();
        for b in &candidates {
            if b.guid == query_guid {
                continue;
            }
            if let Some(p) = self.y_match_profile(b, tree.as_ref()).await? {
                profiles.push(p);
            }
        }
        let mut ranked = navigator_domain::ymatch::rank(&query, &profiles);
        ranked.truncate(200);
        Ok(ranked)
    }

    /// **Genome-level mtDNA placement**: the mt counterpart to [`place_y_consensus`]. Pools every
    /// source's rCRS-coordinate genotype (each alignment's chrM placement calls, each imported mtDNA
    /// FASTA sequence, the chip mt panel) by [`pool_bases`] vote keyed by **position** (rCRS is the
    /// only mt coordinate system → no name indirection), then places the pooled set on the FTDNA mt
    /// tree once. `Ok(None)` when the subject has no mt-bearing source.
    pub async fn place_mt_consensus(&self, biosample_guid: SampleGuid) -> Result<Option<HaploAssignment>, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;

        let mut sources: Vec<(SourceType, HashMap<i64, char>)> = Vec::new();
        let mut has_wgs = false;

        // Each alignment's chrM genotype (rCRS coordinates; base_calls maps a CHM13 chrM back to rCRS).
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            let Ok((_, _, calls)) = self.assign_haplogroup_detail(a.id, "chrM", &tree_json).await else {
                continue;
            };
            if !calls.is_empty() {
                sources.push((SourceType::WgsShortRead, calls));
                has_wgs = true;
            }
        }

        // Each imported mtDNA FASTA — the full sequence sampled at every rCRS position.
        for s in &self.list_mtdna_sequences(biosample_guid).await? {
            let Some(seq) = mtdna_store::get(self.store.pool(), s.id).await? else {
                continue;
            };
            let calls: HashMap<i64, char> = seq
                .sequence
                .bytes()
                .enumerate()
                .filter_map(|(i, b)| {
                    let u = b.to_ascii_uppercase();
                    matches!(u, b'A' | b'C' | b'G' | b'T').then_some(((i + 1) as i64, u as char))
                })
                .collect();
            if !calls.is_empty() {
                sources.push((SourceType::Imported, calls));
            }
        }

        // The chip mt panel (consumer arrays carry a sparse rCRS MT panel), strand-reconciled.
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let chip_mt: HashMap<i64, char> = sets
            .iter()
            .filter(|s| s.source_type == SourceType::Chip)
            .flat_map(|s| s.calls.iter())
            .filter(|c| {
                c.contig.eq_ignore_ascii_case("chrM")
                    || c.contig.eq_ignore_ascii_case("chrMT")
                    || c.contig.eq_ignore_ascii_case("mt")
                    || c.contig.eq_ignore_ascii_case("m")
            })
            .filter_map(|c| c.alternate.chars().next().map(|b| (c.position, b.to_ascii_uppercase())))
            .collect();
        if !chip_mt.is_empty() {
            sources.push((SourceType::Chip, strand_reconcile_to_tree(&tree, chip_mt)));
        }

        if sources.is_empty() {
            return Ok(None);
        }

        // mt is single-coordinate (rCRS) across all sources, so a base vote is strand-safe here.
        let pooled = pool_votes(&sources);
        let assignment = if has_wgs {
            assemble_assignment(&tree, &pooled)
        } else {
            assemble_assignment_robust(&tree, &pooled)
        };
        Ok(Some(assignment))
    }

    /// Build a YFull-style [`DescentReport`] for a subject's Y or mtDNA lineage from the **already
    /// persisted** variant profile — no re-genotyping. Reads the cached profile for its terminal +
    /// per-SNP states (keyed by build-independent SNP name), then walks the FTDNA tree from the
    /// terminal to the root, attaching each node's defining SNPs with the sample's call (`NoCall` for
    /// an untested equivalent). `Ok(None)` when the profile isn't built yet or has no terminal — the
    /// UI then offers to build it (one expensive, persisted step that also powers the variant tabs).
    pub async fn descent_report(
        &self,
        biosample_guid: SampleGuid,
        dna: DnaType,
    ) -> Result<Option<DescentReport>, AppError> {
        use navigator_domain::consensus::ConsensusState;

        // Cheap first: the persisted profile. No profile / no terminal → nothing to draw, and we
        // skip the (multi-MB) tree fetch + parse entirely.
        let profile = match dna {
            DnaType::Y => self.cached_y_profile(biosample_guid).await?,
            DnaType::Mt => self.cached_mt_profile(biosample_guid).await?,
        };
        let Some(profile) = profile else { return Ok(None) };
        let Some(terminal) = profile.terminal.clone() else { return Ok(None) };

        let tree_json = match dna {
            DnaType::Y => self.fetch_ftdna_y_tree().await?,
            DnaType::Mt => self.fetch_ftdna_mt_tree().await?,
        };
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;
        let Some(terminal_id) = tree.nodes.iter().find(|(_, n)| n.name == terminal).map(|(id, _)| *id) else {
            return Ok(None); // terminal not in this tree (provider/build skew) — nothing to draw
        };

        let state_by_name: std::collections::HashMap<String, CallState> = profile
            .variants
            .iter()
            .map(|v| {
                let state = match v.consensus {
                    ConsensusState::Derived => CallState::Derived,
                    ConsensusState::Ancestral => CallState::Ancestral,
                    ConsensusState::NoCall => CallState::NoCall,
                };
                (v.name.clone(), state)
            })
            .collect();

        let nodes = navigator_analysis::haplo::descent_by_node(&tree, terminal_id, &state_by_name);
        Ok(Some(DescentReport { dna, terminal, nodes }))
    }

    /// The persisted autosomal consensus-profile snapshot for a subject, if built — cheap (no
    /// genotyping). `None` until [`build_autosomal_profile`](Self::build_autosomal_profile) runs.
    pub async fn cached_autosomal_profile(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<DiploidProfile>, AppError> {
        match navigator_store::consensus_profile::get(self.store.pool(), biosample_guid, "Auto").await? {
            Some(row) => Ok(Some(serde_json::from_str(&row.payload)?)),
            None => Ok(None),
        }
    }

    /// Build (and persist) the multi-source **autosomal** consensus profile — the diploid (0/1/2)
    /// adapter over the generic [`navigator_domain::consensus`] engine. Genotypes every WGS alignment
    /// and imported chip over the canonical CHM13 **IBD panel** ([`ibd_panel_dosages`](Self::ibd_panel_dosages))
    /// and reconciles the per-site dosages into a voted genotype (confirmed where sources agree,
    /// conflict where they don't), keyed by rsID. Persisted with `dna_type='Auto'`. Requires the IBD
    /// panel asset (built with `panelbuild ibd-panel`); errors if it's missing.
    pub async fn build_autosomal_profile(&self, biosample_guid: SampleGuid) -> Result<DiploidProfile, AppError> {
        use navigator_domain::consensus::{reconcile_diploid, summarize_diploid, DiploidObs};

        let to_obs = |gts: Vec<SiteGenotype>| -> Vec<DiploidObs> {
            gts.into_iter()
                .map(|g| DiploidObs {
                    name: g.name,
                    contig: g.contig,
                    position: g.position,
                    reference: g.reference_allele,
                    alternate: g.alternate_allele,
                    dosage: g.dosage as i8,
                    depth: (g.depth > 0).then_some(g.depth),
                })
                .collect()
        };

        let mut sources: Vec<(String, SourceType, Vec<DiploidObs>)> = Vec::new();
        // Remember the last source error: if *every* source fails (e.g. the panel asset is missing),
        // surface it rather than silently returning an empty profile; a one-off per-source failure
        // (a chip with no stored raw file, an alignment lacking a BAM) is just skipped.
        let mut last_err: Option<AppError> = None;

        // One source per WGS alignment (panel-genotyped, cached per alignment). The IBD panel is
        // CHM13-coordinate and the alignment path has no liftover, so only CHM13 alignments can be
        // genotyped directly — non-CHM13 builds reach the panel via the chip path (multi-build
        // coordinates) or a future lift, and are skipped here rather than yielding wrong-locus calls.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            if !matches!(
                canonical_build(&a.reference_build),
                Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs)
            ) {
                continue;
            }
            match self.ibd_panel_dosages(IbdSource::Alignment(a.id)).await {
                Ok(gts) => {
                    let obs = to_obs(gts);
                    if !obs.is_empty() {
                        sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
                    }
                }
                Err(e) => last_err = Some(e),
            }
        }

        // One source per imported chip (resolved to canonical panel dosages, no alignment needed).
        let chips = self.list_chip_profiles(biosample_guid).await?;
        for c in &chips {
            match self.ibd_panel_dosages(IbdSource::Chip(c.id)).await {
                Ok(gts) => {
                    let obs = to_obs(gts);
                    if !obs.is_empty() {
                        sources.push((format!("{} (chip #{})", c.provider, c.id), SourceType::Chip, obs));
                    }
                }
                Err(e) => last_err = Some(e),
            }
        }

        if sources.is_empty() {
            if let Some(e) = last_err {
                return Err(e); // e.g. the IBD panel asset isn't built yet
            }
        }

        let source_summaries: Vec<YSourceSummary> = sources
            .iter()
            .map(|(label, st, obs)| YSourceSummary {
                label: label.clone(),
                source_type: *st,
                variant_count: obs.len(),
            })
            .collect();
        let variants = reconcile_diploid(&sources);
        let summary = summarize_diploid(&variants);
        let profile = DiploidProfile {
            variants,
            summary,
            sources: source_summaries,
        };

        // Persist (keyed dna_type='Auto'; no lineage label, no tree provider).
        self.persist_consensus_row(
            biosample_guid,
            "Auto",
            None,
            &profile.summary,
            profile.sources.len(),
            None,
            serde_json::to_string(&profile)?,
        )
        .await?;
        Ok(profile)
    }

    /// Build the `com.decodingus.atmosphere.haplogroupReconciliation` record JSON for a
    /// subject + DNA type from the stored consensus, per-run calls, manual override, and
    /// audit log. mtDNA heteroplasmy observations and an optional identity-verification
    /// result are passed in (the caller computes them from the relevant alignments).
    async fn reconciliation_record(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: &[HeteroplasmySite],
        identity: Option<&IdentityVerification>,
    ) -> Result<serde_json::Value, AppError> {
        let consensus = self
            .haplogroup_consensus(biosample_guid, dna_type)
            .await?
            .ok_or_else(|| {
                AppError::Store(StoreError::NotFound(format!(
                    "no {} haplogroup calls for {}",
                    dna_type.as_str(),
                    biosample_guid.0
                )))
            })?;
        let calls = haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?;

        let run_calls = calls
            .iter()
            .map(|c| RunHaplogroupCallRecord {
                source_ref: c.source_label.clone(),
                haplogroup: c.haplogroup.clone(),
                confidence: c.score.to_string(),
                call_method: "SNP_PHYLOGENETIC".into(),
                score: Some(c.score.to_string()),
                supporting_snps: Some(c.matched),
                conflicting_snps: Some((c.expected - c.matched).max(0)),
            })
            .collect();

        let status = ReconciliationStatusRecord {
            compatibility_level: compat_lexicon(consensus.compatibility).into(),
            consensus_haplogroup: consensus.haplogroup.clone(),
            confidence: Some(consensus.confidence.to_string()),
            divergence_point: consensus.divergence_point.clone(),
            branch_compatibility_score: None,
            snp_concordance: identity.and_then(|i| i.snp_concordance).map(|f| f.to_string()),
            run_count: consensus.run_count as i64,
            warnings: consensus.warnings.clone(),
        };

        // Heteroplasmy is mtDNA-only; major frequency is 1 − minor fraction.
        let heteroplasmy_observations = if dna_type == DnaType::Mt {
            heteroplasmy
                .iter()
                .map(|h| HeteroplasmyObservationRecord {
                    position: h.position,
                    major_allele: h.major_base.to_string(),
                    minor_allele: h.minor_base.to_string(),
                    major_allele_frequency: (1.0 - h.minor_fraction).to_string(),
                    depth: Some(h.depth as i64),
                    is_defining_snp: None,
                    affected_haplogroup: None,
                })
                .collect()
        } else {
            Vec::new()
        };

        let identity_verification = identity.map(|i| IdentityVerificationRecord {
            kinship_coefficient: None,
            fingerprint_snp_concordance: i.snp_concordance.map(|f| f.to_string()),
            y_str_distance: i.y_str_distance,
            verification_status: Some(verification_lexicon(i.status).into()),
            verification_method: Some(i.method.clone()),
        });

        let manual_override = recon_store::get_override(self.store.pool(), biosample_guid, dna_type)
            .await?
            .map(|(hg, reason)| ManualOverrideRecord {
                overridden_haplogroup: hg,
                reason,
                overridden_at: Utc::now().to_rfc3339(),
                overridden_by: self.current_account(),
            });

        let audit_log = self
            .reconciliation_audit(biosample_guid, dna_type)
            .await?
            .into_iter()
            .map(|a| AuditEntryRecord {
                timestamp: a.timestamp,
                action: a.action,
                previous_consensus: None,
                new_consensus: None,
                run_ref: None,
                notes: Some(a.note),
            })
            .collect();

        let record = HaplogroupReconciliationRecord::new(
            biosample_guid.0.to_string(),
            dna_type_lexicon(dna_type),
            Utc::now().to_rfc3339(),
            status,
            run_calls,
            heteroplasmy_observations,
            identity_verification,
            manual_override,
            audit_log,
        );
        Ok(serde_json::to_value(&record)?)
    }

    /// Publish a subject's haplogroup reconciliation using an explicit `client` (the
    /// testable core; production callers use [`publish_reconciliation`](Self::publish_reconciliation)).
    pub async fn publish_reconciliation_with(
        &self,
        client: &PdsClient,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: &[HeteroplasmySite],
        identity: Option<&IdentityVerification>,
    ) -> Result<RecordRef, AppError> {
        let value = self
            .reconciliation_record(biosample_guid, dna_type, heteroplasmy, identity)
            .await?;
        Ok(client
            .create_record(HAPLOGROUP_RECONCILIATION_COLLECTION, value, None)
            .await?)
    }

    /// Publish a subject's haplogroup reconciliation record to the signed-in account's PDS
    /// (with refresh-on-expiry and retry/backoff via [`AsyncSync`]).
    pub async fn publish_reconciliation(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: &[HeteroplasmySite],
        identity: Option<&IdentityVerification>,
    ) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self
            .reconciliation_record(biosample_guid, dna_type, heteroplasmy, identity)
            .await?;
        let entity_ref = format!("reconciliation:{biosample_guid}:{dna_type:?}");
        self.enqueue_publish(
            "reconciliation",
            &entity_ref,
            HAPLOGROUP_RECONCILIATION_COLLECTION,
            None,
            value,
        )
        .await
    }

    /// All recorded per-source calls for a subject + DNA type (for display / audit).
    pub async fn haplogroup_calls(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Vec<RunHaplogroupCall>, AppError> {
        Ok(haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?)
    }

    /// Like [`assign_mtdna_haplogroup`](Self::assign_mtdna_haplogroup) but with the tree
    /// JSON supplied directly (no network) — the testable core.
    pub async fn assign_mtdna_haplogroup_with_tree(
        &self,
        mtdna_id: i64,
        tree_json: &str,
    ) -> Result<HaploAssignment, AppError> {
        let seq = mtdna_store::get(self.store.pool(), mtdna_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("mtDNA sequence {mtdna_id}"))))?;

        // Sample base at each (rCRS-coordinate) position, straight from the full sequence.
        let mut calls: HashMap<i64, char> = HashMap::new();
        for (i, b) in seq.sequence.bytes().enumerate() {
            let u = b.to_ascii_uppercase();
            if matches!(u, b'A' | b'C' | b'G' | b'T') {
                calls.insert((i + 1) as i64, u as char);
            }
        }

        let tree = navigator_analysis::haplo::parse_ftdna_json(tree_json).map_err(AppError::Import)?;
        Ok(assemble_assignment(&tree, &calls))
    }

    /// FTDNA mt-DNA haplotree JSON, from the on-disk cache or freshly downloaded + cached.
    pub(crate) async fn fetch_ftdna_mt_tree(&self) -> Result<String, AppError> {
        self.fetch_tree(
            "https://www.familytreedna.com/public/mt-dna-haplotree/get",
            "ftdna-mttree.json",
        )
        .await
    }

    /// DecodingUs Y-DNA tree-with-variants JSON from our AppView (`/api/v1/y-tree/full`),
    /// host from [`decodingus_appview_url`]. On-disk cached like the FTDNA tree.
    pub(crate) async fn fetch_decodingus_y_tree(&self) -> Result<String, AppError> {
        let url = format!("{}/api/v1/y-tree/full", decodingus_appview_url());
        self.fetch_tree(&url, "decodingus-ytree.json").await
    }

    /// FTDNA Y-DNA haplotree JSON, from the on-disk cache or freshly downloaded + cached.
    pub(crate) async fn fetch_ftdna_y_tree(&self) -> Result<String, AppError> {
        self.fetch_tree(
            "https://www.familytreedna.com/public/y-dna-haplotree/get",
            "ftdna-ytree.json",
        )
        .await
    }

    /// The AppView's full instrument→lab map (`GET /api/v1/sequencer/lab-instruments`), on-disk
    /// cached like the trees (7-day TTL + offline fallback). Looked up locally so a batch import
    /// makes one network call, not one per sample.
    async fn fetch_lab_instruments(&self) -> Result<Vec<SequencerLabInfo>, AppError> {
        let url = format!("{}/api/v1/sequencer/lab-instruments", decodingus_appview_url());
        let json = self.fetch_tree(&url, "sequencer-lab-instruments.json").await?;
        serde_json::from_str(&json).map_err(|e| AppError::Import(format!("parsing lab-instruments: {e}")))
    }

    /// Resolve an instrument id to a lab display name via the AppView (cached). Normalizes the
    /// returned name to the local [`labs`] catalog's canonical display name when it matches.
    /// `None` if the instrument has no association or the AppView is unreachable (best-effort).
    pub async fn lookup_lab_by_instrument(&self, instrument_id: &str) -> Option<String> {
        let id = instrument_id.trim();
        if id.is_empty() {
            return None;
        }
        let list = self.fetch_lab_instruments().await.ok()?;
        let raw = list.into_iter().find(|l| l.instrument_id == id)?.lab_name;
        Some(
            navigator_domain::labs::find(&raw)
                .map(|l| l.display_name.to_string())
                .unwrap_or(raw),
        )
    }

    /// Resolve the FTDNA Big Y **generation** for a generic Targeted-Y run from its callable chrY
    /// footprint: a Big Y-500 covers ≤ ~10 Mb of callable chrY, and only the newer Big Y-700
    /// consistently exceeds it. Only acts on an FTDNA `TARGETED_Y` run — a header `@RG LB` label
    /// already pins the generation at import (those are `BIG_Y_500`/`BIG_Y_700`, never `TARGETED_Y`,
    /// so they're never second-guessed here), and a non-FTDNA targeted-Y stays generic. Idempotent.
    pub(crate) async fn refine_big_y_generation(&self, run: &SequenceRun, callable_chr_y: u64) -> Option<&'static str> {
        const BIG_Y_500_MAX_CALLABLE: u64 = 10_000_000;
        if run.test_type != "TARGETED_Y" {
            return None;
        }
        let is_ftdna = run
            .sequencing_facility
            .as_deref()
            .and_then(navigator_domain::labs::find)
            .map(|l| l.abbreviation)
            == Some("FTDNA");
        if !is_ftdna {
            return None;
        }
        let code = if callable_chr_y > BIG_Y_500_MAX_CALLABLE {
            "BIG_Y_700"
        } else {
            "BIG_Y_500"
        };
        let _ = sequence_run::set_test_type(self.store.pool(), run.id, code).await;
        Some(code)
    }

    /// [`Self::refine_big_y_generation`] keyed off an alignment's freshly computed (or cached)
    /// coverage — the callable-chrY base count is the discriminator. Called after coverage runs.
    /// Returns the new code when the generation changed (so the caller can refresh the run card).
    pub async fn refine_big_y_generation_for_alignment(
        &self,
        alignment_id: i64,
        coverage: &Coverage,
    ) -> Result<Option<&'static str>, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        if let Some(run) = sequence_run::get(self.store.pool(), aln.sequence_run_id).await? {
            return Ok(self.refine_big_y_generation(&run, callable_chr_y_bases(coverage)).await);
        }
        Ok(None)
    }

    /// Resolve the sequencing lab for every run that has an inferred `instrument_id` but no facility
    /// yet, via the AppView (one cached fetch). Best-effort; returns how many were filled. Run after
    /// import and on startup so pre-existing runs pick up newly-seeded associations.
    pub async fn backfill_run_labs(&self) -> Result<usize, AppError> {
        // One network/cache fetch (empty when offline — the FTDNA test-type normalization below is
        // local and still runs for runs whose facility was resolved earlier, e.g. subject 103589).
        let list = self.fetch_lab_instruments().await.unwrap_or_default();
        let by_instrument: HashMap<&str, &str> = list
            .iter()
            .map(|l| (l.instrument_id.as_str(), l.lab_name.as_str()))
            .collect();
        let mut filled = 0usize;
        for biosample in biosample::list_all(self.store.pool()).await? {
            for run in sequence_run::list_for_biosample(self.store.pool(), biosample.guid).await? {
                // Resolve + record the facility from the instrument→lab map when not already known.
                let facility = match run.sequencing_facility.clone() {
                    Some(f) => Some(f),
                    None => match run.instrument_id.as_deref().and_then(|i| by_instrument.get(i.trim())) {
                        Some(raw) => {
                            let lab = navigator_domain::labs::find(raw).map(|l| l.display_name).unwrap_or(raw);
                            if sequence_run::set_facility(self.store.pool(), run.id, lab)
                                .await
                                .unwrap_or(false)
                            {
                                filled += 1;
                            }
                            Some(lab.to_string())
                        }
                        None => None,
                    },
                };
                // A run we now know is FTDNA but typed as the generic TARGETED_Y is a Big Y — pick
                // its generation (500/700) from cached coverage when it's already been analyzed, so
                // pre-existing runs get corrected on startup without a re-analysis.
                let _ = facility; // (resolved above; the refine reads facility off the run record)
                if run.test_type == "TARGETED_Y" {
                    if let Ok(Some(cov)) = self.cached_coverage_for_run(run.id).await {
                        // Re-read the run so the just-set facility is visible to the FTDNA gate.
                        if let Ok(Some(fresh)) = sequence_run::get(self.store.pool(), run.id).await {
                            self.refine_big_y_generation(&fresh, callable_chr_y_bases(&cov)).await;
                        }
                    }
                }
            }
        }
        Ok(filled)
    }

    /// Cached coverage for a run, via its first alignment that has one (Big Y runs have a single
    /// alignment). `None` when the run hasn't been analyzed yet.
    async fn cached_coverage_for_run(&self, run_id: i64) -> Result<Option<Coverage>, AppError> {
        for aln in alignment::list_for_run(self.store.pool(), run_id).await? {
            if let Some(cov) = self.cached_coverage(aln.id).await? {
                return Ok(Some(cov));
            }
        }
        Ok(None)
    }

    /// A cached-or-downloaded haplotree JSON. The on-disk cache has a **7-day life** (see
    /// [`TREE_CACHE_TTL`]): a fresh cache short-circuits the network; a stale or missing cache
    /// triggers a re-download (and refresh). If the re-download fails (e.g. the AppView is
    /// unreachable) but a stale copy exists, the stale copy is used rather than failing — so the
    /// app keeps working offline, just on an older tree. (A server-side ETag/version would let us
    /// revalidate without a full re-download; tracked as an AppView backlog item.)
    async fn fetch_tree(&self, url: &str, cache_file: &str) -> Result<String, AppError> {
        let path = tree_cache_path(cache_file);
        let cached = std::fs::read_to_string(&path).ok().filter(|c| !c.trim().is_empty());
        if let Some(cached) = &cached {
            if tree_cache_is_fresh(&path) {
                return Ok(cached.clone());
            }
        }
        // Stale or absent → (re)download, falling back to a stale copy on network failure.
        let downloaded = self
            .auth
            .http
            .get(url)
            .send()
            .await
            .and_then(|r| r.error_for_status())
            .map_err(|e| AppError::Import(format!("downloading {url}: {e}")));
        match downloaded {
            Ok(resp) => {
                let body = resp
                    .text()
                    .await
                    .map_err(|e| AppError::Import(format!("reading {url}: {e}")))?;
                if let Some(parent) = path.parent() {
                    let _ = std::fs::create_dir_all(parent);
                }
                let _ = std::fs::write(&path, &body);
                Ok(body)
            }
            Err(e) => match cached {
                Some(stale) => {
                    eprintln!("tree refresh failed ({e}); using the cached copy at {}", path.display());
                    Ok(stale)
                }
                None => Err(e),
            },
        }
    }

    /// Assign an mtDNA haplogroup directly from an alignment's chrM reads (FTDNA mt tree),
    /// the BAM-based counterpart to [`assign_mtdna_haplogroup`]. Requires a GRCh38/rCRS
    /// chrM (the tree is in rCRS coordinates).
    pub async fn assign_mtdna_haplogroup_from_alignment(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();
        let source_key = format!("aln:{alignment_id}:mt");
        let tree_json = self.fetch_ftdna_mt_tree().await?;

        // Cache: skip re-scoring when the file and the mt tree are unchanged.
        let fingerprint = self
            .alignment_content_hash(alignment_id)
            .await
            .ok()
            .map(|file_hash| format!("f:{}|mt:{}", &file_hash[..16], &sha256_str(&tree_json)[..16]));
        if let (Some(bio), Some(fp)) = (bio, fingerprint.as_deref()) {
            if haplogroup_call::stored_fingerprint(self.store.pool(), bio, DnaType::Mt, &source_key)
                .await?
                .as_deref()
                == Some(fp)
            {
                if let Some(call) = haplogroup_call::get_one(self.store.pool(), bio, DnaType::Mt, &source_key).await? {
                    return Ok(assignment_from_call(&call));
                }
            }
        }

        let assignment = self
            .assign_haplogroup_from_alignment(alignment_id, "chrM", &tree_json)
            .await?;
        if let Some(bio) = bio {
            self.record_call_fp(
                bio,
                DnaType::Mt,
                &source_key,
                format!("aln #{alignment_id} mtDNA"),
                &assignment,
                fingerprint.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// mtDNA assignment + per-SNP lineage evidence (for exact GRCh38-vs-CHM13 comparison).
    pub async fn assign_mtdna_haplogroup_detail(
        &self,
        alignment_id: i64,
    ) -> Result<(HaploAssignment, Vec<SnpEvidence>, HashMap<i64, char>), AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        self.assign_haplogroup_detail(alignment_id, "chrM", &tree_json).await
    }

    /// Scan an alignment's chrM pileup for heteroplasmic positions — sites where a second
    /// mitochondrial allele coexists above the noise floor. A screening pass for the
    /// reconciliation view (a curator judges real heteroplasmy vs. artefacts); ascending
    /// by position. Requires a chrM-bearing BAM.
    pub async fn mtdna_heteroplasmy(&self, alignment_id: i64) -> Result<Vec<HeteroplasmySite>, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = aln.reference_path.map(PathBuf::from);
        tokio::task::spawn_blocking(move || {
            heteroplasmy::detect_heteroplasmy(&bam, "chrM", &HeteroplasmyParams::default(), reference.as_deref())
        })
        .await?
        .map_err(Into::into)
    }

    /// Estimate the donor's ancestry for an alignment by the allele-frequency likelihood: load
    /// the (build-matched) AIMs panel, genotype the sample at its sites with the GL caller, and
    /// score each super-population's binomial likelihood. Persists the result; returns it for
    /// display. Requires a recorded BAM/CRAM and a resolvable reference (CRAM/genotyping).
    /// Estimate autosomal ancestry from the subject's **consensus** — no BAM genotyping. Reads the
    /// cached autosomal [`DiploidProfile`] (reconciled 0/1/2 dosages over the probe panel, pooled
    /// across all WGS + chip sources), bridges it to genotypes, and runs the same estimators as the
    /// per-alignment path used to. Persisted under the consensus pseudo-source
    /// ([`CONSENSUS_SOURCE_ID`]). Errors if the autosomal consensus hasn't been built yet.
    pub async fn estimate_ancestry_from_consensus(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<AncestryResult, AppError> {
        let profile = self.cached_autosomal_profile(biosample_guid).await?.ok_or_else(|| {
            AppError::Import("build the autosomal consensus first (Autosomal tab) before estimating ancestry".into())
        })?;
        let genotypes = consensus_genotypes(&profile);

        // The consensus is canonical CHM13; the AIM freq / PCA assets are keyed by (contig,pos) there.
        let build = ReferenceBuild::Chm13v2;
        let reference_version = "chm13v2.0".to_string();
        let panel_path = ancestry_panel_path(build);
        let panel_bytes = read_verified_asset(build, &panel_path)?
            .ok_or_else(|| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        let optional = |path: PathBuf| {
            read_verified_asset(build, &path).unwrap_or_else(|e| {
                eprintln!("{e}");
                None
            })
        };
        let pca_bytes = optional(ancestry_pca_path(build));
        let ancient_pca_bytes = optional(ancestry_pca_ancient_path(build));
        let fine_bytes = optional(ancestry_freq_global_path(build));

        let (result, pca_gmm, nmonte, fine) = tokio::task::spawn_blocking(move || {
            let mut result = ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            let fine = fine_bytes
                .and_then(|b| ancestry_analysis::AncestryPanel::from_bytes(&b).ok())
                .map(|fp| ancestry_analysis::estimate_fine_admixture(&genotypes, &fp, &reference_version));
            let modern_pca = pca_bytes.and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok());
            if let Some(pca) = &modern_pca {
                result.pca_coordinates = Some(ancestry_analysis::project_pca(&genotypes, pca));
            }
            let gmm_pca = ancient_pca_bytes
                .and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok())
                .or(modern_pca);
            let (pca_gmm, nmonte) = match &gmm_pca {
                Some(pca) => (
                    Some(ancestry_analysis::estimate_pca_gmm(&genotypes, pca, &reference_version)),
                    Some(ancestry_analysis::estimate_nmonte(&genotypes, pca, &reference_version)),
                ),
                None => (None, None),
            };
            (result, pca_gmm, nmonte, fine)
        })
        .await?;

        let required = ancestry_min_snps();
        if result.snps_with_genotype < required {
            return Err(AppError::InsufficientAncestryData {
                genotyped: result.snps_with_genotype,
                required,
            });
        }
        ancestry_result::upsert(self.store.pool(), biosample_guid, CONSENSUS_SOURCE_ID, &result).await?;
        for extra in [pca_gmm.as_ref(), nmonte.as_ref(), fine.as_ref()].into_iter().flatten() {
            ancestry_result::upsert(self.store.pool(), biosample_guid, CONSENSUS_SOURCE_ID, extra).await?;
        }
        Ok(result)
    }

    /// The persisted ancestry estimate for an alignment, if one has been computed.
    pub async fn ancestry_for_alignment(&self, alignment_id: i64) -> Result<Option<AncestryResult>, AppError> {
        Ok(ancestry_result::get_for_alignment(self.store.pool(), alignment_id).await?)
    }

    /// The persisted **fine-population** admixture estimate for an alignment, if one was computed
    /// (the `ancestry_freq_global` asset was present at estimation time). Drives the super→fine
    /// hierarchy rows; the super-pop donut keeps using the primary ([`ancestry_for_alignment`]).
    pub async fn fine_ancestry_for_alignment(&self, alignment_id: i64) -> Result<Option<AncestryResult>, AppError> {
        Ok(ancestry_result::get_for_alignment_method(self.store.pool(), alignment_id, "FINE_ADMIXTURE").await?)
    }

    /// Reference population centroids on (PC1, PC2) for the alignment's build — the backdrop
    /// for the PCA scatter. `(population_code, pc1, pc2)`; empty if no PCA loadings are present.
    /// Reference population centroids in the **consensus** PC frame for the PCA scatter. The donor's
    /// projected coordinate (`AncestryResult::pca_coordinates`) is always computed against the CHM13
    /// PCA asset (the canonical consensus frame — see [`estimate_ancestry_from_consensus`]), so the
    /// backdrop centroids must come from that same asset regardless of which source is selected.
    /// Returns an empty vec when the asset isn't installed (the caller shows "reference not built").
    pub async fn ancestry_pca_reference(&self) -> Result<Vec<(String, f64, f64)>, AppError> {
        let build = ReferenceBuild::Chm13v2;
        let Ok(bytes) = std::fs::read(ancestry_pca_path(build)) else {
            return Ok(Vec::new());
        };
        let pca = navigator_analysis::ancestry::PcaLoadings::from_bytes(&bytes)?;
        Ok(pca
            .populations
            .iter()
            .enumerate()
            .map(|(p, code)| {
                let c = pca.centroid(p);
                (
                    code.clone(),
                    c.first().copied().unwrap_or(0.0) as f64,
                    c.get(1).copied().unwrap_or(0.0) as f64,
                )
            })
            .collect())
    }

    /// The cached chromosome painting for a subject, if one was painted from the **current** autosomal
    /// consensus (signature = the consensus's `last_reconciled_at`). `None` if absent or stale (the
    /// consensus was rebuilt since). Cheap — a cache read, no genotyping or HMM.
    pub async fn cached_painting(&self, biosample_guid: SampleGuid) -> Result<Option<Vec<AncestrySegment>>, AppError> {
        let Some(row) = consensus_profile::get(self.store.pool(), biosample_guid, "Auto").await? else {
            return Ok(None);
        };
        let Some(p) = consensus_painting::get(self.store.pool(), biosample_guid).await? else {
            return Ok(None);
        };
        if p.consensus_sig == row.last_reconciled_at {
            Ok(Some(serde_json::from_str(&p.segments)?))
        } else {
            Ok(None) // painted from an older consensus — stale
        }
    }

    /// Paint each chromosome with diploid local ancestry from the subject's **consensus** — no BAM
    /// walk. Returns the cached painting when it matches the current consensus signature; otherwise
    /// runs the diploid pair-state HMM over the consensus genotypes (anchored on the admixture prior)
    /// and caches it keyed to the consensus's `last_reconciled_at`.
    pub async fn paint_local_ancestry_from_consensus(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Vec<AncestrySegment>, AppError> {
        let row = consensus_profile::get(self.store.pool(), biosample_guid, "Auto")
            .await?
            .ok_or_else(|| {
                AppError::Import("build the autosomal consensus first (Autosomal tab) before painting".into())
            })?;
        let sig = row.last_reconciled_at.clone();

        // Cache hit (same consensus signature) → return without recomputing.
        if let Some(p) = consensus_painting::get(self.store.pool(), biosample_guid).await? {
            if p.consensus_sig == sig {
                return Ok(serde_json::from_str(&p.segments)?);
            }
        }

        let profile: DiploidProfile = serde_json::from_str(&row.payload)?;
        let genotypes = consensus_genotypes(&profile);
        let build = ReferenceBuild::Chm13v2;
        let reference_version = "chm13v2.0".to_string();
        let panel_path = ancestry_panel_path(build);
        let panel_bytes = read_verified_asset(build, &panel_path)?
            .ok_or_else(|| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        let segments = tokio::task::spawn_blocking(move || {
            let composition = ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            let prior: Vec<(String, f64)> = composition
                .components
                .iter()
                .map(|c| (c.population_code.clone(), c.percentage / 100.0))
                .collect();
            ancestry_analysis::paint_local_ancestry(
                &genotypes,
                &panel,
                &prior,
                &ancestry_analysis::PaintParams::default(),
            )
        })
        .await?;

        // Cache keyed to the consensus signature so it's reused until the consensus is rebuilt.
        consensus_painting::upsert(
            self.store.pool(),
            biosample_guid,
            &sig,
            &serde_json::to_string(&segments)?,
            &Utc::now().to_rfc3339(),
        )
        .await?;
        Ok(segments)
    }

    /// An alignment's content SHA-256, computed once at import. Read from the record if present,
    /// else computed now (hashing the file) and stored — so batch-imported alignments are hashed
    /// lazily on first analysis, then cached on the row.
    async fn alignment_content_hash(&self, alignment_id: i64) -> Result<String, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = aln.bam_path.clone();
        let hash = match aln.content_sha256.clone() {
            Some(h) => h,
            None => {
                let path = bam.clone().ok_or(AppError::MissingPaths(alignment_id))?;
                let h = sha256_file_async(PathBuf::from(path)).await?;
                let _ = alignment::set_content_hash(self.store.pool(), alignment_id, &h).await;
                h
            }
        };
        // Register the file by its content hash (gap §5-p2): stable identity across moves, and the
        // dedup/accessibility registry. Idempotent — a moved file just updates its path here.
        if let Some(path) = bam {
            let size = std::fs::metadata(&path).ok().map(|m| m.len() as i64);
            let now = Utc::now().to_rfc3339();
            if source_file::upsert_by_checksum(self.store.pool(), &hash, Some(&path), size, Some("BAM"), &now)
                .await
                .is_ok()
            {
                let _ = source_file::link_to_alignment(self.store.pool(), &hash, alignment_id, &now).await;
            }
        }
        Ok(hash)
    }

    /// All tracked source files (content-hash identity) — for the Data Sources view.
    pub async fn list_source_files(&self) -> Result<Vec<SourceFile>, AppError> {
        Ok(source_file::list(self.store.pool()).await?)
    }

    /// Re-check each tracked file's path on disk and update its accessibility flag. Returns how many
    /// are now missing (moved/deleted) — surfaced as a "file missing" marker in the UI.
    pub async fn verify_source_files(&self) -> Result<usize, AppError> {
        let now = Utc::now().to_rfc3339();
        let mut missing = 0;
        for f in source_file::list(self.store.pool()).await? {
            let ok = f
                .file_path
                .as_deref()
                .map(|p| std::path::Path::new(p).exists())
                .unwrap_or(false);
            if !ok {
                missing += 1;
            }
            if ok != f.is_accessible {
                let _ = source_file::set_accessible(self.store.pool(), f.id, ok, &now).await;
            }
        }
        Ok(missing)
    }

    /// Fingerprint of the inputs to a Y-haplogroup score: the alignment's content hash + the
    /// active Y-tree's content hash. Unchanged inputs → a re-score is unnecessary. Errors (e.g.
    /// the tree is unreachable and uncached) disable caching for this run rather than failing.
    async fn y_score_fingerprint(&self, alignment_id: i64) -> Result<String, AppError> {
        let file_hash = self.alignment_content_hash(alignment_id).await?;
        let tree_json = match y_tree_provider() {
            YTreeProvider::DecodingUs => self.fetch_decodingus_y_tree().await?,
            YTreeProvider::Ftdna => self.fetch_ftdna_y_tree().await?,
        };
        let tree_hash = sha256_str(&tree_json);
        Ok(format!("f:{}|yt:{}", &file_hash[..16], &tree_hash[..16]))
    }

    /// Assign a Y haplogroup to an alignment: place the sample against the configured Y tree
    /// (DecodingUs by default — our tree, native CHM13 coords, no liftover — falling back to
    /// FTDNA if the AppView is unreachable), call the sample's base at each tree position on
    /// chrY, and rank by Kulczynski. Requires a recorded BAM/CRAM path. Skips re-scoring when
    /// the alignment file and tree are unchanged since the last run (see [`Self::y_score_fingerprint`]).
    pub async fn assign_y_haplogroup(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();
        let source_key = format!("aln:{alignment_id}");

        // Input fingerprint = alignment content hash + active Y-tree content hash. If it matches
        // the recorded call's stamp, neither the file nor the tree changed → return the recorded
        // call without re-scoring (the expensive BAM genotyping).
        let fingerprint = self.y_score_fingerprint(alignment_id).await.ok();
        if let (Some(bio), Some(fp)) = (bio, fingerprint.as_deref()) {
            if haplogroup_call::stored_fingerprint(self.store.pool(), bio, DnaType::Y, &source_key)
                .await?
                .as_deref()
                == Some(fp)
            {
                if let Some(call) = haplogroup_call::get_one(self.store.pool(), bio, DnaType::Y, &source_key).await? {
                    return Ok(assignment_from_call(&call));
                }
            }
        }

        let assignment = self.y_assignment_full(alignment_id).await?;
        if let Some(bio) = bio {
            self.record_call_fp(
                bio,
                DnaType::Y,
                &source_key,
                format!("aln #{alignment_id} Y"),
                &assignment,
                fingerprint.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Freshly place an alignment against the configured Y tree, returning the **full** assignment
    /// **including per-branch SNP evidence** (the cached [`assign_y_haplogroup`] path returns only
    /// the terminal). Expensive (genotypes chrY tree sites in the BAM) — used by the Y-variant
    /// profile, which the user builds explicitly.
    async fn y_assignment_full(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        match y_tree_provider() {
            YTreeProvider::DecodingUs => match self.assign_y_decodingus(alignment_id).await {
                Ok(a) => Ok(a),
                Err(e) => {
                    // AppView unreachable / build unsupported / parse failure → FTDNA fallback.
                    eprintln!("DecodingUs Y tree unavailable ({e}); falling back to FTDNA");
                    let tree_json = self.fetch_ftdna_y_tree().await?;
                    self.assign_haplogroup_from_alignment(alignment_id, "chrY", &tree_json)
                        .await
                }
            },
            YTreeProvider::Ftdna => {
                let tree_json = self.fetch_ftdna_y_tree().await?;
                self.assign_haplogroup_from_alignment(alignment_id, "chrY", &tree_json)
                    .await
            }
        }
    }

    /// Place against the DecodingUs Y tree from our AppView, using the alignment's **native**
    /// build coordinates (`hs1` for CHM13, `GRCh38`, `GRCh37`) — queried directly, **no
    /// liftover**. This is the intended architecture (the AppView owns multi-build coordinates;
    /// Navigator stays liftover-free). Today the AppView's `hs1` coords cover the decoding-us
    /// backbone but not the FTDNA-grafted tips, so deep CHM13 placement is limited until the
    /// AppView enriches `hs1` for every variant (lift GRCh38→hs1 at ingest or on the fly).
    async fn assign_y_decodingus(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let (tree, calls) = self.y_decodingus_tree_calls(alignment_id).await?;
        Ok(assemble_assignment(&tree, &calls))
    }

    /// The (DecodingUs tree at the alignment's **native** build, full tree-locus base calls) for one
    /// alignment — the genotype [`assign_y_decodingus`] scores. Factored so the consensus pool can
    /// re-key it by SNP name and merge it with other sources.
    async fn y_decodingus_tree_calls(
        &self,
        alignment_id: i64,
    ) -> Result<(navigator_analysis::haplo::HaploTree, HashMap<i64, char>), AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let build_key = decodingus_build_key(&aln.reference_build).ok_or_else(|| {
            AppError::Import(format!(
                "no DecodingUs tree coordinates for build {}",
                aln.reference_build
            ))
        })?;
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_decodingus_json(&tree_json, build_key).map_err(AppError::Import)?;
        // Native build → no liftover (tree_source_build = None → direct query).
        let calls = self.base_calls(alignment_id, "chrY", &tree, None).await?;
        Ok((tree, calls))
    }

    /// Assign a Y haplogroup from the subject's imported **BISDNA / Y-SNP-panel** calls — no
    /// alignment required. Builds a derived-allele call map from the subject's `Chip`-sourced
    /// variant sets (the panel's positive calls, each `position → derived base`) and scores it
    /// against the Y tree on `build` (the subject's alignment build, else `"hs1"`). Uses the
    /// DecodingUs tree at the native build (FTDNA fallback only on GRCh38, where positions
    /// match), and the chip-robust terminal selection ([`assemble_assignment_robust`]). The
    /// call is recorded as a reconciliation source. Only derived (positive) calls drive the
    /// Kulczynski ranking, so the stored positives-only variant set is sufficient.
    pub async fn assign_y_bisdna(
        &self,
        biosample_guid: SampleGuid,
        build: Option<&str>,
    ) -> Result<HaploAssignment, AppError> {
        // Derived-allele calls from the subject's chip-sourced variant sets (BISDNA positives).
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;

        // Placement build: explicit override, else the build stored on a chip set at import,
        // else (pre-migration sets with no stored build) re-derive from the subject's alignment.
        let build = match build {
            Some(b) => b.to_string(),
            None => match sets
                .iter()
                .filter(|s| s.source_type == SourceType::Chip)
                .find_map(|s| s.reference_build.clone())
            {
                Some(b) => b,
                None => self.bisdna_target_build(biosample_guid).await,
            },
        };

        let mut calls: HashMap<i64, char> = HashMap::new();
        for s in &sets {
            if s.source_type != SourceType::Chip {
                continue;
            }
            for c in &s.calls {
                if !c.contig.eq_ignore_ascii_case("chrY") && !c.contig.eq_ignore_ascii_case("y") {
                    continue;
                }
                if let Some(b) = c.alternate.chars().next() {
                    calls.insert(c.position, b.to_ascii_uppercase());
                }
            }
        }
        if calls.is_empty() {
            return Err(AppError::Import(
                "no Y-SNP panel calls to place — import a BISDNA file for this subject first".into(),
            ));
        }

        let tree = self.chip_y_tree(&build).await?;
        // Chip alleles (BISDNA + consumer arrays) are plus-strand; flip the minority recorded on
        // the tree's opposite strand so they score against the right allele. No-op for BISDNA.
        let calls = strand_reconcile_to_tree(&tree, calls);
        let assignment = assemble_assignment_robust(&tree, &calls);
        self.record_call(
            biosample_guid,
            DnaType::Y,
            "bisdna",
            "Chip Y-SNP panel".into(),
            &assignment,
        )
        .await?;
        Ok(assignment)
    }

    /// chrY genotype calls (`position → uppercase ALT base`) from a variant set — the shared
    /// extractor behind the chip / vendor-VCF Y placements.
    fn vset_chr_y_calls(set: &VariantSet) -> HashMap<i64, char> {
        set.calls
            .iter()
            .filter(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"))
            .filter_map(|c| c.alternate.chars().next().map(|b| (c.position, b.to_ascii_uppercase())))
            .collect()
    }

    /// Place the subject's vendor **Y-NGS VCF** variant sets — FTDNA Big Y / Full Genomes Y Elite /
    /// YSEQ / Nebula / Dante, i.e. anything imported as a non-[`Chip`](SourceType::Chip)
    /// [`VariantSet`] carrying chrY calls — and record a per-source donor call for each. These are
    /// direct Y-SNP genotype calls (the gold-standard placement input), placed against the configured
    /// tree on each set's stored build (FTDNA Big Y is GRCh38, the FTDNA tree's native build → no
    /// liftover). Best-effort per set: one that errors or lacks chrY calls is skipped. Returns the
    /// number of sets placed. Called on import (so a Big Y VCF places without a manual Refresh) and
    /// re-runnable. The vendor-VCF counterpart to [`assign_y_bisdna`](Self::assign_y_bisdna).
    pub async fn assign_y_vendor_vcfs(&self, biosample_guid: SampleGuid) -> Result<usize, AppError> {
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut tree_cache: HashMap<String, navigator_analysis::haplo::HaploTree> = HashMap::new();
        let mut placed = 0;
        for set in &sets {
            if set.source_type == SourceType::Chip {
                continue; // chips place via assign_y_bisdna / the chip-panel path
            }
            let calls = Self::vset_chr_y_calls(set);
            if calls.is_empty() {
                continue;
            }
            let build = set.reference_build.clone().unwrap_or_else(|| "GRCh38".to_string());
            if !tree_cache.contains_key(&build) {
                match self.chip_y_tree(&build).await {
                    Ok(t) => {
                        tree_cache.insert(build.clone(), t);
                    }
                    Err(e) => {
                        eprintln!("vendor Y-VCF placement deferred for build {build} ({e})");
                        continue;
                    }
                }
            }
            let assignment = Self::place_chip_panel(&tree_cache[&build], calls);
            self.record_call(
                biosample_guid,
                DnaType::Y,
                &format!("vcf#{}", set.id),
                set.source_label.clone(),
                &assignment,
            )
            .await?;
            placed += 1;
        }
        Ok(placed)
    }

    /// Fetch + parse the Y haplotree for a chip placement on `build`. DecodingUs is native multi-build
    /// (no liftover); the FTDNA tree is GRCh38-only, so it's a fallback only when the calls are GRCh38.
    /// Shared by the combined [`assign_y_bisdna`](Self::assign_y_bisdna) placement and the per-panel
    /// Y-profile sources, so the tree is fetched once.
    pub(crate) async fn chip_y_tree(&self, build: &str) -> Result<navigator_analysis::haplo::HaploTree, AppError> {
        // Honor the configured Y-tree provider (the alignment placement path does too). With the
        // FTDNA provider, place against the FTDNA tree directly — no DecodingUs call. The default
        // (DecodingUs) keeps the prior behavior, with an FTDNA fallback for a GRCh38 chip when the
        // DecodingUs tree is unavailable.
        if matches!(y_tree_provider(), YTreeProvider::Ftdna) {
            let json = self.fetch_ftdna_y_tree().await?;
            let mut tree = navigator_analysis::haplo::parse_ftdna_json(&json).map_err(AppError::Import)?;
            // Repair FTDNA's reference-as-ancestral polarity against DecodingUs (best effort).
            if let Some(pol) = self.decodingus_y_polarity().await {
                navigator_analysis::haplo::normalize_polarity(&mut tree, &pol);
            }
            return Ok(tree);
        }
        match self.fetch_decodingus_y_tree().await {
            Ok(json) => navigator_analysis::haplo::parse_decodingus_json(&json, build).map_err(AppError::Import),
            Err(e) if build == "GRCh38" => {
                eprintln!("DecodingUs Y tree unavailable ({e}); falling back to FTDNA (GRCh38)");
                let json = self.fetch_ftdna_y_tree().await?;
                navigator_analysis::haplo::parse_ftdna_json(&json).map_err(AppError::Import)
            }
            Err(e) => Err(e),
        }
    }

    /// Place one chip/BISDNA panel's chrY calls on `tree` (strand-reconciled), without persisting —
    /// for assembling the per-panel sources of the Y-variant profile.
    fn place_chip_panel(tree: &navigator_analysis::haplo::HaploTree, calls: HashMap<i64, char>) -> HaploAssignment {
        let calls = strand_reconcile_to_tree(tree, calls);
        assemble_assignment_robust(tree, &calls)
    }

    /// Place an mtDNA haplogroup from the subject's chip-sourced MT genotype calls (e.g. 23andMe
    /// `MT` rows) against the FTDNA mt tree. Consumer-array MT positions are rCRS coordinates,
    /// which the tree uses directly (no liftover). Reads every `Chip`-source variant set's chrM
    /// calls, reconciles strand, and uses the robust (sparse-chip) terminal selection. Records a
    /// donor call. The counterpart to [`assign_y_bisdna`](Self::assign_y_bisdna) for mtDNA.
    pub async fn assign_mt_chip(&self, biosample_guid: SampleGuid) -> Result<HaploAssignment, AppError> {
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut calls: HashMap<i64, char> = HashMap::new();
        for s in &sets {
            if s.source_type != SourceType::Chip {
                continue;
            }
            for c in &s.calls {
                let mt = c.contig.eq_ignore_ascii_case("chrM")
                    || c.contig.eq_ignore_ascii_case("chrMT")
                    || c.contig.eq_ignore_ascii_case("mt")
                    || c.contig.eq_ignore_ascii_case("m");
                if !mt {
                    continue;
                }
                if let Some(b) = c.alternate.chars().next() {
                    calls.insert(c.position, b.to_ascii_uppercase());
                }
            }
        }
        if calls.is_empty() {
            return Err(AppError::Import(
                "no chip mtDNA calls to place — import a 23andMe raw-data file for this subject first".into(),
            ));
        }
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;
        let calls = strand_reconcile_to_tree(&tree, calls);
        let assignment = assemble_assignment_robust(&tree, &calls);
        self.record_call(
            biosample_guid,
            DnaType::Mt,
            "chip-mt",
            "Chip mtDNA panel".into(),
            &assignment,
        )
        .await?;
        Ok(assignment)
    }

    /// Genotype an alignment at a haplotree's positions on `contig` and rank haplogroups by
    /// the Kulczynski measure. The networkless core shared by [`assign_y_haplogroup`] (also
    /// directly testable with a local tree + contig).
    pub async fn assign_haplogroup_from_alignment(
        &self,
        alignment_id: i64,
        contig: &str,
        tree_json: &str,
    ) -> Result<HaploAssignment, AppError> {
        let (tree, calls) = self.tree_base_calls(alignment_id, contig, tree_json).await?;
        Ok(assemble_assignment(&tree, &calls))
    }

    /// Like [`assign_haplogroup_from_alignment`], but also returns the per-SNP evidence along
    /// the called terminal's lineage (each defining mutation's Derived/Ancestral/NoCall state).
    /// For exact comparisons (e.g. GRCh38 vs a lifted CHM13 call).
    /// Full Y-haplogroup placement **report** for an alignment (gap §8): the ranked candidate
    /// haplogroups (with score / matched-vs-expected) + the defining-SNP evidence along the reported
    /// lineage (each SNP's derived / ancestral / no-call state). A fresh placement against the
    /// configured provider tree — heavier than the cached terminal label, so it's button-driven.
    pub async fn y_haplogroup_report(
        &self,
        alignment_id: i64,
    ) -> Result<(HaploAssignment, Vec<SnpEvidence>), AppError> {
        let tree_json = match y_tree_provider() {
            YTreeProvider::DecodingUs => self.fetch_decodingus_y_tree().await?,
            YTreeProvider::Ftdna => self.fetch_ftdna_y_tree().await?,
        };
        let (assignment, lineage, _calls) = self.assign_haplogroup_detail(alignment_id, "chrY", &tree_json).await?;
        Ok((assignment, lineage))
    }

    pub async fn assign_haplogroup_detail(
        &self,
        alignment_id: i64,
        contig: &str,
        tree_json: &str,
    ) -> Result<(HaploAssignment, Vec<SnpEvidence>, HashMap<i64, char>), AppError> {
        let (tree, calls) = self.tree_base_calls(alignment_id, contig, tree_json).await?;
        let assignment = assemble_assignment(&tree, &calls);
        let lineage = match assignment.ranked.first() {
            Some(top) => navigator_analysis::haplo::lineage_evidence(&tree, &calls, top.id),
            None => Vec::new(),
        };
        Ok((assignment, lineage, calls))
    }

    /// Parse the tree, build the per-position base calls (lifting onto the alignment's build
    /// when needed), and return both. Shared by the assignment + detail entry points.
    async fn tree_base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        tree_json: &str,
    ) -> Result<(navigator_analysis::haplo::HaploTree, HashMap<i64, char>), AppError> {
        let mut tree = navigator_analysis::haplo::parse_ftdna_json(tree_json).map_err(AppError::Import)?;
        // Harden the FTDNA tree: it records the GRCh38 *reference* base as "ancestral", so at the
        // sites where the reference carries the derived allele its polarity is inverted (the
        // CT-M168 amber artifact). Normalize against the DecodingUs tree's true polarity (best
        // effort — offline FTDNA mode keeps the raw FTDNA polarity). Y only (no mt polarity source).
        if contig.eq_ignore_ascii_case("chrY") {
            if let Some(pol) = self.decodingus_y_polarity().await {
                let flipped = navigator_analysis::haplo::normalize_polarity(&mut tree, &pol);
                if flipped > 0 {
                    eprintln!("normalized {flipped} FTDNA Y loci to DecodingUs (true) polarity");
                }
            }
        }
        // FTDNA tree positions are in the tree's own build (Y → GRCh38, mt → rCRS/direct).
        let source_build = tree_build_for_contig(contig);
        let calls = self.base_calls(alignment_id, contig, &tree, source_build).await?;
        Ok((tree, calls))
    }

    /// Best-effort true-polarity map (SNP name → ancestral/derived) from the DecodingUs Y tree,
    /// used to repair an FTDNA tree's reference-as-ancestral polarity ([`normalize_polarity`]).
    /// `None` when the DecodingUs tree is unavailable (offline FTDNA mode keeps raw FTDNA polarity).
    async fn decodingus_y_polarity(&self) -> Option<HashMap<String, (String, String)>> {
        let json = self.fetch_decodingus_y_tree().await.ok()?;
        navigator_analysis::haplo::decodingus_polarity_map(&json).ok()
    }

    /// Resolve a canonical contig name (`chrY`, `chrM`) to the name actually present in the
    /// alignment header, tolerating naming conventions: the `chr` prefix (GRCh37/hg19 drop it)
    /// and the `M`/`MT` mitochondrial spelling. Returns `None` when no equivalent contig is in
    /// the header (the caller then queries the requested name and surfaces the original error).
    async fn resolve_header_contig(
        &self,
        bam: &Path,
        reference: Option<&Path>,
        contig: &str,
    ) -> Result<Option<String>, AppError> {
        let bam = bam.to_path_buf();
        let reference = reference.map(|p| p.to_path_buf());
        let names = tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, reference.as_deref()))
            .await??;
        // Candidate spellings for the requested contig, in preference order.
        let bare = contig.strip_prefix("chr").unwrap_or(contig);
        let mut candidates: Vec<String> = vec![contig.to_string(), bare.to_string()];
        if bare.eq_ignore_ascii_case("M") || bare.eq_ignore_ascii_case("MT") {
            for alt in ["chrM", "chrMT", "M", "MT"] {
                candidates.push(alt.to_string());
            }
        }
        Ok(candidates
            .into_iter()
            .find(|cand| names.iter().any(|n| n.eq_ignore_ascii_case(cand)))
            .and_then(|cand| {
                // Return the header's exact casing/spelling so the region query matches.
                names.iter().find(|n| n.eq_ignore_ascii_case(&cand)).cloned()
            }))
    }

    /// Base-call an alignment at a parsed tree's positions on `contig`. `tree_source_build` is
    /// the build the tree's positions are in: when it differs from the alignment build the
    /// positions are lifted (chrY chain), queried there, and mapped back; `None` (e.g. a
    /// DecodingUs tree already in the alignment's build, or mt/rCRS-direct) queries directly.
    async fn base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        // Resolve the reference even when none was stored at import. A CRAM can't be decoded
        // without it, so resolve (download on a miss) via the gateway from the alignment's build —
        // e.g. the already-cached `chm13v2.0.fa` for a CHM13 CRAM. A BAM needs no reference to
        // read reads, so only adopt a *cached* FASTA (never force a multi-GB download just to
        // supply the chrM liftover map).
        let is_cram = bam.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        let reference = match aln.reference_path {
            Some(p) => Some(PathBuf::from(p)),
            None if is_cram => Some(
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?,
            ),
            None => self.gateway.cached_reference(&aln.reference_build),
        };

        let targets: HashSet<i64> = tree
            .nodes
            .values()
            .flat_map(|n| n.loci.iter().map(|l| l.position))
            .collect();

        let lifted = self
            .lifted_targets(
                &aln.reference_build,
                reference.as_deref(),
                contig,
                &targets,
                tree_source_build,
            )
            .await?;

        let calls = match lifted {
            Some(lifted) => self.build_calls_from_lifted(&bam, reference.as_deref(), lifted).await?,
            None => {
                // Match the requested contig to the header's naming convention: GRCh37/hg19
                // (still the medical-space default) use bare `Y`/`MT`; CHM13/GRCh38 use
                // `chrY`/`chrM`. Without this a `chrY` query against a `Y`-named header errors
                // out and the placement falls back to a worse tree.
                let resolved = self
                    .resolve_header_contig(&bam, reference.as_deref(), contig)
                    .await?
                    .unwrap_or_else(|| contig.to_string());
                let bam = bam.clone();
                let reference = reference.clone();
                let targets = targets.clone();
                tokio::task::spawn_blocking(move || {
                    let params = adaptive_haploid_params(&bam, reference.as_deref()); // HiFi -> lower min_depth
                    caller::call_bases_at(&bam, &resolved, &targets, &params, reference.as_deref())
                })
                .await??
            }
        };
        Ok(calls)
    }

    /// Lift the haplotree's positions onto the alignment's build, or `None` to query the tree
    /// positions directly. **chrY**: uses the (auto-downloaded) GRCh38→build liftover chain.
    /// **chrM**: a self-generated rCRS↔`chrM` map — bundled rCRS aligned to *this* reference's
    /// `chrM` (CHM13 builds only; GRCh38/rCRS `chrM` is already rCRS → direct).
    pub(crate) async fn lifted_targets(
        &self,
        reference_build: &str,
        reference: Option<&Path>,
        contig: &str,
        targets: &HashSet<i64>,
        tree_source_build: Option<&str>,
    ) -> Result<Option<Vec<LiftedPos>>, AppError> {
        if targets.is_empty() {
            return Ok(None);
        }

        // chrY: downloaded nuclear chain (when the tree build differs from the alignment).
        if let Some(src) = tree_source_build {
            let differ =
                matches!((canonical_build(src), canonical_build(reference_build)), (Some(s), Some(t)) if s != t);
            if differ && self.gateway.chain_available(src, reference_build) {
                self.gateway.resolve_chain(src, reference_build, &mut |_, _| {}).await?;
                let targets_vec: Vec<i64> = targets.iter().copied().collect();
                return Ok(Some(self.gateway.lift_positions(
                    src,
                    reference_build,
                    contig,
                    &targets_vec,
                )?));
            }
            return Ok(None);
        }

        // chrM on CHM13: self-generated rCRS↔chrM alignment map (no chain exists).
        if contig.eq_ignore_ascii_case("chrM") && canonical_build(reference_build) == Some(ReferenceBuild::Chm13v2) {
            let Some(reference) = reference else { return Ok(None) };
            let reference = reference.to_path_buf();
            // Align bundled rCRS to this reference's chrM (cheap, ~16.5 kb) → (rcrs, chrM) pairs.
            let map = tokio::task::spawn_blocking(move || {
                navigator_analysis::reader::read_contig_sequence(&reference, "chrM").map(|chrm| {
                    let chrm = String::from_utf8_lossy(&chrm).into_owned();
                    // Rotation-aware: CHM13's chrM is a circular permutation of rCRS.
                    navigator_analysis::mtvariants::mt_position_map(navigator_analysis::mtvariants::rcrs(), &chrm)
                })
            })
            .await?;
            let Ok(pairs) = map else { return Ok(None) }; // chrM absent/unreadable → direct fallback
                                                          // rcrs_idx/chrm_idx are 0-based; tree + query positions are 1-based.
            let by_rcrs: HashMap<i64, i64> = pairs.into_iter().map(|(r, c)| (r as i64 + 1, c as i64 + 1)).collect();
            let lifted = targets
                .iter()
                .filter_map(|&t| {
                    by_rcrs.get(&t).map(|&c| LiftedPos {
                        tree_pos: t,
                        contig: "chrM".to_string(),
                        pos: c,
                        reverse: false,
                    })
                })
                .collect();
            return Ok(Some(lifted));
        }

        Ok(None)
    }

    /// Query the already-lifted positions and map observed bases back to the original tree
    /// positions so [`assemble_assignment`] (which keys on tree positions) scores unchanged.
    /// Queries each lifted contig present in the BAM header; minus-strand lifts are
    /// reverse-complemented.
    async fn build_calls_from_lifted(
        &self,
        bam: &Path,
        reference: Option<&Path>,
        lifted: Vec<LiftedPos>,
    ) -> Result<HashMap<i64, char>, AppError> {
        // Group lifted positions by their target contig + a back-map (lifted → tree position,
        // plus whether the lift was to the minus strand → the base needs complementing).
        let mut by_contig: HashMap<String, HashSet<i64>> = HashMap::new();
        let mut back: HashMap<(String, i64), (i64, bool)> = HashMap::new();
        for lp in lifted {
            by_contig.entry(lp.contig.clone()).or_default().insert(lp.pos);
            back.insert((lp.contig, lp.pos), (lp.tree_pos, lp.reverse));
        }

        // Only query contigs the alignment actually has (drop off-target lifts).
        let header_contigs: HashSet<String> = {
            let bam = bam.to_path_buf();
            let reference = reference.map(|p| p.to_path_buf());
            tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, reference.as_deref()))
                .await??
                .into_iter()
                .collect()
        };

        let mut calls: HashMap<i64, char> = HashMap::new();
        for (qcontig, set) in by_contig {
            // Tolerate naming conventions between the lift target and the header (e.g. a `chrY`
            // lift against a GRCh37 `Y`-named header): query the header's actual spelling.
            let bare = qcontig.strip_prefix("chr").unwrap_or(&qcontig);
            let Some(query_contig) = header_contigs
                .iter()
                .find(|n| n.eq_ignore_ascii_case(&qcontig) || n.eq_ignore_ascii_case(bare))
                .cloned()
            else {
                continue; // off-target lift the alignment lacks
            };
            let bam = bam.to_path_buf();
            let reference = reference.map(|p| p.to_path_buf());
            let qc = query_contig;
            let lifted_calls = tokio::task::spawn_blocking(move || {
                let params = adaptive_haploid_params(&bam, reference.as_deref());
                caller::call_bases_at(&bam, &qc, &set, &params, reference.as_deref())
            })
            .await??;
            for (lpos, base) in lifted_calls {
                if let Some(&(tree_pos, reverse)) = back.get(&(qcontig.clone(), lpos)) {
                    // Inverted tracts (common on the CHM13 Y): the tree allele is GRCh38-forward,
                    // so complement the base read off the minus-strand-lifted CHM13 position.
                    calls.insert(tree_pos, if reverse { complement_base(base) } else { base });
                }
            }
        }
        Ok(calls)
    }
}
