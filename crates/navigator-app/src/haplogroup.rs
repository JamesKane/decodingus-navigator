//! `impl App` methods extracted from `lib.rs` (the `haplogroup` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;
use crate::fastpath::{chr_m_gvcf_for_alignment, chr_y_gvcf_for_alignment};

/// Analysis-artifact `kind` for cached per-alignment tree-genotype base calls (see
/// [`App::base_calls`]). The `algorithm_version` carries the site-set hash, so distinct trees /
/// contigs / lift paths get distinct cache rows.
const GENOTYPE_KIND: &str = "tree-genotype";

/// Session memo of fetched haplotree JSON (resolved cache path → body), shared by [`App::fetch_tree`] and
/// cleared by [`App::refresh_trees`]. Trees are 4–121 MB and consulted many times per placement, so
/// they're resolved at most once per process; a corrected AppView tree is picked up by clearing this
/// + the on-disk cache and re-fetching.
static TREE_MEMO: std::sync::OnceLock<std::sync::Mutex<HashMap<String, String>>> = std::sync::OnceLock::new();
fn tree_memo() -> &'static std::sync::Mutex<HashMap<String, String>> {
    TREE_MEMO.get_or_init(|| std::sync::Mutex::new(HashMap::new()))
}

/// Per-source weighted genotype calls (position → base) for one build's pooling group in
/// [`App::place_y_consensus`]: one `(source type, position→base)` entry per contributing source.
type YSourceCalls = Vec<(SourceType, HashMap<i64, char>)>;

/// A stable cache key (used as the artifact `algorithm_version`) for a tree-genotype base call:
/// the queried `contig`, the lift source build, and an FNV-1a hash of the **sorted target
/// positions** + their count. A changed tree (added/removed/moved positions) changes the hash →
/// cache miss → fresh walk; the BAM `source_sig` handles a changed alignment file separately.
fn genotype_cache_key(contig: &str, source_build: Option<&str>, targets: &HashSet<i64>) -> String {
    let mut sorted: Vec<i64> = targets.iter().copied().collect();
    sorted.sort_unstable();
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |bytes: &[u8]| {
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(0x100_0000_01b3);
        }
    };
    feed(contig.as_bytes());
    feed(b"|");
    feed(source_build.unwrap_or("native").as_bytes());
    feed(b"|");
    for p in &sorted {
        feed(&p.to_le_bytes());
    }
    // `g3`: chrY native genotyping now also resolves indel loci (additive derived sentinels), so the
    // cached result differs from the SNP-only `g1` payload — bump on any genotyping-logic change so a
    // stale payload isn't reused (the site-set hash alone doesn't capture logic changes).
    format!("g3:{contig}:{}:{h:016x}", sorted.len())
}

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
            ExportRequest::DiploidVcf(id) => self.diploid_vcf_genome(*id, navigator_analysis::CancelToken::none()).await,
            ExportRequest::ConsensusDiploidVcf(guid) => self.consensus_diploid_vcf(*guid).await,
            ExportRequest::SubjectBriefHtml(guid) => {
                let brief = self.subject_brief(*guid).await?;
                // Fold in the AI story if one is already cached (no generation during an export).
                let narration = self.cached_narration(&brief);
                Ok(export::subject_brief_html(&brief, narration.as_ref()))
            }
            ExportRequest::DescentTsv(guid, dna) => {
                let report = self.descent_report(*guid, *dna).await?.ok_or_else(|| {
                    AppError::Store(StoreError::NotFound(
                        "no descent report yet — build the variant profile first".into(),
                    ))
                })?;
                Ok(export::descent_tsv(&report))
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

    /// Record (upsert) a source's haplogroup call for donor-level reconciliation. Defaults to the
    /// internal `NavigatorWalk` provenance tier (the external fast path records via `record_call_fp`).
    pub async fn record_haplogroup_call(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        call: &RunHaplogroupCall,
    ) -> Result<(), AppError> {
        self.record_haplogroup_call_fp(biosample_guid, dna_type, source_key, call, CallProvenance::NavigatorWalk, None)
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
        provenance: CallProvenance,
        fingerprint: Option<&str>,
    ) -> Result<(), AppError> {
        haplogroup_call::upsert(
            self.store.pool(),
            biosample_guid,
            dna_type,
            source_key,
            call,
            provenance,
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
        self.record_call_fp(
            biosample_guid,
            dna_type,
            source_key,
            source_label,
            assignment,
            CallProvenance::NavigatorWalk,
            None,
        )
        .await
    }

    /// Like [`record_call`](Self::record_call) but stamps the input fingerprint and the provenance
    /// tier (the sidecar fast path passes [`CallProvenance::External`]; internal genotyping passes
    /// [`CallProvenance::NavigatorWalk`]).
    #[allow(clippy::too_many_arguments)]
    pub(crate) async fn record_call_fp(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        source_label: String,
        assignment: &HaploAssignment,
        provenance: CallProvenance,
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
            self.record_haplogroup_call_fp(biosample_guid, dna_type, source_key, &call, provenance, fingerprint)
                .await?;
        }
        Ok(())
    }

    /// The preferred external (sidecar fast-path) call for an alignment's DNA type — present only
    /// when the "prefer external caller" policy is on **and** such a call exists. When present,
    /// Navigator's internal caller must not re-walk the CRAM: returning this call instead is what
    /// protects an external GATK4/1240K placement from being diluted or overwritten (the
    /// PRJEB37976 ancient-DNA fix). See `docs/design/external-caller-precedence.md`.
    pub(crate) async fn preferred_external_call(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        alignment_id: i64,
    ) -> Result<Option<RunHaplogroupCall>, AppError> {
        if !prefer_external_calls() {
            return Ok(None);
        }
        let key = match dna_type {
            DnaType::Y => external_y_source_key(alignment_id),
            DnaType::Mt => external_mt_source_key(alignment_id),
        };
        Ok(haplogroup_call::get_one(self.store.pool(), biosample_guid, dna_type, &key).await?)
    }

    /// Whether an alignment already carries a preferred external call for a DNA type — the gate the
    /// UI worker uses to skip enqueuing the internal Y/mt genotyping in "Full Analysis".
    pub async fn has_preferred_external_call(&self, alignment_id: i64, dna_type: DnaType) -> Result<bool, AppError> {
        let Ok(bio) = self.biosample_of_alignment(alignment_id).await else {
            return Ok(false);
        };
        Ok(self.preferred_external_call(bio, dna_type, alignment_id).await?.is_some())
    }

    /// "Compare callers": the trusted external caller vs Navigator's internal caller for one
    /// alignment. **Forces** the internal walk regardless of the prefer-external policy — it records
    /// its own `aln:{id}` / `aln:{id}:mt` (`NavigatorWalk`) rows and never touches the external `:ext`
    /// row, so the comparison is non-destructive to the external call. Returns Y (for Y-bearing
    /// subjects) and mtDNA, each with both terminals; a divergence is the ancient-DNA-damage signal
    /// the "skip the internal walk" default is protecting against. See external-caller-precedence §6.
    pub async fn compare_callers(&self, alignment_id: i64) -> Result<Vec<CallerComparison>, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();
        let mut out = Vec::new();

        let y_bearing = match bio {
            Some(g) => self.subject_has_y_dna(g).await.unwrap_or(true),
            None => true,
        };
        if y_bearing {
            let external = match bio {
                Some(g) => haplogroup_call::get_one(self.store.pool(), g, DnaType::Y, &external_y_source_key(alignment_id))
                    .await?
                    .map(|c| c.haplogroup),
                None => None,
            };
            let navigator = self
                .assign_y_haplogroup_walk(alignment_id, bio)
                .await
                .ok()
                .and_then(|a| a.ranked.first().map(|r| r.name.clone()));
            out.push(CallerComparison {
                dna_type: DnaType::Y,
                external,
                navigator,
            });
        }

        let external_mt = match bio {
            Some(g) => haplogroup_call::get_one(self.store.pool(), g, DnaType::Mt, &external_mt_source_key(alignment_id))
                .await?
                .map(|c| c.haplogroup),
            None => None,
        };
        let navigator_mt = self
            .assign_mtdna_haplogroup_walk(alignment_id, bio)
            .await
            .ok()
            .and_then(|a| a.ranked.first().map(|r| r.name.clone()));
        out.push(CallerComparison {
            dna_type: DnaType::Mt,
            external: external_mt,
            navigator: navigator_mt,
        });

        Ok(out)
    }

    /// The reconciled donor-level haplogroup consensus across all recorded sources. A user
    /// manual override, when set, replaces the computed terminal (flagged `overridden`).
    pub async fn haplogroup_consensus(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<Consensus>, AppError> {
        let calls = haplogroup_call::list_for_with_provenance(self.store.pool(), biosample_guid, dna_type).await?;
        let prefer_external = prefer_external_calls();
        let has_external = calls.iter().any(|(p, _)| *p == CallProvenance::External);
        // Per-run label reconciliation supplies the lineage / compatibility / divergence warnings —
        // honoring provenance: when the user prefers the external caller and one placed this subject,
        // it wins the vote (a damaged ancient-DNA CRAM walk cannot out-score it).
        let mut consensus = reconciliation::reconcile_with_provenance(&calls, prefer_external);

        // …the genome-level PLACED call (consensus_profile.consensus_label, from build_{y,mt}_profile)
        // is normally authoritative. Phase 2 makes that placement GVCF-sourced on preferred-external
        // subjects (place_{y,mt}_consensus → consensus_base_calls, no CRAM walk), so a freshly built
        // label already agrees with the external call. We still skip it here so a *stale* label left
        // by a pre-Phase-2 (CRAM-pooled) build cannot resurface before the profile is rebuilt — the
        // external reconcile is the safe authority for these subjects.
        let use_placed_label = !(prefer_external && has_external);
        if use_placed_label && matches!(dna_type, DnaType::Y | DnaType::Mt) {
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
        let mut groups: HashMap<(SampleGuid, DnaType), Vec<(CallProvenance, RunHaplogroupCall)>> = HashMap::new();
        for (guid, dna_type, prov, call) in haplogroup_call::list_all(self.store.pool()).await? {
            groups.entry((guid, dna_type)).or_default().push((prov, call));
        }
        let prefer_external = prefer_external_calls();
        let mut has_external: std::collections::HashSet<(SampleGuid, DnaType)> = std::collections::HashSet::new();
        let mut out: HashMap<SampleGuid, (Option<String>, Option<String>)> = HashMap::new();
        for ((guid, dna_type), calls) in groups {
            if calls.iter().any(|(p, _)| *p == CallProvenance::External) {
                has_external.insert((guid, dna_type));
            }
            if let Some(c) = reconciliation::reconcile_with_provenance(&calls, prefer_external) {
                let entry = out.entry(guid).or_default();
                match dna_type {
                    DnaType::Y => entry.0 = Some(c.haplogroup),
                    DnaType::Mt => entry.1 = Some(c.haplogroup),
                }
            }
        }
        // The genome-level placed terminal (build_{y,mt}_profile) wins over the per-run label vote,
        // so the subjects table matches the detail tab — except on a preferred-external subject, where
        // the CRAM-pooled placement is skipped in favor of the external call (as in haplogroup_consensus).
        for (guid_s, dna_type_s, label) in navigator_store::consensus_profile::list_labels(self.store.pool()).await? {
            let Ok(uuid) = guid_s.parse::<uuid::Uuid>() else {
                continue;
            };
            let guid = SampleGuid(uuid);
            let dna = match dna_type_s.as_str() {
                "Y" => DnaType::Y,
                "Mt" => DnaType::Mt,
                _ => continue,
            };
            if prefer_external && has_external.contains(&(guid, dna)) {
                continue;
            }
            let entry = out.entry(guid).or_default();
            match dna {
                DnaType::Y => entry.0 = Some(label),
                DnaType::Mt => entry.1 = Some(label),
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

    /// Load the persisted **observations** for a subject + DNA type — the raw genotype snapshot with
    /// no interpretation. Cheap (no genotyping). The shared loader behind [`cached_y_profile`] /
    /// [`cached_mt_profile`]; those interpret it against the current tree. `None` until a build runs.
    ///
    /// Backward-compat: a payload written before the observation-first switch is a baked
    /// [`ConsensusProfile`] (no `schema_version`); it is normalized to an [`ObservedProfile`] using
    /// its stored per-source bases (a source with no base becomes a no-call until the next rebuild).
    async fn load_observed_profile(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<navigator_domain::consensus::ObservedProfile>, AppError> {
        let Some(row) =
            navigator_store::consensus_profile::get(self.store.pool(), biosample_guid, dna_type.as_str()).await?
        else {
            return Ok(None);
        };
        let value: serde_json::Value = serde_json::from_str(&row.payload)?;
        // New payloads carry `schema_version`; legacy baked profiles don't.
        if value.get("schema_version").is_some() {
            Ok(Some(serde_json::from_value(value)?))
        } else {
            let legacy: ConsensusProfile = serde_json::from_value(value)?;
            Ok(Some(observed_from_legacy(legacy)))
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

    /// Persist a Y/mt **observation** snapshot (the payload) plus the interpreted `summary` header for
    /// quick listing. Only observations are stored; state/status are re-derived on load by
    /// [`interpret_y_profile`](Self::interpret_y_profile) / [`interpret_mt_profile`].
    async fn persist_observed_profile(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        observed: &navigator_domain::consensus::ObservedProfile,
        summary: &navigator_domain::consensus::ConsensusSummary,
        tree_provider: Option<String>,
    ) -> Result<(), AppError> {
        self.persist_consensus_row(
            biosample_guid,
            dna_type.as_str(),
            observed.terminal_hint.clone(),
            summary,
            observed.sources.len(),
            tree_provider,
            serde_json::to_string(observed)?,
        )
        .await
    }

    /// The Y/mt polarity map (SNP name → ancestral/derived) from the **current** tree for the
    /// configured provider — the input to [`navigator_domain::consensus::interpret`]. DecodingUs uses
    /// the tree's true phylogenetic polarity; FTDNA the parsed FTDNA tree's polarity. Empty when the
    /// tree is unavailable (interpret then falls back to each variant's stored ref/alt).
    async fn current_y_polarity(&self) -> std::collections::BTreeMap<String, (String, String)> {
        match y_tree_provider() {
            YTreeProvider::DecodingUs => self.decodingus_y_polarity().await.unwrap_or_default().into_iter().collect(),
            YTreeProvider::Ftdna => self
                .fetch_ftdna_y_tree()
                .await
                .ok()
                .and_then(|j| navigator_analysis::haplo::parse_ftdna_json(&j).ok())
                .map(|tree| navigator_analysis::haplo::polarity_from_tree(&tree))
                .unwrap_or_default(),
        }
    }

    /// The mtDNA polarity map from the current rCRS tree (DecodingUs remapped, FTDNA fallback).
    async fn current_mt_polarity(&self) -> std::collections::BTreeMap<String, (String, String)> {
        match self.mt_tree_rcrs().await {
            Ok((tree, _)) => navigator_analysis::haplo::polarity_from_tree(&tree),
            Err(_) => std::collections::BTreeMap::new(),
        }
    }

    /// Interpret stored Y observations against the current Y tree polarity into the display profile.
    async fn interpret_y_profile(&self, observed: navigator_domain::consensus::ObservedProfile) -> ConsensusProfile {
        let pol = self.current_y_polarity().await;
        interpret_observed(observed, &pol)
    }

    /// Interpret stored mtDNA observations against the current rCRS tree polarity.
    async fn interpret_mt_profile(&self, observed: navigator_domain::consensus::ObservedProfile) -> ConsensusProfile {
        let pol = self.current_mt_polarity().await;
        interpret_observed(observed, &pol)
    }

    /// The Y-profile for a subject, if one has been built — cheap (no genotyping). `None` until
    /// [`build_y_profile`](Self::build_y_profile) runs.
    ///
    /// Loads the stored **observations** and interprets them against the **current** Y tree polarity
    /// on every read — so a corrected/updated tree (or a provider switch) flips the derived/ancestral
    /// states with no rebuild and no BAM re-read. Legacy profiles are normalized to observations on
    /// load; those persisted before bases were stored show no-calls until one rebuild.
    pub async fn cached_y_profile(&self, biosample_guid: SampleGuid) -> Result<Option<YProfile>, AppError> {
        match self.load_observed_profile(biosample_guid, DnaType::Y).await? {
            Some(observed) => Ok(Some(self.interpret_y_profile(observed).await)),
            None => Ok(None),
        }
    }

    /// Build (and persist) the multi-source Y-variant profile: reconcile each Y-bearing source's
    /// per-SNP calls — every alignment's haplogroup placement, the combined chip/BISDNA placement,
    /// and the private-Y bucket — into one concordance view (confirmed / novel / conflict /
    /// single-source per SNP, with per-source provenance + per-observation quality weighting).
    /// Expensive (re-genotypes each alignment), so it's an explicit action; the result is persisted
    /// so [`cached_y_profile`](Self::cached_y_profile) reloads it instantly. Sources without Y data
    /// are skipped.
    pub async fn build_y_profile(&self, biosample_guid: SampleGuid) -> Result<YProfile, AppError> {
        // Females have no Y chromosome — don't build or persist a Y variant profile for them.
        if !self.subject_has_y_dna(biosample_guid).await? {
            return Ok(YProfile {
                variants: Vec::new(),
                summary: yprofile::summarize(&[]),
                terminal: None,
                sources: Vec::new(),
            });
        }

        let mut sources: Vec<(String, SourceType, Vec<YObsInput>)> = Vec::new();

        // WGS / Y-NGS evidence: the **genome-consensus** deep placement — every alignment's chrY
        // calls pooled on ONE tree+coordinate space and placed once ([`place_y_consensus`]). Its
        // root→terminal lineage is the SNP set the sample carries all the way down to the deep
        // terminal, so the descent report renders a populated backbone.
        //
        // Previously this looped each alignment through `y_assignment_full`, whose *per-alignment*
        // placement is shallow on lifted CHM13 Big Y data (it stops a few clades down): the profile
        // then carried only the root→shallow-terminal SNPs while the profile's terminal came from
        // the deeper pooled consensus — so the descent walked terminal→root over SNPs the profile
        // never recorded, rendering every node below the shallow terminal as no-call (the "all
        // no-call below F / reversed SNPs" bug). Pooling first keeps the variants and the terminal
        // on the same deep placement. Genotypes are cached, so this reuses the Y walk already paid.
        let consensus_assignment = self.place_y_consensus(biosample_guid).await?;
        if let Some(asg) = &consensus_assignment {
            let obs = snp_obs_from_assignment(asg, true);
            if !obs.is_empty() {
                sources.push(("genome consensus".to_string(), SourceType::WgsShortRead, obs));
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
                    // Carry the observed base (= the called alt) so interpret re-derives Derived from
                    // the call's own ref/alt — no baked state, consistent with every other source.
                    let mut o = YObsInput::observed(
                        name,
                        v.position,
                        v.reference.to_string(),
                        v.alternate.to_string(),
                        Some(v.alternate),
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

        // Group the sources into an observation-only snapshot — no baked state. State/status are
        // interpreted against the current tree on read (`interpret_y_profile`), so a corrected tree
        // (incl. an FTDNA tree whose reference-as-ancestral polarity is inverted at some sites) flips
        // the display without a rebuild.
        let mut observed = yprofile::to_observed(&sources);
        // Genome-level placement: the pooled call set placed once (computed above — not a vote among
        // the per-run terminal labels). Falls back to the label reconciliation when nothing places.
        observed.terminal_hint = match &consensus_assignment {
            Some(a) => a.ranked.first().map(|r| r.name.clone()),
            None => self
                .haplogroup_consensus(biosample_guid, DnaType::Y)
                .await?
                .map(|c| c.haplogroup),
        };

        // Interpret once for the return value + the persisted summary header.
        let profile = self.interpret_y_profile(observed.clone()).await;
        let provider = Some(match y_tree_provider() {
            YTreeProvider::DecodingUs => "decodingus".to_string(),
            YTreeProvider::Ftdna => "ftdna".to_string(),
        });
        self.persist_observed_profile(biosample_guid, DnaType::Y, &observed, &profile.summary, provider)
            .await?;
        Ok(profile)
    }

    /// The mtDNA consensus profile for a subject, if built — cheap (no genotyping). Loads stored
    /// observations and interprets them against the current rCRS tree polarity on every read (the
    /// mtDNA half of the observation-first fix — previously mt states could not re-interpret at all).
    pub async fn cached_mt_profile(&self, biosample_guid: SampleGuid) -> Result<Option<ConsensusProfile>, AppError> {
        match self.load_observed_profile(biosample_guid, DnaType::Mt).await? {
            Some(observed) => Ok(Some(self.interpret_mt_profile(observed).await)),
            None => Ok(None),
        }
    }

    /// Build (and persist) the multi-source mtDNA consensus profile — the mtDNA adapter over the
    /// generic [`navigator_domain::consensus`] engine. Reconciles each mt-bearing source's
    /// defining-mutation calls (every alignment's chrM placement, each imported mtDNA FASTA
    /// sequence's placement, and the combined chip mtDNA placement) into one concordance view,
    /// keyed by phylotree **mutation name** (rCRS-coordinate, build-independent). Persisted with
    /// `dna_type='Mt'` so [`cached_mt_profile`](Self::cached_mt_profile) reloads it instantly.
    /// Expensive (re-places each alignment's chrM), so it's an explicit action; mt-less sources skip.
    pub async fn build_mt_profile(&self, biosample_guid: SampleGuid) -> Result<ConsensusProfile, AppError> {
        // One mt tree in rCRS coordinates (DecodingUs remapped from hs1, FTDNA fallback), shared by
        // the per-source placements below and the pooled terminal — so the variants and the terminal
        // sit on the same tree + coordinate space (the Y-profile fix, applied to mtDNA).
        let (tree, provider) = self.mt_tree_rcrs().await?;
        let source_calls = self.mt_source_calls(biosample_guid, &tree).await?;

        // One source per contributing test — each alignment's chrM, each imported FASTA, the chip mt
        // panel — placed individually on the shared tree so the profile shows which test confirmed
        // each mutation (name-keyed reconcile across sources). Sparse sources (chip) use the robust
        // assembler; dense ones (WGS/FASTA) the exact one.
        let mut sources: Vec<(String, SourceType, Vec<YObsInput>)> = Vec::new();
        for (label, st, calls) in &source_calls {
            let assignment = if *st == SourceType::Chip {
                assemble_assignment_robust(&tree, calls)
            } else {
                assemble_assignment(&tree, calls)
            };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push((label.clone(), *st, obs));
            }
        }

        // Observation-only snapshot; state/status interpreted against the rCRS tree polarity on read.
        let mut observed = yprofile::to_observed(&sources);
        // Interpret once (against the current mt polarity) for the return value + summary header.
        let mut profile = self.interpret_mt_profile(observed.clone()).await;
        // Genome-level placement of the pooled chrM call set on the same tree. A subject with no
        // derived mutations carries no real placement — an alignment with a handful of off-target
        // chrM reads (a Big Y) genotypes to nothing below the root. Report no terminal rather than a
        // root label, and never resurrect a stale persisted root call (the old "very few mt reads →
        // RSRS" artifact); the profile is meaningful only once some mutation is derived (checked on
        // the *interpreted* variants).
        let terminal = if profile
            .variants
            .iter()
            .all(|v| v.consensus != navigator_domain::consensus::ConsensusState::Derived)
        {
            None
        } else {
            let has_wgs = source_calls.iter().any(|(_, st, _)| *st == SourceType::WgsShortRead);
            let pooled_input: Vec<(SourceType, HashMap<i64, char>)> =
                source_calls.iter().map(|(_, st, calls)| (*st, calls.clone())).collect();
            let pooled = pool_votes(&pooled_input);
            let assignment = if has_wgs {
                assemble_assignment(&tree, &pooled)
            } else {
                assemble_assignment_robust(&tree, &pooled)
            };
            assignment.ranked.first().map(|r| r.name.clone())
        };
        observed.terminal_hint = terminal.clone();
        profile.terminal = terminal;

        // Persist observations (keyed dna_type='Mt') with the tree provider actually used.
        self.persist_observed_profile(biosample_guid, DnaType::Mt, &observed, &profile.summary, Some(provider.to_string()))
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
        // Females have no Y chromosome — no genome consensus to place.
        if !self.subject_has_y_dna(biosample_guid).await? {
            return Ok(None);
        }
        // The consensus follows the user's configured tree provider (Preferences /
        // NAVIGATOR_Y_TREE_PROVIDER), same as the per-alignment placement.
        match y_tree_provider() {
            YTreeProvider::DecodingUs => self.place_y_consensus_decodingus(biosample_guid).await,
            YTreeProvider::Ftdna => self.place_y_consensus_ftdna(biosample_guid).await,
        }
    }

    /// FTDNA-provider genome consensus: pool every WGS alignment + GRCh38 vendor Y-VCF on the FTDNA
    /// GRCh38 tree (`base_calls` lifts CHM13/GRCh37 sources into GRCh38) and place once. One tree +
    /// one coordinate space keeps polarity/coverage consistent; chips (sparse, various builds) stay in
    /// the variant profile's name-keyed reconcile.
    async fn place_y_consensus_ftdna(&self, biosample_guid: SampleGuid) -> Result<Option<HaploAssignment>, AppError> {
        let tree_json = self.fetch_ftdna_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;

        let mut sources: Vec<(SourceType, HashMap<i64, char>)> = Vec::new();
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            // Lifted GRCh38-coordinate calls; sources lacking chrY / a reference are skipped. A
            // preferred-external alignment is genotyped from its chrY GVCF (lifted native→GRCh38 by
            // `gvcf_base_calls`) instead of walking the CRAM.
            let calls = match (prefer_external_calls(), chr_y_gvcf_for_alignment(a)) {
                (true, Some(gvcf)) => self
                    .gvcf_base_calls(a.id, "chrY", &gvcf, &tree, tree_build_for_contig("chrY"))
                    .await
                    .ok(),
                _ => self.assign_haplogroup_detail(a.id, "chrY", &tree_json).await.ok().map(|(_, _, c)| c),
            };
            let Some(calls) = calls else { continue };
            if !calls.is_empty() {
                sources.push((SourceType::WgsShortRead, calls));
            }
        }
        // Dense GRCh38 vendor Y-NGS VCFs pool alongside the WGS; non-GRCh38 sets wouldn't match.
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

    /// DecodingUs-provider genome consensus (the default): genotype every WGS alignment against the
    /// DecodingUs Y tree in each source's *native* build, group by build, pool by position, and place
    /// on the build carrying the most evidence.
    async fn place_y_consensus_decodingus(&self, biosample_guid: SampleGuid) -> Result<Option<HaploAssignment>, AppError> {
        // Genotype every WGS alignment against the **DecodingUs** Y tree — the workspace's configured
        // provider, served from the local cache — in each source's *native* build (`hs1` for CHM13,
        // `GRCh38`, `GRCh37`). No liftover and no FTDNA dependency: the per-alignment genotype is
        // exactly the one the Y assignment already cached, so this reuses that walk rather than paying
        // a second, FTDNA-coordinate one. Sources are grouped by build and pooled by **position**
        // within one coordinate space (a single build per subject is the norm); the build carrying the
        // most evidence is placed once. Pooling across builds by position would mix coordinate systems
        // — cross-build merging lives in the variant profile's name-keyed reconcile, not here.
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let vsets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;

        // Parse the DecodingUs tree once per distinct build the sources use (cheap — the JSON is
        // memoized). Built up front so the async genotyping loop holds only shared borrows of `trees`.
        let mut builds: std::collections::HashSet<&'static str> =
            alignments.iter().filter_map(|a| decodingus_build_key(&a.reference_build)).collect();
        for set in &vsets {
            if set.source_type != SourceType::Chip {
                // Unknown vendor build → GRCh38 (the vendor-Y-VCF import default).
                if let Some(bk) = set.reference_build.as_deref().map_or(Some("GRCh38"), decodingus_build_key) {
                    builds.insert(bk);
                }
            }
        }
        let mut trees: HashMap<&'static str, navigator_analysis::haplo::HaploTree> = HashMap::new();
        for bk in builds {
            if let Ok(t) = navigator_analysis::haplo::parse_decodingus_json(&tree_json, bk) {
                trees.insert(bk, t);
            }
        }

        let mut by_build: HashMap<&'static str, YSourceCalls> = HashMap::new();
        for a in &alignments {
            let Some(bk) = decodingus_build_key(&a.reference_build) else { continue };
            let Some(tree) = trees.get(bk) else { continue };
            // Native build → no liftover; the cache-key matches the Y assignment's, so a CRAM walk is
            // a hit — but a preferred-external alignment is genotyped from its GVCF instead (no decode).
            let Ok(calls) = self.consensus_base_calls(a, "chrY", tree, None).await else { continue };
            if !calls.is_empty() {
                by_build.entry(bk).or_default().push((SourceType::WgsShortRead, calls));
            }
        }

        // Vendor Y-NGS VCFs (FTDNA Big Y / YSEQ / Full Genomes / Nebula) are dense direct Y-SNP calls;
        // fold each into its own build's group (strand-reconciled to that build's tree). Chips stay in
        // the variant profile's name-keyed reconcile.
        for set in &vsets {
            if set.source_type == SourceType::Chip {
                continue;
            }
            let Some(bk) = set.reference_build.as_deref().map_or(Some("GRCh38"), decodingus_build_key) else { continue };
            let Some(tree) = trees.get(bk) else { continue };
            let calls = Self::vset_chr_y_calls(set);
            if !calls.is_empty() {
                by_build.entry(bk).or_default().push((set.source_type, strand_reconcile_to_tree(tree, calls)));
            }
        }

        // Place on the build carrying the most evidence (the subject's primary coordinate space).
        let Some(bk) = by_build
            .iter()
            .max_by_key(|(_, s)| s.iter().map(|(_, c)| c.len()).sum::<usize>())
            .map(|(bk, _)| *bk)
        else {
            return Ok(None);
        };
        let pooled = pool_votes(&by_build[bk]);
        Ok(Some(assemble_assignment(&trees[bk], &pooled)))
    }

    /// Diagnostic: dump the Y **descent** for one subject SNP-by-SNP — the reported state + observed
    /// base against the **incoming tree's** polarity in every DecodingUs build (hs1 / GRCh38 / GRCh37).
    /// This is the "compare the tree vs the calls" log: a backbone SNP the sample must carry that reads
    /// ancestral shows here as `state=Ancestral base=<der allele>` with a build whose polarity is
    /// flipped (`hs1: A>G  GRCh38: G>A`), pinpointing a tree-polarity problem vs a genotyping one.
    /// Read-only. TSV: `node  snp  pos  state  base  hs1  GRCh38  GRCh37`.
    pub async fn debug_y_descent(&self, biosample_guid: SampleGuid) -> Result<String, AppError> {
        use navigator_analysis::haplo;
        let Some(report) = self.descent_report(biosample_guid, DnaType::Y).await? else {
            return Ok("no Y descent (build the variant profile first)".into());
        };
        let json = self.fetch_decodingus_y_tree().await?;
        let pol = |bk: &str| -> std::collections::BTreeMap<String, (String, String)> {
            haplo::parse_decodingus_json(&json, bk)
                .ok()
                .map(|t| haplo::polarity_from_tree(&t))
                .unwrap_or_default()
        };
        let (hs1, g38, g37) = (pol("hs1"), pol("GRCh38"), pol("GRCh37"));
        let show = |m: &std::collections::BTreeMap<String, (String, String)>, name: &str| -> String {
            m.get(&name.trim().to_uppercase())
                .map(|(a, d)| format!("{a}>{d}"))
                .unwrap_or_else(|| "-".into())
        };
        let mut out = format!(
            "descent terminal: {}  nodes-on-path: {}\nnode\tsnp\tpos\tstate\tbase\ths1\tGRCh38\tGRCh37\n",
            report.terminal,
            report.nodes.len()
        );
        let (mut derived, mut ancestral, mut nocall) = (0usize, 0usize, 0usize);
        for node in &report.nodes {
            for snp in &node.snps {
                match snp.state {
                    navigator_analysis::haplo::CallState::Derived => derived += 1,
                    navigator_analysis::haplo::CallState::Ancestral => ancestral += 1,
                    navigator_analysis::haplo::CallState::NoCall => nocall += 1,
                }
                out.push_str(&format!(
                    "{}\t{}\t{}\t{:?}\t{}\t{}\t{}\t{}\n",
                    node.name,
                    snp.name,
                    snp.position,
                    snp.state,
                    snp.base.map(|c| c.to_string()).unwrap_or_else(|| "-".into()),
                    show(&hs1, &snp.name),
                    show(&g38, &snp.name),
                    show(&g37, &snp.name),
                ));
            }
        }
        out.push_str(&format!("\ntotals: derived={derived} ancestral={ancestral} nocall={nocall}\n"));
        Ok(out)
    }

    /// Diagnostic: genotype a **single alignment** against the DecodingUs Y tree in its native build
    /// and dump, per SNP down the placed lineage, the raw read pileup **behind** each call — the
    /// reference base, the A/C/G/T passing-read tally, the consensus base, the tree's ancestral/derived
    /// alleles, and the resulting state. This is the "calls generated" log: it shows whether a backbone
    /// SNP that reads ancestral is (a) genuinely ancestral in the reads, (b) a coordinate/position
    /// mismatch (reads a different base than the tree allele), or (c) a low-depth artifact. Read-only.
    /// TSV: `node  snp  pos  tree(anc>der)  ref  A  C  G  T  depth  called  state`.
    pub async fn debug_y_calls(&self, alignment_id: i64) -> Result<String, AppError> {
        use navigator_analysis::{caller, haplo, reader};
        let aln = self.alignment_or_err(alignment_id).await?;
        let Some(bk) = decodingus_build_key(&aln.reference_build) else {
            return Ok(format!(
                "alignment {alignment_id}: build {} has no DecodingUs tree key",
                aln.reference_build
            ));
        };
        let json = self.fetch_decodingus_y_tree().await?;
        let tree = haplo::parse_decodingus_json(&json, bk).map_err(AppError::Import)?;
        // Native-build genotyping (no liftover) — the same walk place_y_consensus uses; the cache hit
        // means `base_calls` returns the identical winning bases we're auditing here.
        let calls = self.base_calls(alignment_id, "chrY", &tree, None).await?;
        let assignment = assemble_assignment(&tree, &calls);
        if assignment.lineage.is_empty() {
            return Ok(format!("alignment {alignment_id} ({bk}): no Y placement\n"));
        }

        // Localize the BAM/CRAM and resolve its reference exactly as `base_calls` does, then tally the
        // raw reads at the lineage positions and read the reference base there.
        let bam = self
            .localize(Path::new(&aln.bam_path.clone().ok_or(AppError::MissingPaths(alignment_id))?))
            .await;
        let is_cram = bam.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        let reference = match aln.reference_path.clone() {
            Some(p) => Some(PathBuf::from(p)),
            None if is_cram => Some(self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?),
            None => self.gateway.cached_reference(&aln.reference_build),
        };
        let targets: HashSet<i64> = assignment.lineage.iter().map(|e| e.position).collect();
        let resolved = self
            .resolve_header_contig(&bam, reference.as_deref(), "chrY")
            .await?
            .unwrap_or_else(|| "chrY".to_string());
        let (bam2, ref2, contig2, targets2) = (bam.clone(), reference.clone(), resolved.clone(), targets.clone());
        let counts = tokio::task::spawn_blocking(move || {
            let params = adaptive_haploid_params(&bam2, ref2.as_deref());
            caller::tally_at(&bam2, &contig2, &targets2, &params, ref2.as_deref())
        })
        .await??;
        // Reference base per lineage position (0-based index into the contig), best-effort.
        let refseq: Option<Vec<u8>> = match reference.as_deref() {
            Some(r) => {
                let (r, c) = (r.to_path_buf(), resolved.clone());
                tokio::task::spawn_blocking(move || reader::read_contig_sequence(&r, &c).ok()).await?
            }
            None => None,
        };
        let ref_at = |pos: i64| -> char {
            refseq
                .as_ref()
                .and_then(|s| s.get((pos - 1) as usize))
                .map(|b| (*b as char).to_ascii_uppercase())
                .unwrap_or('?')
        };

        let mut out = format!(
            "alignment {alignment_id}  build {bk}  contig {resolved}  terminal {}\nsnp\tpos\ttree\tref\tA\tC\tG\tT\tdepth\tcalled\tstate\n",
            assignment.ranked.first().map(|r| r.name.as_str()).unwrap_or("?"),
        );
        let (mut derived, mut ancestral, mut nocall) = (0usize, 0usize, 0usize);
        for e in &assignment.lineage {
            match e.state {
                navigator_analysis::haplo::CallState::Derived => derived += 1,
                navigator_analysis::haplo::CallState::Ancestral => ancestral += 1,
                navigator_analysis::haplo::CallState::NoCall => nocall += 1,
            }
            let c = counts.get(&e.position).copied().unwrap_or([0; 4]);
            let depth: u32 = c.iter().sum();
            out.push_str(&format!(
                "{}\t{}\t{}>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\n",
                e.name,
                e.position,
                e.ancestral,
                e.derived,
                ref_at(e.position),
                c[0],
                c[1],
                c[2],
                c[3],
                depth,
                e.base.map(|b| b.to_string()).unwrap_or_else(|| "-".into()),
                e.state,
            ));
        }
        out.push_str(&format!("\ntotals: derived={derived} ancestral={ancestral} nocall={nocall}\n"));
        Ok(out)
    }

    /// Pick a single alignment to target for [`Self::debug_y_calls`] when only a subject is given:
    /// prefer a CHM13/HiFi alignment (native tree, no liftover — the cleanest to audit), else the
    /// first. Returns `None` when the subject has no alignments.
    pub async fn pick_y_debug_alignment(&self, biosample_guid: SampleGuid) -> Result<Option<i64>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let pick = alignments
            .iter()
            .find(|a| {
                decodingus_build_key(&a.reference_build) == Some("hs1")
                    && a.aligner.to_ascii_lowercase().contains("pbmm2")
            })
            .or_else(|| alignments.iter().find(|a| decodingus_build_key(&a.reference_build) == Some("hs1")))
            .or_else(|| alignments.first());
        Ok(pick.map(|a| a.id))
    }

    /// Pick a single alignment to genotype **mtDNA** against: skip Y-only runs (an FTDNA Big-Y
    /// carries no `chrM` reads, so it would yield an all-no-call report while a usable WGS
    /// alignment sat unselected), then prefer CHM13, else the first survivor. If *every* run is
    /// Y-only, fall back to the first alignment rather than reporting "no alignment".
    pub async fn pick_mt_alignment(&self, biosample_guid: SampleGuid) -> Result<Option<i64>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut mt_capable: Vec<&Alignment> = Vec::new();
        for a in &alignments {
            let y_only = match sequence_run::get(self.store.pool(), a.sequence_run_id).await? {
                // `target_of` is tolerant of unknown codes (→ None), which stay eligible.
                Some(run) => navigator_domain::testtype::target_of(&run.test_type)
                    == Some(navigator_domain::testtype::TargetType::YChromosome),
                None => false,
            };
            if !y_only {
                mt_capable.push(a);
            }
        }
        let pick = mt_capable
            .iter()
            .copied()
            .find(|a| decodingus_build_key(&a.reference_build) == Some("hs1"))
            .or_else(|| mt_capable.first().copied())
            .or_else(|| alignments.first());
        Ok(pick.map(|a| a.id))
    }

    /// The alignment to genotype `dna` against for a subject-keyed query.
    pub async fn pick_alignment_for(&self, guid: SampleGuid, dna: DnaType) -> Result<Option<i64>, AppError> {
        match dna {
            DnaType::Y => self.pick_y_debug_alignment(guid).await,
            DnaType::Mt => self.pick_mt_alignment(guid).await,
        }
    }

    /// Diagnostic: trace the DecodingUs genome-consensus Y placement for one subject — the pooled
    /// build, the Kulczynski top candidates + admissibility, the assembled terminal, a strict
    /// root→tip descent, and the per-node derived/ancestral tally down the assembled lineage (so an
    /// over-deepening tunnel through ancestral branches is visible). Read-only.
    pub async fn debug_y_placement(&self, biosample_guid: SampleGuid) -> Result<String, AppError> {
        use navigator_analysis::haplo;
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let vsets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut builds: std::collections::HashSet<&'static str> =
            alignments.iter().filter_map(|a| decodingus_build_key(&a.reference_build)).collect();
        for set in &vsets {
            if set.source_type != SourceType::Chip {
                if let Some(bk) = set.reference_build.as_deref().map_or(Some("GRCh38"), decodingus_build_key) {
                    builds.insert(bk);
                }
            }
        }
        let mut trees: HashMap<&'static str, haplo::HaploTree> = HashMap::new();
        for bk in builds {
            if let Ok(t) = haplo::parse_decodingus_json(&tree_json, bk) {
                trees.insert(bk, t);
            }
        }
        let mut by_build: HashMap<&'static str, YSourceCalls> = HashMap::new();
        for a in &alignments {
            let Some(bk) = decodingus_build_key(&a.reference_build) else { continue };
            let Some(tree) = trees.get(bk) else { continue };
            if let Ok(calls) = self.base_calls(a.id, "chrY", tree, None).await {
                if !calls.is_empty() {
                    by_build.entry(bk).or_default().push((SourceType::WgsShortRead, calls));
                }
            }
        }
        for set in &vsets {
            if set.source_type == SourceType::Chip {
                continue;
            }
            let Some(bk) = set.reference_build.as_deref().map_or(Some("GRCh38"), decodingus_build_key) else {
                continue;
            };
            let Some(tree) = trees.get(bk) else { continue };
            let calls = Self::vset_chr_y_calls(set);
            if !calls.is_empty() {
                by_build.entry(bk).or_default().push((set.source_type, strand_reconcile_to_tree(tree, calls)));
            }
        }
        let Some(bk) = by_build
            .iter()
            .max_by_key(|(_, s)| s.iter().map(|(_, c)| c.len()).sum::<usize>())
            .map(|(bk, _)| *bk)
        else {
            return Ok("no Y calls".into());
        };
        let pooled = pool_votes(&by_build[bk]);
        let tree = &trees[bk];

        let mut out = format!("build={bk} pooled_calls={}\n", pooled.len());
        let ranked = haplo::score(tree, &pooled);
        out.push_str("Kulczynski top-10 (name score matched/expected admissible):\n");
        for r in ranked.iter().take(10) {
            out.push_str(&format!(
                "  {:<28} {:.3} {}/{} adm={}\n",
                r.name,
                r.score,
                r.matched,
                r.expected,
                haplo::path_admissible(tree, &pooled, r.id)
            ));
        }
        let asg = assemble_assignment(tree, &pooled);
        let asm_id = asg.ranked.first().map(|r| r.id);
        out.push_str(&format!(
            "assemble_assignment terminal: {}\n",
            asg.ranked.first().map(|r| r.name.as_str()).unwrap_or("?")
        ));
        let mut roots: Vec<i64> = tree.nodes.values().filter(|n| n.is_root).map(|n| n.id).collect();
        roots.sort_unstable();
        for &r in &roots {
            let term = haplo::deepen_terminal(tree, &pooled, r);
            out.push_str(&format!(
                "deepen_terminal(root={}) -> {}\n",
                tree.nodes.get(&r).map(|n| n.name.as_str()).unwrap_or("?"),
                tree.nodes.get(&term).map(|n| n.name.as_str()).unwrap_or("?")
            ));
        }
        if let Some(tid) = asm_id {
            out.push_str("assembled lineage (root->terminal) per-node (derived/ancestral/nocall):\n");
            for id in lineage_ids(tree, tid) {
                let (d, a, n) = haplo::node_call_counts(tree, &pooled, id);
                let name = tree.nodes.get(&id).map(|x| x.name.clone()).unwrap_or_default();
                out.push_str(&format!("  {:<32} d={d} a={a} n={n}\n", name));
            }
        }
        Ok(out)
    }

    /// The DecodingUs coordinate key (`hs1` / `GRCh38` / `GRCh37`) for a subject's Y data, taken from
    /// its first alignment's reference build. Defaults to `hs1` (CHM13, the DecodingUs native build)
    /// when the subject has no build-resolvable alignment. Used to parse the DecodingUs tree in the
    /// subject's own coordinate space for the descent report.
    async fn subject_y_build_key(&self, biosample_guid: SampleGuid) -> &'static str {
        alignment::list_for_biosample(self.store.pool(), biosample_guid)
            .await
            .ok()
            .and_then(|alns| alns.iter().find_map(|a| decodingus_build_key(&a.reference_build)))
            .unwrap_or("hs1")
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

    /// Per-source rCRS-coordinate mtDNA calls for a subject — `(label, type, calls)` keyed by rCRS
    /// position — shared by [`place_mt_consensus`] (pooled placement) and [`build_mt_profile`]
    /// (per-source concordance). `tree` must be in rCRS coordinates: each alignment's `chrM` is
    /// genotyped against it (cached; `base_calls` maps a CHM13 `chrM` back to rCRS, GRCh38/rCRS
    /// direct), each imported FASTA is sampled at every rCRS position, and the chip mt panel is
    /// strand-reconciled to it.
    async fn mt_source_calls(
        &self,
        biosample_guid: SampleGuid,
        tree: &navigator_analysis::haplo::HaploTree,
    ) -> Result<Vec<(String, SourceType, HashMap<i64, char>)>, AppError> {
        let mut sources: Vec<(String, SourceType, HashMap<i64, char>)> = Vec::new();

        // Each alignment's chrM genotype. `None` source-build → rCRS-direct / CHM13-chrM lift. A
        // preferred-external alignment is genotyped from its chrM GVCF instead of the CRAM.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            let Ok(calls) = self.consensus_base_calls(a, "chrM", tree, None).await else {
                continue;
            };
            if !calls.is_empty() {
                sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, calls));
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
                let vendor = mt_vendor_label(s.source_file_name.as_deref(), s.defline.as_deref());
                sources.push((format!("{vendor} (mt seq #{})", s.id), SourceType::Imported, calls));
            }
        }

        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;

        // Each imported **non-chip** variant set carrying chrM SNPs — a whole-genome VCF or a
        // CompleteGenomics masterVar. These report forward-strand ref/alt on rCRS coordinates
        // (GRCh37/GRCh38 chrM = rCRS), so the alt base is used raw like an alignment's chrM call
        // (no TOP-strand reconciliation, unlike the chip panel below).
        for set in sets.iter().filter(|s| s.source_type != SourceType::Chip) {
            let calls: HashMap<i64, char> = set
                .calls
                .iter()
                .filter(|c| {
                    c.contig.eq_ignore_ascii_case("chrM")
                        || c.contig.eq_ignore_ascii_case("chrMT")
                        || c.contig.eq_ignore_ascii_case("mt")
                        || c.contig.eq_ignore_ascii_case("m")
                })
                .filter_map(|c| c.alternate.chars().next().map(|b| (c.position, b.to_ascii_uppercase())))
                .collect();
            if !calls.is_empty() {
                sources.push((set.source_label.clone(), set.source_type, calls));
            }
        }

        // The chip mt panel (consumer arrays carry a sparse rCRS MT panel), strand-reconciled.
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
            sources.push(("Chip mtDNA panel".to_string(), SourceType::Chip, strand_reconcile_to_tree(tree, chip_mt)));
        }

        Ok(sources)
    }

    /// **Genome-level mtDNA placement**: the mt counterpart to [`place_y_consensus`]. Pools every
    /// source's rCRS-coordinate genotype ([`mt_source_calls`]) by [`pool_votes`] vote keyed by
    /// **position** (rCRS is the only mt coordinate system → no name indirection), then places the
    /// pooled set on the mt tree once. The tree is the **DecodingUs** mt tree (the configured
    /// provider) remapped onto rCRS, with the FTDNA mt tree as fallback ([`mt_tree_rcrs`]). `Ok(None)`
    /// when the subject has no mt-bearing source.
    pub async fn place_mt_consensus(&self, biosample_guid: SampleGuid) -> Result<Option<HaploAssignment>, AppError> {
        let (tree, _provider) = self.mt_tree_rcrs().await?;
        let sources = self.mt_source_calls(biosample_guid, &tree).await?;
        if sources.is_empty() {
            return Ok(None);
        }
        let has_wgs = sources.iter().any(|(_, st, _)| *st == SourceType::WgsShortRead);
        // mt is single-coordinate (rCRS) across all sources, so a base vote is strand-safe here.
        let pooled_input: Vec<(SourceType, HashMap<i64, char>)> =
            sources.into_iter().map(|(_, st, calls)| (st, calls)).collect();
        let pooled = pool_votes(&pooled_input);
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

        // Render on the configured provider's tree so the node names + defining SNPs line up with the
        // profile's placement (which followed the same provider). Y: DecodingUs in the subject's
        // native build, or the FTDNA GRCh38 tree; mtDNA: DecodingUs remapped hs1→rCRS, or FTDNA — via
        // `mt_tree_rcrs`, which already honors the provider.
        let tree = match dna {
            DnaType::Y => match y_tree_provider() {
                YTreeProvider::DecodingUs => {
                    let json = self.fetch_decodingus_y_tree().await?;
                    let build_key = self.subject_y_build_key(biosample_guid).await;
                    navigator_analysis::haplo::parse_decodingus_json(&json, build_key).map_err(AppError::Import)?
                }
                YTreeProvider::Ftdna => {
                    let json = self.fetch_ftdna_y_tree().await?;
                    navigator_analysis::haplo::parse_ftdna_json(&json).map_err(AppError::Import)?
                }
            },
            DnaType::Mt => self.mt_tree_rcrs().await?.0,
        };
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
        // The actual consensus nucleotide per SNP, so the descent shows/exports the observed allele.
        let base_by_name: std::collections::HashMap<String, char> = profile
            .variants
            .iter()
            .filter_map(|v| {
                v.consensus_base
                    .as_deref()
                    .and_then(|b| b.chars().next())
                    .map(|c| (v.name.clone(), c))
            })
            .collect();

        let mut nodes = navigator_analysis::haplo::descent_by_node(&tree, terminal_id, &state_by_name);
        for node in &mut nodes {
            for snp in &mut node.snps {
                if snp.base.is_none() {
                    snp.base = base_by_name.get(&snp.name).copied();
                }
            }
        }
        Ok(Some(DescentReport { dna, terminal, nodes }))
    }

    /// Build a [`BranchReport`]: the sample's genotype at every defining marker of `node_query`'s
    /// descendant subtree (Y or mtDNA), with per-marker evidence — for spot-checking placement and
    /// exchanging observations. `node_query` matches a haplogroup name (`R-FGC29071`) or a defining
    /// marker (`FGC29071`); `max_depth` bounds descent (`None` = the whole subtree).
    ///
    /// Genotypes the subtree **fresh** over the tree's loci (not the placement profile), so branches
    /// the sample is *ancestral* for are reported too. Observed bases + evidence come from a per-
    /// sample chrY GVCF sidecar when present (rich DP/AD/GQ, ref blocks), else the pileup caller.
    pub async fn branch_report(
        &self,
        alignment_id: i64,
        dna: DnaType,
        node_query: &str,
        max_depth: Option<usize>,
    ) -> Result<BranchReport, AppError> {
        use std::collections::HashSet;

        let aln = self.alignment_or_err(alignment_id).await?;

        // Tree + observed base calls over ALL tree loci (covers off-path descendant branches).
        let (tree, calls, contig, gvcf) = match dna {
            DnaType::Y => {
                // Probe the sidecar *before* genotyping: with one present the tree comes straight
                // from JSON and the calls from the GVCF, so the per-locus pileup walk is skipped
                // entirely (the point of the fast path — cf. `assign_y_from_gvcf`).
                match crate::fastpath::chr_y_gvcf_for_alignment(&aln) {
                    Some(gvcf) => {
                        let build_key = decodingus_build_key(&aln.reference_build).ok_or_else(|| {
                            AppError::Import(format!(
                                "no DecodingUs tree coordinates for build {}",
                                aln.reference_build
                            ))
                        })?;
                        let tree_json = self.fetch_decodingus_y_tree().await?;
                        let tree = navigator_analysis::haplo::parse_decodingus_json(&tree_json, build_key)
                            .map_err(AppError::Import)?;
                        let calls = self.gvcf_base_calls(alignment_id, "chrY", &gvcf, &tree, None).await?;
                        (tree, calls, "chrY", Some(gvcf))
                    }
                    None => {
                        let (tree, calls) = self.y_decodingus_tree_calls(alignment_id).await?;
                        (tree, calls, "chrY", None)
                    }
                }
            }
            DnaType::Mt => {
                // `mt_tree_rcrs` hands back the *provider* (decodingus/ftdna), not a reference
                // build. `tree_source_build` must stay `None`: a non-build string there makes
                // `lifted_targets` return early, skipping the rCRS↔chrM map a CHM13 alignment
                // needs (its chrM is a circular permutation of rCRS).
                let (tree, _provider) = self.mt_tree_rcrs().await?;
                let calls = self.base_calls(alignment_id, "chrM", &tree, None).await?;
                (tree, calls, "chrM", None)
            }
        };
        let gvcf_backed = gvcf.is_some();

        let root_id = navigator_analysis::haplo::find_node(&tree, node_query).ok_or_else(|| {
            let t = match dna {
                DnaType::Y => "Y",
                DnaType::Mt => "mtDNA",
            };
            AppError::Import(format!("node '{node_query}' not found in the {t} tree"))
        })?;
        let root = tree.nodes.get(&root_id).map(|n| n.name.clone()).unwrap_or_default();
        let subtree = navigator_analysis::haplo::subtree_report(&tree, &calls, root_id, max_depth);

        // GVCF per-marker evidence (Y with a sidecar): ungated DP/AD/GQ, off-thread.
        let evidence = match gvcf {
            Some(gvcf) => {
                let positions: HashSet<i64> = subtree.iter().map(|r| r.snp.position).collect();
                let c = contig.to_string();
                tokio::task::spawn_blocking(move || navigator_analysis::gvcf::read_site_evidence(&gvcf, &c, &positions))
                    .await??
            }
            None => HashMap::new(),
        };

        let rows = subtree
            .into_iter()
            .map(|r| {
                let pos = r.snp.position;
                // Either allele being multi-base (or empty) means this isn't a clean SNV.
                let is_indel = r.snp.derived.chars().count() != 1 || r.snp.ancestral.chars().count() != 1;
                let (source, dp, ad, gq) = match evidence.get(&pos).copied() {
                    Some(e) if e.refblock => ("gvcf_refblock", None, None, e.gq),
                    Some(e) => ("gvcf_variant", e.dp, e.ad, e.gq),
                    None if gvcf_backed => ("gvcf", None, None, None),
                    None => ("pileup", None, None, None),
                };
                // The conditions are orthogonal — an uncalled indel is both — so compose the tags
                // rather than report only the first. Picking one let "indel/MNV" hide a no-call.
                let mut tags: Vec<&str> = Vec::new();
                if is_indel {
                    tags.push("indel/MNV");
                }
                if source == "gvcf_refblock" {
                    tags.push("hom-ref block");
                }
                if r.snp.state == CallState::NoCall {
                    tags.push("no call");
                }
                let note = tags.join("; ");
                BranchRow {
                    node: r.node,
                    parent: r.parent,
                    marker: r.snp.name,
                    position: pos,
                    ancestral: r.snp.ancestral,
                    derived: r.snp.derived,
                    observed_base: r.snp.base,
                    state: r.snp.state,
                    ad,
                    dp,
                    gq,
                    source,
                    note,
                }
            })
            .collect();

        Ok(BranchReport {
            dna,
            root,
            contig: contig.to_string(),
            gvcf_backed,
            rows,
        })
    }

    /// Subject-keyed [`branch_report`](Self::branch_report) for the UI: resolves the alignment to
    /// genotype `dna` against, then builds the report. `Ok(None)` when the subject has no alignment.
    pub async fn branch_report_for_subject(
        &self,
        guid: SampleGuid,
        dna: DnaType,
        node_query: &str,
        max_depth: Option<usize>,
    ) -> Result<Option<BranchReport>, AppError> {
        let Some(alignment_id) = self.pick_alignment_for(guid, dna).await? else {
            return Ok(None);
        };
        Ok(Some(self.branch_report(alignment_id, dna, node_query, max_depth).await?))
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
        // Full build: genotype any alignment whose panel dosages aren't cached yet.
        self.build_autosomal_profile_inner(biosample_guid, false).await
    }

    /// **Progressive refresh** of the autosomal consensus (progressive-consensus, docs §7.17): reduce
    /// over the per-source dosages that are **already available** — every chip / WGS-VCF (which resolve
    /// cheaply with no decode) plus any alignment whose panel dosages are *cached*
    /// ([`Self::cached_alignment_panel_dosages`]) — **without** decoding an uncached alignment. Cheap
    /// and safe to call after every import; alignments get their dosages populated separately by the
    /// panel batch-process mode, and the next refresh folds them in. Returns the refreshed profile, or
    /// `Ok(None)` when the subject has no available autosomal source yet.
    pub async fn refresh_autosomal_consensus(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<DiploidProfile>, AppError> {
        self.build_autosomal_profile_inner(biosample_guid, true).await.map(Some).or_else(|e| match e {
            // "no source" isn't an error for a refresh — the subject just has nothing cached yet.
            AppError::Import(_) => Ok(None),
            other => Err(other),
        })
    }

    /// **Panel batch-process mode** (progressive-consensus, docs §7.17): genotype one alignment at
    /// the full-1240k IBD panel and **cache** the dosages ([`Self::ibd_panel_dosages`]) — the
    /// expensive per-source step (a whole-genome decode) that populates the consensus progressively.
    /// Returns the number of panel sites genotyped. **Does not** refresh the consensus — the caller
    /// refreshes **once** after a batch (reconciling millions of observations per source is wasted
    /// work if repeated per alignment); use [`Self::refresh_autosomal_consensus`] at the batch
    /// boundary. If the dosages are already cached this is a cheap read.
    pub async fn genotype_panel_for_alignment(&self, alignment_id: i64) -> Result<usize, AppError> {
        Ok(self.ibd_panel_dosages(IbdSource::Alignment(alignment_id)).await?.len())
    }

    /// The subject's best alignment for panel genotyping, by **callable quality**:
    /// `genome_territory × pct_10x × (1 − pct_exc_mapq)` — well-mapped, diploid-callable bases. Build-
    /// agnostic (the IBD panel re-keys GRCh37/38 as well as CHM13), so it picks the cleanest
    /// whole-genome WGS over a deep-but-targeted or too-shallow test. Requires a recorded BAM/CRAM and
    /// a cached coverage artifact; `None` when the subject has neither.
    async fn best_callable_alignment(&self, biosample_guid: SampleGuid) -> Result<Option<i64>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut best: Option<(f64, i64)> = None;
        for a in &alignments {
            if a.bam_path.is_none() {
                continue;
            }
            let Some(cov) = self.cached_coverage(a.id).await? else {
                continue;
            };
            let score = cov.genome_territory as f64 * cov.pct_10x * (1.0 - cov.pct_exc_mapq);
            if best.as_ref().map_or(true, |(s, _)| score > *s) {
                best = Some((score, a.id));
            }
        }
        Ok(best.map(|(_, id)| id))
    }

    /// **Panel batch-process mode, subject-level** (progressive-consensus, docs §7.17): genotype the
    /// subject's single **best-callable** alignment ([`Self::best_callable_alignment`]) at the 1240k
    /// panel and refresh the autosomal consensus **once**. Chips and WGS-VCFs need no genotyping —
    /// they resolve into the consensus during the refresh — so this pays at most one whole-genome
    /// decode per subject (vs one per redundant same-person alignment). Returns
    /// `(alignment_id, sites)`, or `None` when the subject has no callable alignment (its chips/VCFs
    /// still get folded into the consensus).
    pub async fn genotype_panel_for_subject(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<(i64, usize)>, AppError> {
        let picked = if let Some(aln) = self.best_callable_alignment(biosample_guid).await? {
            let sites = self.genotype_panel_for_alignment(aln).await?;
            Some((aln, sites))
        } else {
            None
        };
        // Reconcile once — folds the freshly-cached alignment (if any) plus every chip / WGS-VCF.
        let _ = self.refresh_autosomal_consensus(biosample_guid).await?;
        Ok(picked)
    }

    async fn build_autosomal_profile_inner(
        &self,
        biosample_guid: SampleGuid,
        cached_alignments_only: bool,
    ) -> Result<DiploidProfile, AppError> {
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

        // One source per WGS alignment (panel-genotyped, cached per alignment). The IBD panel carries
        // every build's coordinates, so `ibd_panel_dosages` genotypes a CHM13 alignment at its native
        // loci and a GRCh37/GRCh38 alignment at that build's loci, re-keying the result to canonical
        // CHM13. A build the panel doesn't cover yields no genotypes and is skipped downstream.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            // Progressive refresh reduces over cached dosages only — an uncached alignment is skipped
            // (its dosages get populated by the panel batch mode), never decoded inline here.
            let dosages = if cached_alignments_only {
                self.cached_alignment_panel_dosages(a.id).await?
            } else {
                match self.ibd_panel_dosages(IbdSource::Alignment(a.id)).await {
                    Ok(g) => Some(g),
                    Err(e) => {
                        last_err = Some(e);
                        None
                    }
                }
            };
            if let Some(gts) = dosages {
                let obs = to_obs(gts);
                if !obs.is_empty() {
                    sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
                }
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

        // One source per **genome-wide** imported variant set — a WGS VCF or a CompleteGenomics
        // masterVar (no alignment needed; resolved to panel dosages with unlisted sites taken as
        // hom-reference). Only `WgsShortRead`/`WgsLongRead`: that hom-ref default is valid solely
        // for a source that genotyped the whole genome — a targeted Big Y (`TargetedNgs`) or Sanger
        // panel lists only a handful of sites and must NOT imply hom-ref everywhere else.
        let vsets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for set in &vsets {
            if !matches!(set.source_type, SourceType::WgsShortRead | SourceType::WgsLongRead) {
                continue;
            }
            match self.variant_set_panel_dosages(set).await {
                Ok(gts) => {
                    let obs = to_obs(gts);
                    if !obs.is_empty() {
                        sources.push((set.source_label.clone(), set.source_type, obs));
                    }
                }
                Err(e) => last_err = Some(e),
            }
        }

        // One source per imported **external autosomal call set** (a trusted 1240K EIGENSTRAT set —
        // GATK4 / pileupCaller). Resolved to CHM13 panel dosages at import and stored, so it pools in
        // with no CRAM decode (available to both the full build and the progressive refresh).
        for row in navigator_store::external_panel_dosage::list_for_biosample(self.store.pool(), biosample_guid).await? {
            match serde_json::from_str::<Vec<SiteGenotype>>(&row.dosages) {
                Ok(gts) => {
                    let obs = to_obs(gts);
                    if !obs.is_empty() {
                        sources.push((row.source_label, SourceType::Imported, obs));
                    }
                }
                Err(e) => last_err = Some(AppError::Import(format!("decoding external panel dosages: {e}"))),
            }
        }

        if sources.is_empty() {
            if let Some(e) = last_err {
                return Err(e); // e.g. the IBD panel asset isn't built yet
            }
            if cached_alignments_only {
                // A progressive refresh with nothing available yet: don't persist an empty consensus
                // (refresh_autosomal_consensus maps this to Ok(None)).
                return Err(AppError::Import("no cached autosomal sources yet".into()));
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

    /// DecodingUs mtDNA tree-with-variants JSON from our AppView (`/api/v1/mt-tree/full`), host
    /// from [`decodingus_appview_url`]. Same schema as the Y tree; coordinates are keyed by build,
    /// but the mt tree currently carries only `hs1` (CHM13 `chrM`) positions — a *rotation* of rCRS
    /// (~577, plus local indels), so callers must remap onto rCRS via [`mt_tree_rcrs`]. On-disk
    /// cached like the other trees.
    pub(crate) async fn fetch_decodingus_mt_tree(&self) -> Result<String, AppError> {
        let url = format!("{}/api/v1/mt-tree/full", decodingus_appview_url());
        self.fetch_tree(&url, "decodingus-mttree.json").await
    }

    /// The **mtDNA placement tree in rCRS coordinates**, with the provider tag. Honors the
    /// configured Y-tree provider (the Preferences toggle / `NAVIGATOR_Y_TREE_PROVIDER`): when it's
    /// set to FTDNA, use the FTDNA mt tree (already rCRS) directly. Otherwise prefer the DecodingUs
    /// mt tree remapped from its native `hs1` (CHM13 `chrM`) positions onto rCRS — so it drops
    /// straight into the existing rCRS mt pipeline (FASTA/chip sources and the `chrM` genotyper all
    /// speak rCRS) — and still fall back to FTDNA when the DecodingUs tree or the CHM13 `chrM` needed
    /// to build the remap is unavailable.
    async fn mt_tree_rcrs(&self) -> Result<(navigator_analysis::haplo::HaploTree, &'static str), AppError> {
        if !matches!(y_tree_provider(), YTreeProvider::Ftdna) {
            if let Some(tree) = self.decodingus_mt_tree_rcrs().await {
                return Ok((tree, "decodingus"));
            }
        }
        let json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&json).map_err(AppError::Import)?;
        Ok((tree, "ftdna"))
    }

    /// The DecodingUs mt tree parsed and remapped from `hs1` (CHM13 `chrM`) coordinates onto rCRS.
    /// `None` (→ FTDNA fallback) when the tree can't be fetched, or the CHM13 reference isn't cached
    /// to build the `hs1`↔rCRS map. Best-effort so an offline / reference-less workspace still works.
    async fn decodingus_mt_tree_rcrs(&self) -> Option<navigator_analysis::haplo::HaploTree> {
        let json = self.fetch_decodingus_mt_tree().await.ok()?;
        let mut tree = navigator_analysis::haplo::parse_decodingus_json(&json, "hs1").ok()?;
        let hs1_to_rcrs = self.hs1_to_rcrs_mt_map().await?;
        // Remap each defining locus from hs1 (CHM13 chrM) to rCRS; drop any that don't map (indel
        // regions near the rotation wrap). An emptied node still exists in the topology.
        for node in tree.nodes.values_mut() {
            node.loci.retain_mut(|l| match hs1_to_rcrs.get(&l.position) {
                Some(&r) => {
                    l.position = r;
                    true
                }
                None => false,
            });
        }
        Some(tree)
    }

    /// The `hs1` (CHM13 `chrM`, 1-based) → rCRS (1-based) position map, memoized for the process.
    /// Built by aligning the bundled rCRS to the cached CHM13 reference's `chrM` (rotation-aware).
    /// `None` when the CHM13 reference isn't cached (never forces a multi-GB download for this).
    async fn hs1_to_rcrs_mt_map(&self) -> Option<HashMap<i64, i64>> {
        static MAP: std::sync::OnceLock<Option<HashMap<i64, i64>>> = std::sync::OnceLock::new();
        if let Some(m) = MAP.get() {
            return m.clone();
        }
        let reference = self.gateway.cached_reference("chm13v2.0")?;
        let pairs = tokio::task::spawn_blocking(move || {
            navigator_analysis::reader::read_contig_sequence(&reference, "chrM").ok().map(|chrm| {
                let chrm = String::from_utf8_lossy(&chrm).into_owned();
                navigator_analysis::mtvariants::mt_position_map(navigator_analysis::mtvariants::rcrs(), &chrm)
            })
        })
        .await
        .ok()
        .flatten();
        // Invert to hs1(chrM)→rCRS; both stored 1-based (mt_position_map yields 0-based pairs).
        let map = pairs.map(|p| p.into_iter().map(|(r, c)| (c as i64 + 1, r as i64 + 1)).collect());
        MAP.get_or_init(|| map).clone()
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
    /// Force a fresh pull of the haplotrees on the next placement: clear the session memo AND delete
    /// the on-disk tree caches, so a corrected AppView tree (e.g. a polarity fix) is picked up without
    /// an app restart. Observation-first profiles then re-interpret against the new tree on read — no
    /// re-genotyping. Returns the number of cache files removed.
    pub async fn refresh_trees(&self) -> Result<usize, AppError> {
        tree_memo().lock().unwrap().clear();
        let mut removed = 0usize;
        if let Some(dir) = tree_cache_path("_").parent().map(|p| p.to_path_buf()) {
            if let Ok(entries) = std::fs::read_dir(&dir) {
                for e in entries.flatten() {
                    let p = e.path();
                    if p.extension().and_then(|x| x.to_str()) == Some("json")
                        && std::fs::remove_file(&p).is_ok()
                    {
                        removed += 1;
                    }
                }
            }
        }
        Ok(removed)
    }

    async fn fetch_tree(&self, url: &str, cache_file: &str) -> Result<String, AppError> {
        // Session memo: the Y/mt haplotrees are 4–121 MB and each placement consults them several
        // times (per alignment, per vendor set, and for the polarity map). A single genome-consensus
        // build alone would otherwise re-read/re-validate them repeatedly, and a *stale*-cache refresh
        // blocks on the network. Resolve each tree at most once per process and serve every later call
        // from memory — trees are effectively static within a session, so this is the batch's biggest
        // win (a project pass was spending minutes per subject re-fetching the 121 MB FTDNA tree).
        // Keyed by the *resolved* path, not the bare file name: `NAVIGATOR_TREE_DIR` can point the
        // same `cache_file` at different trees, and a name-keyed memo would serve the first one for
        // the rest of the process.
        let path = tree_cache_path(cache_file);
        let key = path.to_string_lossy().into_owned();
        let memo = tree_memo();
        if let Some(hit) = memo.lock().unwrap().get(&key).cloned() {
            return Ok(hit);
        }

        let cached = std::fs::read_to_string(&path).ok().filter(|c| !c.trim().is_empty());
        // A fresh (within-TTL) on-disk cache short-circuits the network entirely.
        let fresh = cached.is_some() && tree_cache_is_fresh(&path);
        let json = if fresh {
            cached.expect("fresh implies present")
        } else {
            // Stale or absent → *conditional*, time-bounded refresh. When we have a cached copy and
            // its stored ETag we send `If-None-Match`: an unchanged tree comes back as a tiny `304`
            // (a few bytes) instead of re-streaming the full ~60–127 MB body — the curated tree
            // changes only every week or so, so most refreshes are 304s. Any failure (connect/
            // timeout, a non-2xx/304 status, or a body read cut short — the whole-request timeout
            // also covers streaming the body, see [`TREE_DOWNLOAD_TIMEOUT`]) falls back to the cached
            // copy when present; only a first-ever fetch with no cache errors.
            enum TreeFetch {
                NotModified,
                Modified { body: String, etag: Option<String> },
            }
            let etag_path = tree_etag_path(&path);
            let prior_etag = cached
                .as_ref()
                .and_then(|_| std::fs::read_to_string(&etag_path).ok())
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty());

            let fetched = async {
                let mut req = self.auth.http.get(url).timeout(TREE_DOWNLOAD_TIMEOUT);
                if let Some(etag) = &prior_etag {
                    req = req.header(reqwest::header::IF_NONE_MATCH, etag.as_str());
                }
                let resp = req
                    .send()
                    .await
                    .map_err(|e| AppError::Import(format!("downloading {url}: {e}")))?;
                if resp.status() == reqwest::StatusCode::NOT_MODIFIED {
                    return Ok(TreeFetch::NotModified);
                }
                let resp = resp
                    .error_for_status()
                    .map_err(|e| AppError::Import(format!("downloading {url}: {e}")))?;
                let etag = resp
                    .headers()
                    .get(reqwest::header::ETAG)
                    .and_then(|v| v.to_str().ok())
                    .map(str::to_owned);
                let body = resp
                    .text()
                    .await
                    .map_err(|e| AppError::Import(format!("reading {url}: {e}")))?;
                Ok(TreeFetch::Modified { body, etag })
            }
            .await;

            match fetched {
                Ok(TreeFetch::NotModified) => match cached {
                    // Still current. Re-write the identical bytes to bump the cache mtime so the TTL
                    // resets and later runs short-circuit without even a conditional request.
                    Some(body) => {
                        let _ = std::fs::write(&path, &body);
                        body
                    }
                    None => return Err(AppError::Import(format!("{url}: 304 with no cached tree"))),
                },
                Ok(TreeFetch::Modified { body, etag }) => {
                    if let Some(parent) = path.parent() {
                        let _ = std::fs::create_dir_all(parent);
                    }
                    let _ = std::fs::write(&path, &body);
                    // Persist the new ETag (or clear a stale one) for the next conditional refresh.
                    match etag {
                        Some(e) => {
                            let _ = std::fs::write(&etag_path, e);
                        }
                        None => {
                            let _ = std::fs::remove_file(&etag_path);
                        }
                    }
                    body
                }
                Err(e) => match cached {
                    Some(stale) => {
                        eprintln!("tree refresh failed ({e}); using the cached copy at {}", path.display());
                        stale
                    }
                    None => return Err(e),
                },
            }
        };
        memo.lock().unwrap().insert(key, json.clone());
        Ok(json)
    }

    /// Assign an mtDNA haplogroup directly from an alignment's chrM reads (FTDNA mt tree),
    /// the BAM-based counterpart to [`assign_mtdna_haplogroup`]. Requires a GRCh38/rCRS
    /// chrM (the tree is in rCRS coordinates).
    pub async fn assign_mtdna_haplogroup_from_alignment(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();

        // Prefer an external (sidecar-GVCF) mt call over re-walking the CRAM — same rationale as the
        // Y path (see `assign_y_haplogroup`); this is the guard the unguarded single-alignment
        // "Full Analysis" was missing, so an internal re-run no longer overwrites the GATK4 mt call.
        if let Some(guid) = bio {
            if let Some(call) = self.preferred_external_call(guid, DnaType::Mt, alignment_id).await? {
                return Ok(assignment_from_call(&call));
            }
        }

        self.assign_mtdna_haplogroup_walk(alignment_id, bio).await
    }

    /// The internal-caller mtDNA placement: place chrM against the FTDNA mt tree and record under the
    /// walk key (`aln:{id}:mt`, `NavigatorWalk`), skipping the re-score when the fingerprint is
    /// unchanged. Split out of [`assign_mtdna_haplogroup_from_alignment`] so [`compare_callers`] can
    /// force the internal walk even when an external call is preferred.
    pub(crate) async fn assign_mtdna_haplogroup_walk(
        &self,
        alignment_id: i64,
        bio: Option<SampleGuid>,
    ) -> Result<HaploAssignment, AppError> {
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
                CallProvenance::NavigatorWalk,
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
        // Resolve the reference for decode (see alignment_reference_for_decode): required for a CRAM,
        // None for a BAM. The chrM pileup finds a second allele from reads; it needs no reference base.
        let (bam, reference) = self.alignment_reference_for_decode(alignment_id).await?;
        tokio::task::spawn_blocking(move || {
            navigator_analysis::guard_walk("heteroplasmy", || {
                heteroplasmy::detect_heteroplasmy(&bam, "chrM", &HeteroplasmyParams::default(), reference.as_deref())
            })
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
        // Auto-download the prebuilt panels on first use (no `panelbuild`). The super-pop panel is
        // required; PCA + fine frequencies are optional (best-effort — the feature degrades if absent).
        self.ensure_ancestry_asset(build, &ancestry_panel_path(build)).await?;
        let _ = self.ensure_ancestry_asset(build, &ancestry_pca_path(build)).await;
        let _ = self.ensure_ancestry_asset(build, &ancestry_freq_global_path(build)).await;
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
        let fine_bytes = optional(ancestry_freq_global_path(build));
        let ancient_bytes = optional(ancestry_freq_ancient_path(build));

        let (result, ancient, fine) = tokio::task::spawn_blocking(move || {
            let mut result = ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            let fine = fine_bytes
                .and_then(|b| ancestry_analysis::AncestryPanel::from_bytes(&b).ok())
                .map(|fp| ancestry_analysis::estimate_fine_admixture(&genotypes, &fp, &reference_version));
            if let Some(pca) = pca_bytes.and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok()) {
                result.pca_coordinates = Some(ancestry_analysis::project_pca(&genotypes, &pca));
            }
            // Deep (ancient) ancestry is NOT computed here — it is the separate `estimate_deep_ancestry`
            // path. NB it is *not* a heavier genotyping pass: it reads the SAME cached autosomal
            // consensus these modern/fine estimators use (the consensus is the full ~1.15M-site 1240k
            // IBD-panel union, not a 20k subset), and just intersects the larger qpAdm f4 panel against
            // it (docs/design/ancient-ancestry-rebuild.md §7.14). Kept out of this hot path because the
            // qpAdm fit is a distinct, on-demand model.
            let ancient: Option<AncestryResult> = None;
            let _ = &ancient_bytes;
            (result, ancient, fine)
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
        for extra in [ancient.as_ref(), fine.as_ref()].into_iter().flatten() {
            ancestry_result::upsert(self.store.pool(), biosample_guid, CONSENSUS_SOURCE_ID, extra).await?;
        }
        Ok(result)
    }

    /// **Deep (ancient) ancestry via qpAdm** (docs/design/ancient-ancestry-rebuild.md §7.14, §7.16) —
    /// the validated WHG / EEF / Steppe breakdown.
    ///
    /// Consumes the subject's **autosomal consensus** ([`Self::build_autosomal_profile`]) — the same
    /// multi-source object modern/fine ancestry uses. The consensus is built by the IBD panel's
    /// per-build resolver, so it already pools every source **re-keyed to canonical CHM13**: WGS on
    /// any reference (GRCh37/38 as well as CHM13) *and* consumer chips, with no alignment required.
    /// That is why deep ancestry works multi-reference and chip-only — it inherits the frontend the
    /// other estimates share. It fits `target = Σ wᵢ·sourcesᵢ` by qpAdm f4 over the (now
    /// CHM13-canonical, §7.16) qpAdm panel and persists the result under the consensus pseudo-source.
    ///
    /// **Requires** the autosomal consensus (errors if absent — build it via the Autosomal tab, same
    /// contract as modern ancestry). Given it, this is a fast fit over cached genotypes; the heavy
    /// full-1240k genotyping happens once, in the shared consensus build.
    ///
    /// Returns `Ok(None)` when the feature is gated off, the asset is not installed, the consensus has
    /// no autosomal calls, or the deep model does not apply (non-European / model rejected / infeasible
    /// weights). `None` persists nothing — keeping an inapplicable breakdown off the UI *and* out of
    /// the PDS.
    pub async fn estimate_deep_ancestry(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<AncestryResult>, AppError> {
        if !crate::ANCIENT_ANCESTRY_ENABLED {
            return Ok(None);
        }
        let build = ReferenceBuild::Chm13v2;
        let reference_version = "chm13v2.0".to_string();
        // Auto-download the prebuilt qpAdm + super-pop panels on first use (no `panelbuild`). qpAdm is
        // best-effort (deep ancestry is simply unavailable if it can't be fetched); the super panel is required.
        let _ = self.ensure_ancestry_asset(build, &ancestry_qpadm_path(build)).await;
        self.ensure_ancestry_asset(build, &ancestry_panel_path(build)).await?;
        let qpadm_path = ancestry_qpadm_path(build);
        let Some(qpadm_bytes) = read_verified_asset(build, &qpadm_path)? else {
            return Ok(None); // asset not installed → deep ancestry unavailable
        };
        let panel = AncestryPanel::from_bytes(&qpadm_bytes)?;
        let super_path = ancestry_panel_path(build);
        let super_bytes = read_verified_asset(build, &super_path)?
            .ok_or_else(|| AppError::AncestryPanelMissing(super_path.clone()))?;
        let super_panel = AncestryPanel::from_bytes(&super_bytes)?;

        // The pooled autosomal consensus (all sources, any build + chips, canonical CHM13, full 0/1/2
        // dosages). **Required**, not built on demand — same contract as modern ancestry: the heavy
        // build runs through the Autosomal-tab flow (with progress), and this is a fast read over it.
        // Both the scope gate and the qpAdm fit read the *same* genotypes; every panel here is
        // CHM13-canonical (§7.16), so no per-site re-keying is needed.
        let profile = self.cached_autosomal_profile(biosample_guid).await?.ok_or_else(|| {
            AppError::Import("build the autosomal consensus first (Autosomal tab) before estimating deep ancestry".into())
        })?;
        let genotypes = consensus_genotypes(&profile);
        if genotypes.is_empty() {
            return Ok(None); // consensus exists but has no autosomal calls
        }

        let n_pops = panel.populations.len();
        let result = tokio::task::spawn_blocking(move || {
            // Scope gate: a deep three-way is a West-Eurasian model, so refuse non-Europeans. Uses
            // the same consensus genotypes, scored against the (canonical) super-pop AIM panel.
            let modern = ancestry_analysis::estimate_admixture(&genotypes, &super_panel, &reference_version);
            if ancestry_analysis::west_eurasian_share(&modern) < 50.0 {
                return None;
            }
            // Committed Patterson layout: the first three populations are the sources, the rest the
            // sister outgroups.
            let sources: Vec<usize> = (0..3).collect();
            let outgroups: Vec<usize> = (3..n_pops).collect();
            ancestry_analysis::estimate_qpadm_ancestry(
                &genotypes,
                &panel,
                &sources,
                &outgroups,
                &modern,
                &reference_version,
            )
        })
        .await?;

        if let Some(r) = &result {
            ancestry_result::upsert(self.store.pool(), biosample_guid, CONSENSUS_SOURCE_ID, r).await?;
        }
        Ok(result)
    }

    /// **Deep-ancestry stability diagnostic** — the §3.4 validation gates, on a real subject.
    ///
    /// Fits the ancient mixture repeatedly over different *views* of the same person: the pooled
    /// consensus, each contributing source on its own (a 30× WGS and a consumer chip are genotyped
    /// by completely different means, so agreeing across them is the strongest evidence the estimate
    /// tracks the donor and not the assay), and random subsets of the sites.
    ///
    /// This is the test the previous implementation failed most spectacularly — the same person came
    /// out WHG 72.6% from the consensus and WHG 7.4% from their own 28× BAM — so it is the one worth
    /// being able to re-run on demand. Rows are diagnostics, never persisted or published; the
    /// `reported` flag records whether the shipping estimator would have accepted that fit.
    ///
    /// Every row reports its dispersion even when the applicability gate rejects it, so a rejection
    /// can be read as a magnitude rather than taken on faith.
    pub async fn ancient_ancestry_stability(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Vec<AncientFitRow>, AppError> {
        // Build the consensus on demand — this is a diagnostic, and requiring the caller to have
        // clicked through the GUI first would make it useless from the CLI.
        let profile = match self.cached_autosomal_profile(biosample_guid).await? {
            Some(p) => p,
            None => self.build_autosomal_profile(biosample_guid).await?,
        };
        let build = ReferenceBuild::Chm13v2;
        let path = ancestry_freq_ancient_path(build);
        let bytes =
            read_verified_asset(build, &path)?.ok_or_else(|| AppError::AncestryPanelMissing(path.clone()))?;
        let panel = AncestryPanel::from_bytes(&bytes)?;
        // The super-pop panel too: deep ancestry is scoped by the modern estimate, so each view has
        // to be scored by both models or the diagnostic wouldn't be reproducing the shipped policy.
        let super_path = ancestry_panel_path(build);
        let super_bytes = read_verified_asset(build, &super_path)?
            .ok_or_else(|| AppError::AncestryPanelMissing(super_path.clone()))?;
        let super_panel = AncestryPanel::from_bytes(&super_bytes)?;

        // The distinct sources that contributed a call anywhere in the consensus.
        let mut source_labels: Vec<String> = profile
            .variants
            .iter()
            .flat_map(|v| v.sources.iter().map(|s| s.label.clone()))
            .collect();
        source_labels.sort();
        source_labels.dedup();

        tokio::task::spawn_blocking(move || {
            let mut rows = Vec::new();
            let mut fit = |label: String, genotypes: &[SiteGenotype]| {
                let modern = ancestry_analysis::estimate_admixture(genotypes, &super_panel, "chm13v2.0");
                if let Some(r) = ancestry_analysis::ancient_admixture_fit(genotypes, &panel, "chm13v2.0") {
                    rows.push(AncientFitRow {
                        label,
                        sites: r.snps_with_genotype,
                        dispersion: r.fit_distance.unwrap_or(f64::NAN),
                        european: ancestry_analysis::west_eurasian_share(&modern),
                        // The shipping estimator's own verdict — not a re-derivation of it.
                        reported: ancestry_analysis::estimate_ancient_admixture(
                            genotypes,
                            &panel,
                            &modern,
                            "chm13v2.0",
                        )
                        .is_some(),
                        components: r
                            .components
                            .iter()
                            .map(|c| (c.population_code.clone(), c.percentage))
                            .collect(),
                    });
                }
            };

            let consensus = consensus_genotypes(&profile);
            fit("consensus (pooled)".to_string(), &consensus);

            // Sites a chip actually called — the intersection target for the refit below. The
            // stability failure has WGS reporting ~80% Steppe where the chips report ~58%; this asks
            // whether the split is *which sites* each technology reaches (WGS scores ~19.7k, a chip
            // ~5–8k) or the calls themselves. If a WGS source restricted to chip-covered sites moves
            // toward the chip answer, the extra WGS-only sites carry the bias; if it stays put, the
            // WGS dosages do.
            let chip_sites: std::collections::HashSet<String> = profile
                .variants
                .iter()
                .filter(|v| {
                    v.sources
                        .iter()
                        .any(|s| matches!(s.source_type, SourceType::Chip) && s.dosage >= 0)
                })
                .map(|v| v.name.clone())
                .collect();

            // Diagnostic dump (NAVIGATOR_ANCIENT_DUMP=<path>): per consensus site, whether a chip
            // covers it, the pooled dosage, and the three source frequencies. Lets us see directly
            // what makes the non-chip sites favour Steppe once every intrinsic site property
            // (MAF/strand/polarity/ts-tv) has been ruled out.
            if let Ok(dump_path) = std::env::var("NAVIGATOR_ANCIENT_DUMP") {
                let want = std::env::var("NAVIGATOR_ANCIENT_ALN").unwrap_or_else(|_| "#9".into());
                let freq: std::collections::HashMap<(&str, i64), &Vec<f32>> =
                    panel.sites.iter().map(|s| ((s.contig.as_str(), s.position), &s.freqs)).collect();
                let mut out = String::from("chip\taln_dosage\twhg\tanf\tsteppe\n");
                for v in &profile.variants {
                    let Some(f) = freq.get(&(v.contig.as_str(), v.position)) else { continue };
                    if f.len() != 3 {
                        continue;
                    }
                    // One clean alignment's own genotype at this ancient-panel site.
                    let d = v
                        .sources
                        .iter()
                        .find(|s| {
                            matches!(s.source_type, SourceType::WgsShortRead | SourceType::WgsLongRead)
                                && s.dosage >= 0
                                && s.label.contains(&want)
                        })
                        .map_or(-1, |s| s.dosage as i32);
                    let chip = u8::from(chip_sites.contains(&v.name));
                    out.push_str(&format!("{}\t{}\t{:.4}\t{:.4}\t{:.4}\n", chip, d, f[0], f[1], f[2]));
                }
                let _ = std::fs::write(&dump_path, out);
            }

            // Strand-ambiguous SNPs (A/T, C/G): ref and alt are Watson–Crick complements, so which
            // allele the panel counted as "alt" can't be recovered from the alleles alone. When a
            // panel built from one dataset is genotyped against reads oriented by another, these are
            // the sites that silently invert — the classic merge bias, and one chips routinely drop.
            let is_ambiguous = |v: &navigator_domain::consensus::DiploidVariant| -> bool {
                matches!(
                    (v.reference.as_str(), v.alternate.as_str()),
                    ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C")
                )
            };

            // One source's own observed dosages, restricted to the variants the predicate keeps.
            let build_single =
                |label: &str, keep: &dyn Fn(&navigator_domain::consensus::DiploidVariant) -> bool| -> Vec<SiteGenotype> {
                    profile
                        .variants
                        .iter()
                        .filter(|v| keep(v))
                        .filter_map(|v| {
                            let obs = v.sources.iter().find(|s| s.label.as_str() == label)?;
                            (obs.dosage >= 0).then(|| SiteGenotype {
                                name: v.name.clone(),
                                contig: v.contig.clone(),
                                position: v.position,
                                reference_allele: v.reference.clone(),
                                alternate_allele: v.alternate.clone(),
                                ploidy: 2,
                                dosage: obs.dosage as i32,
                                gq: 0,
                                depth: 0,
                                ref_depth: 0,
                                alt_depth: 0,
                                pls: Vec::new(),
                                gt: None,
                                allele_depths: None,
                            })
                        })
                        .collect()
                };

            // Each source alone: take that source's own observed dosage at each site. For a WGS
            // source, also refit it three ways to localize the stability bias:
            //   ∩chip   — sites the chips cover (does WGS match the chip answer there?)
            //   ∁chip   — the WGS-only complement (do those sites carry the bias?)
            //   ¬ambig  — all sites minus strand-ambiguous A/T,C/G (does dropping them fix it?)
            for label in &source_labels {
                let is_chip = profile.variants.iter().any(|v| {
                    v.sources
                        .iter()
                        .any(|s| s.label.as_str() == label && matches!(s.source_type, SourceType::Chip))
                });
                fit(format!("source: {label}"), &build_single(label, &|_| true));
                if !is_chip {
                    fit(format!("source: {label} ∩chip"), &build_single(label, &|v| chip_sites.contains(&v.name)));
                    fit(format!("source: {label} ∁chip"), &build_single(label, &|v| !chip_sites.contains(&v.name)));
                    fit(format!("source: {label} ¬ambig"), &build_single(label, &|v| !is_ambiguous(v)));
                }
            }

            // Density: deterministic thinning of the pooled consensus. A well-conditioned fit barely
            // moves when half the evidence is removed; an over-fit one lurches.
            for (keep, label) in [(2usize, "consensus ÷2 sites"), (4, "consensus ÷4 sites")] {
                let thinned: Vec<SiteGenotype> =
                    consensus.iter().step_by(keep).cloned().collect();
                fit(label.to_string(), &thinned);
            }
            rows
        })
        .await
        .map_err(AppError::from)
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
        self.ensure_ancestry_asset(build, &ancestry_panel_path(build)).await?;
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

    /// The cached ROH result for a subject, if one was computed from the **current** autosomal
    /// consensus (signature = the consensus's `last_reconciled_at`). `None` if absent or stale (the
    /// consensus was rebuilt since). Cheap — a cache read, no genotyping or HMM.
    pub async fn cached_roh(&self, biosample_guid: SampleGuid) -> Result<Option<RohResult>, AppError> {
        let Some(row) = consensus_profile::get(self.store.pool(), biosample_guid, "Auto").await? else {
            return Ok(None);
        };
        let Some(r) = consensus_roh::get(self.store.pool(), biosample_guid).await? else {
            return Ok(None);
        };
        if r.consensus_sig == row.last_reconciled_at {
            Ok(Some(serde_json::from_str(&r.roh)?))
        } else {
            Ok(None) // computed from an older consensus — stale
        }
    }

    /// Detect runs of homozygosity from the subject's **consensus** — no BAM walk. Returns the cached
    /// result when it matches the current consensus signature; otherwise runs the 2-state autozygosity
    /// HMM over the consensus genotypes and caches it keyed to the consensus's `last_reconciled_at`.
    /// The genome-wide F_ROH and length-class breakdown are the endogamy / consanguinity signal.
    pub async fn compute_roh_from_consensus(&self, biosample_guid: SampleGuid) -> Result<RohResult, AppError> {
        let row = consensus_profile::get(self.store.pool(), biosample_guid, "Auto")
            .await?
            .ok_or_else(|| {
                AppError::Import("build the autosomal consensus first (Autosomal tab) before computing ROH".into())
            })?;
        let sig = row.last_reconciled_at.clone();

        // Cache hit (same consensus signature) → return without recomputing.
        if let Some(r) = consensus_roh::get(self.store.pool(), biosample_guid).await? {
            if r.consensus_sig == sig {
                return Ok(serde_json::from_str(&r.roh)?);
            }
        }

        let profile: DiploidProfile = serde_json::from_str(&row.payload)?;
        let genotypes = consensus_genotypes(&profile);
        let result = tokio::task::spawn_blocking(move || {
            // Per-contig max position → genetic-map lengths (CHM13 consensus space, uniform fallback).
            let mut lengths: std::collections::BTreeMap<String, i32> = std::collections::BTreeMap::new();
            for g in &genotypes {
                let e = lengths.entry(g.contig.clone()).or_insert(1);
                *e = (*e).max(g.position as i32);
            }
            let pairs: Vec<(&str, i32)> = lengths.iter().map(|(k, v)| (k.as_str(), *v)).collect();
            let gmap = crate::load_genetic_map(ReferenceBuild::Chm13v2, &pairs);
            navigator_analysis::roh::detect_roh(&genotypes, &gmap, &navigator_analysis::roh::RohConfig::default())
        })
        .await?;

        // Cache keyed to the consensus signature so it's reused until the consensus is rebuilt.
        consensus_roh::upsert(
            self.store.pool(),
            biosample_guid,
            &sig,
            &serde_json::to_string(&result)?,
            &Utc::now().to_rfc3339(),
        )
        .await?;
        Ok(result)
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
    /// Whether a subject should be scored for Y-DNA. Females have no Y chromosome, so Y placement,
    /// consensus, and the Y variant profile only produce an empty/degenerate call from mismapped
    /// chrY reads — skip them. Sex comes from `biosample.sex` (user-provided, or written back by the
    /// sex walker, which runs before the Y step). Male / Unknown / unrecorded → scored (`true`), so a
    /// low-confidence or missing inference, or an XXY subject, is never silently dropped.
    pub(crate) async fn subject_has_y_dna(&self, biosample_guid: SampleGuid) -> Result<bool, AppError> {
        let sex = biosample::get(self.store.pool(), biosample_guid)
            .await?
            .and_then(|b| b.sex);
        Ok(!matches!(sex.as_deref().map(str::trim), Some(s) if s.eq_ignore_ascii_case("female")))
    }

    pub async fn assign_y_haplogroup(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();

        // Females have no Y chromosome — don't genotype chrY or record a Y call for them.
        if let Some(guid) = bio {
            if !self.subject_has_y_dna(guid).await? {
                return Ok(HaploAssignment {
                    ranked: Vec::new(),
                    branches: Vec::new(),
                    lineage: Vec::new(),
                });
            }
        }

        // A trusted external caller (GATK4 GVCF) already placed this alignment via the sidecar fast
        // path and the user prefers it: return that call instead of re-walking the CRAM. On damaged
        // ancient DNA the walk would place a different, wrong terminal and clobber the external one.
        if let Some(guid) = bio {
            if let Some(call) = self.preferred_external_call(guid, DnaType::Y, alignment_id).await? {
                return Ok(assignment_from_call(&call));
            }
        }

        self.assign_y_haplogroup_walk(alignment_id, bio).await
    }

    /// The internal-caller Y placement: genotype chrY against the configured tree and record the call
    /// under the walk key (`aln:{id}`, `NavigatorWalk` provenance), skipping the re-score when the
    /// alignment + tree fingerprint is unchanged. Split out of [`assign_y_haplogroup`] so
    /// [`compare_callers`] can force the internal walk even when an external call is preferred.
    pub(crate) async fn assign_y_haplogroup_walk(
        &self,
        alignment_id: i64,
        bio: Option<SampleGuid>,
    ) -> Result<HaploAssignment, AppError> {
        let source_key = format!("aln:{alignment_id}");
        // Input fingerprint = alignment content hash + active Y-tree content hash. If it matches the
        // recorded call's stamp, neither the file nor the tree changed → return the recorded call
        // without re-scoring (the expensive BAM genotyping).
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
                CallProvenance::NavigatorWalk,
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
    pub(crate) async fn y_decodingus_tree_calls(
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

    /// Autosomal variant sites of a set as `(bare-contig, position, a1, a2)` reference-forward
    /// allele pairs — the input to the whole-genome IBD-panel resolve. Only chr1–chr22 (the panel
    /// is autosomal). The genotype string is turned back into an allele pair: `1/1` → `(alt, alt)`;
    /// het (`0/1`, `1/.`, or an absent genotype — one listed alt means at least one copy) →
    /// `(ref, alt)`; tri-allelic `1/2` and haploid `1` (Y/mt, never autosomal) are dropped as
    /// ambiguous for a diploid dosage.
    fn vset_autosomal_calls(set: &VariantSet) -> Vec<(String, i64, char, char)> {
        set.calls
            .iter()
            .filter_map(|c| {
                let bare = c.contig.strip_prefix("chr").unwrap_or(&c.contig);
                let n: u32 = bare.parse().ok()?;
                if !(1..=22).contains(&n) {
                    return None;
                }
                let r = c.reference.chars().next()?;
                let a = c.alternate.chars().next()?;
                let gt = c.genotype.as_deref().unwrap_or("");
                let (a1, a2) = if gt == "1/1" || gt == "1|1" {
                    (a, a)
                } else if matches!(gt, "0/1" | "1/0" | "0|1" | "1|0" | "1/." | "./1") || gt.is_empty() {
                    (r, a)
                } else {
                    return None; // 1/2 tri-allelic, haploid "1", or an unrecognized genotype
                };
                Some((bare.to_string(), c.position, a1, a2))
            })
            .collect()
    }

    /// Resolve a **genome-wide** variant set to canonical CHM13 IBD-panel dosages (unlisted panel
    /// sites ⇒ hom-reference — see [`IbdPanel::resolve_whole_genome`]). Needs the IBD panel asset.
    /// Returns an empty vec for a set with no autosomal calls (e.g. a Y-only VCF). Not cached — the
    /// resolve is cheap and, like the chip path, recomputed when the autosomal consensus is rebuilt.
    pub(crate) async fn variant_set_panel_dosages(&self, set: &VariantSet) -> Result<Vec<SiteGenotype>, AppError> {
        let calls = Self::vset_autosomal_calls(set);
        if calls.is_empty() {
            return Ok(Vec::new());
        }
        let build = set.reference_build.clone().unwrap_or_else(|| "GRCh37".to_string());
        let panel = self.load_ibd_panel().await?;
        let dosages = tokio::task::spawn_blocking(move || panel.resolve_whole_genome(&build, &calls)).await?;
        Ok(dosages)
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
            // A set with no tree-defining SNP matched carries no Y signal — e.g. an off-haplotree
            // FTDNA "Private Variants" report (novel loci only), or an autosomal/mt VCF that slipped
            // through. Recording its placeholder placement would conflict with a real per-source
            // call and collapse the donor consensus to root, so skip it.
            if assignment.ranked.first().map_or(true, |t| t.matched == 0) {
                continue;
            }
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
        // Route through the provider-correct placement (DecodingUs native multi-build, FTDNA
        // fallback with polarity normalization) — not a raw `parse_ftdna_json` of whatever tree the
        // provider returns, which fails on the DecodingUs schema. The assignment already carries the
        // root→terminal lineage evidence.
        let assignment = self.y_assignment_full(alignment_id).await?;
        let lineage = assignment.lineage.clone();
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
    pub(crate) async fn tree_base_calls(
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
    ///
    /// The result (tree-position → base) is cached as a versioned analysis artifact keyed by the
    /// queried **site set** (a hash of the tree's positions) + contig + lift source, and
    /// invalidated by the alignment's `source_sig` (BAM/CRAM mtime:size). This is the BAM-walk
    /// chokepoint for *every* genotyping path (Y/mt placement, the variant profile, genome
    /// consensus), so a profile **rebuild** reuses the cached genotypes instead of re-walking the
    /// reads — only a changed file or a changed tree site set forces a fresh walk.
    /// Tree-locus base calls for one alignment for the **genome-consensus placement**, preferring the
    /// alignment's external sidecar GVCF (no CRAM decode) when the "prefer external caller" policy is
    /// on and the GVCF is present; otherwise the cached CRAM walk ([`base_calls`]). A drop-in for the
    /// per-alignment genotype in `place_{y,mt}_consensus`, so a preferred-external (e.g. ancient-DNA)
    /// subject's damaged CRAM is not re-walked and cannot dilute the pooled placement (Phase 2 of
    /// `docs/design/external-caller-precedence.md` §4.5). `tree_source_build` matches what `base_calls`
    /// receives — `None` for a native-build tree (DecodingUs Y, rCRS mt), the tree's build for a lift.
    pub(crate) async fn consensus_base_calls(
        &self,
        aln: &Alignment,
        contig: &str,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        if prefer_external_calls() {
            let gvcf = if contig.eq_ignore_ascii_case("chrM") {
                chr_m_gvcf_for_alignment(aln)
            } else {
                chr_y_gvcf_for_alignment(aln)
            };
            if let Some(gvcf) = gvcf {
                return self.gvcf_base_calls(aln.id, contig, &gvcf, tree, tree_source_build).await;
            }
        }
        self.base_calls(aln.id, contig, tree, tree_source_build).await
    }

    async fn base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        let targets: HashSet<i64> = tree
            .nodes
            .values()
            .flat_map(|n| n.loci.iter().map(|l| l.position))
            .collect();
        if targets.is_empty() {
            return Ok(HashMap::new());
        }

        // Cache hit → skip both the reference resolution (a CRAM would otherwise download) and the
        // read walk entirely.
        let cache_key = genotype_cache_key(contig, tree_source_build, &targets);
        if let Some(pairs) = self
            .load_analysis::<Vec<(i64, char)>>(alignment_id, GENOTYPE_KIND, &cache_key)
            .await?
        {
            return Ok(pairs.into_iter().collect());
        }

        let aln = self.alignment_or_err(alignment_id).await?;
        // Copy off a slow/removable volume to local disk first — the per-locus genotyping read is a
        // network round-trip per record otherwise (see App::localize).
        let bam = self
            .localize(Path::new(&aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?))
            .await;
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

        let lifted = self
            .lifted_targets(
                &aln.reference_build,
                reference.as_deref(),
                contig,
                &targets,
                tree_source_build,
            )
            .await?;

        // Indel loci (multi-base ancestral/derived) on the tree — genotyped separately on the native
        // chrY path (VCF left-anchored; needs the reference to normalize + know deleted bases). Their
        // resolved sentinel overlays the (meaningless) base call at the anchor. Liftover of indel
        // coordinates isn't handled, so only the native path (no lift) contributes them.
        let indel_targets: Vec<(i64, String, String)> = if contig.eq_ignore_ascii_case("chrY") {
            tree.nodes
                .values()
                .flat_map(|n| n.loci.iter())
                .filter(|l| l.derived.chars().count() > 1 || l.ancestral.chars().count() > 1)
                .map(|l| (l.position, l.ancestral.clone(), l.derived.clone()))
                .collect()
        } else {
            Vec::new()
        };

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
                tokio::task::spawn_blocking(move || {
                    let params = adaptive_haploid_params(&bam, reference.as_deref()); // HiFi -> lower min_depth
                    navigator_analysis::guard_walk("haplogroup genotyping", || {
                        let mut calls =
                            caller::call_bases_at(&bam, &resolved, &targets, &params, reference.as_deref())?;
                        if !indel_targets.is_empty() {
                            let indels = caller::call_indels_at(
                                &bam,
                                &resolved,
                                &indel_targets,
                                &params,
                                reference.as_deref(),
                            )?;
                            calls.extend(indels); // sentinel overlays the anchor's base call
                        }
                        Ok(calls)
                    })
                })
                .await??
            }
        };

        // Cache the genotypes (stamped with the BAM source_sig) so a rebuild skips the walk.
        let pairs: Vec<(i64, char)> = calls.iter().map(|(&p, &b)| (p, b)).collect();
        self.save_analysis(alignment_id, GENOTYPE_KIND, &cache_key, &pairs).await?;
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
            // `src` must name a reference build. A tree *provider* ("decodingus"/"ftdna") would
            // fall through the `differ` test to the `return Ok(None)` below, silently disabling
            // liftover — for chrM that skips the rCRS↔chrM map and miscalls every marker. Refuse
            // it here rather than answer with wrong coordinates.
            let Some(src_build) = canonical_build(src) else {
                return Err(AppError::Import(format!(
                    "lifted_targets: tree_source_build {src:?} is not a reference build \
                     (pass None when the tree is already in the query coordinate space)"
                )));
            };
            let differ = matches!(canonical_build(reference_build), Some(t) if src_build != t);
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
                navigator_analysis::guard_walk("haplogroup genotyping", || {
                    caller::call_bases_at(&bam, &qc, &set, &params, reference.as_deref())
                })
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

/// Interpret an [`navigator_domain::consensus::ObservedProfile`] against a polarity map into the
/// app's display [`ConsensusProfile`] (carrying provenance + terminal). The single place observations
/// become the interpreted view.
fn interpret_observed(
    observed: navigator_domain::consensus::ObservedProfile,
    polarity: &std::collections::BTreeMap<String, (String, String)>,
) -> ConsensusProfile {
    let (variants, summary) = navigator_domain::consensus::interpret(&observed, polarity);
    ConsensusProfile {
        variants,
        summary,
        terminal: observed.terminal_hint,
        sources: observed
            .sources
            .into_iter()
            .map(|s| YSourceSummary {
                label: s.label,
                source_type: s.source_type,
                variant_count: s.variant_count,
            })
            .collect(),
    }
}

/// Normalize a legacy baked [`ConsensusProfile`] payload into an
/// [`navigator_domain::consensus::ObservedProfile`], preserving each source's stored observed base (a
/// base-less legacy source becomes a no-call on interpret until the next rebuild). Load-time
/// backward-compat only.
fn observed_from_legacy(legacy: ConsensusProfile) -> navigator_domain::consensus::ObservedProfile {
    use navigator_domain::consensus::{ObservedProfile, ObservedSource, ObservedVariant, SourceSummary};
    ObservedProfile {
        schema_version: 1,
        variants: legacy
            .variants
            .into_iter()
            .map(|v| ObservedVariant {
                name: v.name,
                position: v.position,
                in_tree: v.in_tree,
                ref_allele: Some(v.ancestral),
                alt_allele: Some(v.derived),
                sources: v
                    .sources
                    .into_iter()
                    .map(|s| ObservedSource {
                        label: s.label,
                        source_type: s.source_type,
                        base: s.base,
                        depth: None,
                        mapq: None,
                        callable: None,
                        region_modifier: 1.0,
                    })
                    .collect(),
            })
            .collect(),
        sources: legacy
            .sources
            .into_iter()
            .map(|s| SourceSummary {
                label: s.label,
                source_type: s.source_type,
                variant_count: s.variant_count,
            })
            .collect(),
        terminal_hint: legacy.terminal,
    }
}

#[cfg(test)]
mod lifted_targets_tests {
    use super::*;
    use navigator_store::Store;

    /// CHM13's `chrM` is a circular permutation of rCRS with its origin near rCRS 577. Build one
    /// synthetically: `chrm[i] = rcrs[(i + K) % n]`, i.e. `rcrs[K..] ++ rcrs[..K]`.
    const K: usize = 576;

    /// Write a one-line `chrM` FASTA + its `.fai` (noodles' indexed reader needs the index).
    fn rotated_chrm_fasta(dir: &Path) -> PathBuf {
        let rcrs = navigator_analysis::mtvariants::rcrs();
        let n = rcrs.len();
        let chrm: String = format!("{}{}", &rcrs[K..], &rcrs[..K]);
        assert_eq!(chrm.len(), n);

        let fa = dir.join("rotated-chrM.fa");
        std::fs::write(&fa, format!(">chrM\n{chrm}\n")).unwrap();
        // name, length, offset-of-first-base, bases-per-line, bytes-per-line
        std::fs::write(fa.with_extension("fa.fai"), format!("chrM\t{n}\t6\t{n}\t{}\n", n + 1)).unwrap();
        fa
    }

    /// rCRS 1-based `p` sits at chrM 1-based `((p - 1 - K) mod n) + 1`.
    fn expected_chrm_pos(p: i64, n: i64) -> i64 {
        (p - 1 - K as i64).rem_euclid(n) + 1
    }

    /// The mt regression: with `tree_source_build = None`, `lifted_targets` reaches the rCRS↔chrM
    /// map and rewrites rCRS tree positions into this reference's rotated `chrM` frame. Covers the
    /// wrap (263 lands past the origin) and a plain interior site (750).
    #[tokio::test]
    async fn mt_targets_lift_onto_a_rotated_chm13_chrm() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let dir = std::env::temp_dir().join("nav-lifted-targets-mt");
        std::fs::create_dir_all(&dir).unwrap();
        let fa = rotated_chrm_fasta(&dir);
        let n = navigator_analysis::mtvariants::rcrs().len() as i64;

        let targets: HashSet<i64> = [263, 750].into_iter().collect();
        let lifted = app
            .lifted_targets("chm13v2.0", Some(&fa), "chrM", &targets, None)
            .await
            .expect("lift must succeed")
            .expect("a CHM13 chrM reference must produce a map");

        for want in [263, 750] {
            let got = lifted
                .iter()
                .find(|l| l.tree_pos == want)
                .unwrap_or_else(|| panic!("tree position {want} was not lifted"));
            assert_eq!(
                got.pos,
                expected_chrm_pos(want, n),
                "rCRS {want} must land at its rotated chrM coordinate"
            );
        }
        // 263 is before the rotation origin, so it wraps to the tail of chrM.
        assert!(lifted.iter().any(|l| l.tree_pos == 263 && l.pos > 16_000));
    }

    /// The bug this guards: `mt_tree_rcrs` returns a tree *provider*, not a build. Passing it as
    /// `tree_source_build` used to fall through to `Ok(None)`, silently skipping the chrM map above
    /// and miscalling every marker. It must be refused instead.
    #[tokio::test]
    async fn a_tree_provider_is_refused_as_a_reference_build() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let dir = std::env::temp_dir().join("nav-lifted-targets-provider");
        std::fs::create_dir_all(&dir).unwrap();
        let fa = rotated_chrm_fasta(&dir);
        let targets: HashSet<i64> = [263, 750].into_iter().collect();

        for provider in ["decodingus", "ftdna"] {
            let err = app
                .lifted_targets("chm13v2.0", Some(&fa), "chrM", &targets, Some(provider))
                .await
                .expect_err("a provider name is not a reference build");
            let msg = err.to_string();
            assert!(msg.contains(provider), "error must name the offender: {msg}");
            assert!(msg.contains("not a reference build"), "{msg}");
        }
    }

    /// A real build that matches the alignment's build needs no lift — still `Ok(None)`, not an
    /// error. Pins that the new guard didn't narrow the chrY path.
    #[tokio::test]
    async fn a_matching_reference_build_still_means_no_lift() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let targets: HashSet<i64> = [2_781_000].into_iter().collect();
        let lifted = app
            .lifted_targets("GRCh38", None, "chrY", &targets, Some("GRCh38"))
            .await
            .expect("same build is not an error");
        assert!(lifted.is_none(), "no chain needed when the builds agree");
    }
}

#[cfg(test)]
mod vset_autosomal_calls_tests {
    use super::*;
    use navigator_domain::variants::{SourceType, VariantCall, VariantSet};

    fn call(contig: &str, pos: i64, r: &str, a: &str, gt: &str) -> VariantCall {
        VariantCall {
            contig: contig.into(),
            position: pos,
            reference: r.into(),
            alternate: a.into(),
            rs_id: None,
            genotype: (!gt.is_empty()).then(|| gt.to_string()),
        }
    }

    fn set(calls: Vec<VariantCall>) -> VariantSet {
        VariantSet {
            id: 1,
            biosample_guid: du_domain::ids::SampleGuid(uuid::Uuid::nil()),
            source_label: "cg".into(),
            source_type: SourceType::WgsShortRead,
            reference_build: Some("GRCh37".into()),
            calls,
        }
    }

    #[test]
    fn genotype_becomes_reference_forward_allele_pair() {
        let s = set(vec![
            call("chr1", 100, "C", "T", "1/1"),  // hom-alt → (T, T)
            call("chr1", 200, "A", "G", "0/1"),  // het → (A, G)
            call("chr1", 300, "G", "A", "1/."),  // het w/ no-call partner → (G, A)
            call("chr7", 400, "A", "C", ""),     // no genotype → assume het → (A, C)
            call("chr2", 500, "A", "G", "1/2"),  // tri-allelic → dropped
            call("chrY", 600, "A", "G", "1"),    // not autosomal → dropped
            call("chrM", 700, "A", "G", "1"),    // not autosomal → dropped
        ]);
        let mut got = App::vset_autosomal_calls(&s);
        got.sort_by_key(|(_, p, _, _)| *p);
        assert_eq!(
            got,
            vec![
                ("1".to_string(), 100, 'T', 'T'),
                ("1".to_string(), 200, 'A', 'G'),
                ("1".to_string(), 300, 'G', 'A'),
                ("7".to_string(), 400, 'A', 'C'),
            ]
        );
    }
}
