//! `impl App` methods extracted from `lib.rs` (the `gvcf` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// Load a bundled chrY position mask/blocklist BED (best-effort). `env_var` overrides the path;
/// otherwise `<cache base>/masks/<file>`. Returns `None` if the file is absent, unparseable, or
/// empty — so a missing cohort asset simply skips that filter rather than blocking the analysis.
fn load_y_position_bed(env_var: &str, file: &str) -> Option<navigator_analysis::mask::RegionMask> {
    let path = std::env::var(env_var)
        .map(PathBuf::from)
        .unwrap_or_else(|_| refgenome_cache::base_dir().join("masks").join(file));
    navigator_analysis::mask::RegionMask::from_bed(&path, "chrY")
        .ok()
        .filter(|m| !m.is_empty())
}

impl App {
    // ---- fast path: place haplogroups from precomputed pipeline GVCFs ---------

    /// Build a tree's per-position base calls for an alignment from a **precomputed GVCF**
    /// (the fast path — no CRAM pileup). Lifts tree positions onto the GVCF's build when the
    /// tree's coordinates differ (mt rCRS-tree vs CHM13 `chrM`), exactly as the CRAM path does,
    /// then reads the GVCF instead of walking reads.
    async fn gvcf_base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        gvcf: &Path,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        // The reference is required: a GVCF hom-ref site means "the sample's base == the
        // reference base" — and the reference (e.g. CHM13 = HG002/J1 Y) is itself deep in the
        // tree, so its base there is often the *derived* allele, not the ancestral. We read the
        // reference base at every callable tree position (exactly what call_bases_at observes).
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        let targets: HashSet<i64> = tree
            .nodes
            .values()
            .flat_map(|n| n.loci.iter().map(|l| l.position))
            .collect();
        if targets.is_empty() {
            return Ok(HashMap::new());
        }
        let params = gvcf::GvcfReadParams::default();

        let lifted = self
            .lifted_targets(
                &aln.reference_build,
                Some(&reference),
                contig,
                &targets,
                tree_source_build,
            )
            .await?;

        match lifted {
            // Native: tree positions are already in the GVCF's coordinates → direct read, then
            // resolve hom-ref bases from the reference at the same positions.
            None => {
                let gvcf = gvcf.to_path_buf();
                let contig_s = contig.to_string();
                let targets2 = targets.clone();
                let called =
                    tokio::task::spawn_blocking(move || gvcf::read_called_bases(&gvcf, &contig_s, &targets2, &params))
                        .await??;
                let ref_base = self.reference_bases(&reference, contig, &called.callable).await?;
                Ok(gvcf::assemble_calls(&called, &ref_base))
            }
            // Lifted: read the GVCF at each lifted contig + the reference bases there, then map
            // observations back to tree positions (reverse-complementing minus-strand lifts).
            Some(lifted) => {
                let mut by_contig: HashMap<String, HashSet<i64>> = HashMap::new();
                for lp in &lifted {
                    by_contig.entry(lp.contig.clone()).or_default().insert(lp.pos);
                }
                let mut all = gvcf::CalledBases::default();
                let mut ref_base: HashMap<i64, char> = HashMap::new();
                for (qcontig, set) in by_contig {
                    let gvcf = gvcf.to_path_buf();
                    let qc = qcontig.clone();
                    let set2 = set.clone();
                    let called =
                        tokio::task::spawn_blocking(move || gvcf::read_called_bases(&gvcf, &qc, &set2, &params))
                            .await??;
                    ref_base.extend(self.reference_bases(&reference, &qcontig, &called.callable).await?);
                    all.variant_bases.extend(called.variant_bases);
                    all.callable.extend(called.callable);
                }
                Ok(assemble_calls_lifted(&all, &lifted, &ref_base))
            }
        }
    }

    /// Reference genome bases (uppercase A/C/G/T) at `positions` on `contig`. Reads the contig
    /// sequence once off-thread; positions are 1-based. Non-ACGT / out-of-range positions are
    /// omitted. Used by the GVCF fast path to resolve hom-ref tree sites to the actual base.
    async fn reference_bases(
        &self,
        reference: &Path,
        contig: &str,
        positions: &HashSet<i64>,
    ) -> Result<HashMap<i64, char>, AppError> {
        if positions.is_empty() {
            return Ok(HashMap::new());
        }
        let reference = reference.to_path_buf();
        let contig = contig.to_string();
        let positions: Vec<i64> = positions.iter().copied().collect();
        let map = tokio::task::spawn_blocking(
            move || -> Result<HashMap<i64, char>, navigator_analysis::AnalysisError> {
                let seq = navigator_analysis::reader::read_contig_sequence(&reference, &contig)?;
                let mut m = HashMap::with_capacity(positions.len());
                for p in positions {
                    if p >= 1 && (p as usize) <= seq.len() {
                        let b = seq[p as usize - 1].to_ascii_uppercase();
                        if matches!(b, b'A' | b'C' | b'G' | b'T') {
                            m.insert(p, b as char);
                        }
                    }
                }
                Ok(m)
            },
        )
        .await??;
        Ok(map)
    }

    /// Fingerprint of a GVCF-sourced placement: the GVCF's content hash ⊕ the tree's hash.
    /// Distinct from the CRAM-based [`Self::y_score_fingerprint`] (`gv:` vs `f:` prefix) so a
    /// later deep analyze can tell the call came from a sidecar (phase: deep-pass skip logic).
    async fn gvcf_fingerprint(&self, gvcf: &Path, tree_json: &str, tag: &str) -> Result<String, AppError> {
        let h = sha256_file_async(gvcf.to_path_buf()).await?;
        Ok(format!("gv:{}|{}:{}", &h[..16], tag, &sha256_str(tree_json)[..16]))
    }

    /// Assign a Y haplogroup from a precomputed chrY GVCF — no CRAM walk. Places against the
    /// DecodingUs tree at the alignment's native build (liftover-free), records the call under
    /// the same source key as the CRAM path (`aln:{id}`) with a `gv:`-prefixed fingerprint.
    /// Errors if the build has no DecodingUs coordinates or the tree is unreachable; the caller
    /// (`ingest_sidecars`) treats that as "leave Y for the deep pass".
    pub async fn assign_y_from_gvcf(&self, alignment_id: i64, gvcf: &Path) -> Result<HaploAssignment, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let build_key = decodingus_build_key(&aln.reference_build).ok_or_else(|| {
            AppError::Import(format!(
                "no DecodingUs tree coordinates for build {}",
                aln.reference_build
            ))
        })?;
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_decodingus_json(&tree_json, build_key).map_err(AppError::Import)?;
        let calls = self.gvcf_base_calls(alignment_id, "chrY", gvcf, &tree, None).await?;
        // Robust (proportional-top) selection, not the strict alignment-tuned guard. A
        // joint-genotyped GVCF gives confident calls that include a few stray ancestral
        // contradictions on the deep backbone (recurrent sites, the CHM13=J1 reference, joint
        // hard-filters); strict `path_admissible` then vetoes the genuine deep lineage and
        // drops to a shallow node (HG00096 → A1b instead of its true R1b1a1b1a1a, which `score`
        // ranks top at 344/364). This is the same confident-but-sparse-contradiction regime as
        // BISDNA chip data — see [`assemble_assignment_robust`].
        let assignment = assemble_assignment_robust(&tree, &calls);
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            let fp = self.gvcf_fingerprint(gvcf, &tree_json, "yt").await.ok();
            self.record_call_fp(
                bio,
                DnaType::Y,
                &format!("aln:{alignment_id}"),
                format!("aln #{alignment_id} Y (pipeline GVCF)"),
                &assignment,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Assign an mtDNA haplogroup from a precomputed chrM GVCF — no CRAM walk. Places against
    /// the FTDNA mt tree; on CHM13 the tree's rCRS positions are lifted onto `chrM` (the cheap
    /// self-generated rCRS↔chrM map), on GRCh38 they're read directly. Recorded under the CRAM
    /// path's mt source key (`aln:{id}:mt`) with a `gv:`-prefixed fingerprint.
    pub async fn assign_mt_from_gvcf(&self, alignment_id: i64, gvcf: &Path) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;
        let source_build = tree_build_for_contig("chrM"); // None → rCRS-direct / chrM lift
        let calls = self
            .gvcf_base_calls(alignment_id, "chrM", gvcf, &tree, source_build)
            .await?;
        // Robust selection, as for Y (see assign_y_from_gvcf) — the GVCF's confident calls fit
        // the proportional-top regime better than the strict alignment guard.
        let assignment = assemble_assignment_robust(&tree, &calls);
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            let fp = self.gvcf_fingerprint(gvcf, &tree_json, "mt").await.ok();
            self.record_call_fp(
                bio,
                DnaType::Mt,
                &format!("aln:{alignment_id}:mt"),
                format!("aln #{alignment_id} mtDNA (pipeline GVCF)"),
                &assignment,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Fast-path ingest of a sample's pipeline sidecars onto one alignment: place Y + mt from
    /// the GVCFs, and fill sex / read-metrics / lite-coverage from the text sidecars — all
    /// without touching the CRAM. Each step is independent and best-effort: a failure is
    /// recorded in the returned report and the rest proceed (a missing/!matching sidecar just
    /// leaves that result for the deep pass). Returns what it managed to fill.
    pub async fn ingest_sidecars(
        &self,
        alignment_id: i64,
        sidecars: &SampleSidecars,
    ) -> Result<SidecarIngest, AppError> {
        let mut out = SidecarIngest::default();

        if let Some(gvcf) = &sidecars.chr_y_gvcf {
            match self.assign_y_from_gvcf(alignment_id, gvcf).await {
                Ok(a) => out.y_haplogroup = a.ranked.first().map(|r| r.name.clone()),
                Err(e) => out.errors.push(format!("Y from GVCF: {e}")),
            }
        }
        if let Some(gvcf) = &sidecars.chr_m_gvcf {
            match self.assign_mt_from_gvcf(alignment_id, gvcf).await {
                Ok(a) => out.mt_haplogroup = a.ranked.first().map(|r| r.name.clone()),
                Err(e) => out.errors.push(format!("mt from GVCF: {e}")),
            }
        }
        if let Some(path) = &sidecars.sex {
            match self.ingest_sex_sidecar(alignment_id, path).await {
                Ok(Some(s)) => out.sex = Some(s),
                Ok(None) => {} // kept an existing full result (reimport) — not overwritten
                Err(e) => out.errors.push(format!("sex: {e}")),
            }
        }
        // Read metrics: richest source wins — samtools `stats` (full, with histograms) > Picard
        // AlignmentSummaryMetrics > samtools `flagstat` (counts only).
        match self.ingest_read_metrics(alignment_id, sidecars).await {
            Ok(true) => out.read_metrics = true,
            Ok(false) => {}
            Err(e) => out.errors.push(format!("read metrics: {e}")),
        }
        // Coverage: samtools `coverage` gives per-contig stats; Picard CollectWgsMetrics gives the
        // genome-wide depth distribution (median/sd/MAD, exclusion fractions, pct_Nx). Use whichever
        // are present, overlaying the distribution onto the per-contig breakdown.
        if sidecars.coverage.is_some() || sidecars.wgs_metrics.is_some() {
            match self.ingest_coverage_sidecar(alignment_id, sidecars).await {
                Ok(wrote) => out.lite_coverage = wrote,
                Err(e) => out.errors.push(format!("coverage: {e}")),
            }
        }
        Ok(out)
    }

    /// Ingest inferred sex from the sidecar. `Ok(None)` when an equal-or-fuller result is already
    /// stored (reimport) so we neither overwrite the artifact nor re-stamp the sequence run.
    async fn ingest_sex_sidecar(&self, alignment_id: i64, path: &Path) -> Result<Option<String>, AppError> {
        let text = tokio::fs::read_to_string(path)
            .await
            .map_err(|e| AppError::Import(format!("{}: {e}", path.display())))?;
        let result = sidecar::parse_sex(&text);
        let wrote = self
            .save_analysis_no_downgrade(alignment_id, "sex", "1", &result, "pipeline-sidecar", "full")
            .await?;
        if !wrote {
            return Ok(None);
        }
        self.write_back_inferred_sex(alignment_id, &result).await?;
        Ok(Some(
            match result.inferred_sex {
                InferredSex::Male => "M",
                InferredSex::Female => "F",
                InferredSex::Unknown => "U",
            }
            .to_string(),
        ))
    }

    /// Ingest read metrics from the best available sidecar (priority: samtools `stats` →
    /// Picard AlignmentSummaryMetrics → samtools `flagstat`). Returns whether one was found.
    async fn ingest_read_metrics(&self, alignment_id: i64, sidecars: &SampleSidecars) -> Result<bool, AppError> {
        let read = |p: &Path| {
            let p = p.to_path_buf();
            async move {
                tokio::fs::read_to_string(&p)
                    .await
                    .map_err(|e| AppError::Import(format!("{}: {e}", p.display())))
            }
        };
        // (metrics, completeness): samtools stats is full (carries histograms); the others are
        // counts/scalars only, so `partial` lets a deep read-metrics walk upgrade them later.
        let (metrics, completeness) = if let Some(p) = &sidecars.stats {
            (sidecar::parse_samtools_stats(&read(p).await?), "full")
        } else if let Some(p) = &sidecars.alignment_summary {
            match sidecar::parse_alignment_summary(&read(p).await?) {
                Some(m) => (m, "partial"),
                None => return Ok(false),
            }
        } else if let Some(p) = &sidecars.flagstat {
            (sidecar::parse_flagstat(&read(p).await?), "partial")
        } else {
            return Ok(false);
        };
        // Don't downgrade a full deep walk on reimport — keep it if it's already equal-or-fuller.
        let wrote = self
            .save_analysis_no_downgrade(alignment_id, "read_metrics", "1", &metrics, "pipeline-sidecar", completeness)
            .await?;
        Ok(wrote)
    }

    /// Ingest lite coverage from the sidecar(s). Returns whether it was written (`false` = an
    /// equal-or-fuller coverage artifact already exists, e.g. a deep walk on reimport).
    async fn ingest_coverage_sidecar(&self, alignment_id: i64, sidecars: &SampleSidecars) -> Result<bool, AppError> {
        let read = |p: &Path| {
            let p = p.to_path_buf();
            async move {
                tokio::fs::read_to_string(&p)
                    .await
                    .map_err(|e| AppError::Import(format!("{}: {e}", p.display())))
            }
        };
        // Per-contig stats + callable counts from samtools coverage (empty base if absent).
        let lite = match &sidecars.coverage {
            Some(cp) => {
                let cov = read(cp).await?;
                let summary = match &sidecars.callable_summary {
                    Some(p) => Some(read(p).await?),
                    None => None,
                };
                sidecar::lite_coverage(&cov, summary.as_deref())
            }
            None => CoverageResult::default(),
        };
        // Overlay Picard's genome-wide depth distribution onto the per-contig breakdown: start from
        // the Picard result (median/sd/MAD, exclusion fractions, pct_Nx) and graft the contig stats.
        let result = match &sidecars.wgs_metrics {
            Some(wp) => match sidecar::parse_wgs_metrics(&read(wp).await?) {
                Some(mut w) => {
                    w.contig_coverage_stats = lite.contig_coverage_stats;
                    w.contig_callable = lite.contig_callable;
                    w.callable_bases = lite.callable_bases;
                    if w.genome_territory == 0 {
                        w.genome_territory = lite.genome_territory;
                    }
                    if w.mean_coverage == 0.0 {
                        w.mean_coverage = lite.mean_coverage;
                    }
                    w
                }
                None => lite,
            },
            None => lite,
        };
        // Still `partial`: no per-base depth histogram (only the deep walk produces that), so the
        // deep pass still upgrades this. Stored under the standard coverage key. Never downgrade a
        // full deep-walk coverage on reimport — keep it if one is already present.
        let wrote = self
            .save_analysis_no_downgrade(
                alignment_id,
                "coverage",
                coverage::COVERAGE_VERSION,
                &result,
                "pipeline-sidecar",
                "partial",
            )
            .await?;
        Ok(wrote)
    }

    /// Self-referential callable intervals (BED 0-based half-open) for `contig` from the
    /// alignment's own reads. Parameters adapt to the sample: long reads (HiFi) earn
    /// callability at lower depth, and the CALLABLE-run gate scales with molecule length
    /// (`f`·fragment), so long molecules clear it over far more of chrY. Requires the BAM.
    pub async fn callable_chr_intervals(&self, alignment_id: i64, contig: &str) -> Result<Vec<(i64, i64)>, AppError> {
        // Resolve the reference via the gateway when the alignment has no stored path — a CRAM can't
        // be decoded without one, and most imported alignments leave `reference_path` null (the build
        // alone is recorded). Same resolution the de-novo caller uses.
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let reference = Some(reference);
        let contig = contig.to_string();
        tokio::task::spawn_blocking(move || {
            let (read_len, frag_len) = coverage::estimate_molecule_lengths(&bam, reference.as_deref())?;
            let molecule = frag_len.max(read_len);
            let mut params = CallableLociParams::default();
            // Long, accurate reads (HiFi) are callable from a single read (see adaptive_min_depth).
            params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            let min_run_len = molecule.round().max(1.0) as u32; // f = 1.0
            coverage::callable_intervals(&bam, &contig, &params, min_run_len, reference.as_deref())
        })
        .await?
        .map_err(Into::into)
    }

    /// The **private bucket**: de-novo SNP calls on chrY that the Y placement doesn't
    /// explain (not on the assigned backbone), classified as off-path-known (a finer/
    /// sibling FTDNA branch) or novel (a new-branch candidate). With `callable_bed` (e.g.
    /// the Poznik/1KG `b38_sites.bed`), calls outside reliable regions are dropped.
    pub async fn private_y_variants(
        &self,
        alignment_id: i64,
        callable_bed: Option<&Path>,
    ) -> Result<PrivateBucket, AppError> {
        let mask = match callable_bed {
            Some(p) => Some(navigator_analysis::mask::RegionMask::from_bed(p, "chrY")?),
            None => None,
        };
        self.private_y_core(alignment_id, mask).await
    }

    /// [`private_y_variants`] using the sample's **own** callable-Y BED as the mask
    /// (self-referential — adapts to the sample's depth and read tech; no external file).
    pub async fn private_y_variants_self_masked(&self, alignment_id: i64) -> Result<PrivateBucket, AppError> {
        let intervals = self.callable_chr_intervals(alignment_id, "chrY").await?;
        let mask = navigator_analysis::mask::RegionMask::from_intervals(intervals);
        let bucket = self.private_y_core(alignment_id, Some(mask)).await?;
        // Persist the self-masked bucket so it reloads instead of recomputing next session. Version
        // "2": adds `alt_depth` + the cohort callable-mask / recurrent-exclude filters, so v1 blobs
        // (unfiltered, no alt_depth) must recompute rather than reload.
        self.save_analysis(alignment_id, "private_y", "2", &bucket).await?;
        Ok(bucket)
    }

    /// Cached self-masked private-Y bucket for an alignment, if previously computed.
    pub async fn cached_private_y(&self, alignment_id: i64) -> Result<Option<PrivateBucket>, AppError> {
        self.load_analysis(alignment_id, "private_y", "2").await
    }

    /// Shared core: assign Y, de-novo chrY, subtract the backbone, optionally mask, classify.
    /// The curated CHM13 chrY structural regions (palindrome/amplicon/AZF-DYZ), resolving +
    /// caching the three BEDs on first use. Best-effort: any download/parse failure yields
    /// `None` so the annotation never blocks the analysis.
    /// Genome-region metadata (centromere/telomere/cytoband/PAR) for a build, via the gateway's
    /// 2-layer cache (fetches the UCSC cytoBand table on a cold miss). For QC / display context.
    pub async fn genome_regions(&self, build: &str) -> Result<std::sync::Arc<GenomeRegions>, AppError> {
        Ok(self.gateway.genome_regions(build, &mut |_, _| {}).await?)
    }

    /// Region annotation for a 1-based `position` on `contig` in `build` (centromere/telomere/PAR
    /// membership + cytoband name). Uses the cached regions only — `None` if not yet fetched.
    pub fn region_annotation(&self, build: &str, contig: &str, position: i64) -> Option<RegionAnnotation> {
        self.gateway
            .cached_genome_regions(build)
            .map(|r| r.annotate(contig, position))
    }

    async fn y_structural_regions(&self) -> Option<navigator_analysis::mask::YStructuralRegions> {
        let amplicon = self
            .gateway
            .resolve_mask("chm13v2.0Y_amplicons_v1", &mut |_, _| {})
            .await
            .ok()?;
        let palindrome = self
            .gateway
            .resolve_mask("chm13v2.0Y_inverted_repeats_v1", &mut |_, _| {})
            .await
            .ok()?;
        let azf_dyz = self
            .gateway
            .resolve_mask("chm13v2.0Y_AZF_DYZ_v1", &mut |_, _| {})
            .await
            .ok()?;
        navigator_analysis::mask::YStructuralRegions::from_beds(&amplicon, &palindrome, &azf_dyz).ok()
    }

    async fn private_y_core(
        &self,
        alignment_id: i64,
        mask: Option<navigator_analysis::mask::RegionMask>,
    ) -> Result<PrivateBucket, AppError> {
        // Classify novels against the **DecodingUs** tree — the app's placement authority, which
        // folds in the cohort-derived branches (from the de-novo tree pipeline). A shared lineage
        // variant is named there, so it reads as OffPathKnown, not a false "novel"; a variant absent
        // from this tree yet shared across the cohort is genuinely suspect. FTDNA fallback keeps the
        // report working when the AppView tree is unavailable or the build has no DecodingUs coords.
        let (tree, tree_calls) = match self.y_decodingus_tree_calls(alignment_id).await {
            Ok(tc) => tc,
            Err(e) => {
                eprintln!("DecodingUs Y tree unavailable ({e}); private-Y classifying against FTDNA");
                let tree_json = self.fetch_ftdna_y_tree().await?;
                self.tree_base_calls(alignment_id, "chrY", &tree_json).await?
            }
        };
        let assignment = assemble_assignment(&tree, &tree_calls);
        let terminal = assignment
            .ranked
            .first()
            .ok_or_else(|| AppError::Import("no Y haplogroup match".into()))?;
        let path = navigator_analysis::haplo::path_positions(&tree, terminal.id);
        let known = navigator_analysis::haplo::tree_positions(&tree);

        // The structural BEDs + cohort masks are in CHM13 chrY coordinates, so they only apply to a
        // CHM13 alignment (the de-novo positions are in the alignment's build). Best-effort.
        let aln = self.alignment_or_err(alignment_id).await?;
        let is_chm13 = matches!(
            canonical_build(&aln.reference_build),
            Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs)
        );
        let regions = if is_chm13 { self.y_structural_regions().await } else { None };
        // L3: the cohort **callable mask** (Poznik-style, CALLABLE in ≥90% of a ~3k-male cohort) —
        // only ~25% of non-PAR chrY is reliably callable cohort-wide. L4: a **cohort-shared-sites**
        // blocklist — every position that varies with ≥2 carriers across the cohort (plus homoplasy
        // hotspots). A real shared lineage variant belongs in the DecodingUs tree (and so classifies
        // as off-path-known above); one that is cohort-shared yet *absent* from the tree is a suspect
        // recurrent artifact, not a private SNP. A truly private variant has a single cohort carrier,
        // so it survives this filter. This is the single-sample stand-in for the de-novo pipeline's
        // cohort carrier filter. Both bundled, CHM13-only; absent ⇒ that filter is skipped.
        let cohort_mask =
            if is_chm13 { load_y_position_bed("NAVIGATOR_Y_CALLABLE_MASK", "chrY_callable_mask.chm13v2.bed") } else { None };
        let cohort_shared = if is_chm13 {
            load_y_position_bed("NAVIGATOR_Y_COHORT_SHARED", "chrY_cohort_shared_sites.chm13v2.bed")
        } else {
            None
        };

        // De-novo chrY (cached as an artifact), then keep only off-backbone, callable, non-shared calls.
        let denovo = self.run_denovo_for_alignment(alignment_id, "chrY".to_string()).await?;
        let mut variants: Vec<PrivateVariant> = denovo
            .iter()
            .filter(|c| !path.contains(&c.position))
            .filter(|c| mask.as_ref().map_or(true, |m| m.contains(c.position))) // self-callable
            .filter(|c| cohort_mask.as_ref().map_or(true, |m| m.contains(c.position))) // L3 cohort callable mask
            .filter(|c| cohort_shared.as_ref().map_or(true, |m| !m.contains(c.position))) // L4 cohort-shared exclude
            .map(|c| PrivateVariant {
                position: c.position,
                reference: c.reference_allele,
                alternate: c.alternate_allele,
                depth: c.depth,
                alt_depth: c.alt_depth,
                allele_fraction: c.allele_fraction,
                class: match known.get(&c.position) {
                    Some(name) => PrivateClass::OffPathKnown(name.clone()),
                    None => PrivateClass::Novel,
                },
                region: regions.as_ref().and_then(|r| r.classify(c.position)),
            })
            .collect();
        variants.sort_by_key(|v| v.position);
        Ok(PrivateBucket {
            terminal: terminal.name.clone(),
            variants,
        })
    }
}
