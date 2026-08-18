//! `impl App` methods extracted from `lib.rs` (the `gvcf` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// Load a bundled chrY position mask/blocklist BED for a build (best-effort). `env_var` overrides the
/// path; otherwise the seeded `<cache base>/masks/<stem>.<build_token>.bed`, trying the gzipped
/// `.bed.gz` first (how the bundled assets ship) then a plain `.bed`. Returns `None` if absent,
/// unparseable, or empty — so a missing cohort asset simply skips that filter rather than blocking.
fn load_y_position_bed(env_var: &str, stem: &str, build_token: &str) -> Option<navigator_analysis::mask::RegionMask> {
    let candidates: Vec<PathBuf> = if let Ok(p) = std::env::var(env_var) {
        vec![PathBuf::from(p)]
    } else {
        let dir = refgenome_cache::base_dir().join("masks");
        let file = format!("{stem}.{build_token}.bed");
        vec![dir.join(format!("{file}.gz")), dir.join(file)]
    };
    candidates.into_iter().find_map(|path| {
        navigator_analysis::mask::RegionMask::from_bed(&path, "chrY")
            .ok()
            .filter(|m| !m.is_empty())
    })
}

/// Caller output directories a per-sample GVCF is commonly filed under, beside the alignment.
/// A pipeline that runs several callers keeps each one's output in its own directory rather than
/// beside the CRAM, so looking only at the alignment's own directory misses them.
const CALLER_SUBDIRS: [&str; 3] = ["gatk4", "gatk3", "gvcf"];

/// Locate a per-sample GVCF for `contig` beside an alignment: the alignment's own directory first
/// (a `*.chrY.g.vcf.gz` sidecar, the ytree layout), then the known caller subdirectories, where a
/// bare `chrY.g.vcf.gz` is the usual name.
///
/// Both spellings matter. The ytree flat layout emits `<sample>.chrY.g.vcf.gz` next to the CRAM; a
/// per-run pipeline emits `gatk4/chrY.g.vcf.gz`, whose name has no sample prefix at all. Matching
/// only the dotted suffix missed every file of the second kind — and since finding the GVCF is what
/// lets placement skip decoding the CRAM, missing it silently turns a seconds-long read into a
/// minutes-long whole-chromosome walk.
fn gvcf_beside_alignment(aln: &Alignment, contig_token: &str) -> Option<PathBuf> {
    let dotted = format!(".{contig_token}.g.vcf.gz");
    let bare = format!("{contig_token}.g.vcf.gz");
    let matches = |name: &str| {
        let n = name.to_ascii_lowercase();
        n.ends_with(&dotted) || n == bare
    };
    let dir = Path::new(aln.bam_path.as_ref()?).parent()?;
    let scan = |d: &Path| -> Option<PathBuf> {
        std::fs::read_dir(d)
            .ok()?
            .flatten()
            .find_map(|e| matches(&e.file_name().to_string_lossy()).then(|| e.path()))
    };
    scan(dir).or_else(|| CALLER_SUBDIRS.iter().find_map(|sub| scan(&dir.join(sub))))
}

/// Locate a per-sample chrY GVCF for an alignment: the `NAVIGATOR_Y_GVCF` path override, else
/// [`gvcf_beside_alignment`]. `None` when absent — the private-Y path then falls back to the pileup
/// caller, and placement to a full CRAM walk.
pub(crate) fn chr_y_gvcf_for_alignment(aln: &Alignment) -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_Y_GVCF") {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Some(p);
        }
    }
    gvcf_beside_alignment(aln, "chry")
}

/// Locate a per-sample chrM GVCF for an alignment: the `NAVIGATOR_M_GVCF` path override, else
/// [`gvcf_beside_alignment`]. The mtDNA counterpart to [`chr_y_gvcf_for_alignment`].
pub(crate) fn chr_m_gvcf_for_alignment(aln: &Alignment) -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_M_GVCF") {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Some(p);
        }
    }
    gvcf_beside_alignment(aln, "chrm")
}

/// The bundled-mask filename token for an alignment's reference build, or `None` when no chrY masks
/// ship for it. CHM13 masks are native (hs1); the GRCh38 masks are lifted from them (CrossMap
/// hs1→hg38). GRCh37 has no masks yet (bare-`Y` contig naming + no lifted set).
fn y_mask_build_token(build: &str) -> Option<&'static str> {
    match canonical_build(build) {
        Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs) => Some("chm13v2"),
        Some(ReferenceBuild::Grch38) => Some("grch38"),
        _ => None,
    }
}

/// A shared, process-wide copy of one build's chrY structural regions (see the memo in
/// [`App::y_structural_regions_for`]).
type YRegionsHandle = std::sync::Arc<navigator_analysis::mask::YStructuralRegions>;

impl App {
    // ---- fast path: place haplogroups from precomputed pipeline GVCFs ---------

    /// Build a tree's per-position base calls for an alignment from a **precomputed GVCF**
    /// (the fast path — no CRAM pileup). Lifts tree positions onto the GVCF's build when the
    /// tree's coordinates differ (mt rCRS-tree vs CHM13 `chrM`), exactly as the CRAM path does,
    /// then reads the GVCF instead of walking reads.
    pub(crate) async fn gvcf_base_calls(
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
                &external_y_source_key(alignment_id),
                format!("aln #{alignment_id} Y (pipeline GVCF)"),
                &assignment,
                CallProvenance::External,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Assign an mtDNA haplogroup from a precomputed chrM GVCF — no CRAM walk. Places against
    /// the FTDNA mt tree; on CHM13 the tree's rCRS positions are lifted onto `chrM` (the cheap
    /// self-generated rCRS↔chrM map), on GRCh38 they are read directly. Recorded under the CRAM
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
                &external_mt_source_key(alignment_id),
                format!("aln #{alignment_id} mtDNA (pipeline GVCF)"),
                &assignment,
                CallProvenance::External,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// The sidecar paths this alignment was ingested from, as recorded by [`Self::ingest_sidecars`].
    ///
    /// Read directly, with no source-mtime freshness check: this is a record of *what was used*,
    /// not a derived result that a changed CRAM invalidates. `None` for an alignment that never
    /// went through the fast path (imported before this was recorded, or with no sidecars at all).
    pub async fn recorded_sidecars(&self, alignment_id: i64) -> Result<Option<SampleSidecars>, AppError> {
        match artifact::get(self.store.pool(), alignment_id, SIDECARS_KIND, SIDECARS_VERSION).await? {
            Some(a) => Ok(serde_json::from_str(&a.payload).ok()),
            None => Ok(None),
        }
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

        // Record which files this alignment was ingested from, before using them. Discovery is a
        // directory scan done once at import, so without this the fast path is a one-shot: a Y
        // placement made from a GVCF against the tree of the day could never be re-derived, and the
        // resulting `haplogroup_call` row outlived every tree it was placed against. See
        // `App::replace_against_current_tree`, which replays this. Best-effort — a workspace that
        // can not record the paths should still get the ingest.
        let _ = self
            .save_analysis_with_provenance(
                alignment_id,
                SIDECARS_KIND,
                SIDECARS_VERSION,
                sidecars,
                "sidecar",
                "full",
            )
            .await;

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
        // Do not downgrade a full deep walk on reimport — keep it if it is already equal-or-fuller.
        let wrote = self
            .save_analysis_no_downgrade(
                alignment_id,
                "read_metrics",
                "1",
                &metrics,
                "pipeline-sidecar",
                completeness,
            )
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

    /// The **private bucket**: de-novo SNP calls on chrY that the Y placement does not
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
    ///
    /// With a per-sample GVCF sidecar the self-mask is **skipped**: the GVCF's own confidence
    /// gating is the callable evidence (re-imposing Navigator's callable-loci depth threshold would
    /// discard GATK calls the whole point is to trust), and skipping it avoids a CRAM walk — so the
    /// GVCF fast path stays fast. Reliability then comes from the cohort callable mask + GVCF GQ.
    pub async fn private_y_variants_self_masked(&self, alignment_id: i64) -> Result<PrivateBucket, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let mask = if chr_y_gvcf_for_alignment(&aln).is_some() {
            None
        } else {
            let intervals = self.callable_chr_intervals(alignment_id, "chrY").await?;
            Some(navigator_analysis::mask::RegionMask::from_intervals(intervals))
        };
        let bucket = self.private_y_core(alignment_id, mask).await?;
        // Persist the self-masked bucket so it reloads instead of recomputing next session. Version
        // "3": prefers a per-sample GVCF sidecar as the derived-call source (was pileup-only in v2),
        // so v2 blobs must recompute rather than reload.
        // Version 4: private variants are now classified against structural masks lifted to the
        // alignment's own build. A v3 bucket on a GRCh38 alignment saw no mask at all.
        self.save_analysis(alignment_id, "private_y", "4", &bucket).await?;
        Ok(bucket)
    }

    /// Cached self-masked private-Y bucket for an alignment, if previously computed.
    pub async fn cached_private_y(&self, alignment_id: i64) -> Result<Option<PrivateBucket>, AppError> {
        self.load_analysis(alignment_id, "private_y", "4").await
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

    /// Memo for [`y_structural_regions_for`]: lifting parses the whole chain file, and a project
    /// pass over thousands of subjects would otherwise repeat that per subject — the same trap the
    /// tree fetch fell into. The masks are static within a process, so resolve each build once.
    fn y_regions_memo() -> &'static std::sync::Mutex<HashMap<String, Option<YRegionsHandle>>> {
        static MEMO: std::sync::OnceLock<std::sync::Mutex<HashMap<String, Option<YRegionsHandle>>>> =
            std::sync::OnceLock::new();
        MEMO.get_or_init(|| std::sync::Mutex::new(HashMap::new()))
    }

    /// The curated chrY structural regions **in `build`'s coordinates**.
    ///
    /// The three BEDs are CHM13-native, so anything else is lifted. That matters more than it
    /// sounds: without it, a GRCh38 or GRCh37 source has no structural mask at all, every
    /// palindromic and amplicon call counts as unique sequence, and private-variant counts inflate
    /// into the hundreds — the difference between a donor averaging 4 and one averaging 661.
    ///
    /// Best-effort throughout: any download / chain / parse failure yields `None` so the annotation
    /// never blocks the analysis, exactly as before.
    async fn y_structural_regions_for(&self, build: &str) -> Option<YRegionsHandle> {
        // Keyed by the *canonical* build: `hs1`, `CHM13v2.0` and the masked variant share coordinates
        // and must share one entry rather than lifting three times.
        let key = canonical_build(build)?.as_str().to_string();
        if let Some(hit) = Self::y_regions_memo().lock().unwrap().get(&key) {
            return hit.clone();
        }
        let built = self.build_y_structural_regions(build).await.map(std::sync::Arc::new);
        if built.is_none() {
            // Cached so a batch does not retry a failing download per subject — but *said*, because
            // "no structural mask" is the condition that inflated private-variant counts into the
            // hundreds in the first place, and it must never be reached silently again.
            eprintln!(
                "no chrY structural mask available for {key} — private-variant counts will include \
                 palindromic and amplicon calls"
            );
        }
        Self::y_regions_memo().lock().unwrap().insert(key, built.clone());
        built
    }

    async fn build_y_structural_regions(&self, build: &str) -> Option<navigator_analysis::mask::YStructuralRegions> {
        use navigator_analysis::mask::{RegionMask, YStructuralRegions};

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
        let native = YStructuralRegions::from_beds(&amplicon, &palindrome, &azf_dyz).ok()?;

        let target = canonical_build(build)?;
        if matches!(target, ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs) {
            return Some(native);
        }
        // PAR and heterochromatin are taken natively per build rather than lifted — a chain is least
        // trustworthy in exactly those places (PAR is shared with chrX, Yq12 is satellite).
        let landmarks = navigator_analysis::mask::y_landmarks(build)?;

        self.gateway
            .resolve_chain(ReferenceBuild::Chm13v2.as_str(), target.as_str(), &mut |_, _| {})
            .await
            .ok()?;
        let lift = |m: &RegionMask, what: &str| -> Option<RegionMask> {
            let (iv, dropped) = self
                .gateway
                .lift_intervals(ReferenceBuild::Chm13v2.as_str(), target.as_str(), "chrY", m.intervals())
                .ok()?;
            if iv.is_empty() {
                // A mask that lifted to nothing is not a mask; better to annotate nothing than to
                // report "no structural regions here" as though it had been checked.
                eprintln!("chrY {what} mask lifted CHM13→{} to nothing; skipping", target.as_str());
                return None;
            }
            if dropped > 0 {
                eprintln!(
                    "chrY {what} mask CHM13→{}: {} of {} intervals dropped as unliftable",
                    target.as_str(),
                    dropped,
                    m.intervals().len()
                );
            }
            Some(RegionMask::from_intervals(iv))
        };
        let (palindrome_m, amplicon_m) = native.structural_masks();
        Some(YStructuralRegions::from_masks(
            RegionMask::from_intervals(landmarks.par.to_vec()),
            lift(palindrome_m, "palindrome")?,
            lift(amplicon_m, "amplicon")?,
            // The satellite arrays rarely survive a chain, so the build's heterochromatin bound is
            // the load-bearing part here and the lifted AZF/DYZ intervals only refine it.
            lift(native.heterochromatin_mask(), "AZF/DYZ")
                .unwrap_or_else(|| RegionMask::from_intervals(vec![]))
                .union(&[landmarks.heterochromatin]),
        ))
    }

    /// Derive chrY private-variant candidates from a per-sample GVCF, returning the same
    /// [`VariantCall`] shape the pileup de-novo path produces so `private_y_core`'s downstream
    /// classification is identical. GATK's reassembly recovers SNVs the pileup caller misses.
    async fn run_denovo_from_gvcf(&self, gvcf: &Path) -> Result<Vec<VariantCall>, AppError> {
        let gvcf = gvcf.to_path_buf();
        let snvs = tokio::task::spawn_blocking(move || {
            // min_dp 4 (over the reader's permissive default 2): a real private SNV is covered by ≥4
            // reads, whereas a misaligned-read cluster smears 2–3 reads across many nearby false SNVs
            // — so the depth floor removes those artifact clusters without touching the (DP≥4) truth.
            let params = navigator_analysis::gvcf::GvcfReadParams { min_dp: 4, min_gq: 20 };
            navigator_analysis::gvcf::read_derived_snvs(&gvcf, "chrY", &params)
        })
        .await??;
        Ok(snvs
            .into_iter()
            .map(|s| VariantCall {
                contig: "chrY".to_string(),
                position: s.position,
                reference_allele: s.reference,
                alternate_allele: s.alternate,
                depth: s.depth,
                alt_depth: s.alt_depth,
                allele_fraction: s.allele_fraction,
                quality: None,
            })
            .collect())
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
            // A gone alignment file is not a tree problem: the fallback reads the same absent file,
            // so it can only fail again while logging a tree provider that was never at fault.
            Err(e) if e.is_missing_alignment_file() => return Err(e),
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

        // The structural BEDs are in CHM13 chrY coordinates, so they only annotate a CHM13 alignment.
        // The cohort masks apply per build: native for CHM13, CrossMap-lifted (hs1→hg38) for GRCh38.
        let aln = self.alignment_or_err(alignment_id).await?;
        let regions = self.y_structural_regions_for(&aln.reference_build).await;
        // L2: the cohort **callable mask** (Poznik-style, CALLABLE in ≥90% of a ~3k-male cohort) —
        // only ~25% of non-PAR chrY is reliably callable cohort-wide. L3: a **cohort-shared-sites**
        // blocklist — every position that varies with ≥2 carriers across the cohort (plus homoplasy
        // hotspots). A real shared lineage variant belongs in the DecodingUs tree (and so classifies
        // as off-path-known above); one that is cohort-shared yet *absent* from the tree is a suspect
        // recurrent artifact, not a private SNP. A truly private variant has a single cohort carrier,
        // so it survives this filter. This is the single-sample stand-in for the de-novo pipeline's
        // cohort carrier filter. Bundled per build (CHM13 native, GRCh38 lifted); absent ⇒ skipped.
        let mask_token = y_mask_build_token(&aln.reference_build);
        let cohort_mask =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_CALLABLE_MASK", "chrY_callable_mask", t));
        let cohort_shared =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_COHORT_SHARED", "chrY_cohort_shared_sites", t));

        // Derived-call source. Prefer a per-sample chrY GVCF (GATK HaplotypeCaller's local haplotype
        // reassembly resolves misaligned-ref ~50/50 sites the pileup caller drops — see the WGS229
        // recall gap); fall back to Navigator's de-novo pileup caller when no sidecar is present.
        let denovo = match chr_y_gvcf_for_alignment(&aln) {
            Some(gvcf) => {
                eprintln!("private-Y: sourcing chrY calls from GVCF sidecar {}", gvcf.display());
                self.run_denovo_from_gvcf(&gvcf).await?
            }
            None => self.run_denovo_for_alignment(alignment_id, "chrY".to_string()).await?,
        };
        // Then keep only off-backbone, callable, non-shared calls.
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

#[cfg(test)]
mod gvcf_discovery_tests {
    use super::*;

    fn alignment_at(bam: &Path) -> Alignment {
        Alignment {
            id: 1,
            sequence_run_id: 1,
            bam_path: Some(bam.to_string_lossy().into_owned()),
            reference_path: None,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            content_sha256: None,
            derived_from_alignment_id: None,
            derivation: None,
        }
    }

    /// Each case gets its own directory — these run in parallel.
    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nav-gvcf-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn finds_the_ytree_sibling_sidecar() {
        let d = scratch("sibling");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::write(d.join("HG002.chrY.g.vcf.gz"), "").unwrap();
        let got = chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).unwrap();
        assert_eq!(got.file_name().unwrap(), "HG002.chrY.g.vcf.gz");
    }

    #[test]
    fn finds_a_bare_named_gvcf_in_a_caller_subdirectory() {
        // The D2C per-run layout: `<run>/CP086569.2/gatk4/chrY.g.vcf.gz`, with no sample prefix and
        // one directory down. Matching only `*.chry.g.vcf.gz` beside the CRAM found none of these,
        // so every subject fell back to decoding the whole chromosome.
        let d = scratch("subdir");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::create_dir_all(d.join("gatk4")).unwrap();
        std::fs::write(d.join("gatk4").join("chrY.g.vcf.gz"), "").unwrap();
        let got = chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).unwrap();
        assert_eq!(got.file_name().unwrap(), "chrY.g.vcf.gz");
        assert_eq!(got.parent().unwrap().file_name().unwrap(), "gatk4");
    }

    #[test]
    fn the_alignments_own_directory_wins_over_a_subdirectory() {
        let d = scratch("precedence");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::write(d.join("HG002.chrY.g.vcf.gz"), "").unwrap();
        std::fs::create_dir_all(d.join("gatk4")).unwrap();
        std::fs::write(d.join("gatk4").join("chrY.g.vcf.gz"), "").unwrap();
        let got = chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).unwrap();
        assert_eq!(got.file_name().unwrap(), "HG002.chrY.g.vcf.gz");
    }

    #[test]
    fn a_chr_m_gvcf_is_not_mistaken_for_chr_y() {
        let d = scratch("contig");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::create_dir_all(d.join("gatk4")).unwrap();
        std::fs::write(d.join("gatk4").join("chrM.g.vcf.gz"), "").unwrap();
        let aln = alignment_at(&d.join("chrYM.cram"));
        assert!(
            chr_y_gvcf_for_alignment(&aln).is_none(),
            "chrM must not satisfy a chrY lookup"
        );
        assert!(chr_m_gvcf_for_alignment(&aln).is_some());
    }

    #[test]
    fn absent_sidecars_yield_none() {
        let d = scratch("absent");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        assert!(chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).is_none());
    }
}

/// Read depth a VCF call must show before it can be a private-variant candidate. Matches the GVCF
/// path's floor: below this, a call is far more likely a misaligned-read cluster than a real SNV.
const VCF_PRIVATE_MIN_DP: u32 = 4;

/// Genotype-quality floor, matching the GVCF path's `min_gq`. Aligns the two sources' gates as far
/// as the evidence allows — though it does not make them comparable: see
/// [`App::private_y_from_variant_set`] on why a vendor caller's call set is a different instrument.
const VCF_PRIVATE_MIN_GQ: u32 = 20;

/// Derived-allele fraction a call must reach to count as **deterministic** on a haploid chromosome.
/// chrY carries one copy, so a genuine call is essentially all-alt; a middling fraction is an
/// ambiguous locus, and an ambiguous call can not support a private-variant claim.
const VCF_PRIVATE_MIN_AF: f64 = 0.95;

/// Depth ceiling, as a multiple of the donor's own typical depth at good calls.
///
/// chrY carries one copy, so a locus drawing far more reads than the rest of the chromosome is
/// collecting them from somewhere else — a collapsed repeat. This was found by reviewing a candidate
/// branch whose two carriers sat at DP 413 and 504 against a median of 57, each holding a stubborn
/// ~5% reference allele: the shape of a paralogous pile-up, and it cleared every other gate. Three
/// times the median keeps ~91% of quality calls while removing that tail.
const VCF_PRIVATE_MAX_DEPTH_RATIO: u32 = 3;

/// Minimum quality-passing calls before a depth ratio is trustworthy. Below this the median is not a
/// description of the donor's coverage, and the rule abstains rather than judging against noise.
const VCF_PRIVATE_MIN_CALLS_FOR_RATIO: usize = 20;

impl App {
    /// The **private bucket for a variant set** — the VCF counterpart of [`Self::private_y_variants`].
    ///
    /// Private-Y has always been keyed on an alignment: it walks a BAM/CRAM (or its GVCF sidecar) and
    /// caches against `alignment_id`. A subject whose Y data arrived as an externally processed VCF
    /// has no alignment, so the option was never offered — on R1b-CTS4466Plus that is ~1,600 of 1,881
    /// members, and it is why cohort features that depend on private variants had almost nothing to
    /// work with.
    ///
    /// The classification is deliberately the same as the alignment path's: subtract the placed
    /// backbone, drop anything outside the cohort callable mask or on the cohort-shared blocklist,
    /// then split off-path-known from novel. What differs is where the evidence comes from and how
    /// the donor's own reliability is judged:
    ///
    /// - **Placement** uses [`Self::vset_base_calls`], so the terminal is derived from tree-position
    ///   genotypes (including hom-ref) rather than the handful of derived calls alone.
    /// - **There is no self-callable mask** — a VCF carries no coverage track. Its place is taken by
    ///   the source's own per-call evidence: a `FILTER`-flagged call is dropped, as is one below
    ///   [`VCF_PRIVATE_MIN_DP`], and a chrY heterozygote is dropped outright — on a haploid
    ///   chromosome that is a paralog or mismapping artefact, and it is ~2/3 of a Big Y's chrY rows.
    /// - **Depth and allele fraction are the source's**, so [`PublishGate`] judges these calls on real
    ///   read evidence. A set imported before evidence capture (`call_schema` 1) therefore yields
    ///   nothing publishable, which is the honest outcome rather than a fabricated one.
    ///
    /// **Not comparable to the alignment path's counts, and not yet fit for branch inference.** This
    /// yields a median ~175 novel calls per donor against the GVCF path's 3–13. The gap is the
    /// instrument, not a defect here: the alignment path reads GATK HaplotypeCaller at ploidy 1,
    /// while a vendor export is a diploid caller emitting far more chrY calls, and only ~10% of the
    /// difference is reachable by matching DP/GQ gates. Note also that `Novel` means "not
    /// branch-defining in *this* tree" — the tree is FTDNA's supported branches plus splits solved
    /// from the cohort, not a catalogue of known Y variation — so a real, well-known variant that
    /// defines no branch classifies as novel here. Feeding these buckets to the block tree's
    /// candidate detection took CTS4466 from 3 candidates to 20 (39 conflicts, 105 recurrent
    /// positions dropped) on only 111 of ~1,600 sets, which is why
    /// [`Self::private_y_for_biosamples`] does not union them yet.
    pub async fn private_y_from_variant_set(&self, set: &VariantSet) -> Result<PrivateBucket, AppError> {
        use navigator_analysis::haplo;

        // Without per-call evidence every quality gate below is a no-op, and the result is a list of
        // whatever the vendor's caller emitted — on a real set that is 400-550 "novel" calls against
        // ~70 for the same donor's evidence-bearing set. A call we can not judge is the most
        // non-deterministic kind there is, so refuse rather than publish a number that looks like a
        // finding. Re-importing the source populates `CallEvidence` (migration 0042).
        if !set.has_evidence() {
            return Err(AppError::Import(format!(
                "variant set {} carries no per-call evidence (call_schema {}); re-import it to enable private-Y",
                set.id, set.call_schema
            )));
        }
        let build = set.reference_build.clone().unwrap_or_else(|| "GRCh38".to_string());
        let tree = self.chip_y_tree(&build).await?;

        let cache_key = {
            let targets: HashSet<i64> = tree
                .nodes
                .values()
                .flat_map(|n| n.loci.iter().map(|l| l.position))
                .collect();
            // `pv2`: the chrY structural masks are now lifted to the set's own build, so a `pv1`
            // bucket was classified with **no** structural mask on anything but CHM13 and its counts
            // are inflated. The version is the invalidation — `--force` can not reach this cache.
            format!("pv2:{}", crate::haplogroup::genotype_cache_key("chrY", None, &targets))
        };
        if let Ok(Some(json)) = variant_set_private_y::get(self.store.pool(), set.id, &cache_key).await {
            if let Ok(bucket) = serde_json::from_str::<PrivateBucket>(&json) {
                return Ok(bucket);
            }
        }

        // Place the set to get its backbone. Tree-position genotypes (with ancestral evidence) place
        // far deeper than the derived calls alone would.
        let genotypes = self.vset_base_calls(set, "chrY", &tree).await;
        let assignment = assemble_assignment(&tree, &genotypes);
        let Some(terminal) = assignment.ranked.first() else {
            return Err(AppError::Import("no Y haplogroup match for this call set".into()));
        };
        let path = haplo::path_positions(&tree, terminal.id);
        let known = haplo::tree_positions(&tree);

        let mask_token = y_mask_build_token(&build);
        let cohort_mask =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_CALLABLE_MASK", "chrY_callable_mask", t));
        let cohort_shared =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_COHORT_SHARED", "chrY_cohort_shared_sites", t));
        // The structural BEDs are CHM13-native and lifted to whatever this set is in — without that
        // a GRCh38 set has no structural mask and its private counts inflate into the hundreds.
        let regions = self.y_structural_regions_for(&build).await;

        // Quality-passing calls first, so the depth ceiling below is measured against the donor's own
        // good coverage rather than against a median dragged down by the junk we are about to drop.
        let passing: Vec<&navigator_domain::variants::VariantCall> = set
            .calls
            .iter()
            .filter(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"))
            .filter(|c| !path.contains(&c.position))
            .filter(|c| cohort_mask.as_ref().map_or(true, |m| m.contains(c.position)))
            .filter(|c| cohort_shared.as_ref().map_or(true, |m| !m.contains(c.position)))
            .filter(|c| is_hemizygous(c.genotype.as_deref()))
            .filter(|c| !c.evidence.is_filtered())
            .filter(|c| !c.evidence.dp.is_some_and(|dp| dp < VCF_PRIVATE_MIN_DP))
            .filter(|c| !c.evidence.gq.is_some_and(|gq| gq < VCF_PRIVATE_MIN_GQ))
            .filter(|c| !c.evidence.allele_fraction().is_some_and(|af| af < VCF_PRIVATE_MIN_AF))
            .collect();
        let depth_ceiling = median_depth(&passing).map(|m| m * VCF_PRIVATE_MAX_DEPTH_RATIO);

        let mut variants: Vec<PrivateVariant> = passing
            .into_iter()
            .filter(|c| match (depth_ceiling, c.evidence.dp) {
                (Some(ceiling), Some(dp)) => dp <= ceiling,
                _ => true,
            })
            .filter_map(|c| {
                let (reference, alternate) = (c.reference.chars().next()?, c.alternate.chars().next()?);
                Some(PrivateVariant {
                    position: c.position,
                    reference,
                    alternate,
                    // The source's own numbers; absent when it gave none, which the publish gate
                    // then (correctly) refuses rather than treating as evidence.
                    depth: c.evidence.dp.unwrap_or(0),
                    alt_depth: c.evidence.ad_alt.unwrap_or(0),
                    allele_fraction: c.evidence.allele_fraction().unwrap_or(0.0),
                    class: match known.get(&c.position) {
                        Some(name) => PrivateClass::OffPathKnown(name.clone()),
                        None => PrivateClass::Novel,
                    },
                    region: regions.as_ref().and_then(|r| r.classify(c.position)),
                })
            })
            .collect();
        variants.sort_by_key(|v| v.position);

        let bucket = PrivateBucket {
            terminal: terminal.name.clone(),
            variants,
        };
        if let Ok(json) = serde_json::to_string(&bucket) {
            let _ = variant_set_private_y::upsert(self.store.pool(), set.id, &cache_key, &json).await;
        }
        Ok(bucket)
    }
}

/// Whether a genotype is a single-allele (hemizygous / homozygous-alt) call.
///
/// chrY is haploid, so a heterozygous call there has no biological reading: it is a paralogous or
/// mismapped locus. In a real Big Y export those are ~2/3 of the chrY rows, and admitting them would
/// make the private set mostly artefact. An absent genotype is admitted — a source that reports no GT
/// is not asserting heterozygosity.
fn is_hemizygous(gt: Option<&str>) -> bool {
    let Some(gt) = gt else { return true };
    let alleles: Vec<&str> = gt.split(['/', '|']).filter(|a| *a != ".").collect();
    !alleles.is_empty() && alleles.windows(2).all(|w| w[0] == w[1])
}

#[cfg(test)]
mod vcf_private_y_tests {
    use super::is_hemizygous;

    #[test]
    fn a_chr_y_heterozygote_is_rejected() {
        // chrY is haploid: a het call is a paralog or a mismapping, and it is ~2/3 of a Big Y's
        // chrY rows — admitting them would make the private set mostly artefact.
        assert!(!is_hemizygous(Some("0/1")));
        assert!(!is_hemizygous(Some("1|2")));
        assert!(!is_hemizygous(Some("1/2")));
    }

    #[test]
    fn hemizygous_and_homozygous_alt_calls_are_kept() {
        assert!(is_hemizygous(Some("1")));
        assert!(is_hemizygous(Some("1/1")));
        assert!(is_hemizygous(Some("1|1")));
        assert!(is_hemizygous(Some("2/2")));
    }

    #[test]
    fn a_source_that_reports_no_genotype_is_not_treated_as_heterozygous() {
        // Absence of a GT is not an assertion about ploidy; rejecting it would silently discard
        // every sites-only or CSV-derived set.
        assert!(is_hemizygous(None));
        assert!(is_hemizygous(Some("1/.")), "a partial call still carries one allele");
    }
}

/// Median read depth across `calls`, or `None` when too few carry one to describe the donor's
/// coverage. Median rather than mean: the pile-ups this exists to find would drag a mean upward and
/// hide themselves behind it.
fn median_depth(calls: &[&navigator_domain::variants::VariantCall]) -> Option<u32> {
    let mut depths: Vec<u32> = calls.iter().filter_map(|c| c.evidence.dp).collect();
    if depths.len() < VCF_PRIVATE_MIN_CALLS_FOR_RATIO {
        return None;
    }
    depths.sort_unstable();
    Some(depths[depths.len() / 2]).filter(|&m| m > 0)
}

#[cfg(test)]
mod depth_ratio_tests {
    use super::{median_depth, VCF_PRIVATE_MIN_CALLS_FOR_RATIO};
    use navigator_domain::variants::{CallEvidence, VariantCall};

    fn call(dp: Option<u32>) -> VariantCall {
        VariantCall {
            contig: "chrY".into(),
            position: 1,
            reference: "A".into(),
            alternate: "G".into(),
            rs_id: None,
            genotype: Some("1/1".into()),
            evidence: CallEvidence {
                dp,
                ..Default::default()
            },
        }
    }

    #[test]
    fn the_median_ignores_the_pile_ups_it_exists_to_find() {
        // A mean would be dragged up by the outliers and hide them behind itself.
        let mut calls: Vec<VariantCall> = (0..VCF_PRIVATE_MIN_CALLS_FOR_RATIO).map(|_| call(Some(50))).collect();
        calls.push(call(Some(2584)));
        calls.push(call(Some(1191)));
        let refs: Vec<&VariantCall> = calls.iter().collect();
        assert_eq!(median_depth(&refs), Some(50));
    }

    #[test]
    fn it_abstains_on_too_few_calls_to_describe_coverage() {
        let calls: Vec<VariantCall> = (0..VCF_PRIVATE_MIN_CALLS_FOR_RATIO - 1)
            .map(|_| call(Some(50)))
            .collect();
        let refs: Vec<&VariantCall> = calls.iter().collect();
        assert_eq!(median_depth(&refs), None, "no ceiling rather than one built on noise");
    }

    #[test]
    fn a_source_reporting_no_depth_yields_no_ceiling() {
        // Absent depth must not read as zero and reject everything.
        let calls: Vec<VariantCall> = (0..40).map(|_| call(None)).collect();
        let refs: Vec<&VariantCall> = calls.iter().collect();
        assert_eq!(median_depth(&refs), None);
    }
}
