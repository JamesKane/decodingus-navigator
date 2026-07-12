//! `impl App` methods extracted from `lib.rs` (the `analysis` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- analysis (compute + persist) --------------------------------------

    /// Run the coverage + callable walker on an alignment's BAM and persist the result
    /// as a versioned `coverage` artifact. The blocking noodles I/O runs on a blocking
    /// thread so the async runtime is not stalled.
    pub async fn run_coverage(
        &self,
        alignment_id: i64,
        bam: PathBuf,
        reference: PathBuf,
        contig_allowlist: Option<HashSet<String>>,
        params: CallableLociParams,
    ) -> Result<CoverageResult, AppError> {
        let result = tokio::task::spawn_blocking(move || {
            navigator_analysis::guard_walk("coverage", || {
                coverage::collect_coverage_callable(&bam, &reference, &params, contig_allowlist.as_ref())
            })
        })
        .await??;
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result)
            .await?;
        Ok(result)
    }

    /// Cached `coverage` result for the current algorithm version, if present.
    pub async fn cached_coverage(&self, alignment_id: i64) -> Result<Option<CoverageResult>, AppError> {
        self.load_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION)
            .await
    }

    /// Run coverage using the alignment's own stored BAM/reference paths, then persist.
    /// Errors if the alignment is unknown or has no paths recorded.
    pub async fn run_coverage_for_alignment(&self, alignment_id: i64) -> Result<CoverageResult, AppError> {
        self.run_coverage_for_alignment_with_progress(alignment_id, |_, _| {})
            .await
    }

    /// Like [`run_coverage_for_alignment`], reporting `progress(contigs_done, contigs_total)` as
    /// the whole-genome pass walks each contig (the slow step — minutes on a real WGS BAM — so a
    /// progress bar can advance instead of sitting frozen). The callback runs on the blocking
    /// thread.
    pub async fn run_coverage_for_alignment_with_progress(
        &self,
        alignment_id: i64,
        mut progress: impl FnMut(usize, usize) + Send + 'static,
    ) -> Result<CoverageResult, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        // The reference isn't asked for at import — resolve the alignment's build via the gateway
        // (cached, else download) when no FASTA was stored.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        // For a targeted test (Big Y, etc.) restrict the walk to the target chromosome(s) so the
        // headline depth reflects the target rather than being diluted to ~0 by the empty genome.
        let allowlist = self.coverage_target_allowlist(alignment_id).await?;
        let mut params = CallableLociParams::default();
        let result = tokio::task::spawn_blocking(move || {
            // Adapt the callable threshold to read tech (HiFi → 1×; see adaptive_min_depth).
            if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(&bam, Some(&reference)) {
                params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            }
            navigator_analysis::guard_walk("coverage", || {
                coverage::collect_coverage_callable_with_progress(
                    &bam,
                    &reference,
                    &params,
                    allowlist.as_ref(),
                    &mut progress,
                )
            })
        })
        .await??;
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result)
            .await?;
        // A generic FTDNA Targeted-Y run now has a callable footprint → pin Big Y-500 vs -700.
        self.refine_big_y_generation_for_alignment(alignment_id, &result)
            .await?;
        Ok(result)
    }

    /// The coverage contig allowlist for a targeted test, or `None` (whole genome) for WGS/autosomal.
    /// A Y-targeted test (FTDNA Big Y, Y Elite, …) walks chrY only — plus chrM so the "has mtDNA
    /// reads" signal survives for the few Big Ys that retained mitochondrial reads (the UI hides the
    /// mtDNA sections when chrM has none). An mtDNA-targeted test walks chrM only. Build-agnostic
    /// (both `chr`-prefixed and bare contig names are listed).
    async fn coverage_target_allowlist(&self, alignment_id: i64) -> Result<Option<HashSet<String>>, AppError> {
        use navigator_domain::testtype::TargetType;
        let aln = self.alignment_or_err(alignment_id).await?;
        let Some(run) = sequence_run::get(self.store.pool(), aln.sequence_run_id).await? else {
            return Ok(None);
        };
        // `target_of` (not bare `by_code`) so a stored human label like "Big Y" — which a bulk
        // import / --test-type override writes instead of BIG_Y_500/700 — still scopes the walk to
        // chrY+chrM. Otherwise coverage walks the whole genome, which on a targeted multi-reference
        // CRAM is the ~1-hour batch-analysis stall.
        let contigs: &[&str] = match navigator_domain::testtype::target_of(&run.test_type) {
            Some(TargetType::YChromosome) => &["chrY", "Y", "chrM", "chrMT", "M", "MT"],
            Some(TargetType::MtDna) => &["chrM", "chrMT", "M", "MT"],
            _ => return Ok(None),
        };
        Ok(Some(contigs.iter().map(|s| s.to_string()).collect()))
    }

    /// Whether a cached coverage result was computed at the right scope for the alignment's test.
    /// A targeted test (Big Y, mtFull) must cover only its target contig(s); a whole-genome cached
    /// result for it is stale — the headline depth was diluted across the empty genome — and must
    /// be recomputed. Whole-genome tests (no allowlist) are always in scope.
    pub(crate) async fn coverage_is_correctly_scoped(
        &self,
        alignment_id: i64,
        cov: &CoverageResult,
    ) -> Result<bool, AppError> {
        match self.coverage_target_allowlist(alignment_id).await? {
            None => Ok(true),
            Some(allow) => Ok(cov.contig_coverage_stats.iter().all(|s| allow.contains(&s.contig))),
        }
    }

    /// Cached coverage for analysis reuse: the stored result, but only when it was computed at the
    /// right scope for the test (see [`Self::coverage_is_correctly_scoped`]). A stale whole-genome
    /// result for a targeted test reads as a cache miss so the caller recomputes it correctly.
    pub async fn cached_coverage_for_analysis(&self, alignment_id: i64) -> Result<Option<CoverageResult>, AppError> {
        match self.cached_coverage(alignment_id).await? {
            Some(cov) if self.coverage_is_correctly_scoped(alignment_id, &cov).await? => Ok(Some(cov)),
            _ => Ok(None),
        }
    }

    /// Infer biological sex from the alignment's chrX:autosome read-density ratio, persisting
    /// the result as a `sex` artifact. Cheap (BAI fast-path for BAM). `reference` is used only
    /// for CRAM decode.
    pub async fn run_sex(&self, alignment_id: i64) -> Result<navigator_analysis::sex::SexInferenceResult, AppError> {
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let result =
            tokio::task::spawn_blocking(move || navigator_analysis::sex::infer_from_bam(&bam, reference.as_deref()))
                .await??;
        self.save_analysis(alignment_id, "sex", "1", &result).await?;
        self.write_back_inferred_sex(alignment_id, &result).await?;
        Ok(result)
    }

    /// Write the inferred sex back to the biosample when the user didn't provide one, so it
    /// shows in the subjects table + header instead of "Unknown". No-op for Unknown sex or
    /// when the biosample already carries a sex.
    pub(crate) async fn write_back_inferred_sex(
        &self,
        alignment_id: i64,
        result: &navigator_analysis::sex::SexInferenceResult,
    ) -> Result<(), AppError> {
        let label = match result.inferred_sex {
            InferredSex::Male => Some("Male"),
            InferredSex::Female => Some("Female"),
            InferredSex::Unknown => None,
        };
        if let (Some(label), Ok(guid)) = (label, self.biosample_of_alignment(alignment_id).await) {
            if let Ok(Some(bio)) = biosample::get(self.store.pool(), guid).await {
                if bio.sex.as_deref().map(str::trim).unwrap_or("").is_empty() {
                    biosample::set_sex(self.store.pool(), guid, label).await?;
                }
            }
        }
        Ok(())
    }

    /// Cached `sex` inference, if present.
    pub async fn cached_sex(
        &self,
        alignment_id: i64,
    ) -> Result<Option<navigator_analysis::sex::SexInferenceResult>, AppError> {
        self.load_analysis(alignment_id, "sex", "1").await
    }

    /// Collect read-level QC metrics (alignment summary + read-length/insert-size distributions,
    /// pair orientation, mean MAPQ) and persist as a `read_metrics` artifact.
    pub async fn run_read_metrics(
        &self,
        alignment_id: i64,
    ) -> Result<navigator_analysis::read_metrics::ReadMetrics, AppError> {
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let result = tokio::task::spawn_blocking(move || {
            navigator_analysis::read_metrics::collect_read_metrics(&bam, reference.as_deref())
        })
        .await??;
        self.save_analysis(alignment_id, "read_metrics", "1", &result).await?;
        self.write_back_read_stats(alignment_id, &result).await?;
        Ok(result)
    }

    /// Mirror an alignment's library-level read stats onto its owning sequence run (`total_reads`,
    /// `mean_read_length`, `mean_insert_size`) so the Data Sources run card shows them without
    /// re-walking. Best-effort: a missing alignment/run is ignored. When a run has several
    /// alignments the last write wins — these are per-library properties, so any pass is
    /// representative.
    pub(crate) async fn write_back_read_stats(
        &self,
        alignment_id: i64,
        m: &navigator_analysis::read_metrics::ReadMetrics,
    ) -> Result<(), AppError> {
        if let Some(aln) = alignment::get(self.store.pool(), alignment_id).await? {
            // Paired-end evidence: any reads aligned in pairs ⇒ PAIRED. Only overrides the stored
            // layout when we have aligned reads to judge (else leave the import-time flag value).
            let layout = (m.pf_reads_aligned > 0).then_some(if m.reads_aligned_in_pairs > 0 {
                "PAIRED"
            } else {
                "SINGLE"
            });
            sequence_run::set_read_stats(
                self.store.pool(),
                aln.sequence_run_id,
                Some(m.total_reads as i64),
                (m.mean_read_length > 0.0).then_some(m.mean_read_length),
                (m.mean_insert_size > 0.0).then_some(m.mean_insert_size),
                layout,
                // Exact sequenced yield (Σ read_length_histogram) → the "Gbases" figure of the
                // standardized test label. `None` (empty histogram, no fallback) leaves the column.
                m.total_bases(),
            )
            .await?;
        }
        Ok(())
    }

    /// Cached `read_metrics`, if present.
    pub async fn cached_read_metrics(
        &self,
        alignment_id: i64,
    ) -> Result<Option<navigator_analysis::read_metrics::ReadMetrics>, AppError> {
        self.load_analysis(alignment_id, "read_metrics", "1").await
    }

    /// Local on-disk cache for alignments copied off a slow/removable volume (see [`localize`]).
    pub(crate) fn align_cache_dir() -> std::path::PathBuf {
        navigator_refgenome::cache::base_dir().join("cache").join("aln")
    }

    /// If `remote` lives on a slow/removable volume (a `/Volumes/…` mount), copy it — and its `.crai`
    /// / `.bai` index — into the local cache and return the *local* path; otherwise return `remote`
    /// unchanged. The analysis walkers do random-access record iteration (region seeks, per-read
    /// decode), which is pathologically slow over a network/USB mount even though a plain sequential
    /// **copy** of the same file is fast — so we pay one fast bulk copy up front and let every
    /// subsequent pass read from local disk. The copy is reused across a subject's passes and cleared
    /// per subject by [`clear_align_cache`]. A copy failure falls back to the remote path (slow, but
    /// still works).
    pub(crate) async fn localize(&self, remote: &Path) -> PathBuf {
        if std::env::var_os("NAVIGATOR_NO_LOCALIZE").is_some() || !is_removable_volume(remote) {
            return remote.to_path_buf();
        }
        let cache = Self::align_cache_dir();
        let local = cache.join(local_cache_name(remote));
        if local.is_file() {
            return local;
        }
        let (remote_owned, local2) = (remote.to_path_buf(), local.clone());
        match tokio::task::spawn_blocking(move || copy_with_index(&remote_owned, &local2)).await {
            Ok(Ok(())) => local,
            Ok(Err(e)) => {
                eprintln!("localize: copy failed ({e}); reading from the original (slow)");
                remote.to_path_buf()
            }
            Err(e) => {
                eprintln!("localize: copy task failed ({e}); reading from the original (slow)");
                remote.to_path_buf()
            }
        }
    }

    /// Drop the local alignment cache (called per subject so the batch holds at most one file).
    pub(crate) fn clear_align_cache() {
        let _ = std::fs::remove_dir_all(Self::align_cache_dir());
    }

    /// Run the unified quality-metrics walker — coverage + callable, read-level QC metrics, and
    /// sex inference in **one pass** over the alignment's BAM/CRAM (vs. the separate passes
    /// `run_coverage` + `run_read_metrics` + `run_sex` cost: 2 reads for BAM, 3 for CRAM). All
    /// three sub-results are persisted under their existing artifact keys (`coverage`/
    /// `COVERAGE_VERSION`, `read_metrics`/`"1"`, `sex`/`"1"`), so `cached_coverage`/
    /// `cached_read_metrics`/`cached_sex` and the SV step's reuse logic keep working unchanged.
    pub async fn run_unified_metrics(&self, alignment_id: i64) -> Result<UnifiedMetricsResult, AppError> {
        self.run_unified_metrics_with_progress(alignment_id, |_, _| {}).await
    }

    /// Like [`run_unified_metrics`], reporting `progress(contigs_done, contigs_total)` as the
    /// (slow) whole-genome coverage portion finalizes each contig. Uses the per-contig parallel
    /// walker (falling back to a sequential pass for CRAM / unindexed BAM); the callback is
    /// `Fn + Sync` because it's invoked concurrently from the fan-out's worker threads.
    pub async fn run_unified_metrics_with_progress(
        &self,
        alignment_id: i64,
        progress: impl Fn(usize, usize) + Send + Sync + 'static,
    ) -> Result<UnifiedMetricsResult, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let run_id = aln.sequence_run_id;
        // Copy off a slow/removable volume to local disk first — the walker's random-access record
        // iteration is far slower over a network/USB mount than a one-shot bulk copy.
        let bam = self
            .localize(Path::new(&aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?))
            .await;
        // The walker requires a reference (CRAM decode + reference-N detection); resolve the
        // build via the gateway when no FASTA was stored at import.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        // Restrict a targeted test (Big Y, mtFull) to its target contig(s), exactly like the
        // standalone coverage walker — otherwise the headline depth is diluted across the empty
        // genome (a Big Y reads as ~0.2× instead of ~50× on chrY). WGS keeps the whole-genome walk.
        let allowlist = self.coverage_target_allowlist(alignment_id).await?;
        let mut params = CallableLociParams::default();
        let result = tokio::task::spawn_blocking(move || {
            // Adapt the callable threshold to read tech (HiFi → 1×; see adaptive_min_depth).
            if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(&bam, Some(&reference)) {
                params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            }
            navigator_analysis::guard_walk("metrics", || {
                navigator_analysis::unified::collect_unified_metrics_parallel_with_progress(
                    &bam,
                    &reference,
                    &params,
                    allowlist.as_ref(),
                    &progress,
                )
            })
        })
        .await??;

        // Persist each sub-result under its own existing cache key.
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result.coverage)
            .await?;
        // A generic FTDNA Targeted-Y run now has a callable footprint → pin Big Y-500 vs -700.
        self.refine_big_y_generation_for_alignment(alignment_id, &result.coverage)
            .await?;
        self.save_analysis(alignment_id, "read_metrics", "1", &result.read_metrics)
            .await?;
        self.write_back_read_stats(alignment_id, &result.read_metrics).await?;
        // Sex: a Y-targeted test (Big Y, Y Elite, …) sequences the donor's Y chromosome — he is male
        // by definition. The chrX/autosome ratio the inference needs isn't present in a chrY-scoped
        // walk, and is unreliable even whole-genome (a Big Y's off-target chrX ≈ autosome ≈ 0.4×
        // reads as *female*). So force Male for a Y-targeted test, overriding the inference + any
        // prior auto-assignment; WGS / mt-targeted keep the walk's result.
        let y_targeted = matches!(
            sequence_run::get(self.store.pool(), run_id)
                .await?
                .as_ref()
                .and_then(|r| navigator_domain::testtype::target_of(&r.test_type)),
            Some(navigator_domain::testtype::TargetType::YChromosome)
        );
        // A Y-scoped alignment reads as male the same way a Y-targeted test does — chrY carries
        // essentially all the reads while the autosomes hold only a few dozen mismapped ones (a
        // Y-only extract, e.g. GRCh38 chrY reads realigned to hs1, or a Y-Elite/Big Y capture that
        // came in mislabeled WGS). The ratio walk can then read it as *female*, which silently
        // disables the whole Y pipeline (assign_y_haplogroup skips females before it ever fetches
        // the tree). Detect it from the per-contig read counts and force male, exactly like a Y test.
        let y_scoped = navigator_analysis::sex::is_y_scoped(
            result
                .coverage
                .contig_coverage_stats
                .iter()
                .map(|s| (s.contig.as_str(), s.num_reads)),
        );
        let male_by_scope = y_targeted || y_scoped;
        let sex = if male_by_scope {
            Some(navigator_analysis::sex::SexInferenceResult {
                inferred_sex: navigator_analysis::sex::InferredSex::Male,
                x_autosome_ratio: 0.0,
                autosome_mean_coverage: 0.0,
                x_coverage: 0.0,
                confidence: navigator_analysis::sex::Confidence::High,
            })
        } else {
            result.sex
        };
        if let Some(sex) = &sex {
            self.save_analysis(alignment_id, "sex", "1", sex).await?;
            if male_by_scope {
                // Definitive (Y test / Y-scoped ⇒ male): override any prior auto-inferred sex —
                // including a stale false "Female" — rather than write-if-empty.
                if let Ok(guid) = self.biosample_of_alignment(alignment_id).await {
                    biosample::set_sex(self.store.pool(), guid, "Male").await?;
                }
            } else {
                self.write_back_inferred_sex(alignment_id, sex).await?;
            }
        }
        Ok(result)
    }

    /// Call structural variants (depth-segmentation + paired-end/split-read evidence) and
    /// persist as an `sv` artifact. Needs coverage + insert-size inputs (computed/loaded here)
    /// and **≥10× mean coverage** (the caller errors below that).
    pub async fn run_sv(&self, alignment_id: i64) -> Result<navigator_analysis::sv::types::SvAnalysisResult, AppError> {
        // Resume: a fresh cached SV result (source unchanged) is reused rather than recomputed.
        if let Some(c) = self.cached_sv(alignment_id).await? {
            return Ok(c);
        }
        let aln = self.alignment_or_err(alignment_id).await?;
        let reference_build = aln.reference_build.clone();
        // Resolve the reference for decode (see alignment_reference_for_decode): required for a CRAM,
        // None for a BAM. SV reads records + header lengths; it doesn't consult reference bases.
        let (bam, reference) = self.alignment_reference_for_decode(alignment_id).await?;

        let cov = match self.cached_coverage(alignment_id).await? {
            Some(c) => c,
            None => self.run_coverage_for_alignment(alignment_id).await?,
        };
        let rm = match self.cached_read_metrics(alignment_id).await? {
            Some(m) => m,
            None => self.run_read_metrics(alignment_id).await?,
        };
        let (mean_cov, mean_ins, sd_ins, mean_rl) = (
            cov.mean_coverage,
            rm.mean_insert_size,
            rm.std_insert_size,
            rm.mean_read_length,
        );

        let result = tokio::task::spawn_blocking(move || {
            let lengths = caller::header_contig_lengths(&bam, reference.as_deref())?;
            navigator_analysis::guard_walk("structural variants", || {
                navigator_analysis::sv::caller::call_structural_variants(
                    &bam,
                    &lengths,
                    &reference_build,
                    mean_cov,
                    mean_ins,
                    sd_ins,
                    mean_rl,
                    &navigator_analysis::sv::types::SvCallerConfig::default(),
                )
            })
        })
        .await??;
        self.save_analysis(alignment_id, "sv", "1", &result).await?;
        Ok(result)
    }

    /// Cached `sv` result, if present.
    pub async fn cached_sv(
        &self,
        alignment_id: i64,
    ) -> Result<Option<navigator_analysis::sv::types::SvAnalysisResult>, AppError> {
        self.load_analysis(alignment_id, "sv", "1").await
    }

    /// The HipSTR-format STR reference BED for `reference_build`, if available: the explicit
    /// `NAVIGATOR_STR_REFERENCE` path (env override), else `~/.decodingus/str/{build}.hipstr_reference.bed.gz`.
    /// `None` → the caller surfaces a "configure the STR reference" error.
    fn str_reference_path(reference_build: &str) -> Option<PathBuf> {
        if let Ok(p) = std::env::var("NAVIGATOR_STR_REFERENCE") {
            let p = PathBuf::from(p);
            return p.exists().then_some(p);
        }
        // Use the shared cache base (honors `NAVIGATOR_REFGENOME_DIR`) so this matches where
        // `seed_bundled_str` places the bundled reference.
        let p = navigator_refgenome::cache::base_dir()
            .join("str")
            .join(format!("{reference_build}.hipstr_reference.bed.gz"));
        p.exists().then_some(p)
    }

    /// Genotype short tandem repeats on `contig` from the alignment, via the enclosing-read caller
    /// over the HipSTR reference tracts (haploid for chrY/chrM, diploid elsewhere). Persisted as a
    /// `str:{contig}` artifact (so it's cached + source-invalidated like other analyses). Errors if
    /// no STR reference is configured for the alignment's build (the tracts are build-specific —
    /// CHM13/GRCh37 need their own reference or liftover, not yet wired).
    pub async fn run_str_calls(
        &self,
        alignment_id: i64,
        contig: String,
    ) -> Result<Vec<navigator_analysis::strcaller::StrGenotype>, AppError> {
        let kind = format!("str:{contig}");
        if let Some(c) = self.load_analysis(alignment_id, &kind, "str-1").await? {
            return Ok(c);
        }
        let aln = self.alignment_or_err(alignment_id).await?;
        let build = aln.reference_build.clone();
        let bed = Self::str_reference_path(&build).ok_or_else(|| {
            AppError::Import(format!(
                "no STR reference for build {build} — set NAVIGATOR_STR_REFERENCE to a HipSTR BED, \
                 or place it at ~/.decodingus/str/{build}.hipstr_reference.bed.gz"
            ))
        })?;
        // Resolve the reference for decode (see alignment_reference_for_decode): required for a CRAM,
        // None for a BAM. STR region-genotyping reads the alignment; it doesn't consult reference bases.
        let (bam, reference) = self.alignment_reference_for_decode(alignment_id).await?;
        // chrY / chrM are haploid (one allele); autosomes + chrX (in a female) are diploid. We
        // genotype chrY/chrM haploid and everything else diploid — sex-aware chrX is a refinement.
        let ploidy: u8 = {
            let c = contig.strip_prefix("chr").unwrap_or(&contig).to_ascii_uppercase();
            if c == "Y" || c == "M" || c == "MT" {
                1
            } else {
                2
            }
        };
        let params = navigator_analysis::strcaller::StrCallerParams::default();
        let genos = tokio::task::spawn_blocking(move || {
            let loci = navigator_analysis::strref::load_hipstr_contig(&bed, &contig, 2)?;
            navigator_analysis::strcaller::genotype_str_loci(
                &bam,
                &contig,
                &loci,
                ploidy,
                &params,
                reference.as_deref(),
            )
        })
        .await??;
        self.save_analysis(alignment_id, &kind, "str-1", &genos).await?;
        Ok(genos)
    }

    /// Compare the STR markers called from sequence (mapped to the FTDNA convention via the
    /// corpus-calibrated [`navigator_analysis::strmarker`] table) against the subject's imported
    /// vendor Y-STR profile — the By-Panel concordance view. One row per marker present in either
    /// source: the called value + its calibration status, the imported value, and whether they agree.
    /// `contig` is typically `chrY`. Reuses the cached `str:{contig}` calls.
    pub async fn str_concordance(&self, alignment_id: i64, contig: String) -> Result<Vec<StrConcordanceRow>, AppError> {
        use navigator_analysis::strmarker::{called_markers_build, normalize_marker, MarkerStatus, StrBuild};

        // The FTDNA convention offset is build-dependent for a few markers (the CHM13 liftover shifted
        // some tract boundaries) — select the offsets for this alignment's build.
        let build = alignment::get(self.store.pool(), alignment_id)
            .await?
            .map(|a| StrBuild::from_build_str(&a.reference_build))
            .unwrap_or_default();
        let genos = self.run_str_calls(alignment_id, contig).await?;
        let called = called_markers_build(&genos, build);

        // Imported vendor markers (FTDNA preferred, else the first profile), keyed by normalized name.
        let biosample = self.biosample_of_alignment(alignment_id).await?;
        let profiles = self.list_str_profiles(biosample).await?;
        let chosen = profiles
            .iter()
            .find(|p| p.provider.as_deref().is_some_and(|v| v.eq_ignore_ascii_case("FTDNA")))
            .or_else(|| profiles.first());
        let imported: HashMap<String, String> = chosen
            .map(|p| {
                p.markers
                    .iter()
                    .map(|m| (normalize_marker(&m.marker), m.value.clone()))
                    .collect()
            })
            .unwrap_or_default();

        let mut rows: HashMap<String, StrConcordanceRow> = HashMap::new();
        for c in &called {
            rows.insert(
                c.marker.clone(),
                StrConcordanceRow {
                    marker: c.marker.clone(),
                    called: Some(c.value),
                    status: format!("{:?}", c.status),
                    calibrated: matches!(c.status, MarkerStatus::Reliable | MarkerStatus::ConventionOffset),
                    imported: imported.get(&c.marker).cloned(),
                    depth: c.depth,
                    agree: false,
                },
            );
        }
        for (m, v) in &imported {
            rows.entry(m.clone()).or_insert_with(|| StrConcordanceRow {
                marker: m.clone(),
                called: None,
                status: "NotCalled".into(),
                calibrated: false,
                imported: Some(v.clone()),
                depth: 0,
                agree: false,
            });
        }
        // Agreement: a calibrated call whose value matches the imported single value.
        let mut out: Vec<StrConcordanceRow> = rows
            .into_values()
            .map(|mut r| {
                r.agree =
                    r.calibrated && matches!((&r.called, &r.imported), (Some(c), Some(i)) if i.trim() == c.to_string());
                r
            })
            .collect();
        out.sort_by(|a, b| a.marker.cmp(&b.marker));
        Ok(out)
    }

    /// Pick the subject's best STR-capable alignment and run the Y-STR concordance on chrY — the
    /// entry point the UI calls. "STR-capable" = an alignment whose reference build has a HipSTR
    /// reference present ([`str_reference_path`](Self::str_reference_path)); highest mean coverage
    /// wins. A CRAM needs no stored reference here — [`run_str_calls`](Self::run_str_calls) resolves
    /// it for decode. Errors with guidance when none qualifies (no HipSTR reference / no alignment).
    pub async fn str_concordance_for_subject(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<(i64, Vec<StrConcordanceRow>), AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut best: Option<(i64, f64)> = None;
        for a in &alns {
            if Self::str_reference_path(&a.reference_build).is_none() {
                continue; // no HipSTR reference for this build
            }
            // A CRAM with no stored reference is fine — run_str_calls resolves it via the gateway.
            let cov = self
                .cached_coverage(a.id)
                .await
                .ok()
                .flatten()
                .map(|c| c.mean_coverage)
                .unwrap_or(0.0);
            if best.as_ref().map_or(true, |(_, bc)| cov > *bc) {
                best = Some((a.id, cov));
            }
        }
        let (aln_id, _) = best.ok_or_else(|| {
            AppError::Import(
                "no STR-capable alignment — need a GRCh38/CHM13 BAM or CRAM and the HipSTR \
                 reference at ~/.decodingus/str/{build}.hipstr_reference.bed.gz (or \
                 NAVIGATOR_STR_REFERENCE)"
                    .into(),
            )
        })?;
        let rows = self.str_concordance(aln_id, "chrY".into()).await?;
        Ok((aln_id, rows))
    }

    /// The alignment's BAM (required) + a reference for decoding it (see
    /// [`alignment_reference_for_decode`](Self::alignment_reference_for_decode)): resolved for a CRAM,
    /// `None` for a BAM. Coverage / read-metrics / callable read records but never consult reference
    /// bases, so a BAM needs none.
    pub(crate) async fn alignment_paths(&self, alignment_id: i64) -> Result<(PathBuf, Option<PathBuf>), AppError> {
        self.alignment_reference_for_decode(alignment_id).await
    }

    /// Run de-novo haploid calling on a contig and persist the SNP calls as a versioned
    /// `denovo_snps` artifact.
    pub async fn run_denovo_caller(
        &self,
        alignment_id: i64,
        bam: PathBuf,
        reference: PathBuf,
        contig: String,
        params: HaploidCallerParams,
    ) -> Result<Vec<VariantCall>, AppError> {
        // Resume: reuse a fresh cached de-novo result for this contig (source unchanged).
        if let Some(c) = self.cached_denovo(alignment_id, &contig).await? {
            return Ok(c);
        }
        let kind = denovo_kind(&contig);
        let calls = tokio::task::spawn_blocking(move || {
            navigator_analysis::guard_walk("de-novo calling", || caller::call_denovo(&bam, &reference, &contig, &params))
        })
        .await??;
        self.save_analysis(alignment_id, &kind, caller::DENOVO_VERSION, &calls)
            .await?;
        Ok(calls)
    }

    /// Cached de-novo calls for `contig` at the current caller version, if present.
    pub async fn cached_denovo(&self, alignment_id: i64, contig: &str) -> Result<Option<Vec<VariantCall>>, AppError> {
        self.load_analysis(alignment_id, &denovo_kind(contig), caller::DENOVO_VERSION)
            .await
    }

    /// Whole-contig **de-novo diploid** SNV calling (het 0/1 + hom-alt 1/1) on `contig`, cached per
    /// alignment+contig. Reuses the alignment's BAM + reference (resolved from the build). Returns
    /// [`SiteGenotype`]s in position order — feed to [`Self::diploid_vcf`].
    pub async fn run_diploid_calls(&self, alignment_id: i64, contig: String) -> Result<Vec<SiteGenotype>, AppError> {
        let kind = format!("diploid_denovo:{contig}");
        if let Some(c) = self
            .load_analysis(alignment_id, &kind, caller::GENOTYPE_VERSION)
            .await?
        {
            return Ok(c);
        }
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let params = adaptive_haploid_params(&bam, Some(&reference));
        let calls = tokio::task::spawn_blocking(move || {
            navigator_analysis::guard_walk("diploid calling", || {
                caller::call_denovo_diploid(&bam, &reference, &contig, &params)
            })
        })
        .await??;
        self.save_analysis(alignment_id, &kind, caller::GENOTYPE_VERSION, &calls)
            .await?;
        Ok(calls)
    }

    /// A diploid VCF (VCFv4.2, `GT:AD:DP:GQ:PL`) of the de-novo diploid SNV calls for `contig`
    /// (computing + caching them if needed). The sample column is `aln<id>`.
    pub async fn diploid_vcf(&self, alignment_id: i64, contig: String) -> Result<String, AppError> {
        let calls = self.run_diploid_calls(alignment_id, contig).await?;
        Ok(navigator_analysis::vcf::write_diploid_vcf(
            &format!("aln{alignment_id}"),
            &calls,
        ))
    }

    /// A **whole-genome** diploid VCF: de-novo SNV + indel calls over the diploid primary
    /// chromosomes (1–22, X) of the alignment, per-contig cached. chrY and chrM are **excluded** —
    /// they're haploid, so the diploid (het 0/1) model is wrong for them; their variants come from
    /// the haploid caller and the Y/mt haplogroup + mtDNA-mutation features. Heavy (a real WGS
    /// calling pass); the caller runs it off the UI thread (the export path).
    pub async fn diploid_vcf_genome(&self, alignment_id: i64) -> Result<String, AppError> {
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let contigs =
            tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, Some(&reference))).await??;
        let mut all = Vec::new();
        for contig in contigs
            .into_iter()
            .filter(|c| is_primary_contig(c) && !is_haploid_contig(c))
        {
            all.extend(self.run_diploid_calls(alignment_id, contig).await?);
        }
        Ok(navigator_analysis::vcf::write_diploid_vcf(
            &format!("aln{alignment_id}"),
            &all,
        ))
    }

    /// The subject's alignments on the **dominant reference build** (the build the most alignments
    /// share, compared on the canonical build so `chm13v2`/`hs1` agree). The consensus diploid
    /// genotype pools only same-build alignments — de-novo variant coordinates can't be merged
    /// across builds by position without genome-wide liftover (out of scope). `None` if no alignments.
    async fn consensus_diploid_alignments(&self, biosample_guid: SampleGuid) -> Result<Vec<i64>, AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        if alns.is_empty() {
            return Ok(Vec::new());
        }
        let mut counts: HashMap<Option<ReferenceBuild>, usize> = HashMap::new();
        for a in &alns {
            *counts.entry(canonical_build(&a.reference_build)).or_default() += 1;
        }
        let dominant = counts
            .into_iter()
            .max_by_key(|(_, n)| *n)
            .map(|(b, _)| b)
            .unwrap_or(None);
        Ok(alns
            .into_iter()
            .filter(|a| canonical_build(&a.reference_build) == dominant)
            .map(|a| a.id)
            .collect())
    }

    /// **Subject-level consensus** diploid genotype across the subject's same-build WGS alignments —
    /// the joint genotype (opportunity #3). Per [`reconcile_site_genotypes`]: call each alignment's
    /// variants (cached [`run_diploid_calls`]), union the SNV sites, force-genotype **every**
    /// alignment at the union (so a site absent from one run is its real hom-ref / no-call), and vote
    /// a depth-weighted 0/1/2 dosage per site. Returns the variant (het/hom-alt) consensus sites.
    /// `contigs` limits the scan (None = all primary chromosomes). Heavy (a call pass + a force-call
    /// pass per alignment) — an explicit export action; nothing is persisted.
    pub async fn consensus_diploid_calls(
        &self,
        biosample_guid: SampleGuid,
        contigs: Option<Vec<String>>,
    ) -> Result<Vec<SiteGenotype>, AppError> {
        let aln_ids = self.consensus_diploid_alignments(biosample_guid).await?;
        if aln_ids.is_empty() {
            return Ok(Vec::new());
        }

        // (bam, reference) per same-build alignment, resolved once.
        let mut paths = Vec::new();
        for id in &aln_ids {
            paths.push((*id, self.alignment_bam_reference(*id).await?));
        }

        // 1–2. Call each alignment's variants and union the SNV sites (force-call is SNP-only).
        let mut union: HashMap<(String, i64, String), navigator_analysis::caller::Site> = HashMap::new();
        for (id, (bam, reference)) in &paths {
            let clist = match &contigs {
                Some(c) => c.clone(),
                None => {
                    let bam = bam.clone();
                    let reference = reference.clone();
                    tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, Some(&reference)))
                        .await??
                        .into_iter()
                        .filter(|c| is_primary_contig(c))
                        .collect()
                }
            };
            for contig in clist {
                // Tolerate a contig absent from this alignment's header (heterogeneous inputs) —
                // skip it for this source rather than aborting the whole consensus.
                let Ok(variants) = self.run_diploid_calls(*id, contig).await else {
                    continue;
                };
                for v in variants {
                    if v.reference_allele.len() == 1 && v.alternate_allele.len() == 1 {
                        union
                            .entry((v.contig.clone(), v.position, v.alternate_allele.clone()))
                            .or_insert(navigator_analysis::caller::Site {
                                name: String::new(),
                                contig: v.contig,
                                position: v.position,
                                reference_allele: v.reference_allele,
                                alternate_allele: v.alternate_allele,
                            });
                    }
                }
            }
        }
        if union.is_empty() {
            return Ok(Vec::new());
        }
        let sites: Vec<navigator_analysis::caller::Site> = union.into_values().collect();

        // 3. Force-genotype every alignment at the union (each emits hom-ref / no-call too).
        let mut per_aln: Vec<Vec<SiteGenotype>> = Vec::new();
        for (_, (bam, reference)) in &paths {
            let params = adaptive_haploid_params(bam, Some(reference));
            let (bam, reference, sites) = (bam.clone(), reference.clone(), sites.clone());
            let g = tokio::task::spawn_blocking(move || {
                caller::genotype_sites_all_contigs(&bam, &sites, 2, &params, Some(&reference))
            })
            .await??;
            per_aln.push(g);
        }

        // 4. Vote per site → consensus. min_depth = 2: a run abstains only when essentially
        // uncovered; depth-weighting lets deep runs dominate the rest.
        Ok(caller::reconcile_site_genotypes(&per_aln, 2))
    }

    /// A **consensus** diploid VCF (VCFv4.2) for the subject — the joint genotype across same-build
    /// alignments (see [`consensus_diploid_calls`]), sample column `consensus`. Heavy; the export
    /// path runs it off the UI thread.
    pub async fn consensus_diploid_vcf(&self, biosample_guid: SampleGuid) -> Result<String, AppError> {
        let calls = self.consensus_diploid_calls(biosample_guid, None).await?;
        Ok(navigator_analysis::vcf::write_diploid_vcf("consensus", &calls))
    }

    /// Run de-novo calling on `contig` using the alignment's own stored paths.
    /// The alignment's BAM + a usable reference FASTA: the stored path, else resolved from the
    /// alignment's build via the gateway (cached, else downloaded). Errors only if no BAM is
    /// recorded. Use this in steps that *require* the reference, so the user never has to supply
    /// one (it follows from the header-detected build).
    pub(crate) async fn alignment_bam_reference(&self, alignment_id: i64) -> Result<(PathBuf, PathBuf), AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        Ok((bam, reference))
    }

    /// The alignment's path and a reference suitable for **decoding** it: a CRAM can't be read
    /// without the reference, so resolve it (stored path, else from the build via the gateway,
    /// cache-first); a BAM decodes without one, so return the stored path as-is (usually `None`) and
    /// never force a reference download. Use this for record/pileup reads and SNP-site genotyping
    /// that don't consult reference bases; use [`alignment_bam_reference`](Self::alignment_bam_reference)
    /// for calling paths (de-novo SNV/indel) that need the reference even on a BAM.
    pub(crate) async fn alignment_reference_for_decode(
        &self,
        alignment_id: i64,
    ) -> Result<(PathBuf, Option<PathBuf>), AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let is_cram = bam.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        let reference = match aln.reference_path {
            Some(p) => Some(PathBuf::from(p)),
            None if is_cram => Some(
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?,
            ),
            None => None,
        };
        Ok((bam, reference))
    }

    /// Whether the reference FASTA for `build` is already on disk (no download needed). Lets the UI
    /// worker decide when a reference resolution would trigger a visible download vs. a cache hit.
    pub fn reference_cached(&self, build: &str) -> bool {
        self.gateway.cached_reference(build).is_some()
    }

    /// The distinct reference builds across a subject's alignments — the builds whose FASTA an
    /// analysis of this subject may need. Used to pre-resolve references (with a progress bar) after
    /// import and before a subject-level analysis, so on-demand downloads aren't silent.
    pub async fn reference_builds_for_subject(&self, biosample_guid: SampleGuid) -> Result<Vec<String>, AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut builds: Vec<String> = alns.into_iter().map(|a| a.reference_build).collect();
        builds.sort();
        builds.dedup();
        Ok(builds)
    }

    /// The reference build of a single alignment (`None` if it no longer exists) — for pre-resolving
    /// that alignment's reference before a per-alignment analysis.
    pub async fn reference_build_of_alignment(&self, alignment_id: i64) -> Result<Option<String>, AppError> {
        Ok(alignment::get(self.store.pool(), alignment_id)
            .await?
            .map(|a| a.reference_build))
    }

    /// The alignment IDs (BAM/CRAM only) across a subject's alignments — for pre-building each one's
    /// coordinate index (with a progress bar) after import and before a subject-level analysis.
    pub async fn alignment_ids_for_subject(&self, biosample_guid: SampleGuid) -> Result<Vec<i64>, AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        Ok(alns
            .into_iter()
            .filter(|a| a.bam_path.is_some())
            .map(|a| a.id)
            .collect())
    }

    /// Ensure the alignment's coordinate index (`.bai`/`.crai`) exists, **building it if missing** so
    /// the query-driven analyses (the per-contig walker, callable intervals, the de-novo / STR
    /// callers) can seek by region instead of erroring or degrading to a whole-file linear scan.
    /// Returns the index path if one was built, `None` if it was already present. `progress(done,
    /// total)` reports a byte fraction for a BAM and indeterminate progress (`total = None`) for a
    /// CRAM. The build is a single sequential pass, run on a decode-safe blocking thread.
    pub async fn ensure_alignment_index(
        &self,
        alignment_id: i64,
        mut progress: impl FnMut(u64, Option<u64>) + Send + 'static,
    ) -> Result<Option<PathBuf>, AppError> {
        let (bam, reference) = self.alignment_reference_for_decode(alignment_id).await?;
        let built = tokio::task::spawn_blocking(move || {
            navigator_analysis::guard_walk("build index", || {
                navigator_analysis::index::ensure_index(&bam, reference.as_deref(), &mut progress)
            })
        })
        .await??;
        Ok(built)
    }

    pub async fn run_denovo_for_alignment(
        &self,
        alignment_id: i64,
        contig: String,
    ) -> Result<Vec<VariantCall>, AppError> {
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let probe = bam.clone();
        let probe_ref = reference.clone();
        let params = tokio::task::spawn_blocking(move || adaptive_haploid_params(&probe, Some(&probe_ref))).await?; // HiFi -> lower min_depth
        self.run_denovo_caller(alignment_id, bam, reference, contig, params)
            .await
    }

    /// The [`PublishGate`] for an alignment, adapted to its mean read length (HiFi relaxes the
    /// supporting-read floor — see [`PublishGate::for_read_len`]). Samples the BAM head; any error
    /// falls back to the short-read default.
    pub async fn publish_gate_for_alignment(&self, alignment_id: i64) -> Result<PublishGate, AppError> {
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let read_len = tokio::task::spawn_blocking(move || {
            navigator_analysis::coverage::estimate_molecule_lengths(&bam, Some(&reference))
                .map(|(rl, _)| rl)
                .unwrap_or(0.0)
        })
        .await?;
        Ok(PublishGate::for_read_len(read_len))
    }
}

/// True for a path on a removable/network mount (macOS `/Volumes/…`), where per-record random access
/// is slow but a bulk sequential copy is fast — the case [`App::localize`] copies to local disk.
fn is_removable_volume(p: &Path) -> bool {
    p.starts_with("/Volumes/")
}

/// A collision-free local filename for a remote alignment. Every kit's file is named `chrYM.cram`,
/// so the basename alone collides; hash the full remote path and keep the extension so the reader
/// still finds the sibling index at `<local>.crai` / `<local>.bai`.
fn local_cache_name(remote: &Path) -> String {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    remote.to_string_lossy().hash(&mut h);
    let ext = remote.extension().and_then(|e| e.to_str()).unwrap_or("bam");
    format!("{:016x}.{ext}", h.finish())
}

/// Copy `remote` → `local` plus its index sibling. The index is copied **first** and the main file
/// last (via a temp + rename), so a present `local` always implies its index is present too — the
/// `is_file` cache check in [`App::localize`] can't see a half-copied pair.
fn copy_with_index(remote: &Path, local: &Path) -> std::io::Result<()> {
    if let Some(parent) = local.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let (rstr, lstr) = (remote.to_string_lossy(), local.to_string_lossy());
    for suffix in [".crai", ".bai"] {
        let ri = PathBuf::from(format!("{rstr}{suffix}"));
        if ri.is_file() {
            std::fs::copy(&ri, PathBuf::from(format!("{lstr}{suffix}")))?;
        }
    }
    // BAM index sometimes drops the .bam: `<stem>.bai`.
    let ri = remote.with_extension("bai");
    if ri.is_file() {
        std::fs::copy(&ri, local.with_extension("bai"))?;
    }
    let tmp = local.with_extension("partial");
    std::fs::copy(remote, &tmp)?;
    std::fs::rename(&tmp, local)?;
    Ok(())
}
