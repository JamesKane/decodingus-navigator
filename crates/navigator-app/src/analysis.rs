//! `impl App` methods extracted from `lib.rs` (the `analysis` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;
use navigator_analysis::{contig, CancelToken};

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
        self.run_coverage_for_alignment_with_progress(alignment_id, |_, _| {}, CancelToken::none())
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
        cancel: CancelToken,
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
                    &cancel,
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

    /// Scratch directory for alignments copied off a slow/removable volume (see [`localize`]).
    /// Entries are owned by a [`LocalAlignment`] and removed when the last holder drops.
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
    pub(crate) async fn localize(&self, remote: &Path) -> LocalAlignment {
        if std::env::var_os("NAVIGATOR_NO_LOCALIZE").is_some() || !is_removable_volume(remote) {
            return LocalAlignment::borrowed(remote);
        }
        let local = Self::align_cache_dir().join(local_cache_name(remote));

        // Serialize the cache-check-then-copy per destination. The worker `tokio::spawn`s every
        // command, so a batch walk and a per-alignment command genuinely overlap on one alignment;
        // without this both miss the cache and both copy, writing a second full 40 GB pull over the
        // network for nothing, after which the loser of the rename reads from the remote anyway.
        //
        // The gate is held across the copy, so it must not be taken re-entrantly — no `localize`
        // may be called while another is outstanding *on the same path in the same task*. The three
        // call sites are sequential today (`debug_y_calls` awaits `base_calls` to completion before
        // localizing itself); keep it that way.
        let gate = copy_gate(&local);
        let _copying = gate.lock().await;

        // Size the remote once: it decides both whether an existing copy can be trusted and whether
        // the one we make arrived whole.
        let remote_len = tokio::fs::metadata(remote).await.ok().map(|m| m.len());

        // Another holder is already using this copy — share it and bump the count.
        if LocalAlignment::retain(&local, remote_len) {
            return LocalAlignment::owned(local);
        }
        let (remote_owned, local2) = (remote.to_path_buf(), local.clone());
        match tokio::task::spawn_blocking(move || copy_with_index(&remote_owned, &local2, remote_len)).await {
            // Registering can still fail if the copy was removed in the gap (a concurrent holder
            // finishing and dropping to zero). Returning an `owned` handle to a missing path would
            // fail the walk with a confusing ENOENT and then "clean up" a file that isn't there.
            Ok(Ok(())) if LocalAlignment::retain(&local, remote_len) => LocalAlignment::owned(local),
            Ok(Ok(())) => {
                eprintln!(
                    "localize: the copy at {} went away before it could be used; reading from the original (slow)",
                    local.display()
                );
                LocalAlignment::borrowed(remote)
            }
            Ok(Err(e)) => {
                eprintln!("localize: copy failed ({e}); reading from the original (slow)");
                LocalAlignment::borrowed(remote)
            }
            Err(e) => {
                eprintln!("localize: copy task failed ({e}); reading from the original (slow)");
                LocalAlignment::borrowed(remote)
            }
        }
    }

    /// Run the unified quality-metrics walker — coverage + callable, read-level QC metrics, and
    /// sex inference in **one pass** over the alignment's BAM/CRAM (vs. the separate passes
    /// `run_coverage` + `run_read_metrics` + `run_sex` cost: 2 reads for BAM, 3 for CRAM). All
    /// three sub-results are persisted under their existing artifact keys (`coverage`/
    /// `COVERAGE_VERSION`, `read_metrics`/`"1"`, `sex`/`"1"`), so `cached_coverage`/
    /// `cached_read_metrics`/`cached_sex` and the SV step's reuse logic keep working unchanged.
    pub async fn run_unified_metrics(&self, alignment_id: i64) -> Result<UnifiedMetricsResult, AppError> {
        self.run_unified_metrics_with_progress(alignment_id, |_, _| {}, CancelToken::none())
            .await
    }

    /// Like [`run_unified_metrics`], reporting `progress(contigs_done, contigs_total)` as the
    /// (slow) whole-genome coverage portion finalizes each contig. Uses the per-contig parallel
    /// walker (falling back to a sequential pass for CRAM / unindexed BAM); the callback is
    /// `Fn + Sync` because it's invoked concurrently from the fan-out's worker threads.
    pub async fn run_unified_metrics_with_progress(
        &self,
        alignment_id: i64,
        progress: impl Fn(usize, usize) + Send + Sync + 'static,
        cancel: CancelToken,
    ) -> Result<UnifiedMetricsResult, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let run_id = aln.sequence_run_id;
        // Copy off a slow/removable volume to local disk first — the walker's random-access record
        // iteration is far slower over a network/USB mount than a one-shot bulk copy.
        // Held for the whole walk: dropping it removes the local copy.
        let bam = self
            .localize(Path::new(&aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?))
            .await;
        let bam = bam.path().to_path_buf();
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
                    &cancel,
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
    pub async fn run_sv(
        &self,
        alignment_id: i64,
        cancel: CancelToken,
    ) -> Result<navigator_analysis::sv::types::SvAnalysisResult, AppError> {
        // Resume: a fresh cached SV result (source unchanged) is reused rather than recomputed.
        if let Some(c) = self.cached_sv(alignment_id).await? {
            return Ok(c);
        }
        let aln = self.alignment_or_err(alignment_id).await?;
        let reference_build = aln.reference_build.clone();
        // Resolve the reference for decode (see alignment_reference_for_decode): required for a CRAM,
        // None for a BAM. SV never consults reference *bases* — but decoding a CRAM record does, so
        // the walker needs it too, not just the header-lengths probe.
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
                    reference.as_deref(),
                    &lengths,
                    &reference_build,
                    mean_cov,
                    mean_ins,
                    sd_ins,
                    mean_rl,
                    &navigator_analysis::sv::types::SvCallerConfig::default(),
                    &cancel,
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
        let ploidy: u8 = if contig::is_haploid(&contig) { 1 } else { 2 };
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
        cancel: CancelToken,
    ) -> Result<Vec<VariantCall>, AppError> {
        // Resume: reuse a fresh cached de-novo result for this contig (source unchanged).
        if let Some(c) = self.cached_denovo(alignment_id, &contig).await? {
            return Ok(c);
        }
        let kind = denovo_kind(&contig);
        let calls = tokio::task::spawn_blocking(move || {
            navigator_analysis::guard_walk("de-novo calling", || {
                caller::call_denovo(&bam, &reference, &contig, &params, &cancel)
            })
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
    pub async fn run_diploid_calls(
        &self,
        alignment_id: i64,
        contig: String,
        cancel: CancelToken,
    ) -> Result<Vec<SiteGenotype>, AppError> {
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
                caller::call_denovo_diploid(&bam, &reference, &contig, &params, &cancel)
            })
        })
        .await??;
        self.save_analysis(alignment_id, &kind, caller::GENOTYPE_VERSION, &calls)
            .await?;
        Ok(calls)
    }

    /// A diploid VCF (VCFv4.2, `GT:AD:DP:GQ:PL`) of the de-novo diploid SNV calls for `contig`
    /// (computing + caching them if needed). The sample column is `aln<id>`.
    pub async fn diploid_vcf(
        &self,
        alignment_id: i64,
        contig: String,
        cancel: CancelToken,
    ) -> Result<String, AppError> {
        let calls = self.run_diploid_calls(alignment_id, contig, cancel).await?;
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
    pub async fn diploid_vcf_genome(&self, alignment_id: i64, cancel: CancelToken) -> Result<String, AppError> {
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let contigs =
            tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, Some(&reference))).await??;
        let mut all = Vec::new();
        for contig in contigs
            .into_iter()
            .filter(|c| contig::is_main_assembly(c) && !contig::is_haploid(c))
        {
            all.extend(self.run_diploid_calls(alignment_id, contig, cancel.clone()).await?);
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
    pub(crate) async fn consensus_diploid_alignments(&self, biosample_guid: SampleGuid) -> Result<Vec<i64>, AppError> {
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
        cancel: CancelToken,
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
                        .filter(|c| contig::is_main_assembly(c))
                        .collect()
                }
            };
            for contig in clist {
                // Tolerate a contig absent from this alignment's header (heterogeneous inputs) —
                // skip it for this source rather than aborting the whole consensus.
                let Ok(variants) = self.run_diploid_calls(*id, contig, cancel.clone()).await else {
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
            let cancel = cancel.clone();
            let g = tokio::task::spawn_blocking(move || {
                caller::genotype_sites_all_contigs(&bam, &sites, 2, &params, Some(&reference), &cancel)
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
        let calls = self
            .consensus_diploid_calls(biosample_guid, None, CancelToken::none())
            .await?;
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

    /// Diagnose why an alignment can't be read, naming the **exact file** at fault rather than the
    /// one the failing call happened to be handed. See [`navigator_analysis::preflight`] for why
    /// that distinction is the whole point: an unreadable `.crai` and an unreadable CRAM produce
    /// the same `io error on …cram` message today, and on macOS a privacy (TCC) denial and a Unix
    /// permission denial are told apart only by the raw errno.
    ///
    /// Deliberately **cache-only** for the reference: a diagnostic has to describe the machine as
    /// it is, so resolving (and silently downloading) a missing FASTA here would paper over exactly
    /// the state we were asked to report. A CRAM with no cached reference is a finding, not a task.
    pub async fn diagnose_alignment(
        &self,
        alignment_id: i64,
    ) -> Result<navigator_analysis::preflight::Report, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = match aln.reference_path {
            Some(p) => Some(PathBuf::from(p)),
            None => self.gateway.cached_reference(&aln.reference_build),
        };
        Ok(
            tokio::task::spawn_blocking(move || navigator_analysis::preflight::diagnose(&bam, reference.as_deref()))
                .await?,
        )
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
        self.run_denovo_caller(alignment_id, bam, reference, contig, params, CancelToken::none())
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

/// A scratch name for one caller's in-progress copy of `local`. **Unique per call**: a temp path
/// derived from the destination alone (`<dest>.partial`) is shared by every concurrent copier of
/// that alignment, which lets two of them open the same inode — one truncating what the other is
/// writing, and continuing to write into it after the other has renamed it into place and started
/// reading. Uniqueness makes that unrepresentable rather than merely unlikely.
fn partial_path(local: &Path) -> PathBuf {
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let stem = local.file_name().unwrap_or_default().to_string_lossy().into_owned();
    local.with_file_name(format!("{stem}.partial.{}.{n}", std::process::id()))
}

/// Copy `remote` → `local` plus its index sibling. The index is copied **first** and the main file
/// last (via a temp + rename), so a present `local` always implies its index is present too — the
/// cache check in [`App::localize`] can't see a half-copied pair.
///
/// `expect_len` is the remote's size; when known, a copy that doesn't match it is rejected rather
/// than published. A short copy is otherwise indistinguishable from corrupt data: it surfaces as a
/// decode error ("unexpected end of file", a bad container checksum) tens of gigabytes into a walk,
/// naming the cache path, and the copy is deleted on drop before anyone can look at it.
///
/// Nothing partial survives a failure — neither the temp nor the index copied ahead of it. A
/// leftover `.partial` used to sit in the cache indefinitely, occupying the disk while satisfying
/// no one.
fn copy_with_index(remote: &Path, local: &Path, expect_len: Option<u64>) -> std::io::Result<()> {
    if let Some(parent) = local.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let (rstr, lstr) = (remote.to_string_lossy(), local.to_string_lossy());
    let mut indexes: Vec<PathBuf> = Vec::new();
    let mut copy_index = |from: PathBuf, to: PathBuf| -> std::io::Result<()> {
        if from.is_file() {
            std::fs::copy(&from, &to)?;
            indexes.push(to);
        }
        Ok(())
    };
    let mut copied_indexes = || -> std::io::Result<()> {
        for suffix in [".crai", ".bai"] {
            copy_index(
                PathBuf::from(format!("{rstr}{suffix}")),
                PathBuf::from(format!("{lstr}{suffix}")),
            )?;
        }
        // BAM index sometimes drops the .bam: `<stem>.bai`.
        copy_index(remote.with_extension("bai"), local.with_extension("bai"))
    };

    let tmp = partial_path(local);
    let result = copied_indexes().and_then(|()| {
        std::fs::copy(remote, &tmp)?;
        let got = std::fs::metadata(&tmp)?.len();
        match expect_len {
            Some(want) if got != want => Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!("copy of {} is {got} bytes, expected {want}", remote.display()),
            )),
            _ => std::fs::rename(&tmp, local),
        }
    });
    if result.is_err() {
        let _ = std::fs::remove_file(&tmp);
        for i in &indexes {
            let _ = std::fs::remove_file(i);
        }
    }
    result
}

/// One step of a full analysis of a single alignment, in the order [`App::plan_full_analysis`]
/// returns them. The variants carry whatever the step needs, so a caller's dispatch is total and it
/// cannot silently run a step the plan excluded.
#[derive(Debug, Clone, PartialEq)]
pub enum AnalysisStep {
    /// Coverage + callable, read-level QC, and sex inference in one pass over the alignment.
    QualityMetrics,
    /// CNV + discordant pairs. Needs ≥10× — the step itself reports when the depth is too low.
    ///
    /// **Opt-in only.** SV is experimental and, alone among the steps, walks every read in the file
    /// for a result nothing else consumes — hours per whole-genome sample. It is planned only when
    /// a caller asks for it (`include_sv`); the GUI's "Call SV" button and `analyze --sv` are the
    /// ways in. Nothing runs it unattended.
    StructuralVariants,
    /// De-novo calling on the mitochondrial contig (small and fully callable, unlike whole chrY).
    MitoDenovo { contig: String },
    /// Place the alignment on the Y tree.
    YHaplogroup,
    /// Place the alignment on the mtDNA tree.
    MtHaplogroup,
    /// Genome-consensus Y signature (deep placement + variant profile → descent report).
    YSignature { biosample_guid: SampleGuid },
    /// Genotype the ancestry markers into the autosomal consensus profile.
    AutosomalProfile { biosample_guid: SampleGuid },
    /// Estimate admixture/PCA *from* the autosomal profile — must follow [`Self::AutosomalProfile`].
    Ancestry { biosample_guid: SampleGuid },
}

impl AnalysisStep {
    /// Short name for a progress indicator.
    pub fn label(&self) -> &'static str {
        match self {
            Self::QualityMetrics => "Quality metrics",
            Self::StructuralVariants => "Structural variants",
            Self::MitoDenovo { .. } => "Variant calling",
            Self::YHaplogroup => "Y haplogroup",
            Self::MtHaplogroup => "mtDNA haplogroup",
            Self::YSignature { .. } => "Y signature",
            Self::AutosomalProfile { .. } => "Autosomal profile",
            Self::Ancestry { .. } => "Ancestry",
        }
    }

    /// One line of detail under the label.
    pub fn detail(&self) -> &'static str {
        match self {
            Self::QualityMetrics => "scanning contigs…",
            Self::StructuralVariants => "CNV + discordant pairs (needs ≥10×)",
            Self::MitoDenovo { .. } => "chrM de-novo (haploid)",
            Self::YHaplogroup => "placing on the Y tree",
            Self::MtHaplogroup => "placing on the mt tree",
            Self::YSignature { .. } => "building the descent report",
            Self::AutosomalProfile { .. } => "genotyping ancestry markers",
            Self::Ancestry { .. } => "estimating admixture + ancient components",
        }
    }
}

impl App {
    /// The ordered steps a full analysis of `alignment_id` should run.
    ///
    /// This is the **one** definition of that pipeline. The GUI streams progress events and the CLI
    /// prints a log, so how a step is reported differs — but which steps run, and the conditions
    /// under which one is skipped, must not: the two had drifted, and the CLI's copy re-genotyped Y
    /// unconditionally, overwriting a trusted external call (on ancient DNA, with a worse one).
    ///
    /// `coverage` is the just-computed result when [`AnalysisStep::QualityMetrics`] has already run,
    /// which makes the mitochondrial decision authoritative; pass `None` before that and the cached
    /// coverage (or, with none, the assumption that chrM is present) is used instead. Callers that
    /// show a step count should re-plan after the metrics step, as the count can drop.
    ///
    /// `include_sv` adds the experimental [`AnalysisStep::StructuralVariants`]; see that variant for
    /// why it is off by default. `include_ancestry` likewise gates the two heaviest ancestry steps.
    pub async fn plan_full_analysis(
        &self,
        alignment_id: i64,
        include_ancestry: bool,
        include_sv: bool,
        coverage: Option<&CoverageResult>,
    ) -> Result<Vec<AnalysisStep>, AppError> {
        // Skip the mitochondrial steps when the alignment has no chrM reads (e.g. an FTDNA Big Y):
        // scoring zero chrM data just records a meaningless RSRS root. Unknown coverage (never run)
        // keeps them rather than silently skipping.
        let has_mtdna = match coverage {
            Some(c) => chrm_has_reads(c),
            None => self
                .cached_coverage(alignment_id)
                .await
                .ok()
                .flatten()
                .map(|c| chrm_has_reads(&c))
                .unwrap_or(true),
        };
        // Subject-level steps need the owning subject; an unattached alignment simply skips them.
        let guid = self.biosample_of_alignment(alignment_id).await.ok();

        let mut steps = vec![AnalysisStep::QualityMetrics];
        if include_sv {
            steps.push(AnalysisStep::StructuralVariants);
        }
        if has_mtdna {
            steps.push(AnalysisStep::MitoDenovo { contig: "chrM".into() });
        }
        // Skip the internal Y/mt genotyping when a trusted external caller (GATK4 GVCF, sidecar fast
        // path) already placed this alignment and the user prefers it — re-walking would only produce
        // a secondary call that loses the vote, and on ancient DNA a wrong one. The assign_* commands
        // guard this too; planning around it also avoids the wasted decode. See external-caller-precedence.
        if !self
            .has_preferred_external_call(alignment_id, DnaType::Y)
            .await
            .unwrap_or(false)
        {
            steps.push(AnalysisStep::YHaplogroup);
        }
        if has_mtdna
            && !self
                .has_preferred_external_call(alignment_id, DnaType::Mt)
                .await
                .unwrap_or(false)
        {
            steps.push(AnalysisStep::MtHaplogroup);
        }
        // The Y signature makes the descent report ready without an explicit click (the button stays
        // for an on-demand rebuild). Built once — an existing profile is left alone — and it needs no
        // extra read of the file, since the Y assignment above just cached the chrY genotypes.
        if let Some(guid) = guid {
            if self.cached_y_profile(guid).await?.is_none() {
                steps.push(AnalysisStep::YSignature { biosample_guid: guid });
            }
        }
        // The one-click Simple flow folds in the autosomal ancestry that drives "Your ancestry".
        // Two ordered steps (profile → estimate), last because they are the heaviest. Advanced runs
        // ancestry as a separate deliberate action and passes `include_ancestry = false`.
        if let (true, Some(guid)) = (include_ancestry, guid) {
            steps.push(AnalysisStep::AutosomalProfile { biosample_guid: guid });
            steps.push(AnalysisStep::Ancestry { biosample_guid: guid });
        }
        Ok(steps)
    }
}

/// Whether a coverage result shows any reads on the mitochondrial contig, under either naming.
fn chrm_has_reads(c: &CoverageResult) -> bool {
    c.contig_coverage_stats
        .iter()
        .any(|s| contig::is_chr_m(&s.contig) && s.num_reads > 0)
}

/// Live localized copies: local path → how many [`LocalAlignment`]s hold it. A copy is removed when
/// the count reaches zero, so its lifetime belongs to the copy rather than to a caller remembering
/// to clean up.
fn localized_registry() -> &'static std::sync::Mutex<HashMap<PathBuf, usize>> {
    static REG: std::sync::OnceLock<std::sync::Mutex<HashMap<PathBuf, usize>>> = std::sync::OnceLock::new();
    REG.get_or_init(|| std::sync::Mutex::new(HashMap::new()))
}

/// The lock serializing cache copies for one destination path — see [`App::localize`], which holds
/// it across the copy so a second caller waits for the first rather than duplicating it.
///
/// Async, because it is held across the copy's `await`. Entries are never removed: one small entry
/// per distinct alignment localized in this process, bounded by the workspace's alignment count,
/// which is cheaper than the bookkeeping needed to retire them safely.
fn copy_gate(local: &Path) -> std::sync::Arc<tokio::sync::Mutex<()>> {
    #[allow(clippy::type_complexity)]
    static GATES: std::sync::OnceLock<std::sync::Mutex<HashMap<PathBuf, std::sync::Arc<tokio::sync::Mutex<()>>>>> =
        std::sync::OnceLock::new();
    let gates = GATES.get_or_init(|| std::sync::Mutex::new(HashMap::new()));
    let mut gates = gates.lock().unwrap();
    std::sync::Arc::clone(gates.entry(local.to_path_buf()).or_default())
}

/// An alignment path to read from, owning any local copy made for it.
///
/// The previous design cached copies in a directory cleared by one caller — `analyze_biosample` —
/// while three call sites created them. Every other path (notably the Y genotyping a batch drives)
/// copied ~400 MB per alignment and never cleaned up; that reached 687 files and 145 GB, and filled
/// the volume mid-run. Tying removal to `Drop` makes that leak unrepresentable, and the refcount
/// keeps the original benefit: a subject's several passes still share one copy instead of re-copying
/// per pass.
pub(crate) struct LocalAlignment {
    path: PathBuf,
    /// False when `path` is the original (no copy was made, so nothing to remove).
    owned: bool,
}

impl LocalAlignment {
    fn borrowed(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            owned: false,
        }
    }

    fn owned(path: PathBuf) -> Self {
        Self { path, owned: true }
    }

    /// Register interest in an existing copy; `true` when one was present and is now retained.
    ///
    /// A file with no registry entry is a **leftover from an earlier process** — `Drop` never ran,
    /// so the run was killed — and is validated against `expect_len` (the remote's size) before
    /// being trusted, then discarded if it doesn't match. Adopting a leftover on its existence
    /// alone is how a truncated copy gets read as though it were the alignment: the failure then
    /// appears as a decode error deep into a walk, pointing at a cache path whose file is deleted
    /// moments later. A wrong-sized copy is worth exactly one re-copy to be rid of.
    ///
    /// An entry that *is* registered belongs to a live holder in this process and was validated
    /// when it was made, so it is shared without re-statting.
    fn retain(local: &Path, expect_len: Option<u64>) -> bool {
        let mut reg = localized_registry().lock().unwrap();
        if let Some(n) = reg.get_mut(local) {
            *n += 1;
            return true;
        }
        let Ok(md) = std::fs::metadata(local) else {
            return false;
        };
        if !md.is_file() {
            return false;
        }
        if let Some(want) = expect_len {
            if md.len() != want {
                eprintln!(
                    "localize: discarding a stale cache copy at {} ({} bytes, expected {want})",
                    local.display(),
                    md.len()
                );
                let _ = std::fs::remove_file(local);
                return false;
            }
        }
        reg.insert(local.to_path_buf(), 1);
        true
    }

    /// The path to read from — the local copy when one was made, else the original.
    pub(crate) fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for LocalAlignment {
    fn drop(&mut self) {
        if !self.owned {
            return;
        }
        let mut reg = localized_registry().lock().unwrap();
        let remaining = match reg.get_mut(&self.path) {
            Some(n) => {
                *n = n.saturating_sub(1);
                *n
            }
            None => 0,
        };
        if remaining > 0 {
            return;
        }
        reg.remove(&self.path);
        // Best-effort: a copy left behind is a disk-space problem, not a correctness one, and a
        // panic here would mask whatever the caller was actually doing.
        let _ = std::fs::remove_file(&self.path);
        let p = self.path.to_string_lossy().into_owned();
        for suffix in [".crai", ".bai"] {
            let _ = std::fs::remove_file(PathBuf::from(format!("{p}{suffix}")));
        }
        let _ = std::fs::remove_file(self.path.with_extension("bai"));
    }
}

#[cfg(test)]
mod local_alignment_tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nav-localaln-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn a_local_copy_is_removed_with_its_index_when_dropped() {
        let d = scratch("drop");
        let (cram, crai) = (d.join("a.cram"), d.join("a.cram.crai"));
        std::fs::write(&cram, "x").unwrap();
        std::fs::write(&crai, "i").unwrap();

        assert!(LocalAlignment::retain(&cram, None), "an existing copy registers");
        drop(LocalAlignment::owned(cram.clone()));
        assert!(!cram.is_file(), "the copy is removed when the last holder drops");
        assert!(!crai.is_file(), "and so is its index");
    }

    #[test]
    fn a_shared_copy_survives_until_the_last_holder_drops() {
        // The reason for the refcount: a subject's passes each localize the same alignment, and
        // removing it when the first finishes would force the rest to re-copy ~400 MB.
        let d = scratch("shared");
        let cram = d.join("b.cram");
        std::fs::write(&cram, "x").unwrap();

        assert!(LocalAlignment::retain(&cram, None));
        let first = LocalAlignment::owned(cram.clone());
        assert!(LocalAlignment::retain(&cram, None), "a second holder shares the copy");
        let second = LocalAlignment::owned(cram.clone());

        drop(first);
        assert!(cram.is_file(), "still held by the second");
        drop(second);
        assert!(!cram.is_file(), "removed once nobody holds it");
    }

    /// A leftover from a killed run must be checked, not trusted. Adopting a short copy is how a
    /// truncated cache entry gets read as though it were the alignment — surfacing as a decode
    /// failure tens of gigabytes into a walk, blamed on a cache file that is deleted moments later.
    #[test]
    fn a_wrong_sized_leftover_is_discarded_rather_than_adopted() {
        let d = scratch("stale");
        let cram = d.join("d.cram");
        std::fs::write(&cram, "truncated").unwrap();

        assert!(
            !LocalAlignment::retain(&cram, Some(9_999)),
            "a copy that doesn't match the remote's size must not be adopted"
        );
        assert!(!cram.is_file(), "and it is removed, so the next attempt re-copies");
    }

    #[test]
    fn a_correctly_sized_leftover_is_adopted() {
        let d = scratch("stale-ok");
        let cram = d.join("e.cram");
        std::fs::write(&cram, "123456789").unwrap();

        assert!(LocalAlignment::retain(&cram, Some(9)), "a matching copy is reusable");
        drop(LocalAlignment::owned(cram.clone()));
    }

    /// Two concurrent copiers must never be able to name the same scratch file: sharing one let
    /// them open a single inode, where one truncates what the other is writing — and keeps writing
    /// into it after the other renames it into place and starts reading.
    #[test]
    fn each_copy_gets_its_own_partial_path() {
        let local = PathBuf::from("/tmp/nav-cache/abc.cram");
        let (a, b) = (partial_path(&local), partial_path(&local));
        assert_ne!(a, b, "two callers must not share a scratch path");
        assert_eq!(a.parent(), local.parent(), "the scratch stays beside its destination");
        for p in [&a, &b] {
            let name = p.file_name().unwrap().to_string_lossy().into_owned();
            assert!(name.starts_with("abc.cram.partial."), "{name}");
        }
    }

    /// A copy that arrives short is rejected instead of published, and leaves nothing behind — not
    /// the scratch file, and not the index copied ahead of it. Both used to accumulate in the cache.
    #[test]
    fn a_short_copy_is_rejected_and_leaves_no_debris() {
        let d = scratch("short");
        let (remote, local) = (d.join("src.cram"), d.join("dst.cram"));
        std::fs::write(&remote, "0123456789").unwrap();
        std::fs::write(d.join("src.cram.crai"), "idx").unwrap();

        // Claim the remote is larger than it is — the same shape as a copy cut short.
        let err = copy_with_index(&remote, &local, Some(64)).expect_err("a short copy must fail");
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof, "{err}");

        assert!(!local.is_file(), "a rejected copy is not published");
        assert!(!d.join("dst.cram.crai").is_file(), "nor is the index that preceded it");
        let debris: Vec<_> = std::fs::read_dir(&d)
            .unwrap()
            .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
            .filter(|n| n.contains(".partial"))
            .collect();
        assert!(debris.is_empty(), "scratch files left behind: {debris:?}");
    }

    /// The gate is what stops two overlapping callers from both copying the same alignment. Every
    /// worker command is `tokio::spawn`ed, so that overlap is real, and the duplicate was a second
    /// full pull of a 40 GB CRAM over the network whose only outcome was a failed rename.
    #[tokio::test]
    async fn one_destination_copies_at_a_time_and_others_are_not_blocked() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let path = PathBuf::from("/tmp/nav-cache/gate.cram");
        let (inside, max_seen) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
        let mut tasks = Vec::new();
        for _ in 0..4 {
            let (path, inside, max_seen) = (path.clone(), inside.clone(), max_seen.clone());
            tasks.push(tokio::spawn(async move {
                let gate = copy_gate(&path);
                let _held = gate.lock().await;
                let now = inside.fetch_add(1, Ordering::SeqCst) + 1;
                max_seen.fetch_max(now, Ordering::SeqCst);
                tokio::time::sleep(std::time::Duration::from_millis(20)).await;
                inside.fetch_sub(1, Ordering::SeqCst);
            }));
        }
        for t in tasks {
            t.await.unwrap();
        }
        assert_eq!(
            max_seen.load(Ordering::SeqCst),
            1,
            "two callers were inside the copy for one destination at once"
        );

        // A different alignment must not queue behind it, or a batch would serialize on the cache.
        let first = copy_gate(&path);
        let _held = first.lock().await;
        let other = copy_gate(&PathBuf::from("/tmp/nav-cache/other.cram"));
        assert!(
            other.try_lock().is_ok(),
            "an unrelated destination must not wait on this one"
        );
    }

    #[test]
    fn a_matching_copy_is_published_with_its_index() {
        let d = scratch("good");
        let (remote, local) = (d.join("src.cram"), d.join("dst.cram"));
        std::fs::write(&remote, "0123456789").unwrap();
        std::fs::write(d.join("src.cram.crai"), "idx").unwrap();

        copy_with_index(&remote, &local, Some(10)).expect("a whole copy succeeds");
        assert_eq!(std::fs::read(&local).unwrap(), b"0123456789");
        assert!(d.join("dst.cram.crai").is_file(), "the index lands beside it");
    }

    #[test]
    fn the_original_is_never_removed() {
        // A path we didn't copy (local disk, or a failed copy falling back to the remote) must be
        // left alone — deleting the user's own alignment would be catastrophic.
        let d = scratch("borrowed");
        let original = d.join("c.cram");
        std::fs::write(&original, "x").unwrap();

        let guard = LocalAlignment::borrowed(&original);
        assert_eq!(guard.path(), original);
        drop(guard);
        assert!(original.is_file(), "the source file must survive");
    }
}
