//! `impl App` methods extracted from `lib.rs` (the `analysis` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;
use navigator_analysis::{contig, CancelToken};

impl App {
    // ---- analysis (compute + persist) --------------------------------------

    /// Run the coverage walker and the callable walker on the BAM file of an alignment. The method
    /// writes the result as a `coverage` artifact with a version.
    ///
    /// The I/O of noodles blocks, so it runs on its own thread. The async runtime then continues.
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

    /// Calculate the coverage with the BAM path and the reference path of the alignment, and then
    /// write the result. The method fails when the store holds no such alignment, and when that
    /// alignment holds no path.
    pub async fn run_coverage_for_alignment(&self, alignment_id: i64) -> Result<CoverageResult, AppError> {
        self.run_coverage_for_alignment_with_progress(alignment_id, |_, _| {}, CancelToken::none())
            .await
    }

    /// The same work as [`run_coverage_for_alignment`], with a progress report. The method calls
    /// `progress(contigs_done, contigs_total)` as the whole-genome pass reads each contig.
    ///
    /// That pass is the slow step, and it needs some minutes on a real WGS BAM file. So a progress
    /// bar can move, and the app does not look stopped. The callback runs on the thread that
    /// blocks.
    pub async fn run_coverage_for_alignment_with_progress(
        &self,
        alignment_id: i64,
        mut progress: impl FnMut(usize, usize) + Send + 'static,
        cancel: CancelToken,
    ) -> Result<CoverageResult, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = Self::alignment_file(&aln)?;
        // The import step does not ask the user for the reference. When the alignment holds no
        // FASTA path, the gateway finds the build of that alignment. It reads the cache, and it
        // downloads the file when the cache holds none.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        // For a targeted test, such as Big Y, read only the target chromosomes. The depth that the
        // app reports then describes the target. Across the full genome, most contigs hold no read,
        // and they make that value almost zero.
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

    /// The list of contigs that the coverage walk reads, for a targeted test. The method returns
    /// `None` for a WGS test and an autosomal test, and the walk then reads the full genome.
    ///
    /// A Y test, such as FTDNA Big Y or Y Elite, reads chrY. It also reads chrM, so the signal "this
    /// test holds mtDNA reads" survives. A few Big Y files hold mitochondrial reads, and the UI
    /// hides the mtDNA sections when chrM holds none.
    ///
    /// An mtDNA test reads chrM only.
    ///
    /// The list does not depend on the build. It holds each contig name with the `chr` prefix and
    /// each name without it.
    async fn coverage_target_allowlist(&self, alignment_id: i64) -> Result<Option<HashSet<String>>, AppError> {
        use navigator_domain::testtype::TargetType;
        let aln = self.alignment_or_err(alignment_id).await?;
        let Some(run) = sequence_run::get(self.store.pool(), aln.sequence_run_id).await? else {
            return Ok(None);
        };
        // The code calls `target_of` and not `by_code`. So a label that a person wrote, such as
        // "Big Y", still limits the walk to chrY and chrM. A bulk import writes such a label, and
        // the `--test-type` option also writes one, in place of BIG_Y_500 or BIG_Y_700.
        //
        // Without this call, the walk reads the full genome. On a targeted CRAM file with many
        // references, that walk is the stop of about one hour in a batch analysis.
        let contigs: &[&str] = match navigator_domain::testtype::target_of(&run.test_type) {
            Some(TargetType::YChromosome) => &["chrY", "Y", "chrM", "chrMT", "M", "MT"],
            Some(TargetType::MtDna) => &["chrM", "chrMT", "M", "MT"],
            _ => return Ok(None),
        };
        Ok(Some(contigs.iter().map(|s| s.to_string()).collect()))
    }

    /// Shows whether a cached coverage result covers the correct contigs for the test of this
    /// alignment.
    ///
    /// A targeted test, such as Big Y or mtFull, must cover its target contigs only. A cached
    /// whole-genome result for such a test is wrong, because the depth is small across the contigs
    /// with no read. The app must calculate that result again.
    ///
    /// A whole-genome test has no list of contigs, and its result is always correct.
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

    /// The cached coverage result, for a later analysis. The method returns the stored result only
    /// when that result covers the correct contigs for the test. See
    /// [`Self::coverage_is_correctly_scoped`].
    ///
    /// A whole-genome result for a targeted test is wrong. The method then returns nothing, and the
    /// caller calculates the correct result.
    pub async fn cached_coverage_for_analysis(&self, alignment_id: i64) -> Result<Option<CoverageResult>, AppError> {
        match self.cached_coverage(alignment_id).await? {
            Some(cov) if self.coverage_is_correctly_scoped(alignment_id, &cov).await? => Ok(Some(cov)),
            _ => Ok(None),
        }
    }

    /// Find the biological sex from the ratio between the read density of chrX and the read density
    /// of the autosomes. The method writes the result as a `sex` artifact.
    ///
    /// The step is fast, because a BAM file has a BAI index. The code uses `reference` only to
    /// decode a CRAM file.
    pub async fn run_sex(&self, alignment_id: i64) -> Result<navigator_analysis::sex::SexInferenceResult, AppError> {
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let result =
            tokio::task::spawn_blocking(move || navigator_analysis::sex::infer_from_bam(&bam, reference.as_deref()))
                .await??;
        self.save_analysis(alignment_id, "sex", "1", &result).await?;
        self.write_back_inferred_sex(alignment_id, &result).await?;
        Ok(result)
    }

    /// Write the sex that the code found to the biosample, when the user gave none. The subjects
    /// table and the header then show that value in place of "Unknown".
    ///
    /// The method does nothing when the code found no sex, and when the biosample already holds
    /// one.
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

    /// Copy the library-level read statistics of an alignment to its sequence run. Those values are
    /// `total_reads`, `mean_read_length`, and `mean_insert_size`. The run card of the Data Sources
    /// tab then shows them, and the app reads no file again.
    ///
    /// The step is optional, and the method ignores an absent alignment and an absent run.
    ///
    /// When a run holds more than one alignment, the last write wins. These values describe the
    /// library, so each alignment gives the same answer.
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

    /// The scratch directory for an alignment that the code copied from a slow volume or a
    /// removable volume. See [`localize`]. A [`LocalAlignment`] value owns each entry, and the code
    /// removes that entry after the last holder drops it.
    pub(crate) fn align_cache_dir() -> std::path::PathBuf {
        navigator_refgenome::cache::base_dir().join("cache").join("aln")
    }

    /// Copy `remote` to the local cache and return the *local* path, when that file sits on a slow
    /// volume or a removable volume. Such a volume has a `/Volumes/…` mount point. The method also
    /// copies the `.crai` index or the `.bai` index. For any other path, it returns `remote` with no
    /// change.
    ///
    /// An analysis walker reads records at random positions. It seeks to a region, and it decodes
    /// each read. That access is very slow over a network mount or a USB mount. A plain sequential
    /// **copy** of the same file is fast.
    ///
    /// So the code pays for one fast copy first, and each later pass reads from the local disk. The
    /// passes of one subject share that copy, and [`clear_align_cache`] removes it for each subject.
    ///
    /// A failed copy gives the remote path. The analysis is then slow, and it still works.
    pub(crate) async fn localize(&self, remote: &Path) -> LocalAlignment {
        if std::env::var_os("NAVIGATOR_NO_LOCALIZE").is_some() || !is_removable_volume(remote) {
            return LocalAlignment::borrowed(remote);
        }
        let local = Self::align_cache_dir().join(local_cache_name(remote));

        // Order the two steps, the cache test and the copy, for each destination. The worker calls
        // `tokio::spawn` for each command. So a batch walk and a command for one alignment do
        // overlap on the same alignment.
        //
        // Without this lock, both find no cache entry and both copy. The second copy reads another
        // 40 GB over the network for no result. The task that loses the rename then reads the
        // remote file.
        //
        // The code holds the lock across the copy. So it must not take that lock a second time. No
        // call of `localize` can run while another call is open *on the same path in the same
        // task*.
        //
        // The three call sites run in sequence today. `debug_y_calls` waits for `base_calls` to
        // complete before it localizes its own file. Keep that order.
        let gate = copy_gate(&local);
        let _copying = gate.lock().await;

        // Read the size of the remote file one time. That value answers two questions. It shows
        // whether the code can trust a copy that exists, and it shows whether the new copy is
        // complete.
        let remote_len = tokio::fs::metadata(remote).await.ok().map(|m| m.len());

        // Another holder already uses this copy. Share it, and add one to the count.
        if LocalAlignment::retain(&local, remote_len) {
            return LocalAlignment::owned(local);
        }
        let (remote_owned, local2) = (remote.to_path_buf(), local.clone());
        match tokio::task::spawn_blocking(move || copy_with_index(&remote_owned, &local2, remote_len)).await {
            // This step can still fail, because another holder can remove the copy in the time
            // between the two steps. That holder completes its work and drops the count to zero.
            //
            // An `owned` handle to an absent path fails the walk with an ENOENT error that no user
            // can read. The code then also tries to remove a file that is not there.
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

    /// Run the unified quality-metrics walker. It makes **one pass** over the BAM file or the CRAM
    /// file of the alignment. That pass gives three results: the coverage with the callable regions,
    /// the quality metrics of each read, and the sex.
    ///
    /// The separate calls `run_coverage`, `run_read_metrics`, and `run_sex` cost more. They read a
    /// BAM file two times, and a CRAM file three times.
    ///
    /// The method writes each of the three results under its existing artifact key. Those keys are
    /// `coverage` with `COVERAGE_VERSION`, `read_metrics` with `"1"`, and `sex` with `"1"`.
    ///
    /// So `cached_coverage`, `cached_read_metrics`, `cached_sex`, and the reuse rule of the SV step
    /// each work with no change.
    pub async fn run_unified_metrics(&self, alignment_id: i64) -> Result<UnifiedMetricsResult, AppError> {
        self.run_unified_metrics_with_progress(alignment_id, |_, _| {}, CancelToken::none())
            .await
    }

    /// The same work as [`run_unified_metrics`], with a progress report. The method calls
    /// `progress(contigs_done, contigs_total)` as the whole-genome coverage step completes each
    /// contig. That step is the slow one.
    ///
    /// The method uses the parallel walker, which works on each contig at the same time. For a CRAM
    /// file, and for a BAM file with no index, it reads the file from start to end instead.
    ///
    /// The callback is `Fn + Sync`, because the worker threads of the parallel walker call it at the
    /// same time.
    pub async fn run_unified_metrics_with_progress(
        &self,
        alignment_id: i64,
        progress: impl Fn(usize, usize) + Send + Sync + 'static,
        cancel: CancelToken,
    ) -> Result<UnifiedMetricsResult, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let run_id = aln.sequence_run_id;
        // Copy the file from a slow volume or a removable volume to the local disk first. The
        // walker reads records at random positions, and that access is much slower over a network
        // mount or a USB mount than one bulk copy.
        //
        // The code holds this value for the full walk. A drop of it removes the local copy.
        let bam = self.localize(&Self::alignment_file(&aln)?).await;
        let bam = bam.path().to_path_buf();
        // The walker needs a reference. It decodes the CRAM file with that reference, and it finds
        // each N base of the reference. When the import stored no FASTA path, the gateway finds the
        // build.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        // Limit a targeted test, such as Big Y or mtFull, to its target contigs. The separate
        // coverage walker uses the same rule.
        //
        // Across the full genome, most contigs hold no read, and they make the depth small. A Big Y
        // test then reads as about 0.2x, and its true depth on chrY is about 50x. A WGS test keeps
        // the whole-genome walk.
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
        // The sex. A Y test, such as Big Y or Y Elite, reads the Y chromosome of the donor. So that
        // donor is male, by definition.
        //
        // A walk of chrY alone holds no ratio between chrX and the autosomes, and the code needs
        // that ratio. The ratio is also wrong across the full genome. In a Big Y file, chrX and the
        // autosomes each hold about 0.4x, and the code then reads the donor as *female*.
        //
        // So the code writes Male for a Y test. That value replaces the result of the ratio, and it
        // replaces a value from an earlier run. A WGS test and an mt test keep the result of the
        // walk.
        let y_targeted = matches!(
            sequence_run::get(self.store.pool(), run_id)
                .await?
                .as_ref()
                .and_then(|r| navigator_domain::testtype::target_of(&r.test_type)),
            Some(navigator_domain::testtype::TargetType::YChromosome)
        );
        // An alignment with reads on chrY only is male, as a Y test is. Its chrY contig holds
        // almost each read, and its autosomes hold a few reads that the mapper placed wrongly.
        //
        // Two files have that shape. One is a chrY extract, such as GRCh38 chrY reads that the app
        // realigned to hs1. The other is a Y Elite capture, or a Big Y capture, that arrived with a
        // WGS label.
        //
        // The ratio can read such a file as *female*. That value stops the full Y pipeline with no
        // message, because `assign_y_haplogroup` skips a female subject before it reads the tree.
        //
        // So the code finds this shape from the read count of each contig, and it writes Male, as it
        // does for a Y test.
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
                // This value is definite: a Y test, or an alignment with reads on chrY only, is
                // male. So the code replaces a sex from an earlier run, and that set holds a wrong
                // "Female" value. It does not only write into an empty field.
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
        // A new SV result in the cache, from a source file that did not change, is correct. The
        // code uses that result and calculates nothing.
        if let Some(c) = self.cached_sv(alignment_id).await? {
            return Ok(c);
        }
        let aln = self.alignment_or_err(alignment_id).await?;
        let reference_build = aln.reference_build.clone();
        // Find the reference for the decoder. See alignment_reference_for_decode. A CRAM file needs
        // it, and a BAM file uses None.
        //
        // The SV step reads no reference *base*. But a decode of a CRAM record does read one. So the
        // walker also needs the reference, and not only the step that reads the contig lengths from
        // the header.
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

    /// Genotype the short tandem repeats on `contig` from the alignment. The caller reads each
    /// record that covers a full tract, and it uses the HipSTR reference tracts. It calls chrY and
    /// chrM as haploid, and each other contig as diploid.
    ///
    /// The method writes the result as a `str:{contig}` artifact. So the cache holds it, and a
    /// change to the source file makes it invalid, as it does for another analysis.
    ///
    /// The method fails when no STR reference exists for the build of the alignment. The tracts
    /// belong to one build. CHM13 and GRCh37 each need their own reference, or a liftover, and no
    /// code does that work yet.
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
        // None for a BAM. STR region-genotyping reads the alignment; it does not consult reference bases.
        let (bam, reference) = self.alignment_reference_for_decode(alignment_id).await?;
        // A cell holds one copy of chrY and one copy of chrM, so each has one allele. It holds two
        // copies of each autosome, and a female cell holds two copies of chrX.
        //
        // So the code calls chrY and chrM as haploid, and each other contig as diploid. A rule for
        // chrX that reads the sex is a later improvement.
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

    /// Compare the STR markers from the sequence data with the vendor Y-STR profile that the user
    /// imported. The By-Panel view shows this comparison.
    ///
    /// The [`navigator_analysis::strmarker`] table changes each called value to the FTDNA
    /// convention. A corpus of real kits calibrated that table.
    ///
    /// The result holds one row for each marker in either source. A row holds the called value with
    /// its calibration state, the imported value, and a flag that shows whether the two agree.
    ///
    /// The `contig` value is usually `chrY`. The method reads the `str:{contig}` calls from the
    /// cache.
    pub async fn str_concordance(&self, alignment_id: i64, contig: String) -> Result<Vec<StrConcordanceRow>, AppError> {
        use navigator_analysis::strmarker::{called_markers_build, normalize_marker, MarkerStatus, StrBuild};

        // For a few markers, the offset of the FTDNA convention changes with the build. The CHM13
        // liftover moved the boundary of some tracts. So the code reads the offsets of the build of
        // this alignment.
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

    /// Select the best alignment of the subject for STR work, and compare the Y-STR markers on
    /// chrY. The UI calls this method.
    ///
    /// An alignment can do STR work when a HipSTR reference exists for its build. See
    /// [`str_reference_path`](Self::str_reference_path). Among those alignments, the one with the
    /// highest mean coverage wins.
    ///
    /// A CRAM file needs no stored reference here, because
    /// [`run_str_calls`](Self::run_str_calls) finds one for the decoder.
    ///
    /// The method fails with a hint when no alignment passes. The two causes are an absent HipSTR
    /// reference and a subject with no alignment.
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
            // A CRAM file with no stored reference is acceptable, because run_str_calls finds one
            // through the gateway.
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

    /// Call the **de-novo diploid** SNVs across the full `contig`. The caller writes a heterozygous
    /// call as 0/1 and a homozygous alternate call as 1/1. The cache key is the alignment with the
    /// contig.
    ///
    /// The method reads the BAM file of the alignment and its reference, which the code finds from
    /// the build. It returns the [`SiteGenotype`] values in the order of their positions. Give them
    /// to [`Self::diploid_vcf`].
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

    /// A diploid VCF file of the de-novo diploid SNV calls of `contig`. The file uses VCFv4.2, and
    /// its format field is `GT:AD:DP:GQ:PL`. The method calculates those calls and writes them to
    /// the cache when the cache holds none. The sample column is `aln<id>`.
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

    /// A **whole-genome** diploid VCF file. It holds the de-novo SNV calls and indel calls across
    /// the diploid primary chromosomes of the alignment, which are 1 to 22 and X. The cache holds
    /// the result of each contig.
    ///
    /// The file holds **no** chrY data and **no** chrM data. A cell holds one copy of each, so the
    /// diploid model, with its 0/1 calls, is wrong for them. Their variants come from the haploid
    /// caller, and from the Y and mt haplogroup features with the mtDNA mutation list.
    ///
    /// This method is a full WGS calling pass, and it costs much. The caller runs it away from the
    /// UI thread, on the export path.
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

    /// The alignments of the subject on the **most frequent reference build**. The code compares the
    /// canonical build, so `chm13v2` and `hs1` count as one build here.
    ///
    /// The consensus diploid genotype pools the alignments of one build only. The position of a
    /// de-novo variant does not compare across two builds, and a join by position needs a liftover
    /// of the full genome. That work is not in this feature.
    ///
    /// The method returns `None` when the subject has no alignment.
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

    /// The **consensus diploid genotype of a subject**, across its WGS alignments on one build. This
    /// value is the joint genotype, which is opportunity #3.
    ///
    /// [`reconcile_site_genotypes`] does the work in four steps.
    ///
    /// It calls the variants of each alignment, and [`run_diploid_calls`] gives those calls from the
    /// cache. It joins the SNV sites of each alignment into one set. It then genotypes **each**
    /// alignment at each site of that set. So a site that one run does not hold gets its real
    /// hom-ref call or no-call. It then votes a dosage of 0, 1, or 2 at each site, and a deeper run
    /// has more weight.
    ///
    /// The method returns the consensus sites with a variant, which are the heterozygous sites and
    /// the homozygous alternate sites. The `contigs` value limits the scan, and `None` reads each
    /// primary chromosome.
    ///
    /// The method costs much: one call pass and one forced-call pass for each alignment. The user
    /// starts it from the export screen, and the method stores nothing.
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

        // The pair (bam, reference) of each alignment on this build. The code finds each pair one
        // time.
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
                // The header of this alignment can hold no such contig, because the inputs differ.
                // Skip that contig for this source, and do not stop the full consensus.
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

        // 4. Vote at each site to get the consensus. The value min_depth = 2 means that a run
        // gives no vote only when it has almost no coverage there. A deeper run has more weight than
        // a shallow one.
        Ok(caller::reconcile_site_genotypes(&per_aln, 2))
    }

    /// A **consensus** diploid VCF file for the subject, in VCFv4.2. It holds the joint genotype
    /// across the alignments on one build. See [`consensus_diploid_calls`]. The sample column is
    /// `consensus`. The method costs much, and the export path runs it away from the UI thread.
    pub async fn consensus_diploid_vcf(&self, biosample_guid: SampleGuid) -> Result<String, AppError> {
        let calls = self
            .consensus_diploid_calls(biosample_guid, None, CancelToken::none())
            .await?;
        Ok(navigator_analysis::vcf::write_diploid_vcf("consensus", &calls))
    }

    /// Call the de-novo variants on `contig` with the stored paths of the alignment.
    ///
    /// The method returns the BAM path of the alignment and a reference FASTA path that the code can
    /// use. That reference is the stored path. When the alignment holds none, the gateway finds one
    /// from the build of the alignment, from the cache or by a download.
    ///
    /// The method fails only when the alignment holds no BAM path. Use it in a step that *needs* the
    /// reference. The user then supplies no reference, because the build in the header gives it.
    pub(crate) async fn alignment_bam_reference(&self, alignment_id: i64) -> Result<(PathBuf, PathBuf), AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = Self::alignment_file(&aln)?;
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

    /// The path of the alignment, and a reference that the code can use to **decode** it.
    ///
    /// No reader can open a CRAM file without its reference. So for a CRAM file the method takes
    /// the stored path first, and then the build through the gateway. It reads the cache before it
    /// starts a download.
    ///
    /// A reader can open a BAM file with no reference. So for a BAM file the method returns the
    /// stored path with no change, and that value is usually `None`. It never starts a download.
    ///
    /// Use this method to read records, to read a pileup, and to genotype a SNP site. None of those
    /// steps reads a reference base.
    ///
    /// Use [`alignment_bam_reference`](Self::alignment_bam_reference) for a caller path, such as a
    /// de-novo SNV call or indel call. Those paths need the reference for a BAM file also.
    pub(crate) async fn alignment_reference_for_decode(
        &self,
        alignment_id: i64,
    ) -> Result<(PathBuf, Option<PathBuf>), AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = Self::alignment_file(&aln)?;
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

    /// Each distinct reference build across the alignments of a subject. An analysis of that subject
    /// can need the FASTA file of any of them.
    ///
    /// The code reads this list after an import, and before an analysis of the subject. It then
    /// downloads each file with a progress bar. So a download during the analysis never surprises
    /// the user.
    pub async fn reference_builds_for_subject(&self, biosample_guid: SampleGuid) -> Result<Vec<String>, AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut builds: Vec<String> = alns.into_iter().map(|a| a.reference_build).collect();
        builds.sort();
        builds.dedup();
        Ok(builds)
    }

    /// The reference build of one alignment. The method returns `None` when the store holds no such
    /// alignment. The code reads this value to find the reference before it analyzes that
    /// alignment.
    pub async fn reference_build_of_alignment(&self, alignment_id: i64) -> Result<Option<String>, AppError> {
        Ok(alignment::get(self.store.pool(), alignment_id)
            .await?
            .map(|a| a.reference_build))
    }

    /// The id of each alignment of a subject that has a BAM file or a CRAM file. The code reads this
    /// list to make the coordinate index of each one, with a progress bar. It does that work after an
    /// import, and before an analysis of the subject.
    pub async fn alignment_ids_for_subject(&self, biosample_guid: SampleGuid) -> Result<Vec<i64>, AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        Ok(alns
            .into_iter()
            .filter(|a| a.bam_path.is_some())
            .map(|a| a.id)
            .collect())
    }

    /// Make sure that the coordinate index of the alignment exists. That index is a `.bai` file or a
    /// `.crai` file, and the method **makes it** when the disk holds none.
    ///
    /// Each analysis that queries a region needs that index. Those analyses are the walker that
    /// works on one contig, the step that finds the callable intervals, the de-novo caller, and the
    /// STR caller. Without an index, such a step fails, or it reads the full file from start to end.
    ///
    /// The method returns the path of the index when it made one. It returns `None` when the index
    /// already existed.
    ///
    /// The method calls `progress(done, total)`. For a BAM file, that call gives a fraction of the
    /// bytes. For a CRAM file, the `total` value is `None`, and the progress has no end value.
    ///
    /// The method reads the file one time, from start to end, on a thread that can decode
    /// safely.
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

    /// Report the reason that the app can not read an alignment. The report names the **exact
    /// file** at fault, and not the file that the failed call received.
    ///
    /// [`navigator_analysis::preflight`] gives the reason for that rule. A `.crai` file that the app
    /// can not read, and a CRAM file that it can not read, give the same `io error on …cram` message
    /// today. On macOS, only the raw errno separates a privacy denial from TCC and a Unix permission
    /// denial.
    ///
    /// For the reference, the method reads the **cache only**, by design. A diagnostic must describe
    /// the machine as it is. A download of an absent FASTA file here hides the exact state that the
    /// user asked about. A CRAM file with no reference in the cache is a result of this check. It
    /// is not a task for this method.
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

    /// The [`PublishGate`] of an alignment, for its mean read length. A HiFi read needs fewer reads
    /// at a site than a short read. See [`PublishGate::for_read_len`].
    ///
    /// The method reads the first records of the BAM file. After any error, it returns the default
    /// gate for a short read.
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

/// True for a path on a removable mount or a network mount. On macOS such a mount is below
/// `/Volumes/…`.
///
/// A read of one record at a random position is slow on that mount, and a bulk copy from start to
/// end is fast. [`App::localize`] copies such a file to the local disk.
fn is_removable_volume(p: &Path) -> bool {
    p.starts_with("/Volumes/")
}

/// A local file name for a remote alignment. Two such names are never the same.
///
/// The file of each kit has the name `chrYM.cram`, so the base name alone gives the same local name
/// for two kits. The function hashes the full remote path and keeps the extension. So the reader
/// still finds the index beside the file, at `<local>.crai` or `<local>.bai`.
fn local_cache_name(remote: &Path) -> String {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    remote.to_string_lossy().hash(&mut h);
    let ext = remote.extension().and_then(|e| e.to_str()).unwrap_or("bam");
    format!("{:016x}.{ext}", h.finish())
}

/// A scratch name for the copy of `local` that one caller is writing. Each call gives a **different
/// name**.
///
/// A name from the destination alone, such as `<dest>.partial`, is the same name for each caller
/// that copies that alignment at the same time. Two callers then open the same inode.
///
/// One caller empties the file that the other one writes. It also continues to write into that file
/// after the other caller renames it into place and starts to read it.
///
/// A different name for each call makes that state impossible. It does not only make it rare.
fn partial_path(local: &Path) -> PathBuf {
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let stem = local.file_name().unwrap_or_default().to_string_lossy().into_owned();
    local.with_file_name(format!("{stem}.partial.{}.{n}", std::process::id()))
}

/// Copy `remote` to `local`, and copy the index beside it.
///
/// The function copies the index **first** and the main file last. For the main file it writes a
/// temporary file and then renames it. So a `local` file that exists always has its index. The cache
/// test in [`App::localize`] can then never see one file of the pair.
///
/// `expect_len` is the size of the remote file. When the caller gives that value, the function
/// refuses a copy with a different size and publishes nothing.
///
/// A short copy looks the same as damaged data. It gives a decode error tens of GB into a walk, such
/// as "unexpected end of file" or a bad container checksum. That error names the cache path, and the
/// code deletes the copy at the drop, before a user can look at it.
///
/// No partial file survives a failure. That rule covers the temporary file and the index that the
/// function copied first. An old `.partial` file stayed in the cache for all time, and it filled the
/// disk with data that no code used.
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

/// One step of a full analysis of one alignment, in the order that [`App::plan_full_analysis`]
/// gives. Each variant carries the values that its step needs. So a caller must handle each variant,
/// and it can not run a step that the plan left out.
#[derive(Debug, Clone, PartialEq)]
pub enum AnalysisStep {
    /// Coverage + callable, read-level QC, and sex inference in one pass over the alignment.
    QualityMetrics,
    /// The CNV calls and the discordant pairs. The step needs a depth of 10x or more, and it
    /// reports a depth below that value.
    ///
    /// **The user must select this step.** It is experimental. It is also the one step that reads
    /// each record in the file for a result that no other step uses. It needs hours for one
    /// whole-genome sample.
    ///
    /// The plan holds this step only when a caller sets `include_sv`. The "Call SV" button of the
    /// GUI and the `analyze --sv` command do that. No code runs it without a user.
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
    /// Estimate the admixture and the PCA values *from* the autosomal profile. This step must come
    /// after [`Self::AutosomalProfile`].
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
    /// The steps of a full analysis of `alignment_id`, in their order.
    ///
    /// This function is the **one** definition of that pipeline. The GUI sends progress events, and
    /// the CLI writes a log. So the two report a step in different ways.
    ///
    /// But the set of steps, and each condition that skips one, must be the same. The two did become
    /// different. The copy in the CLI genotyped the Y chromosome at each run, and it replaced a
    /// trusted external call. On ancient DNA it replaced that call with a worse one.
    ///
    /// The `coverage` value is the result that [`AnalysisStep::QualityMetrics`] calculated, when
    /// that step already ran. The mitochondrial decision then has the correct data. Before that
    /// step, pass `None`. The function then reads the cached coverage. With no cached value, it
    /// assumes that the file holds chrM reads.
    ///
    /// A caller that shows a count of steps must call this function again after the metrics step,
    /// because that count can become smaller.
    ///
    /// `include_sv` adds the experimental step [`AnalysisStep::StructuralVariants`]. That variant
    /// gives the reason for its default. `include_ancestry` adds the two ancestry steps that cost
    /// the most.
    pub async fn plan_full_analysis(
        &self,
        alignment_id: i64,
        include_ancestry: bool,
        include_sv: bool,
        coverage: Option<&CoverageResult>,
    ) -> Result<Vec<AnalysisStep>, AppError> {
        // Leave out the mitochondrial steps when the alignment holds no chrM read. An FTDNA Big Y
        // file is one example. A placement with no chrM data writes the RSRS root, and that result
        // has no value.
        //
        // With no coverage result, the plan keeps those steps. The code must not remove a step with
        // no message.
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
        // A step at the level of a subject needs that subject. The plan leaves out each such step
        // for an alignment with no subject.
        let guid = self.biosample_of_alignment(alignment_id).await.ok();

        let mut steps = vec![AnalysisStep::QualityMetrics];
        if include_sv {
            steps.push(AnalysisStep::StructuralVariants);
        }
        if has_mtdna {
            steps.push(AnalysisStep::MitoDenovo { contig: "chrM".into() });
        }
        // Leave out the internal Y step and the internal mt step under two conditions. A trusted
        // external caller already placed this alignment, and the user prefers that caller. Such a
        // caller is a GATK4 GVCF file through the sidecar fast path.
        //
        // A second walk gives a call that loses the vote. On ancient DNA it also gives a wrong call.
        //
        // Each `assign_*` command has the same guard. A plan without the step also saves the decode.
        // See external-caller-precedence.
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
        // The Y signature makes the descent report ready, and the user presses no button. That
        // button stays, and it rebuilds the report at any time.
        //
        // The code builds the signature one time, and it does not change a profile that exists. It
        // reads the file no more times, because the Y step above wrote the chrY genotypes to the
        // cache.
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

/// Shows whether a coverage result holds a read on the mitochondrial contig. The function accepts
/// both names of that contig.
fn chrm_has_reads(c: &CoverageResult) -> bool {
    c.contig_coverage_stats
        .iter()
        .any(|s| contig::is_chr_m(&s.contig) && s.num_reads > 0)
}

/// The local copies that exist now. The map takes a local path and gives the count of
/// [`LocalAlignment`] values that hold it.
///
/// The code removes a copy when its count reaches zero. So the copy itself controls its life, and no
/// caller must remember to remove it.
fn localized_registry() -> &'static std::sync::Mutex<HashMap<PathBuf, usize>> {
    static REG: std::sync::OnceLock<std::sync::Mutex<HashMap<PathBuf, usize>>> = std::sync::OnceLock::new();
    REG.get_or_init(|| std::sync::Mutex::new(HashMap::new()))
}

/// The lock that orders the cache copies for one destination path. [`App::localize`] holds it across
/// the copy, so a second caller waits for the first one and makes no second copy.
///
/// The lock is async, because the code holds it across the `await` of the copy.
///
/// The map never removes an entry. It holds one small entry for each alignment that this process
/// copies, and the count of alignments in the workspace limits that number. The code to remove an
/// entry safely costs more than those entries.
fn copy_gate(local: &Path) -> std::sync::Arc<tokio::sync::Mutex<()>> {
    #[allow(clippy::type_complexity)]
    static GATES: std::sync::OnceLock<std::sync::Mutex<HashMap<PathBuf, std::sync::Arc<tokio::sync::Mutex<()>>>>> =
        std::sync::OnceLock::new();
    let gates = GATES.get_or_init(|| std::sync::Mutex::new(HashMap::new()));
    let mut gates = gates.lock().unwrap();
    std::sync::Arc::clone(gates.entry(local.to_path_buf()).or_default())
}

/// A path to read an alignment from. This value owns the local copy of that alignment, when one
/// exists.
///
/// The earlier design kept each copy in a directory that one caller emptied, and that caller was
/// `analyze_biosample`. But three call sites made a copy.
///
/// Each other path copied about 400 MB for one alignment and removed nothing. The Y genotype step of
/// a batch is the main example. That fault reached 687 files and 145 GB, and it filled the volume
/// during a run.
///
/// A removal at the `Drop` call makes that fault impossible. The count of holders keeps the first
/// advantage. The passes of one subject still share one copy, and the code copies that file one
/// time.
pub(crate) struct LocalAlignment {
    path: PathBuf,
    /// False when `path` is the original file. The code made no copy, so it removes nothing.
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

    /// Take a share of a copy that exists. The method returns `true` when such a copy was there and
    /// now has one more holder.
    ///
    /// A file with no entry in the registry is a **file from an earlier process**. No `Drop` call
    /// ran there, so something stopped that run.
    ///
    /// The method compares such a file with `expect_len`, which is the size of the remote file. It
    /// removes the file when the two differ.
    ///
    /// A method that took such a file on its existence alone would read a short copy as the
    /// alignment. The fault then appears as a decode error deep in a walk. It names a cache path,
    /// and the code deletes that file a moment later. One more copy is a small price to remove a
    /// file of the wrong size.
    ///
    /// A file *with* an entry belongs to a holder in this process, and the code checked it when it
    /// made that copy. So the method shares it and reads no metadata.
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

    /// The path to read from. The value is the local copy when the code made one, and the original
    /// path when it made none.
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
        // The step is optional. A copy that stays on the disk costs space, and it gives no wrong
        // result. A panic here also hides the work of the caller.
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
        // The reason for the count of holders. Each pass of one subject localizes the same
        // alignment. A removal after the first pass makes each later pass copy about 400 MB
        // again.
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

    /// The code must check a file from a run that stopped. It must not trust that file.
    ///
    /// A method that took a short copy would read that copy as the alignment. The fault then appears
    /// as a decode failure tens of GB into a walk. The message names a cache file, and the code
    /// deletes that file a moment later.
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

    /// Two callers that copy at the same time must never use the same scratch file name.
    ///
    /// With one name, the two open the same inode. One caller then empties the file that the other
    /// one writes. It also continues to write into that file after the other caller renames it into
    /// place and starts to read it.
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

    /// The code refuses a copy that is too short, and it publishes nothing. It also leaves no file
    /// behind. That rule covers the scratch file and the index that it copied first. The cache used
    /// to collect both of them.
    #[test]
    fn a_short_copy_is_rejected_and_leaves_no_debris() {
        let d = scratch("short");
        let (remote, local) = (d.join("src.cram"), d.join("dst.cram"));
        std::fs::write(&remote, "0123456789").unwrap();
        std::fs::write(d.join("src.cram.crai"), "idx").unwrap();

        // Give a remote size that is larger than the real size. The result has the same shape as a
        // copy that stopped early.
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

    /// The lock stops two callers that overlap, so only one of them copies the alignment. Each
    /// worker command
    /// runs under `tokio::spawn`, so that overlap does occur. The second copy read a 40 GB CRAM file
    /// over the network again, and its only result was a failed rename.
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
        // The code must not change a path that it did not copy. Such a path is a file on the local
        // disk, or the remote file after a failed copy. A delete of the alignment of the user is a
        // very bad fault.
        let d = scratch("borrowed");
        let original = d.join("c.cram");
        std::fs::write(&original, "x").unwrap();

        let guard = LocalAlignment::borrowed(&original);
        assert_eq!(guard.path(), original);
        drop(guard);
        assert!(original.is_file(), "the source file must survive");
    }
}
