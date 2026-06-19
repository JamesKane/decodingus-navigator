//! `impl App` methods extracted from `lib.rs` (the `import_unified` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- unified import ----------------------------------------------------

    /// Detect a file's type and route it to the right subject importer (STR / variants /
    /// chip / mtDNA), using sensible defaults. Returns the detected type. Alignment files
    /// are rejected here — they attach to a sequencing test, not directly to a subject.
    /// Probe a BAM/CRAM header for the build/aligner/platform/test-type (best-effort).
    pub async fn probe_alignment(&self, path: PathBuf) -> Result<AlignmentProbe, AppError> {
        tokio::task::spawn_blocking(move || navigator_analysis::probe::probe_alignment(&path))
            .await?
            .map_err(AppError::from)
    }

    /// Scan a bounded prefix of an alignment's reads to infer the instrument/library identity —
    /// the `@RG SM/LB/PU` tags plus the most-frequent instrument/flowcell/platform from read names
    /// (the crowd-source input for resolving the lab). Off-thread (blocking IO + CRAM decode);
    /// `reference` is required for CRAM. Best-effort — callers tolerate an error.
    pub async fn library_stats(
        &self,
        path: PathBuf,
        reference: Option<PathBuf>,
    ) -> Result<navigator_analysis::library_stats::LibraryStats, AppError> {
        tokio::task::spawn_blocking(move || {
            navigator_analysis::library_stats::scan_library_stats(
                &path,
                reference.as_deref(),
                navigator_analysis::library_stats::DEFAULT_MAX_READS,
            )
        })
        .await?
        .map_err(AppError::from)
    }

    /// Auto-import an alignment file by probing its header: create the sequencing run (test type,
    /// platform, instrument) and the alignment (reference build + aligner) with no questions
    /// asked. The reference FASTA is **not** required — it's resolved from the build on demand;
    /// if already cached it's stored so every analysis step has it immediately.
    async fn import_alignment_file(&self, biosample_guid: SampleGuid, path: &Path) -> Result<(), AppError> {
        // Idempotent: skip if this exact BAM/CRAM is already recorded as an alignment.
        let path_str = path.to_string_lossy().into_owned();
        if alignment::list_all(self.store.pool())
            .await?
            .iter()
            .any(|a| a.bam_path.as_deref() == Some(path_str.as_str()))
        {
            return Ok(());
        }
        // Best-effort: a probe failure falls back to filename/defaults rather than aborting.
        let probe = self.probe_alignment(path.to_path_buf()).await.unwrap_or_default();

        // Resolve the reference first — the read-name scan needs it to decode a CRAM.
        let reference_build = probe
            .reference_build
            .clone()
            .unwrap_or_else(|| reference_build_for(path));
        // Store the cached reference path if we have it; otherwise leave it unset (resolved on
        // demand) — never block import on a download.
        let reference_path = self
            .gateway
            .cached_reference(&reference_build)
            .map(|p| p.to_string_lossy().into_owned());

        // Read-name scan → instrument/library identity (the lab crowd-source input). Best-effort:
        // it fills the platform/model the header `@RG` left blank, and the instrument/flowcell that
        // never live in the header. Skipped silently if the file can't be read (e.g. CRAM with no
        // resolved reference yet).
        let stats = self
            .library_stats(path.to_path_buf(), reference_path.as_deref().map(PathBuf::from))
            .await
            .ok();

        // Platform/model: prefer the header `@RG` (PL/PM); fall back to the read-name inference.
        let platform_name = probe
            .platform
            .clone()
            .or_else(|| {
                stats
                    .as_ref()
                    .and_then(|s| s.platform.clone())
                    .map(|p| p.to_uppercase())
            })
            .unwrap_or_else(|| "UNKNOWN".into());
        let instrument_model = probe
            .instrument_model
            .clone()
            .or_else(|| stats.as_ref().and_then(|s| s.instrument_model.clone()));

        // Test type: refine the header/platform guess with coverage *shape* from the BAI index —
        // a targeted-Y pile-up (autosomes empty) → Big Y / Y Elite / YSEQ; an mtDNA pile-up →
        // mtFull. Best-effort and cheap (O(contigs), no read scan); CRAM / unindexed BAMs have no
        // profile and keep the platform-based guess.
        let test_type = {
            let p = path.to_path_buf();
            let profile =
                tokio::task::spawn_blocking(move || navigator_analysis::testtype::coverage_profile_from_bai(&p, None))
                    .await
                    .ok()
                    .flatten();
            navigator_analysis::testtype::infer_test_type(
                profile.as_ref(),
                probe.platform.as_deref(),
                probe.vendor_hint.as_deref(),
                None,
            )
            .or_else(|| probe.test_type.clone())
            .unwrap_or_else(|| "WGS".into())
        };

        let run = self
            .record_sequence_run(NewSequenceRun {
                biosample_guid,
                platform_name,
                instrument_model,
                test_type,
                library_layout: stats.as_ref().and_then(|s| s.library_layout.clone()),
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await?;

        // Persist the inferred lab/instrument identity block (the crowd-source key). The lab
        // (`sequencing_facility`) stays unset — set manually, or resolved from `instrument_id`
        // once the AppView lookup ships (roadmap D8).
        if let Some(s) = &stats {
            let _ = sequence_run::set_library_stats(
                self.store.pool(),
                run.id,
                s.instrument_id.as_deref(),
                s.sample_name.as_deref(),
                s.library_id.as_deref(),
                s.platform_unit.as_deref(),
                s.flowcell_id.as_deref(),
            )
            .await;
            // Resolve the lab from the instrument id via the AppView (best-effort, cached).
            if let Some(inst) = s.instrument_id.as_deref() {
                if let Some(lab) = self.lookup_lab_by_instrument(inst).await {
                    let _ = sequence_run::set_facility(self.store.pool(), run.id, &lab).await;
                }
            }
        }

        // Defer the content hash (the file's identity, used to invalidate cached analyses): a
        // whole-file SHA-256 of a multi-GB alignment would block this import for minutes with no
        // feedback. Like the batch path, leave it `None` — `alignment_content_hash` computes and
        // caches it lazily on the first analysis that needs it.
        self.record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build,
            aligner: probe.aligner.clone().unwrap_or_else(|| "unknown".into()),
            variant_caller: None,
            bam_path: Some(path.to_string_lossy().into_owned()),
            reference_path,
            content_sha256: None,
        })
        .await?;
        Ok(())
    }

    pub async fn add_data(&self, biosample_guid: SampleGuid, path: &Path) -> Result<DetectedData, AppError> {
        let name = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        let lower = name.to_ascii_lowercase();
        // Binary/structured formats are detected by extension; only text needs a sniff.
        let by_ext = lower.ends_with(".bam")
            || lower.ends_with(".cram")
            || lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
            || [".fasta", ".fa", ".fna", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz"]
                .iter()
                .any(|e| lower.ends_with(e));
        let head = if by_ext { String::new() } else { read_head(path)? };
        let detected = filetype::detect(&name, &head);

        match detected {
            DetectedData::Variants => {
                self.import_variants_from_file(biosample_guid, path, variants::SourceType::Imported)
                    .await?;
            }
            DetectedData::StrProfile => {
                self.import_str_profile_from_csv(biosample_guid, "CUSTOM", None, Some("IMPORTED".into()), path)
                    .await?;
            }
            DetectedData::YSnpPanel => {
                // Build resolved from the subject's alignment, else "hs1" (project default).
                self.import_bisdna_from_file(biosample_guid, path, None).await?;
            }
            DetectedData::ChipData => {
                self.import_chip_profile_from_csv(biosample_guid, None, None, path)
                    .await?;
            }
            DetectedData::MtdnaFasta => {
                self.import_mtdna_from_fasta(biosample_guid, path).await?;
            }
            DetectedData::Alignment => {
                self.import_alignment_file(biosample_guid, path).await?;
            }
            DetectedData::Unknown => {
                return Err(AppError::Import(format!("could not recognize the data in {name}")));
            }
        }
        Ok(detected)
    }

    /// Batch [`add_data`]: expand any directories among `paths` into their recognized data files,
    /// then auto-detect + import each into the subject, collecting a [`BatchImportSummary`]. A
    /// failed/unrecognized file is recorded (not propagated) so one bad file doesn't abort the
    /// batch. `progress(done, total)` ticks per file. The unified multi-file / folder importer
    /// behind the GUI's Add Data button + drag-and-drop. (Distinct from [`import_project_dir`],
    /// which builds a *new* multi-subject project from a NAS layout; this adds to *this* subject.)
    pub async fn add_data_batch(
        &self,
        biosample_guid: SampleGuid,
        paths: Vec<PathBuf>,
        progress: impl Fn(usize, usize),
    ) -> Result<BatchImportSummary, AppError> {
        let mut files = Vec::new();
        for p in &paths {
            collect_data_files(p, &mut files, 0);
        }
        files.dedup();
        let total = files.len();
        let mut summary = BatchImportSummary::default();
        for (i, f) in files.iter().enumerate() {
            progress(i, total);
            let name = f
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_default();
            match self.add_data(biosample_guid, f).await {
                Ok(d) => summary.imported.push((name, d.description().to_string())),
                Err(e) => summary.skipped.push((name, e.to_string())),
            }
        }
        progress(total, total);
        Ok(summary)
    }

    /// Batch-import a NAS project directory: scan `{dir}/{sample}/…` and create the Project
    /// plus its Biosample → SequenceRun → Alignment rows. The reference is resolved per
    /// alignment: pass `Some(fasta)` to use a specific FASTA (validated with its `.fai`) for
    /// every alignment, or `None` to let the gateway resolve each file's inferred build from
    /// the cache. If a needed build isn't cached, returns [`AppError::ReferenceNeeded`]
    /// **before any DB writes** so the UI can prompt + download, then retry. Idempotent: an
    /// existing project (by name), biosample (by donor id), or alignment (by path) is reused.
    /// Coverage is NOT computed here — run it per alignment or via the project report.
    pub async fn import_project_dir(
        &self,
        dir: &Path,
        reference: Option<PathBuf>,
        administrator: String,
        fast_path: bool,
        progress: impl Fn(usize, usize, &str),
    ) -> Result<ProjectImportSummary, AppError> {
        // An explicit FASTA must exist and be indexed; it applies to every alignment.
        if let Some(path) = &reference {
            if !path.exists() {
                return Err(AppError::Import(format!(
                    "reference FASTA not found: {}",
                    path.display()
                )));
            }
            let fai = PathBuf::from(format!("{}.fai", path.display()));
            if !fai.exists() {
                return Err(AppError::Import(format!(
                    "reference FASTA index (.fai) not found: {}",
                    fai.display()
                )));
            }
        }

        let scan_dir = dir.to_path_buf();
        let discovered = tokio::task::spawn_blocking(move || navigator_analysis::scan::scan(&scan_dir)).await??;

        // Resolve each alignment's reference build to a path (explicit FASTA, else the cache).
        // Collect any builds that need downloading and bail before writing anything.
        let explicit = reference.as_ref().map(|p| p.to_string_lossy().into_owned());
        let mut resolved: HashMap<String, String> = HashMap::new();
        let mut needs: Vec<BuildNeed> = Vec::new();
        for sample in &discovered.samples {
            for aln_path in &sample.alignment_files {
                let build = reference_build_for(aln_path);
                if resolved.contains_key(&build) || needs.iter().any(|n| n.build == build) {
                    continue;
                }
                if let Some(ref path) = explicit {
                    resolved.insert(build, path.clone());
                } else if let Some(p) = self.gateway.cached_reference(&build) {
                    resolved.insert(build, p.to_string_lossy().into_owned());
                } else {
                    match self.gateway.reference_status(&build) {
                        RefStatus::NeedsDownload { url, est_bytes } => needs.push(BuildNeed { build, url, est_bytes }),
                        RefStatus::Unknown => {
                            return Err(AppError::Import(format!(
                                "unknown reference build '{build}' — supply a reference FASTA explicitly"
                            )))
                        }
                        RefStatus::Cached(p) | RefStatus::LocalOverride(p) => {
                            resolved.insert(build, p.to_string_lossy().into_owned());
                        }
                    }
                }
            }
        }
        if !needs.is_empty() {
            return Err(AppError::ReferenceNeeded(needs));
        }

        // Project: reuse an existing one with the same name.
        let project = match project::list(self.store.pool())
            .await?
            .into_iter()
            .find(|p| p.name == discovered.project_id)
        {
            Some(p) => p,
            None => {
                self.create_project(NewProject {
                    name: discovered.project_id.clone(),
                    description: None,
                    administrator,
                })
                .await?
            }
        };

        let mut summary = ProjectImportSummary {
            project: project.clone(),
            samples_total: discovered.samples.len(),
            samples_created: 0,
            alignments_created: 0,
            alignments_skipped: 0,
            missing_index: Vec::new(),
            fast_path: FastPathSummary::default(),
        };

        // First tick carries the discovered sample count so the UI meter knows the total up front.
        let total = discovered.samples.len();
        progress(0, total, "");
        for (sample_idx, sample) in discovered.samples.iter().enumerate() {
            // Biosample: reuse by donor identifier within the project.
            let biosample = match biosample::list_for_project(self.store.pool(), project.id)
                .await?
                .into_iter()
                .find(|b| b.donor_identifier == sample.sample_id)
            {
                Some(b) => b,
                None => {
                    summary.samples_created += 1;
                    self.add_biosample(
                        Some(project.id),
                        sample.sample_id.clone(),
                        Some(sample.sample_id.clone()),
                        None,
                    )
                    .await?
                }
            };

            // SequenceRun: reuse the first existing run, else create one (defaults to WGS).
            let run = match sequence_run::list_for_biosample(self.store.pool(), biosample.guid)
                .await?
                .into_iter()
                .next()
            {
                Some(r) => r,
                None => {
                    self.record_sequence_run(NewSequenceRun {
                        biosample_guid: biosample.guid,
                        platform_name: "UNKNOWN".into(),
                        instrument_model: None,
                        test_type: "WGS".into(),
                        library_layout: None,
                        total_reads: None,
                        pf_reads_aligned: None,
                        mean_read_length: None,
                        mean_insert_size: None,
                    })
                    .await?
                }
            };

            let existing = alignment::list_for_run(self.store.pool(), run.id).await?;
            for aln_path in &sample.alignment_files {
                let path_str = aln_path.to_string_lossy().into_owned();
                if existing
                    .iter()
                    .any(|a| a.bam_path.as_deref() == Some(path_str.as_str()))
                {
                    summary.alignments_skipped += 1;
                    continue;
                }
                if !has_sibling_index(aln_path, &sample.index_files) {
                    summary.missing_index.push(sample.sample_id.clone());
                }
                let build = reference_build_for(aln_path);
                let reference_path = resolved.get(&build).cloned();
                self.record_alignment(NewAlignment {
                    sequence_run_id: run.id,
                    reference_build: build,
                    aligner: "unknown".into(),
                    variant_caller: None,
                    bam_path: Some(path_str),
                    reference_path,
                    // Batch import: hash lazily on first analysis (don't stall a bulk NAS import
                    // hashing every multi-GB file up front).
                    content_sha256: None,
                })
                .await?;
                summary.alignments_created += 1;
            }

            // Fast path: ingest the pipeline sidecars onto the build-matching alignment —
            // places Y + mt from the GVCFs and fills sex/metrics/lite-coverage from the text
            // sidecars, no CRAM walk. Best-effort; a failure is tallied and import continues.
            if fast_path && sample.sidecars.has_haplogroup_gvcf() {
                let alns = alignment::list_for_run(self.store.pool(), run.id).await?;
                let chosen = sample
                    .sidecars
                    .build_hint
                    .as_deref()
                    .and_then(|hint| alns.iter().find(|a| build_hint_matches(&a.reference_build, hint)))
                    .or_else(|| alns.iter().find(|a| a.bam_path.is_some()))
                    .or_else(|| alns.first());
                if let Some(a) = chosen {
                    summary.fast_path.samples_with_sidecars += 1;
                    match self.ingest_sidecars(a.id, &sample.sidecars).await {
                        Ok(ing) => {
                            summary.fast_path.y_placed += ing.y_haplogroup.is_some() as usize;
                            summary.fast_path.mt_placed += ing.mt_haplogroup.is_some() as usize;
                            summary.fast_path.sex_filled += ing.sex.is_some() as usize;
                            summary.fast_path.metrics_filled += ing.read_metrics as usize;
                            summary.fast_path.coverage_filled += ing.lite_coverage as usize;
                            for e in ing.errors {
                                summary.fast_path.errors.push(format!("{}: {e}", sample.sample_id));
                            }
                        }
                        Err(e) => summary.fast_path.errors.push(format!("{}: {e}", sample.sample_id)),
                    }
                }
            }
            progress(sample_idx + 1, total, &sample.sample_id);
        }
        Ok(summary)
    }

    /// Cache/override status of a reference build (no network).
    pub fn reference_status(&self, build: &str) -> RefStatus {
        self.gateway.reference_status(build)
    }

    /// Resolve a reference build to a cached, indexed `.fa`, downloading on a miss.
    /// `progress(received, total)` is invoked as bytes arrive.
    pub async fn resolve_reference(
        &self,
        build: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, AppError> {
        Ok(self.gateway.resolve_reference(build, progress).await?)
    }

    /// Resolve (and cache) a liftover chain for a build pair, downloading on a miss. The
    /// cached `.chain` is then available for the haplogroup/liftover path.
    pub async fn resolve_chain(
        &self,
        from: &str,
        to: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, AppError> {
        Ok(self.gateway.resolve_chain(from, to, progress).await?)
    }

    /// Re-hash a cached reference against its integrity sidecar (gap §7) — detects on-disk
    /// corruption of the cached `.fa`. Runs on a blocking thread (re-reads the whole FASTA), so it's
    /// an explicit, user-triggered check (Settings), not the hot path.
    pub async fn verify_reference(&self, build: &str) -> Result<navigator_refgenome::VerifyOutcome, AppError> {
        let gw = self.gateway.clone();
        let build = build.to_string();
        Ok(tokio::task::spawn_blocking(move || gw.verify_reference(&build)).await??)
    }

    /// Lift a whole VCF from `source` build to `target` build (gap §7 — the GATK `LiftoverVcf`
    /// replacement). Ensures the source→target chain and the target reference are resolved
    /// (downloading on a miss, with progress), then runs the line-level lift on a blocking thread.
    /// Returns lift/drop counts.
    pub async fn lift_vcf(
        &self,
        source: &str,
        target: &str,
        in_vcf: PathBuf,
        out_vcf: PathBuf,
        opts: navigator_refgenome::VcfLiftOpts,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<navigator_refgenome::VcfLiftStats, AppError> {
        // Resolve the inputs (chain + target FASTA), downloading on a miss.
        self.gateway.resolve_chain(source, target, progress).await?;
        let target_fa = self.gateway.resolve_reference(target, progress).await?;
        let lo = self.gateway.load_liftover(source, target)?;

        // Target chrY PAR intervals (only needed when filtering them out).
        let target_par: Vec<(i64, i64)> = if opts.filter_par {
            let regions = self.gateway.genome_regions(target, progress).await?;
            regions
                .chromosomes
                .get("chrY")
                .or_else(|| regions.chromosomes.get("Y"))
                .map(|c| c.par.clone())
                .unwrap_or_default()
        } else {
            Vec::new()
        };

        let (src_label, tgt_label) = (source.to_string(), target.to_string());
        let stats = tokio::task::spawn_blocking(move || {
            navigator_refgenome::vcf_lift::lift_vcf(
                &lo,
                &target_fa,
                &target_par,
                &src_label,
                &tgt_label,
                &in_vcf,
                &out_vcf,
                opts,
            )
        })
        .await??;
        Ok(stats)
    }

    pub async fn panel_site_count(&self, panel_id: i64) -> Result<i64, AppError> {
        Ok(panel::site_count(self.store.pool(), panel_id).await?)
    }

    /// Genotype an alignment against a panel at the given ploidy and persist the dosages
    /// (one artifact per alignment+panel+ploidy). Runs the blocking caller off-thread.
    pub async fn genotype_panel(
        &self,
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
    ) -> Result<Vec<SiteGenotype>, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let bam = aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?;
        let sites: Vec<Site> = panel::sites(self.store.pool(), panel_id)
            .await?
            .into_iter()
            .map(|s| Site {
                name: s.name,
                contig: s.chrom,
                position: s.position,
                reference_allele: s.reference_allele,
                alternate_allele: s.alternate_allele,
            })
            .collect();

        let bam_pb = PathBuf::from(bam);
        let reference = aln.reference_path.map(PathBuf::from);
        let params = HaploidCallerParams::default();
        let genotypes = tokio::task::spawn_blocking(move || {
            caller::genotype_sites_all_contigs(&bam_pb, &sites, ploidy, &params, reference.as_deref())
        })
        .await??;

        self.save_analysis(
            alignment_id,
            &panel_kind(panel_id, ploidy),
            caller::GENOTYPE_VERSION,
            &genotypes,
        )
        .await?;
        Ok(genotypes)
    }

    /// Cached panel genotypes for an alignment, if present.
    pub async fn cached_panel_genotypes(
        &self,
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
    ) -> Result<Option<Vec<SiteGenotype>>, AppError> {
        self.load_analysis(alignment_id, &panel_kind(panel_id, ploidy), caller::GENOTYPE_VERSION)
            .await
    }

    /// Resolve an imported chip's genotypes to canonical CHM13 **IBD-panel** dosages — the chip→IBD
    /// path (no alignment, no runtime liftover: the multi-build panel pre-computes coordinates). The
    /// output [`SiteGenotype`]s are over the same CHM13 sites a WGS caller would hit, so a chip and a
    /// WGS sample compare uniformly. Errors if the IBD panel asset isn't built yet.
    pub async fn chip_ibd_dosages(&self, chip_profile_id: i64) -> Result<Vec<SiteGenotype>, AppError> {
        let chip = chip_profile::get(self.store.pool(), chip_profile_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("chip profile {chip_profile_id}"))))?;
        let path = chip.source_path.clone().ok_or_else(|| {
            AppError::Import("this chip has no stored raw-data file — re-import it to enable IBD".into())
        })?;
        let text = std::fs::read_to_string(&path).map_err(|e| AppError::Import(format!("chip file {path}: {e}")))?;
        let from_build = chipprofile::detect_build(&text);
        let calls = chipprofile::autosomal_calls(&text);

        let panel_path = ibd_panel_path(ReferenceBuild::Chm13v2);
        let bytes = read_verified_asset(ReferenceBuild::Chm13v2, &panel_path)?.ok_or_else(|| {
            AppError::Import(format!(
                "IBD panel asset not found at {} — build it with `panelbuild ibd-panel`",
                panel_path.display()
            ))
        })?;
        let panel = navigator_analysis::ibd_panel::IbdPanel::from_bytes(&bytes)?;

        let tuples: Vec<(String, i64, char, char)> =
            calls.into_iter().map(|c| (c.contig, c.position, c.a1, c.a2)).collect();
        Ok(panel.resolve_chip(&from_build, &tuples))
    }

    /// Compare two alignments for IBD, using each one's cached panel genotypes. Both must
    /// have been genotyped against `panel_id` at `ploidy` first.
    pub async fn compare_ibd(
        &self,
        alignment_a: i64,
        alignment_b: i64,
        panel_id: i64,
        ploidy: u8,
        config: IbdDetectorConfig,
    ) -> Result<IbdComparison, AppError> {
        let ga = self
            .cached_panel_genotypes(alignment_a, panel_id, ploidy)
            .await?
            .ok_or_else(|| AppError::NotGenotyped(alignment_a))?;
        let gb = self
            .cached_panel_genotypes(alignment_b, panel_id, ploidy)
            .await?
            .ok_or_else(|| AppError::NotGenotyped(alignment_b))?;

        let build = alignment::get(self.store.pool(), alignment_a)
            .await?
            .and_then(|a| canonical_build(&a.reference_build))
            .unwrap_or(ReferenceBuild::Chm13v2);
        Ok(detect_ibd(&ga, &gb, build, config))
    }

    /// IBD comparison over the **chip-compatible IBD panel** for two samples that may each be a
    /// WGS alignment *or* an imported chip (the volume case). Each source resolves to dosages over
    /// the canonical CHM13 IBD-panel sites ([`Self::ibd_panel_dosages`]); the comparison is then
    /// data-type-agnostic. Requires the IBD panel asset (for the WGS-genotyping / chip-resolve path).
    pub async fn compare_ibd_sources(
        &self,
        a: IbdSource,
        b: IbdSource,
        config: IbdDetectorConfig,
    ) -> Result<IbdComparison, AppError> {
        let ga = self.ibd_panel_dosages(a).await?;
        let gb = self.ibd_panel_dosages(b).await?;
        // The IBD panel is CHM13-coordinate, so the CHM13 genetic map applies to both sources.
        Ok(detect_ibd(&ga, &gb, ReferenceBuild::Chm13v2, config))
    }

    /// IBD comparison between two **subjects** from their autosomal consensuses — each subject's
    /// pooled best genotype per site (across all its WGS + chip sources), no per-source genotyping.
    /// This is the subject-level IBD path (consensus-driven); both subjects must have a built
    /// autosomal consensus. A near-complete genome-wide match is the cross-subject identity (dedup)
    /// signal — read it off the returned [`MatchSummary`]'s relationship estimate.
    pub async fn compare_ibd_consensus(
        &self,
        a: SampleGuid,
        b: SampleGuid,
        config: IbdDetectorConfig,
    ) -> Result<IbdComparison, AppError> {
        let pa = self.cached_autosomal_profile(a).await?.ok_or_else(|| {
            AppError::Import("the first subject has no autosomal consensus yet — build it (Autosomal tab) first".into())
        })?;
        let pb = self.cached_autosomal_profile(b).await?.ok_or_else(|| {
            AppError::Import(
                "the second subject has no autosomal consensus yet — build it (Autosomal tab) first".into(),
            )
        })?;
        let ga = consensus_genotypes(&pa);
        let gb = consensus_genotypes(&pb);
        Ok(detect_ibd(&ga, &gb, ReferenceBuild::Chm13v2, config))
    }

    /// Dosages over the canonical CHM13 IBD-panel sites for a comparison source. A chip resolves
    /// directly ([`Self::chip_ibd_dosages`]); an alignment genotypes the panel's CHM13 sites from
    /// its BAM (cached per alignment, ploidy-2 autosomal).
    pub async fn ibd_panel_dosages(&self, source: IbdSource) -> Result<Vec<SiteGenotype>, AppError> {
        match source {
            IbdSource::Chip(id) => self.chip_ibd_dosages(id).await,
            IbdSource::Alignment(id) => {
                // Salt the cache key with the panel asset's manifest hash, so regenerating the panel
                // (e.g. the probe superset) auto-invalidates stale per-alignment genotypes instead of
                // silently serving genotypes taken over an older site set.
                let kind = ibd_panel_cache_kind();
                if let Some(g) = self.load_analysis(id, &kind, caller::GENOTYPE_VERSION).await? {
                    return Ok(g);
                }
                let aln = self.alignment_or_err(id).await?;
                let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(id))?);
                let reference = aln.reference_path.map(PathBuf::from);
                let panel_path = ibd_panel_path(ReferenceBuild::Chm13v2);
                let bytes = read_verified_asset(ReferenceBuild::Chm13v2, &panel_path)?.ok_or_else(|| {
                    AppError::Import(format!(
                        "IBD panel asset not found at {} — build it with `panelbuild ibd-panel`",
                        panel_path.display()
                    ))
                })?;
                let panel = navigator_analysis::ibd_panel::IbdPanel::from_bytes(&bytes)?;
                let sites: Vec<Site> = panel
                    .sites
                    .iter()
                    .map(|s| Site {
                        name: s.rsid.clone(),
                        contig: s.chm13.contig.clone(),
                        position: s.chm13.position,
                        reference_allele: s.chm13.reference.to_string(),
                        alternate_allele: s.chm13.alternate.to_string(),
                    })
                    .collect();
                let genotypes = tokio::task::spawn_blocking(move || {
                    let params = HaploidCallerParams::default();
                    caller::genotype_sites_all_contigs(&bam, &sites, 2, &params, reference.as_deref())
                })
                .await??;
                self.save_analysis(id, &kind, caller::GENOTYPE_VERSION, &genotypes)
                    .await?;
                Ok(genotypes)
            }
        }
    }

    /// Identity verification — are two alignments the same individual? Autosomal genotype
    /// concordance at the panel sites (primary), corroborated by Y-STR distance when both
    /// subjects have an STR profile. Both alignments must be genotyped against the panel.
    pub async fn verify_identity(
        &self,
        alignment_a: i64,
        alignment_b: i64,
        panel_id: i64,
        ploidy: u8,
    ) -> Result<IdentityVerification, AppError> {
        let ga = self
            .cached_panel_genotypes(alignment_a, panel_id, ploidy)
            .await?
            .ok_or(AppError::NotGenotyped(alignment_a))?;
        let gb = self
            .cached_panel_genotypes(alignment_b, panel_id, ploidy)
            .await?
            .ok_or(AppError::NotGenotyped(alignment_b))?;
        let (matched, sites) = genotype_concordance(&ga, &gb);
        let concordance = (sites > 0).then(|| matched as f64 / sites as f64);

        // Optional Y-STR corroboration from each subject's first STR profile.
        let (mut y_dist, mut y_markers) = (None, 0i64);
        if let (Ok(ba), Ok(bb)) = (
            self.biosample_of_alignment(alignment_a).await,
            self.biosample_of_alignment(alignment_b).await,
        ) {
            let (pa, pb) = (self.list_str_profiles(ba).await?, self.list_str_profiles(bb).await?);
            if let (Some(a), Some(b)) = (pa.first(), pb.first()) {
                let (d, c) = strprofile::str_distance(&a.markers, &b.markers);
                if c > 0 {
                    y_dist = Some(d);
                    y_markers = c;
                }
            }
        }
        Ok(reconciliation::classify_identity(concordance, sites, y_dist, y_markers))
    }

    /// Subject-level identity verification (gap §8) — "are these two subjects the same individual?"
    /// (duplicate detection). The consensus counterpart to [`verify_identity`]: pooled autosomal
    /// consensus genotype concordance (no panel selection), corroborated by Y-STR distance. Both
    /// subjects need a built autosomal consensus.
    pub async fn verify_identity_consensus(
        &self,
        a: SampleGuid,
        b: SampleGuid,
    ) -> Result<IdentityVerification, AppError> {
        let pa = self.cached_autosomal_profile(a).await?.ok_or_else(|| {
            AppError::Import("the first subject has no autosomal consensus yet — build it (Autosomal tab) first".into())
        })?;
        let pb = self.cached_autosomal_profile(b).await?.ok_or_else(|| {
            AppError::Import(
                "the second subject has no autosomal consensus yet — build it (Autosomal tab) first".into(),
            )
        })?;
        let (ga, gb) = (consensus_genotypes(&pa), consensus_genotypes(&pb));
        let (matched, sites) = genotype_concordance(&ga, &gb);
        let concordance = (sites > 0).then(|| matched as f64 / sites as f64);

        // Y-STR corroboration from each subject's first STR profile.
        let (mut y_dist, mut y_markers) = (None, 0i64);
        let (sa, sb) = (self.list_str_profiles(a).await?, self.list_str_profiles(b).await?);
        if let (Some(x), Some(y)) = (sa.first(), sb.first()) {
            let (d, c) = strprofile::str_distance(&x.markers, &y.markers);
            if c > 0 {
                y_dist = Some(d);
                y_markers = c;
            }
        }
        Ok(reconciliation::classify_identity(concordance, sites, y_dist, y_markers))
    }
}
