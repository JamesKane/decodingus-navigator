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
    async fn import_alignment_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        test_type_override: Option<&str>,
    ) -> Result<(), AppError> {
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
        // An explicit override (e.g. a Big_Y-700/500 directory the caller recognized) wins over
        // inference — CRAMs ship no `.bai`, so the coverage-shape detector below can't see the
        // targeted-Y pile-up and would otherwise fall back to the platform default (WGS).
        let test_type = match test_type_override {
            Some(t) => t.to_string(),
            None => {
                let p = path.to_path_buf();
                let profile = tokio::task::spawn_blocking(move || {
                    navigator_analysis::testtype::coverage_profile_from_bai(&p, None)
                })
                .await
                .ok()
                .flatten();
                navigator_analysis::testtype::infer_test_type(
                    profile.as_ref(),
                    probe.platform.as_deref(),
                    probe.vendor_hint.as_deref(),
                    None,
                    probe.big_y_code.as_deref(),
                )
                .or_else(|| probe.test_type.clone())
                .unwrap_or_else(|| "WGS".into())
            }
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
                s.read_type.as_deref(),
            )
            .await;
            // Resolve the lab from the instrument id via the AppView (best-effort, cached). The
            // FTDNA Big Y generation comes from the header `@RG LB` label (already in `test_type`
            // above) or, on older headers that omit it, from the callable-chrY footprint after
            // analysis ([`Self::refine_big_y_generation`]) — not guessed from the lab here.
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
        self.add_data_with_test_type(biosample_guid, path, None).await
    }

    /// Like [`add_data`], but forces the sequencing-run `test_type` for an alignment file instead
    /// of inferring it (e.g. a bulk Big Y import where the directory layout names the test). The
    /// override is ignored for non-alignment inputs (their type is intrinsic to the file).
    pub async fn add_data_with_test_type(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        test_type: Option<&str>,
    ) -> Result<DetectedData, AppError> {
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
            || lower.ends_with(".vcf.bgz")
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
            DetectedData::FtdnaCsvVariants => {
                self.import_ftdna_csv_variants(biosample_guid, path).await?;
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
                self.import_alignment_file(biosample_guid, path, test_type).await?;
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
            // Guard against a single picked folder that's really a *parent* of several per-sample
            // folders (e.g. an FTDNA download root): recursing it would silently merge sibling
            // samples into this one subject. Refuse with guidance rather than import the wrong data.
            if p.is_dir() {
                let mut these = Vec::new();
                collect_data_files(p, &mut these, 0);
                let subdirs = contributing_subdirs(p, &these);
                if subdirs.len() >= 2 {
                    let sample = subdirs.iter().take(3).cloned().collect::<Vec<_>>().join(", ");
                    return Err(AppError::Import(format!(
                        "{} holds data for {} separate samples ({sample}…) — import one sample's \
                         folder at a time, or use Project Import for a multi-sample directory.",
                        p.display(),
                        subdirs.len(),
                    )));
                }
                files.extend(these);
            } else {
                collect_data_files(p, &mut files, 0);
            }
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

    /// Ingest one staged **sample directory** onto an existing subject (the CLI `ingest` fast path
    /// for the D2C bulk side-load). Scans `dir` into a single sample, records its alignment(s) onto
    /// `biosample_guid` from the **header only** (no read decode / library scan), imports any variant
    /// files, then — when `fast_path` and a haplogroup GVCF is present — runs [`Self::ingest_sidecars`]
    /// to place Y + mt from the BGZF GVCFs and fill sex / read-metrics / lite-coverage from the text
    /// sidecars, still **without decoding the CRAM**. Per-file [`Self::add_data`] can't do this: it
    /// can't group `*.callable.bed` / `coverage.txt` / `stats.txt` to their alignment, and it would
    /// route a `*.g.vcf.gz` through the plain-VCF importer instead of the GVCF haplogroup fast path.
    ///
    /// A directory that holds no alignment, variant, or haplogroup GVCF falls back to a best-effort
    /// per-file [`Self::add_data`] of its contents — so a plain folder of chip/STR/mtDNA exports
    /// still imports as before.
    pub async fn add_sample_dir(
        &self,
        biosample_guid: SampleGuid,
        dir: &Path,
        fast_path: bool,
    ) -> Result<SampleDirSummary, AppError> {
        let scan_dir = dir.to_path_buf();
        let sample = tokio::task::spawn_blocking(move || navigator_analysis::scan::scan_sample(&scan_dir)).await?;
        let mut summary = SampleDirSummary::default();

        // No primary sequencing data (no alignment, no variant, no haplogroup GVCF): treat the
        // directory as a loose bundle of subject files and import each as add_data would.
        let has_primary = !sample.alignment_files.is_empty()
            || !sample.variant_files.is_empty()
            || sample.sidecars.has_haplogroup_gvcf();
        if !has_primary {
            for f in sample
                .all_files
                .iter()
                .filter(|f| f.kind != navigator_analysis::scan::DiscoveredFileType::Index)
            {
                let name = f.path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
                match self.add_data(biosample_guid, &f.path).await {
                    Ok(d) => summary.imported.push((name, d.description().to_string())),
                    Err(e) => summary.skipped.push((name, e.to_string())),
                }
            }
            return Ok(summary);
        }

        // Sequence run: reuse the subject's first run, else create one (WGS default).
        let run = match sequence_run::list_for_biosample(self.store.pool(), biosample_guid)
            .await?
            .into_iter()
            .next()
        {
            Some(r) => r,
            None => {
                self.record_sequence_run(NewSequenceRun {
                    biosample_guid,
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

        // Record each alignment from the header only — cheap, no read decode (the whole point of the
        // fast path). Idempotent on the alignment's stored path.
        let existing = alignment::list_for_run(self.store.pool(), run.id).await?;
        for aln_path in &sample.alignment_files {
            let path_str = aln_path.to_string_lossy().into_owned();
            if existing.iter().any(|a| a.bam_path.as_deref() == Some(path_str.as_str())) {
                summary.alignments_skipped += 1;
                continue;
            }
            let probe_path = aln_path.clone();
            let (build, _source) = tokio::task::spawn_blocking(move || detect_build_for(&probe_path)).await?;
            let reference_path = self.gateway.cached_reference(&build).map(|p| p.to_string_lossy().into_owned());
            self.record_alignment(NewAlignment {
                sequence_run_id: run.id,
                reference_build: build,
                aligner: "unknown".into(),
                variant_caller: None,
                bam_path: Some(path_str),
                reference_path,
                content_sha256: None,
            })
            .await?;
            summary.alignments_created += 1;
        }

        // Import bundled variant files ONLY when there is no haplogroup GVCF. When a GVCF is present
        // the fast path below is the authoritative Y/mt source, so a called `chrY.vcf.gz` sitting
        // beside it (the GATK repo layout ships both) is redundant — importing it would fire a second
        // Y placement and, because variant-set import isn't content-idempotent, would duplicate the
        // set on a resumable re-run. Non-GVCF tiers (e.g. the b38 aengine `variants.vcf.gz`) still
        // import here: there the VCF *is* the Y source. GVCFs themselves are `.g.vcf.gz`, which `scan`
        // also lists as variant files — the guard keeps them out of this loop too.
        if !sample.sidecars.has_haplogroup_gvcf() {
            for vcf in &sample.variant_files {
                let name = vcf.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
                match self
                    .import_variants_from_file(biosample_guid, vcf, variants::SourceType::Imported)
                    .await
                {
                    Ok(_) => {
                        summary.variants_imported += 1;
                        summary.imported.push((name, DetectedData::Variants.description().to_string()));
                    }
                    Err(e) => summary.skipped.push((name, e.to_string())),
                }
            }
        }

        // Fast path: place Y + mt from the GVCFs and fill sex / read-metrics / lite-coverage from the
        // text sidecars onto the build-matching alignment — no CRAM walk. Best-effort (mirrors the
        // project-import chooser at import_project_sample).
        if fast_path && sample.sidecars.has_haplogroup_gvcf() {
            let alns = alignment::list_for_run(self.store.pool(), run.id).await?;
            let chosen = sample
                .sidecars
                .build_hint
                .as_deref()
                .and_then(|hint| alns.iter().find(|a| build_hint_matches(&a.reference_build, hint)))
                .or_else(|| alns.iter().find(|a| a.bam_path.is_some()))
                .or_else(|| alns.first());
            match chosen {
                Some(a) => match self.ingest_sidecars(a.id, &sample.sidecars).await {
                    Ok(ing) => {
                        summary.sidecars_ingested = true;
                        summary.y_haplogroup = ing.y_haplogroup;
                        summary.mt_haplogroup = ing.mt_haplogroup;
                        summary.sex = ing.sex;
                        summary.read_metrics = ing.read_metrics;
                        summary.lite_coverage = ing.lite_coverage;
                        summary.errors.extend(ing.errors);
                    }
                    Err(e) => summary.errors.push(format!("sidecar ingest: {e}")),
                },
                None => summary
                    .errors
                    .push("haplogroup GVCF present but no alignment to attach it to".into()),
            }
        }

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
    ) -> Result<ProjectImportSummary, AppError> {
        self.import_project_dir_with_progress(dir, reference, administrator, fast_path, |_, _, _| {})
            .await
    }

    /// [`Self::import_project_dir`] with a per-sample progress callback `progress(done, total,
    /// sample_id)`, invoked before each sample so a large NAS import (thousands of samples) can
    /// stream a status bar instead of appearing frozen. `done` is the 0-based index about to
    /// process; the first call fires only after the (potentially slow) header-probe pre-flight.
    pub async fn import_project_dir_with_progress(
        &self,
        dir: &Path,
        reference: Option<PathBuf>,
        administrator: String,
        fast_path: bool,
        mut progress: impl FnMut(usize, usize, &str),
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

        // Detect each alignment's reference build from its **header** (only the header, so it's
        // cheap and needs no reference FASTA). The filename is an unreliable signal — most NAS
        // project layouts don't put the build in the name — so probe first, fall back to the
        // filename, and record how each build was decided for the import diagnostics.
        let all_paths: Vec<PathBuf> = discovered
            .samples
            .iter()
            .flat_map(|s| s.alignment_files.iter().cloned())
            .collect();
        let detected: HashMap<PathBuf, (String, &'static str)> = tokio::task::spawn_blocking(move || {
            all_paths.into_iter().map(|p| { let d = detect_build_for(&p); (p, d) }).collect()
        })
        .await?;

        // Resolve each *distinct* detected build to a reference path. A build the gateway can't
        // canonicalize falls back to the CHM13v2.0 default rather than aborting the whole batch;
        // a known build that isn't cached is surfaced as a recoverable download need. `effective_of`
        // maps a detected build to the one actually stored on the alignment (after any fallback).
        let explicit = reference.as_ref().map(|p| p.to_string_lossy().into_owned());
        let mut resolved: HashMap<String, String> = HashMap::new(); // effective build -> FASTA path
        let mut effective_of: HashMap<String, String> = HashMap::new(); // detected build -> effective build
        let mut needs: Vec<BuildNeed> = Vec::new();
        let mut reference_notes: Vec<String> = Vec::new();

        // Alignment count + a representative detection source, per distinct detected build.
        let mut per_build: BTreeMap<String, (usize, &'static str)> = BTreeMap::new();
        for (build, source) in detected.values() {
            let e = per_build.entry(build.clone()).or_insert((0, *source));
            e.0 += 1;
        }

        for (detected_build, (count, source)) in &per_build {
            let count = *count;
            // Effective build: keep the detected one when the gateway recognizes it (or an explicit
            // FASTA overrides everything); otherwise fall back to the default so unlabeled files
            // still import instead of killing the batch.
            let (effective, defaulted) = if explicit.is_some()
                || !matches!(self.gateway.reference_status(detected_build), RefStatus::Unknown)
            {
                (detected_build.clone(), false)
            } else {
                (DEFAULT_IMPORT_BUILD.to_string(), true)
            };
            effective_of.insert(detected_build.clone(), effective.clone());

            // Resolve the effective build to a FASTA once (explicit > already-resolved > cache >
            // gateway status). A download need is collected; an unresolvable build is recorded
            // without a FASTA (resolved on demand at analysis time) rather than aborting.
            let path: Option<String> = if let Some(ref p) = explicit {
                Some(p.clone())
            } else if let Some(p) = resolved.get(&effective) {
                Some(p.clone())
            } else if let Some(p) = self.gateway.cached_reference(&effective) {
                Some(p.to_string_lossy().into_owned())
            } else {
                match self.gateway.reference_status(&effective) {
                    RefStatus::Cached(p) | RefStatus::LocalOverride(p) => Some(p.to_string_lossy().into_owned()),
                    RefStatus::NeedsDownload { url, est_bytes } => {
                        if !needs.iter().any(|n| n.build == effective) {
                            needs.push(BuildNeed { build: effective.clone(), url, est_bytes });
                        }
                        None
                    }
                    RefStatus::Unknown => None,
                }
            };
            if let Some(ref p) = path {
                resolved.entry(effective.clone()).or_insert_with(|| p.clone());
            }

            let note = match (&path, defaulted) {
                (Some(p), false) => format!("{detected_build}: {count} alignment(s) → {p} ({source})"),
                (Some(p), true) => format!(
                    "{detected_build}: {count} alignment(s) → {effective} default → {p} ({source}; build undetectable from header/filename)"
                ),
                (None, _) if needs.iter().any(|n| n.build == effective) => {
                    format!("{detected_build}: {count} alignment(s) → {effective} (needs download)")
                }
                (None, _) => format!(
                    "{detected_build}: {count} alignment(s) → {effective} (no reference available; resolved on demand)"
                ),
            };
            eprintln!("project import: {note}");
            reference_notes.push(note);
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
            sample_errors: Vec::new(),
            reference_notes,
            fast_path: FastPathSummary::default(),
        };

        // Import each sample independently: a single sample's failure (unreadable file, DB hiccup)
        // is logged + tallied into `sample_errors` and the batch continues with the rest, rather
        // than one bad sample aborting the whole import.
        let total = discovered.samples.len();
        for (i, sample) in discovered.samples.iter().enumerate() {
            progress(i, total, &sample.sample_id);
            if let Err(e) = self
                .import_project_sample(sample, &project, fast_path, &detected, &effective_of, &resolved, &mut summary)
                .await
            {
                eprintln!(
                    "project import: sample {} failed ({e}); skipping and continuing with the rest",
                    sample.sample_id
                );
                summary.sample_errors.push(format!("{}: {e}", sample.sample_id));
            }
        }
        Ok(summary)
    }

    /// Import one sample's subject, run, alignments, and fast-path sidecars. Extracted so a failure
    /// here bubbles up as this sample's error (caught by [`Self::import_project_dir`]) instead of
    /// aborting the whole batch. `detected`/`effective_of`/`resolved` are the pre-flight reference
    /// maps from the caller; `summary` is updated in place with what this sample contributed.
    #[allow(clippy::too_many_arguments)]
    async fn import_project_sample(
        &self,
        sample: &navigator_analysis::scan::DiscoveredSample,
        project: &Project,
        fast_path: bool,
        detected: &HashMap<PathBuf, (String, &'static str)>,
        effective_of: &HashMap<String, String>,
        resolved: &HashMap<String, String>,
        summary: &mut ProjectImportSummary,
    ) -> Result<(), AppError> {
        // Biosample: reuse an existing subject with this donor identifier **anywhere in the
        // workspace** — a person is one subject across projects. Scoping the lookup to the target
        // project duplicated everyone when the same folder was re-imported under a different
        // project name (a person then existed once per project). Create only when truly new.
        let biosample = match biosample::find_by_donor(self.store.pool(), &sample.sample_id).await? {
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
        // Ensure the subject is a member of this project (idempotent on the (guid, project) PK).
        // A reused subject whose *home* project is another one still joins this project's roster.
        biosample_project::add(self.store.pool(), biosample.guid, project.id, None, &Utc::now().to_rfc3339())
            .await?;

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
            // Store the *effective* build (the detected one, or the CHM13v2.0 fallback) so every
            // downstream analysis step reads the same reference the pre-flight resolved.
            let detected_build = detected
                .get(aln_path)
                .map(|(b, _)| b.clone())
                .unwrap_or_else(|| reference_build_for(aln_path));
            let build = effective_of.get(&detected_build).cloned().unwrap_or(detected_build);
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
        Ok(())
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
