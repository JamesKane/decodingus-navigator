//! `impl App` methods extracted from `lib.rs` (the `import_unified` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- unified import ----------------------------------------------------

    /// Find the type of a file and send it to the correct importer for a subject. The importers
    /// cover STR data, variants, a chip export, and mtDNA data. The method uses a default value
    /// where it needs one, and it returns the type that it found.
    ///
    /// This method refuses an alignment file. Such a file belongs to a sequence test, and it does
    /// not attach to a subject directly.
    ///
    /// The method also reads the header of a BAM file or a CRAM file. From that header it finds
    /// the build, the aligner, the platform, and the test type. That step is optional.
    pub async fn probe_alignment(&self, path: PathBuf) -> Result<AlignmentProbe, AppError> {
        tokio::task::spawn_blocking(move || navigator_analysis::probe::probe_alignment(&path))
            .await?
            .map_err(AppError::from)
    }

    /// Read a limited count of records from the start of an alignment, and find the identity of the
    /// instrument and the library.
    ///
    /// That identity is the `@RG` tags `SM`, `LB`, and `PU`. It also holds the most frequent
    /// instrument, flowcell, and platform in the read names. The AppView uses those values to find
    /// the laboratory.
    ///
    /// The method runs on another thread, because it blocks on I/O and decodes a CRAM file. A CRAM
    /// file also needs the `reference` value. The step is optional, and each caller continues after
    /// an error.
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

    /// Import an alignment file with no question to the user. The method reads the header of that
    /// file.
    ///
    /// It then makes the sequence run, with the test type, the platform, and the instrument. It also
    /// makes the alignment, with the reference build and the aligner.
    ///
    /// The method does **not** need the reference FASTA file. It finds that file from the build when
    /// a step needs it. When the cache already holds the file, the method stores its path, and each
    /// analysis step then has it at once.
    async fn import_alignment_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        test_type_override: Option<&str>,
    ) -> Result<(), AppError> {
        // A second import is safe for one subject. The code skips the file only when *this*
        // subject already holds the alignment.
        //
        // An earlier version compared across each subject. So the code skipped a file for a new
        // subject when another subject already held it, and it gave no message. The new subject
        // stayed empty, and the app showed an "imported" message that was not true.
        //
        // One case is a second import of a file after the user deleted its earlier subject, when
        // another subject also holds that file.
        let path_str = path.to_string_lossy().into_owned();
        if alignment::list_for_biosample(self.store.pool(), biosample_guid)
            .await?
            .iter()
            .any(|a| a.bam_path.as_deref() == Some(path_str.as_str()))
        {
            return Ok(());
        }
        // The step is optional. After a failed read of the header, the code uses the file name and
        // its default values. It does not stop the import.
        let probe = self.probe_alignment(path.to_path_buf()).await.unwrap_or_default();

        // Find the reference first. The scan of the read names needs it to decode a CRAM file.
        let reference_build = probe
            .reference_build
            .clone()
            .unwrap_or_else(|| reference_build_for(path));
        // Store the path of the reference when the cache holds that file. If not, leave the field
        // empty, and the code finds the file when a step needs it. An import must never wait for a
        // download.
        let reference_path = self
            .gateway
            .cached_reference(&reference_build)
            .map(|p| p.to_string_lossy().into_owned());

        // The scan of the read names gives the identity of the instrument and the library, and the
        // AppView uses those values to find the laboratory.
        //
        // The step is optional. It fills the platform and the model when the `@RG` header holds
        // neither. It also fills the instrument and the flowcell, which no header holds.
        //
        // The code skips this step when it can not read the file. One case is a CRAM file with no
        // reference yet.
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

        // The test type. The code reads the *shape* of the coverage from the BAI index, and that
        // shape corrects the value from the header and the platform.
        //
        // Many reads on chrY, with no read on an autosome, mark a Big Y test, a Y Elite test, or a
        // YSEQ test. Many reads on chrM mark an mtFull test.
        //
        // The step is optional and fast. It costs O(contigs) and reads no record. A CRAM file and a
        // BAM file with no index hold no such profile, and they keep the value from the platform.
        //
        // A value from the caller wins over each value above. One example is a Big_Y-700 directory
        // or a Big_Y-500 directory that the caller recognized. A CRAM file has no `.bai` file, so
        // the detector below can not see the reads on chrY. Without the value from the caller, the
        // code would use the default of the platform, which is WGS.
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
                instrument_model,
                library_layout: stats.as_ref().and_then(|s| s.library_layout.clone()),
                ..NewSequenceRun::new(biosample_guid, platform_name, test_type)
            })
            .await?;

        // Write the identity of the laboratory and the instrument that the code found. The AppView
        // uses those values as its key.
        //
        // The `sequencing_facility` field stays empty. The user sets it, or the AppView lookup
        // gives it from `instrument_id` after that feature ships. See roadmap D8.
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
            // Find the laboratory from the instrument id, through the AppView. The step is
            // optional, and the result goes into the cache.
            //
            // The generation of an FTDNA Big Y test comes from the `@RG LB` label of the header,
            // and the step above already put it in `test_type`. An older header holds no such
            // label. For such a file, the callable area of chrY gives the generation after the
            // analysis, in [`Self::refine_big_y_generation`].
            //
            // This code does not estimate the generation from the laboratory.
            if let Some(inst) = s.instrument_id.as_deref() {
                if let Some(lab) = self.lookup_lab_by_instrument(inst).await {
                    let _ = sequence_run::set_facility(self.store.pool(), run.id, &lab).await;
                }
            }
        }

        // Do not calculate the content hash here. That hash is the identity of the file, and the
        // app uses it to find an old cache entry.
        //
        // A SHA-256 hash of a full alignment of many GB stops this import for some minutes, and the
        // user sees nothing. The batch path also leaves the field `None`. The function
        // `alignment_content_hash` calculates the hash at the first analysis that needs it, and it
        // writes the value to the cache.
        self.record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build,
            aligner: probe.aligner.clone().unwrap_or_else(|| "unknown".into()),
            variant_caller: None,
            bam_path: Some(path.to_string_lossy().into_owned()),
            reference_path,
            content_sha256: None,
            // An imported alignment is an original alignment, and no other row made it. Only a
            // realignment writes these fields, and it adds its own row.
            derived_from_alignment_id: None,
            derivation: None,
        })
        .await?;
        Ok(())
    }

    pub async fn add_data(&self, biosample_guid: SampleGuid, path: &Path) -> Result<DetectedData, AppError> {
        self.add_data_with_test_type(biosample_guid, path, None).await
    }

    /// The same work as [`add_data`], but the caller gives the `test_type` of the sequence run for
    /// an alignment file. The code does not find that value itself. One case is a bulk Big Y import,
    /// where the layout of the directories names the test.
    ///
    /// The method ignores that value for each other kind of file, because the file itself gives the
    /// type.
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
        // The extension of a file gives its type for each binary format and each structured
        // format. Only a text file needs a look at its content.
        //
        // The code now also looks inside a VCF file. A VCF with each site is a 1240K call set. A VCF
        // with the variants only is a plain variant set. See
        // `filetype::looks_like_genotyped_callset_vcf`. So this list holds no `.vcf` pattern.
        let by_ext = lower.ends_with(".bam")
            || lower.ends_with(".cram")
            || lower.ends_with(".geno")
            || lower.ends_with(".snp")
            || lower.ends_with(".ind")
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
            DetectedData::CompleteGenomicsVar => {
                self.import_mastervar_from_file(biosample_guid, path).await?;
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
            DetectedData::EigenstratCallSet => {
                self.import_callset_from_file(biosample_guid, path).await?;
            }
            DetectedData::GvcfCallSet => {
                self.import_gvcf_callset_from_file(biosample_guid, path).await?;
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

    /// The batch form of [`add_data`].
    ///
    /// The method expands each directory in `paths` into the data files that it recognizes. It then
    /// finds the type of each file and imports it into the subject. It collects the results in a
    /// [`BatchImportSummary`] value.
    ///
    /// The summary records a file that failed, and a file that the code does not recognize. The
    /// method does not return an error for such a file, so one bad file does not stop the batch.
    ///
    /// The method calls `progress(done, total)` after each file.
    ///
    /// The Add Data button of the GUI calls this method, and a drag-and-drop action also calls it.
    ///
    /// This method is not [`import_project_dir`]. That method makes a *new* project with many
    /// subjects from a NAS layout. This method adds files to *this* subject.
    pub async fn add_data_batch(
        &self,
        biosample_guid: SampleGuid,
        paths: Vec<PathBuf>,
        progress: impl Fn(usize, usize),
    ) -> Result<BatchImportSummary, AppError> {
        let mut files = Vec::new();
        for p in &paths {
            // Guard against one folder that the user picked and that is the *parent* of the
            // folders of many samples. An FTDNA download root is one example.
            //
            // A read of each folder below it would add the samples of each folder to this one
            // subject, with no message. The method refuses and gives the user a hint. It must not
            // import the wrong data.
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

    /// Import one **sample directory** onto a subject that exists. The CLI `ingest` command uses
    /// this fast path for the D2C bulk load.
    ///
    /// The method reads `dir` as one sample. It records each alignment onto `biosample_guid` from
    /// the **header only**. It decodes no read and scans no library. It then imports each variant
    /// file.
    ///
    /// When the caller sets `fast_path` and the directory holds a haplogroup GVCF file, the method
    /// calls [`Self::ingest_sidecars`]. That method places the Y haplogroup and the mt haplogroup
    /// from the BGZF GVCF files. It also fills the sex, the read metrics, and a small coverage
    /// result from the text sidecar files. It decodes **no CRAM file**.
    ///
    /// A call of [`Self::add_data`] for each file can not do this work. That method can not join a
    /// `*.callable.bed` file, a `coverage.txt` file, or a `stats.txt` file to its alignment. It also
    /// sends a `*.g.vcf.gz` file to the plain-VCF importer, and not to the GVCF fast path.
    ///
    /// A directory with no alignment, no variant file, and no haplogroup GVCF takes another path.
    /// The method then calls [`Self::add_data`] for each file, and each call is optional. So a plain
    /// folder of chip exports, STR exports, and mtDNA exports imports as it did before.
    pub async fn add_sample_dir(
        &self,
        biosample_guid: SampleGuid,
        dir: &Path,
        fast_path: bool,
    ) -> Result<SampleDirSummary, AppError> {
        let scan_dir = dir.to_path_buf();
        let sample = tokio::task::spawn_blocking(move || navigator_analysis::scan::scan_sample(&scan_dir)).await?;
        let mut summary = SampleDirSummary::default();

        // The directory holds no primary sequence data: no alignment, no variant file, and no
        // haplogroup GVCF file. So the code reads it as a set of separate subject files, and it
        // imports each file as add_data does.
        let has_primary = !sample.alignment_files.is_empty()
            || !sample.variant_files.is_empty()
            || sample.sidecars.has_haplogroup_gvcf();
        if !has_primary {
            for f in sample
                .all_files
                .iter()
                .filter(|f| f.kind != navigator_analysis::scan::DiscoveredFileType::Index)
            {
                let name = f
                    .path
                    .file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default();
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
                self.record_sequence_run(NewSequenceRun::new(biosample_guid, "UNKNOWN", "WGS"))
                    .await?
            }
        };

        // Record each alignment from the header only. That read is fast and decodes no record,
        // which is the purpose of the fast path. A second call with the same stored path is safe.
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
            let probe_path = aln_path.clone();
            let (build, _source) = tokio::task::spawn_blocking(move || detect_build_for(&probe_path)).await?;
            let reference_path = self
                .gateway
                .cached_reference(&build)
                .map(|p| p.to_string_lossy().into_owned());
            self.record_alignment(NewAlignment {
                sequence_run_id: run.id,
                reference_build: build,
                aligner: "unknown".into(),
                variant_caller: None,
                bam_path: Some(path_str),
                reference_path,
                content_sha256: None,
                // An imported alignment is an original; see above.
                derived_from_alignment_id: None,
                derivation: None,
            })
            .await?;
            summary.alignments_created += 1;
        }

        // Import the variant files of this directory ONLY when it holds no haplogroup GVCF file.
        //
        // With a GVCF file, the fast path below is the source of the Y value and the mt value. A
        // called `chrY.vcf.gz` file beside it holds the same data, and the GATK layout ships both
        // files.
        //
        // An import of that file starts a second Y placement. A second import of a variant set also
        // adds a second copy, because that import does not compare the content. So a run that
        // continues an earlier run would duplicate the set.
        //
        // A directory with no GVCF file still imports its variant files here. In the b38 aengine
        // layout, for example, the `variants.vcf.gz` file *is* the Y source.
        //
        // A GVCF file has the name `*.g.vcf.gz`, and `scan` also lists it as a variant file. This
        // guard keeps it out of this loop.
        if !sample.sidecars.has_haplogroup_gvcf() {
            for vcf in &sample.variant_files {
                let name = vcf
                    .file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default();
                match self
                    .import_variants_from_file(biosample_guid, vcf, variants::SourceType::Imported)
                    .await
                {
                    Ok(_) => {
                        summary.variants_imported += 1;
                        summary
                            .imported
                            .push((name, DetectedData::Variants.description().to_string()));
                    }
                    Err(e) => summary.skipped.push((name, e.to_string())),
                }
            }
        }

        // The fast path. It places the Y haplogroup and the mt haplogroup from the GVCF files. It
        // fills the sex, the read metrics, and a small coverage result from the text sidecar files,
        // onto the alignment with the same build. It walks no CRAM file. The step is optional, and
        // the chooser in import_project_sample works in the same way.
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

        // The progressive consensus, in docs §7.17. The code adds each autosomal dosage that the
        // store now holds to the consensus of the subject.
        //
        // The step is fast. A chip and a WGS VCF resolve with no decode. A WGS alignment from a
        // recent import has no dosage in the cache. The code skips such an alignment until the
        // batch mode of the panel genotypes it.
        //
        // The step is optional. A fault in the consensus must not fail the import.
        if let Err(e) = self.refresh_autosomal_consensus(biosample_guid).await {
            summary.errors.push(format!("consensus refresh: {e}"));
        }

        Ok(summary)
    }

    /// Import a NAS project directory as a batch. The method reads `{dir}/{sample}/…` and makes
    /// the project with its Biosample, SequenceRun, and Alignment rows.
    ///
    /// The method finds the reference of each alignment. With `Some(fasta)`, it uses that one FASTA
    /// file for each alignment, and it checks the `.fai` file. With `None`, the gateway finds the
    /// build of each file in the cache.
    ///
    /// When the cache holds no file for a build that the import needs, the method returns
    /// [`AppError::ReferenceNeeded`] **before it writes to the database**. The UI can then ask the
    /// user, download the file, and call the method again.
    ///
    /// A second call is safe. The method uses a project with the same name, a biosample with the
    /// same donor id, and an alignment with the same path.
    ///
    /// The method does NOT calculate the coverage. Run that step for one alignment, or from the
    /// project report.
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

    /// Run the sidecar fast path again, for each alignment of a subject whose source directory
    /// still holds the GVCF files of the pipeline.
    ///
    /// The method returns the external Y calls and mt calls from GATK4. An older build ran its
    /// internal walk and replaced those calls, before the app recorded a provenance.
    ///
    /// The method is fast. It reads the small GVCF files and never the CRAM file.
    ///
    /// Each external call goes to its own `:ext` key. So it can replace no other call, and the
    /// "prefer external caller" policy makes it win the consensus.
    ///
    /// The method returns `(y_placed, mt_placed)`. It is the correction for a workspace that a user
    /// imported before the app had external-caller precedence.
    pub async fn reingest_external_for_biosample(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<(usize, usize), AppError> {
        let alns = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let (mut y_placed, mut mt_placed) = (0usize, 0usize);
        for a in &alns {
            let Some(dir) = a.bam_path.as_deref().map(Path::new).and_then(Path::parent) else {
                continue;
            };
            let sample = navigator_analysis::scan::scan_sample(dir);
            if !sample.sidecars.has_haplogroup_gvcf() {
                continue;
            }
            let ingest = self.ingest_sidecars(a.id, &sample.sidecars).await?;
            if ingest.y_haplogroup.is_some() {
                y_placed += 1;
            }
            if ingest.mt_haplogroup.is_some() {
                mt_placed += 1;
            }
        }
        Ok((y_placed, mt_placed))
    }

    /// The work of [`Self::import_project_dir`], with a progress callback for each sample. The
    /// method calls `progress(done, total, sample_id)` before each sample.
    ///
    /// So a large NAS import of some thousands of samples can move a status bar. Without it, the app
    /// looks stopped.
    ///
    /// The `done` value is the 0-based index of the next sample. The first call comes after the
    /// header probe of the preflight, and that step can be slow.
    pub async fn import_project_dir_with_progress(
        &self,
        dir: &Path,
        reference: Option<PathBuf>,
        administrator: String,
        fast_path: bool,
        mut progress: impl FnMut(usize, usize, &str),
    ) -> Result<ProjectImportSummary, AppError> {
        // A FASTA file from the caller must exist and must have an index. It applies to each
        // alignment.
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

        // Find the reference build of each alignment from its **header**. The code reads the
        // header only, so this step is fast and needs no reference FASTA file.
        //
        // The file name is not a reliable source, because most NAS layouts do not put the build in
        // that name. So the code reads the header first and uses the file name second. It also
        // records the source of each build, for the import report.
        let all_paths: Vec<PathBuf> = discovered
            .samples
            .iter()
            .flat_map(|s| s.alignment_files.iter().cloned())
            .collect();
        let detected: HashMap<PathBuf, (String, &'static str)> = tokio::task::spawn_blocking(move || {
            all_paths
                .into_iter()
                .map(|p| {
                    let d = detect_build_for(&p);
                    (p, d)
                })
                .collect()
        })
        .await?;

        // Find a reference path for each *distinct* build that the code detected.
        //
        // A build that the gateway does not recognize takes the CHM13v2.0 default, so the batch
        // continues. A known build with no file in the cache becomes a download that the UI can
        // start.
        //
        // The `effective_of` map takes a build that the code detected and gives the build that the
        // alignment row holds, after each default above.
        let explicit = reference.as_ref().map(|p| p.to_string_lossy().into_owned());
        let mut resolved: HashMap<String, String> = HashMap::new(); // effective build -> FASTA path
        let mut effective_of: HashMap<String, String> = HashMap::new(); // detected build -> effective build
        let mut needs: Vec<BuildNeed> = Vec::new();
        let mut reference_notes: Vec<String> = Vec::new();

        // The count of alignments, and one example of the detection source, for each distinct
        // build.
        let mut per_build: BTreeMap<String, (usize, &'static str)> = BTreeMap::new();
        for (build, source) in detected.values() {
            let e = per_build.entry(build.clone()).or_insert((0, *source));
            e.0 += 1;
        }

        for (detected_build, (count, source)) in &per_build {
            let count = *count;
            // The build that the row holds. Keep the build that the code detected when the gateway
            // recognizes it. A FASTA file from the caller replaces each such value. If neither
            // applies, use the default build. A file with no label then still imports, and the
            // batch continues.
            let (effective, defaulted) =
                if explicit.is_some() || !matches!(self.gateway.reference_status(detected_build), RefStatus::Unknown) {
                    (detected_build.clone(), false)
                } else {
                    (DEFAULT_IMPORT_BUILD.to_string(), true)
                };
            effective_of.insert(detected_build.clone(), effective.clone());

            // Find one FASTA file for that build. The order is the file from the caller, a file
            // that the code already found, the cache, and then the status of the gateway.
            //
            // The method collects each download that the import needs. It records a build with no
            // file and no FASTA path, and the analysis finds that file later. It does not stop the
            // import.
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
                            needs.push(BuildNeed {
                                build: effective.clone(),
                                url,
                                est_bytes,
                            });
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
                (Some(p), false) => format!("{detected_build}: {count} alignment(s) › {p} ({source})"),
                (Some(p), true) => format!(
                    "{detected_build}: {count} alignment(s) › {effective} default › {p} ({source}; build undetectable from header/filename)"
                ),
                (None, _) if needs.iter().any(|n| n.build == effective) => {
                    format!("{detected_build}: {count} alignment(s) › {effective} (needs download)")
                }
                (None, _) => format!(
                    "{detected_build}: {count} alignment(s) › {effective} (no reference available; resolved on demand)"
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

        // Import each sample on its own. A failure in one sample goes into the log and into the
        // `sample_errors` count, and the batch continues with the other samples. The causes are a
        // file that the code can not read, and a fault in the database. One bad sample must not
        // stop the full import.
        let total = discovered.samples.len();
        for (i, sample) in discovered.samples.iter().enumerate() {
            progress(i, total, &sample.sample_id);
            if let Err(e) = self
                .import_project_sample(
                    sample,
                    &project,
                    fast_path,
                    &detected,
                    &effective_of,
                    &resolved,
                    &mut summary,
                )
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

    /// Import the subject, the run, the alignments, and the fast-path sidecar files of one sample.
    ///
    /// This method is separate, so a failure here becomes the error of this sample.
    /// [`Self::import_project_dir`] catches that error, and the batch continues.
    ///
    /// The `detected`, `effective_of`, and `resolved` maps come from the preflight of the caller.
    /// The method writes what this sample added to `summary`.
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
        // The biosample. Use a subject with this donor identifier from **any place in the
        // workspace**. One person is one subject in each project.
        //
        // An earlier version looked in the target project only. So a second import of the same
        // folder, under another project name, made a second subject for each person. A person then
        // had one subject in each project.
        //
        // Make a subject only when the workspace holds none.
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
        // Make sure that the subject is a member of this project. A second call is safe, because
        // the primary key is the pair (guid, project). A subject whose *home* project is another
        // project also joins the list of this project.
        biosample_project::add(
            self.store.pool(),
            biosample.guid,
            project.id,
            None,
            &Utc::now().to_rfc3339(),
        )
        .await?;

        // SequenceRun: reuse the first existing run, else create one (defaults to WGS).
        let run = match sequence_run::list_for_biosample(self.store.pool(), biosample.guid)
            .await?
            .into_iter()
            .next()
        {
            Some(r) => r,
            None => {
                self.record_sequence_run(NewSequenceRun::new(biosample.guid, "UNKNOWN", "WGS"))
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
                // This is a batch import. Calculate the hash at the first analysis. A bulk NAS
                // import must not stop while it hashes each file of many GB.
                content_sha256: None,
                // An imported alignment is an original; see above.
                derived_from_alignment_id: None,
                derivation: None,
            })
            .await?;
            summary.alignments_created += 1;
        }

        // The fast path. It reads the pipeline sidecar files onto the alignment with the same
        // build.
        //
        // It places the Y haplogroup and the mt haplogroup from the GVCF files. It also fills the
        // sex, the metrics, and a small coverage result from the text sidecar files. It walks no
        // CRAM file.
        //
        // The step is optional. A failure goes into a count, and the import continues.
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

    /// Find the indexed `.fa` file of a reference build in the cache. The method downloads that
    /// file when the cache holds none. It calls `progress(received, total)` as each part of the file
    /// arrives.
    pub async fn resolve_reference(
        &self,
        build: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, AppError> {
        Ok(self.gateway.resolve_reference(build, progress).await?)
    }

    /// Find the liftover chain of a build pair, and write it to the cache. The method downloads
    /// that file when the cache holds none. The haplogroup path and the liftover path then read the
    /// `.chain` file from the cache.
    pub async fn resolve_chain(
        &self,
        from: &str,
        to: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, AppError> {
        Ok(self.gateway.resolve_chain(from, to, progress).await?)
    }

    /// Calculate the hash of a reference in the cache again, and compare it with the sidecar file
    /// that holds the correct value. See gap §7. The method finds a `.fa` file that the disk
    /// damaged.
    ///
    /// The method reads the full FASTA file, so it runs on its own thread. The user starts it from
    /// the Settings screen, and no analysis calls it.
    pub async fn verify_reference(&self, build: &str) -> Result<navigator_refgenome::VerifyOutcome, AppError> {
        let gw = self.gateway.clone();
        let build = build.to_string();
        Ok(tokio::task::spawn_blocking(move || gw.verify_reference(&build)).await??)
    }

    /// Move a full VCF file from the `source` build to the `target` build. See gap §7. This method
    /// takes the place of the GATK `LiftoverVcf` tool.
    ///
    /// The method first finds the chain from the source to the target, and the reference of the
    /// target. It downloads each file that the cache does not hold, and it reports the progress.
    ///
    /// It then moves each line on its own thread. It returns the count of the lines that it moved
    /// and the count of the lines that it removed.
    pub async fn lift_vcf(
        &self,
        source: &str,
        target: &str,
        in_vcf: PathBuf,
        out_vcf: PathBuf,
        opts: navigator_refgenome::VcfLiftOpts,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<navigator_refgenome::VcfLiftStats, AppError> {
        // Find the two input files, which are the chain and the FASTA file of the target. Download
        // each file that the cache does not hold.
        self.gateway.resolve_chain(source, target, progress).await?;
        let target_fa = self.gateway.resolve_reference(target, progress).await?;
        let lo = self.gateway.load_liftover(source, target)?;

        // The PAR intervals of chrY on the target build. The code needs them only when it removes
        // those intervals.
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

    /// Make sure that the ancestry asset or IBD asset at `path` is present **and current**.
    /// [`asset_action`] makes the decision that this method acts on.
    ///
    /// The method downloads the asset, and the manifest that it checks the asset against, from the
    /// published GitHub release. A user receives each panel in this way, and no user runs the
    /// offline `panelbuild` tool.
    ///
    /// There are three cases, and the manifest decides each one.
    ///
    /// * The asset is **absent**. Download it when the manifest lists it. An optional asset that
    ///   the team did not publish stays absent, and its feature gives less data.
    /// * The **manifest does not list the asset**. Read the manifest again one time, then test
    ///   again. The code refreshes a cached manifest at no other time. So an installation from
    ///   before the publication of an asset would never learn that the asset exists. That fault
    ///   occurred with `ancestry_haps`.
    /// * The asset is **present with the wrong size**. The team published a new version, so replace
    ///   the file. A test for the file alone can not see a new version. So without this test, an
    ///   installation with the asset keeps the old file for all time.
    ///
    ///   The code moves the old file to another name and does not delete it. It puts that file back
    ///   when the download fails, because an old asset is better than no asset.
    ///
    /// The method never downloads over a path from a `$NAVIGATOR_*` variable, and it never repairs
    /// such a file. That file belongs to the user. Each step is optional, and a network failure
    /// changes nothing on the disk.
    pub(crate) async fn ensure_ancestry_asset(&self, build: ReferenceBuild, path: &Path) -> Result<(), AppError> {
        let Some(name) = path.file_name().and_then(|n| n.to_str()).map(str::to_string) else {
            return Ok(());
        };
        // Download only to the default place in the cache. A path from a variable names a file of
        // the user.
        let default = refgenome_cache::base_dir().join("ancestry").join(&name);
        if path != default {
            return Ok(());
        }
        let base = asset_release_base_url(build);
        let manifest_name = format!("ancestry_manifest_{}.json", build.as_str());
        let manifest_path = default.with_file_name(&manifest_name);

        // (1) The manifest. Download it when the cache holds none. Download it again when it does
        //     not list this asset, because the cache can hold a manifest from before the team
        //     published that asset.
        //
        //     Keep the old copy in memory. A failed download then does not remove the check values
        //     that the app already has.
        let listed = |m: &Option<navigator_analysis::manifest::AssetManifest>| {
            m.as_ref().is_some_and(|m| m.assets.contains_key(&name))
        };
        let mut manifest = load_asset_manifest(build);
        if !listed(&manifest) {
            let previous = std::fs::read(&manifest_path).ok();
            let _ = std::fs::remove_file(&manifest_path); // else the gateway serves the cached copy
            let url = format!("{base}/{manifest_name}");
            match self
                .gateway
                .resolve_ancestry_asset(&manifest_name, &url, &mut |_, _| {})
                .await
            {
                Ok(_) => manifest = load_asset_manifest(build),
                Err(e) => {
                    if let Some(bytes) = previous {
                        let _ = std::fs::write(&manifest_path, bytes);
                    }
                    eprintln!(
                        "ancestry assets: could not fetch {manifest_name} ({e}) — leaving {name} to on-disk state"
                    );
                    return Ok(());
                }
            }
        }
        // Fetch only assets the manifest publishes for this build.
        let Some(entry) = manifest.as_ref().and_then(|m| m.assets.get(&name)).cloned() else {
            return Ok(());
        };

        // (2) The decision about the file on disk. `read_verified_asset` checks the content at
        //     each read. A hash of each asset here costs some seconds at each paint.
        let on_disk = std::fs::metadata(&default).ok().map(|m| m.len());
        let action = asset_action(Some(&entry), on_disk);
        if action == AssetAction::Ready {
            return Ok(());
        }
        // Wrong size → a revision (or a truncated download). Quarantine, then fetch.
        let stale = if action == AssetAction::Replace {
            let aside = default.with_extension("stale");
            match std::fs::rename(&default, &aside) {
                Ok(()) => {
                    eprintln!(
                        "ancestry assets: {name} is {} bytes, the published asset is {} — replacing it",
                        on_disk.unwrap_or(0),
                        entry.bytes
                    );
                    Some(aside)
                }
                Err(e) => {
                    eprintln!("ancestry assets: cannot replace stale {name} ({e}) — keeping the existing file");
                    return Ok(());
                }
            }
        } else {
            eprintln!("ancestry assets: downloading {name} (first use — no build needed) …");
            None
        };

        let mut last = 0u64;
        let fetched = self
            .gateway
            .resolve_ancestry_asset(&name, &format!("{base}/{name}"), &mut |done, total| {
                if done.saturating_sub(last) >= 16 * 1024 * 1024 {
                    last = done;
                    match total {
                        Some(t) if t > 0 => eprintln!("  {name}: {} / {} MB", done / 1_048_576, t / 1_048_576),
                        _ => eprintln!("  {name}: {} MB", done / 1_048_576),
                    }
                }
            })
            .await;
        match fetched {
            Ok(_) => {
                if let Some(aside) = stale {
                    let _ = std::fs::remove_file(aside);
                }
                eprintln!("ancestry assets: {name} ready");
                Ok(())
            }
            Err(e) => {
                // Put the old asset back. An analysis against an old panel is better than no
                // analysis.
                if let Some(aside) = stale {
                    let _ = std::fs::rename(&aside, &default);
                    eprintln!("ancestry assets: {name} download failed ({e}) — kept the existing copy");
                    return Ok(());
                }
                Err(AppError::Import(format!("downloading ancestry asset {name}: {e}")))
            }
        }
    }

    /// Read the CHM13 IBD panel. At the first use, the method downloads the asset from the release,
    /// and no user runs `panelbuild`.
    ///
    /// This method is the one entry point for that panel. It takes the place of five call sites.
    /// Each of those sites gave the error "build it with `panelbuild ibd-panel`" when the asset was
    /// absent.
    pub(crate) async fn load_ibd_panel(&self) -> Result<navigator_analysis::ibd_panel::IbdPanel, AppError> {
        let build = ReferenceBuild::Chm13v2;
        let path = ibd_panel_path(build);
        self.ensure_ancestry_asset(build, &path).await?;
        let bytes = read_verified_asset(build, &path)?.ok_or_else(|| {
            AppError::Import(format!(
                "IBD panel asset not found at {} and could not be downloaded — check your network, or \
                 set NAVIGATOR_IBD_PANEL to a local copy",
                path.display()
            ))
        })?;
        Ok(navigator_analysis::ibd_panel::IbdPanel::from_bytes(&bytes)?)
    }

    /// Change the genotypes of an imported chip into dosages at the canonical CHM13 **IBD panel**
    /// sites. This method is the path from a chip to the IBD data.
    ///
    /// It needs no alignment and does no liftover at run time, because the panel holds the
    /// coordinates of each build.
    ///
    /// The [`SiteGenotype`] values that the method returns cover the same CHM13 sites that a WGS
    /// caller reaches. So a chip sample and a WGS sample compare in the same way.
    ///
    /// The method fails when the IBD panel asset does not exist.
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

        let panel = self.load_ibd_panel().await?;

        let tuples: Vec<(String, i64, char, char)> =
            calls.into_iter().map(|c| (c.contig, c.position, c.a1, c.a2)).collect();
        Ok(panel.resolve_chip(&from_build, &tuples))
    }

    /// Import the autosomal **1240K EIGENSTRAT call set** of a trusted external caller for a
    /// subject. The files are `.geno`, `.snp`, and `.ind`. This method is the autosomal form of the
    /// Y and mt sidecar fast path.
    ///
    /// The method changes the genotypes of the target individual into dosages at the canonical CHM13
    /// panel sites. It decodes no CRAM file, and `resolve_chip` orients each call against the CHM13
    /// alleles itself.
    ///
    /// It stores those dosages as an `external` source, and it builds the autosomal consensus again.
    /// The modern ancestry, the fine ancestry, the deep ancestry, and the IBD steps then read them.
    ///
    /// The `path` value can name any of the three files, because the code finds the other two from
    /// the base name.
    ///
    /// The build of the `.snp` file is GRCh37, which the AADR 1240K set uses. The
    /// `NAVIGATOR_CALLSET_BUILD` variable replaces that value.
    ///
    /// The method returns the count of the panel sites that it resolved.
    pub async fn import_callset_from_file(&self, biosample_guid: SampleGuid, path: &Path) -> Result<usize, AppError> {
        // Resolve the triplet from any member by shared basename.
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| AppError::Import(format!("{}: not an EIGENSTRAT file name", path.display())))?;
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        let geno = dir.join(format!("{stem}.geno"));
        let snp = dir.join(format!("{stem}.snp"));
        let ind = dir.join(format!("{stem}.ind"));
        for (ext, p) in [("geno", &geno), ("snp", &snp), ("ind", &ind)] {
            if !p.is_file() {
                return Err(AppError::Import(format!(
                    "EIGENSTRAT .{ext} not found for {} (expected {})",
                    path.display(),
                    p.display()
                )));
            }
        }

        // The AADR 1240K set uses GRCh37, which is also hg19. The variable lets a user import a
        // call set on GRCh38.
        let build = std::env::var("NAVIGATOR_CALLSET_BUILD").unwrap_or_else(|_| "GRCh37".to_string());
        let (g, s, i, b) = (geno.clone(), snp.clone(), ind.clone(), build.clone());
        let callset =
            tokio::task::spawn_blocking(move || navigator_analysis::callset::read_eigenstrat(&g, &s, &i, None, &b))
                .await??;

        // Resolve to canonical CHM13 panel dosages (same path as a chip; self-orients to CHM13).
        let panel = self.load_ibd_panel().await?;
        let dosages = panel.resolve_chip(&callset.build, &callset.calls);
        let label = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "external call set".into());
        self.store_external_dosages(biosample_guid, &format!("{label} (1240K call set)"), dosages, || {
            format!(
                "the call set resolved to 0 panel sites — check the build (got {} genotypes on {}; \
                 set NAVIGATOR_CALLSET_BUILD if it is not GRCh37)",
                callset.calls.len(),
                callset.build
            )
        })
        .await
    }

    /// Import the autosomal genotypes of a trusted external caller for a subject, from a **diploid
    /// VCF file or gVCF file**. This method is the VCF path of the autosomal fast path, in phases 4
    /// and 5.
    ///
    /// The method reads a GATK4 gVCF file, which holds variant records and hom-ref blocks. It also
    /// reads a VCF file with each site genotyped, such as the output of `bcftools mpileup` and
    /// `bcftools call` across the 1240K sites. In that second file, each site holds a `GT` value.
    ///
    /// The method genotypes the panel loci directly and decodes **no CRAM file**. It then re-keys
    /// each call to the canonical CHM13 sites with `resolve_chip`, stores the dosages as an
    /// `external` source, and builds the autosomal consensus again.
    ///
    /// The code finds the build in the header of the VCF file, and the `NAVIGATOR_CALLSET_BUILD`
    /// variable replaces that value. The method returns the count of the panel sites that it
    /// resolved.
    pub async fn import_gvcf_callset_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
    ) -> Result<usize, AppError> {
        let build = callset_build_for(path);

        let panel = self.load_ibd_panel().await?;

        // The panel loci in the build of the gVCF file. The code groups them by contig and sorts
        // each group. It keeps the reference allele of each site, so a hom-ref block gives the pair
        // (ref, ref).
        let mut targets_by_contig: std::collections::HashMap<String, Vec<i64>> = std::collections::HashMap::new();
        let mut ref_allele: std::collections::HashMap<(String, i64), char> = std::collections::HashMap::new();
        for site in &panel.sites {
            if let Some(l) = site.locus(&build) {
                targets_by_contig.entry(l.contig.clone()).or_default().push(l.position);
                ref_allele.insert((l.contig.clone(), l.position), l.reference);
            }
        }
        for v in targets_by_contig.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        if targets_by_contig.is_empty() {
            return Err(AppError::Import(format!(
                "the IBD panel has no {build} loci — set NAVIGATOR_CALLSET_BUILD to the gVCF's build"
            )));
        }

        let gvcf = path.to_path_buf();
        let params = navigator_analysis::gvcf::GvcfReadParams::default();
        let calls = tokio::task::spawn_blocking(move || {
            navigator_analysis::gvcf::read_diploid_calls(&gvcf, &targets_by_contig, &params)
        })
        .await??;

        // Build reference-forward allele pairs → resolve to CHM13.
        let tuples: Vec<(String, i64, char, char)> = calls
            .into_iter()
            .map(|((contig, pos), call)| {
                let (a1, a2) = match call {
                    navigator_analysis::gvcf::GvcfDiploid::Genotype(a, b) => (a, b),
                    navigator_analysis::gvcf::GvcfDiploid::HomRef => {
                        let r = ref_allele.get(&(contig.clone(), pos)).copied().unwrap_or('N');
                        (r, r)
                    }
                };
                (contig, pos, a1, a2)
            })
            .collect();
        let called = tuples.len();
        let dosages = panel.resolve_chip(&build, &tuples);
        let label = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "external VCF".into());
        self.store_external_dosages(
            biosample_guid,
            &format!("{label} (1240K call set, {build})"),
            dosages,
            || format!("the VCF genotyped 0 panel sites on {build} ({called} calls) — check the build/VCF"),
        )
        .await
    }

    /// Persist a resolved external autosomal call set (CHM13 panel dosages) as an `external` source
    /// and fold it into the autosomal consensus (no decode). Shared by the EIGENSTRAT and gVCF paths.
    /// `on_empty` supplies the error message when nothing resolved.
    async fn store_external_dosages(
        &self,
        biosample_guid: SampleGuid,
        source_label: &str,
        dosages: Vec<navigator_analysis::caller::SiteGenotype>,
        on_empty: impl FnOnce() -> String,
    ) -> Result<usize, AppError> {
        let site_count = dosages.len();
        if site_count == 0 {
            return Err(AppError::Import(on_empty()));
        }
        let json =
            serde_json::to_string(&dosages).map_err(|e| AppError::Import(format!("serializing dosages: {e}")))?;
        let row = navigator_store::external_panel_dosage::StoredPanelDosage {
            biosample_guid: biosample_guid.0.to_string(),
            source_label: source_label.to_string(),
            provenance: navigator_domain::reconciliation::CallProvenance::External
                .as_str()
                .to_string(),
            panel_sig: Some(ibd_panel_cache_kind()),
            site_count: site_count as i64,
            dosages: json,
            created_at: chrono::Utc::now().to_rfc3339(),
        };
        navigator_store::external_panel_dosage::upsert(self.store.pool(), &row).await?;
        // Add the dosages to the autosomal consensus at once. The step is fast, it decodes
        // nothing, and it is optional.
        let _ = self.refresh_autosomal_consensus(biosample_guid).await;
        Ok(site_count)
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

    /// Compare two **subjects** for IBD segments, from their autosomal consensus values.
    ///
    /// The consensus of a subject holds the best genotype at each site, from each of its WGS sources
    /// and chip sources. The method genotypes no source again.
    ///
    /// This path works at the level of a subject, and the consensus drives it. Both subjects must
    /// have an autosomal consensus.
    ///
    /// A match across almost the full genome shows that the two subjects are the same person. Read
    /// that result from the relationship estimate of the [`MatchSummary`] value.
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

    /// The dosages at the canonical CHM13 IBD-panel sites, for one comparison source.
    ///
    /// A chip resolves directly, in [`Self::chip_ibd_dosages`]. An alignment genotypes the CHM13
    /// sites of the panel from its BAM file. The code caches that result for each alignment, and it
    /// uses ploidy 2 on an autosome.
    pub async fn ibd_panel_dosages(&self, source: IbdSource) -> Result<Vec<SiteGenotype>, AppError> {
        match source {
            IbdSource::Chip(id) => self.chip_ibd_dosages(id).await,
            IbdSource::VariantSet(id) => {
                let set = variant_set::get(self.store.pool(), id)
                    .await?
                    .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("variant set {id}"))))?;
                self.variant_set_panel_dosages(&set).await
            }
            IbdSource::Alignment(id) => {
                // Add the manifest hash of the panel asset to the cache key. A new panel, such as
                // one with more probes, then makes each stored genotype of an alignment invalid. So
                // the app does not return a genotype from an older set of sites with no message.
                let kind = ibd_panel_cache_kind();
                if let Some(g) = self.load_analysis(id, &kind, caller::GENOTYPE_VERSION).await? {
                    return Ok(g);
                }
                // Find the reference for the decoder. See alignment_reference_for_decode. A CRAM
                // file needs it, and a BAM file uses None.
                //
                // The panel genotype step counts the reads at each SNP site, and the panel gives the
                // reference allele and the alternate allele. So the code reads no reference base for
                // a BAM file, and it must not start a download for one.
                let build = self.alignment_or_err(id).await?.reference_build;
                let (bam, reference) = self.alignment_reference_for_decode(id).await?;
                let panel = self.load_ibd_panel().await?;
                let is_chm13 = matches!(
                    canonical_build(&build),
                    Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs)
                );

                let genotypes = if is_chm13 {
                    // A native CHM13 alignment. Genotype it directly at the canonical CHM13 loci
                    // of the panel. Each dosage is then already in the CHM13 space, and the code
                    // changes no key.
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
                    tokio::task::spawn_blocking(move || {
                        let params = HaploidCallerParams::default();
                        caller::genotype_sites_all_contigs(
                            &bam,
                            &sites,
                            2,
                            &params,
                            reference.as_deref(),
                            &navigator_analysis::CancelToken::none(),
                        )
                    })
                    .await??
                } else if panel.sites.iter().any(|s| s.locus(&build).is_some()) {
                    // A GRCh37 alignment or a GRCh38 alignment. The panel already holds the
                    // coordinates of that build, from an offline liftover that read the alleles.
                    //
                    // Genotype the sample at the loci of that build. Then re-key each dosage to the
                    // canonical CHM13 sites, in [`IbdPanel::resolve_alignment`]. The code does no
                    // liftover at run time.
                    //
                    // Match the contig names of the panel for that build to the names in the file.
                    // One file uses `chr1`, and another file uses `1`.
                    let (bam_h, ref_h) = (bam.clone(), reference.clone());
                    let file_contigs = tokio::task::spawn_blocking(move || {
                        navigator_analysis::reader::contig_names(&bam_h, ref_h.as_deref())
                    })
                    .await??;
                    let index: HashMap<String, String> = file_contigs
                        .into_iter()
                        .map(|c| (navigator_analysis::contig::bare_upper(&c), c))
                        .collect();
                    let sites: Vec<Site> = panel
                        .sites
                        .iter()
                        .filter_map(|s| {
                            let l = s.locus(&build)?;
                            let key = navigator_analysis::contig::bare_upper(&l.contig);
                            let contig = index.get(&key)?.clone();
                            Some(Site {
                                name: s.rsid.clone(),
                                contig,
                                position: l.position,
                                reference_allele: l.reference.to_string(),
                                alternate_allele: l.alternate.to_string(),
                            })
                        })
                        .collect();
                    let raw = tokio::task::spawn_blocking(move || {
                        let params = HaploidCallerParams::default();
                        caller::genotype_sites_all_contigs(
                            &bam,
                            &sites,
                            2,
                            &params,
                            reference.as_deref(),
                            &navigator_analysis::CancelToken::none(),
                        )
                    })
                    .await??;
                    panel.resolve_alignment(&build, &raw)
                } else {
                    // The panel holds no coordinates for this build, so there is nothing to
                    // genotype. The code gives a smaller result. It must not read the wrong
                    // loci.
                    Vec::new()
                };
                self.save_analysis(id, &kind, caller::GENOTYPE_VERSION, &genotypes)
                    .await?;
                Ok(genotypes)
            }
        }
    }

    /// The IBD-panel dosages of an alignment from the cache. The method genotypes **nothing**.
    ///
    /// It returns `Ok(None)` when no earlier run calculated those dosages. So a caller can work with
    /// the data that the store holds, and it does not start a decode of the full genome.
    ///
    /// [`Self::ibd_panel_dosages`] calculates the dosages and writes them to the cache. This method
    /// only reads them, and the progressive-consensus refresh calls it.
    pub async fn cached_alignment_panel_dosages(
        &self,
        alignment_id: i64,
    ) -> Result<Option<Vec<SiteGenotype>>, AppError> {
        self.load_analysis(alignment_id, &ibd_panel_cache_kind(), caller::GENOTYPE_VERSION)
            .await
    }

    /// A test of identity at the level of a subject, in gap §8. It answers the question "are these
    /// two subjects the same person?", and the app uses it to find a duplicate.
    ///
    /// The method compares the genotypes of the two pooled autosomal consensus values, and it
    /// selects no panel. The distance between the Y-STR values supports that comparison.
    ///
    /// Both subjects need an autosomal consensus.
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

/// The action that [`App::ensure_ancestry_asset`] must take for one asset.
///
/// The decision reads two values. The first is the manifest entry, which is `None` when the manifest
/// does not list the asset. The second is the size of the file on disk, which is `None` when the
/// file is absent.
///
/// The decision uses the size and not a hash. The team makes a new version of an asset when it
/// builds that asset again, and the new file has a different length. The read step already checks
/// the content, and that check has authority.
///
/// A hash of a panel of 133 MB at each call costs some seconds at each paint. It finds only one more
/// case, where two files have the same size and different bytes, and the read step finds that case
/// also.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AssetAction {
    /// The file is present, and its size equals the published size. Use it.
    Ready,
    /// The file is absent. Download it.
    Download,
    /// The file is present, but a newer version exists, or the file is not complete. Move it to
    /// another name and download the new file.
    Replace,
    /// The team published no file of this asset for this build. Leave it absent, and let its
    /// feature give a smaller result.
    Skip,
}

pub(crate) fn asset_action(
    entry: Option<&navigator_analysis::manifest::AssetEntry>,
    on_disk: Option<u64>,
) -> AssetAction {
    match (entry, on_disk) {
        (None, _) => AssetAction::Skip,
        (Some(_), None) => AssetAction::Download,
        (Some(e), Some(len)) if len == e.bytes => AssetAction::Ready,
        (Some(_), Some(_)) => AssetAction::Replace,
    }
}

#[cfg(test)]
mod asset_tests {
    use super::*;
    use navigator_analysis::manifest::AssetEntry;

    fn entry(bytes: u64) -> AssetEntry {
        AssetEntry {
            sha256: "deadbeef".into(),
            bytes,
        }
    }

    #[test]
    fn asset_action_covers_the_present_absent_and_superseded_cases() {
        // Unpublished asset: nothing to fetch whatever is on disk (the feature degrades instead).
        assert_eq!(asset_action(None, None), AssetAction::Skip);
        assert_eq!(asset_action(None, Some(10)), AssetAction::Skip);
        // Published + absent → download; published + right size → use it.
        assert_eq!(asset_action(Some(&entry(100)), None), AssetAction::Download);
        assert_eq!(asset_action(Some(&entry(100)), Some(100)), AssetAction::Ready);
        // The case a plain existence check misses: a locally-present asset the release has revised.
        assert_eq!(
            asset_action(Some(&entry(139_815_581)), Some(13_774_065)),
            AssetAction::Replace
        );
        // …and a truncated download.
        assert_eq!(asset_action(Some(&entry(100)), Some(41)), AssetAction::Replace);
    }
}
