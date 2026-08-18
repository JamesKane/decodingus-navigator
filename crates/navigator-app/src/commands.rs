//! `impl App` methods extracted from `lib.rs` (the `commands` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- commands ----------------------------------------------------------

    pub async fn create_project(&self, new: NewProject) -> Result<Project, AppError> {
        Ok(project::create(self.store.pool(), &new).await?)
    }

    /// Update a project's editable fields (name required; description optional; administrator
    /// defaults to "unknown" when blank). Returns the updated record.
    pub async fn update_project(
        &self,
        id: i64,
        name: String,
        description: Option<String>,
        administrator: String,
    ) -> Result<Project, AppError> {
        let name = name.trim();
        if name.is_empty() {
            return Err(AppError::Conflict("project name cannot be empty".into()));
        }
        let desc = description.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let admin = administrator.trim();
        let admin = if admin.is_empty() { "unknown" } else { admin };
        let updated = project::update(self.store.pool(), id, name, desc.as_deref(), admin).await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("project {id}"))));
        }
        project::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("project {id}"))))
    }

    /// Delete a project, and remove its members from it first.
    ///
    /// A subject is an independent record, and many projects can hold the same subject. So a delete
    /// of the project keeps each subject. It removes only the membership of that subject in this
    /// project. It also clears the old home column of a subject whose home is this project.
    pub async fn delete_project(&self, id: i64) -> Result<(), AppError> {
        // A project is only a group. A subject is an independent record, and many projects can
        // hold the same subject.
        //
        // A delete does three steps. It removes each membership from the M:N table. It clears the
        // old home column of a subject whose home is this project. It then removes the project.
        // Each subject stays in the workspace.
        //
        // So a user can undo an import that went to the wrong project. The user deletes the project
        // and imports again. Without this behaviour, the message "N subjects still belong to it"
        // stops the user.
        biosample_project::remove_all_for_project(self.store.pool(), id).await?;
        biosample::clear_home_project(self.store.pool(), id).await?;
        if !project::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("project {id}"))));
        }
        Ok(())
    }

    /// Add a biosample and give it a stable `SampleGuid` here. The app layer decides the identity,
    /// and the UI does not. The method checks that the target project exists first. So the caller
    /// receives a clear `NotFound` error and not a raw foreign-key error.
    pub async fn add_biosample(
        &self,
        project_id: Option<i64>,
        donor_identifier: impl Into<String>,
        sample_accession: Option<String>,
        sex: Option<String>,
    ) -> Result<Biosample, AppError> {
        if let Some(pid) = project_id {
            if project::get(self.store.pool(), pid).await?.is_none() {
                return Err(AppError::Store(StoreError::NotFound(format!("project {pid}"))));
            }
        }
        let b = Biosample {
            sample_accession,
            sex,
            project_id,
            ..Biosample::new(SampleGuid(Uuid::new_v4()), donor_identifier)
        };
        biosample::create(self.store.pool(), &b).await?;
        Ok(b)
    }

    /// Change the fields of a subject that the user can edit. They are the identity, the
    /// accession, the description, the center, and the sex. The method changes an empty string to
    /// NULL. It returns the new record.
    pub async fn update_biosample(
        &self,
        guid: SampleGuid,
        donor_identifier: String,
        sample_accession: Option<String>,
        description: Option<String>,
        center_name: Option<String>,
        sex: Option<String>,
    ) -> Result<Biosample, AppError> {
        let donor = donor_identifier.trim();
        if donor.is_empty() {
            return Err(AppError::Conflict("subject identifier cannot be empty".into()));
        }
        let norm = |o: Option<String>| o.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let (acc, desc, center, sex) = (norm(sample_accession), norm(description), norm(center_name), norm(sex));
        let updated = biosample::update(
            self.store.pool(),
            guid,
            donor,
            acc.as_deref(),
            desc.as_deref(),
            center.as_deref(),
            sex.as_deref(),
        )
        .await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))));
        }
        biosample::get(self.store.pool(), guid)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))))
    }

    /// Add a subject to a project. The method checks that the project exists. A value of `None`
    /// removes the subject from its project.
    pub async fn add_biosample_to_project(&self, guid: SampleGuid, project_id: Option<i64>) -> Result<(), AppError> {
        if let Some(pid) = project_id {
            if project::get(self.store.pool(), pid).await?.is_none() {
                return Err(AppError::Store(StoreError::NotFound(format!("project {pid}"))));
            }
        }
        if !biosample::set_project(self.store.pool(), guid, project_id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))));
        }
        Ok(())
    }

    /// Delete a subject.
    ///
    /// The method refuses, and gives a clear message, when the subject still has data. That data is
    /// a sequence run or an imported profile. So the user removes the data first. Without this
    /// guard, the delete leaves rows with no subject and gives no message.
    pub async fn delete_biosample(&self, guid: SampleGuid) -> Result<(), AppError> {
        let runs = self.list_sequence_runs(guid).await?.len();
        let strs = self.list_str_profiles(guid).await?.len();
        let variants = self.list_variant_sets(guid).await?.len();
        let chips = self.list_chip_profiles(guid).await?.len();
        let mt = self.list_mtdna_sequences(guid).await?.len();
        let total = runs + strs + variants + chips + mt;
        if total > 0 {
            return Err(AppError::Conflict(format!(
                "cannot delete subject: it still has {runs} sequencing run(s), {strs} STR, \
                 {variants} variant-set, {chips} chip, {mt} mtDNA record(s) — remove its data first"
            )));
        }
        // The guard above makes sure that no run and no profile stays. This step then removes
        // each derived row with no owner. Such a row is an old haplogroup, consensus,
        // reconciliation, ancestry, or IBD row from a delete that did not complete. So the delete
        // of the subject can never leave a row behind.
        biosample::clear_data(self.store.pool(), guid).await?;
        if !biosample::delete(self.store.pool(), guid).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))));
        }
        Ok(())
    }

    pub async fn record_sequence_run(&self, run: NewSequenceRun) -> Result<SequenceRun, AppError> {
        let guid = run.biosample_guid;
        let created = sequence_run::create(self.store.pool(), &run).await?;
        self.assign_male_for_y_evidence(guid).await?;
        Ok(created)
    }

    /// A Y test or a Y-STR profile is proof that the subject is male. A Y test is a Big Y test, a
    /// Targeted Y test, or a Y-SNP pack.
    ///
    /// The method sets the sex of the biosample to "Male" when such data exists and the record does
    /// not already hold that value.
    ///
    /// The step is optional, and a second call is safe. Call it after any run import or STR-profile
    /// import. It reads the stored data and decides again at each call.
    pub(crate) async fn assign_male_for_y_evidence(&self, guid: SampleGuid) -> Result<(), AppError> {
        use navigator_domain::testtype::{by_code, TargetType};
        let has_y_test = self
            .list_sequence_runs(guid)
            .await?
            .iter()
            .any(|r| by_code(&r.test_type).map(|t| t.target) == Some(TargetType::YChromosome));
        let has_ystr = !self.list_str_profiles(guid).await?.is_empty();
        if !(has_y_test || has_ystr) {
            return Ok(());
        }
        let already_male = biosample::get(self.store.pool(), guid)
            .await?
            .and_then(|b| b.sex)
            .is_some_and(|s| s.trim().eq_ignore_ascii_case("male"));
        if !already_male {
            biosample::set_sex(self.store.pool(), guid, "Male").await?;
        }
        Ok(())
    }

    pub async fn record_alignment(&self, aln: NewAlignment) -> Result<Alignment, AppError> {
        Ok(alignment::create(self.store.pool(), &aln).await?)
    }

    /// Change the descriptive fields of a sequence run. The test type is necessary. A blank
    /// platform becomes "UNKNOWN". The instrument and the layout are optional. The method keeps the
    /// read metrics and returns the new record.
    pub async fn update_sequence_run(
        &self,
        id: i64,
        platform_name: String,
        instrument_model: Option<String>,
        test_type: String,
        library_layout: Option<String>,
        sequencing_facility: Option<String>,
    ) -> Result<SequenceRun, AppError> {
        let test_type = test_type.trim();
        if test_type.is_empty() {
            return Err(AppError::Conflict("test type cannot be empty".into()));
        }
        let norm = |o: Option<String>| o.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let platform = platform_name.trim();
        let platform = if platform.is_empty() { "UNKNOWN" } else { platform };
        let updated = sequence_run::update(
            self.store.pool(),
            id,
            platform,
            norm(instrument_model).as_deref(),
            test_type,
            norm(library_layout).as_deref(),
            norm(sequencing_facility).as_deref(),
        )
        .await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("sequence run {id}"))));
        }
        sequence_run::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("sequence run {id}"))))
    }

    /// Change the descriptive fields of an alignment. The reference build and the aligner are
    /// necessary, and the variant caller is optional. The import step and the probe step control
    /// the file paths. The method returns the new record.
    pub async fn update_alignment(
        &self,
        id: i64,
        reference_build: String,
        aligner: String,
        variant_caller: Option<String>,
    ) -> Result<Alignment, AppError> {
        let build = reference_build.trim();
        let aligner = aligner.trim();
        if build.is_empty() || aligner.is_empty() {
            return Err(AppError::Conflict("reference build and aligner are required".into()));
        }
        let caller = variant_caller.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let updated = alignment::update(self.store.pool(), id, build, aligner, caller.as_deref()).await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("alignment {id}"))));
        }
        self.alignment_or_err(id).await
    }

    /// The alignment with this id, or `None` when the store holds no such row.
    ///
    /// This method is public because a caller now asks about one alignment directly. The provenance
    /// feature caused that change. The UI needs the row to write "realigned to hs1 from alignment
    /// #N". Before, it only listed the files.
    pub async fn alignment(&self, id: i64) -> Result<Option<Alignment>, AppError> {
        Ok(alignment::get(self.store.pool(), id).await?)
    }

    /// Read the alignment with this id, and change an absent row into a `NotFound` error. Each
    /// analysis method and each query method uses this method to resolve an `alignment_id` before
    /// it opens the BAM file or the CRAM file.
    pub(crate) async fn alignment_or_err(&self, id: i64) -> Result<Alignment, AppError> {
        alignment::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {id}"))))
    }

    /// The alignment's BAM/CRAM path, confirmed to still resolve on disk. The standard way a read
    /// path turns an [`Alignment`] into a path to open.
    ///
    /// The two checks belong together, and they belong *early*.
    ///
    /// A recorded path that no longer points to a file is normal in a workspace with a long life. A
    /// user removes old vendor downloads, and a user disconnects a volume.
    ///
    /// No code checked for that state. So the fault appeared as a plain `No such file or directory`
    /// error from inside the reader. By that point the caller had already downloaded a haplotree of
    /// many MB.
    ///
    /// The result was worse than one unclear message. One caller read that io error as an absent
    /// *tree*. It then used the FTDNA tree, and that read failed on the same absent file. One
    /// deleted BAM file gave an incorrect log line, a download with no purpose, and a change of
    /// tree provider with no message.
    ///
    /// Another process can delete the file directly after this check. That result is acceptable.
    /// The check names the common case correctly and at a low cost. It does not make the open
    /// operation safe, and each caller still handles a read error from that operation.
    pub(crate) fn alignment_file(aln: &Alignment) -> Result<PathBuf, AppError> {
        let path = aln.bam_path.clone().ok_or(AppError::MissingPaths(aln.id))?;
        let p = PathBuf::from(&path);
        if !p.exists() {
            return Err(AppError::AlignmentFileMissing { id: aln.id, path });
        }
        Ok(p)
    }

    /// Delete a sequence run and everything beneath it (its alignments + cached analysis
    /// artifacts). This is how a mistaken BAM/CRAM import is undone.
    pub async fn delete_sequence_run(&self, id: i64) -> Result<(), AppError> {
        // Read the subject and the alignments of the run before the cascade. The code then
        // removes each derived haplogroup row and consensus row with a key on those alignments.
        // Without this step, those rows stay and become incorrect.
        let biosample = sequence_run::get(self.store.pool(), id)
            .await?
            .map(|r| r.biosample_guid);
        let alignment_ids: Vec<i64> = alignment::list_for_run(self.store.pool(), id)
            .await?
            .into_iter()
            .map(|a| a.id)
            .collect();
        if !sequence_run::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("sequence run {id}"))));
        }
        if let Some(guid) = biosample {
            self.purge_alignment_derived(guid, &alignment_ids).await?;
        }
        Ok(())
    }

    /// Join the `secondary` sequence run to the `primary` run. Both runs must belong to
    /// `biosample_guid`.
    ///
    /// The method moves each alignment of the secondary run to the primary run. It then deletes the
    /// secondary run, which is now empty. Each analysis artifact moves with its alignment, because
    /// the key of an artifact is the alignment.
    ///
    /// This method destroys data, and the user can not undo it. It returns the count of the
    /// alignments that it moved.
    pub async fn merge_sequence_runs(
        &self,
        biosample_guid: SampleGuid,
        primary: i64,
        secondary: i64,
    ) -> Result<usize, AppError> {
        if primary == secondary {
            return Err(AppError::Import("cannot merge a run into itself".into()));
        }
        // Both runs must exist and belong to this subject (guards a cross-subject merge).
        let runs = sequence_run::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for id in [primary, secondary] {
            if !runs.iter().any(|r| r.id == id) {
                return Err(AppError::Store(StoreError::NotFound(format!(
                    "sequence run {id} for this subject"
                ))));
            }
        }
        let moved = alignment::list_for_run(self.store.pool(), secondary).await?;
        let mut count = 0usize;
        for a in &moved {
            if alignment::set_sequence_run(self.store.pool(), a.id, primary).await? {
                count += 1;
            }
        }
        // The secondary run is now empty, so delete it. The cascade does nothing for the
        // alignments, because the code already moved them.
        sequence_run::delete(self.store.pool(), secondary).await?;
        Ok(count)
    }

    /// Delete one alignment and each analysis artifact in its cache. The method keeps the parent
    /// run.
    pub async fn delete_alignment(&self, id: i64) -> Result<(), AppError> {
        // Find the subject through the run before the delete. The code then removes each derived
        // haplogroup row and consensus row.
        let biosample = match alignment::get(self.store.pool(), id).await? {
            Some(a) => sequence_run::get(self.store.pool(), a.sequence_run_id)
                .await?
                .map(|r| r.biosample_guid),
            None => None,
        };
        if !alignment::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("alignment {id}"))));
        }
        if let Some(guid) = biosample {
            self.purge_alignment_derived(guid, &[id]).await?;
        }
        Ok(())
    }

    /// Remove the derived data whose key is an alignment that the app deleted.
    ///
    /// That data is the Y haplogroup call and the mt haplogroup call of each alignment, which use
    /// the keys `aln:<id>` and `aln:<id>:mt`. It is also the genome-level consensus profiles and
    /// the painting of the subject, for Y, mt, and Auto. The app pooled those results from sources
    /// that can now be absent.
    ///
    /// The app can calculate a consensus again at any time. After this method clears it, the
    /// displayed haplogroup comes from the cached calls that remain, or from nothing. Without this
    /// step, the app shows an old placement.
    ///
    /// The method keeps a value that the user set.
    async fn purge_alignment_derived(&self, biosample: SampleGuid, alignment_ids: &[i64]) -> Result<(), AppError> {
        let pool = self.store.pool();
        for &aln in alignment_ids {
            haplogroup_call::delete_one(pool, biosample, DnaType::Y, &format!("aln:{aln}")).await?;
            haplogroup_call::delete_one(pool, biosample, DnaType::Mt, &format!("aln:{aln}:mt")).await?;
            // The ancestry estimate of each alignment goes with that alignment.
            ancestry_result::delete_for_alignment(pool, aln).await?;
        }
        for dna in ["Y", "Mt", "Auto"] {
            consensus_profile::delete(pool, biosample, dna).await?;
        }
        // Each signature-keyed cache, from the one list. Before this list, the code named three
        // of the four caches by hand. It left the Tier-B archaic segments in the store, with a key
        // on an alignment that the app had deleted.
        for cache in sig_cache::ALL {
            cache.delete(pool, biosample).await?;
        }
        // The audit log describes the consensus that this method removed. Clear the log also.
        // Without that step, a delete of the last run leaves an old RUN_RECORDED entry that names
        // absent alignments. The app writes the log again at the next rebuild of the consensus,
        // from the calls that remain.
        recon_store::clear_audit(pool, biosample, DnaType::Y).await?;
        recon_store::clear_audit(pool, biosample, DnaType::Mt).await?;
        Ok(())
    }

    /// Reset the analysis of a subject. The method clears **all** sequence data and each derived
    /// or imported result.
    ///
    /// It removes the runs, the alignments, and the cached artifacts. It removes the Y and mt
    /// haplogroups with their consensus and reconciliation rows. It also removes the ancestry, the
    /// IBD results, and the chip, STR, variant, and mtDNA profiles.
    ///
    /// It keeps the subject. It also keeps the identity of that subject, which is the name, the
    /// sex, and the center. It keeps the vendor IDs, the project memberships, and the MDKA
    /// genealogy.
    ///
    /// This method is the recovery tool for an import that went wrong. It removes each row with no
    /// owner, so the user can import the subject again. The work is atomic, because
    /// [`biosample::clear_data`] runs in one transaction.
    pub async fn clear_biosample_data(&self, guid: SampleGuid) -> Result<(), AppError> {
        biosample::clear_data(self.store.pool(), guid).await?;
        // The dosages of an imported external autosomal call set are in their own table, outside
        // the cascade of the biosample. Remove them also, so a subject that the user clears holds
        // nothing.
        navigator_store::external_panel_dosage::delete_for_biosample(self.store.pool(), guid).await?;
        Ok(())
    }

    /// Reset only the haplogroup placement of the subject. That placement is the calls, the
    /// consensus, and the override and audit rows, for Y and for mt.
    ///
    /// The method keeps the coverage, the ancestry, and the imported data. It removes an old
    /// lineage, so the next analysis places the subject again. The placement returns at the next
    /// full analysis of a WGS sample, or at the next import of vendor data.
    pub async fn clear_haplogroup_data(&self, guid: SampleGuid) -> Result<(), AppError> {
        biosample::clear_haplogroup_data(self.store.pool(), guid).await?;
        Ok(())
    }

    /// Delete an imported STR profile (and its markers).
    pub async fn delete_str_profile(&self, id: i64) -> Result<(), AppError> {
        if !str_profile::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("STR profile {id}"))));
        }
        Ok(())
    }

    /// Delete an imported variant set (and its calls).
    pub async fn delete_variant_set(&self, id: i64) -> Result<(), AppError> {
        if !variant_set::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("variant set {id}"))));
        }
        Ok(())
    }

    /// Delete an imported chip/array profile.
    pub async fn delete_chip_profile(&self, id: i64) -> Result<(), AppError> {
        if !chip_profile::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("chip profile {id}"))));
        }
        Ok(())
    }

    /// Delete an imported mtDNA sequence.
    pub async fn delete_mtdna_sequence(&self, id: i64) -> Result<(), AppError> {
        if !mtdna_store::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("mtDNA sequence {id}"))));
        }
        Ok(())
    }

    /// Persist a typed analysis result as a versioned artifact (JSON payload). The
    /// `algorithm_version` is part of the cache key, so a newer version supersedes the
    /// old entry. Pair with [`App::load_analysis`].
    pub async fn save_analysis<T: Serialize>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
        result: &T,
    ) -> Result<AnalysisArtifact, AppError> {
        // Default provenance: a full result from a Navigator CRAM walk.
        self.save_analysis_with_provenance(alignment_id, kind, algorithm_version, result, "navigator-walk", "full")
            .await
    }

    /// The same work as [`save_analysis`], but the method also writes the provenance. The
    /// provenance is `source`, which is `navigator-walk` or `pipeline-sidecar`, and `completeness`,
    /// which is `full` or `partial`.
    ///
    /// The fast-path sidecar import uses this method. The manual deep pass can then see the
    /// difference between a partial sidecar result and a full walk. It replaces the partial result
    /// and does not skip it.
    pub async fn save_analysis_with_provenance<T: Serialize>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
        result: &T,
        source: &str,
        completeness: &str,
    ) -> Result<AnalysisArtifact, AppError> {
        let payload = serde_json::to_string(result)?;
        // Stamp the source file's current signature so a later re-align (same path, new content)
        // invalidates this cached result (see `load_analysis`).
        let sig = self.bam_source_sig(alignment_id).await;
        Ok(artifact::upsert(
            self.store.pool(),
            alignment_id,
            kind,
            algorithm_version,
            Utc::now(),
            &payload,
            source,
            completeness,
            sig.as_deref(),
        )
        .await?)
    }

    /// The same work as [`save_analysis_with_provenance`], but the method never replaces a better
    /// artifact with a worse one.
    ///
    /// The store can already hold a result for this `(kind, version)` pair. When the completeness
    /// of that result is the same as the new one, or higher, the method keeps the stored artifact.
    /// One example is a full `navigator-walk` scan against a new `partial` sidecar result.
    ///
    /// The fast-path sidecar import uses this method. So a second import of a project folder can
    /// not replace a real deep scan with the smaller statistics of a sidecar.
    ///
    /// The method returns `true` when it wrote the artifact. It returns `false` when it kept the
    /// stored result, which was the same or better.
    pub async fn save_analysis_no_downgrade<T: Serialize>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
        result: &T,
        source: &str,
        completeness: &str,
    ) -> Result<bool, AppError> {
        if let Some((_src, existing)) = self.analysis_provenance(alignment_id, kind, algorithm_version).await? {
            if completeness_rank(&existing) >= completeness_rank(completeness) {
                return Ok(false);
            }
        }
        self.save_analysis_with_provenance(alignment_id, kind, algorithm_version, result, source, completeness)
            .await?;
        Ok(true)
    }

    /// Write a mark that shows a failed Navigator walk for this alignment. One cause is a CRAM
    /// file that the reader can not decode.
    ///
    /// The store holds the mark as the `error` artifact with the value `"1"`. The project report
    /// then shows a "Failed" cell and not an empty cell. [`clear_analysis_error`] removes the mark
    /// after the next good walk.
    ///
    /// The step is optional. The code hides a failure to write the mark, because the mark is only a
    /// diagnostic.
    pub async fn record_analysis_error(&self, alignment_id: i64, step: &str, message: &str) {
        let mut message = message.to_string();
        message.truncate(500); // keep the payload small; the head carries the cause
        let marker = AnalysisError {
            step: step.to_string(),
            message,
        };
        let _ = self
            .save_analysis(alignment_id, ERROR_KIND, ERROR_VERSION, &marker)
            .await;
    }

    /// Clear any persisted [`record_analysis_error`] marker for this alignment (no-op when absent).
    pub async fn clear_analysis_error(&self, alignment_id: i64) {
        let _ = artifact::delete(self.store.pool(), alignment_id, ERROR_KIND, ERROR_VERSION).await;
    }

    /// The persisted analysis-failure marker for this alignment, if any. Read directly (no
    /// source-mtime freshness check) so the marker shows until an explicit success clears it.
    pub async fn analysis_error(&self, alignment_id: i64) -> Result<Option<AnalysisError>, AppError> {
        match artifact::get(self.store.pool(), alignment_id, ERROR_KIND, ERROR_VERSION).await? {
            Some(a) => Ok(serde_json::from_str(&a.payload).ok()),
            None => Ok(None),
        }
    }

    /// The signature of the source file of the alignment, as `mtime:size`. The code uses it to
    /// find an old cache entry.
    ///
    /// The value is `None` when the alignment is absent, when its path is absent, or when the
    /// operating system can not read the metadata. The code then trusts the cache, because it has
    /// no value to compare.
    ///
    /// The call is fast. It reads the metadata and does not read the file. The content hash is a
    /// separate path for the federation identity, and it runs later.
    async fn bam_source_sig(&self, alignment_id: i64) -> Option<String> {
        let aln = alignment::get(self.store.pool(), alignment_id).await.ok().flatten()?;
        file_signature(Path::new(&aln.bam_path?))
    }

    /// The `(source, completeness)` pair of a cached artifact. A `None` column becomes
    /// `("navigator-walk", "full")`, because a row from before the provenance change holds no
    /// value. The method returns `None` when no artifact exists.
    pub async fn analysis_provenance(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
    ) -> Result<Option<(String, String)>, AppError> {
        Ok(artifact::get(self.store.pool(), alignment_id, kind, algorithm_version)
            .await?
            .map(|a| {
                (
                    a.source.unwrap_or_else(|| "navigator-walk".into()),
                    a.completeness.unwrap_or_else(|| "full".into()),
                )
            }))
    }

    /// Load and deserialize a stored analysis result, if present for this version.
    pub async fn load_analysis<T: DeserializeOwned>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
    ) -> Result<Option<T>, AppError> {
        match artifact::get(self.store.pool(), alignment_id, kind, algorithm_version).await? {
            Some(a) => {
                // Treat a cached result as absent when the source file changed after the
                // calculation. The mtime of the BAM file shows that change. The caller then
                // calculates the result again and writes a new signature.
                let current = self.bam_source_sig(alignment_id).await;
                if !artifact_is_fresh(a.source_sig.as_deref(), current.as_deref()) {
                    return Ok(None);
                }
                Ok(Some(serde_json::from_str(&a.payload)?))
            }
            None => Ok(None),
        }
    }
}

#[cfg(test)]
mod alignment_file_tests {
    use super::*;

    fn aln(id: i64, bam_path: Option<String>) -> Alignment {
        Alignment {
            id,
            sequence_run_id: 1,
            reference_build: "GRCh38".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path,
            reference_path: None,
            content_sha256: None,
            derived_from_alignment_id: None,
            derivation: None,
        }
    }

    #[test]
    fn a_path_that_no_longer_resolves_is_told_apart_from_one_never_recorded() {
        // The distinction the callers act on: `MissingPaths` means the import never recorded a
        // file, which no sweep can do anything about; `AlignmentFileMissing` means it did and the
        // file has since gone, which a sweep skips past.
        assert!(matches!(
            App::alignment_file(&aln(7, None)),
            Err(AppError::MissingPaths(7))
        ));

        let gone = "/Users/nobody/Downloads/FTDNA/23771/2461/2461-gone.bam";
        let err = App::alignment_file(&aln(9, Some(gone.into()))).unwrap_err();
        assert!(err.is_missing_alignment_file(), "got {err:?}");
        // The message names the path, since "which file?" is the user's first question.
        assert!(err.to_string().contains(gone), "{err}");
        assert!(!AppError::MissingPaths(9).is_missing_alignment_file());
    }

    #[test]
    fn a_present_file_resolves() {
        let f = std::env::temp_dir().join("navigator-alignment-file-test.bam");
        std::fs::write(&f, b"x").unwrap();
        let got = App::alignment_file(&aln(3, Some(f.to_string_lossy().into_owned()))).unwrap();
        assert_eq!(got, f);
        let _ = std::fs::remove_file(&f);
    }
}
