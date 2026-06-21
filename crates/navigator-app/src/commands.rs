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

    /// Delete a project. Refused (with a clear message) while subjects still belong to it, so
    /// the user reassigns them first rather than orphaning the rows.
    pub async fn delete_project(&self, id: i64) -> Result<(), AppError> {
        let members = biosample::count_members_for_project(self.store.pool(), id).await?;
        if members > 0 {
            return Err(AppError::Conflict(format!(
                "cannot delete project: {members} subject(s) still belong to it — reassign them first"
            )));
        }
        if !project::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("project {id}"))));
        }
        Ok(())
    }

    /// Register a biosample, assigning its stable `SampleGuid` here (identity is an
    /// app-layer decision, not the UI's). Verifies the target project exists first so
    /// the caller gets a clear `NotFound` rather than a raw foreign-key error.
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
            guid: SampleGuid(Uuid::new_v4()),
            sample_accession,
            donor_identifier: donor_identifier.into(),
            description: None,
            center_name: None,
            sex,
            project_id,
        };
        biosample::create(self.store.pool(), &b).await?;
        Ok(b)
    }

    /// Update a subject's editable fields (identity, accession, description, center, sex).
    /// Empty strings are normalized to NULL. Returns the updated record.
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

    /// Assign a subject to a project (validating the project exists). `None` clears it.
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

    /// Delete a subject. Refused (with a clear message) when it still has dependent data —
    /// sequencing runs or any imported profile — so the user removes data first rather than
    /// silently orphaning rows.
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
        // The guard above ensures no runs/profiles remain; sweep any derived-only orphans
        // (stale haplogroup/consensus/reconciliation/ancestry/IBD rows from an earlier
        // incomplete delete) so removing the subject can never leave dangling rows.
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

    /// A Y-targeted test (Big Y, Targeted Y, a Y-SNP pack, …) or any Y-STR profile is definitive
    /// evidence of a male subject. Set the biosample's sex to "Male" when such data is present and
    /// it isn't already recorded as male. Best-effort and idempotent — safe to call after any run
    /// or STR-profile import (it re-derives the verdict from the stored data each time).
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

    /// Update a sequence run's descriptive fields (test type required; platform defaults to
    /// "UNKNOWN" when blank; instrument/layout optional). Read metrics are preserved. Returns
    /// the updated record.
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

    /// Update an alignment's descriptive fields (reference build + aligner required; variant
    /// caller optional). File paths are managed by import/probe. Returns the updated record.
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

    /// Fetch an alignment by id, mapping a missing row to a `NotFound` error. The standard way
    /// the analysis/query methods resolve an `alignment_id` before touching its BAM/CRAM.
    pub(crate) async fn alignment_or_err(&self, id: i64) -> Result<Alignment, AppError> {
        alignment::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {id}"))))
    }

    /// Delete a sequence run and everything beneath it (its alignments + cached analysis
    /// artifacts). This is how a mistaken BAM/CRAM import is undone.
    pub async fn delete_sequence_run(&self, id: i64) -> Result<(), AppError> {
        // Capture the run's subject + alignments before the cascade so we can purge any derived
        // haplogroup/consensus data keyed on those alignments (it would otherwise go stale).
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

    /// Merge `secondary` sequence run into `primary` (both must belong to `biosample_guid`):
    /// reparent the secondary run's alignments onto the primary, then delete the now-empty secondary
    /// (its analysis artifacts travel with the alignments — they're alignment-keyed). Destructive +
    /// irreversible. Returns the number of alignments moved.
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
        // The secondary is now empty; delete it (cascade is a no-op for alignments — already moved).
        sequence_run::delete(self.store.pool(), secondary).await?;
        Ok(count)
    }

    /// Delete a single alignment and its cached analysis artifacts (the parent run is kept).
    pub async fn delete_alignment(&self, id: i64) -> Result<(), AppError> {
        // Resolve the subject (via run) before deleting, to purge derived haplogroup/consensus data.
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

    /// Remove derived data keyed on now-deleted alignments: each alignment's Y + mt haplogroup calls
    /// (`aln:<id>` / `aln:<id>:mt`), and the subject's genome-level consensus profiles + painting
    /// (Y/mt/Auto), which were pooled from sources that may no longer exist. The consensus is
    /// recomputable on demand; clearing it makes the displayed haplogroup fall back to reconciling the
    /// remaining cached calls (or nothing), rather than showing a stale placement. A user manual
    /// override is left intact.
    async fn purge_alignment_derived(&self, biosample: SampleGuid, alignment_ids: &[i64]) -> Result<(), AppError> {
        let pool = self.store.pool();
        for &aln in alignment_ids {
            haplogroup_call::delete_one(pool, biosample, DnaType::Y, &format!("aln:{aln}")).await?;
            haplogroup_call::delete_one(pool, biosample, DnaType::Mt, &format!("aln:{aln}:mt")).await?;
            // The per-alignment ancestry estimates die with the alignment.
            ancestry_result::delete_for_alignment(pool, aln).await?;
        }
        for dna in ["Y", "Mt", "Auto"] {
            consensus_profile::delete(pool, biosample, dna).await?;
        }
        consensus_painting::delete(pool, biosample).await?;
        // The audit log describes the consensus we just wiped; clear it so deleting the last run
        // can't leave a stale RUN_RECORDED history pointing at gone alignments. It is re-appended
        // when the consensus is next rebuilt from any remaining calls.
        recon_store::clear_audit(pool, biosample, DnaType::Y).await?;
        recon_store::clear_audit(pool, biosample, DnaType::Mt).await?;
        Ok(())
    }

    /// Reset a subject's analysis: clear **all** sequencing + derived/imported data (runs,
    /// alignments, cached artifacts, Y/mt haplogroups + consensus + reconciliation, ancestry, IBD
    /// results, and chip/STR/variant/mtDNA profiles) while keeping the subject itself — its
    /// identity (name/sex/center), vendor IDs, project memberships, and MDKA genealogy. The
    /// recovery tool for a botched import: clears orphaned/garbage rows so the subject can be
    /// re-imported cleanly. Atomic ([`biosample::clear_data`] runs in one transaction).
    pub async fn clear_biosample_data(&self, guid: SampleGuid) -> Result<(), AppError> {
        biosample::clear_data(self.store.pool(), guid).await?;
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

    /// Like [`save_analysis`] but stamps provenance: `source` (`navigator-walk` |
    /// `pipeline-sidecar`) and `completeness` (`full` | `partial`). The fast-path sidecar
    /// ingest uses this so the manual deep pass can tell a sidecar/partial result apart from a
    /// full walk and upgrade it rather than skip it.
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

    /// The alignment's source-file signature (`mtime:size`) for cache staleness. `None` when the
    /// alignment / its path is gone or unstattable — then the cache is trusted (nothing to
    /// recompute against). Cheap: a metadata stat, no file read (content hashing is the separate,
    /// deferred federation-identity path).
    async fn bam_source_sig(&self, alignment_id: i64) -> Option<String> {
        let aln = alignment::get(self.store.pool(), alignment_id).await.ok().flatten()?;
        file_signature(Path::new(&aln.bam_path?))
    }

    /// `(source, completeness)` of a cached artifact, defaulting `None` columns to
    /// `("navigator-walk", "full")` (pre-provenance rows). `None` when no artifact exists.
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
                // Treat a cached result as a miss when the source file changed since it was computed
                // (BAM-mtime invalidation) — the caller then recomputes + re-stamps it.
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
