//! Navigator application/command layer — the single API the UI dispatches to, and the
//! antidote to the `WorkbenchViewModel` god object. Orchestrates `navigator-store` (and
//! later analysis/sync) behind commands and queries; holds policy the old dialogs
//! embedded (identity assignment, existence checks, result (de)serialization). The UI
//! holds only view-state and dispatch — no DB calls or domain decisions in widgets.

use std::collections::HashSet;
use std::path::PathBuf;

use chrono::Utc;
use du_domain::ids::SampleGuid;
use navigator_analysis::caller::{self, HaploidCallerParams, VariantCall};
use navigator_analysis::coverage::{self, CallableLociParams, CoverageResult};

// Re-export the analysis result types the command API returns, so the UI depends only
// on navigator-app (ui -> app), not directly on navigator-analysis.
pub use navigator_analysis::caller::VariantCall as DenovoCall;
pub use navigator_analysis::coverage::CoverageResult as Coverage;
use navigator_domain::workspace::{
    Alignment, AnalysisArtifact, Biosample, NewAlignment, NewProject, NewSequenceRun, Project,
    SequenceRun,
};
use navigator_store::{alignment, artifact, biosample, project, sequence_run, Store, StoreError};
use serde::de::DeserializeOwned;
use serde::Serialize;
use uuid::Uuid;

pub mod error;
pub use error::AppError;

/// Artifact kind for de-novo calls, keyed per contig so different contigs don't
/// overwrite each other in the cache.
fn denovo_kind(contig: &str) -> String {
    format!("denovo_snps:{contig}")
}

/// A project plus a rolled-up count for list/dashboard views.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectOverview {
    pub project: Project,
    pub sample_count: i64,
}

/// The application. Cheap to clone (the store wraps a connection pool).
#[derive(Clone)]
pub struct App {
    store: Store,
}

impl App {
    pub fn new(store: Store) -> Self {
        App { store }
    }

    /// Open/create the workspace database and build the app.
    pub async fn open(path: &std::path::Path) -> Result<Self, AppError> {
        Ok(App::new(Store::open(path).await?))
    }

    // ---- commands ----------------------------------------------------------

    pub async fn create_project(&self, new: NewProject) -> Result<Project, AppError> {
        Ok(project::create(self.store.pool(), &new).await?)
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

    pub async fn record_sequence_run(&self, run: NewSequenceRun) -> Result<SequenceRun, AppError> {
        Ok(sequence_run::create(self.store.pool(), &run).await?)
    }

    pub async fn record_alignment(&self, aln: NewAlignment) -> Result<Alignment, AppError> {
        Ok(alignment::create(self.store.pool(), &aln).await?)
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
        let payload = serde_json::to_string(result)?;
        Ok(artifact::upsert(self.store.pool(), alignment_id, kind, algorithm_version, Utc::now(), &payload).await?)
    }

    /// Load and deserialize a stored analysis result, if present for this version.
    pub async fn load_analysis<T: DeserializeOwned>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
    ) -> Result<Option<T>, AppError> {
        match artifact::get(self.store.pool(), alignment_id, kind, algorithm_version).await? {
            Some(a) => Ok(Some(serde_json::from_str(&a.payload)?)),
            None => Ok(None),
        }
    }

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
            coverage::collect_coverage_callable(&bam, &reference, &params, contig_allowlist.as_ref())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result).await?;
        Ok(result)
    }

    /// Cached `coverage` result for the current algorithm version, if present.
    pub async fn cached_coverage(&self, alignment_id: i64) -> Result<Option<CoverageResult>, AppError> {
        self.load_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION).await
    }

    /// Run coverage using the alignment's own stored BAM/reference paths, then persist.
    /// Errors if the alignment is unknown or has no paths recorded.
    pub async fn run_coverage_for_alignment(&self, alignment_id: i64) -> Result<CoverageResult, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let (Some(bam), Some(reference)) = (aln.bam_path, aln.reference_path) else {
            return Err(AppError::MissingPaths(alignment_id));
        };
        self.run_coverage(
            alignment_id,
            PathBuf::from(bam),
            PathBuf::from(reference),
            None,
            CallableLociParams::default(),
        )
        .await
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
        let kind = denovo_kind(&contig);
        let calls = tokio::task::spawn_blocking(move || {
            caller::call_denovo(&bam, &reference, &contig, &params)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, &kind, caller::DENOVO_VERSION, &calls).await?;
        Ok(calls)
    }

    /// Cached de-novo calls for `contig` at the current caller version, if present.
    pub async fn cached_denovo(&self, alignment_id: i64, contig: &str) -> Result<Option<Vec<VariantCall>>, AppError> {
        self.load_analysis(alignment_id, &denovo_kind(contig), caller::DENOVO_VERSION).await
    }

    /// Run de-novo calling on `contig` using the alignment's own stored paths.
    pub async fn run_denovo_for_alignment(&self, alignment_id: i64, contig: String) -> Result<Vec<VariantCall>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let (Some(bam), Some(reference)) = (aln.bam_path, aln.reference_path) else {
            return Err(AppError::MissingPaths(alignment_id));
        };
        self.run_denovo_caller(
            alignment_id,
            PathBuf::from(bam),
            PathBuf::from(reference),
            contig,
            HaploidCallerParams::default(),
        )
        .await
    }

    // ---- queries -----------------------------------------------------------

    /// Biosamples belonging to a project.
    pub async fn list_biosamples(&self, project_id: i64) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_for_project(self.store.pool(), project_id).await?)
    }

    /// Sequence runs for a biosample.
    pub async fn list_sequence_runs(&self, biosample_guid: SampleGuid) -> Result<Vec<SequenceRun>, AppError> {
        Ok(sequence_run::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    /// Alignments for a sequence run.
    pub async fn list_alignments(&self, sequence_run_id: i64) -> Result<Vec<Alignment>, AppError> {
        Ok(alignment::list_for_run(self.store.pool(), sequence_run_id).await?)
    }

    /// Projects with their sample counts, for a dashboard/list view.
    pub async fn project_overview(&self) -> Result<Vec<ProjectOverview>, AppError> {
        let mut out = Vec::new();
        for project in project::list(self.store.pool()).await? {
            let sample_count = biosample::count_for_project(self.store.pool(), project.id).await?;
            out.push(ProjectOverview { project, sample_count });
        }
        Ok(out)
    }
}
