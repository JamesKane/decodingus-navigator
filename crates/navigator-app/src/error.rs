//! Application-layer error: store failures plus artifact (de)serialization.

#[derive(Debug, thiserror::Error)]
pub enum AppError {
    #[error(transparent)]
    Store(#[from] navigator_store::StoreError),

    #[error(transparent)]
    Analysis(#[from] navigator_analysis::AnalysisError),

    #[error("serialization error: {0}")]
    Serde(#[from] serde_json::Error),

    /// A blocking analysis task failed to join (panicked or was cancelled).
    #[error("analysis task failed: {0}")]
    Join(String),

    #[error("alignment {0} has no BAM/reference path recorded")]
    MissingPaths(i64),

    #[error("alignment {0} has not been genotyped against this panel")]
    NotGenotyped(i64),

    /// The ancestry reference panel file is missing — build it with `navigator-panelbuild`
    /// and install it (or set `$NAVIGATOR_ANCESTRY_PANEL`).
    #[error("ancestry panel not found at {0} — build it with navigator-panelbuild")]
    AncestryPanelMissing(std::path::PathBuf),

    /// The bundled panel is for a different reference build than the alignment.
    #[error("ancestry panel build {panel} does not match alignment build {alignment}")]
    AncestryPanelBuildMismatch { panel: String, alignment: String },

    /// Too few sites genotyped for a reliable ancestry estimate.
    #[error("insufficient data for ancestry: {genotyped} SNPs genotyped, {required} required")]
    InsufficientAncestryData { genotyped: usize, required: usize },

    #[error("not signed in — log in to a PDS account first")]
    NotAuthenticated,

    /// An AppView API call failed (e.g. federated IBD). 403 → the device key isn't
    /// registered/verified yet; 422 → clock skew; otherwise the server's reason.
    #[error("appview error: {0}")]
    AppView(String),

    #[error("could not read file: {0}")]
    Io(#[from] std::io::Error),

    #[error("import error: {0}")]
    Import(String),

    #[error(transparent)]
    Sync(#[from] navigator_sync::SyncError),

    #[error(transparent)]
    Refgenome(#[from] navigator_refgenome::RefgenomeError),

    /// Import needs reference build(s) that aren't cached — the UI prompts, downloads via
    /// the gateway, then retries. No DB writes happened.
    #[error("reference download required: {0:?}")]
    ReferenceNeeded(Vec<crate::BuildNeed>),

    /// A requested mutation is refused because of current state (e.g. deleting a subject that
    /// still has sequencing data or profiles).
    #[error("{0}")]
    Conflict(String),
}
