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

    /// A local-LLM operation failed (server unreachable, bad response, etc.). The message is
    /// plain-language for the Settings "Test connection" UI.
    #[error("{0}")]
    Llm(String),

    /// An installer-update check failed (GitHub Releases unreachable, bad response, etc.). Surfaced
    /// as a plain-language status line; a failed check is non-fatal (the app runs regardless).
    #[error("update check failed: {0}")]
    Update(String),
}

impl AppError {
    /// Whether this is a user-requested cancellation rather than a genuine failure.
    ///
    /// Cancellation travels as an error so it unwinds the walk from wherever it was, but callers
    /// must be able to tell the two apart: a cancelled run holds a partial result that must not be
    /// persisted, and the UI has to say "cancelled" instead of showing an error. Lives here rather
    /// than in the UI so the layers above never have to reach past `navigator-app` for it.
    pub fn is_cancellation(&self) -> bool {
        matches!(self, AppError::Analysis(navigator_analysis::AnalysisError::Cancelled))
    }
}

impl From<tokio::task::JoinError> for AppError {
    fn from(e: tokio::task::JoinError) -> Self {
        AppError::Join(e.to_string())
    }
}
