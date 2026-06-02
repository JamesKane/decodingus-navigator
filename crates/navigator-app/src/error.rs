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
}
