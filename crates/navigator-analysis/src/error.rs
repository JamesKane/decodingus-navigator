//! Error type for the analysis layer (plan §6: one `thiserror` enum per layer).

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum AnalysisError {
    #[error("io error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("{0}")]
    Message(String),
}

impl AnalysisError {
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        AnalysisError::Io { path: path.into(), source }
    }
}
