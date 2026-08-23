//! Error type for the mapping layer (one `thiserror` enum for each layer, as elsewhere).

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum AlignError {
    #[error("io error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// The read technology does not resolve to a mapper preset. This is an error and not a guess.
    /// A short-read preset that maps long reads makes bad alignments, and so does the reverse. The
    /// mapper gives no warning when it does this. An unknown technology must stop the job and ask.
    #[error("cannot choose a mapper preset for {what} — pass one explicitly")]
    UnknownTechnology { what: String },

    #[error("{0}")]
    Message(String),

    /// The user cancelled the job. This is a different variant, so that a caller can tell a stop
    /// that the user asked for from a failure. The contract is the same as `AnalysisError::Cancelled`.
    #[error("cancelled")]
    Cancelled,
}

impl AlignError {
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        AlignError::Io {
            path: path.into(),
            source,
        }
    }
}
