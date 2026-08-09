//! Error type for the mapping layer (one `thiserror` enum per layer, as elsewhere).

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum AlignError {
    #[error("io error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// The read technology could not be resolved to a mapper preset. Deliberately an error and not
    /// a guess: mapping long reads with a short-read preset (or the reverse) silently produces bad
    /// alignments rather than failing, so an unknown technology has to stop the job and ask.
    #[error("cannot choose a mapper preset for {what} — pass one explicitly")]
    UnknownTechnology { what: String },

    #[error("{0}")]
    Message(String),

    /// The job stopped because cancellation was requested. A distinct variant so callers can tell
    /// a user-requested stop from a failure — same contract as `AnalysisError::Cancelled`.
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
