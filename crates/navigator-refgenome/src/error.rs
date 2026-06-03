//! Error type for the reference-genome layer.

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum RefgenomeError {
    #[error("http error fetching {url}: {source}")]
    Http {
        url: String,
        #[source]
        source: reqwest::Error,
    },

    #[error("io error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("unknown reference build: {0}")]
    UnknownBuild(String),

    #[error("no liftover chain registered for {from} -> {to}")]
    NoChain { from: String, to: String },

    #[error("{0}")]
    Message(String),
}

impl RefgenomeError {
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        RefgenomeError::Io { path: path.into(), source }
    }
}
