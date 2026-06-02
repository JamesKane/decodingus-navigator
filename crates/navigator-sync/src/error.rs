//! Sync-layer error (plan §6: one `thiserror` enum per layer).

#[derive(Debug, thiserror::Error)]
pub enum SyncError {
    #[error(transparent)]
    Atproto(#[from] du_atproto::error::AtprotoError),

    #[error(transparent)]
    Http(#[from] reqwest::Error),

    #[error(transparent)]
    Keyring(#[from] keyring::Error),

    #[error(transparent)]
    Json(#[from] serde_json::Error),

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),

    #[error("oauth error: {0}")]
    Oauth(String),
}
