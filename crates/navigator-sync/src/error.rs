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

    /// A device-key crypto/encoding fault (bad seed length, corrupt base64).
    #[error("device key error: {0}")]
    Crypto(String),

    /// The access token was rejected (HTTP 401) — refresh and retry.
    #[error("unauthorized (token expired or revoked)")]
    Unauthorized,

    /// A 5xx from the PDS/auth server — transient, worth retrying.
    #[error("server error {0}: {1}")]
    Server(u16, String),
}

impl SyncError {
    /// Whether retrying the same request later might succeed: transport failures
    /// (offline, timeout) and 5xx server errors. 4xx/validation/auth errors are not.
    pub fn is_transient(&self) -> bool {
        match self {
            SyncError::Http(e) => e.is_connect() || e.is_timeout() || e.is_request(),
            SyncError::Server(code, _) => *code >= 500,
            _ => false,
        }
    }
}
