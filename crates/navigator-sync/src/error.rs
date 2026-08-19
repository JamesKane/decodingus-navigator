//! Sync-layer error (plan §6: one `thiserror` enum for each layer).

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

    /// The server refused the access token (HTTP 401). Refresh it and try again.
    #[error("unauthorized (token expired or revoked)")]
    Unauthorized,

    /// A 5xx from the PDS or the auth server. It is transient, so a second try can work.
    #[error("server error {0}: {1}")]
    Server(u16, String),
}

impl SyncError {
    /// Whether a second try of the same request, later, can succeed. A transport failure, such as
    /// offline or a timeout, and a 5xx server error, can. A 4xx, a validation error, and an auth
    /// error can not.
    pub fn is_transient(&self) -> bool {
        match self {
            SyncError::Http(e) => e.is_connect() || e.is_timeout() || e.is_request(),
            SyncError::Server(code, _) => *code >= 500,
            _ => false,
        }
    }
}
