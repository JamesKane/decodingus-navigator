//! Store-layer error (plan §6: one `thiserror` enum per layer; propagate with `?`).

#[derive(Debug, thiserror::Error)]
pub enum StoreError {
    #[error(transparent)]
    Sqlx(#[from] sqlx::Error),

    #[error("migration failed: {0}")]
    Migrate(#[from] sqlx::migrate::MigrateError),

    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),

    /// A column held data the domain can't decode (e.g. a malformed GUID).
    #[error("decode error: {0}")]
    Decode(String),

    #[error("not found: {0}")]
    NotFound(String),
}

/// Parse a stored GUID string into a [`SampleGuid`], tagging any decode error with `context`
/// (the originating table/column) so a malformed value is traceable.
pub(crate) fn parse_sample_guid(guid: &str, context: &str) -> Result<du_domain::ids::SampleGuid, StoreError> {
    uuid::Uuid::parse_str(guid)
        .map(du_domain::ids::SampleGuid)
        .map_err(|e| StoreError::Decode(format!("{context} guid {guid:?}: {e}")))
}
