//! The store-layer error (plan §6: one `thiserror` enum for each layer, sent up with `?`).

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

/// Parse a stored GUID string into a [`SampleGuid`]. It marks any decode error with `context`,
/// which names the table and column that the value came from. So a reader can trace a bad value.
pub(crate) fn parse_sample_guid(guid: &str, context: &str) -> Result<du_domain::ids::SampleGuid, StoreError> {
    uuid::Uuid::parse_str(guid)
        .map(du_domain::ids::SampleGuid)
        .map_err(|e| StoreError::Decode(format!("{context} guid {guid:?}: {e}")))
}
