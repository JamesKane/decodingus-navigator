//! Store-layer error (plan §6: one `thiserror` enum per layer; propagate with `?`).

#[derive(Debug, thiserror::Error)]
pub enum StoreError {
    #[error(transparent)]
    Sqlx(#[from] sqlx::Error),

    #[error("migration failed: {0}")]
    Migrate(#[from] sqlx::migrate::MigrateError),

    /// A column held data the domain can't decode (e.g. a malformed GUID).
    #[error("decode error: {0}")]
    Decode(String),

    #[error("not found: {0}")]
    NotFound(String),
}
