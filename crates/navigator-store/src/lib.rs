//! Navigator local persistence — SQLite via `sqlx`, replacing the H2/Slick layer.
//! Query-module-per-aggregate over a `SqlitePool`; complex children modelled as proper
//! rows, with versioned, reversible migrations. Persisted state is authoritative.

use std::path::Path;
use std::str::FromStr;

use sqlx::sqlite::{SqliteConnectOptions, SqlitePoolOptions};
use sqlx::SqlitePool;

pub mod alignment;
pub mod artifact;
pub mod biosample;
pub mod chip_profile;
pub mod mtdna;
pub mod error;
pub mod panel;
pub mod project;
pub mod sequence_run;
pub mod str_profile;
pub mod variant_set;

pub use error::StoreError;

static MIGRATOR: sqlx::migrate::Migrator = sqlx::migrate!("./migrations");

/// Handle to the workspace database.
#[derive(Clone, Debug)]
pub struct Store {
    pool: SqlitePool,
}

impl Store {
    /// Open (creating if absent) the SQLite database at `path` and run migrations.
    pub async fn open(path: &Path) -> Result<Self, StoreError> {
        let opts = SqliteConnectOptions::new()
            .filename(path)
            .create_if_missing(true)
            .foreign_keys(true);
        let pool = SqlitePoolOptions::new().connect_with(opts).await?;
        MIGRATOR.run(&pool).await?;
        Ok(Store { pool })
    }

    /// Open an in-memory database (one connection, so all ops share it) for tests.
    pub async fn open_in_memory() -> Result<Self, StoreError> {
        let opts = SqliteConnectOptions::from_str("sqlite::memory:")?.foreign_keys(true);
        let pool = SqlitePoolOptions::new()
            .max_connections(1)
            .connect_with(opts)
            .await?;
        MIGRATOR.run(&pool).await?;
        Ok(Store { pool })
    }

    pub fn pool(&self) -> &SqlitePool {
        &self.pool
    }
}
