//! Navigator local persistence — SQLite via `sqlx`, replacing the H2/Slick layer.
//! Query-module-per-aggregate over a `SqlitePool`; complex children modelled as proper
//! rows, with versioned, reversible migrations. Persisted state is authoritative.

use std::path::Path;
use std::str::FromStr;

use sqlx::sqlite::{SqliteConnectOptions, SqlitePoolOptions};
use sqlx::SqlitePool;

pub mod alignment;
pub mod ancestry_result;
pub mod artifact;
pub mod biosample;
pub mod biosample_project;
pub mod chip_profile;
pub mod consensus_painting;
pub mod consensus_profile;
pub mod consensus_roh;
pub mod dm;
pub mod error;
pub mod external_id;
pub mod external_panel_dosage;
pub mod ftdna_member;
pub mod haplogroup_call;
pub mod ibd_exchange;
pub mod mdka;
pub mod mtdna;
pub mod project;
pub mod reconciliation;
pub mod sequence_run;
pub mod source_file;
pub mod str_profile;
pub mod sync_history;
pub mod sync_outbox;
pub mod sync_state;
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
            .foreign_keys(true)
            // WAL + a generous busy timeout so the GUI and a concurrent CLI (`navigator analyze`)
            // can share the one workspace file without immediate "database is locked" failures:
            // WAL lets readers run alongside a single writer, and the timeout makes a contended
            // writer wait for the other to finish instead of erroring out.
            .journal_mode(sqlx::sqlite::SqliteJournalMode::Wal)
            .synchronous(sqlx::sqlite::SqliteSynchronous::Normal)
            .busy_timeout(std::time::Duration::from_secs(30));
        let pool = SqlitePoolOptions::new().connect_with(opts).await?;
        MIGRATOR.run(&pool).await?;
        Ok(Store { pool })
    }

    /// Open an in-memory database (one connection, so all ops share it) for tests.
    pub async fn open_in_memory() -> Result<Self, StoreError> {
        let opts = SqliteConnectOptions::from_str("sqlite::memory:")?.foreign_keys(true);
        let pool = SqlitePoolOptions::new().max_connections(1).connect_with(opts).await?;
        MIGRATOR.run(&pool).await?;
        Ok(Store { pool })
    }

    pub fn pool(&self) -> &SqlitePool {
        &self.pool
    }
}
