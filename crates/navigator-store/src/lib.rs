//! Navigator local persistence: SQLite through `sqlx`, in place of the H2/Slick layer.
//!
//! Each aggregate has one query module over a `SqlitePool`. A complex child is a real row, and the
//! migrations carry a version and can go backwards. The stored state has authority.

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
pub mod consensus_profile;
pub mod dm;
pub mod error;
pub mod external_id;
pub mod external_panel_dosage;
pub mod ftdna_member;
pub mod haplogroup_call;
pub mod ibd_exchange;
pub mod ibd_request;
pub mod mdka;
pub mod mtdna;
pub mod project;
pub mod reconciliation;
pub mod sequence_run;
pub mod sig_cache;
pub mod source_file;
pub mod str_profile;
pub mod sync_history;
pub mod sync_outbox;
pub mod sync_state;
pub mod variant_set;
pub mod variant_set_genotype;
pub mod variant_set_private_y;

pub use error::StoreError;

static MIGRATOR: sqlx::migrate::Migrator = sqlx::migrate!("./migrations");

/// Handle to the workspace database.
#[derive(Clone, Debug)]
pub struct Store {
    pool: SqlitePool,
}

impl Store {
    /// Open the SQLite database at `path`, make it when it is absent, and run the migrations.
    pub async fn open(path: &Path) -> Result<Self, StoreError> {
        let opts = SqliteConnectOptions::new()
            .filename(path)
            .create_if_missing(true)
            .foreign_keys(true)
            // WAL, with a long busy timeout. The GUI and a CLI can then run at the same time, such
            // as `navigator analyze`, and share the one workspace file. Neither one then sees an
            // immediate `database is locked` failure.
            //
            // WAL lets a reader run beside the one writer. The timeout makes a second writer wait
            // for the first to finish, and it gives no error.
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
