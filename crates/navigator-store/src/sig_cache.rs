//! Result caches with a signature key. Each row holds one result for one biosample. The result is
//! JSON text. The row also holds the signature of the input data that made the result.
//!
//! Four tables have this shape. The app uses each of the four tables in the same sequence:
//!
//! 1. Read the row for the biosample.
//! 2. Compare the signature in the row with the signature of the current input data.
//! 3. If the two signatures are different, calculate the result again.
//!
//! Before this module, there were four modules with the same code. The four copies became
//! different. Each of the two purge paths in the app forgot a different table. One module prevents
//! this fault.
//!
//! Each table gives different names to its columns. One table has `consensus_sig` and another table
//! has `source_sig`. So each cache keeps its own column names, and `get` changes these names to the
//! names in the [`Cached`] structure. This module does not change the database schema. It changes
//! only the Rust code.

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// One cached result. `sig` is the signature of the input data. `payload` is the result as JSON
/// text. `computed_at` is the time of the calculation. The caller compares [`Cached::sig`] with the
/// signature of the current input data. If the two are different, the payload is out of date.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct Cached {
    pub biosample_guid: String,
    pub sig: String,
    pub payload: String,
    pub computed_at: String,
}

/// One cache table. The table name and the three columns that are not the key define it.
#[derive(Debug, Clone, Copy)]
pub struct SigCache {
    /// The table name. This value is always a constant in the code. The caller never supplies it.
    /// For this reason, it is safe to put the value into the SQL text.
    table: &'static str,
    sig_col: &'static str,
    payload_col: &'static str,
    at_col: &'static str,
}

/// The cached chromosome painting. The painting holds the local ancestry segments. The key is the
/// `last_reconciled_at` value of the autosomal consensus.
pub const PAINTING: SigCache = SigCache::new("consensus_painting", "consensus_sig", "segments", "painted_at");

/// The cached runs-of-homozygosity result. The result holds the segments and a summary. The key is
/// the autosomal consensus.
pub const ROH: SigCache = SigCache::new("consensus_roh", "consensus_sig", "roh", "computed_at");

/// The cached archaic **Tier A** marker count for Neanderthal and Denisovan. The key is the
/// autosomal consensus.
pub const ARCHAIC: SigCache = SigCache::new("consensus_archaic", "consensus_sig", "archaic", "computed_at");

/// The cached archaic **Tier B** segment calls. The key is the alignment, not the consensus. The
/// caller finds these segments from de-novo diploid calls on one alignment across the genome. The
/// consensus holds only the 1240k panel loci. The signature is the alignment id and the genotype
/// version of the caller. So a newer caller makes the cache out of date.
pub const ARCHAIC_SEGMENTS: SigCache =
    SigCache::new("consensus_archaic_segments", "source_sig", "segments", "computed_at");

/// All of the caches, in one list. A purge that must remove the derived results of a subject
/// removes all of them. Before this list, the two purge paths named the tables one by one. Each
/// path did not name all of the tables.
pub const ALL: [SigCache; 4] = [PAINTING, ROH, ARCHAIC, ARCHAIC_SEGMENTS];

impl SigCache {
    const fn new(table: &'static str, sig_col: &'static str, payload_col: &'static str, at_col: &'static str) -> Self {
        SigCache {
            table,
            sig_col,
            payload_col,
            at_col,
        }
    }

    /// The table of this cache. A caller needs the name when it deletes many tables in its own
    /// transaction. The [`SigCache::delete`] function uses the pool, so that function can not join
    /// such a transaction.
    pub const fn table(&self) -> &'static str {
        self.table
    }

    /// Insert or replace the cached result for this biosample.
    pub async fn upsert(
        &self,
        pool: &SqlitePool,
        guid: SampleGuid,
        sig: &str,
        payload: &str,
        computed_at: &str,
    ) -> Result<(), StoreError> {
        let (table, s, p, a) = (self.table, self.sig_col, self.payload_col, self.at_col);
        sqlx::query(&format!(
            "INSERT INTO {table} (biosample_guid, {s}, {p}, {a}) VALUES (?, ?, ?, ?) \
             ON CONFLICT(biosample_guid) DO UPDATE SET \
             {s} = excluded.{s}, {p} = excluded.{p}, {a} = excluded.{a}"
        ))
        .bind(guid.0.to_string())
        .bind(sig)
        .bind(payload)
        .bind(computed_at)
        .execute(pool)
        .await?;
        Ok(())
    }

    /// The cached result for this biosample, if a result exists. The caller must check
    /// [`Cached::sig`]. A row is not a usable result until that check passes.
    pub async fn get(&self, pool: &SqlitePool, guid: SampleGuid) -> Result<Option<Cached>, StoreError> {
        let (table, s, p, a) = (self.table, self.sig_col, self.payload_col, self.at_col);
        let row: Option<Cached> = sqlx::query_as(&format!(
            "SELECT biosample_guid, {s} AS sig, {p} AS payload, {a} AS computed_at \
             FROM {table} WHERE biosample_guid = ?"
        ))
        .bind(guid.0.to_string())
        .fetch_optional(pool)
        .await?;
        Ok(row)
    }

    /// Remove the cached result for this biosample. `false` shows that there was no result.
    pub async fn delete(&self, pool: &SqlitePool, guid: SampleGuid) -> Result<bool, StoreError> {
        let affected = sqlx::query(&format!("DELETE FROM {} WHERE biosample_guid = ?", self.table))
            .bind(guid.0.to_string())
            .execute(pool)
            .await?
            .rows_affected();
        Ok(affected > 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    /// Each cache keeps and returns the same data. `upsert` replaces a row and does not add a
    /// second row. This occurs when the input data changes and the app calculates the result again.
    /// The test uses [`ALL`], so a new table always gets a test.
    #[tokio::test]
    async fn every_cache_round_trips_and_upsert_replaces() {
        let pool = crate::Store::open_in_memory().await.unwrap();
        let g = SampleGuid(Uuid::new_v4());
        crate::biosample::create(pool.pool(), &navigator_domain::workspace::Biosample::new(g, "S1"))
            .await
            .unwrap();

        for cache in ALL {
            assert!(cache.get(pool.pool(), g).await.unwrap().is_none(), "{cache:?}");
            cache
                .upsert(pool.pool(), g, "sig-1", "{}", "2026-07-22T01:00:00Z")
                .await
                .unwrap();
            let got = cache.get(pool.pool(), g).await.unwrap().unwrap();
            assert_eq!(got.sig, "sig-1", "{cache:?}");
            assert_eq!(got.payload, "{}", "{cache:?}");

            cache
                .upsert(pool.pool(), g, "sig-2", r#"{"segments":[]}"#, "2026-07-23T01:00:00Z")
                .await
                .unwrap();
            let got = cache.get(pool.pool(), g).await.unwrap().unwrap();
            assert_eq!(got.sig, "sig-2", "{cache:?}");
            assert_eq!(got.payload, r#"{"segments":[]}"#, "{cache:?}");
            assert_eq!(got.computed_at, "2026-07-23T01:00:00Z", "{cache:?}");

            assert!(cache.delete(pool.pool(), g).await.unwrap(), "{cache:?}");
            assert!(cache.get(pool.pool(), g).await.unwrap().is_none(), "{cache:?}");
        }
    }

    /// The caches are independent. A write to one cache must not change a different cache. Two of
    /// the tables use the column name `segments`: `consensus_painting` and
    /// `consensus_archaic_segments`.
    #[tokio::test]
    async fn caches_do_not_alias_each_other() {
        let pool = crate::Store::open_in_memory().await.unwrap();
        let g = SampleGuid(Uuid::new_v4());
        crate::biosample::create(pool.pool(), &navigator_domain::workspace::Biosample::new(g, "S1"))
            .await
            .unwrap();

        PAINTING.upsert(pool.pool(), g, "sig-p", "painted", "t").await.unwrap();
        ARCHAIC_SEGMENTS
            .upsert(pool.pool(), g, "sig-a", "called", "t")
            .await
            .unwrap();

        assert_eq!(PAINTING.get(pool.pool(), g).await.unwrap().unwrap().payload, "painted");
        assert_eq!(
            ARCHAIC_SEGMENTS.get(pool.pool(), g).await.unwrap().unwrap().payload,
            "called"
        );
    }
}
