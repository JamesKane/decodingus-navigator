//! Signature-keyed result caches — one row per biosample, holding an opaque JSON result plus the
//! signature of the input it was computed from.
//!
//! Four tables share this exact shape, and the app uses them the same way every time: read the row,
//! compare its signature against what the current inputs hash to, and recompute on a mismatch. They
//! were four hand-copied modules until this one replaced them; the copies had already drifted (the
//! two purge paths in the app each forgot a different table), which is the argument for having one.
//!
//! The columns are *named* differently per table for historical reasons — `consensus_sig` vs
//! `source_sig`, `roh` vs `archaic` vs `segments`, `computed_at` vs `painted_at` — so each cache
//! carries its column names and `get` aliases them back to the common [`Cached`] shape. The schema
//! is untouched; only the Rust side is unified.

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A cached result: the signature of the inputs it came from, the result itself as opaque JSON,
/// and when it was computed. The caller compares [`Cached::sig`] against the current inputs to
/// decide whether the payload is still good.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct Cached {
    pub biosample_guid: String,
    pub sig: String,
    pub payload: String,
    pub computed_at: String,
}

/// One signature-keyed cache table, identified by its name and its three non-key columns.
#[derive(Debug, Clone, Copy)]
pub struct SigCache {
    /// The table name. A compile-time constant in every case — never caller input — so
    /// interpolating it into SQL is safe.
    table: &'static str,
    sig_col: &'static str,
    payload_col: &'static str,
    at_col: &'static str,
}

/// Cached chromosome painting (local-ancestry segments), keyed to the autosomal consensus's
/// `last_reconciled_at`.
pub const PAINTING: SigCache = SigCache::new("consensus_painting", "consensus_sig", "segments", "painted_at");

/// Cached runs-of-homozygosity result (segments + summary), keyed to the autosomal consensus.
pub const ROH: SigCache = SigCache::new("consensus_roh", "consensus_sig", "roh", "computed_at");

/// Cached archaic (Neanderthal / Denisovan) **Tier A** marker count, keyed to the autosomal
/// consensus.
pub const ARCHAIC: SigCache = SigCache::new("consensus_archaic", "consensus_sig", "archaic", "computed_at");

/// Cached archaic **Tier B** segment calls. Keyed to the *alignment* they were called from rather
/// than to the consensus: segments come from genome-wide de-novo diploid calls on one alignment,
/// whereas the consensus only carries the 1240k panel loci. The signature is the alignment id plus
/// the caller's genotype version, so re-calling with a newer caller invalidates the cache.
pub const ARCHAIC_SEGMENTS: SigCache =
    SigCache::new("consensus_archaic_segments", "source_sig", "segments", "computed_at");

/// Every signature-keyed cache, in one list — so a purge that means "drop this subject's derived
/// results" drops *all* of them. Both purge paths used to enumerate tables by hand and both had
/// fallen behind the set.
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

    /// The table this cache lives in — for callers that must fold it into a wider delete inside
    /// their own transaction, where [`SigCache::delete`]'s pool-level call would not enlist.
    pub const fn table(&self) -> &'static str {
        self.table
    }

    /// Insert or replace this biosample's cached result.
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

    /// This biosample's cached result, if one exists. The caller checks [`Cached::sig`] for
    /// staleness — a row here is not by itself a usable result.
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

    /// Remove this biosample's cached result. `false` means there was nothing to remove.
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

    /// Every cache round-trips, and upsert replaces rather than duplicating (a recompute after the
    /// inputs changed). Running the same body over [`ALL`] is what keeps a newly added table from
    /// silently going untested.
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

    /// The caches are independent: writing one must not disturb another that happens to share a
    /// column name (`segments` is `consensus_painting`'s *and* `consensus_archaic_segments`'s).
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
