//! Cached tree-position genotypes for a variant set. This is the VCF counterpart of the
//! `tree-genotype` analysis artifact that each alignment has.
//!
//! The key is `(variant_set_id, cache_key)`. The `cache_key` holds a hash of the target positions,
//! exactly as the `algorithm_version` of the alignment path does. A changed tree gives a different
//! key, so a stale genotype misses the cache. Without that, the code would place the donor against
//! sites that moved, and say nothing.

use sqlx::SqlitePool;

use crate::StoreError;

/// Store (or replace) the genotypes for one `(set, cache key)`.
pub async fn upsert(pool: &SqlitePool, set_id: i64, cache_key: &str, calls_json: &str) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO variant_set_genotype (variant_set_id, cache_key, calls, computed_at) \
         VALUES (?, ?, ?, datetime('now')) \
         ON CONFLICT(variant_set_id, cache_key) DO UPDATE SET calls = excluded.calls, computed_at = excluded.computed_at",
    )
    .bind(set_id)
    .bind(cache_key)
    .bind(calls_json)
    .execute(pool)
    .await?;
    Ok(())
}

/// The cached genotypes, or `None` on a miss (absent set, or a different target site-set).
pub async fn get(pool: &SqlitePool, set_id: i64, cache_key: &str) -> Result<Option<String>, StoreError> {
    Ok(
        sqlx::query_scalar("SELECT calls FROM variant_set_genotype WHERE variant_set_id = ? AND cache_key = ?")
            .bind(set_id)
            .bind(cache_key)
            .fetch_optional(pool)
            .await?,
    )
}

/// Drop every cached genotype for a set. The code calls this when new calls replace the set's
/// calls.
pub async fn delete_for_set(pool: &SqlitePool, set_id: i64) -> Result<u64, StoreError> {
    Ok(sqlx::query("DELETE FROM variant_set_genotype WHERE variant_set_id = ?")
        .bind(set_id)
        .execute(pool)
        .await?
        .rows_affected())
}
