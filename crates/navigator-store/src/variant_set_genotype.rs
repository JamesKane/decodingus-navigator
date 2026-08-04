//! Cached tree-position genotypes for a variant set — the VCF counterpart of the per-alignment
//! `tree-genotype` analysis artifact.
//!
//! Keyed by `(variant_set_id, cache_key)`, where `cache_key` carries a hash of the target positions
//! exactly as the alignment path's `algorithm_version` does: a changed tree yields a different key,
//! so a stale genotype misses rather than quietly placing the donor against sites that moved.

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

/// Drop every cached genotype for a set — used when the set's calls are replaced.
pub async fn delete_for_set(pool: &SqlitePool, set_id: i64) -> Result<u64, StoreError> {
    Ok(sqlx::query("DELETE FROM variant_set_genotype WHERE variant_set_id = ?")
        .bind(set_id)
        .execute(pool)
        .await?
        .rows_affected())
}
