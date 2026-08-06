//! Cached private-Y buckets for a variant set — the VCF counterpart of the per-alignment `private_y`
//! artifact. Same `(set, cache_key)` shape as [`crate::variant_set_genotype`], so a changed tree is a
//! miss rather than a bucket classified against sites that have moved.

use sqlx::SqlitePool;

use crate::StoreError;

pub async fn upsert(pool: &SqlitePool, set_id: i64, cache_key: &str, bucket_json: &str) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO variant_set_private_y (variant_set_id, cache_key, bucket, computed_at) \
         VALUES (?, ?, ?, datetime('now')) \
         ON CONFLICT(variant_set_id, cache_key) DO UPDATE SET bucket = excluded.bucket, computed_at = excluded.computed_at",
    )
    .bind(set_id)
    .bind(cache_key)
    .bind(bucket_json)
    .execute(pool)
    .await?;
    Ok(())
}

pub async fn get(pool: &SqlitePool, set_id: i64, cache_key: &str) -> Result<Option<String>, StoreError> {
    Ok(
        sqlx::query_scalar("SELECT bucket FROM variant_set_private_y WHERE variant_set_id = ? AND cache_key = ?")
            .bind(set_id)
            .bind(cache_key)
            .fetch_optional(pool)
            .await?,
    )
}

/// Every cached bucket for a subject's sets, newest cache first — for the donor-level union.
pub async fn list_for_biosample(pool: &SqlitePool, guid: &str) -> Result<Vec<(i64, String)>, StoreError> {
    Ok(sqlx::query_as(
        "SELECT p.variant_set_id, p.bucket FROM variant_set_private_y p \
         JOIN variant_set vs ON vs.id = p.variant_set_id \
         WHERE vs.biosample_guid = ? ORDER BY p.variant_set_id",
    )
    .bind(guid)
    .fetch_all(pool)
    .await?)
}
