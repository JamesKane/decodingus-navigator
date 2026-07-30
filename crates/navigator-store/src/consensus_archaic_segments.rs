//! Cached **Tier B** archaic SEGMENT calls per subject.
//!
//! Keyed to the alignment they were called from rather than to the autosomal consensus: segments
//! come from genome-wide de-novo diploid calls on one alignment, whereas the consensus only carries
//! the 1240k panel loci. `source_sig` is the alignment id plus the caller's genotype version, so
//! re-calling with a newer caller invalidates the cache. Mirrors [`crate::consensus_archaic`].

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A stored segments marker count result: the consensus signature it was computed from + the full result (opaque JSON).
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredArchaicSegments {
    pub biosample_guid: String,
    pub source_sig: String,
    pub segments: String,
    pub computed_at: String,
}

/// Insert or replace the cached segments marker count result for a biosample.
pub async fn upsert(
    pool: &SqlitePool,
    guid: SampleGuid,
    source_sig: &str,
    segments: &str,
    computed_at: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO consensus_archaic_segments (biosample_guid, source_sig, segments, computed_at) \
         VALUES (?, ?, ?, ?) \
         ON CONFLICT(biosample_guid) DO UPDATE SET \
         source_sig = excluded.source_sig, segments = excluded.segments, computed_at = excluded.computed_at",
    )
    .bind(guid.0.to_string())
    .bind(source_sig)
    .bind(segments)
    .bind(computed_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// The cached segments marker count result for a biosample, if one exists (caller checks the signature for staleness).
pub async fn get(pool: &SqlitePool, guid: SampleGuid) -> Result<Option<StoredArchaicSegments>, StoreError> {
    let row: Option<StoredArchaicSegments> = sqlx::query_as("SELECT * FROM consensus_archaic_segments WHERE biosample_guid = ?")
        .bind(guid.0.to_string())
        .fetch_optional(pool)
        .await?;
    Ok(row)
}

/// Remove a biosample's cached segments marker count result.
pub async fn delete(pool: &SqlitePool, guid: SampleGuid) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM consensus_archaic_segments WHERE biosample_guid = ?")
        .bind(guid.0.to_string())
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    #[tokio::test]
    async fn upsert_get_delete_round_trip() {
        let pool = crate::Store::open_in_memory().await.unwrap();
        let g = SampleGuid(Uuid::new_v4());
        let bio = navigator_domain::workspace::Biosample {
            guid: g,
            sample_accession: None,
            donor_identifier: "S1".into(),
            description: None,
            center_name: None,
            sex: None,
            project_id: None,
        };
        crate::biosample::create(pool.pool(), &bio).await.unwrap();

        assert!(get(pool.pool(), g).await.unwrap().is_none());
        upsert(pool.pool(), g, "2026-07-22T00:00:00Z", "{}", "2026-07-22T01:00:00Z")
            .await
            .unwrap();
        let got = get(pool.pool(), g).await.unwrap().unwrap();
        assert_eq!(got.source_sig, "2026-07-22T00:00:00Z");
        // Upsert replaces (a recompute after a consensus rebuild).
        upsert(
            pool.pool(),
            g,
            "2026-07-23T00:00:00Z",
            r#"{"segments":[]}"#,
            "2026-07-23T01:00:00Z",
        )
        .await
        .unwrap();
        let got = get(pool.pool(), g).await.unwrap().unwrap();
        assert_eq!(got.source_sig, "2026-07-23T00:00:00Z");
        assert_eq!(got.segments, r#"{"segments":[]}"#);
        assert!(delete(pool.pool(), g).await.unwrap());
        assert!(get(pool.pool(), g).await.unwrap().is_none());
    }
}
