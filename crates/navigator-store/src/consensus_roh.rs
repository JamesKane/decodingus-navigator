//! Cached runs-of-homozygosity (ROH) result per subject, keyed to the autosomal consensus signature
//! it was computed from. The app upserts on (re)compute and reads back when the signature still
//! matches the current consensus — otherwise it recomputes. One row per biosample. Mirrors
//! [`crate::consensus_painting`].

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A stored ROH result: the consensus signature it was computed from + the full result (opaque JSON).
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredRoh {
    pub biosample_guid: String,
    pub consensus_sig: String,
    pub roh: String,
    pub computed_at: String,
}

/// Insert or replace the cached ROH result for a biosample.
pub async fn upsert(
    pool: &SqlitePool,
    guid: SampleGuid,
    consensus_sig: &str,
    roh: &str,
    computed_at: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO consensus_roh (biosample_guid, consensus_sig, roh, computed_at) \
         VALUES (?, ?, ?, ?) \
         ON CONFLICT(biosample_guid) DO UPDATE SET \
         consensus_sig = excluded.consensus_sig, roh = excluded.roh, computed_at = excluded.computed_at",
    )
    .bind(guid.0.to_string())
    .bind(consensus_sig)
    .bind(roh)
    .bind(computed_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// The cached ROH result for a biosample, if one exists (caller checks the signature for staleness).
pub async fn get(pool: &SqlitePool, guid: SampleGuid) -> Result<Option<StoredRoh>, StoreError> {
    let row: Option<StoredRoh> = sqlx::query_as("SELECT * FROM consensus_roh WHERE biosample_guid = ?")
        .bind(guid.0.to_string())
        .fetch_optional(pool)
        .await?;
    Ok(row)
}

/// Remove a biosample's cached ROH result.
pub async fn delete(pool: &SqlitePool, guid: SampleGuid) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM consensus_roh WHERE biosample_guid = ?")
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
        assert_eq!(got.consensus_sig, "2026-07-22T00:00:00Z");
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
        assert_eq!(got.consensus_sig, "2026-07-23T00:00:00Z");
        assert_eq!(got.roh, r#"{"segments":[]}"#);
        assert!(delete(pool.pool(), g).await.unwrap());
        assert!(get(pool.pool(), g).await.unwrap().is_none());
    }
}
