//! Persisted Y-profile snapshots (one per biosample). The full reconciled profile is opaque JSON in
//! `payload` (the app's `YProfile`); the scalar columns mirror the header for quick listing. Upsert
//! on rebuild, load on tab open.

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A stored Y-profile row: the scalar header + the opaque profile `payload` (app `YProfile` JSON).
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredYProfile {
    pub biosample_guid: String,
    pub consensus_haplogroup: Option<String>,
    pub overall_confidence: f64,
    pub source_count: i64,
    pub total: i64,
    pub confirmed: i64,
    pub novel: i64,
    pub conflict: i64,
    pub single_source: i64,
    pub tree_provider: Option<String>,
    pub payload: String,
    pub last_reconciled_at: String,
}

/// Insert or replace the snapshot for a biosample.
pub async fn upsert(pool: &SqlitePool, guid: SampleGuid, p: &StoredYProfile) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO y_profile (biosample_guid, consensus_haplogroup, overall_confidence, source_count, \
         total, confirmed, novel, conflict, single_source, tree_provider, payload, last_reconciled_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid) DO UPDATE SET \
         consensus_haplogroup = excluded.consensus_haplogroup, overall_confidence = excluded.overall_confidence, \
         source_count = excluded.source_count, total = excluded.total, confirmed = excluded.confirmed, \
         novel = excluded.novel, conflict = excluded.conflict, single_source = excluded.single_source, \
         tree_provider = excluded.tree_provider, payload = excluded.payload, \
         last_reconciled_at = excluded.last_reconciled_at",
    )
    .bind(guid.0.to_string())
    .bind(&p.consensus_haplogroup)
    .bind(p.overall_confidence)
    .bind(p.source_count)
    .bind(p.total)
    .bind(p.confirmed)
    .bind(p.novel)
    .bind(p.conflict)
    .bind(p.single_source)
    .bind(&p.tree_provider)
    .bind(&p.payload)
    .bind(&p.last_reconciled_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// Load the snapshot for a biosample, if one has been built.
pub async fn get(pool: &SqlitePool, guid: SampleGuid) -> Result<Option<StoredYProfile>, StoreError> {
    let row: Option<StoredYProfile> = sqlx::query_as("SELECT * FROM y_profile WHERE biosample_guid = ?")
        .bind(guid.0.to_string())
        .fetch_optional(pool)
        .await?;
    Ok(row)
}

/// Remove a biosample's snapshot (e.g. on subject delete). Returns whether a row was removed.
pub async fn delete(pool: &SqlitePool, guid: SampleGuid) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM y_profile WHERE biosample_guid = ?")
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

    fn sample(guid: &str) -> StoredYProfile {
        StoredYProfile {
            biosample_guid: guid.to_string(),
            consensus_haplogroup: Some("R-FGC29071".into()),
            overall_confidence: 0.9,
            source_count: 2,
            total: 10,
            confirmed: 8,
            novel: 1,
            conflict: 1,
            single_source: 0,
            tree_provider: Some("decodingus".into()),
            payload: r#"{"variants":[],"summary":{},"terminal":"R-FGC29071","sources":[]}"#.into(),
            last_reconciled_at: "2026-06-13T00:00:00Z".into(),
        }
    }

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
        upsert(pool.pool(), g, &sample(&g.0.to_string())).await.unwrap();
        let got = get(pool.pool(), g).await.unwrap().unwrap();
        assert_eq!(got.consensus_haplogroup.as_deref(), Some("R-FGC29071"));
        assert_eq!(got.confirmed, 8);
        // Upsert replaces.
        let mut updated = sample(&g.0.to_string());
        updated.confirmed = 9;
        upsert(pool.pool(), g, &updated).await.unwrap();
        assert_eq!(get(pool.pool(), g).await.unwrap().unwrap().confirmed, 9);
        assert!(delete(pool.pool(), g).await.unwrap());
        assert!(get(pool.pool(), g).await.unwrap().is_none());
    }
}
