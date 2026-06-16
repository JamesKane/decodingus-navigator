//! Persisted consensus-profile snapshots — one per (biosample, DNA type). The full reconciled
//! profile is opaque JSON in `payload` (the app's `ConsensusProfile`); the scalar columns mirror the
//! header for quick listing. Upsert on rebuild, load on tab open. DNA-type-agnostic: `dna_type` keys
//! the row ('Y' today; 'Mt' / autosomal adapters reuse this table).

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A stored consensus-profile row: the scalar header + the opaque profile `payload`
/// (app `ConsensusProfile` JSON). `dna_type` is the lineage key ("Y" / "Mt" / …).
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredConsensusProfile {
    pub biosample_guid: String,
    pub dna_type: String,
    pub consensus_label: Option<String>,
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

/// Insert or replace the snapshot for a (biosample, dna_type).
pub async fn upsert(pool: &SqlitePool, p: &StoredConsensusProfile) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO consensus_profile (biosample_guid, dna_type, consensus_label, overall_confidence, \
         source_count, total, confirmed, novel, conflict, single_source, tree_provider, payload, last_reconciled_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, dna_type) DO UPDATE SET \
         consensus_label = excluded.consensus_label, overall_confidence = excluded.overall_confidence, \
         source_count = excluded.source_count, total = excluded.total, confirmed = excluded.confirmed, \
         novel = excluded.novel, conflict = excluded.conflict, single_source = excluded.single_source, \
         tree_provider = excluded.tree_provider, payload = excluded.payload, \
         last_reconciled_at = excluded.last_reconciled_at",
    )
    .bind(&p.biosample_guid)
    .bind(&p.dna_type)
    .bind(&p.consensus_label)
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

/// Load the snapshot for a (biosample, dna_type), if one has been built.
pub async fn get(pool: &SqlitePool, guid: SampleGuid, dna_type: &str) -> Result<Option<StoredConsensusProfile>, StoreError> {
    let row: Option<StoredConsensusProfile> =
        sqlx::query_as("SELECT * FROM consensus_profile WHERE biosample_guid = ? AND dna_type = ?")
            .bind(guid.0.to_string())
            .bind(dna_type)
            .fetch_optional(pool)
            .await?;
    Ok(row)
}

/// All persisted `(biosample_guid, dna_type, consensus_label)` rows that carry a non-null label —
/// the genome-level placed terminal per subject + arm, for the subjects-list bulk path (one query
/// instead of a `get` per subject).
pub async fn list_labels(pool: &SqlitePool) -> Result<Vec<(String, String, String)>, StoreError> {
    let rows: Vec<(String, String, String)> = sqlx::query_as(
        "SELECT biosample_guid, dna_type, consensus_label FROM consensus_profile \
         WHERE consensus_label IS NOT NULL AND consensus_label <> ''",
    )
    .fetch_all(pool)
    .await?;
    Ok(rows)
}

/// Remove a biosample's snapshot for one DNA type. Returns whether a row was removed.
pub async fn delete(pool: &SqlitePool, guid: SampleGuid, dna_type: &str) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM consensus_profile WHERE biosample_guid = ? AND dna_type = ?")
        .bind(guid.0.to_string())
        .bind(dna_type)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    fn sample(guid: &str) -> StoredConsensusProfile {
        StoredConsensusProfile {
            biosample_guid: guid.to_string(),
            dna_type: "Y".into(),
            consensus_label: Some("R-FGC29071".into()),
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

        assert!(get(pool.pool(), g, "Y").await.unwrap().is_none());
        upsert(pool.pool(), &sample(&g.0.to_string())).await.unwrap();
        let got = get(pool.pool(), g, "Y").await.unwrap().unwrap();
        assert_eq!(got.consensus_label.as_deref(), Some("R-FGC29071"));
        assert_eq!(got.confirmed, 8);
        // A different DNA type is a distinct row (no collision).
        assert!(get(pool.pool(), g, "Mt").await.unwrap().is_none());
        // An autosomal ('Auto') snapshot coexists with the Y row for the same biosample.
        let mut auto = sample(&g.0.to_string());
        auto.dna_type = "Auto".into();
        auto.consensus_label = None; // autosomes have no lineage label
        auto.confirmed = 42;
        upsert(pool.pool(), &auto).await.unwrap();
        assert_eq!(get(pool.pool(), g, "Auto").await.unwrap().unwrap().confirmed, 42);
        assert_eq!(get(pool.pool(), g, "Y").await.unwrap().unwrap().confirmed, 8); // Y row untouched
        // Upsert replaces (same biosample + dna_type).
        let mut updated = sample(&g.0.to_string());
        updated.confirmed = 9;
        upsert(pool.pool(), &updated).await.unwrap();
        assert_eq!(get(pool.pool(), g, "Y").await.unwrap().unwrap().confirmed, 9);
        assert!(delete(pool.pool(), g, "Y").await.unwrap());
        assert!(get(pool.pool(), g, "Y").await.unwrap().is_none());
        assert!(get(pool.pool(), g, "Auto").await.unwrap().is_some()); // deleting Y leaves Auto
    }
}
