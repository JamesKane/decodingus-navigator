//! Persisted federated-IBD exchange results (gap §4) — one row per completed exchange session. The
//! match summary scalars are columns for quick listing; the segment list + both signed attestations
//! are opaque JSON. Upsert by `session_id` (re-running an exchange overwrites). PII-free: DIDs +
//! opaque sample refs + cM only.

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A stored IBD-exchange result row. `segments` / `my_attestation` / `partner_attestation` are opaque
/// JSON (the app's `IbdSegment[]` / `IbdAttestation`); the rest mirror the match header for listing.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredIbdExchange {
    pub session_id: String,
    pub request_uri: String,
    pub my_did: String,
    pub partner_did: String,
    pub biosample_guid: String,
    pub partner_sample_ref: Option<String>,
    pub total_shared_cm: f64,
    pub segment_count: i64,
    pub longest_segment_cm: f64,
    pub relationship: String,
    pub agreed: bool,
    pub segments: String,
    pub my_attestation: String,
    pub partner_attestation: String,
    pub created_at: String,
}

/// Insert or replace the result for a session.
pub async fn upsert(pool: &SqlitePool, r: &StoredIbdExchange) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO ibd_exchange_result (session_id, request_uri, my_did, partner_did, biosample_guid, \
         partner_sample_ref, total_shared_cm, segment_count, longest_segment_cm, relationship, agreed, \
         segments, my_attestation, partner_attestation, created_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(session_id) DO UPDATE SET \
         request_uri = excluded.request_uri, my_did = excluded.my_did, partner_did = excluded.partner_did, \
         biosample_guid = excluded.biosample_guid, partner_sample_ref = excluded.partner_sample_ref, \
         total_shared_cm = excluded.total_shared_cm, segment_count = excluded.segment_count, \
         longest_segment_cm = excluded.longest_segment_cm, relationship = excluded.relationship, \
         agreed = excluded.agreed, segments = excluded.segments, my_attestation = excluded.my_attestation, \
         partner_attestation = excluded.partner_attestation, created_at = excluded.created_at",
    )
    .bind(&r.session_id)
    .bind(&r.request_uri)
    .bind(&r.my_did)
    .bind(&r.partner_did)
    .bind(&r.biosample_guid)
    .bind(&r.partner_sample_ref)
    .bind(r.total_shared_cm)
    .bind(r.segment_count)
    .bind(r.longest_segment_cm)
    .bind(&r.relationship)
    .bind(r.agreed)
    .bind(&r.segments)
    .bind(&r.my_attestation)
    .bind(&r.partner_attestation)
    .bind(&r.created_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// All exchange results, newest first.
pub async fn list(pool: &SqlitePool) -> Result<Vec<StoredIbdExchange>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM ibd_exchange_result ORDER BY created_at DESC")
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

/// Exchange results for one local subject, newest first.
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<StoredIbdExchange>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM ibd_exchange_result WHERE biosample_guid = ? ORDER BY created_at DESC")
        .bind(guid.0.to_string())
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    fn row(session: &str, guid: &str, cm: f64) -> StoredIbdExchange {
        StoredIbdExchange {
            session_id: session.into(),
            request_uri: "exchange:r1".into(),
            my_did: "did:key:zA".into(),
            partner_did: "did:key:zB".into(),
            biosample_guid: guid.into(),
            partner_sample_ref: Some("bio-B".into()),
            total_shared_cm: cm,
            segment_count: 1,
            longest_segment_cm: cm,
            relationship: "ThirdCousin".into(),
            agreed: true,
            segments: "[]".into(),
            my_attestation: "{}".into(),
            partner_attestation: "{}".into(),
            created_at: "2026-06-17T00:00:00Z".into(),
        }
    }

    #[tokio::test]
    async fn upsert_and_list_round_trip() {
        let store = crate::Store::open_in_memory().await.unwrap();
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
        crate::biosample::create(store.pool(), &bio).await.unwrap();

        assert!(list_for_biosample(store.pool(), g).await.unwrap().is_empty());
        upsert(store.pool(), &row("sess-1", &g.0.to_string(), 75.0))
            .await
            .unwrap();
        let got = list_for_biosample(store.pool(), g).await.unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].total_shared_cm, 75.0);
        assert!(got[0].agreed);
        // Upsert on the same session replaces.
        upsert(store.pool(), &row("sess-1", &g.0.to_string(), 40.0))
            .await
            .unwrap();
        let got = list_for_biosample(store.pool(), g).await.unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].total_shared_cm, 40.0);
        assert_eq!(list(store.pool()).await.unwrap().len(), 1);
    }
}
