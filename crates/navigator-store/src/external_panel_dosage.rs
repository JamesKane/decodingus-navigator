//! Per-source autosomal 1240K panel dosages from a trusted external caller (a GATK4 / 1240K
//! EIGENSTRAT call set). One row per (biosample, source_label). `dosages` is opaque JSON (the app's
//! resolved `Vec<SiteGenotype>`, CHM13-oriented); the store just persists and lists it so the
//! autosomal consensus can pool it with no CRAM decode. See migration `0037`.

use du_domain::ids::SampleGuid;
use sqlx::SqlitePool;

use crate::StoreError;

/// A stored external panel-dosage row: provenance/label header + the opaque `dosages` payload.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct StoredPanelDosage {
    pub biosample_guid: String,
    pub source_label: String,
    pub provenance: String,
    pub panel_sig: Option<String>,
    pub site_count: i64,
    pub dosages: String,
    pub created_at: String,
}

/// Insert or replace the dosages for a (biosample, source_label).
pub async fn upsert(pool: &SqlitePool, row: &StoredPanelDosage) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO external_panel_dosage \
         (biosample_guid, source_label, provenance, panel_sig, site_count, dosages, created_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, source_label) DO UPDATE SET \
         provenance = excluded.provenance, panel_sig = excluded.panel_sig, \
         site_count = excluded.site_count, dosages = excluded.dosages, created_at = excluded.created_at",
    )
    .bind(&row.biosample_guid)
    .bind(&row.source_label)
    .bind(&row.provenance)
    .bind(&row.panel_sig)
    .bind(row.site_count)
    .bind(&row.dosages)
    .bind(&row.created_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// All external panel-dosage rows for a biosample (each a distinct source).
pub async fn list_for_biosample(
    pool: &SqlitePool,
    guid: SampleGuid,
) -> Result<Vec<StoredPanelDosage>, StoreError> {
    let rows: Vec<StoredPanelDosage> = sqlx::query_as(
        "SELECT biosample_guid, source_label, provenance, panel_sig, site_count, dosages, created_at \
         FROM external_panel_dosage WHERE biosample_guid = ? ORDER BY id",
    )
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;
    Ok(rows)
}

/// Remove all external panel-dosage rows for a biosample (e.g. clearing a subject's analysis).
/// Returns the number removed.
pub async fn delete_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<u64, StoreError> {
    let affected = sqlx::query("DELETE FROM external_panel_dosage WHERE biosample_guid = ?")
        .bind(guid.0.to_string())
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    #[tokio::test]
    async fn upsert_list_delete_round_trip() {
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
        let row = StoredPanelDosage {
            biosample_guid: g.0.to_string(),
            source_label: "aadr 1240K (EIGENSTRAT)".into(),
            provenance: "external".into(),
            panel_sig: Some("abc123".into()),
            site_count: 3,
            dosages: r#"[{"name":"rs1","dosage":2}]"#.into(),
            created_at: "2026-07-21T00:00:00Z".into(),
        };
        upsert(store.pool(), &row).await.unwrap();
        let got = list_for_biosample(store.pool(), g).await.unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].site_count, 3);
        assert_eq!(got[0].provenance, "external");

        // Upsert on the same (biosample, label) replaces, not duplicates.
        let mut row2 = row.clone();
        row2.site_count = 5;
        upsert(store.pool(), &row2).await.unwrap();
        let got = list_for_biosample(store.pool(), g).await.unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].site_count, 5);

        assert_eq!(delete_for_biosample(store.pool(), g).await.unwrap(), 1);
        assert!(list_for_biosample(store.pool(), g).await.unwrap().is_empty());
    }
}
