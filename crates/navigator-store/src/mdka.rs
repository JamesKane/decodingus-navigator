//! Most Distant Known Ancestor (FTDNA project-import design §4.3). One row per Subject per lineage,
//! upserted on `(biosample_guid, lineage)`.
//!
//! PII / project-shared-private: the most sensitive data in the importer. See the migration
//! `0030_mdka` header — never federated, never stored in AppView.

use du_domain::ids::SampleGuid;
use navigator_domain::identity::{Mdka, NewMdka};
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    biosample_guid: String,
    lineage: String,
    ancestor_name: Option<String>,
    birth_year: Option<i64>,
    death_year: Option<i64>,
    origin_place: Option<String>,
    origin_country: Option<String>,
    latitude: Option<f64>,
    longitude: Option<f64>,
    source: Option<String>,
    notes: Option<String>,
    updated_at: String,
}

impl Row {
    fn into_domain(self) -> Result<Mdka, StoreError> {
        Ok(Mdka {
            id: self.id,
            biosample_guid: parse_sample_guid(&self.biosample_guid, "mdka")?,
            lineage: self.lineage,
            ancestor_name: self.ancestor_name,
            birth_year: self.birth_year.map(|y| y as i32),
            death_year: self.death_year.map(|y| y as i32),
            origin_place: self.origin_place,
            origin_country: self.origin_country,
            latitude: self.latitude,
            longitude: self.longitude,
            source: self.source,
            notes: self.notes,
            updated_at: self.updated_at,
        })
    }
}

const COLS: &str = "id, biosample_guid, lineage, ancestor_name, birth_year, death_year, origin_place, \
                    origin_country, latitude, longitude, source, notes, updated_at";

/// Insert or replace the MDKA for a Subject's lineage. `updated_at` is the caller's ISO-8601 stamp.
pub async fn upsert(pool: &SqlitePool, guid: SampleGuid, m: &NewMdka, updated_at: &str) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO mdka \
         (biosample_guid, lineage, ancestor_name, birth_year, death_year, origin_place, origin_country, \
          latitude, longitude, source, notes, updated_at) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, lineage) DO UPDATE SET \
           ancestor_name = excluded.ancestor_name, \
           birth_year = excluded.birth_year, \
           death_year = excluded.death_year, \
           origin_place = excluded.origin_place, \
           origin_country = excluded.origin_country, \
           latitude = excluded.latitude, \
           longitude = excluded.longitude, \
           source = excluded.source, \
           notes = excluded.notes, \
           updated_at = excluded.updated_at",
    )
    .bind(guid.0.to_string())
    .bind(&m.lineage)
    .bind(&m.ancestor_name)
    .bind(m.birth_year)
    .bind(m.death_year)
    .bind(&m.origin_place)
    .bind(&m.origin_country)
    .bind(m.latitude)
    .bind(m.longitude)
    .bind(&m.source)
    .bind(&m.notes)
    .bind(updated_at)
    .execute(pool)
    .await?;
    Ok(())
}

pub async fn list_for(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<Mdka>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM mdka WHERE biosample_guid = ? ORDER BY lineage"
    ))
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_domain::identity::Lineage;
    use navigator_domain::workspace::Biosample;

    async fn seed(pool: &SqlitePool) -> SampleGuid {
        let guid = SampleGuid(uuid::Uuid::new_v4());
        crate::biosample::create(
            pool,
            &Biosample {
                guid,
                sample_accession: None,
                donor_identifier: "GFX".into(),
                description: None,
                center_name: None,
                sex: None,
                project_id: None,
            },
        )
        .await
        .unwrap();
        guid
    }

    #[tokio::test]
    async fn one_per_lineage_upsert() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        let g = seed(pool).await;
        let t = "2026-06-19T00:00:00Z";

        // B5163's paternal MDKA (the real fixture row).
        upsert(
            pool,
            g,
            &NewMdka {
                lineage: Lineage::Y.as_str().into(),
                ancestor_name: Some("Thomas Michael Kane".into()),
                birth_year: Some(1830),
                death_year: Some(1908),
                origin_place: Some("Creegh South, Co. Clare, Ireland".into()),
                origin_country: Some("Ireland".into()),
                latitude: Some(52.75),
                longitude: Some(-9.43),
                source: Some("FTDNA".into()),
                notes: None,
            },
            t,
        )
        .await
        .unwrap();
        upsert(
            pool,
            g,
            &NewMdka {
                lineage: Lineage::Mt.as_str().into(),
                ancestor_name: Some("Maternal line".into()),
                ..Default::default()
            },
            t,
        )
        .await
        .unwrap();

        let rows = list_for(pool, g).await.unwrap();
        assert_eq!(rows.len(), 2, "one per lineage");
        let y = rows.iter().find(|m| m.lineage == "Y").unwrap();
        assert_eq!(y.ancestor_name.as_deref(), Some("Thomas Michael Kane"));
        assert_eq!(y.birth_year, Some(1830));
        assert_eq!(y.latitude, Some(52.75));

        // Re-upsert the Y lineage replaces in place (no duplicate row).
        upsert(
            pool,
            g,
            &NewMdka {
                lineage: Lineage::Y.as_str().into(),
                ancestor_name: Some("Thomas M. Kane".into()),
                origin_country: Some("Ireland".into()),
                ..Default::default()
            },
            t,
        )
        .await
        .unwrap();
        let rows = list_for(pool, g).await.unwrap();
        assert_eq!(rows.len(), 2);
        let y = rows.iter().find(|m| m.lineage == "Y").unwrap();
        assert_eq!(y.ancestor_name.as_deref(), Some("Thomas M. Kane"));
        assert_eq!(y.birth_year, None, "replaced, not merged");
    }
}
