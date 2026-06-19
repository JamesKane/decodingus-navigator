//! FTDNA-reported member labels (FTDNA project-import design §4.2). One row per Subject, upserted on
//! the `biosample_guid` PK. These are FTDNA's *reported* labels — our own computed haplogroups live
//! in `haplogroup_call` (different provenance).
//!
//! PII / never-federated: see the migration `0029_subject_identity` header.

use du_domain::ids::SampleGuid;
use navigator_domain::identity::FtdnaMember;
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    biosample_guid: String,
    member_name: Option<String>,
    y_haplogroup_ftdna: Option<String>,
    mt_haplogroup_ftdna: Option<String>,
    haplo_status: Option<String>,
    access_granted: Option<String>,
    publicly_shares: Option<i64>,
}

impl Row {
    fn into_domain(self) -> Result<FtdnaMember, StoreError> {
        Ok(FtdnaMember {
            biosample_guid: parse_sample_guid(&self.biosample_guid, "ftdna_member")?,
            member_name: self.member_name,
            y_haplogroup_ftdna: self.y_haplogroup_ftdna,
            mt_haplogroup_ftdna: self.mt_haplogroup_ftdna,
            haplo_status: self.haplo_status,
            access_granted: self.access_granted,
            publicly_shares: self.publicly_shares.map(|n| n != 0),
        })
    }
}

const COLS: &str = "biosample_guid, member_name, y_haplogroup_ftdna, mt_haplogroup_ftdna, haplo_status, \
                    access_granted, publicly_shares";

/// Insert or update the FTDNA member labels for a Subject.
pub async fn upsert(pool: &SqlitePool, m: &FtdnaMember) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO ftdna_member \
         (biosample_guid, member_name, y_haplogroup_ftdna, mt_haplogroup_ftdna, haplo_status, access_granted, publicly_shares) \
         VALUES (?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid) DO UPDATE SET \
           member_name = excluded.member_name, \
           y_haplogroup_ftdna = excluded.y_haplogroup_ftdna, \
           mt_haplogroup_ftdna = excluded.mt_haplogroup_ftdna, \
           haplo_status = excluded.haplo_status, \
           access_granted = excluded.access_granted, \
           publicly_shares = excluded.publicly_shares",
    )
    .bind(m.biosample_guid.0.to_string())
    .bind(&m.member_name)
    .bind(&m.y_haplogroup_ftdna)
    .bind(&m.mt_haplogroup_ftdna)
    .bind(&m.haplo_status)
    .bind(&m.access_granted)
    .bind(m.publicly_shares.map(|b| b as i64))
    .execute(pool)
    .await?;
    Ok(())
}

pub async fn get(pool: &SqlitePool, guid: SampleGuid) -> Result<Option<FtdnaMember>, StoreError> {
    let row: Option<Row> = sqlx::query_as(&format!("SELECT {COLS} FROM ftdna_member WHERE biosample_guid = ?"))
        .bind(guid.0.to_string())
        .fetch_optional(pool)
        .await?;
    row.map(Row::into_domain).transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
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
    async fn upsert_replaces_and_roundtrips_bool() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        let g = seed(pool).await;

        upsert(
            pool,
            &FtdnaMember {
                biosample_guid: g,
                member_name: Some("REDACTED".into()),
                y_haplogroup_ftdna: Some("R-FGC29071".into()),
                mt_haplogroup_ftdna: None,
                haplo_status: Some("confirmed".into()),
                access_granted: Some("Limited".into()),
                publicly_shares: Some(true),
            },
        )
        .await
        .unwrap();
        let got = get(pool, g).await.unwrap().unwrap();
        assert_eq!(got.y_haplogroup_ftdna.as_deref(), Some("R-FGC29071"));
        assert_eq!(got.publicly_shares, Some(true));

        // Upsert overwrites in place (PK conflict).
        upsert(
            pool,
            &FtdnaMember {
                biosample_guid: g,
                member_name: Some("REDACTED".into()),
                y_haplogroup_ftdna: Some("R-FGC29071".into()),
                mt_haplogroup_ftdna: Some("U5a1b1g".into()),
                haplo_status: Some("confirmed".into()),
                access_granted: Some("Advanced".into()),
                publicly_shares: Some(false),
            },
        )
        .await
        .unwrap();
        let got = get(pool, g).await.unwrap().unwrap();
        assert_eq!(got.mt_haplogroup_ftdna.as_deref(), Some("U5a1b1g"));
        assert_eq!(got.access_granted.as_deref(), Some("Advanced"));
        assert_eq!(got.publicly_shares, Some(false));
    }
}
