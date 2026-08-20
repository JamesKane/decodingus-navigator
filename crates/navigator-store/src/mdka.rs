//! Most Distant Known Ancestor (FTDNA project-import design §4.3). Each Subject has one row for
//! each lineage, and the upsert key is `(biosample_guid, lineage)`.
//!
//! This is PII, and a project shares it in private. It is the most sensitive data in the importer.
//! See the header of the migration `0030_mdka`. It never goes to the federation, and the AppView
//! never stores it.

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

/// Remove a Subject's MDKA for one lineage. It returns `true` when it removed a row.
pub async fn delete(pool: &SqlitePool, guid: SampleGuid, lineage: &str) -> Result<bool, StoreError> {
    let res = sqlx::query("DELETE FROM mdka WHERE biosample_guid = ? AND lineage = ?")
        .bind(guid.0.to_string())
        .bind(lineage)
        .execute(pool)
        .await?;
    Ok(res.rows_affected() > 0)
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

/// MDKA rows for one lineage that the workspace may publish.
///
/// The predicate is the whole consent story, and it is deliberately narrow:
///
/// - **the workspace holds primary data for the subject**, which is a `variant_set`, or a
///   `sequence_run` with an `alignment`. A roster row on its own is another person's kit that we
///   happen to know of. To publish its genealogy would publish again what the tester gave to a
///   vendor, and not to us.
/// - **and the tester has not opted out.** `ftdna_member.publicly_shares` is the member's own FTDNA
///   share setting, which arrives with the roster. A `0` there is an explicit "do not show me", and
///   it wins over everything else. With no roster row there is no opt-out to obey.
///
/// Measured on the reference workspace: 583 Y rows sit on subjects with primary data, and 558 of
/// those share publicly. The other 25 must never publish.
pub async fn publishable(pool: &SqlitePool, lineage: &str) -> Result<Vec<Mdka>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(&format!(
        "SELECT {COLS} FROM mdka m \
         WHERE m.lineage = ? \
           AND ( EXISTS (SELECT 1 FROM variant_set v WHERE v.biosample_guid = m.biosample_guid) \
              OR EXISTS (SELECT 1 FROM sequence_run r JOIN alignment a ON a.sequence_run_id = r.id \
                          WHERE r.biosample_guid = m.biosample_guid) ) \
           AND NOT EXISTS (SELECT 1 FROM ftdna_member f \
                            WHERE f.biosample_guid = m.biosample_guid \
                              AND COALESCE(f.publicly_shares, 0) = 0) \
         ORDER BY m.biosample_guid"
    ))
    .bind(lineage)
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

/// How many MDKA rows a lineage has, whether the app may publish them or not. It is the
/// denominator that makes a count of the publishable rows readable.
pub async fn count_for_lineage(pool: &SqlitePool, lineage: &str) -> Result<usize, StoreError> {
    let n: i64 = sqlx::query_scalar("SELECT count(*) FROM mdka WHERE lineage = ?")
        .bind(lineage)
        .fetch_one(pool)
        .await?;
    Ok(n as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_domain::identity::Lineage;
    use navigator_domain::workspace::Biosample;

    async fn seed(pool: &SqlitePool) -> SampleGuid {
        let guid = SampleGuid(uuid::Uuid::new_v4());
        crate::biosample::create(pool, &Biosample::new(guid, "GFX"))
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

        // Delete removes only the named lineage; a second delete is a no-op.
        assert!(delete(pool, g, "Y").await.unwrap(), "Y removed");
        assert!(!delete(pool, g, "Y").await.unwrap(), "second delete is a no-op");
        let rows = list_for(pool, g).await.unwrap();
        assert_eq!(rows.len(), 1, "only Mt remains");
        assert_eq!(rows[0].lineage, "Mt");
    }

    /// The consent predicate is the whole story of what may leave the workspace, so this test pins
    /// it case by case. Against the reference workspace it selects 558 of 583 Y rows. The 25 that
    /// it drops are testers who told FTDNA not to share them publicly.
    #[tokio::test]
    async fn publishable_requires_primary_data_and_no_opt_out() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        let t = "2026-08-07T00:00:00Z";

        let mk = |name: &str| NewMdka {
            lineage: Lineage::Y.as_str().into(),
            ancestor_name: Some(name.into()),
            origin_country: Some("Ireland".into()),
            ..Default::default()
        };

        // (a) primary data, and no roster row, so there is nothing to opt out of.
        let owned = seed(pool).await;
        upsert(pool, owned, &mk("Thomas Kane"), t).await.unwrap();
        sqlx::query("INSERT INTO variant_set (biosample_guid, source_label, source_type, reference_build, call_schema, source_path) VALUES (?,?,?,?,?,?)")
            .bind(owned.0.to_string()).bind("v").bind("VCF").bind("hs1").bind("s").bind("/tmp/v.vcf")
            .execute(pool).await.unwrap();

        // (b) primary data, and the tester shares publicly.
        let sharing = seed(pool).await;
        upsert(pool, sharing, &mk("Mary Sullivan"), t).await.unwrap();
        sqlx::query("INSERT INTO variant_set (biosample_guid, source_label, source_type, reference_build, call_schema, source_path) VALUES (?,?,?,?,?,?)")
            .bind(sharing.0.to_string()).bind("v").bind("VCF").bind("hs1").bind("s").bind("/tmp/w.vcf")
            .execute(pool).await.unwrap();
        sqlx::query("INSERT INTO ftdna_member (biosample_guid, member_name, publicly_shares) VALUES (?,?,1)")
            .bind(sharing.0.to_string())
            .bind("M S")
            .execute(pool)
            .await
            .unwrap();

        // (c) primary data, but the tester opted OUT. Must never publish.
        let opted_out = seed(pool).await;
        upsert(pool, opted_out, &mk("Patrick Walsh"), t).await.unwrap();
        sqlx::query("INSERT INTO variant_set (biosample_guid, source_label, source_type, reference_build, call_schema, source_path) VALUES (?,?,?,?,?,?)")
            .bind(opted_out.0.to_string()).bind("v").bind("VCF").bind("hs1").bind("s").bind("/tmp/x.vcf")
            .execute(pool).await.unwrap();
        sqlx::query("INSERT INTO ftdna_member (biosample_guid, member_name, publicly_shares) VALUES (?,?,0)")
            .bind(opted_out.0.to_string())
            .bind("P W")
            .execute(pool)
            .await
            .unwrap();

        // (d) a roster row and an MDKA, but the workspace holds no data of its own. This is
        // another person's kit that we only know of.
        let roster_only = seed(pool).await;
        upsert(pool, roster_only, &mk("Bridget Moore"), t).await.unwrap();
        sqlx::query("INSERT INTO ftdna_member (biosample_guid, member_name, publicly_shares) VALUES (?,?,1)")
            .bind(roster_only.0.to_string())
            .bind("B M")
            .execute(pool)
            .await
            .unwrap();

        let out = publishable(pool, "Y").await.unwrap();
        let guids: Vec<_> = out.iter().map(|m| m.biosample_guid).collect();
        assert!(guids.contains(&owned), "primary data, no roster row");
        assert!(guids.contains(&sharing), "primary data, sharing publicly");
        assert!(!guids.contains(&opted_out), "an explicit opt-out wins over everything");
        assert!(
            !guids.contains(&roster_only),
            "a roster row alone is not ours to publish"
        );
        assert_eq!(out.len(), 2);

        // The query obeys the lineage: an Mt row on a publishable subject is not a Y row.
        assert!(publishable(pool, "Mt").await.unwrap().is_empty());
    }
}
