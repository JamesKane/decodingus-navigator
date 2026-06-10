//! Y-STR profile queries: profiles and their marker values. Profiles attach to a
//! biosample (the `SampleGuid` is stored as its hyphenated TEXT form, like elsewhere).

use du_domain::ids::SampleGuid;
use navigator_domain::strprofile::{NewStrProfile, StrMarker, StrProfile};
use sqlx::SqlitePool;
use uuid::Uuid;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct ProfileRow {
    id: i64,
    biosample_guid: String,
    panel_name: String,
    provider: Option<String>,
    source: Option<String>,
}

#[derive(sqlx::FromRow)]
struct MarkerRow {
    marker: String,
    value: String,
}

impl MarkerRow {
    fn into_domain(self) -> StrMarker {
        StrMarker { marker: self.marker, value: self.value }
    }
}

/// Create an STR profile and bulk-insert its markers in one transaction.
pub async fn create(pool: &SqlitePool, new: &NewStrProfile) -> Result<StrProfile, StoreError> {
    let mut tx = pool.begin().await?;
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO str_profile (biosample_guid, panel_name, provider, source) \
         VALUES (?, ?, ?, ?) RETURNING id",
    )
    .bind(new.biosample_guid.0.to_string())
    .bind(&new.panel_name)
    .bind(&new.provider)
    .bind(&new.source)
    .fetch_one(&mut *tx)
    .await?;
    for m in &new.markers {
        sqlx::query("INSERT INTO str_marker (str_profile_id, marker, value) VALUES (?, ?, ?)")
            .bind(id)
            .bind(&m.marker)
            .bind(&m.value)
            .execute(&mut *tx)
            .await?;
    }
    tx.commit().await?;
    Ok(StrProfile {
        id,
        biosample_guid: new.biosample_guid,
        panel_name: new.panel_name.clone(),
        provider: new.provider.clone(),
        source: new.source.clone(),
        markers: new.markers.clone(),
    })
}

async fn markers_for(pool: &SqlitePool, profile_id: i64) -> Result<Vec<StrMarker>, StoreError> {
    let rows: Vec<MarkerRow> =
        sqlx::query_as("SELECT marker, value FROM str_marker WHERE str_profile_id = ? ORDER BY id")
            .bind(profile_id)
            .fetch_all(pool)
            .await?;
    Ok(rows.into_iter().map(MarkerRow::into_domain).collect())
}

/// Delete an STR profile and its markers (children-first; FKs are enforced). Returns whether
/// the profile row was removed.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let mut tx = pool.begin().await?;
    sqlx::query("DELETE FROM str_marker WHERE str_profile_id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?;
    let affected = sqlx::query("DELETE FROM str_profile WHERE id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?
        .rows_affected();
    tx.commit().await?;
    Ok(affected > 0)
}

/// All STR profiles for a biosample, with their markers.
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<StrProfile>, StoreError> {
    let rows: Vec<ProfileRow> = sqlx::query_as(
        "SELECT id, biosample_guid, panel_name, provider, source FROM str_profile \
         WHERE biosample_guid = ? ORDER BY id",
    )
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;

    let mut profiles = Vec::with_capacity(rows.len());
    for r in rows {
        let uuid = Uuid::parse_str(&r.biosample_guid)
            .map_err(|e| StoreError::Decode(format!("str_profile guid {:?}: {e}", r.biosample_guid)))?;
        let markers = markers_for(pool, r.id).await?;
        profiles.push(StrProfile {
            id: r.id,
            biosample_guid: SampleGuid(uuid),
            panel_name: r.panel_name,
            provider: r.provider,
            source: r.source,
            markers,
        });
    }
    Ok(profiles)
}
