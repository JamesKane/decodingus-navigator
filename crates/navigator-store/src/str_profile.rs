//! Y-STR profile queries: the profiles and their marker values. A profile attaches to a biosample.
//! The `SampleGuid` goes in as its hyphenated TEXT form, as it does everywhere else.

use du_domain::ids::SampleGuid;
use navigator_domain::strprofile::{NewStrProfile, StrMarker, StrProfile};
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
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
        StrMarker {
            marker: self.marker,
            value: self.value,
        }
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

/// The first STR profile of a biosample that matches `panel_name`, with its markers, if there is
/// one. The import uses it to merge a panel that arrives a second time, such as a Big Y CUSTOM set,
/// into the profile that exists. Without it the import would add a duplicate.
pub async fn find_by_panel(
    pool: &SqlitePool,
    guid: SampleGuid,
    panel_name: &str,
) -> Result<Option<StrProfile>, StoreError> {
    let row: Option<ProfileRow> = sqlx::query_as(
        "SELECT id, biosample_guid, panel_name, provider, source FROM str_profile \
         WHERE biosample_guid = ? AND panel_name = ? ORDER BY id LIMIT 1",
    )
    .bind(guid.0.to_string())
    .bind(panel_name)
    .fetch_optional(pool)
    .await?;
    let Some(r) = row else { return Ok(None) };
    let biosample_guid = parse_sample_guid(&r.biosample_guid, "str_profile")?;
    let markers = markers_for(pool, r.id).await?;
    Ok(Some(StrProfile {
        id: r.id,
        biosample_guid,
        panel_name: r.panel_name,
        provider: r.provider,
        source: r.source,
        markers,
    }))
}

/// Replace every marker of a profile, with a delete and then an insert, in one transaction. The
/// merge of a panel that arrives a second time into an existing profile uses it.
pub async fn replace_markers(pool: &SqlitePool, profile_id: i64, markers: &[StrMarker]) -> Result<(), StoreError> {
    let mut tx = pool.begin().await?;
    sqlx::query("DELETE FROM str_marker WHERE str_profile_id = ?")
        .bind(profile_id)
        .execute(&mut *tx)
        .await?;
    for m in markers {
        sqlx::query("INSERT INTO str_marker (str_profile_id, marker, value) VALUES (?, ?, ?)")
            .bind(profile_id)
            .bind(&m.marker)
            .bind(&m.value)
            .execute(&mut *tx)
            .await?;
    }
    tx.commit().await?;
    Ok(())
}

/// Delete an STR profile and its markers. The children go first, because the database enforces the
/// FKs. It returns whether it removed the profile row.
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
        let biosample_guid = parse_sample_guid(&r.biosample_guid, "str_profile")?;
        let markers = markers_for(pool, r.id).await?;
        profiles.push(StrProfile {
            id: r.id,
            biosample_guid,
            panel_name: r.panel_name,
            provider: r.provider,
            source: r.source,
            markers,
        });
    }
    Ok(profiles)
}
