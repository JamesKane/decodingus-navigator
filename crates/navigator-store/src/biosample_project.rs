//! M:N Subject↔Project membership (FTDNA project-import design §4.1). `biosample_project` is the
//! source of truth for which projects a Subject belongs to; `biosample.project_id` lingers as a
//! nullable "home project" pointer.

use du_domain::ids::SampleGuid;
use navigator_domain::identity::ProjectMembership;
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    biosample_guid: String,
    project_id: i64,
    role: Option<String>,
    added_at: String,
}

impl Row {
    fn into_domain(self) -> Result<ProjectMembership, StoreError> {
        Ok(ProjectMembership {
            biosample_guid: parse_sample_guid(&self.biosample_guid, "biosample_project")?,
            project_id: self.project_id,
            role: self.role,
            added_at: self.added_at,
        })
    }
}

/// Add a membership (idempotent on the `(guid, project_id)` PK). On re-add, the `role` is updated.
pub async fn add(
    pool: &SqlitePool,
    guid: SampleGuid,
    project_id: i64,
    role: Option<&str>,
    added_at: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO biosample_project (biosample_guid, project_id, role, added_at) VALUES (?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, project_id) DO UPDATE SET role = excluded.role",
    )
    .bind(guid.0.to_string())
    .bind(project_id)
    .bind(role)
    .bind(added_at)
    .execute(pool)
    .await?;
    Ok(())
}

/// Remove a membership. Returns whether a row was removed.
pub async fn remove(pool: &SqlitePool, guid: SampleGuid, project_id: i64) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM biosample_project WHERE biosample_guid = ? AND project_id = ?")
        .bind(guid.0.to_string())
        .bind(project_id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Detach every subject from a project (used when deleting the project — the subjects themselves
/// are first-class and survive). Returns the number of memberships removed.
pub async fn remove_all_for_project(pool: &SqlitePool, project_id: i64) -> Result<u64, StoreError> {
    let affected = sqlx::query("DELETE FROM biosample_project WHERE project_id = ?")
        .bind(project_id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected)
}

/// Project ids a Subject belongs to.
pub async fn list_projects_for(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<i64>, StoreError> {
    let ids: Vec<i64> =
        sqlx::query_scalar("SELECT project_id FROM biosample_project WHERE biosample_guid = ? ORDER BY project_id")
            .bind(guid.0.to_string())
            .fetch_all(pool)
            .await?;
    Ok(ids)
}

/// Subject guids that belong to a project (via the membership table).
pub async fn list_biosamples_for(pool: &SqlitePool, project_id: i64) -> Result<Vec<SampleGuid>, StoreError> {
    let rows: Vec<String> =
        sqlx::query_scalar("SELECT biosample_guid FROM biosample_project WHERE project_id = ? ORDER BY biosample_guid")
            .bind(project_id)
            .fetch_all(pool)
            .await?;
    rows.iter().map(|g| parse_sample_guid(g, "biosample_project")).collect()
}

/// Full membership rows for a Subject (with role + timestamp).
pub async fn memberships_for(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<ProjectMembership>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(
        "SELECT biosample_guid, project_id, role, added_at FROM biosample_project \
         WHERE biosample_guid = ? ORDER BY project_id",
    )
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;
    rows.into_iter().map(Row::into_domain).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_domain::workspace::{Biosample, NewProject};

    async fn seed_project(pool: &SqlitePool, name: &str) -> i64 {
        crate::project::create(
            pool,
            &NewProject {
                name: name.into(),
                description: None,
                administrator: "admin".into(),
            },
        )
        .await
        .unwrap()
        .id
    }

    async fn seed_biosample(pool: &SqlitePool, donor: &str) -> SampleGuid {
        let guid = SampleGuid(uuid::Uuid::new_v4());
        crate::biosample::create(
            pool,
            &Biosample {
                guid,
                sample_accession: None,
                donor_identifier: donor.into(),
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
    async fn membership_is_mn_and_idempotent() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        let t = "2026-06-19T00:00:00Z";
        let p1 = seed_project(pool, "P1").await;
        let p2 = seed_project(pool, "P2").await;
        let s = seed_biosample(pool, "kitA").await;

        add(pool, s, p1, Some("subgroupX"), t).await.unwrap();
        add(pool, s, p2, None, t).await.unwrap();
        // Re-add the same pair updates role, doesn't duplicate.
        add(pool, s, p1, Some("subgroupY"), t).await.unwrap();

        let projects = list_projects_for(pool, s).await.unwrap();
        assert_eq!(projects, vec![p1, p2]);
        let members = memberships_for(pool, s).await.unwrap();
        assert_eq!(members.len(), 2);
        assert_eq!(members[0].role.as_deref(), Some("subgroupY"));

        assert_eq!(list_biosamples_for(pool, p1).await.unwrap(), vec![s]);
        assert!(remove(pool, s, p2).await.unwrap());
        assert_eq!(list_projects_for(pool, s).await.unwrap(), vec![p1]);
    }
}
