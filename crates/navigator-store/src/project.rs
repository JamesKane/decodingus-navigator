//! Project queries.

use navigator_domain::workspace::{NewProject, Project};
use sqlx::SqlitePool;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    id: i64,
    name: String,
    description: Option<String>,
    administrator: String,
}

impl Row {
    fn into_domain(self) -> Project {
        Project {
            id: self.id,
            name: self.name,
            description: self.description,
            administrator: self.administrator,
        }
    }
}

pub async fn create(pool: &SqlitePool, p: &NewProject) -> Result<Project, StoreError> {
    let id: i64 =
        sqlx::query_scalar("INSERT INTO project (name, description, administrator) VALUES (?, ?, ?) RETURNING id")
            .bind(&p.name)
            .bind(&p.description)
            .bind(&p.administrator)
            .fetch_one(pool)
            .await?;
    Ok(Project {
        id,
        name: p.name.clone(),
        description: p.description.clone(),
        administrator: p.administrator.clone(),
    })
}

/// Update the fields of a project that a user can edit. It returns whether it changed a row.
pub async fn update(
    pool: &SqlitePool,
    id: i64,
    name: &str,
    description: Option<&str>,
    administrator: &str,
) -> Result<bool, StoreError> {
    let affected = sqlx::query("UPDATE project SET name = ?, description = ?, administrator = ? WHERE id = ?")
        .bind(name)
        .bind(description)
        .bind(administrator)
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

/// Delete a project row. It returns whether it removed a row. A caller must make sure that no
/// biosample still names the project, because the database enforces the FKs. The app layer guards
/// that.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let affected = sqlx::query("DELETE FROM project WHERE id = ?")
        .bind(id)
        .execute(pool)
        .await?
        .rows_affected();
    Ok(affected > 0)
}

pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<Project>, StoreError> {
    let row: Option<Row> = sqlx::query_as("SELECT id, name, description, administrator FROM project WHERE id = ?")
        .bind(id)
        .fetch_optional(pool)
        .await?;
    Ok(row.map(Row::into_domain))
}

pub async fn list(pool: &SqlitePool) -> Result<Vec<Project>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as("SELECT id, name, description, administrator FROM project ORDER BY id")
        .fetch_all(pool)
        .await?;
    Ok(rows.into_iter().map(Row::into_domain).collect())
}
