//! Genotyping-panel queries: panels and their SNP sites.

use navigator_domain::workspace::{Panel, PanelSite};
use sqlx::SqlitePool;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct PanelRow {
    id: i64,
    name: String,
}

impl PanelRow {
    fn into_domain(self) -> Panel {
        Panel { id: self.id, name: self.name }
    }
}

#[derive(sqlx::FromRow)]
struct SiteRow {
    chrom: String,
    position: i64,
    reference_allele: String,
    alternate_allele: String,
    name: String,
}

impl SiteRow {
    fn into_domain(self) -> PanelSite {
        PanelSite {
            chrom: self.chrom,
            position: self.position,
            reference_allele: self.reference_allele,
            alternate_allele: self.alternate_allele,
            name: self.name,
        }
    }
}

/// Create a panel and bulk-insert its sites in one transaction.
pub async fn create(pool: &SqlitePool, name: &str, sites: &[PanelSite]) -> Result<Panel, StoreError> {
    let mut tx = pool.begin().await?;
    let id: i64 = sqlx::query_scalar("INSERT INTO panel (name) VALUES (?) RETURNING id")
        .bind(name)
        .fetch_one(&mut *tx)
        .await?;
    for s in sites {
        sqlx::query(
            "INSERT INTO panel_site (panel_id, chrom, position, reference_allele, alternate_allele, name) \
             VALUES (?, ?, ?, ?, ?, ?)",
        )
        .bind(id)
        .bind(&s.chrom)
        .bind(s.position)
        .bind(&s.reference_allele)
        .bind(&s.alternate_allele)
        .bind(&s.name)
        .execute(&mut *tx)
        .await?;
    }
    tx.commit().await?;
    Ok(Panel { id, name: name.to_string() })
}

pub async fn list(pool: &SqlitePool) -> Result<Vec<Panel>, StoreError> {
    let rows: Vec<PanelRow> = sqlx::query_as("SELECT id, name FROM panel ORDER BY id")
        .fetch_all(pool)
        .await?;
    Ok(rows.into_iter().map(PanelRow::into_domain).collect())
}

pub async fn sites(pool: &SqlitePool, panel_id: i64) -> Result<Vec<PanelSite>, StoreError> {
    let rows: Vec<SiteRow> = sqlx::query_as(
        "SELECT chrom, position, reference_allele, alternate_allele, name FROM panel_site \
         WHERE panel_id = ? ORDER BY chrom, position",
    )
    .bind(panel_id)
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(SiteRow::into_domain).collect())
}

pub async fn site_count(pool: &SqlitePool, panel_id: i64) -> Result<i64, StoreError> {
    let n: i64 = sqlx::query_scalar("SELECT COUNT(*) FROM panel_site WHERE panel_id = ?")
        .bind(panel_id)
        .fetch_one(pool)
        .await?;
    Ok(n)
}
