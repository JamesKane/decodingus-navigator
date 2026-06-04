//! Persisted ancestry estimates — one row per (biosample, alignment, panel). Upsert
//! replaces a re-run of the same panel on the same alignment. The ranked components and
//! super-population summary are stored as JSON blobs.

use du_domain::ids::SampleGuid;
use navigator_domain::ancestry::{AncestryResult, PopulationComponent, SuperPopulationSummary};
use sqlx::SqlitePool;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    alignment_id: i64,
    panel_type: String,
    reference_version: String,
    confidence_level: f64,
    snps_analyzed: i64,
    snps_with_genotype: i64,
    components_json: String,
    super_pop_json: String,
    pca_json: Option<String>,
}

impl Row {
    fn into_domain(self) -> Result<AncestryResult, StoreError> {
        let components: Vec<PopulationComponent> = serde_json::from_str(&self.components_json)?;
        let super_population_summary: Vec<SuperPopulationSummary> =
            serde_json::from_str(&self.super_pop_json)?;
        let pca_coordinates: Option<Vec<f64>> = match self.pca_json {
            Some(s) => serde_json::from_str(&s)?,
            None => None,
        };
        let snps_analyzed = self.snps_analyzed.max(0) as usize;
        let snps_with_genotype = self.snps_with_genotype.max(0) as usize;
        Ok(AncestryResult {
            panel_type: self.panel_type,
            snps_analyzed,
            snps_with_genotype,
            snps_missing: snps_analyzed.saturating_sub(snps_with_genotype),
            components,
            super_population_summary,
            confidence_level: self.confidence_level,
            pipeline_version: String::new(),
            reference_version: self.reference_version,
            pca_coordinates,
        })
    }
}

/// Insert or replace the ancestry estimate for `alignment_id` under this panel type.
pub async fn upsert(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    alignment_id: i64,
    result: &AncestryResult,
) -> Result<(), StoreError> {
    let components_json = serde_json::to_string(&result.components)?;
    let super_pop_json = serde_json::to_string(&result.super_population_summary)?;
    let pca_json = match &result.pca_coordinates {
        Some(c) => Some(serde_json::to_string(c)?),
        None => None,
    };
    sqlx::query(
        "INSERT INTO ancestry_result \
         (biosample_guid, alignment_id, panel_type, reference_version, confidence_level, \
          snps_analyzed, snps_with_genotype, components_json, super_pop_json, pca_json) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, alignment_id, panel_type) DO UPDATE SET \
         reference_version = excluded.reference_version, confidence_level = excluded.confidence_level, \
         snps_analyzed = excluded.snps_analyzed, snps_with_genotype = excluded.snps_with_genotype, \
         components_json = excluded.components_json, super_pop_json = excluded.super_pop_json, \
         pca_json = excluded.pca_json",
    )
    .bind(biosample_guid.0.to_string())
    .bind(alignment_id)
    .bind(&result.panel_type)
    .bind(&result.reference_version)
    .bind(result.confidence_level)
    .bind(result.snps_analyzed as i64)
    .bind(result.snps_with_genotype as i64)
    .bind(components_json)
    .bind(super_pop_json)
    .bind(pca_json)
    .execute(pool)
    .await?;
    Ok(())
}

/// The most recent ancestry estimate recorded for `alignment_id`, if any.
pub async fn get_for_alignment(
    pool: &SqlitePool,
    alignment_id: i64,
) -> Result<Option<AncestryResult>, StoreError> {
    let row: Option<Row> = sqlx::query_as(
        "SELECT alignment_id, panel_type, reference_version, confidence_level, snps_analyzed, \
         snps_with_genotype, components_json, super_pop_json, pca_json FROM ancestry_result \
         WHERE alignment_id = ? ORDER BY id DESC LIMIT 1",
    )
    .bind(alignment_id)
    .fetch_optional(pool)
    .await?;
    row.map(Row::into_domain).transpose()
}

/// All ancestry estimates recorded for a biosample, newest first.
pub async fn for_biosample(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
) -> Result<Vec<(i64, AncestryResult)>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(
        "SELECT alignment_id, panel_type, reference_version, confidence_level, snps_analyzed, \
         snps_with_genotype, components_json, super_pop_json, pca_json FROM ancestry_result \
         WHERE biosample_guid = ? ORDER BY id DESC",
    )
    .bind(biosample_guid.0.to_string())
    .fetch_all(pool)
    .await?;
    rows.into_iter()
        .map(|r| {
            let aln = r.alignment_id;
            r.into_domain().map(|d| (aln, d))
        })
        .collect()
}
