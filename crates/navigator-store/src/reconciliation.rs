//! Reconciliation persistence: the manual override of the consensus haplogroup and the
//! audit log of reconciliation actions, keyed by (biosample, DNA type).

use du_domain::ids::SampleGuid;
use navigator_domain::reconciliation::{AuditEntry, DnaType};
use sqlx::SqlitePool;

use crate::StoreError;

/// Set (upsert) a manual override of the consensus haplogroup.
pub async fn set_override(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    haplogroup: &str,
    reason: Option<&str>,
) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO reconciliation_override (biosample_guid, dna_type, haplogroup, reason) \
         VALUES (?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, dna_type) DO UPDATE SET haplogroup = excluded.haplogroup, reason = excluded.reason",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .bind(haplogroup)
    .bind(reason)
    .execute(pool)
    .await?;
    Ok(())
}

/// The manual override (haplogroup, reason) for a subject + DNA type, if any.
pub async fn get_override(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
) -> Result<Option<(String, Option<String>)>, StoreError> {
    let row: Option<(String, Option<String>)> =
        sqlx::query_as("SELECT haplogroup, reason FROM reconciliation_override WHERE biosample_guid = ? AND dna_type = ?")
            .bind(biosample_guid.0.to_string())
            .bind(dna_type.as_str())
            .fetch_optional(pool)
            .await?;
    Ok(row)
}

/// Every manual override as `(guid, dna_type, haplogroup)` — for the subjects-list summary.
pub async fn list_all_overrides(pool: &SqlitePool) -> Result<Vec<(SampleGuid, DnaType, String)>, StoreError> {
    let rows: Vec<(String, String, String)> =
        sqlx::query_as("SELECT biosample_guid, dna_type, haplogroup FROM reconciliation_override")
            .fetch_all(pool)
            .await?;
    let mut out = Vec::with_capacity(rows.len());
    for (g, dt, hg) in rows {
        let guid = uuid::Uuid::parse_str(&g)
            .map_err(|e| StoreError::Decode(format!("override guid {g:?}: {e}")))?;
        let dna_type = match dt.as_str() {
            "Y" => DnaType::Y,
            "Mt" => DnaType::Mt,
            other => return Err(StoreError::Decode(format!("override dna_type {other:?}"))),
        };
        out.push((SampleGuid(guid), dna_type, hg));
    }
    Ok(out)
}

/// Remove a manual override.
pub async fn clear_override(pool: &SqlitePool, biosample_guid: SampleGuid, dna_type: DnaType) -> Result<(), StoreError> {
    sqlx::query("DELETE FROM reconciliation_override WHERE biosample_guid = ? AND dna_type = ?")
        .bind(biosample_guid.0.to_string())
        .bind(dna_type.as_str())
        .execute(pool)
        .await?;
    Ok(())
}

/// Append an audit-log entry.
pub async fn append_audit(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    entry: &AuditEntry,
) -> Result<(), StoreError> {
    sqlx::query("INSERT INTO reconciliation_audit (biosample_guid, dna_type, ts, action, note) VALUES (?, ?, ?, ?, ?)")
        .bind(biosample_guid.0.to_string())
        .bind(dna_type.as_str())
        .bind(&entry.timestamp)
        .bind(&entry.action)
        .bind(&entry.note)
        .execute(pool)
        .await?;
    Ok(())
}

/// The audit log for a subject + DNA type, oldest first.
pub async fn list_audit(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
) -> Result<Vec<AuditEntry>, StoreError> {
    let rows: Vec<(String, String, String)> = sqlx::query_as(
        "SELECT ts, action, note FROM reconciliation_audit WHERE biosample_guid = ? AND dna_type = ? ORDER BY id",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(|(timestamp, action, note)| AuditEntry { timestamp, action, note }).collect())
}
