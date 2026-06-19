//! Per-source Y/mtDNA haplogroup calls — the inputs to donor-level reconciliation. One
//! row per (biosample, dna_type, source); upsert replaces a re-run from the same source.

use du_domain::ids::SampleGuid;
use navigator_domain::reconciliation::{DnaType, RunHaplogroupCall};
use sqlx::SqlitePool;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct Row {
    source_label: String,
    haplogroup: String,
    lineage: String,
    score: f64,
    matched: i64,
    expected: i64,
}

impl Row {
    fn into_domain(self) -> RunHaplogroupCall {
        let lineage = if self.lineage.is_empty() {
            Vec::new()
        } else {
            self.lineage.split('\t').map(str::to_string).collect()
        };
        RunHaplogroupCall {
            source_label: self.source_label,
            haplogroup: self.haplogroup,
            lineage,
            score: self.score,
            matched: self.matched,
            expected: self.expected,
        }
    }
}

/// Insert or replace the call from `source_key` for this biosample + DNA type. `fingerprint`
/// stamps the inputs (file + tree content hashes) so a later run can skip re-scoring.
pub async fn upsert(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    source_key: &str,
    call: &RunHaplogroupCall,
    fingerprint: Option<&str>,
) -> Result<(), StoreError> {
    let lineage = call.lineage.join("\t");
    sqlx::query(
        "INSERT INTO haplogroup_call \
         (biosample_guid, dna_type, source_key, source_label, haplogroup, lineage, score, matched, expected, source_fingerprint) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, dna_type, source_key) DO UPDATE SET \
         source_label = excluded.source_label, haplogroup = excluded.haplogroup, \
         lineage = excluded.lineage, score = excluded.score, matched = excluded.matched, \
         expected = excluded.expected, source_fingerprint = excluded.source_fingerprint",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .bind(source_key)
    .bind(&call.source_label)
    .bind(&call.haplogroup)
    .bind(lineage)
    .bind(call.score)
    .bind(call.matched)
    .bind(call.expected)
    .bind(fingerprint)
    .execute(pool)
    .await?;
    Ok(())
}

/// The stored input-fingerprint for one source's call, if recorded. Used to decide whether a
/// re-score is needed (the inputs are unchanged when this matches the current fingerprint).
pub async fn stored_fingerprint(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    source_key: &str,
) -> Result<Option<String>, StoreError> {
    let fp: Option<Option<String>> = sqlx::query_scalar(
        "SELECT source_fingerprint FROM haplogroup_call \
         WHERE biosample_guid = ? AND dna_type = ? AND source_key = ?",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .bind(source_key)
    .fetch_optional(pool)
    .await?;
    Ok(fp.flatten())
}

/// One source's recorded call (for returning a cached result without re-scoring).
pub async fn get_one(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    source_key: &str,
) -> Result<Option<RunHaplogroupCall>, StoreError> {
    let row: Option<Row> = sqlx::query_as(
        "SELECT source_label, haplogroup, lineage, score, matched, expected FROM haplogroup_call \
         WHERE biosample_guid = ? AND dna_type = ? AND source_key = ?",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .bind(source_key)
    .fetch_optional(pool)
    .await?;
    Ok(row.map(Row::into_domain))
}

/// All recorded calls for a biosample + DNA type.
pub async fn list_for(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
) -> Result<Vec<RunHaplogroupCall>, StoreError> {
    let rows: Vec<Row> = sqlx::query_as(
        "SELECT source_label, haplogroup, lineage, score, matched, expected FROM haplogroup_call \
         WHERE biosample_guid = ? AND dna_type = ? ORDER BY id",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(Row::into_domain).collect())
}

#[derive(sqlx::FromRow)]
struct AllRow {
    biosample_guid: String,
    dna_type: String,
    source_label: String,
    haplogroup: String,
    lineage: String,
    score: f64,
    matched: i64,
    expected: i64,
}

/// Every recorded call across all subjects, as `(guid, dna_type, call)` — for building a
/// donor-level haplogroup summary (the subjects list) in one query.
pub async fn list_all(pool: &SqlitePool) -> Result<Vec<(SampleGuid, DnaType, RunHaplogroupCall)>, StoreError> {
    let rows: Vec<AllRow> = sqlx::query_as(
        "SELECT biosample_guid, dna_type, source_label, haplogroup, lineage, score, matched, expected \
         FROM haplogroup_call ORDER BY id",
    )
    .fetch_all(pool)
    .await?;
    let mut out = Vec::with_capacity(rows.len());
    for r in rows {
        let guid = crate::error::parse_sample_guid(&r.biosample_guid, "haplogroup_call")?;
        let dna_type = match r.dna_type.as_str() {
            "Y" => DnaType::Y,
            "Mt" => DnaType::Mt,
            other => return Err(StoreError::Decode(format!("haplogroup_call dna_type {other:?}"))),
        };
        let lineage = if r.lineage.is_empty() {
            Vec::new()
        } else {
            r.lineage.split('\t').map(str::to_string).collect()
        };
        out.push((
            guid,
            dna_type,
            RunHaplogroupCall {
                source_label: r.source_label,
                haplogroup: r.haplogroup,
                lineage,
                score: r.score,
                matched: r.matched,
                expected: r.expected,
            },
        ));
    }
    Ok(out)
}
