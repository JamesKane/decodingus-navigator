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

/// Insert or replace the call from `source_key` for this biosample + DNA type.
pub async fn upsert(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    source_key: &str,
    call: &RunHaplogroupCall,
) -> Result<(), StoreError> {
    let lineage = call.lineage.join("\t");
    sqlx::query(
        "INSERT INTO haplogroup_call \
         (biosample_guid, dna_type, source_key, source_label, haplogroup, lineage, score, matched, expected) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, dna_type, source_key) DO UPDATE SET \
         source_label = excluded.source_label, haplogroup = excluded.haplogroup, \
         lineage = excluded.lineage, score = excluded.score, matched = excluded.matched, \
         expected = excluded.expected",
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
    .execute(pool)
    .await?;
    Ok(())
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
