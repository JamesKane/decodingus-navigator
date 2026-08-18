//! Per-source Y/mtDNA haplogroup calls — the inputs to donor-level reconciliation. One
//! row per (biosample, dna_type, source); upsert replaces a re-run from the same source.

use du_domain::ids::SampleGuid;
use navigator_domain::reconciliation::{CallProvenance, DnaType, RunHaplogroupCall};
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

/// Insert or replace the call from `source_key` for this biosample + DNA type. `provenance` records
/// which caller produced it (external / navigator-walk / manual) — the precedence tier used at
/// reconcile. `fingerprint` stamps the inputs (file + tree content hashes) so a later run can skip
/// re-scoring. External and internal calls use *distinct* `source_key`s, so this upsert never lets
/// one overwrite the other.
pub async fn upsert(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    source_key: &str,
    call: &RunHaplogroupCall,
    provenance: CallProvenance,
    fingerprint: Option<&str>,
) -> Result<(), StoreError> {
    let lineage = call.lineage.join("\t");
    sqlx::query(
        "INSERT INTO haplogroup_call \
         (biosample_guid, dna_type, source_key, source_label, haplogroup, lineage, score, matched, expected, provenance, source_fingerprint) \
         VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) \
         ON CONFLICT(biosample_guid, dna_type, source_key) DO UPDATE SET \
         source_label = excluded.source_label, haplogroup = excluded.haplogroup, \
         lineage = excluded.lineage, score = excluded.score, matched = excluded.matched, \
         expected = excluded.expected, provenance = excluded.provenance, \
         source_fingerprint = excluded.source_fingerprint",
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
    .bind(provenance.as_str())
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

/// Subjects whose recorded calls were placed against **a different haplotree** than the one now
/// active — the selector for a workspace re-evaluation after a new tree lands.
///
/// A fingerprint is `<source tag>:<hash>|<tree tag>:<hash>` (`f:…|yt:…` for an alignment-placed Y
/// call, `gv:…|mt:…` for a GVCF-placed mt call), so `tree_tag` is `"yt:"` or `"mt:"` and `tree_hash`
/// is the first 16 hex of the current tree's SHA-256.
///
/// `include_unknown` decides what to do with a call carrying **no tree tag** — placed before the
/// fingerprint existed, so which tree it used can not be established. These are a different job from
/// a tree change, and much more expensive: in this workspace 3,780 Y calls are *provably* on a
/// superseded tree while 15,648 have no fingerprint at all, and most of the latter own a BAM that a
/// re-placement would re-walk. Default them out so the routine "a new tree landed" sweep stays
/// runnable, and let a caller opt into the backfill explicitly.
pub async fn biosamples_placed_against_another_tree(
    pool: &SqlitePool,
    dna_type: DnaType,
    tree_tag: &str,
    tree_hash: &str,
    include_unknown: bool,
) -> Result<Vec<SampleGuid>, StoreError> {
    // `instr(...) = 0` catches a fingerprint that carries no tree tag at all; without it `substr`
    // would slice from offset 1 and could, in principle, compare equal to the hash. Whether that
    // case *selects* is `include_unknown`; it must never fall through to the hash comparison.
    let unknown = "(source_fingerprint IS NULL OR source_fingerprint = '' OR instr(source_fingerprint, ?2) = 0)";
    let sql = format!(
        "SELECT DISTINCT biosample_guid FROM haplogroup_call \
         WHERE dna_type = ?1 AND ( \
             ({unknown} AND ?4) \
             OR (NOT {unknown} \
                 AND substr(source_fingerprint, instr(source_fingerprint, ?2) + length(?2)) <> ?3))"
    );
    let rows: Vec<String> = sqlx::query_scalar(&sql)
        .bind(dna_type.as_str())
        .bind(tree_tag)
        .bind(tree_hash)
        .bind(include_unknown)
        .fetch_all(pool)
        .await?;
    Ok(rows.iter().filter_map(|g| g.parse().ok().map(SampleGuid)).collect())
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

/// Delete a single recorded call (one source) for a biosample + DNA type. Returns whether a row was
/// removed. Used to drop alignment-derived calls when their sequencing run/alignment is deleted.
pub async fn delete_one(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
    source_key: &str,
) -> Result<bool, StoreError> {
    let affected =
        sqlx::query("DELETE FROM haplogroup_call WHERE biosample_guid = ? AND dna_type = ? AND source_key = ?")
            .bind(biosample_guid.0.to_string())
            .bind(dna_type.as_str())
            .bind(source_key)
            .execute(pool)
            .await?
            .rows_affected();
    Ok(affected > 0)
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
struct ProvRow {
    source_label: String,
    haplogroup: String,
    lineage: String,
    score: f64,
    matched: i64,
    expected: i64,
    provenance: String,
}

impl ProvRow {
    fn into_domain(self) -> (CallProvenance, RunHaplogroupCall) {
        let provenance = CallProvenance::from_token(&self.provenance);
        let lineage = if self.lineage.is_empty() {
            Vec::new()
        } else {
            self.lineage.split('\t').map(str::to_string).collect()
        };
        (
            provenance,
            RunHaplogroupCall {
                source_label: self.source_label,
                haplogroup: self.haplogroup,
                lineage,
                score: self.score,
                matched: self.matched,
                expected: self.expected,
            },
        )
    }
}

/// All recorded calls for a biosample + DNA type, each tagged with its provenance tier — the input
/// to [`navigator_domain::reconciliation::reconcile_with_provenance`].
pub async fn list_for_with_provenance(
    pool: &SqlitePool,
    biosample_guid: SampleGuid,
    dna_type: DnaType,
) -> Result<Vec<(CallProvenance, RunHaplogroupCall)>, StoreError> {
    let rows: Vec<ProvRow> = sqlx::query_as(
        "SELECT source_label, haplogroup, lineage, score, matched, expected, provenance FROM haplogroup_call \
         WHERE biosample_guid = ? AND dna_type = ? ORDER BY id",
    )
    .bind(biosample_guid.0.to_string())
    .bind(dna_type.as_str())
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(ProvRow::into_domain).collect())
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
    provenance: String,
}

/// Every recorded call across all subjects, as `(guid, dna_type, provenance, call)` — for building a
/// donor-level haplogroup summary (the subjects list) in one query.
pub async fn list_all(
    pool: &SqlitePool,
) -> Result<Vec<(SampleGuid, DnaType, CallProvenance, RunHaplogroupCall)>, StoreError> {
    let rows: Vec<AllRow> = sqlx::query_as(
        "SELECT biosample_guid, dna_type, source_label, haplogroup, lineage, score, matched, expected, provenance \
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
            CallProvenance::from_token(&r.provenance),
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
