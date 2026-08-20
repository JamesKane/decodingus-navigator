//! The Y and mtDNA haplogroup call of each source, which are the inputs to the donor-level
//! reconcile. Each (biosample, dna_type, source) triple has one row, and an upsert replaces a
//! second run from the same source.

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

/// Insert or replace the call from `source_key`, for this biosample and DNA type.
///
/// `provenance` records which caller produced it: external, navigator-walk, or manual. That is the
/// tier that the reconcile uses to decide precedence. `fingerprint` stamps the inputs, which are
/// the content hashes of the file and the tree, so a later run can skip the score step.
///
/// An external call and an internal call use *different* `source_key`s. So this upsert can never
/// let one of them overwrite the other.
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

/// The stored input fingerprint of one source's call, if the store holds one. The caller reads it
/// to decide whether it must score again. The inputs are the same when this equals the current
/// fingerprint.
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

/// Subjects whose recorded calls sit on **a different haplotree** from the one now active. This is
/// the selector for a fresh look at the workspace, after a new tree arrives.
///
/// A fingerprint reads `<source tag>:<hash>|<tree tag>:<hash>`. An alignment-placed Y call gives
/// `f:…|yt:…`, and a GVCF-placed mt call gives `gv:…|mt:…`. So `tree_tag` is `"yt:"` or `"mt:"`, and
/// `tree_hash` is the first 16 hex characters of the current tree's SHA-256.
///
/// `include_unknown` decides what to do with a call that carries **no tree tag**. Such a call comes
/// from before the fingerprint existed, so nobody can say which tree it used.
///
/// Those calls are a different job from a tree change, and they cost far more. In this workspace
/// 3,780 Y calls are *provably* on a tree that a newer one replaced, and 15,648 carry no
/// fingerprint at all. Most of that second group owns a BAM that a second placement would walk
/// again. So they stay out by default, and the routine "a new tree arrived" sweep stays quick. A
/// caller must ask for the backfill.
pub async fn biosamples_placed_against_another_tree(
    pool: &SqlitePool,
    dna_type: DnaType,
    tree_tag: &str,
    tree_hash: &str,
    include_unknown: bool,
) -> Result<Vec<SampleGuid>, StoreError> {
    // `instr(...) = 0` catches a fingerprint with no tree tag at all. Without it, `substr` would
    // cut from offset 1, and the result could in theory equal the hash. `include_unknown` decides
    // whether that case *selects*. It must never reach the comparison against the hash.
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

/// The recorded call of one source. The caller returns it as a cached result, and does not score
/// again.
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

/// Delete one recorded call, from one source, for a biosample and DNA type. It returns whether it
/// removed a row. A delete of a sequencing run or an alignment uses it to drop the calls that came
/// from that alignment.
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

/// Every recorded call for a biosample and DNA type, each one with its provenance tier. This is
/// the input to [`navigator_domain::reconciliation::reconcile_with_provenance`].
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

/// Every recorded call of every subject, as `(guid, dna_type, provenance, call)`. One query then
/// builds the donor-level haplogroup summary for the subjects list.
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
