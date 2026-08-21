//! Subject SNP-variant queries: the variant sets and their calls. A set attaches to a biosample.
//! The `SampleGuid` goes in as its hyphenated TEXT form, as it does everywhere else.

use du_domain::ids::SampleGuid;
use navigator_domain::variants::{self, CallEvidence, NewVariantSet, SourceType, VariantCall, VariantSet};
use sqlx::SqlitePool;

use crate::error::parse_sample_guid;
use crate::StoreError;

#[derive(sqlx::FromRow)]
struct SetRow {
    id: i64,
    biosample_guid: String,
    source_label: String,
    source_type: String,
    reference_build: Option<String>,
    call_schema: i64,
    source_path: Option<String>,
}

#[derive(sqlx::FromRow)]
struct CallRow {
    contig: String,
    position: i64,
    reference: String,
    alternate: String,
    rs_id: Option<String>,
    genotype: Option<String>,
    qual: Option<f64>,
    filter: Option<String>,
    dp: Option<i64>,
    gq: Option<i64>,
    ad_ref: Option<i64>,
    ad_alt: Option<i64>,
}

/// Columns read back for a call, in `CallRow` order.
const CALL_COLS: &str = "contig, position, reference, alternate, rs_id, genotype, qual, filter, dp, gq, ad_ref, ad_alt";

impl CallRow {
    fn into_domain(self) -> VariantCall {
        // The column is an INTEGER, because SQLite has no unsigned type. A negative value would be
        // corrupt data, so the read gives absent. It must not wrap into a huge count.
        let count = |v: Option<i64>| v.and_then(|n| u32::try_from(n).ok());
        VariantCall {
            contig: self.contig,
            position: self.position,
            reference: self.reference,
            alternate: self.alternate,
            rs_id: self.rs_id,
            genotype: self.genotype,
            evidence: CallEvidence {
                qual: self.qual,
                filter: self.filter,
                dp: count(self.dp),
                gq: count(self.gq),
                ad_ref: count(self.ad_ref),
                ad_alt: count(self.ad_alt),
            },
        }
    }
}

/// Create a variant set and bulk-insert its calls in one transaction.
pub async fn create(pool: &SqlitePool, new: &NewVariantSet) -> Result<VariantSet, StoreError> {
    let mut tx = pool.begin().await?;
    // The schema tag comes from what the import captured, and not from which importer ran.
    //
    // A sites-only VCF, and a CSV marker table, truly hold no evidence. A consumer that asks "can I
    // gate on quality here?" needs the true answer. A tag keyed to the importer version would claim
    // evidence that those sets can never supply.
    let schema = if new.calls.iter().any(|c| !c.evidence.is_empty()) {
        variants::CALL_SCHEMA_EVIDENCE
    } else {
        variants::CALL_SCHEMA_BASIC
    };
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO variant_set (biosample_guid, source_label, source_type, reference_build, call_schema, source_path) \
         VALUES (?, ?, ?, ?, ?, ?) RETURNING id",
    )
    .bind(new.biosample_guid.0.to_string())
    .bind(&new.source_label)
    .bind(new.source_type.as_str())
    .bind(&new.reference_build)
    .bind(schema)
    .bind(&new.source_path)
    .fetch_one(&mut *tx)
    .await?;
    for c in &new.calls {
        let e = &c.evidence;
        sqlx::query(
            "INSERT INTO variant_call \
             (variant_set_id, contig, position, reference, alternate, rs_id, genotype, \
              qual, filter, dp, gq, ad_ref, ad_alt) \
             VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        )
        .bind(id)
        .bind(&c.contig)
        .bind(c.position)
        .bind(&c.reference)
        .bind(&c.alternate)
        .bind(&c.rs_id)
        .bind(&c.genotype)
        .bind(e.qual)
        .bind(&e.filter)
        .bind(e.dp.map(i64::from))
        .bind(e.gq.map(i64::from))
        .bind(e.ad_ref.map(i64::from))
        .bind(e.ad_alt.map(i64::from))
        .execute(&mut *tx)
        .await?;
    }
    tx.commit().await?;
    Ok(VariantSet {
        id,
        biosample_guid: new.biosample_guid,
        source_label: new.source_label.clone(),
        source_type: new.source_type,
        reference_build: new.reference_build.clone(),
        calls: new.calls.clone(),
        call_schema: schema,
        source_path: new.source_path.clone(),
    })
}

/// One variant set (with its calls) by id, or `None` if it does not exist.
pub async fn get(pool: &SqlitePool, id: i64) -> Result<Option<VariantSet>, StoreError> {
    let Some(r) = sqlx::query_as::<_, SetRow>(
        "SELECT id, biosample_guid, source_label, source_type, reference_build, call_schema, source_path FROM variant_set WHERE id = ?",
    )
    .bind(id)
    .fetch_optional(pool)
    .await?
    else {
        return Ok(None);
    };
    let biosample_guid = parse_sample_guid(&r.biosample_guid, "variant_set")?;
    let calls = calls_for(pool, r.id).await?;
    Ok(Some(VariantSet {
        id: r.id,
        biosample_guid,
        source_label: r.source_label,
        source_type: SourceType::from_code(&r.source_type),
        reference_build: r.reference_build,
        calls,
        call_schema: r.call_schema,
        source_path: r.source_path,
    }))
}

async fn calls_for(pool: &SqlitePool, set_id: i64) -> Result<Vec<VariantCall>, StoreError> {
    let rows: Vec<CallRow> = sqlx::query_as(&format!(
        "SELECT {CALL_COLS} FROM variant_call WHERE variant_set_id = ? ORDER BY id"
    ))
    .bind(set_id)
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(CallRow::into_domain).collect())
}

/// Delete a variant set and its calls. The children go first, because the database enforces the
/// FKs. It returns whether it removed the set row.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let mut tx = pool.begin().await?;
    sqlx::query("DELETE FROM variant_set_genotype WHERE variant_set_id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?;
    sqlx::query("DELETE FROM variant_call WHERE variant_set_id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?;
    let affected = sqlx::query("DELETE FROM variant_set WHERE id = ?")
        .bind(id)
        .execute(&mut *tx)
        .await?
        .rows_affected();
    tx.commit().await?;
    Ok(affected > 0)
}

/// All variant sets for a biosample, with their calls.
pub async fn list_for_biosample(pool: &SqlitePool, guid: SampleGuid) -> Result<Vec<VariantSet>, StoreError> {
    let rows: Vec<SetRow> = sqlx::query_as(
        "SELECT id, biosample_guid, source_label, source_type, reference_build, call_schema, source_path FROM variant_set \
         WHERE biosample_guid = ? ORDER BY id",
    )
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;

    let mut sets = Vec::with_capacity(rows.len());
    for r in rows {
        let biosample_guid = parse_sample_guid(&r.biosample_guid, "variant_set")?;
        let calls = calls_for(pool, r.id).await?;
        sets.push(VariantSet {
            id: r.id,
            biosample_guid,
            source_label: r.source_label,
            source_type: SourceType::from_code(&r.source_type),
            reference_build: r.reference_build,
            calls,
            call_schema: r.call_schema,
            source_path: r.source_path,
        });
    }
    Ok(sets)
}
