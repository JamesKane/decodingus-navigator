//! Subject SNP-variant queries: variant sets and their calls. Sets attach to a biosample
//! (the `SampleGuid` is stored as its hyphenated TEXT form, like elsewhere).

use du_domain::ids::SampleGuid;
use navigator_domain::variants::{NewVariantSet, SourceType, VariantCall, VariantSet};
use sqlx::SqlitePool;
use uuid::Uuid;

use crate::StoreError;

#[derive(sqlx::FromRow)]
struct SetRow {
    id: i64,
    biosample_guid: String,
    source_label: String,
    source_type: String,
}

#[derive(sqlx::FromRow)]
struct CallRow {
    contig: String,
    position: i64,
    reference: String,
    alternate: String,
    rs_id: Option<String>,
    genotype: Option<String>,
}

impl CallRow {
    fn into_domain(self) -> VariantCall {
        VariantCall {
            contig: self.contig,
            position: self.position,
            reference: self.reference,
            alternate: self.alternate,
            rs_id: self.rs_id,
            genotype: self.genotype,
        }
    }
}

/// Create a variant set and bulk-insert its calls in one transaction.
pub async fn create(pool: &SqlitePool, new: &NewVariantSet) -> Result<VariantSet, StoreError> {
    let mut tx = pool.begin().await?;
    let id: i64 = sqlx::query_scalar(
        "INSERT INTO variant_set (biosample_guid, source_label, source_type) VALUES (?, ?, ?) RETURNING id",
    )
    .bind(new.biosample_guid.0.to_string())
    .bind(&new.source_label)
    .bind(new.source_type.as_str())
    .fetch_one(&mut *tx)
    .await?;
    for c in &new.calls {
        sqlx::query(
            "INSERT INTO variant_call (variant_set_id, contig, position, reference, alternate, rs_id, genotype) \
             VALUES (?, ?, ?, ?, ?, ?, ?)",
        )
        .bind(id)
        .bind(&c.contig)
        .bind(c.position)
        .bind(&c.reference)
        .bind(&c.alternate)
        .bind(&c.rs_id)
        .bind(&c.genotype)
        .execute(&mut *tx)
        .await?;
    }
    tx.commit().await?;
    Ok(VariantSet {
        id,
        biosample_guid: new.biosample_guid,
        source_label: new.source_label.clone(),
        source_type: new.source_type,
        calls: new.calls.clone(),
    })
}

async fn calls_for(pool: &SqlitePool, set_id: i64) -> Result<Vec<VariantCall>, StoreError> {
    let rows: Vec<CallRow> = sqlx::query_as(
        "SELECT contig, position, reference, alternate, rs_id, genotype FROM variant_call \
         WHERE variant_set_id = ? ORDER BY id",
    )
    .bind(set_id)
    .fetch_all(pool)
    .await?;
    Ok(rows.into_iter().map(CallRow::into_domain).collect())
}

/// Delete a variant set and its calls (children-first; FKs are enforced). Returns whether the
/// set row was removed.
pub async fn delete(pool: &SqlitePool, id: i64) -> Result<bool, StoreError> {
    let mut tx = pool.begin().await?;
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
        "SELECT id, biosample_guid, source_label, source_type FROM variant_set WHERE biosample_guid = ? ORDER BY id",
    )
    .bind(guid.0.to_string())
    .fetch_all(pool)
    .await?;

    let mut sets = Vec::with_capacity(rows.len());
    for r in rows {
        let uuid = Uuid::parse_str(&r.biosample_guid)
            .map_err(|e| StoreError::Decode(format!("variant_set guid {:?}: {e}", r.biosample_guid)))?;
        let calls = calls_for(pool, r.id).await?;
        sets.push(VariantSet {
            id: r.id,
            biosample_guid: SampleGuid(uuid),
            source_label: r.source_label,
            source_type: SourceType::from_code(&r.source_type),
            calls,
        });
    }
    Ok(sets)
}
