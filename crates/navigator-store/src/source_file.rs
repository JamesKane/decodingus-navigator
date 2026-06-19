//! Content-hash file identity (gap §5-p2). Imported files are keyed by their SHA-256, not their path —
//! so moving/renaming a BAM updates the path in place instead of orphaning its analyses, and a
//! re-import of the same content is recognised as a duplicate. Ports the legacy `SourceFileRepository`.

use sqlx::SqlitePool;

use crate::StoreError;

/// A tracked source file, identified by its content hash.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct SourceFile {
    pub id: i64,
    pub content_sha256: String,
    pub file_path: Option<String>,
    pub file_size: Option<i64>,
    pub file_format: Option<String>,
    pub alignment_id: Option<i64>,
    pub is_accessible: bool,
    pub last_verified_at: Option<String>,
    pub created_at: String,
    pub updated_at: String,
}

/// Register a file by content hash, or update the recorded path/size/format if its content is already
/// known (the file moved). Returns the row. Marks it accessible + verified now.
pub async fn upsert_by_checksum(
    pool: &SqlitePool,
    content_sha256: &str,
    file_path: Option<&str>,
    file_size: Option<i64>,
    file_format: Option<&str>,
    now: &str,
) -> Result<SourceFile, StoreError> {
    sqlx::query(
        "INSERT INTO source_file (content_sha256, file_path, file_size, file_format, is_accessible, last_verified_at, created_at, updated_at) \
         VALUES (?, ?, ?, ?, 1, ?, ?, ?) \
         ON CONFLICT(content_sha256) DO UPDATE SET \
         file_path = excluded.file_path, file_size = excluded.file_size, file_format = excluded.file_format, \
         is_accessible = 1, last_verified_at = excluded.last_verified_at, updated_at = excluded.updated_at",
    )
    .bind(content_sha256)
    .bind(file_path)
    .bind(file_size)
    .bind(file_format)
    .bind(now)
    .bind(now)
    .bind(now)
    .execute(pool)
    .await?;
    Ok(find_by_checksum(pool, content_sha256).await?.expect("just upserted"))
}

/// Look up a file by its content hash (the dedup probe on import).
pub async fn find_by_checksum(pool: &SqlitePool, content_sha256: &str) -> Result<Option<SourceFile>, StoreError> {
    let row = sqlx::query_as("SELECT * FROM source_file WHERE content_sha256 = ?")
        .bind(content_sha256)
        .fetch_optional(pool)
        .await?;
    Ok(row)
}

/// Link a tracked file to the alignment produced from it.
pub async fn link_to_alignment(
    pool: &SqlitePool,
    content_sha256: &str,
    alignment_id: i64,
    now: &str,
) -> Result<(), StoreError> {
    sqlx::query("UPDATE source_file SET alignment_id = ?, updated_at = ? WHERE content_sha256 = ?")
        .bind(alignment_id)
        .bind(now)
        .bind(content_sha256)
        .execute(pool)
        .await?;
    Ok(())
}

/// Update a file's accessibility (path still resolves on disk?) after a re-verify pass.
pub async fn set_accessible(pool: &SqlitePool, id: i64, accessible: bool, now: &str) -> Result<(), StoreError> {
    sqlx::query("UPDATE source_file SET is_accessible = ?, last_verified_at = ?, updated_at = ? WHERE id = ?")
        .bind(accessible)
        .bind(now)
        .bind(now)
        .bind(id)
        .execute(pool)
        .await?;
    Ok(())
}

/// All tracked source files.
pub async fn list(pool: &SqlitePool) -> Result<Vec<SourceFile>, StoreError> {
    let rows = sqlx::query_as("SELECT * FROM source_file ORDER BY created_at DESC")
        .fetch_all(pool)
        .await?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn upsert_moves_path_without_duplicating() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let t = "2026-06-17T00:00:00Z";
        let a = upsert_by_checksum(store.pool(), "deadbeef", Some("/data/x.bam"), Some(100), Some("BAM"), t)
            .await
            .unwrap();
        assert_eq!(a.file_path.as_deref(), Some("/data/x.bam"));
        assert!(a.is_accessible);
        // Same content, new path → same row, path updated (no duplicate).
        let b = upsert_by_checksum(
            store.pool(),
            "deadbeef",
            Some("/archive/x.bam"),
            Some(100),
            Some("BAM"),
            t,
        )
        .await
        .unwrap();
        assert_eq!(a.id, b.id);
        assert_eq!(b.file_path.as_deref(), Some("/archive/x.bam"));
        assert_eq!(list(store.pool()).await.unwrap().len(), 1);
        // Dedup probe + accessibility flip.
        assert!(find_by_checksum(store.pool(), "deadbeef").await.unwrap().is_some());
        set_accessible(store.pool(), a.id, false, t).await.unwrap();
        assert!(
            !find_by_checksum(store.pool(), "deadbeef")
                .await
                .unwrap()
                .unwrap()
                .is_accessible
        );
    }
}
