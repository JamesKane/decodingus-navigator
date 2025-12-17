package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Abstract base class for repositories with sync status tracking.
 *
 * Provides default implementations for common sync operations, reducing
 * code duplication across repositories. Subclasses only need to provide:
 *   - tableName: The database table name
 *   - mapRow: Entity mapping from ResultSet
 *
 * @tparam E Entity type extending Entity[UUID]
 */
abstract class SyncableRepositoryBase[E <: Entity[UUID]]
  extends SyncableRepository[E, UUID]:

  /**
   * The database table name for this entity type.
   * Used in SQL queries for sync operations.
   */
  protected def tableName: String

  /**
   * Map a ResultSet row to an entity instance.
   * Each subclass implements entity-specific mapping.
   */
  protected def mapRow(rs: ResultSet): E

  // ============================================
  // Default Sync Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[E] =
    queryList(
      s"SELECT * FROM $tableName WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      s"UPDATE $tableName SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      s"""UPDATE $tableName SET
         |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
         |WHERE id = ?
       """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  override def count()(using conn: Connection): Long =
    queryOne(s"SELECT COUNT(*) FROM $tableName") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  /**
   * Find entities pending sync (Local or Modified status).
   * Commonly used for sync queue processing.
   */
  def findPendingSync()(using conn: Connection): List[E] =
    queryList(
      s"""SELECT * FROM $tableName
         |WHERE sync_status IN ('Local', 'Modified')
         |ORDER BY updated_at ASC
       """.stripMargin
    )(mapRow)

  // ============================================
  // Entity Metadata Helpers
  // ============================================

  /**
   * Extract EntityMeta from a ResultSet.
   *
   * Shared helper for mapRow implementations. Expects standard column names:
   * sync_status, at_uri, at_cid, version, created_at, updated_at
   */
  protected def readEntityMeta(rs: ResultSet): EntityMeta = EntityMeta(
    syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
    atUri = getOptString(rs, "at_uri"),
    atCid = getOptString(rs, "at_cid"),
    version = rs.getInt("version"),
    createdAt = getDateTime(rs, "created_at"),
    updatedAt = getDateTime(rs, "updated_at")
  )
