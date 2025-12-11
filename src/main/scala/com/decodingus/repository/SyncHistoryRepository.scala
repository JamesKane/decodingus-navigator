package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Sync direction.
 */
enum SyncDirection:
  case Push
  case Pull

object SyncDirection:
  def fromString(s: String): SyncDirection = s match
    case "PUSH" => Push
    case "PULL" => Pull
    case other => throw new IllegalArgumentException(s"Unknown direction: $other")

  def toDbString(d: SyncDirection): String = d match
    case Push => "PUSH"
    case Pull => "PULL"

/**
 * Sync result status.
 */
enum SyncResultStatus:
  case Success
  case Failed
  case Conflict
  case Skipped

object SyncResultStatus:
  def fromString(s: String): SyncResultStatus = s match
    case "SUCCESS" => Success
    case "FAILED" => Failed
    case "CONFLICT" => Conflict
    case "SKIPPED" => Skipped
    case other => throw new IllegalArgumentException(s"Unknown sync result status: $other")

  def toDbString(s: SyncResultStatus): String = s match
    case Success => "SUCCESS"
    case Failed => "FAILED"
    case Conflict => "CONFLICT"
    case Skipped => "SKIPPED"

/**
 * Sync history entry entity - audit trail of sync operations.
 */
case class SyncHistoryEntity(
  id: UUID,
  entityType: SyncEntityType,
  entityId: UUID,
  atUri: Option[String],
  operation: SyncOperation,
  direction: SyncDirection,
  status: SyncResultStatus,
  errorMessage: Option[String],
  startedAt: LocalDateTime,
  completedAt: LocalDateTime,
  durationMs: Option[Long],
  localVersionBefore: Option[Int],
  localVersionAfter: Option[Int],
  remoteVersionBefore: Option[Int],
  remoteVersionAfter: Option[Int],
  localCidBefore: Option[String],
  localCidAfter: Option[String],
  remoteCid: Option[String],
  createdAt: LocalDateTime
) extends Entity[UUID]

object SyncHistoryEntity:
  /**
   * Create a new sync history entry.
   */
  def create(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    direction: SyncDirection,
    status: SyncResultStatus,
    startedAt: LocalDateTime,
    completedAt: LocalDateTime,
    atUri: Option[String] = None,
    errorMessage: Option[String] = None,
    localVersionBefore: Option[Int] = None,
    localVersionAfter: Option[Int] = None,
    remoteVersionBefore: Option[Int] = None,
    remoteVersionAfter: Option[Int] = None,
    localCidBefore: Option[String] = None,
    localCidAfter: Option[String] = None,
    remoteCid: Option[String] = None
  ): SyncHistoryEntity =
    val now = LocalDateTime.now()
    val durationMs = java.time.Duration.between(startedAt, completedAt).toMillis
    SyncHistoryEntity(
      id = UUID.randomUUID(),
      entityType = entityType,
      entityId = entityId,
      atUri = atUri,
      operation = operation,
      direction = direction,
      status = status,
      errorMessage = errorMessage,
      startedAt = startedAt,
      completedAt = completedAt,
      durationMs = Some(durationMs),
      localVersionBefore = localVersionBefore,
      localVersionAfter = localVersionAfter,
      remoteVersionBefore = remoteVersionBefore,
      remoteVersionAfter = remoteVersionAfter,
      localCidBefore = localCidBefore,
      localCidAfter = localCidAfter,
      remoteCid = remoteCid,
      createdAt = now
    )

/**
 * Repository for sync history (audit trail).
 */
class SyncHistoryRepository:

  // ============================================
  // Core Operations
  // ============================================

  def findById(id: UUID)(using conn: Connection): Option[SyncHistoryEntity] =
    queryOne(
      "SELECT * FROM sync_history WHERE id = ?",
      Seq(id)
    )(mapRow)

  def findAll(limit: Int = 100)(using conn: Connection): List[SyncHistoryEntity] =
    queryList(
      "SELECT * FROM sync_history ORDER BY completed_at DESC LIMIT ?",
      Seq(limit)
    )(mapRow)

  def insert(entity: SyncHistoryEntity)(using conn: Connection): SyncHistoryEntity =
    executeUpdate(
      """INSERT INTO sync_history (
        |  id, entity_type, entity_id, at_uri, operation, direction, status,
        |  error_message, started_at, completed_at, duration_ms,
        |  local_version_before, local_version_after, remote_version_before, remote_version_after,
        |  local_cid_before, local_cid_after, remote_cid, created_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        SyncEntityType.toDbString(entity.entityType),
        entity.entityId,
        entity.atUri,
        SyncOperation.toDbString(entity.operation),
        SyncDirection.toDbString(entity.direction),
        SyncResultStatus.toDbString(entity.status),
        entity.errorMessage,
        entity.startedAt,
        entity.completedAt,
        entity.durationMs,
        entity.localVersionBefore,
        entity.localVersionAfter,
        entity.remoteVersionBefore,
        entity.remoteVersionAfter,
        entity.localCidBefore,
        entity.localCidAfter,
        entity.remoteCid,
        entity.createdAt
      )
    )
    entity

  def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM sync_history WHERE id = ?", Seq(id)) > 0

  // ============================================
  // Query Operations
  // ============================================

  /**
   * Find history for a specific entity.
   */
  def findByEntity(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): List[SyncHistoryEntity] =
    queryList(
      "SELECT * FROM sync_history WHERE entity_type = ? AND entity_id = ? ORDER BY completed_at DESC",
      Seq(SyncEntityType.toDbString(entityType), entityId)
    )(mapRow)

  /**
   * Find history by status.
   */
  def findByStatus(status: SyncResultStatus, limit: Int = 100)(using conn: Connection): List[SyncHistoryEntity] =
    queryList(
      "SELECT * FROM sync_history WHERE status = ? ORDER BY completed_at DESC LIMIT ?",
      Seq(SyncResultStatus.toDbString(status), limit)
    )(mapRow)

  /**
   * Find history by direction (push/pull).
   */
  def findByDirection(direction: SyncDirection, limit: Int = 100)(using conn: Connection): List[SyncHistoryEntity] =
    queryList(
      "SELECT * FROM sync_history WHERE direction = ? ORDER BY completed_at DESC LIMIT ?",
      Seq(SyncDirection.toDbString(direction), limit)
    )(mapRow)

  /**
   * Find recent failures.
   */
  def findRecentFailures(limit: Int = 50)(using conn: Connection): List[SyncHistoryEntity] =
    queryList(
      "SELECT * FROM sync_history WHERE status = 'FAILED' ORDER BY completed_at DESC LIMIT ?",
      Seq(limit)
    )(mapRow)

  /**
   * Find history within a time range.
   */
  def findInTimeRange(start: LocalDateTime, end: LocalDateTime)(using conn: Connection): List[SyncHistoryEntity] =
    queryList(
      "SELECT * FROM sync_history WHERE completed_at >= ? AND completed_at <= ? ORDER BY completed_at DESC",
      Seq(start, end)
    )(mapRow)

  /**
   * Get last sync for an entity.
   */
  def getLastSyncForEntity(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): Option[SyncHistoryEntity] =
    queryOne(
      "SELECT * FROM sync_history WHERE entity_type = ? AND entity_id = ? ORDER BY completed_at DESC LIMIT 1",
      Seq(SyncEntityType.toDbString(entityType), entityId)
    )(mapRow)

  /**
   * Get last successful sync for an entity.
   */
  def getLastSuccessfulSync(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): Option[SyncHistoryEntity] =
    queryOne(
      "SELECT * FROM sync_history WHERE entity_type = ? AND entity_id = ? AND status = 'SUCCESS' ORDER BY completed_at DESC LIMIT 1",
      Seq(SyncEntityType.toDbString(entityType), entityId)
    )(mapRow)

  // ============================================
  // Convenience Recording Methods
  // ============================================

  /**
   * Record a successful sync.
   */
  def recordSuccess(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    direction: SyncDirection,
    startedAt: LocalDateTime,
    atUri: Option[String] = None,
    localVersionBefore: Option[Int] = None,
    localVersionAfter: Option[Int] = None,
    remoteCid: Option[String] = None
  )(using conn: Connection): SyncHistoryEntity =
    val entity = SyncHistoryEntity.create(
      entityType = entityType,
      entityId = entityId,
      operation = operation,
      direction = direction,
      status = SyncResultStatus.Success,
      startedAt = startedAt,
      completedAt = LocalDateTime.now(),
      atUri = atUri,
      localVersionBefore = localVersionBefore,
      localVersionAfter = localVersionAfter,
      remoteCid = remoteCid
    )
    insert(entity)

  /**
   * Record a failed sync.
   */
  def recordFailure(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    direction: SyncDirection,
    startedAt: LocalDateTime,
    errorMessage: String,
    atUri: Option[String] = None,
    localVersionBefore: Option[Int] = None
  )(using conn: Connection): SyncHistoryEntity =
    val entity = SyncHistoryEntity.create(
      entityType = entityType,
      entityId = entityId,
      operation = operation,
      direction = direction,
      status = SyncResultStatus.Failed,
      startedAt = startedAt,
      completedAt = LocalDateTime.now(),
      atUri = atUri,
      errorMessage = Some(errorMessage),
      localVersionBefore = localVersionBefore
    )
    insert(entity)

  /**
   * Record a conflict.
   */
  def recordConflict(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    direction: SyncDirection,
    startedAt: LocalDateTime,
    localVersionBefore: Int,
    remoteVersionBefore: Int,
    atUri: Option[String] = None
  )(using conn: Connection): SyncHistoryEntity =
    val entity = SyncHistoryEntity.create(
      entityType = entityType,
      entityId = entityId,
      operation = operation,
      direction = direction,
      status = SyncResultStatus.Conflict,
      startedAt = startedAt,
      completedAt = LocalDateTime.now(),
      atUri = atUri,
      localVersionBefore = Some(localVersionBefore),
      remoteVersionBefore = Some(remoteVersionBefore)
    )
    insert(entity)

  // ============================================
  // Statistics
  // ============================================

  /**
   * Count entries by status.
   */
  def countByStatus()(using conn: Connection): Map[SyncResultStatus, Long] =
    queryList(
      "SELECT status, COUNT(*) as cnt FROM sync_history GROUP BY status"
    ) { rs =>
      val status = SyncResultStatus.fromString(rs.getString("status"))
      val count = rs.getLong("cnt")
      (status, count)
    }.toMap

  /**
   * Get sync statistics for a time period.
   */
  def getStatsForPeriod(start: LocalDateTime, end: LocalDateTime)(using conn: Connection): SyncStats =
    val results = findInTimeRange(start, end)
    SyncStats(
      total = results.size,
      successful = results.count(_.status == SyncResultStatus.Success),
      failed = results.count(_.status == SyncResultStatus.Failed),
      conflicts = results.count(_.status == SyncResultStatus.Conflict),
      skipped = results.count(_.status == SyncResultStatus.Skipped),
      pushes = results.count(_.direction == SyncDirection.Push),
      pulls = results.count(_.direction == SyncDirection.Pull),
      avgDurationMs = if results.isEmpty then 0L else results.flatMap(_.durationMs).sum / results.size
    )

  /**
   * Cleanup old history entries.
   */
  def cleanupOlderThan(olderThanDays: Int)(using conn: Connection): Int =
    val cutoff = LocalDateTime.now().minusDays(olderThanDays)
    executeUpdate(
      "DELETE FROM sync_history WHERE created_at < ?",
      Seq(cutoff)
    )

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): SyncHistoryEntity =
    SyncHistoryEntity(
      id = getUUID(rs, "id"),
      entityType = SyncEntityType.fromString(rs.getString("entity_type")),
      entityId = getUUID(rs, "entity_id"),
      atUri = getOptString(rs, "at_uri"),
      operation = SyncOperation.fromString(rs.getString("operation")),
      direction = SyncDirection.fromString(rs.getString("direction")),
      status = SyncResultStatus.fromString(rs.getString("status")),
      errorMessage = getOptString(rs, "error_message"),
      startedAt = getDateTime(rs, "started_at"),
      completedAt = getDateTime(rs, "completed_at"),
      durationMs = getOptLong(rs, "duration_ms"),
      localVersionBefore = getOptInt(rs, "local_version_before"),
      localVersionAfter = getOptInt(rs, "local_version_after"),
      remoteVersionBefore = getOptInt(rs, "remote_version_before"),
      remoteVersionAfter = getOptInt(rs, "remote_version_after"),
      localCidBefore = getOptString(rs, "local_cid_before"),
      localCidAfter = getOptString(rs, "local_cid_after"),
      remoteCid = getOptString(rs, "remote_cid"),
      createdAt = getDateTime(rs, "created_at")
    )

/**
 * Sync statistics summary.
 */
case class SyncStats(
  total: Int,
  successful: Int,
  failed: Int,
  conflicts: Int,
  skipped: Int,
  pushes: Int,
  pulls: Int,
  avgDurationMs: Long
)
