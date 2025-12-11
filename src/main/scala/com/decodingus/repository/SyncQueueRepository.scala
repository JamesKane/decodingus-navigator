package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import io.circe.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Entity types that can be synced to PDS.
 */
enum SyncEntityType:
  case Biosample
  case Project
  case SequenceRun
  case Alignment
  case StrProfile
  case ChipProfile
  case YSnpPanel
  case HaplogroupReconciliation

object SyncEntityType:
  def fromString(s: String): SyncEntityType = s match
    case "BIOSAMPLE" => Biosample
    case "PROJECT" => Project
    case "SEQUENCE_RUN" => SequenceRun
    case "ALIGNMENT" => Alignment
    case "STR_PROFILE" => StrProfile
    case "CHIP_PROFILE" => ChipProfile
    case "Y_SNP_PANEL" => YSnpPanel
    case "HAPLOGROUP_RECONCILIATION" => HaplogroupReconciliation
    case other => throw new IllegalArgumentException(s"Unknown entity type: $other")

  def toDbString(t: SyncEntityType): String = t match
    case Biosample => "BIOSAMPLE"
    case Project => "PROJECT"
    case SequenceRun => "SEQUENCE_RUN"
    case Alignment => "ALIGNMENT"
    case StrProfile => "STR_PROFILE"
    case ChipProfile => "CHIP_PROFILE"
    case YSnpPanel => "Y_SNP_PANEL"
    case HaplogroupReconciliation => "HAPLOGROUP_RECONCILIATION"

/**
 * Sync operations.
 */
enum SyncOperation:
  case Create
  case Update
  case Delete

object SyncOperation:
  def fromString(s: String): SyncOperation = s match
    case "CREATE" => Create
    case "UPDATE" => Update
    case "DELETE" => Delete
    case other => throw new IllegalArgumentException(s"Unknown operation: $other")

  def toDbString(o: SyncOperation): String = o match
    case Create => "CREATE"
    case Update => "UPDATE"
    case Delete => "DELETE"

/**
 * Queue entry status.
 */
enum QueueStatus:
  case Pending
  case InProgress
  case Completed
  case Failed
  case Cancelled

object QueueStatus:
  def fromString(s: String): QueueStatus = s match
    case "PENDING" => Pending
    case "IN_PROGRESS" => InProgress
    case "COMPLETED" => Completed
    case "FAILED" => Failed
    case "CANCELLED" => Cancelled
    case other => throw new IllegalArgumentException(s"Unknown queue status: $other")

  def toDbString(s: QueueStatus): String = s match
    case Pending => "PENDING"
    case InProgress => "IN_PROGRESS"
    case Completed => "COMPLETED"
    case Failed => "FAILED"
    case Cancelled => "CANCELLED"

/**
 * Sync queue entry entity.
 */
case class SyncQueueEntity(
  id: UUID,
  entityType: SyncEntityType,
  entityId: UUID,
  operation: SyncOperation,
  status: QueueStatus,
  priority: Int,
  queuedAt: LocalDateTime,
  startedAt: Option[LocalDateTime],
  completedAt: Option[LocalDateTime],
  attemptCount: Int,
  nextRetryAt: Option[LocalDateTime],
  lastError: Option[String],
  payloadSnapshot: Option[Json],
  createdAt: LocalDateTime,
  updatedAt: LocalDateTime
) extends Entity[UUID]

object SyncQueueEntity:
  /**
   * Create a new sync queue entry.
   */
  def create(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    priority: Int = 5,
    payloadSnapshot: Option[Json] = None
  ): SyncQueueEntity =
    val now = LocalDateTime.now()
    SyncQueueEntity(
      id = UUID.randomUUID(),
      entityType = entityType,
      entityId = entityId,
      operation = operation,
      status = QueueStatus.Pending,
      priority = priority,
      queuedAt = now,
      startedAt = None,
      completedAt = None,
      attemptCount = 0,
      nextRetryAt = None,
      lastError = None,
      payloadSnapshot = payloadSnapshot,
      createdAt = now,
      updatedAt = now
    )

/**
 * Repository for sync queue management.
 */
class SyncQueueRepository:

  // ============================================
  // Core Operations
  // ============================================

  def findById(id: UUID)(using conn: Connection): Option[SyncQueueEntity] =
    queryOne(
      "SELECT * FROM sync_queue WHERE id = ?",
      Seq(id)
    )(mapRow)

  def findAll()(using conn: Connection): List[SyncQueueEntity] =
    queryList("SELECT * FROM sync_queue ORDER BY priority ASC, queued_at ASC")(mapRow)

  def insert(entity: SyncQueueEntity)(using conn: Connection): SyncQueueEntity =
    val payloadJson = entity.payloadSnapshot.map(j => JsonValue(j.noSpaces))

    executeUpdate(
      """INSERT INTO sync_queue (
        |  id, entity_type, entity_id, operation, status, priority, queued_at,
        |  started_at, completed_at, attempt_count, next_retry_at, last_error,
        |  payload_snapshot, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        SyncEntityType.toDbString(entity.entityType),
        entity.entityId,
        SyncOperation.toDbString(entity.operation),
        QueueStatus.toDbString(entity.status),
        entity.priority,
        entity.queuedAt,
        entity.startedAt,
        entity.completedAt,
        entity.attemptCount,
        entity.nextRetryAt,
        entity.lastError,
        payloadJson,
        entity.createdAt,
        entity.updatedAt
      )
    )
    entity

  def update(entity: SyncQueueEntity)(using conn: Connection): SyncQueueEntity =
    val payloadJson = entity.payloadSnapshot.map(j => JsonValue(j.noSpaces))
    val now = LocalDateTime.now()

    executeUpdate(
      """UPDATE sync_queue SET
        |  entity_type = ?, entity_id = ?, operation = ?, status = ?, priority = ?,
        |  queued_at = ?, started_at = ?, completed_at = ?, attempt_count = ?,
        |  next_retry_at = ?, last_error = ?, payload_snapshot = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        SyncEntityType.toDbString(entity.entityType),
        entity.entityId,
        SyncOperation.toDbString(entity.operation),
        QueueStatus.toDbString(entity.status),
        entity.priority,
        entity.queuedAt,
        entity.startedAt,
        entity.completedAt,
        entity.attemptCount,
        entity.nextRetryAt,
        entity.lastError,
        payloadJson,
        now,
        entity.id
      )
    )
    entity.copy(updatedAt = now)

  def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM sync_queue WHERE id = ?", Seq(id)) > 0

  // ============================================
  // Queue Operations
  // ============================================

  /**
   * Enqueue a sync operation.
   * If an entry already exists for this entity/operation, updates it.
   */
  def enqueue(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    priority: Int = 5,
    payloadSnapshot: Option[Json] = None
  )(using conn: Connection): SyncQueueEntity =
    // Check for existing pending entry
    findByEntityAndOperation(entityType, entityId, operation) match
      case Some(existing) if existing.status == QueueStatus.Pending =>
        // Update existing entry
        val updated = existing.copy(
          priority = math.min(existing.priority, priority),
          payloadSnapshot = payloadSnapshot.orElse(existing.payloadSnapshot)
        )
        update(updated)
      case _ =>
        // Create new entry
        val entity = SyncQueueEntity.create(entityType, entityId, operation, priority, payloadSnapshot)
        insert(entity)

  /**
   * Find pending entries ready to process.
   */
  def findPendingBatch(batchSize: Int = 10)(using conn: Connection): List[SyncQueueEntity] =
    queryList(
      """SELECT * FROM sync_queue
        |WHERE status = 'PENDING'
        |AND (next_retry_at IS NULL OR next_retry_at <= ?)
        |ORDER BY priority ASC, queued_at ASC
        |LIMIT ?
      """.stripMargin,
      Seq(LocalDateTime.now(), batchSize)
    )(mapRow)

  /**
   * Find entry by entity and operation.
   */
  def findByEntityAndOperation(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation
  )(using conn: Connection): Option[SyncQueueEntity] =
    queryOne(
      """SELECT * FROM sync_queue
        |WHERE entity_type = ? AND entity_id = ? AND operation = ?
        |AND status IN ('PENDING', 'IN_PROGRESS')
      """.stripMargin,
      Seq(
        SyncEntityType.toDbString(entityType),
        entityId,
        SyncOperation.toDbString(operation)
      )
    )(mapRow)

  /**
   * Find all entries for an entity.
   */
  def findByEntity(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): List[SyncQueueEntity] =
    queryList(
      "SELECT * FROM sync_queue WHERE entity_type = ? AND entity_id = ? ORDER BY queued_at DESC",
      Seq(SyncEntityType.toDbString(entityType), entityId)
    )(mapRow)

  /**
   * Find entries by status.
   */
  def findByStatus(status: QueueStatus)(using conn: Connection): List[SyncQueueEntity] =
    queryList(
      "SELECT * FROM sync_queue WHERE status = ? ORDER BY priority ASC, queued_at ASC",
      Seq(QueueStatus.toDbString(status))
    )(mapRow)

  // ============================================
  // Status Update Operations
  // ============================================

  /**
   * Mark entry as in-progress.
   */
  def markInProgress(id: UUID)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_queue SET
        |  status = 'IN_PROGRESS', started_at = ?, attempt_count = attempt_count + 1, updated_at = ?
        |WHERE id = ? AND status = 'PENDING'
      """.stripMargin,
      Seq(now, now, id)
    ) > 0

  /**
   * Mark entry as completed.
   */
  def markCompleted(id: UUID)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_queue SET
        |  status = 'COMPLETED', completed_at = ?, last_error = NULL, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(now, now, id)
    ) > 0

  /**
   * Mark entry as failed and schedule retry.
   * Uses exponential backoff capped at 1 hour.
   */
  def markFailedWithRetry(id: UUID, error: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    // Get current attempt count for backoff calculation
    findById(id) match
      case Some(entry) =>
        val backoffSeconds = math.min(
          math.pow(2, entry.attemptCount).toLong * 60, // Exponential backoff
          3600 // Cap at 1 hour
        )
        val nextRetry = now.plusSeconds(backoffSeconds)

        executeUpdate(
          """UPDATE sync_queue SET
            |  status = 'PENDING', last_error = ?, next_retry_at = ?, updated_at = ?
            |WHERE id = ?
          """.stripMargin,
          Seq(error, nextRetry, now, id)
        ) > 0
      case None => false

  /**
   * Mark entry as permanently failed.
   */
  def markFailed(id: UUID, error: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_queue SET
        |  status = 'FAILED', completed_at = ?, last_error = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(now, error, now, id)
    ) > 0

  /**
   * Cancel a pending entry.
   */
  def cancel(id: UUID)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_queue SET
        |  status = 'CANCELLED', completed_at = ?, updated_at = ?
        |WHERE id = ? AND status IN ('PENDING', 'FAILED')
      """.stripMargin,
      Seq(now, now, id)
    ) > 0

  /**
   * Cancel all pending entries for an entity.
   */
  def cancelByEntity(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): Int =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_queue SET
        |  status = 'CANCELLED', completed_at = ?, updated_at = ?
        |WHERE entity_type = ? AND entity_id = ? AND status IN ('PENDING', 'FAILED')
      """.stripMargin,
      Seq(now, now, SyncEntityType.toDbString(entityType), entityId)
    )

  // ============================================
  // Statistics
  // ============================================

  /**
   * Count entries by status.
   */
  def countByStatus()(using conn: Connection): Map[QueueStatus, Long] =
    queryList(
      "SELECT status, COUNT(*) as cnt FROM sync_queue GROUP BY status"
    ) { rs =>
      val status = QueueStatus.fromString(rs.getString("status"))
      val count = rs.getLong("cnt")
      (status, count)
    }.toMap

  /**
   * Count pending entries.
   */
  def countPending()(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM sync_queue WHERE status = 'PENDING'"
    )(_.getLong(1)).getOrElse(0L)

  /**
   * Count entries ready to retry now.
   */
  def countReadyToRetry()(using conn: Connection): Long =
    queryOne(
      """SELECT COUNT(*) FROM sync_queue
        |WHERE status = 'PENDING' AND (next_retry_at IS NULL OR next_retry_at <= ?)
      """.stripMargin,
      Seq(LocalDateTime.now())
    )(_.getLong(1)).getOrElse(0L)

  /**
   * Cleanup completed entries older than specified days.
   */
  def cleanupCompleted(olderThanDays: Int)(using conn: Connection): Int =
    val cutoff = LocalDateTime.now().minusDays(olderThanDays)
    executeUpdate(
      "DELETE FROM sync_queue WHERE status IN ('COMPLETED', 'CANCELLED') AND completed_at < ?",
      Seq(cutoff)
    )

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): SyncQueueEntity =
    val payloadJson = getOptJsonString(rs, "payload_snapshot")
    val payload = payloadJson.flatMap(json => parse(json).toOption)

    SyncQueueEntity(
      id = getUUID(rs, "id"),
      entityType = SyncEntityType.fromString(rs.getString("entity_type")),
      entityId = getUUID(rs, "entity_id"),
      operation = SyncOperation.fromString(rs.getString("operation")),
      status = QueueStatus.fromString(rs.getString("status")),
      priority = rs.getInt("priority"),
      queuedAt = getDateTime(rs, "queued_at"),
      startedAt = getOptDateTime(rs, "started_at"),
      completedAt = getOptDateTime(rs, "completed_at"),
      attemptCount = rs.getInt("attempt_count"),
      nextRetryAt = getOptDateTime(rs, "next_retry_at"),
      lastError = getOptString(rs, "last_error"),
      payloadSnapshot = payload,
      createdAt = getDateTime(rs, "created_at"),
      updatedAt = getDateTime(rs, "updated_at")
    )
