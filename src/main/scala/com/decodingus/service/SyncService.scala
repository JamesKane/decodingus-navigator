package com.decodingus.service

import com.decodingus.db.Transactor
import com.decodingus.repository.*
import io.circe.Json

import java.time.LocalDateTime
import java.util.UUID

/**
 * Service for managing PDS synchronization.
 *
 * Provides high-level operations for:
 * - Enqueueing sync operations (event-driven outgoing)
 * - Processing the sync queue
 * - Managing conflicts
 * - Tracking sync history
 * - Sync statistics and status
 */
trait SyncService:

  // ============================================
  // Queue Operations
  // ============================================

  /**
   * Enqueue an entity for sync to PDS.
   * Called when local changes are made.
   */
  def enqueuePush(
                   entityType: SyncEntityType,
                   entityId: UUID,
                   operation: SyncOperation,
                   priority: Int = 5,
                   payloadSnapshot: Option[Json] = None
                 ): Either[String, SyncQueueEntity]

  /**
   * Get pending sync count.
   */
  def getPendingCount(): Either[String, Long]

  /**
   * Get next batch of items to sync.
   */
  def getNextBatch(batchSize: Int = 10): Either[String, List[SyncQueueEntity]]

  /**
   * Start processing a queued item.
   */
  def startProcessing(queueId: UUID): Either[String, Boolean]

  /**
   * Mark a queued item as successfully synced.
   */
  def markSynced(
                  queueId: UUID,
                  atUri: Option[String] = None,
                  remoteCid: Option[String] = None
                ): Either[String, Boolean]

  /**
   * Mark a queued item as failed (will retry).
   */
  def markFailed(queueId: UUID, error: String): Either[String, Boolean]

  /**
   * Mark a queued item as permanently failed.
   */
  def markPermanentlyFailed(queueId: UUID, error: String): Either[String, Boolean]

  /**
   * Cancel a pending sync operation.
   */
  def cancelSync(queueId: UUID): Either[String, Boolean]

  /**
   * Cancel all pending syncs for an entity.
   */
  def cancelAllForEntity(entityType: SyncEntityType, entityId: UUID): Either[String, Int]

  // ============================================
  // Conflict Operations
  // ============================================

  /**
   * Record a sync conflict.
   */
  def recordConflict(
                      entityType: SyncEntityType,
                      entityId: UUID,
                      localVersion: Int,
                      remoteVersion: Int,
                      atUri: Option[String] = None,
                      localChanges: Option[Json] = None,
                      remoteChanges: Option[Json] = None,
                      overlappingFields: Option[Json] = None,
                      localSnapshot: Option[Json] = None,
                      remoteSnapshot: Option[Json] = None,
                      suggestedResolution: Option[ConflictResolution] = None,
                      resolutionReason: Option[String] = None
                    ): Either[String, SyncConflictEntity]

  /**
   * Get unresolved conflicts.
   */
  def getUnresolvedConflicts(): Either[String, List[SyncConflictEntity]]

  /**
   * Get unresolved conflict count.
   */
  def getUnresolvedConflictCount(): Either[String, Long]

  /**
   * Get conflict by ID.
   */
  def getConflict(conflictId: UUID): Either[String, Option[SyncConflictEntity]]

  /**
   * Resolve a conflict by keeping local version.
   */
  def resolveKeepLocal(conflictId: UUID, resolvedBy: String = "user"): Either[String, Boolean]

  /**
   * Resolve a conflict by accepting remote version.
   */
  def resolveAcceptRemote(conflictId: UUID, resolvedBy: String = "user"): Either[String, Boolean]

  /**
   * Resolve a conflict by merging.
   */
  def resolveMerge(conflictId: UUID, resolvedBy: String = "user"): Either[String, Boolean]

  /**
   * Dismiss a conflict.
   */
  def dismissConflict(conflictId: UUID, resolvedBy: String = "user"): Either[String, Boolean]

  // ============================================
  // History Operations
  // ============================================

  /**
   * Get recent sync history.
   */
  def getRecentHistory(limit: Int = 100): Either[String, List[SyncHistoryEntity]]

  /**
   * Get sync history for an entity.
   */
  def getEntityHistory(entityType: SyncEntityType, entityId: UUID): Either[String, List[SyncHistoryEntity]]

  /**
   * Get last successful sync for an entity.
   */
  def getLastSuccessfulSync(entityType: SyncEntityType, entityId: UUID): Either[String, Option[SyncHistoryEntity]]

  /**
   * Get recent failures.
   */
  def getRecentFailures(limit: Int = 50): Either[String, List[SyncHistoryEntity]]

  // ============================================
  // Status and Statistics
  // ============================================

  /**
   * Get overall sync status summary.
   */
  def getSyncStatus(): Either[String, SyncStatusInfo]

  /**
   * Get sync statistics for a time period.
   */
  def getStats(start: LocalDateTime, end: LocalDateTime): Either[String, SyncStats]

  /**
   * Check if sync is healthy (no errors, queue not too large).
   */
  def isSyncHealthy(): Either[String, Boolean]

  // ============================================
  // Maintenance Operations
  // ============================================

  /**
   * Cleanup completed queue entries.
   */
  def cleanupQueue(olderThanDays: Int = 7): Either[String, Int]

  /**
   * Cleanup old history entries.
   */
  def cleanupHistory(olderThanDays: Int = 90): Either[String, Int]

  /**
   * Cleanup resolved conflicts.
   */
  def cleanupConflicts(olderThanDays: Int = 30): Either[String, Int]

/**
 * Overall sync status information.
 */
case class SyncStatusInfo(
                           pendingCount: Long,
                           inProgressCount: Long,
                           failedCount: Long,
                           unresolvedConflicts: Long,
                           lastSuccessfulSync: Option[LocalDateTime],
                           lastSyncAttempt: Option[LocalDateTime],
                           isHealthy: Boolean,
                           healthMessage: Option[String]
                         )

/**
 * H2 database-backed implementation of SyncService.
 */
class H2SyncService(
                     transactor: Transactor,
                     queueRepo: SyncQueueRepository,
                     historyRepo: SyncHistoryRepository,
                     conflictRepo: SyncConflictRepository
                   ) extends SyncService:

  // ============================================
  // Queue Operations
  // ============================================

  override def enqueuePush(
                            entityType: SyncEntityType,
                            entityId: UUID,
                            operation: SyncOperation,
                            priority: Int,
                            payloadSnapshot: Option[Json]
                          ): Either[String, SyncQueueEntity] =
    transactor.readWrite {
      queueRepo.enqueue(entityType, entityId, operation, priority, payloadSnapshot)
    }

  override def getPendingCount(): Either[String, Long] =
    transactor.readOnly {
      queueRepo.countPending()
    }

  override def getNextBatch(batchSize: Int): Either[String, List[SyncQueueEntity]] =
    transactor.readOnly {
      queueRepo.findPendingBatch(batchSize)
    }

  override def startProcessing(queueId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      queueRepo.markInProgress(queueId)
    }

  override def markSynced(
                           queueId: UUID,
                           atUri: Option[String],
                           remoteCid: Option[String]
                         ): Either[String, Boolean] =
    transactor.readWrite {
      // Get the queue entry to record in history
      queueRepo.findById(queueId) match
        case Some(entry) =>
          // Record success in history
          historyRepo.recordSuccess(
            entityType = entry.entityType,
            entityId = entry.entityId,
            operation = entry.operation,
            direction = SyncDirection.Push,
            startedAt = entry.startedAt.getOrElse(entry.queuedAt),
            atUri = atUri,
            remoteCid = remoteCid
          )
          // Mark queue entry as completed
          queueRepo.markCompleted(queueId)
        case None =>
          false
    }

  override def markFailed(queueId: UUID, error: String): Either[String, Boolean] =
    transactor.readWrite {
      queueRepo.findById(queueId) match
        case Some(entry) =>
          // Record failure in history
          historyRepo.recordFailure(
            entityType = entry.entityType,
            entityId = entry.entityId,
            operation = entry.operation,
            direction = SyncDirection.Push,
            startedAt = entry.startedAt.getOrElse(entry.queuedAt),
            errorMessage = error
          )
          // Mark for retry
          queueRepo.markFailedWithRetry(queueId, error)
        case None =>
          false
    }

  override def markPermanentlyFailed(queueId: UUID, error: String): Either[String, Boolean] =
    transactor.readWrite {
      queueRepo.findById(queueId) match
        case Some(entry) =>
          // Record failure in history
          historyRepo.recordFailure(
            entityType = entry.entityType,
            entityId = entry.entityId,
            operation = entry.operation,
            direction = SyncDirection.Push,
            startedAt = entry.startedAt.getOrElse(entry.queuedAt),
            errorMessage = error
          )
          // Mark as permanently failed
          queueRepo.markFailed(queueId, error)
        case None =>
          false
    }

  override def cancelSync(queueId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      queueRepo.cancel(queueId)
    }

  override def cancelAllForEntity(entityType: SyncEntityType, entityId: UUID): Either[String, Int] =
    transactor.readWrite {
      queueRepo.cancelByEntity(entityType, entityId)
    }

  // ============================================
  // Conflict Operations
  // ============================================

  override def recordConflict(
                               entityType: SyncEntityType,
                               entityId: UUID,
                               localVersion: Int,
                               remoteVersion: Int,
                               atUri: Option[String],
                               localChanges: Option[Json],
                               remoteChanges: Option[Json],
                               overlappingFields: Option[Json],
                               localSnapshot: Option[Json],
                               remoteSnapshot: Option[Json],
                               suggestedResolution: Option[ConflictResolution],
                               resolutionReason: Option[String]
                             ): Either[String, SyncConflictEntity] =
    transactor.readWrite {
      val entity = SyncConflictEntity.create(
        entityType = entityType,
        entityId = entityId,
        localVersion = localVersion,
        remoteVersion = remoteVersion,
        atUri = atUri,
        localChanges = localChanges,
        remoteChanges = remoteChanges,
        overlappingFields = overlappingFields,
        suggestedResolution = suggestedResolution,
        resolutionReason = resolutionReason,
        localSnapshot = localSnapshot,
        remoteSnapshot = remoteSnapshot
      )
      conflictRepo.insert(entity)
    }

  override def getUnresolvedConflicts(): Either[String, List[SyncConflictEntity]] =
    transactor.readOnly {
      conflictRepo.findUnresolved()
    }

  override def getUnresolvedConflictCount(): Either[String, Long] =
    transactor.readOnly {
      conflictRepo.countUnresolved()
    }

  override def getConflict(conflictId: UUID): Either[String, Option[SyncConflictEntity]] =
    transactor.readOnly {
      conflictRepo.findById(conflictId)
    }

  override def resolveKeepLocal(conflictId: UUID, resolvedBy: String): Either[String, Boolean] =
    transactor.readWrite {
      conflictRepo.resolveKeepLocal(conflictId, resolvedBy)
    }

  override def resolveAcceptRemote(conflictId: UUID, resolvedBy: String): Either[String, Boolean] =
    transactor.readWrite {
      conflictRepo.resolveAcceptRemote(conflictId, resolvedBy)
    }

  override def resolveMerge(conflictId: UUID, resolvedBy: String): Either[String, Boolean] =
    transactor.readWrite {
      conflictRepo.resolveMerge(conflictId, resolvedBy)
    }

  override def dismissConflict(conflictId: UUID, resolvedBy: String): Either[String, Boolean] =
    transactor.readWrite {
      conflictRepo.dismiss(conflictId, resolvedBy)
    }

  // ============================================
  // History Operations
  // ============================================

  override def getRecentHistory(limit: Int): Either[String, List[SyncHistoryEntity]] =
    transactor.readOnly {
      historyRepo.findAll(limit)
    }

  override def getEntityHistory(entityType: SyncEntityType, entityId: UUID): Either[String, List[SyncHistoryEntity]] =
    transactor.readOnly {
      historyRepo.findByEntity(entityType, entityId)
    }

  override def getLastSuccessfulSync(entityType: SyncEntityType, entityId: UUID): Either[String, Option[SyncHistoryEntity]] =
    transactor.readOnly {
      historyRepo.getLastSuccessfulSync(entityType, entityId)
    }

  override def getRecentFailures(limit: Int): Either[String, List[SyncHistoryEntity]] =
    transactor.readOnly {
      historyRepo.findRecentFailures(limit)
    }

  // ============================================
  // Status and Statistics
  // ============================================

  override def getSyncStatus(): Either[String, SyncStatusInfo] =
    transactor.readOnly {
      val queueStatusCounts = queueRepo.countByStatus()
      val pendingCount = queueStatusCounts.getOrElse(QueueStatus.Pending, 0L)
      val inProgressCount = queueStatusCounts.getOrElse(QueueStatus.InProgress, 0L)
      val failedCount = queueStatusCounts.getOrElse(QueueStatus.Failed, 0L)
      val unresolvedConflicts = conflictRepo.countUnresolved()

      // Get last sync times from history
      val recentHistory = historyRepo.findAll(1)
      val lastSyncAttempt = recentHistory.headOption.map(_.completedAt)

      val successfulHistory = historyRepo.findByStatus(SyncResultStatus.Success, 1)
      val lastSuccessfulSync = successfulHistory.headOption.map(_.completedAt)

      // Determine health
      val (isHealthy, healthMessage) = determineHealth(pendingCount, failedCount, unresolvedConflicts)

      SyncStatusInfo(
        pendingCount = pendingCount,
        inProgressCount = inProgressCount,
        failedCount = failedCount,
        unresolvedConflicts = unresolvedConflicts,
        lastSuccessfulSync = lastSuccessfulSync,
        lastSyncAttempt = lastSyncAttempt,
        isHealthy = isHealthy,
        healthMessage = healthMessage
      )
    }

  override def getStats(start: LocalDateTime, end: LocalDateTime): Either[String, SyncStats] =
    transactor.readOnly {
      historyRepo.getStatsForPeriod(start, end)
    }

  override def isSyncHealthy(): Either[String, Boolean] =
    getSyncStatus().map(_.isHealthy)

  private def determineHealth(pending: Long, failed: Long, conflicts: Long): (Boolean, Option[String]) =
    if conflicts > 0 then
      (false, Some(s"$conflicts unresolved conflict(s) require attention"))
    else if failed > 10 then
      (false, Some(s"$failed permanently failed sync operations"))
    else if pending > 100 then
      (false, Some(s"Large sync backlog: $pending pending operations"))
    else
      (true, None)

  // ============================================
  // Maintenance Operations
  // ============================================

  override def cleanupQueue(olderThanDays: Int): Either[String, Int] =
    transactor.readWrite {
      queueRepo.cleanupCompleted(olderThanDays)
    }

  override def cleanupHistory(olderThanDays: Int): Either[String, Int] =
    transactor.readWrite {
      historyRepo.cleanupOlderThan(olderThanDays)
    }

  override def cleanupConflicts(olderThanDays: Int): Either[String, Int] =
    transactor.readWrite {
      conflictRepo.cleanupResolved(olderThanDays)
    }
