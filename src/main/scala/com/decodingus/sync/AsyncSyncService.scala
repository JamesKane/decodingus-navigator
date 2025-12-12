package com.decodingus.sync

import com.decodingus.auth.User
import com.decodingus.config.FeatureToggles
import com.decodingus.db.Transactor
import com.decodingus.repository.*
import com.decodingus.util.Logger
import io.circe.Json

import java.sql.Connection
import java.time.LocalDateTime
import java.util.UUID
import java.util.concurrent.{Executors, ScheduledExecutorService, ScheduledFuture, TimeUnit}
import java.util.concurrent.atomic.AtomicBoolean
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
 * Async sync service for background synchronization with Personal Data Store (PDS).
 *
 * Responsibilities:
 * - Queue entities for async sync (event-driven on edit)
 * - Process outgoing sync queue in background
 * - Optional hourly pull for remote changes
 * - Detect and track conflicts via SyncConflictRepository
 * - Notify UI of sync status via ConflictNotifier
 *
 * Design principles:
 * - Non-blocking: User continues working, sync happens in background
 * - Indefinite offline: Queue persists forever, no timeout warnings
 * - Exponential backoff: Failed syncs retry with increasing delays (capped at 1 hour)
 */
class AsyncSyncService(
  transactor: Transactor,
  syncQueueRepo: SyncQueueRepository,
  syncConflictRepo: SyncConflictRepository,
  syncHistoryRepo: SyncHistoryRepository,
  biosampleRepo: BiosampleRepository,
  projectRepo: ProjectRepository,
  sequenceRunRepo: SequenceRunRepository,
  alignmentRepo: AlignmentRepository,
  conflictNotifier: ConflictNotifier
)(using ec: ExecutionContext):

  private val log = Logger[AsyncSyncService]
  private val scheduler: ScheduledExecutorService = Executors.newScheduledThreadPool(2)
  private val isProcessing = AtomicBoolean(false)
  private val incomingSyncEnabled = AtomicBoolean(true)
  private var currentUser: Option[User] = None

  // Scheduled task handles for cleanup
  private var incomingSyncTask: Option[ScheduledFuture[?]] = None
  private var queueProcessorTask: Option[ScheduledFuture[?]] = None

  // ============================================
  // Lifecycle Management
  // ============================================

  /**
   * Start the sync service with the given user.
   * Begins background queue processing and optional incoming sync polling.
   */
  def start(user: Option[User]): Unit =
    currentUser = user

    // Start queue processor (runs every 30 seconds)
    queueProcessorTask = Some(scheduler.scheduleWithFixedDelay(
      () => processOutgoingQueueSafe(),
      5,    // Initial delay: 5 seconds
      30,   // Period: 30 seconds
      TimeUnit.SECONDS
    ))

    // Start incoming sync poller (runs hourly, if enabled)
    incomingSyncTask = Some(scheduler.scheduleWithFixedDelay(
      () => if incomingSyncEnabled.get() then pullRemoteChangesSafe(),
      60,   // Initial delay: 1 minute
      3600, // Period: 1 hour
      TimeUnit.SECONDS
    ))

    log.info("Started background sync service")

  /**
   * Stop the sync service.
   * Allows in-progress operations to complete gracefully.
   */
  def shutdown(): Unit =
    log.info("Shutting down...")
    incomingSyncTask.foreach(_.cancel(false))
    queueProcessorTask.foreach(_.cancel(false))
    scheduler.shutdown()
    if !scheduler.awaitTermination(30, TimeUnit.SECONDS) then
      scheduler.shutdownNow()
    log.info("Shutdown complete")

  /**
   * Enable or disable incoming sync polling.
   * User preference - can work offline indefinitely.
   */
  def setIncomingSyncEnabled(enabled: Boolean): Unit =
    incomingSyncEnabled.set(enabled)
    log.info(s"Incoming sync ${if enabled then "enabled" else "disabled"}")

  /**
   * Update the current user (e.g., after login/logout).
   */
  def setUser(user: Option[User]): Unit =
    currentUser = user

  // ============================================
  // Queueing Operations
  // ============================================

  /**
   * Queue an entity for async sync to PDS.
   * Called immediately when user makes an edit.
   * Returns immediately; actual sync happens in background.
   *
   * @param entityType The type of entity being synced
   * @param entityId The entity's unique identifier
   * @param operation The sync operation (Create, Update, Delete)
   * @param priority Lower numbers = higher priority (default: 5)
   * @param payload Optional JSON snapshot of entity state
   */
  def queueForSync(
    entityType: SyncEntityType,
    entityId: UUID,
    operation: SyncOperation,
    priority: Int = 5,
    payload: Option[Json] = None
  ): Future[SyncQueueEntity] = Future {
    transactor.readWrite {
      val entry = syncQueueRepo.enqueue(entityType, entityId, operation, priority, payload)
      log.debug(s"Queued $operation for $entityType:$entityId")
      entry
    }.getOrElse {
      throw new RuntimeException("Failed to enqueue sync operation")
    }
  }.andThen {
    case Success(_) =>
      // Trigger immediate processing attempt
      processOutgoingQueueAsync()
      updatePendingCount()
    case Failure(e) =>
      log.error(s"Failed to queue sync: ${e.getMessage}")
  }

  /**
   * Convenience method for queueing biosample sync.
   */
  def queueBiosampleSync(biosampleId: UUID, operation: SyncOperation): Future[SyncQueueEntity] =
    queueForSync(SyncEntityType.Biosample, biosampleId, operation)

  /**
   * Convenience method for queueing project sync.
   */
  def queueProjectSync(projectId: UUID, operation: SyncOperation): Future[SyncQueueEntity] =
    queueForSync(SyncEntityType.Project, projectId, operation)

  /**
   * Convenience method for queueing sequence run sync.
   */
  def queueSequenceRunSync(sequenceRunId: UUID, operation: SyncOperation): Future[SyncQueueEntity] =
    queueForSync(SyncEntityType.SequenceRun, sequenceRunId, operation)

  /**
   * Convenience method for queueing alignment sync.
   */
  def queueAlignmentSync(alignmentId: UUID, operation: SyncOperation): Future[SyncQueueEntity] =
    queueForSync(SyncEntityType.Alignment, alignmentId, operation)

  // ============================================
  // Queue Processing
  // ============================================

  /**
   * Manually trigger queue processing.
   * Used when user explicitly requests sync (e.g., "Sync Now" button).
   */
  def processQueueNow(): Future[Int] =
    if isProcessing.compareAndSet(false, true) then
      Future {
        try
          processOutgoingQueue()
        finally
          isProcessing.set(false)
          updatePendingCount()
      }
    else
      Future.successful(0)

  /**
   * Get current queue statistics.
   */
  def getQueueStats: Future[QueueStats] = Future {
    transactor.readOnly {
      val counts = syncQueueRepo.countByStatus()
      QueueStats(
        pending = counts.getOrElse(QueueStatus.Pending, 0L),
        inProgress = counts.getOrElse(QueueStatus.InProgress, 0L),
        failed = counts.getOrElse(QueueStatus.Failed, 0L),
        completed = counts.getOrElse(QueueStatus.Completed, 0L)
      )
    }.getOrElse(QueueStats(0, 0, 0, 0))
  }

  private def processOutgoingQueueSafe(): Unit =
    try
      if isProcessing.compareAndSet(false, true) then
        try processOutgoingQueue()
        finally isProcessing.set(false)
    catch
      case e: Exception =>
        log.error(s"Queue processing error: ${e.getMessage}")

  private def processOutgoingQueueAsync(): Unit =
    Future(processOutgoingQueueSafe())

  private def processOutgoingQueue(): Int =
    if !FeatureToggles.atProtocolEnabled then
      return 0

    currentUser match
      case None =>
        // Offline - leave queue intact for later
        0
      case Some(user) =>
        var processed = 0
        var continueProcessing = true

        while continueProcessing do
          val batch = transactor.readWrite {
            syncQueueRepo.findPendingBatch(10)
          }.getOrElse(List.empty)

          if batch.isEmpty then
            continueProcessing = false
          else
            batch.foreach { entry =>
              processEntry(entry, user) match
                case Success(_) =>
                  processed += 1
                case Failure(e) =>
                  log.warn(s"Failed to process ${entry.entityType}:${entry.entityId}: ${e.getMessage}")
            }

        if processed > 0 then
          log.debug(s"Processed $processed queue entries")
        processed

  private def processEntry(entry: SyncQueueEntity, user: User): Try[Unit] =
    // Mark as in-progress
    transactor.readWrite {
      syncQueueRepo.markInProgress(entry.id)
    }

    val result = Try {
      entry.operation match
        case SyncOperation.Create => pushCreate(entry, user)
        case SyncOperation.Update => pushUpdate(entry, user)
        case SyncOperation.Delete => pushDelete(entry, user)
    }

    result match
      case Success(_) =>
        transactor.readWrite {
          syncQueueRepo.markCompleted(entry.id)
          recordSyncHistory(entry, SyncDirection.Push, SyncResultStatus.Success, None)
        }
      case Failure(e) =>
        transactor.readWrite {
          syncQueueRepo.markFailedWithRetry(entry.id, e.getMessage)
          recordSyncHistory(entry, SyncDirection.Push, SyncResultStatus.Failed, Some(e.getMessage))
        }

    result

  // ============================================
  // PDS Operations (Stubs - implement with actual PDS client)
  // ============================================

  private def pushCreate(entry: SyncQueueEntity, user: User): Unit =
    // TODO: Implement actual PDS create
    // For now, simulate success
    log.debug(s"Would push CREATE for ${entry.entityType}:${entry.entityId} to PDS")
    // PdsClient.createRecord(user, entry.entityType, getEntityPayload(entry))

  private def pushUpdate(entry: SyncQueueEntity, user: User): Unit =
    // TODO: Implement actual PDS update
    log.debug(s"Would push UPDATE for ${entry.entityType}:${entry.entityId} to PDS")
    // PdsClient.updateRecord(user, entry.entityType, entry.entityId, getEntityPayload(entry))

  private def pushDelete(entry: SyncQueueEntity, user: User): Unit =
    // TODO: Implement actual PDS delete
    log.debug(s"Would push DELETE for ${entry.entityType}:${entry.entityId} to PDS")
    // PdsClient.deleteRecord(user, entry.entityType, entry.entityId)

  // ============================================
  // Incoming Sync (Remote Changes)
  // ============================================

  private def pullRemoteChangesSafe(): Unit =
    try pullRemoteChanges()
    catch
      case e: Exception =>
        log.error(s"Pull remote changes error: ${e.getMessage}")

  private def pullRemoteChanges(): Unit =
    if !FeatureToggles.atProtocolEnabled then return

    currentUser match
      case None => // Offline
      case Some(user) =>
        log.debug("Checking for remote changes...")
        // TODO: Implement actual PDS fetch and conflict detection
        // 1. Fetch remote records since last sync
        // 2. Compare with local versions
        // 3. Detect conflicts
        // 4. Create SyncConflictEntity for each conflict
        // 5. Notify UI via conflictNotifier

        val conflicts = detectRemoteConflicts(user)
        if conflicts.nonEmpty then
          conflictNotifier.notifyConflicts(conflicts)
          log.info(s"Detected ${conflicts.size} conflicts")

  private def detectRemoteConflicts(user: User): List[SyncConflictEntity] =
    // TODO: Implement actual conflict detection
    // For now, return empty list (no conflicts)
    List.empty

  // ============================================
  // History Recording
  // ============================================

  private def recordSyncHistory(
    entry: SyncQueueEntity,
    direction: SyncDirection,
    status: SyncResultStatus,
    error: Option[String]
  )(using conn: Connection): Unit =
    val now = LocalDateTime.now()
    syncHistoryRepo.insert(SyncHistoryEntity.create(
      entityType = entry.entityType,
      entityId = entry.entityId,
      operation = entry.operation,
      direction = direction,
      status = status,
      startedAt = entry.startedAt.getOrElse(now),
      completedAt = now,
      errorMessage = error
    ))

  // ============================================
  // UI Integration
  // ============================================

  private def updatePendingCount(): Unit =
    Future {
      transactor.readOnly {
        val pending = syncQueueRepo.countPending()
        val conflicts = syncConflictRepo.countUnresolved()
        conflictNotifier.updateCounts(pending.toInt, conflicts.toInt)
      }
    }

/**
 * Queue statistics.
 */
case class QueueStats(
  pending: Long,
  inProgress: Long,
  failed: Long,
  completed: Long
):
  def total: Long = pending + inProgress + failed + completed

