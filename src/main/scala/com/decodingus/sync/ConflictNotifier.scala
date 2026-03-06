package com.decodingus.sync

import com.decodingus.repository.SyncConflictEntity
import scalafx.application.Platform
import scalafx.beans.property.{BooleanProperty, IntegerProperty, ObjectProperty}
import scalafx.collections.ObservableBuffer

/**
 * Observable state for sync status UI integration.
 *
 * Provides:
 * - Conflict tracking with observable properties for UI binding
 * - Pending sync count for status bar
 * - Non-blocking notifications (UI thread safe)
 *
 * Usage in ScalaFX:
 * ```
 * val notifier = ConflictNotifier()
 * statusLabel.text <== notifier.pendingCount.asString.concat(" pending")
 * conflictBadge.visible <== notifier.hasConflicts
 * ```
 */
class ConflictNotifier:

  // ============================================
  // Observable Properties for UI Binding
  // ============================================

  /** Whether there are any unresolved conflicts */
  val hasConflicts: BooleanProperty = BooleanProperty(false)

  /** Number of unresolved conflicts */
  val conflictCount: IntegerProperty = IntegerProperty(0)

  /** Number of pending sync operations in queue */
  val pendingCount: IntegerProperty = IntegerProperty(0)

  /** Whether sync is currently in progress */
  val isSyncing: BooleanProperty = BooleanProperty(false)

  /** Whether the service is online (can reach PDS) */
  val isOnline: BooleanProperty = BooleanProperty(true)

  /** Most recent error message (if any) */
  val lastError: ObjectProperty[Option[String]] = ObjectProperty(None)

  /** List of current conflicts for UI display */
  val conflicts: ObservableBuffer[SyncConflictEntity] = ObservableBuffer.empty

  // ============================================
  // Notification Methods
  // ============================================

  /**
   * Update counts for pending syncs and conflicts.
   * Called from background thread, runs on FX thread.
   */
  def updateCounts(pending: Int, conflictCount: Int): Unit =
    runOnFxThread {
      this.pendingCount.value = pending
      this.conflictCount.value = conflictCount
      this.hasConflicts.value = conflictCount > 0
    }

  /**
   * Notify about new conflicts detected from remote.
   * Replaces current conflict list and updates counts.
   */
  def notifyConflicts(newConflicts: List[SyncConflictEntity]): Unit =
    runOnFxThread {
      conflicts.clear()
      conflicts ++= newConflicts
      conflictCount.value = newConflicts.size
      hasConflicts.value = newConflicts.nonEmpty
    }

  /**
   * Add a single conflict to the list.
   */
  def addConflict(conflict: SyncConflictEntity): Unit =
    runOnFxThread {
      // Remove existing conflict for same entity if present
      conflicts.filterInPlace(c =>
        !(c.entityType == conflict.entityType && c.entityId == conflict.entityId)
      )
      conflicts += conflict
      conflictCount.value = conflicts.size
      hasConflicts.value = conflicts.nonEmpty
    }

  /**
   * Remove a resolved conflict from the list.
   */
  def clearConflict(entityId: java.util.UUID): Unit =
    runOnFxThread {
      conflicts.filterInPlace(_.entityId != entityId)
      conflictCount.value = conflicts.size
      hasConflicts.value = conflicts.nonEmpty
    }

  /**
   * Clear all conflicts.
   */
  def clearAllConflicts(): Unit =
    runOnFxThread {
      conflicts.clear()
      conflictCount.value = 0
      hasConflicts.value = false
    }

  /**
   * Update sync status indicators.
   */
  def setSyncStatus(syncing: Boolean, online: Boolean): Unit =
    runOnFxThread {
      isSyncing.value = syncing
      isOnline.value = online
    }

  /**
   * Set an error message.
   */
  def setError(error: String): Unit =
    runOnFxThread {
      lastError.value = Some(error)
    }

  /**
   * Clear error message.
   */
  def clearError(): Unit =
    runOnFxThread {
      lastError.value = None
    }

  /**
   * Increment pending count (when queueing new sync).
   */
  def incrementPending(): Unit =
    runOnFxThread {
      pendingCount.value = pendingCount.value + 1
    }

  /**
   * Decrement pending count (when sync completes).
   */
  def decrementPending(): Unit =
    runOnFxThread {
      pendingCount.value = math.max(0, pendingCount.value - 1)
    }

  // ============================================
  // Helper
  // ============================================

  private def runOnFxThread(action: => Unit): Unit =
    try
      if Platform.isFxApplicationThread then
        action
      else
        Platform.runLater(action)
    catch
      case _: RuntimeException =>
        // No JavaFX toolkit available (headless tests) - run directly
        action

/**
 * Factory for ConflictNotifier.
 */
object ConflictNotifier:
  def apply(): ConflictNotifier = new ConflictNotifier()
