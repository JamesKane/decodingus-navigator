package com.decodingus.ui.components

import com.decodingus.sync.ConflictNotifier
import scalafx.Includes._
import scalafx.application.Platform
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.{Button, Label, Tooltip}
import scalafx.scene.layout.{HBox, Priority, Region}

/**
 * Application status bar displayed at the bottom of the main window.
 *
 * Shows:
 * - Connection status (online/offline)
 * - Pending sync count
 * - Conflict indicator (if any)
 * - Cache status summary
 *
 * Binds to ConflictNotifier observable properties for reactive updates.
 */
class StatusBar(notifier: ConflictNotifier) extends HBox {
  alignment = Pos.CenterLeft
  padding = Insets(4, 10, 4, 10)
  spacing = 15
  style = "-fx-background-color: #2D2D2D; -fx-border-color: #404040; -fx-border-width: 1 0 0 0;"

  // Connection status indicator
  private val connectionLabel = new Label("● Online") {
    style = "-fx-font-size: 11px; -fx-text-fill: #4CAF50;"
    tooltip = Tooltip("Connected and ready to sync")
  }

  // Sync status indicator
  private val syncLabel = new Label("") {
    style = "-fx-font-size: 11px; -fx-text-fill: #9E9E9E;"
    visible = false
    managed <== visible
  }

  // Conflict indicator (hidden when no conflicts)
  private val conflictLabel = new Label("") {
    style = "-fx-font-size: 11px; -fx-text-fill: #FF9800; -fx-cursor: hand;"
    visible = false
    managed <== visible
    tooltip = Tooltip("Click to view conflicts")
  }

  // Spacer to push right-side content
  private val spacer = new Region()
  HBox.setHgrow(spacer, Priority.Always)

  // Right-side status (version, cache size, etc.)
  private val cacheLabel = new Label("") {
    style = "-fx-font-size: 11px; -fx-text-fill: #757575;"
  }

  children = Seq(connectionLabel, syncLabel, conflictLabel, spacer, cacheLabel)

  // Bind to ConflictNotifier properties
  setupBindings()

  private def setupBindings(): Unit = {
    // Online/offline status
    notifier.isOnline.onChange { (_, _, isOnline) =>
      Platform.runLater {
        if isOnline then
          connectionLabel.text = "● Online"
          connectionLabel.style = "-fx-font-size: 11px; -fx-text-fill: #4CAF50;"
          connectionLabel.tooltip = Tooltip("Connected and ready to sync")
        else
          connectionLabel.text = "○ Offline"
          connectionLabel.style = "-fx-font-size: 11px; -fx-text-fill: #9E9E9E;"
          connectionLabel.tooltip = Tooltip("Working offline - changes will sync when connected")
      }
    }

    // Syncing status
    notifier.isSyncing.onChange { (_, _, isSyncing) =>
      Platform.runLater {
        if isSyncing then
          syncLabel.text = "↻ Syncing..."
          syncLabel.style = "-fx-font-size: 11px; -fx-text-fill: #2196F3;"
          syncLabel.visible = true
        else
          syncLabel.visible = false
      }
    }

    // Pending count
    notifier.pendingCount.onChange { (_, _, count) =>
      Platform.runLater {
        val pending = count.intValue()
        if pending > 0 && !notifier.isSyncing.value then
          syncLabel.text = s"$pending pending"
          syncLabel.style = "-fx-font-size: 11px; -fx-text-fill: #FF9800;"
          syncLabel.tooltip = Tooltip(s"$pending changes waiting to sync")
          syncLabel.visible = true
        else if !notifier.isSyncing.value then
          syncLabel.visible = false
      }
    }

    // Conflict indicator
    notifier.hasConflicts.onChange { (_, _, hasConflicts) =>
      Platform.runLater {
        conflictLabel.visible = hasConflicts
      }
    }

    notifier.conflictCount.onChange { (_, _, count) =>
      Platform.runLater {
        val conflicts = count.intValue()
        if conflicts > 0 then
          conflictLabel.text = s"⚠ $conflicts conflict${if conflicts != 1 then "s" else ""}"
          conflictLabel.tooltip = Tooltip(s"$conflicts sync conflict${if conflicts != 1 then "s" else ""} need resolution")
      }
    }

    // Error status
    notifier.lastError.onChange { (_, _, errorOpt) =>
      Platform.runLater {
        errorOpt match
          case Some(error) =>
            syncLabel.text = "⚠ Sync error"
            syncLabel.style = "-fx-font-size: 11px; -fx-text-fill: #F44336;"
            syncLabel.tooltip = Tooltip(error)
            syncLabel.visible = true
          case None =>
            // Handled by other listeners
      }
    }
  }

  /**
   * Update cache statistics display.
   */
  def updateCacheStats(artifactCount: Int, cacheSize: Long): Unit =
    Platform.runLater {
      val sizeStr = formatBytes(cacheSize)
      cacheLabel.text = s"Cache: $artifactCount artifacts ($sizeStr)"
    }

  /**
   * Set a custom right-side message (e.g., version info).
   */
  def setRightMessage(message: String): Unit =
    Platform.runLater {
      cacheLabel.text = message
    }

  private def formatBytes(bytes: Long): String =
    if bytes < 1024 then s"$bytes B"
    else if bytes < 1024 * 1024 then f"${bytes / 1024.0}%.1f KB"
    else if bytes < 1024 * 1024 * 1024 then f"${bytes / (1024.0 * 1024)}%.1f MB"
    else f"${bytes / (1024.0 * 1024 * 1024)}%.1f GB"
}

object StatusBar:
  def apply(notifier: ConflictNotifier): StatusBar = new StatusBar(notifier)
