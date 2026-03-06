package com.decodingus.workspace.model

/** Represents the synchronization status between local and remote (PDS) workspace state */
sealed trait SyncStatus

object SyncStatus {
  /** Local and remote are in sync */
  case object Synced extends SyncStatus

  /** Currently syncing to remote */
  case object Syncing extends SyncStatus

  /** Local changes pending sync to remote */
  case object Pending extends SyncStatus

  /** Sync failed - check lastSyncError for details */
  case object Error extends SyncStatus

  /** Offline mode - remote unavailable, using local only */
  case object Offline extends SyncStatus
}
