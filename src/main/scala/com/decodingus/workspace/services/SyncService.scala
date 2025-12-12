package com.decodingus.workspace.services

import com.decodingus.auth.User
import com.decodingus.config.FeatureToggles
import com.decodingus.pds.PdsClient
import com.decodingus.util.Logger
import com.decodingus.workspace.WorkspaceService
import com.decodingus.workspace.model.{Workspace, SyncStatus}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

/**
 * Result of a sync operation.
 */
sealed trait SyncResult
object SyncResult {
  case object Success extends SyncResult
  case object NoChange extends SyncResult
  case object Offline extends SyncResult
  case object Disabled extends SyncResult
  case class Error(message: String) extends SyncResult
}

/**
 * Callback interface for sync status updates.
 */
trait SyncStatusListener {
  def onSyncStatusChanged(status: SyncStatus): Unit
  def onSyncError(message: String): Unit
}

/**
 * Service for synchronizing workspace data with Personal Data Store (PDS).
 *
 * Handles:
 * - Loading workspace from PDS
 * - Saving workspace to PDS
 * - Sync status management
 * - Offline/online mode detection
 */
class SyncService(
  workspaceService: WorkspaceService
)(implicit ec: ExecutionContext) {

  private val log = Logger[SyncService]

  /**
   * Attempts to load workspace from PDS if AT Protocol is enabled and user is logged in.
   *
   * @param user Current logged-in user (if any)
   * @param currentWorkspace The current local workspace for comparison
   * @param onStatusChange Callback for status updates
   * @return Future containing the workspace (either remote or current) and sync result
   */
  def syncFromPds(
    user: Option[User],
    currentWorkspace: Workspace,
    onStatusChange: SyncStatus => Unit = _ => ()
  ): Future[(Workspace, SyncResult)] = {
    if (!FeatureToggles.atProtocolEnabled) {
      Future.successful((currentWorkspace, SyncResult.Disabled))
    } else {
      user match {
        case Some(u) =>
          onStatusChange(SyncStatus.Syncing)
          log.info(s"Syncing workspace from PDS for user ${u.did}...")

          PdsClient.loadWorkspace(u).map { remoteWorkspace =>
            if (remoteWorkspace != currentWorkspace) {
              log.info("PDS workspace differs from local, updating...")
              // Save remote to local cache
              workspaceService.save(remoteWorkspace)
              onStatusChange(SyncStatus.Synced)
              (remoteWorkspace, SyncResult.Success)
            } else {
              log.debug("PDS workspace matches local, no update needed")
              onStatusChange(SyncStatus.Synced)
              (currentWorkspace, SyncResult.NoChange)
            }
          }.recover {
            case e: Exception =>
              log.error(s"Failed to sync from PDS: ${e.getMessage}")
              onStatusChange(SyncStatus.Error)
              (currentWorkspace, SyncResult.Error(e.getMessage))
          }

        case None =>
          log.debug("AT Protocol enabled but no user logged in, using local workspace only")
          onStatusChange(SyncStatus.Offline)
          Future.successful((currentWorkspace, SyncResult.Offline))
      }
    }
  }

  /**
   * Attempts to sync workspace to PDS if AT Protocol is enabled and user is logged in.
   *
   * @param user Current logged-in user (if any)
   * @param workspace The workspace to sync
   * @param onStatusChange Callback for status updates
   * @return Future containing the sync result
   */
  def syncToPds(
    user: Option[User],
    workspace: Workspace,
    onStatusChange: SyncStatus => Unit = _ => ()
  ): Future[SyncResult] = {
    if (!FeatureToggles.atProtocolEnabled) {
      Future.successful(SyncResult.Disabled)
    } else {
      user match {
        case Some(u) =>
          onStatusChange(SyncStatus.Syncing)
          log.info(s"Syncing workspace to PDS for user ${u.did}...")

          PdsClient.saveWorkspace(u, workspace).map { _ =>
            log.info("Successfully synced workspace to PDS")
            onStatusChange(SyncStatus.Synced)
            SyncResult.Success
          }.recover {
            case e: Exception =>
              log.error(s"Failed to sync to PDS: ${e.getMessage}")
              onStatusChange(SyncStatus.Error)
              SyncResult.Error(e.getMessage)
          }

        case None =>
          onStatusChange(SyncStatus.Offline)
          Future.successful(SyncResult.Offline)
      }
    }
  }

  /**
   * Saves workspace locally and optionally syncs to PDS.
   *
   * @param workspace The workspace to save
   * @param user Current logged-in user (if any)
   * @param onStatusChange Callback for status updates
   * @return Either error message or sync result
   */
  def saveAndSync(
    workspace: Workspace,
    user: Option[User],
    onStatusChange: SyncStatus => Unit = _ => ()
  ): Future[Either[String, SyncResult]] = {
    // Step 1: Save to local JSON (synchronous)
    workspaceService.save(workspace) match {
      case Left(error) =>
        log.error(s"Error saving workspace locally: $error")
        Future.successful(Left(error))

      case Right(_) =>
        log.debug("Workspace saved to local cache")
        // Step 2: Sync to PDS in background if available
        syncToPds(user, workspace, onStatusChange).map(Right(_))
    }
  }

  /**
   * Checks if sync is available (AT Protocol enabled and user logged in).
   */
  def isSyncAvailable(user: Option[User]): Boolean = {
    FeatureToggles.atProtocolEnabled && user.isDefined
  }

  /**
   * Gets the appropriate sync status for the current state.
   */
  def getInitialSyncStatus(user: Option[User]): SyncStatus = {
    if (!FeatureToggles.atProtocolEnabled) {
      SyncStatus.Synced // Local-only mode, always "synced"
    } else if (user.isEmpty) {
      SyncStatus.Offline
    } else {
      SyncStatus.Synced
    }
  }
}
