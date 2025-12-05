package com.decodingus.workspace

import com.decodingus.auth.User
import com.decodingus.config.FeatureToggles
import com.decodingus.pds.PdsClient
import com.decodingus.workspace.model.{Workspace, Project, Biosample, WorkspaceContent, SyncStatus}
import scalafx.beans.property.{ObjectProperty, ReadOnlyObjectProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.application.Platform

import java.time.LocalDateTime
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.{Success, Failure}

class WorkbenchViewModel(val workspaceService: WorkspaceService) {

  // --- Sync State ---
  val syncStatus: ObjectProperty[SyncStatus] = ObjectProperty(SyncStatus.Synced)
  val lastSyncError: StringProperty = StringProperty("")

  // Current authenticated user (set from GenomeNavigatorApp when user logs in)
  val currentUser: ObjectProperty[Option[User]] = ObjectProperty(None)

  // --- Model State ---
  private val _workspace = ObjectProperty(
    Workspace(
      lexicon = 1,
      id = "com.decodingus.atmosphere.workspace",
      main = WorkspaceContent(
        samples = List.empty,
        projects = List.empty
      )
    )
  )
  // Exposed as ReadOnlyProperty for external observation, preventing direct modification
  val workspace: ReadOnlyObjectProperty[Workspace] = _workspace

  // --- UI Observable Collections ---
  // These buffers will be updated when _workspace changes, and views will bind to them
  val projects: ObservableBuffer[Project] = ObservableBuffer[Project]()
  val samples: ObservableBuffer[Biosample] = ObservableBuffer[Biosample]()

  // --- Selected Items ---
  // Properties for currently selected project/subject, allowing two-way binding with UI
  val selectedProject: ObjectProperty[Option[Project]] = ObjectProperty(None)
  val selectedSubject: ObjectProperty[Option[Biosample]] = ObjectProperty(None)

  // Listen to changes in the internal _workspace and update observable buffers
  // NOTE: This listener MUST be registered BEFORE loadWorkspace() is called,
  // otherwise the initial load won't trigger the buffer updates.
  _workspace.onChange { (_, oldWorkspace, newWorkspace) =>
    println(s"[DEBUG] WorkbenchViewModel: _workspace changed from (samples: ${oldWorkspace.main.samples.size}, projects: ${oldWorkspace.main.projects.size}) to (samples: ${newWorkspace.main.samples.size}, projects: ${newWorkspace.main.projects.size})")
    syncBuffers(newWorkspace)
  }

  // --- Initialization ---
  // Load workspace on creation of ViewModel (AFTER onChange listener is registered)
  loadWorkspace()

  /** Synchronizes the observable buffers with the workspace state */
  private def syncBuffers(workspace: Workspace): Unit = {
    println(s"[DEBUG] WorkbenchViewModel: Syncing buffers with workspace...")
    projects.clear()
    projects ++= workspace.main.projects
    samples.clear()
    samples ++= workspace.main.samples
    println(s"[DEBUG] WorkbenchViewModel: Buffers synced. Projects: ${projects.size}, Samples: ${samples.size}")
  }

  // --- Commands (Business Logic) ---

  /**
   * Loads workspace with the following strategy:
   * 1. Load from local JSON immediately (fast startup)
   * 2. If AT Protocol enabled and user logged in, fetch from PDS in background
   * 3. If PDS has newer/different data, update local state and cache
   */
  def loadWorkspace(): Unit = {
    // Step 1: Load from local cache immediately
    workspaceService.load().fold(
      error => {
        println(s"[ViewModel] Error loading local workspace: $error")
        _workspace.value = emptyWorkspace
      },
      loadedWorkspace => {
        println(s"[ViewModel] Loaded workspace from local cache: ${loadedWorkspace.main.samples.size} samples")
        _workspace.value = loadedWorkspace
      }
    )

    // Step 2: If AT Protocol enabled and user is logged in, sync from PDS
    syncFromPdsIfAvailable()
  }

  /**
   * Attempts to sync workspace from PDS if AT Protocol is enabled and user is logged in.
   * Updates local state and cache if remote data differs.
   */
  def syncFromPdsIfAvailable(): Unit = {
    if (FeatureToggles.atProtocolEnabled) {
      currentUser.value match {
        case Some(user) =>
          syncStatus.value = SyncStatus.Syncing
          println(s"[ViewModel] Syncing workspace from PDS for user ${user.did}...")

          PdsClient.loadWorkspace(user).onComplete {
            case Success(remoteWorkspace) =>
              Platform.runLater {
                // Compare and update if different
                if (remoteWorkspace != _workspace.value) {
                  println(s"[ViewModel] PDS workspace differs from local, updating...")
                  _workspace.value = remoteWorkspace
                  // Update local cache
                  workspaceService.save(remoteWorkspace)
                } else {
                  println(s"[ViewModel] PDS workspace matches local, no update needed")
                }
                syncStatus.value = SyncStatus.Synced
                lastSyncError.value = ""
              }
            case Failure(e) =>
              Platform.runLater {
                println(s"[ViewModel] Failed to sync from PDS: ${e.getMessage}")
                syncStatus.value = SyncStatus.Error
                lastSyncError.value = e.getMessage
                // Continue with local data - app remains functional
              }
          }
        case None =>
          println(s"[ViewModel] AT Protocol enabled but no user logged in, using local workspace only")
          syncStatus.value = SyncStatus.Offline
      }
    }
  }

  /**
   * Saves workspace with optimistic update strategy:
   * 1. Update local state immediately (already done by caller)
   * 2. Save to local JSON (fast, ensures durability)
   * 3. Sync to PDS in background if AT Protocol enabled
   */
  def saveWorkspace(): Unit = {
    // Step 1: Save to local JSON (fast, synchronous)
    workspaceService.save(_workspace.value).fold(
      error => {
        println(s"[ViewModel] Error saving workspace locally: $error")
      },
      _ => {
        println(s"[ViewModel] Workspace saved to local cache")
      }
    )

    // Step 2: Sync to PDS in background if available
    syncToPdsIfAvailable()
  }

  /**
   * Attempts to sync workspace to PDS if AT Protocol is enabled and user is logged in.
   */
  private def syncToPdsIfAvailable(): Unit = {
    if (FeatureToggles.atProtocolEnabled) {
      currentUser.value match {
        case Some(user) =>
          syncStatus.value = SyncStatus.Syncing
          println(s"[ViewModel] Syncing workspace to PDS for user ${user.did}...")

          PdsClient.saveWorkspace(user, _workspace.value).onComplete {
            case Success(_) =>
              Platform.runLater {
                println(s"[ViewModel] Successfully synced workspace to PDS")
                syncStatus.value = SyncStatus.Synced
                lastSyncError.value = ""
              }
            case Failure(e) =>
              Platform.runLater {
                println(s"[ViewModel] Failed to sync to PDS: ${e.getMessage}")
                syncStatus.value = SyncStatus.Error
                lastSyncError.value = e.getMessage
                // Local save succeeded, so data is not lost
              }
          }
        case None =>
          // No user logged in, local-only mode
          syncStatus.value = SyncStatus.Offline
      }
    }
  }

  private def emptyWorkspace: Workspace = Workspace(
    lexicon = 1,
    id = "com.decodingus.atmosphere.workspace",
    main = WorkspaceContent(samples = List.empty, projects = List.empty)
  )

  def addSubject(newBiosample: Biosample): Unit = {
    // Update the internal workspace object
    val updatedSamples = _workspace.value.main.samples :+ newBiosample
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))
    
    // Select the newly added subject (this will trigger UI binding)
    selectedSubject.value = Some(newBiosample)
    
    // Save changes
    saveWorkspace()
  }

  def addProject(newProject: Project): Unit = {
    val updatedProjects = _workspace.value.main.projects :+ newProject
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))
    
    selectedProject.value = Some(newProject)
    
    saveWorkspace()
  }

  // --- Other Business Logic (e.g., Analysis) ---
  // These will be fleshed out later and likely operate on selectedSubject/selectedProject
  def performInitialAnalysis(filePath: String): Unit = {
    println(s"Performing initial analysis on $filePath (logic to be implemented in ViewModel).")
    // This will involve calling analysis processors and updating the selectedSubject's sequenceData
  }

  def performDeepCoverageAnalysis(biosample: Biosample): Unit = {
    println(s"Performing deep coverage analysis for ${biosample.donorIdentifier} (logic to be implemented in ViewModel).")
  }
}
