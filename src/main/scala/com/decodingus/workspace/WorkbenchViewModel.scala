package com.decodingus.workspace

import com.decodingus.analysis.*
import com.decodingus.auth.User
import com.decodingus.client.DecodingUsClient
import com.decodingus.config.{FeatureToggles, UserPreferencesService}
import com.decodingus.db.Transactor
import com.decodingus.haplogroup.model.HaplogroupResult as AnalysisHaplogroupResult
import com.decodingus.haplogroup.processor.HaplogroupProcessor
import com.decodingus.haplogroup.tree.{TreeProviderType, TreeType}
import com.decodingus.model.{LibraryStats, WgsMetrics}
import com.decodingus.refgenome.config.ReferenceConfigService
import com.decodingus.refgenome.{ReferenceGateway, ReferenceResolveResult, YRegionAnnotator}
import com.decodingus.repository.BiosampleRepository
import com.decodingus.service.{DatabaseContext, H2WorkspaceService, SequenceDataManager}
import com.decodingus.util.Logger
import com.decodingus.workspace.model.{Alignment, AlignmentMetrics, Biosample, CallMethod, ChipProfile, ContigMetrics, DnaType, FileInfo, HaplogroupAssignments, HaplogroupTechnology, Project, RecordMeta, RunHaplogroupCall, SequenceRun, StrProfile, SyncStatus, Workspace, WorkspaceContent, HaplogroupResult as WorkspaceHaplogroupResult}
import com.decodingus.workspace.services.*
import com.decodingus.yprofile.model.*
import com.decodingus.yprofile.repository.*
import com.decodingus.yprofile.service.YProfileService
import com.decodingus.yprofile.repository.YSnpPanelRepository
import htsjdk.samtools.SamReaderFactory
import scalafx.application.Platform
import scalafx.beans.property.*
import scalafx.collections.ObservableBuffer

import java.io.File
import java.nio.file.{Files, Path}
import java.time.LocalDateTime
import java.util.UUID
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

class WorkbenchViewModel(
                          databaseContext: DatabaseContext
                        ) {
  private val log = Logger[WorkbenchViewModel]

  // H2 service for atomic CRUD operations
  private val h2Service: H2WorkspaceService = databaseContext.workspaceService

  // Legacy adapter for load/sync operations
  private val workspaceService: WorkspaceService = H2WorkspaceAdapter(databaseContext)

  // --- Service Instances ---
  private val workspaceOps = new WorkspaceOperations()
  private val syncService = new SyncService(workspaceService)
  private val fingerprintMatchService = new FingerprintMatchService()

  // Centralized manager for SequenceRun and Alignment CRUD operations
  private val sequenceDataManager: SequenceDataManager = databaseContext.sequenceDataManager

  // Y Profile service
  private val yProfileService: Option[YProfileService] = Some(databaseContext.transactor).map { tx =>
    YProfileService(
      tx,
      YChromosomeProfileRepository(),
      YProfileSourceRepository(),
      YProfileRegionRepository(),
      YProfileVariantRepository(),
      YVariantSourceCallRepository(),
      YVariantAuditRepository(),
      YSourceCallAlignmentRepository()
    )
  }

  // Analysis coordinator (uses YProfileService for auto-populating Y profiles during analysis)
  // Receives h2Service for persisting analysis results after each step
  private val analysisCoordinator = new AnalysisCoordinator(h2Service, yProfileService)

  // --- Sync State ---
  val syncStatus: ObjectProperty[SyncStatus] = ObjectProperty(SyncStatus.Synced)
  val lastSyncError: StringProperty = StringProperty("")

  // Async sync notifier for status bar binding
  // Provides observable properties for pending sync count, conflicts, online status
  val syncNotifier: com.decodingus.sync.ConflictNotifier = com.decodingus.sync.ConflictNotifier()

  // Current authenticated user (set from GenomeNavigatorApp when user logs in)
  val currentUser: ObjectProperty[Option[User]] = ObjectProperty(None)

  // When user logs in, backfill atUri for any samples/projects created while logged out
  currentUser.onChange { (_, _, newUser) =>
    newUser.foreach(user => backfillAtUris(user.did))
  }

  // --- Model State ---
  private val _workspace = ObjectProperty(Workspace.empty)
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

  // --- Filtering ---
  // Filter text for projects and subjects lists
  val projectFilter: StringProperty = StringProperty("")
  val subjectFilter: StringProperty = StringProperty("")

  // Filtered observable collections (updated when filter or source data changes)
  val filteredProjects: ObservableBuffer[Project] = ObservableBuffer[Project]()
  val filteredSamples: ObservableBuffer[Biosample] = ObservableBuffer[Biosample]()

  // Apply filters when filter text changes
  projectFilter.onChange { (_, _, _) => applyFilters() }
  subjectFilter.onChange { (_, _, _) => applyFilters() }

  /** Applies current filter text to projects and samples */
  private def applyFilters(): Unit = {
    val projectQuery = projectFilter.value.toLowerCase.trim
    val subjectQuery = subjectFilter.value.toLowerCase.trim

    // Filter projects
    filteredProjects.clear()
    if (projectQuery.isEmpty) {
      filteredProjects ++= projects
    } else {
      filteredProjects ++= projects.filter { p =>
        p.projectName.toLowerCase.contains(projectQuery) ||
          p.description.exists(_.toLowerCase.contains(projectQuery)) ||
          p.administrator.toLowerCase.contains(projectQuery)
      }
    }

    // Filter subjects
    filteredSamples.clear()
    if (subjectQuery.isEmpty) {
      filteredSamples ++= samples
    } else {
      filteredSamples ++= samples.filter { s =>
        s.donorIdentifier.toLowerCase.contains(subjectQuery) ||
          s.sampleAccession.toLowerCase.contains(subjectQuery) ||
          s.description.exists(_.toLowerCase.contains(subjectQuery)) ||
          s.centerName.exists(_.toLowerCase.contains(subjectQuery))
      }
    }
  }

  // Listen to changes in the internal _workspace and update observable buffers
  // NOTE: This listener MUST be registered BEFORE loadWorkspace() is called,
  // otherwise the initial load won't trigger the buffer updates.
  _workspace.onChange { (_, oldWorkspace, newWorkspace) =>
    log.debug(s" WorkbenchViewModel: _workspace changed from (samples: ${oldWorkspace.main.samples.size}, projects: ${oldWorkspace.main.projects.size}) to (samples: ${newWorkspace.main.samples.size}, projects: ${newWorkspace.main.projects.size})")
    syncBuffers(newWorkspace)
  }

  // --- Initialization ---
  // Set up the SequenceDataManager callbacks to keep in-memory state in sync
  sequenceDataManager.setWorkspaceCallbacks(
    getter = () => _workspace.value,
    updater = (updated: Workspace) => _workspace.value = updated
  )

  // Load workspace on creation of ViewModel (AFTER onChange listener is registered)
  loadWorkspace()

  /** Synchronizes the observable buffers with the workspace state */
  private def syncBuffers(workspace: Workspace): Unit = {
    log.debug(s" WorkbenchViewModel.syncBuffers: Syncing buffers with workspace...")
    log.debug(s"   workspace.main.samples: ${workspace.main.samples.size}")
    workspace.main.samples.foreach { s =>
      log.debug(s"     Sample ${s.sampleAccession}: sequenceRunRefs=${s.sequenceRunRefs.size} ${s.sequenceRunRefs.mkString(", ")}")
    }
    log.debug(s"   workspace.main.sequenceRuns: ${workspace.main.sequenceRuns.size}")
    workspace.main.sequenceRuns.foreach { sr =>
      log.debug(s"     SequenceRun atUri=${sr.atUri}, biosampleRef=${sr.biosampleRef}, alignmentRefs=${sr.alignmentRefs.size}")
    }
    log.debug(s"   workspace.main.alignments: ${workspace.main.alignments.size}")

    // Preserve current selection identifiers before updating buffers
    val selectedProjectName = selectedProject.value.map(_.projectName)
    val selectedSampleAccession = selectedSubject.value.map(_.sampleAccession)

    // Use setAll for atomic update (single onChange event instead of clear + add)
    import scala.jdk.CollectionConverters.*
    projects.delegate.setAll(workspace.main.projects.asJava)
    samples.delegate.setAll(workspace.main.samples.asJava)
    // Also refresh filtered lists
    applyFilters()

    // Restore selection by finding the updated objects with the same identifiers
    selectedProjectName.foreach { name =>
      workspace.main.projects.find(_.projectName == name).foreach { project =>
        selectedProject.value = Some(project)
      }
    }
    selectedSampleAccession.foreach { accession =>
      workspace.main.samples.find(_.sampleAccession == accession).foreach { sample =>
        selectedSubject.value = Some(sample)
      }
    }

    log.debug(s" WorkbenchViewModel.syncBuffers: Buffers synced. Projects: ${projects.size}, Samples: ${samples.size}")
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
        log.error(s"Error loading local workspace: $error")
        _workspace.value = emptyWorkspace
      },
      loadedWorkspace => {
        log.info(s" Loaded workspace from local cache: ${loadedWorkspace.main.samples.size} samples")
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
    syncService.syncFromPds(
      currentUser.value,
      _workspace.value,
      status => Platform.runLater {
        syncStatus.value = status
      }
    ).foreach { case (workspace, result) =>
      Platform.runLater {
        result match {
          case SyncResult.Success =>
            _workspace.value = workspace
            lastSyncError.value = ""
          case SyncResult.NoChange =>
            lastSyncError.value = ""
          case SyncResult.Error(msg) =>
            lastSyncError.value = msg
          case _ =>
            lastSyncError.value = ""
        }
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
    syncService.saveAndSync(
      _workspace.value,
      currentUser.value,
      status => Platform.runLater {
        syncStatus.value = status
      }
    ).foreach {
      case Left(error) =>
        Platform.runLater {
          lastSyncError.value = error
        }
      case Right(SyncResult.Error(msg)) =>
        Platform.runLater {
          lastSyncError.value = msg
        }
      case Right(_) =>
        Platform.runLater {
          lastSyncError.value = ""
        }
    }
  }

  private def emptyWorkspace: Workspace = Workspace.empty

  /** Gets the current workspace state for service calls */
  private def currentState: WorkspaceState = WorkspaceState(_workspace.value)

  /** Applies a new workspace state and saves */

  /**
   * Updates in-memory state from a WorkspaceState.
   * Used by atomic operations that have already persisted to H2.
   * Does NOT trigger a save (data is already persisted).
   */
  private def updateInMemoryState(newState: WorkspaceState): Unit = {
    _workspace.value = newState.workspace
  }

  /**
   * Updates in-memory state from a WorkspaceState.
   * Used by batch operations (like AnalysisCoordinator) that handle their own
   * H2 persistence internally.
   *
   * Note: This no longer calls saveWorkspace() because all CRUD operations now
   * persist atomically to H2. This method only updates the in-memory state.
   */
  private def applyState(newState: WorkspaceState): Unit = {
    _workspace.value = newState.workspace
  }

  // --- Subject CRUD Operations ---

  /**
   * Creates a new subject and adds it to the workspace.
   * ATOMIC: Persists to H2 first, then updates in-memory state.
   */
  def addSubject(newBiosample: Biosample): Unit = {
    val userDid = currentUser.value.map(_.did)
    val (newState, enrichedBiosample) = workspaceOps.addSubject(currentState, newBiosample, userDid)

    // Persist to H2 atomically
    h2Service.createBiosample(enrichedBiosample) match {
      case Right(persistedBiosample) =>
        log.info(s" Biosample persisted to H2: ${persistedBiosample.sampleAccession}")
        // IMPORTANT: Update state with persisted biosample (which has the correct DB-assigned atUri/UUID)
        // The enrichedBiosample may have a different atUri that doesn't contain a valid UUID
        val stateWithPersistedBiosample = newState.copy(
          workspace = newState.workspace.copy(
            main = newState.workspace.main.copy(
              samples = newState.workspace.main.samples.map { s =>
                if s.sampleAccession == persistedBiosample.sampleAccession then persistedBiosample else s
              }
            )
          )
        )
        updateInMemoryState(stateWithPersistedBiosample)
        selectedSubject.value = Some(persistedBiosample)
      case Left(error) =>
        log.error(s"Failed to persist Biosample to H2: $error")
        // Fallback: still update in-memory but note it's not persisted
        updateInMemoryState(newState)
        selectedSubject.value = Some(enrichedBiosample)
    }
  }

  /**
   * Updates an existing subject identified by sampleAccession.
   * ATOMIC: Persists to H2 first, then updates in-memory state.
   */
  def updateSubject(updatedBiosample: Biosample): Unit = {
    val newState = workspaceOps.updateSubject(currentState, updatedBiosample)
    val biosampleWithMeta = newState.workspace.main.samples.find(_.sampleAccession == updatedBiosample.sampleAccession)

    biosampleWithMeta.foreach { bs =>
      h2Service.updateBiosample(bs) match {
        case Right(persisted) =>
          log.info(s" Biosample updated in H2: ${persisted.sampleAccession}")
        case Left(error) =>
          log.error(s"Failed to update Biosample in H2: $error")
      }
    }
    updateInMemoryState(newState)
    selectedSubject.value = biosampleWithMeta
  }

  /**
   * Internal: Updates a subject without modifying meta (used when meta is already updated).
   * ATOMIC: Persists to H2 first, then updates in-memory state.
   */
  private def updateSubjectDirect(updatedBiosample: Biosample): Unit = {
    val newState = workspaceOps.updateSubjectDirect(currentState, updatedBiosample)

    h2Service.updateBiosample(updatedBiosample) match {
      case Right(persisted) =>
        log.info(s" Biosample updated (direct) in H2: ${persisted.sampleAccession}")
      case Left(error) =>
        log.error(s"Failed to update Biosample (direct) in H2: $error")
    }
    updateInMemoryState(newState)
    selectedSubject.value = Some(updatedBiosample)
  }

  /**
   * Deletes a subject identified by sampleAccession.
   * ATOMIC: Deletes from H2 first, then updates in-memory state.
   */
  def deleteSubject(sampleAccession: String): Unit = {
    import com.decodingus.service.EntityConversions.parseIdFromRef

    // Find the biosample to get its ID
    findSubject(sampleAccession) match {
      case None =>
        log.info(s" Cannot delete - biosample not found: $sampleAccession")
        return
      case Some(biosample) =>
        // Extract UUID from atUri
        val biosampleId = biosample.atUri.flatMap(parseIdFromRef)
        biosampleId match {
          case Some(id) =>
            // Delete from H2 first
            h2Service.deleteBiosample(id) match {
              case Right(deleted) if deleted =>
                log.info(s" Biosample deleted from H2: $sampleAccession")
                // Update in-memory state
                val newState = workspaceOps.deleteSubject(currentState, sampleAccession)
                updateInMemoryState(newState)
              case Right(_) =>
                log.info(s" Biosample not found in H2, removing from in-memory: $sampleAccession")
                val newState = workspaceOps.deleteSubject(currentState, sampleAccession)
                updateInMemoryState(newState)
              case Left(error) =>
                log.error(s"Failed to delete Biosample from H2: $error")
                // Still remove from in-memory to keep UI consistent
                val newState = workspaceOps.deleteSubject(currentState, sampleAccession)
                updateInMemoryState(newState)
            }
          case None =>
            log.info(s" Cannot extract ID from biosample atUri, removing from in-memory only")
            val newState = workspaceOps.deleteSubject(currentState, sampleAccession)
            updateInMemoryState(newState)
        }
    }

    // Clear selection if the deleted subject was selected
    selectedSubject.value match {
      case Some(selected) if selected.sampleAccession == sampleAccession =>
        selectedSubject.value = None
      case _ => // Keep current selection
    }
  }

  /** Finds a subject by sampleAccession */
  def findSubject(sampleAccession: String): Option[Biosample] = {
    _workspace.value.main.samples.find(_.sampleAccession == sampleAccession)
  }

  // --- Project CRUD Operations ---

  /**
   * Creates a new project and adds it to the workspace.
   * ATOMIC: Persists to H2 first, then updates in-memory state.
   */
  def addProject(newProject: Project): Unit = {
    val userDid = currentUser.value.map(_.did)
    val (newState, enrichedProject) = workspaceOps.addProject(currentState, newProject, userDid)

    // Persist to H2 atomically
    h2Service.createProject(enrichedProject) match {
      case Right(persistedProject) =>
        log.info(s" Project persisted to H2: ${persistedProject.projectName}")
        updateInMemoryState(newState)
        selectedProject.value = Some(persistedProject)
      case Left(error) =>
        log.error(s"Failed to persist Project to H2: $error")
        // Fallback: still update in-memory but note it's not persisted
        updateInMemoryState(newState)
        selectedProject.value = Some(enrichedProject)
    }
  }

  /**
   * Updates an existing project identified by projectName.
   * ATOMIC: Persists to H2 first, then updates in-memory state.
   */
  def updateProject(updatedProject: Project): Unit = {
    val newState = workspaceOps.updateProject(currentState, updatedProject)
    val projectWithMeta = newState.workspace.main.projects.find(_.projectName == updatedProject.projectName)

    projectWithMeta.foreach { proj =>
      h2Service.updateProject(proj) match {
        case Right(persisted) =>
          log.info(s" Project updated in H2: ${persisted.projectName}")
        case Left(error) =>
          log.error(s"Failed to update Project in H2: $error")
      }
    }
    updateInMemoryState(newState)
    selectedProject.value = projectWithMeta
  }

  /**
   * Deletes a project by projectName.
   * ATOMIC: Deletes from H2 first, then updates in-memory state.
   */
  def deleteProject(projectName: String): Unit = {
    import com.decodingus.service.EntityConversions.parseIdFromRef

    // Find the project to get its ID
    findProject(projectName) match {
      case None =>
        log.info(s" Cannot delete - project not found: $projectName")
        return
      case Some(project) =>
        // Extract UUID from atUri
        val projectId = project.atUri.flatMap(parseIdFromRef)
        projectId match {
          case Some(id) =>
            // Delete from H2 first
            h2Service.deleteProject(id) match {
              case Right(deleted) if deleted =>
                log.info(s" Project deleted from H2: $projectName")
                val newState = workspaceOps.deleteProject(currentState, projectName)
                updateInMemoryState(newState)
              case Right(_) =>
                log.info(s" Project not found in H2, removing from in-memory: $projectName")
                val newState = workspaceOps.deleteProject(currentState, projectName)
                updateInMemoryState(newState)
              case Left(error) =>
                log.error(s"Failed to delete Project from H2: $error")
                val newState = workspaceOps.deleteProject(currentState, projectName)
                updateInMemoryState(newState)
            }
          case None =>
            log.info(s" Cannot extract ID from project atUri, removing from in-memory only")
            val newState = workspaceOps.deleteProject(currentState, projectName)
            updateInMemoryState(newState)
        }
    }

    selectedProject.value match {
      case Some(selected) if selected.projectName == projectName =>
        selectedProject.value = None
      case _ =>
    }
  }

  /** Finds a project by projectName */
  def findProject(projectName: String): Option[Project] = {
    _workspace.value.main.projects.find(_.projectName == projectName)
  }

  /**
   * Backfills atUri for any samples/projects that were created while logged out.
   * ATOMIC: Persists each updated entity to H2, then updates in-memory state.
   */
  private def backfillAtUris(did: String): Unit = {
    val newState = workspaceOps.backfillAtUris(currentState, did)
    if (newState.workspace != currentState.workspace) {
      // Persist updated biosamples to H2
      newState.workspace.main.samples.foreach { sample =>
        h2Service.updateBiosample(sample) match {
          case Right(_) => // Success
          case Left(error) =>
            log.error(s"Failed to backfill biosample atUri in H2: $error")
        }
      }
      // Persist updated projects to H2
      newState.workspace.main.projects.foreach { project =>
        h2Service.updateProject(project) match {
          case Right(_) => // Success
          case Left(error) =>
            log.error(s"Failed to backfill project atUri in H2: $error")
        }
      }
      updateInMemoryState(newState)
      log.info(s" Backfilled atUri for samples/projects after login")
    }
  }

  /**
   * Adds a subject (by accession) to a project's members list.
   * ATOMIC: Persists updated project to H2 first, then updates in-memory state.
   */
  def addSubjectToProject(projectName: String, sampleAccession: String): Boolean = {
    workspaceOps.addSubjectToProject(currentState, projectName, sampleAccession) match {
      case Right(newState) =>
        // Find the updated project and persist to H2
        val updatedProject = newState.workspace.main.projects.find(_.projectName == projectName)
        updatedProject.foreach { proj =>
          h2Service.updateProject(proj) match {
            case Right(persisted) =>
              log.info(s" Project membership updated in H2: ${persisted.projectName}")
            case Left(error) =>
              log.error(s"Failed to update project membership in H2: $error")
          }
        }
        updateInMemoryState(newState)
        selectedProject.value = updatedProject
        true
      case Left(error) =>
        log.info(s" $error")
        false
    }
  }

  /**
   * Removes a subject (by accession) from a project's members list.
   * ATOMIC: Persists updated project to H2 first, then updates in-memory state.
   */
  def removeSubjectFromProject(projectName: String, sampleAccession: String): Boolean = {
    workspaceOps.removeSubjectFromProject(currentState, projectName, sampleAccession) match {
      case Right(newState) =>
        // Find the updated project and persist to H2
        val updatedProject = newState.workspace.main.projects.find(_.projectName == projectName)
        updatedProject.foreach { proj =>
          h2Service.updateProject(proj) match {
            case Right(persisted) =>
              log.info(s" Project membership updated in H2: ${persisted.projectName}")
            case Left(error) =>
              log.error(s"Failed to update project membership in H2: $error")
          }
        }
        updateInMemoryState(newState)
        selectedProject.value = updatedProject
        true
      case Left(error) =>
        log.info(s" $error")
        false
    }
  }

  /**
   * Internal: Updates a project without modifying meta (used when meta is already updated).
   * ATOMIC: Persists to H2 first, then updates in-memory state.
   */
  private def updateProjectDirect(updatedProject: Project): Unit = {
    val newState = workspaceOps.updateProjectDirect(currentState, updatedProject)

    h2Service.updateProject(updatedProject) match {
      case Right(persisted) =>
        log.info(s" Project updated (direct) in H2: ${persisted.projectName}")
      case Left(error) =>
        log.error(s"Failed to update Project (direct) in H2: $error")
    }
    updateInMemoryState(newState)
    selectedProject.value = Some(updatedProject)
  }

  /** Gets subjects that are members of a project */
  def getProjectMembers(projectName: String): List[Biosample] = {
    findProject(projectName) match {
      case Some(project) =>
        project.memberRefs.flatMap(accession => findSubject(accession))
      case None =>
        List.empty
    }
  }

  /** Gets subjects that are NOT members of a project (for adding) */
  def getNonProjectMembers(projectName: String): List[Biosample] = {
    findProject(projectName) match {
      case Some(project) =>
        _workspace.value.main.samples.filterNot(s => project.memberRefs.contains(s.sampleAccession))
      case None =>
        _workspace.value.main.samples
    }
  }

  // --- SequenceRun CRUD Operations (first-class records in workspace) ---

  /**
   * Creates a new SequenceRun from a FileInfo.
   * All metadata (platform, reads, etc.) will be populated during analysis.
   * Returns the index of the new entry, or -1 if duplicate/error.
   *
   * Uses SequenceDataManager for atomic H2 persistence + in-memory state update.
   */
  def addSequenceRunFromFile(sampleAccession: String, fileInfo: FileInfo): Int = {
    // Create initial SequenceRun with placeholder values
    val initialRun = SequenceRun(
      atUri = None, // Will be assigned by SequenceDataManager
      meta = RecordMeta.initial,
      biosampleRef = "", // Will be set by SequenceDataManager
      platformName = "Unknown",
      instrumentModel = None,
      testType = "Unknown",
      libraryLayout = None,
      totalReads = None,
      readLength = None,
      meanInsertSize = None,
      files = List.empty, // fileInfo added by manager
      alignmentRefs = List.empty
    )

    // SequenceDataManager handles: validation, persistence, in-memory state update
    sequenceDataManager.createSequenceRun(sampleAccession, initialRun, fileInfo) match {
      case Right(result) =>
        log.info(s"SequenceRun created via SequenceDataManager: ${result.sequenceRun.atUri}")
        result.index
      case Left(error) =>
        log.error(s"Failed to create SequenceRun: $error")
        -1
    }
  }

  // Backward compatibility alias
  def addSequenceDataFromFile(sampleAccession: String, fileInfo: FileInfo): Int =
    addSequenceRunFromFile(sampleAccession, fileInfo)

  /**
   * Result of adding a file - either a new run or added to existing.
   */
  sealed trait AddFileResult

  case class NewRunCreated(index: Int) extends AddFileResult

  case class AddedToExistingRun(index: Int, referenceBuild: String) extends AddFileResult

  /**
   * Adds a file and immediately runs library stats analysis.
   * This is the primary flow for adding new sequencing data.
   *
   * Pipeline:
   * 1. Check for exact duplicate by checksum
   * 2. Analyze file to get library stats + fingerprint
   * 3. Check for fingerprint match (same run, different reference)
   * 4. If match: add alignment to existing SequenceRun
   * 5. If no match: create new SequenceRun + Alignment
   *
   * @param sampleAccession The subject's accession ID
   * @param fileInfo        The file information
   * @param onProgress      Progress callback (message, percent)
   * @param onComplete      Completion callback with result
   */
  def addFileAndAnalyze(
                         sampleAccession: String,
                         fileInfo: FileInfo,
                         onProgress: (String, Double) => Unit,
                         onComplete: Either[String, (Int, LibraryStats)] => Unit
                       ): Unit = {
    // Step 1: Check for exact duplicate by checksum (including alignment files)
    onProgress("Checking for duplicates...", 0.05)

    findSubject(sampleAccession) match {
      case None =>
        onComplete(Left(s"Subject $sampleAccession not found"))
        return
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        val seqRunChecksums = sequenceRuns.flatMap(_.files.flatMap(_.checksum))
        // Also include alignment file checksums
        val alignmentChecksums = sequenceRuns.flatMap { sr =>
          _workspace.value.main.getAlignmentsForSequenceRun(sr).flatMap(_.files.flatMap(_.checksum))
        }
        val existingChecksums = (seqRunChecksums ++ alignmentChecksums).toSet
        if (fileInfo.checksum.exists(existingChecksums.contains)) {
          val checksumPreview = fileInfo.checksum.map(_.take(12) + "...").getOrElse("unknown")
          log.warn(s"Duplicate file detected: ${fileInfo.fileName} (checksum: $checksumPreview)")
          onComplete(Left(s"Duplicate file - ${fileInfo.fileName} has already been added (checksum: $checksumPreview)"))
          return
        }
    }

    val bamPath = fileInfo.location.getOrElse("")
    log.info(s" Starting add+analyze pipeline for ${fileInfo.fileName}")

    analysisInProgress.value = true
    analysisError.value = ""
    analysisProgress.value = "Initializing analysis..."
    analysisProgressPercent.value = 0.1

    // Run analysis in background
    Future {
      try {
        // Step 2: Detect reference build from header
        onProgress("Reading BAM/CRAM header...", 0.15)
        updateProgress("Reading BAM/CRAM header...", 0.15)

        val header = SamReaderFactory.makeDefault().open(new File(bamPath)).getFileHeader
        val libraryStatsProcessor = new LibraryStatsProcessor()
        val referenceBuild = libraryStatsProcessor.detectReferenceBuild(header)

        if (referenceBuild == "Unknown") {
          throw new IllegalStateException("Could not determine reference build from BAM/CRAM header.")
        }

        // Step 3: Resolve reference genome
        onProgress(s"Resolving reference: $referenceBuild", 0.25)
        updateProgress(s"Resolving reference: $referenceBuild", 0.25)

        val referenceGateway = new ReferenceGateway((done, total) => {
          val pct = if (total > 0) 0.25 + (done.toDouble / total) * 0.25 else 0.25
          onProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct)
          updateProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct)
        })

        val referencePath = referenceGateway.resolve(referenceBuild) match {
          case Right(path) => path.toString
          case Left(error) => throw new Exception(s"Failed to resolve reference: $error")
        }

        // Step 4: Collect library stats
        onProgress("Analyzing library statistics...", 0.55)
        updateProgress("Analyzing library statistics...", 0.55)

        val libraryStats = libraryStatsProcessor.process(bamPath, referencePath, (message, current, total) => {
          val pct = 0.55 + (current.toDouble / total) * 0.35
          onProgress(s"Library Stats: $message", pct)
          updateProgress(s"Library Stats: $message", pct)
        })

        // Step 5: Check fingerprint and create/update SequenceRun + Alignment
        onProgress("Checking for matching sequence runs...", 0.92)
        updateProgress("Checking for matching sequence runs...", 0.92)

        val fingerprint = libraryStats.computeRunFingerprint
        val biosampleRef = findSubject(sampleAccession)
          .flatMap(_.atUri)
          .getOrElse(s"local:biosample:$sampleAccession")

        Platform.runLater {
          findSubject(sampleAccession) match {
            case Some(subject) =>
              // Check for fingerprint match (same run, different reference)
              val matchResult = findMatchingSequenceRun(biosampleRef, fingerprint, libraryStats)

              val runResult: Option[(Int, SequenceRun, Boolean)] = matchResult match {
                case FingerprintMatchResult.MatchFound(existingRun, idx, confidence) =>
                  // For LOW confidence matches, ask user to confirm
                  if (confidence == "LOW") {
                    import com.decodingus.ui.components.{FingerprintMatchDecision, FingerprintMatchDialog, GroupTogether, KeepSeparate}
                    val dialog = new FingerprintMatchDialog(
                      existingRun = existingRun,
                      newReferenceBuild = libraryStats.referenceBuild,
                      matchConfidence = confidence,
                      totalReads = libraryStats.readCount.toLong,
                      sampleName = libraryStats.sampleName,
                      libraryId = libraryStats.libraryId
                    )
                    val dialogResult = dialog.showAndWait()
                    val decision: Option[FingerprintMatchDecision] = dialogResult.asInstanceOf[Option[Option[FingerprintMatchDecision]]].flatten
                    decision match {
                      case Some(GroupTogether) =>
                        log.info(s" User confirmed LOW confidence match - adding ${libraryStats.referenceBuild} alignment to existing run")
                        Some((idx, existingRun, false))
                      case Some(KeepSeparate) =>
                        log.info(s" User chose to keep separate - creating new sequence run")
                        val newIndex = addSequenceRunFromFile(sampleAccession, fileInfo)
                        if (newIndex < 0) {
                          analysisInProgress.value = false
                          onComplete(Left("Failed to create sequence run entry"))
                          None
                        } else {
                          // Re-fetch updated subject to get correct sequence run refs
                          val updatedSubject = findSubject(sampleAccession).getOrElse(subject)
                          val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(updatedSubject)
                          Some((newIndex, sequenceRuns(newIndex), true))
                        }
                      case None =>
                        // User cancelled
                        analysisInProgress.value = false
                        onComplete(Left("User cancelled fingerprint match decision"))
                        None
                    }
                  } else {
                    // HIGH/MEDIUM confidence - auto-group
                    log.info(s" Fingerprint match found ($confidence confidence) - adding ${libraryStats.referenceBuild} alignment to existing run")
                    Some((idx, existingRun, false))
                  }

                case FingerprintMatchResult.NoMatch =>
                  // Create new SequenceRun
                  val newIndex = addSequenceRunFromFile(sampleAccession, fileInfo)
                  if (newIndex < 0) {
                    analysisInProgress.value = false
                    onComplete(Left("Failed to create sequence run entry"))
                    None
                  } else {
                    // Re-fetch updated subject to get correct sequence run refs
                    val updatedSubject = findSubject(sampleAccession).getOrElse(subject)
                    val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(updatedSubject)
                    Some((newIndex, sequenceRuns(newIndex), true))
                  }
              }

              runResult.foreach { case (resultIndex, seqRun, isNewRun) =>

                // Create alignment URI
                val alignUri = s"local:alignment:${subject.sampleAccession}:${libraryStats.referenceBuild}:${java.util.UUID.randomUUID().toString.take(8)}"

                val newAlignment = Alignment(
                  atUri = Some(alignUri),
                  meta = RecordMeta.initial,
                  sequenceRunRef = seqRun.atUri.getOrElse(""),
                  biosampleRef = Some(subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")),
                  referenceBuild = libraryStats.referenceBuild,
                  aligner = libraryStats.aligner,
                  files = List(fileInfo),
                  metrics = None
                )

                // Update the SequenceRun with inferred metadata
                val instrumentId = if (libraryStats.mostFrequentInstrumentId != "Unknown")
                  Some(libraryStats.mostFrequentInstrumentId)
                else
                  None

                // Only update sampleName if not already set (preserve manual edits)
                val sampleNameFromBam = seqRun.sampleName.orElse {
                  if (libraryStats.sampleName != "Unknown") Some(libraryStats.sampleName) else None
                }

                // Extract fingerprint fields (preserve manual edits)
                val libraryIdFromBam = seqRun.libraryId.orElse {
                  if (libraryStats.libraryId != "Unknown") Some(libraryStats.libraryId) else None
                }
                val platformUnitFromBam = seqRun.platformUnit.orElse(libraryStats.platformUnit)
                val runFingerprint = seqRun.runFingerprint.orElse(Some(fingerprint))

                // Add file to sequence run if not already present (for matched runs)
                val updatedFiles = if (seqRun.files.exists(_.checksum == fileInfo.checksum)) {
                  seqRun.files
                } else {
                  seqRun.files :+ fileInfo
                }

                val updatedSeqRun = seqRun.copy(
                  meta = seqRun.meta.updated("analysis"),
                  platformName = libraryStats.inferredPlatform,
                  instrumentModel = Some(libraryStats.mostFrequentInstrument),
                  instrumentId = instrumentId,
                  sampleName = sampleNameFromBam,
                  libraryId = libraryIdFromBam,
                  platformUnit = platformUnitFromBam,
                  runFingerprint = runFingerprint,
                  testType = inferTestType(libraryStats),
                  libraryLayout = Some(if (libraryStats.pairedReads > libraryStats.readCount / 2) "Paired-End" else "Single-End"),
                  totalReads = Some(libraryStats.readCount.toLong),
                  readLength = calculateMeanReadLength(libraryStats.lengthDistribution),
                  maxReadLength = libraryStats.lengthDistribution.keys.maxOption,
                  meanInsertSize = calculateMeanInsertSize(libraryStats.insertSizeDistribution),
                  files = updatedFiles,
                  alignmentRefs = if (seqRun.alignmentRefs.contains(alignUri)) seqRun.alignmentRefs else seqRun.alignmentRefs :+ alignUri
                )

                // Trigger facility lookup only if not already set (preserve manual edits)
                if (seqRun.sequencingFacility.isEmpty) {
                  instrumentId.foreach { id =>
                    lookupAndUpdateFacility(sampleAccession, resultIndex, id)
                  }
                }

                // Persist to H2 - alignment first, then sequence run
                import com.decodingus.service.EntityConversions.parseIdFromRef
                val persistResult = seqRun.atUri.flatMap(parseIdFromRef) match {
                  case Some(seqRunId) =>
                    h2Service.createAlignment(newAlignment, seqRunId) match {
                      case Right(persisted) =>
                        log.info(s" Alignment persisted to H2: ${persisted.atUri}")
                        // Only update sequence run if alignment succeeded
                        h2Service.updateSequenceRun(updatedSeqRun) match {
                          case Right(seqRunPersisted) =>
                            log.info(s" SequenceRun updated in H2: ${seqRunPersisted.atUri}")
                            Right(())
                          case Left(error) =>
                            log.error(s"Failed to update SequenceRun in H2: $error")
                            Left(s"Failed to update sequence run: $error")
                        }
                      case Left(error) =>
                        log.error(s"Failed to persist Alignment to H2: $error")
                        Left(s"Failed to create alignment: $error")
                    }
                  case None =>
                    log.warn(s"Cannot persist alignment: sequence run has no valid URI")
                    Left("Sequence run has no valid URI")
                }

                persistResult match {
                  case Right(_) =>
                    // Update in-memory state only on success
                    val updatedSequenceRuns = _workspace.value.main.sequenceRuns.map { sr =>
                      if (sr.atUri == seqRun.atUri) updatedSeqRun else sr
                    }
                    val updatedAlignments = _workspace.value.main.alignments :+ newAlignment
                    val updatedContent = _workspace.value.main.copy(
                      sequenceRuns = updatedSequenceRuns,
                      alignments = updatedAlignments
                    )
                    _workspace.value = _workspace.value.copy(main = updatedContent)

                    lastLibraryStats.value = Some(libraryStats)
                    analysisInProgress.value = false
                    analysisProgress.value = if (isNewRun) "Analysis complete" else s"Added ${libraryStats.referenceBuild} alignment to existing run"
                    analysisProgressPercent.value = 1.0
                    onProgress("Complete", 1.0)
                    onComplete(Right((resultIndex, libraryStats)))

                  case Left(error) =>
                    analysisInProgress.value = false
                    analysisProgress.value = s"Failed: $error"
                    onComplete(Left(error))
                }
              } // end runResult.foreach

            case None =>
              analysisInProgress.value = false
              onComplete(Left("Subject was removed during analysis"))
          }
        }

        libraryStats
      } catch {
        case e: Exception =>
          Platform.runLater {
            analysisInProgress.value = false
            analysisError.value = e.getMessage
            analysisProgress.value = s"Analysis failed: ${e.getMessage}"
            onComplete(Left(e.getMessage))
          }
          throw e
      }
    }
  }

  /** Infer test type from library stats */
  private def inferTestType(stats: LibraryStats): String = {
    // HiFi reads are typically very long (>10kb average)
    val avgReadLength = if (stats.lengthDistribution.nonEmpty) {
      val total = stats.lengthDistribution.map { case (len, count) => len.toLong * count }.sum
      val count = stats.lengthDistribution.values.sum
      if (count > 0) total / count else 0
    } else 0

    // Use codes that match TestTypes definitions for proper capability lookup
    if (stats.inferredPlatform == "PacBio" && avgReadLength > 10000) "WGS_HIFI"
    else if (stats.inferredPlatform == "PacBio") "WGS_CLR"
    else if (stats.inferredPlatform == "Oxford Nanopore") "WGS_NANOPORE"
    else "WGS" // Default assumption for Illumina/MGI
  }

  /** Calculate mean insert size from distribution */
  private def calculateMeanReadLength(distribution: Map[Int, Int]): Option[Int] = {
    if (distribution.isEmpty) None
    else {
      val totalReads = distribution.values.sum.toDouble
      val weightedSum = distribution.map { case (len, count) => len.toLong * count }.sum
      if (totalReads > 0) Some((weightedSum / totalReads).round.toInt) else None
    }
  }

  private def calculateMeanInsertSize(distribution: Map[Long, Int]): Option[Double] = {
    if (distribution.isEmpty) None
    else {
      val total = distribution.map { case (size, count) => size * count }.sum
      val count = distribution.values.sum
      if (count > 0) Some(total.toDouble / count) else None
    }
  }

  /** Gets all checksums for a subject's sequence run and alignment files */
  def getExistingChecksums(sampleAccession: String): Set[String] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        val seqRunChecksums = sequenceRuns.flatMap(_.files.flatMap(_.checksum))
        // Also include alignment file checksums
        val alignmentChecksums = sequenceRuns.flatMap { sr =>
          _workspace.value.main.getAlignmentsForSequenceRun(sr).flatMap(_.files.flatMap(_.checksum))
        }
        (seqRunChecksums ++ alignmentChecksums).toSet
      case None =>
        Set.empty
    }
  }

  /** Gets a specific sequence run by index for a subject */
  def getSequenceRun(sampleAccession: String, index: Int): Option[SequenceRun] = {
    findSubject(sampleAccession).flatMap { subject =>
      _workspace.value.main.getSequenceRunsForBiosample(subject).lift(index)
    }
  }

  /**
   * Removes a sequence run from a subject by index.
   * Uses SequenceDataManager for atomic H2 deletion + in-memory state update.
   */
  def removeSequenceData(sampleAccession: String, index: Int): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (index >= 0 && index < sequenceRuns.size) {
          val seqRunToRemove = sequenceRuns(index)
          seqRunToRemove.atUri match {
            case Some(uri) =>
              sequenceDataManager.deleteSequenceRun(sampleAccession, uri) match {
                case Right(deleted) =>
                  if (deleted) log.info(s"SequenceRun deleted via SequenceDataManager: $uri")
                  else log.warn(s"SequenceRun not found for deletion: $uri")
                case Left(error) =>
                  log.error(s"Failed to delete SequenceRun: $error")
              }
            case None =>
              log.warn(s"Cannot delete sequence run at index $index: no URI")
          }
        } else {
          log.info(s"Cannot remove sequence run: index $index out of bounds")
        }
      case None =>
        log.info(s"Cannot remove sequence run: subject $sampleAccession not found")
    }
  }

  /**
   * Updates a sequence run at a specific index for a subject.
   * Uses SequenceDataManager for atomic H2 update + in-memory state update.
   */
  def updateSequenceRun(sampleAccession: String, index: Int, updatedRun: SequenceRun): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (index >= 0 && index < sequenceRuns.size) {
          val originalRun = sequenceRuns(index)
          // Ensure the updated run has the original URI and updated metadata
          val runToUpdate = updatedRun.copy(
            atUri = originalRun.atUri,
            meta = originalRun.meta.updated("edit")
          )

          sequenceDataManager.updateSequenceRun(runToUpdate) match {
            case Right(persisted) =>
              log.info(s"SequenceRun updated via SequenceDataManager: ${persisted.atUri}")
            case Left(error) =>
              log.error(s"Failed to update SequenceRun: $error")
          }
        } else {
          log.info(s"Cannot update sequence run: index $index out of bounds")
        }
      case None =>
        log.info(s" Cannot update sequence run: subject $sampleAccession not found")
    }
  }

  /**
   * Looks up the sequencing facility for an instrument ID via the API
   * and updates the sequence run if found.
   *
   * This runs asynchronously in the background and updates the workspace
   * when the lookup completes. Failures are logged but don't interrupt
   * the main workflow.
   *
   * @param sampleAccession The subject's accession ID
   * @param index           The sequence run index
   * @param instrumentId    The instrument identifier to look up
   */
  private def lookupAndUpdateFacility(sampleAccession: String, index: Int, instrumentId: String): Unit = {
    DecodingUsClient.lookupLabByInstrument(instrumentId).onComplete {
      case Success(Some(labInfo)) =>
        Platform.runLater {
          getSequenceRun(sampleAccession, index).foreach { seqRun =>
            // Only update if facility not already set
            if (seqRun.sequencingFacility.isEmpty) {
              val updatedRun = seqRun.copy(
                sequencingFacility = Some(labInfo.labName),
                // Also update instrument model if available from API
                instrumentModel = labInfo.model.orElse(seqRun.instrumentModel)
              )
              updateSequenceRun(sampleAccession, index, updatedRun)
              log.info(s" Updated facility for instrument $instrumentId: ${labInfo.labName}")
            }
          }
        }
      case Success(None) =>
        log.info(s" No facility found for instrument ID: $instrumentId")
      case Failure(error) =>
        log.error(s"Failed to lookup facility for instrument $instrumentId: ${error.getMessage}")
    }
  }

  /**
   * Manually triggers a facility lookup for a sequence run.
   * Useful if the initial lookup failed or needs to be refreshed.
   */
  def refreshFacilityLookup(sampleAccession: String, index: Int): Unit = {
    getSequenceRun(sampleAccession, index).foreach { seqRun =>
      seqRun.instrumentId.foreach { instrumentId =>
        lookupAndUpdateFacility(sampleAccession, index, instrumentId)
      }
    }
  }

  // --- Run Fingerprint Matching ---

  /**
   * Find an existing sequence run that matches the given fingerprint.
   * Delegates to FingerprintMatchService for matching logic.
   *
   * @param biosampleRef The biosample to search within
   * @param fingerprint  The computed fingerprint to match
   * @param libraryStats Full stats for additional matching criteria
   * @return Match result with confidence level
   */
  def findMatchingSequenceRun(
                               biosampleRef: String,
                               fingerprint: String,
                               libraryStats: LibraryStats
                             ): FingerprintMatchResult = {
    val candidateRuns = _workspace.value.main.sequenceRuns
      .filter(_.biosampleRef == biosampleRef)
      .zipWithIndex
      .toList

    fingerprintMatchService.findMatch(candidateRuns, fingerprint, libraryStats)
  }

  /**
   * Add an alignment to an existing sequence run (for multi-reference scenarios).
   * Uses SequenceDataManager for atomic alignment creation + in-memory state update.
   *
   * @param sampleAccession  The biosample's accession ID
   * @param sequenceRunIndex Index of the sequence run to add alignment to
   * @param newAlignment     The new alignment to add
   * @param fileInfo         The file info for the new alignment
   */
  def addAlignmentToExistingRun(
                                 sampleAccession: String,
                                 sequenceRunIndex: Int,
                                 newAlignment: Alignment,
                                 fileInfo: FileInfo
                               ): Unit = {
    findSubject(sampleAccession).foreach { subject =>
      val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
      if (sequenceRunIndex >= 0 && sequenceRunIndex < sequenceRuns.size) {
        val seqRun = sequenceRuns(sequenceRunIndex)

        seqRun.atUri match {
          case Some(seqRunUri) =>
            sequenceDataManager.createAlignment(seqRunUri, newAlignment, Some(fileInfo)) match {
              case Right(result) =>
                log.info(s"Alignment created via SequenceDataManager: ${result.alignment.atUri}")
              case Left(error) =>
                log.error(s"Failed to create Alignment: $error")
            }
          case None =>
            log.warn(s"Cannot add alignment: sequence run at index $sequenceRunIndex has no URI")
        }
      } else {
        log.info(s"Cannot add alignment: index $sequenceRunIndex out of bounds")
      }
    }
  }

  /**
   * Merges two sequence runs that represent the same source sequencing data.
   * The secondary run's alignments are moved to the primary run, then the secondary is deleted.
   *
   * Use case: Same HiFi/Illumina data aligned to different references (GRCh38, CHM13, etc.)
   *
   * @param sampleAccession  The biosample's accession ID
   * @param primaryIndex     Index of the run to keep (alignments will be added here)
   * @param secondaryIndex   Index of the run to merge and delete
   * @return Either error message or success with count of moved alignments
   */
  def mergeSequenceRuns(
    sampleAccession: String,
    primaryIndex: Int,
    secondaryIndex: Int
  ): Either[String, Int] = {
    import com.decodingus.service.EntityConversions.parseIdFromRef

    findSubject(sampleAccession) match {
      case None =>
        Left(s"Subject $sampleAccession not found")

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)

        if (primaryIndex < 0 || primaryIndex >= sequenceRuns.size) {
          return Left(s"Primary run index $primaryIndex out of bounds")
        }
        if (secondaryIndex < 0 || secondaryIndex >= sequenceRuns.size) {
          return Left(s"Secondary run index $secondaryIndex out of bounds")
        }
        if (primaryIndex == secondaryIndex) {
          return Left("Cannot merge a run with itself")
        }

        val primaryRun = sequenceRuns(primaryIndex)
        val secondaryRun = sequenceRuns(secondaryIndex)

        val primaryUri = primaryRun.atUri match {
          case Some(uri) => uri
          case None => return Left("Primary run has no URI")
        }

        val secondaryUri = secondaryRun.atUri match {
          case Some(uri) => uri
          case None => return Left("Secondary run has no URI")
        }

        val secondaryId = parseIdFromRef(secondaryUri) match {
          case Some(id) => id
          case None => return Left(s"Invalid secondary run URI: $secondaryUri")
        }

        log.info(s"Merging sequence runs for $sampleAccession: secondary[$secondaryIndex] -> primary[$primaryIndex]")

        // Get alignments from the secondary run
        val secondaryAlignments = _workspace.value.main.getAlignmentsForSequenceRun(secondaryRun)
        var movedCount = 0
        val movedAlignmentUris = scala.collection.mutable.ArrayBuffer[String]()

        // Step 1: Update each alignment's sequenceRunRef to point to the primary run (DB only)
        // h2Service methods don't trigger in-memory workspace updates
        secondaryAlignments.foreach { alignment =>
          val updatedAlignment = alignment.copy(
            sequenceRunRef = primaryUri,
            meta = alignment.meta.updated("merge")
          )

          // Update alignment in database only
          h2Service.updateAlignment(updatedAlignment) match {
            case Right(_) =>
              movedCount += 1
              alignment.atUri.foreach(uri => movedAlignmentUris += uri)
              log.info(s"Moved alignment ${alignment.atUri} from $secondaryUri to $primaryUri")
            case Left(error) =>
              log.error(s"Failed to update alignment ${alignment.atUri}: $error")
          }
        }

        // Step 2: Update primary run's alignmentRefs to include secondary's alignments
        val updatedAlignmentRefs = primaryRun.alignmentRefs ++ secondaryRun.alignmentRefs.filterNot(primaryRun.alignmentRefs.contains)

        // Move files from secondary to primary
        val updatedFiles = primaryRun.files ++ secondaryRun.files.filterNot { sf =>
          primaryRun.files.exists(pf => pf.checksum == sf.checksum && sf.checksum.isDefined)
        }

        // Merge fingerprint data (keep primary's data, fill in gaps from secondary)
        val updatedPrimaryRun = primaryRun.copy(
          alignmentRefs = updatedAlignmentRefs,
          files = updatedFiles,
          sampleName = primaryRun.sampleName.orElse(secondaryRun.sampleName),
          libraryId = primaryRun.libraryId.orElse(secondaryRun.libraryId),
          platformUnit = primaryRun.platformUnit.orElse(secondaryRun.platformUnit),
          runFingerprint = primaryRun.runFingerprint.orElse(secondaryRun.runFingerprint),
          meta = primaryRun.meta.updated("merge")
        )

        // Step 3: Update the primary run in DB
        h2Service.updateSequenceRun(updatedPrimaryRun) match {
          case Right(_) =>
            log.info(s"Updated primary run with ${updatedAlignmentRefs.size} alignment refs")
          case Left(error) =>
            log.error(s"Failed to update primary run: $error")
            return Left(s"Failed to update primary run: $error")
        }

        // Step 4: Delete the secondary run from DB (takes UUID)
        h2Service.deleteSequenceRun(secondaryId) match {
          case Right(deleted) =>
            if (deleted) {
              log.info(s"Deleted secondary run: $secondaryUri")
            } else {
              log.warn(s"Secondary run not found for deletion: $secondaryUri")
            }
          case Left(error) =>
            log.error(s"Failed to delete secondary run: $error")
            // Continue anyway - alignments have been moved
        }

        // Step 5: Manually update in-memory workspace state in one atomic operation
        // This avoids the multiple syncBuffers calls that can cause race conditions
        val workspace = _workspace.value

        // Update alignments to point to primary
        val updatedAlignments = workspace.main.alignments.map { alignment =>
          if (movedAlignmentUris.contains(alignment.atUri.getOrElse(""))) {
            alignment.copy(sequenceRunRef = primaryUri)
          } else {
            alignment
          }
        }

        // Update the primary sequence run and remove secondary
        val updatedSequenceRuns = workspace.main.sequenceRuns
          .filterNot(_.atUri.contains(secondaryUri)) // Remove secondary
          .map { sr =>
            if (sr.atUri == primaryRun.atUri) updatedPrimaryRun else sr // Update primary
          }

        // Update biosample's sequenceRunRefs to remove secondary
        val updatedSamples = workspace.main.samples.map { sample =>
          if (sample.sampleAccession == sampleAccession) {
            sample.copy(sequenceRunRefs = sample.sequenceRunRefs.filterNot(_ == secondaryUri))
          } else {
            sample
          }
        }

        val updatedContent = workspace.main.copy(
          samples = updatedSamples,
          sequenceRuns = updatedSequenceRuns,
          alignments = updatedAlignments
        )
        _workspace.value = workspace.copy(main = updatedContent)

        log.info(s"Merge complete: moved $movedCount alignments, merged ${secondaryRun.files.size} files")
        Right(movedCount)
    }
  }

  // --- Analysis State ---
  // Observable properties for tracking analysis progress
  val analysisInProgress: BooleanProperty = BooleanProperty(false)
  val analysisProgress: StringProperty = StringProperty("")
  val analysisProgressPercent: DoubleProperty = DoubleProperty(0.0)
  val analysisError: StringProperty = StringProperty("")

  // Store analysis results for UI access
  val lastLibraryStats: ObjectProperty[Option[LibraryStats]] = ObjectProperty(None)
  val lastWgsMetrics: ObjectProperty[Option[WgsMetrics]] = ObjectProperty(None)

  // --- Reference Download State ---
  // These properties allow the UI to show a prompt when a reference download is needed
  sealed trait ReferenceDownloadRequest

  case class PendingDownload(build: String, url: String, sizeMB: Int, onConfirm: () => Unit, onCancel: () => Unit) extends ReferenceDownloadRequest

  case object NoDownloadPending extends ReferenceDownloadRequest

  val pendingReferenceDownload: ObjectProperty[ReferenceDownloadRequest] = ObjectProperty(NoDownloadPending)

  /**
   * Checks if a reference is available and resolves it.
   * If prompting is enabled and download is required, sets pendingReferenceDownload for UI to handle.
   *
   * @param referenceBuild The build to resolve (e.g., "GRCh38")
   * @param onProgress     Progress callback for download
   * @param onResolved     Called with the resolved path when available
   * @param onError        Called if resolution fails or is cancelled
   */
  def resolveReferenceWithPrompt(
                                  referenceBuild: String,
                                  onProgress: (Long, Long) => Unit,
                                  onResolved: String => Unit,
                                  onError: String => Unit
                                ): Unit = {
    val config = ReferenceConfigService.load()
    val referenceGateway = new ReferenceGateway(onProgress)

    referenceGateway.checkAvailability(referenceBuild) match {
      case ReferenceResolveResult.Available(path) =>
        onResolved(path.toString)

      case ReferenceResolveResult.DownloadRequired(build, url, sizeMB) =>
        if (config.promptBeforeDownload) {
          // Set pending download for UI to handle
          Platform.runLater {
            pendingReferenceDownload.value = PendingDownload(
              build = build,
              url = url,
              sizeMB = sizeMB,
              onConfirm = () => {
                pendingReferenceDownload.value = NoDownloadPending
                // User confirmed - proceed with download in background
                Future {
                  referenceGateway.downloadAndResolve(build) match {
                    case Right(path) =>
                      Platform.runLater {
                        onResolved(path.toString)
                      }
                    case Left(error) =>
                      Platform.runLater {
                        onError(error)
                      }
                  }
                }
              },
              onCancel = () => {
                pendingReferenceDownload.value = NoDownloadPending
                onError(s"Reference download cancelled. Configure a local path in Settings or enable auto-download.")
              }
            )
          }
        } else {
          // Auto-download enabled - proceed without prompting
          Future {
            referenceGateway.downloadAndResolve(build) match {
              case Right(path) =>
                Platform.runLater {
                  onResolved(path.toString)
                }
              case Left(error) =>
                Platform.runLater {
                  onError(error)
                }
            }
          }
        }

      case ReferenceResolveResult.Error(message) =>
        onError(message)
    }
  }

  // --- Analysis Operations ---

  /**
   * Runs the full initial analysis pipeline for a sequencing run:
   * 1. Detects reference build from BAM/CRAM header
   * 2. Resolves/downloads reference genome if needed
   * 3. Collects library statistics
   * 4. Updates the SequenceRun with results and creates Alignment
   *
   * @param sampleAccession  The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param onComplete       Callback when analysis completes (success or failure)
   */
  def runInitialAnalysis(
                          sampleAccession: String,
                          sequenceRunIndex: Int,
                          onComplete: Either[String, LibraryStats] => Unit
                        ): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case Some(seqRun) =>
            seqRun.files.headOption match {
              case Some(fileInfo) =>
                val bamPath = fileInfo.location.getOrElse("")
                log.info(s" Starting initial analysis for ${fileInfo.fileName}")

                analysisInProgress.value = true
                analysisError.value = ""
                analysisProgress.value = "Initializing analysis..."
                analysisProgressPercent.value = 0.0

                // Run analysis in background thread
                Future {
                  try {
                    // Step 1: Detect reference build from header
                    updateProgress("Reading BAM/CRAM header...", 0.1)
                    val header = SamReaderFactory.makeDefault().open(new File(bamPath)).getFileHeader
                    val libraryStatsProcessor = new LibraryStatsProcessor()
                    val referenceBuild = libraryStatsProcessor.detectReferenceBuild(header)

                    if (referenceBuild == "Unknown") {
                      throw new IllegalStateException("Could not determine reference build from BAM/CRAM header.")
                    }

                    // Step 2: Resolve reference genome
                    updateProgress(s"Resolving reference: $referenceBuild", 0.2)
                    val referenceGateway = new ReferenceGateway((done, total) => {
                      val pct = if (total > 0) 0.2 + (done.toDouble / total) * 0.3 else 0.2
                      updateProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct)
                    })

                    val referencePath = referenceGateway.resolve(referenceBuild) match {
                      case Right(path) => path.toString
                      case Left(error) => throw new Exception(s"Failed to resolve reference: $error")
                    }

                    // Step 3: Collect library stats
                    updateProgress("Analyzing library statistics...", 0.5)
                    val libraryStats = libraryStatsProcessor.process(bamPath, referencePath, (message, current, total) => {
                      val pct = 0.5 + (current.toDouble / total) * 0.4
                      updateProgress(s"Library Stats: $message", pct)
                    })

                    // Step 4: Update SequenceRun and create/update Alignment
                    updateProgress("Saving results...", 0.95)
                    Platform.runLater {
                      // Create or update alignment - find one matching THIS reference build
                      val existingAlignment = seqRun.alignmentRefs.flatMap { ref =>
                        _workspace.value.main.alignments.find(a => a.atUri.contains(ref) && a.referenceBuild == libraryStats.referenceBuild)
                      }.headOption

                      val alignUri = existingAlignment.flatMap(_.atUri).getOrElse(
                        s"local:alignment:${subject.sampleAccession}:${libraryStats.referenceBuild}:${java.util.UUID.randomUUID().toString.take(8)}"
                      )
                      val newAlignment = Alignment(
                        atUri = Some(alignUri),
                        meta = existingAlignment.map(_.meta.updated("analysis")).getOrElse(RecordMeta.initial),
                        sequenceRunRef = seqRun.atUri.getOrElse(""),
                        biosampleRef = Some(subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")),
                        referenceBuild = libraryStats.referenceBuild,
                        aligner = libraryStats.aligner,
                        files = seqRun.files,
                        metrics = existingAlignment.flatMap(_.metrics)
                      )

                      val updatedSeqRun = seqRun.copy(
                        meta = seqRun.meta.updated("analysis"),
                        platformName = if (seqRun.platformName == "Unknown" || seqRun.platformName == "Other") libraryStats.inferredPlatform else seqRun.platformName,
                        instrumentModel = seqRun.instrumentModel.orElse(Some(libraryStats.mostFrequentInstrument)),
                        testType = inferTestType(libraryStats),
                        libraryId = if (libraryStats.libraryId != "Unknown") Some(libraryStats.libraryId) else seqRun.libraryId,
                        platformUnit = libraryStats.platformUnit.orElse(seqRun.platformUnit),
                        libraryLayout = Some(if (libraryStats.pairedReads > libraryStats.readCount / 2) "Paired-End" else "Single-End"),
                        totalReads = Some(libraryStats.readCount.toLong),
                        readLength = calculateMeanReadLength(libraryStats.lengthDistribution).orElse(seqRun.readLength),
                        maxReadLength = libraryStats.lengthDistribution.keys.maxOption.orElse(seqRun.maxReadLength),
                        meanInsertSize = calculateMeanInsertSize(libraryStats.insertSizeDistribution).orElse(seqRun.meanInsertSize),
                        alignmentRefs = if (seqRun.alignmentRefs.contains(alignUri)) seqRun.alignmentRefs else seqRun.alignmentRefs :+ alignUri
                      )

                      // Persist to H2 atomically
                      import com.decodingus.service.EntityConversions.parseIdFromRef
                      seqRun.atUri.flatMap(parseIdFromRef).foreach { seqRunId =>
                        if (existingAlignment.isDefined) {
                          h2Service.updateAlignment(newAlignment) match {
                            case Right(persisted) =>
                              log.info(s" Alignment updated in H2: ${persisted.atUri}")
                            case Left(error) =>
                              log.error(s"Failed to update Alignment in H2: $error")
                          }
                        } else {
                          h2Service.createAlignment(newAlignment, seqRunId) match {
                            case Right(persisted) =>
                              log.info(s" Alignment created in H2: ${persisted.atUri}")
                            case Left(error) =>
                              log.error(s"Failed to create Alignment in H2: $error")
                          }
                        }
                      }
                      h2Service.updateSequenceRun(updatedSeqRun) match {
                        case Right(persisted) =>
                          log.info(s" SequenceRun updated in H2: ${persisted.atUri}")
                        case Left(error) =>
                          log.error(s"Failed to update SequenceRun in H2: $error")
                      }

                      // Update in-memory state
                      val updatedSequenceRuns = _workspace.value.main.sequenceRuns.map { sr =>
                        if (sr.atUri == seqRun.atUri) updatedSeqRun else sr
                      }
                      val updatedAlignments = if (existingAlignment.isDefined) {
                        _workspace.value.main.alignments.map { a =>
                          if (a.atUri.contains(alignUri)) newAlignment else a
                        }
                      } else {
                        _workspace.value.main.alignments :+ newAlignment
                      }
                      val updatedContent = _workspace.value.main.copy(
                        sequenceRuns = updatedSequenceRuns,
                        alignments = updatedAlignments
                      )
                      _workspace.value = _workspace.value.copy(main = updatedContent)

                      lastLibraryStats.value = Some(libraryStats)
                      analysisInProgress.value = false
                      analysisProgress.value = "Analysis complete"
                      analysisProgressPercent.value = 1.0
                      onComplete(Right(libraryStats))
                    }

                    libraryStats
                  } catch {
                    case e: Exception =>
                      Platform.runLater {
                        analysisInProgress.value = false
                        analysisError.value = e.getMessage
                        analysisProgress.value = s"Analysis failed: ${e.getMessage}"
                        onComplete(Left(e.getMessage))
                      }
                      throw e
                  }
                }

              case None =>
                onComplete(Left("No alignment file associated with this sequence run"))
            }

          case None =>
            onComplete(Left(s"Sequence run not found at index $sequenceRunIndex for subject $sampleAccession"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /**
   * Runs WGS metrics analysis (deep coverage) for a sequencing run.
   * Requires that initial analysis has already been run.
   *
   * @param sampleAccession  The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param onComplete       Callback when analysis completes
   */
  def runWgsMetricsAnalysis(
                             sampleAccession: String,
                             sequenceRunIndex: Int,
                             onComplete: Either[String, WgsMetrics] => Unit
                           ): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case Some(seqRun) =>
            // Need alignment with reference build
            val alignments = _workspace.value.main.getAlignmentsForSequenceRun(seqRun)
            alignments.headOption match {
              case Some(alignment) =>
                seqRun.files.headOption match {
                  case Some(fileInfo) =>
                    val bamPath = fileInfo.location.getOrElse("")
                    log.info(s" Starting WGS metrics analysis for ${fileInfo.fileName}")

                    analysisInProgress.value = true
                    analysisError.value = ""
                    analysisProgress.value = "Starting WGS metrics analysis..."
                    analysisProgressPercent.value = 0.0

                    Future {
                      try {
                        // Resolve reference
                        updateProgress("Resolving reference genome...", 0.1)
                        val referenceGateway = new ReferenceGateway((_, _) => {})
                        val referencePath = referenceGateway.resolve(alignment.referenceBuild) match {
                          case Right(path) => path.toString
                          case Left(error) => throw new Exception(s"Failed to resolve reference: $error")
                        }

                        // Run GATK CollectWgsMetrics
                        updateProgress("Running GATK CollectWgsMetrics (this may take a while)...", 0.2)
                        val wgsProcessor = new WgsMetricsProcessor()
                        val artifactCtx = ArtifactContext(
                          sampleAccession = sampleAccession,
                          sequenceRunUri = seqRun.atUri,
                          alignmentUri = alignment.atUri
                        )
                        // Single-end reads need COUNT_UNPAIRED (long-read platforms + YSEQ WGS400 with 400bp SE reads)
                        val isLongReadPlatform = seqRun.testType.toUpperCase match {
                          case t if t.contains("HIFI") || t.contains("CLR") || t.contains("NANOPORE") => true
                          case _ => false
                        }
                        val isSingleEnd = isLongReadPlatform || seqRun.libraryLayout.exists(_.toLowerCase == "single-end")
                        val wgsMetrics = wgsProcessor.process(
                          bamPath,
                          referencePath,
                          (message, current, total) => {
                            val pct = 0.2 + (current / total) * 0.7
                            updateProgress(message, pct)
                          },
                          seqRun.maxReadLength, // Pass max read length to handle long reads (e.g., PacBio HiFi, NovaSeq 151bp)
                          Some(artifactCtx),
                          seqRun.totalReads,
                          countUnpaired = isSingleEnd
                        ) match {
                          case Right(metrics) => metrics
                          case Left(error) => throw error
                        }

                        // Update alignment metrics
                        updateProgress("Saving results...", 0.95)
                        Platform.runLater {
                          val existingMetrics = alignment.metrics.getOrElse(AlignmentMetrics())
                          val updatedMetrics = existingMetrics.copy(
                            genomeTerritory = Some(wgsMetrics.genomeTerritory),
                            meanCoverage = Some(wgsMetrics.meanCoverage),
                            medianCoverage = Some(wgsMetrics.medianCoverage),
                            sdCoverage = Some(wgsMetrics.sdCoverage),
                            pctExcDupe = Some(wgsMetrics.pctExcDupe),
                            pctExcMapq = Some(wgsMetrics.pctExcMapq),
                            pct10x = Some(wgsMetrics.pct10x),
                            pct20x = Some(wgsMetrics.pct20x),
                            pct30x = Some(wgsMetrics.pct30x),
                            hetSnpSensitivity = Some(wgsMetrics.hetSnpSensitivity)
                          )
                          val updatedAlignment = alignment.copy(
                            meta = alignment.meta.updated("wgsMetrics"),
                            metrics = Some(updatedMetrics)
                          )

                          // Persist to H2 atomically
                          h2Service.updateAlignment(updatedAlignment) match {
                            case Right(persisted) =>
                              log.info(s" Alignment metrics updated in H2: ${persisted.atUri}")
                            case Left(error) =>
                              log.error(s"Failed to update Alignment metrics in H2: $error")
                          }

                          // Update in-memory state
                          val updatedAlignments = _workspace.value.main.alignments.map { a =>
                            if (a.atUri == alignment.atUri) updatedAlignment else a
                          }
                          val updatedContent = _workspace.value.main.copy(alignments = updatedAlignments)
                          _workspace.value = _workspace.value.copy(main = updatedContent)

                          lastWgsMetrics.value = Some(wgsMetrics)
                          analysisInProgress.value = false
                          analysisProgress.value = "WGS metrics analysis complete"
                          analysisProgressPercent.value = 1.0
                          onComplete(Right(wgsMetrics))
                        }

                        wgsMetrics
                      } catch {
                        case e: Exception =>
                          Platform.runLater {
                            analysisInProgress.value = false
                            analysisError.value = e.getMessage
                            analysisProgress.value = s"Analysis failed: ${e.getMessage}"
                            onComplete(Left(e.getMessage))
                          }
                          throw e
                      }
                    }

                  case None =>
                    onComplete(Left("No alignment file associated with this sequence run"))
                }

              case None =>
                onComplete(Left("Please run initial analysis first to detect reference build"))
            }

          case None =>
            onComplete(Left(s"Sequence run not found at index $sequenceRunIndex for subject $sampleAccession"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /** Helper to update progress on the JavaFX thread */
  private def updateProgress(message: String, percent: Double): Unit = {
    Platform.runLater {
      analysisProgress.value = message
      analysisProgressPercent.value = percent
    }
  }

  // --- Callable Loci Analysis ---

  // Store last callable loci result
  val lastCallableLociResult: ObjectProperty[Option[CallableLociResult]] = ObjectProperty(None)

  /**
   * Runs callable loci analysis for a sequencing run.
   * Analyzes each contig to determine callable vs non-callable regions.
   * Requires that initial analysis has already been run.
   *
   * @param sampleAccession  The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param onComplete       Callback when analysis completes, returns result and artifact directory path
   */
  def runCallableLociAnalysis(
                               sampleAccession: String,
                               sequenceRunIndex: Int,
                               onComplete: Either[String, (CallableLociResult, java.nio.file.Path)] => Unit
                             ): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case Some(seqRun) =>
            // Need alignment with reference build
            val alignments = _workspace.value.main.getAlignmentsForSequenceRun(seqRun)
            alignments.headOption match {
              case Some(alignment) =>
                seqRun.files.headOption match {
                  case Some(fileInfo) =>
                    val bamPath = fileInfo.location.getOrElse("")
                    log.info(s" Starting callable loci analysis for ${fileInfo.fileName}")

                    analysisInProgress.value = true
                    analysisError.value = ""
                    analysisProgress.value = "Starting callable loci analysis..."
                    analysisProgressPercent.value = 0.0

                    Future {
                      try {
                        // Resolve reference
                        updateProgress("Resolving reference genome...", 0.05)
                        val referenceGateway = new ReferenceGateway((_, _) => {})
                        val referencePath = referenceGateway.resolve(alignment.referenceBuild) match {
                          case Right(path) => path.toString
                          case Left(error) => throw new Exception(s"Failed to resolve reference: $error")
                        }

                        // Run GATK CallableLoci
                        updateProgress("Running GATK CallableLoci (analyzing contigs)...", 0.1)
                        val callableLociProcessor = new CallableLociProcessor()
                        val artifactCtx = ArtifactContext(
                          sampleAccession = sampleAccession,
                          sequenceRunUri = seqRun.atUri,
                          alignmentUri = alignment.atUri
                        )
                        val artifactDir = artifactCtx.getArtifactDir
                        // HiFi long reads have higher per-base accuracy, so 2x coverage is callable
                        val minDepth = seqRun.testType.toUpperCase match {
                          case t if t.contains("HIFI") => 2
                          case _ => 4 // Default for short reads
                        }
                        val (result, svgStrings) = callableLociProcessor.process(
                          bamPath,
                          referencePath,
                          (message, current, total) => {
                            val pct = 0.1 + (current.toDouble / total.toDouble) * 0.85
                            updateProgress(message, pct)
                          },
                          Some(artifactCtx),
                          minDepth
                        ) match {
                          case Right(r) => r
                          case Left(error) => throw error
                        }

                        // Update alignment metrics with callable bases count
                        updateProgress("Saving results...", 0.98)
                        Platform.runLater {
                          val existingMetrics = alignment.metrics.getOrElse(AlignmentMetrics())
                          val updatedMetrics = existingMetrics.copy(
                            callableBases = Some(result.callableBases),
                            contigs = result.contigAnalysis.map { cs =>
                              ContigMetrics(
                                contigName = cs.contigName,
                                callable = cs.callable,
                                noCoverage = cs.noCoverage,
                                lowCoverage = cs.lowCoverage,
                                excessiveCoverage = cs.excessiveCoverage,
                                poorMappingQuality = cs.poorMappingQuality
                              )
                            }
                          )
                          val updatedAlignment = alignment.copy(
                            meta = alignment.meta.updated("callableLoci"),
                            metrics = Some(updatedMetrics)
                          )

                          // Persist to H2 atomically
                          h2Service.updateAlignment(updatedAlignment) match {
                            case Right(persisted) =>
                              log.info(s" Alignment callable loci updated in H2: ${persisted.atUri}")
                            case Left(error) =>
                              log.error(s"Failed to update Alignment callable loci in H2: $error")
                          }

                          // Update in-memory state
                          val updatedAlignments = _workspace.value.main.alignments.map { a =>
                            if (a.atUri == alignment.atUri) updatedAlignment else a
                          }
                          val updatedContent = _workspace.value.main.copy(alignments = updatedAlignments)
                          _workspace.value = _workspace.value.copy(main = updatedContent)

                          lastCallableLociResult.value = Some(result)
                          analysisInProgress.value = false
                          analysisProgress.value = "Callable loci analysis complete"
                          analysisProgressPercent.value = 1.0
                          onComplete(Right((result, artifactDir)))
                        }

                        result
                      } catch {
                        case e: Exception =>
                          Platform.runLater {
                            analysisInProgress.value = false
                            analysisError.value = e.getMessage
                            analysisProgress.value = s"Analysis failed: ${e.getMessage}"
                            onComplete(Left(e.getMessage))
                          }
                          throw e
                      }
                    }

                  case None =>
                    onComplete(Left("No alignment file associated with this sequence run"))
                }

              case None =>
                onComplete(Left("Please run initial analysis first to detect reference build"))
            }

          case None =>
            onComplete(Left(s"Sequence run not found at index $sequenceRunIndex for subject $sampleAccession"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  // --- Haplogroup Analysis ---

  // Store last haplogroup analysis result
  val lastHaplogroupResult: ObjectProperty[Option[AnalysisHaplogroupResult]] = ObjectProperty(None)

  /**
   * Runs haplogroup analysis for a subject using the specified tree type.
   * Uses FTDNA tree provider for both Y-DNA and MT-DNA.
   *
   * @param sampleAccession  The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param treeType         The type of haplogroup tree (YDNA or MTDNA)
   * @param onComplete       Callback when analysis completes
   */
  def runHaplogroupAnalysis(
                             sampleAccession: String,
                             sequenceRunIndex: Int,
                             treeType: TreeType,
                             onComplete: Either[String, AnalysisHaplogroupResult] => Unit
                           ): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case Some(seqRun) =>
            // Need alignment with reference build for haplogroup analysis
            val alignments = _workspace.value.main.getAlignmentsForSequenceRun(seqRun)
            alignments.headOption match {
              case Some(alignment) =>
                seqRun.files.headOption match {
                  case Some(fileInfo) =>
                    val bamPath = fileInfo.location.getOrElse("")
                    log.info(s" Starting ${treeType} haplogroup analysis for ${fileInfo.fileName}")

                    analysisInProgress.value = true
                    analysisError.value = ""
                    analysisProgress.value = "Starting haplogroup analysis..."
                    analysisProgressPercent.value = 0.0

                    Future {
                      try {
                        // Build LibraryStats from existing data for the processor
                        val libraryStats = LibraryStats(
                          readCount = seqRun.totalReads.map(_.toInt).getOrElse(0),
                          pairedReads = 0,
                          lengthDistribution = Map.empty,
                          insertSizeDistribution = Map.empty,
                          aligner = alignment.aligner,
                          referenceBuild = alignment.referenceBuild,
                          sampleName = subject.donorIdentifier,
                          flowCells = Map.empty,
                          instruments = Map.empty,
                          mostFrequentInstrument = seqRun.instrumentModel.getOrElse("Unknown"),
                          inferredPlatform = seqRun.platformName,
                          platformCounts = Map.empty
                        )

                        updateProgress("Loading haplogroup tree...", 0.1)

                        val processor = new HaplogroupProcessor()
                        val artifactCtx = ArtifactContext(
                          sampleAccession = sampleAccession,
                          sequenceRunUri = seqRun.atUri,
                          alignmentUri = alignment.atUri
                        )
                        // Select tree provider based on user preferences
                        val treeProviderType = treeType match {
                          case TreeType.YDNA =>
                            if (UserPreferencesService.getYdnaTreeProvider.equalsIgnoreCase("decodingus"))
                              TreeProviderType.DECODINGUS
                            else TreeProviderType.FTDNA
                          case TreeType.MTDNA =>
                            if (UserPreferencesService.getMtdnaTreeProvider.equalsIgnoreCase("decodingus"))
                              TreeProviderType.DECODINGUS
                            else TreeProviderType.FTDNA
                        }

                        // Extract biosampleId for YProfile population
                        val biosampleId = subject.atUri.flatMap { uri =>
                          scala.util.Try(UUID.fromString(uri.split("/").last)).toOption
                        }

                        // Infer source type from platform/test type
                        val yProfileSourceType = {
                          val testType = seqRun.testType.toLowerCase
                          val platform = seqRun.platformName.toLowerCase
                          if (testType.contains("hifi") || testType.contains("pacbio") || platform.contains("pacbio")) {
                            YProfileSourceType.WGS_LONG_READ
                          } else if (testType.contains("nanopore") || platform.contains("nanopore")) {
                            YProfileSourceType.WGS_LONG_READ
                          } else if (testType.contains("targeted") || testType.contains("panel")) {
                            YProfileSourceType.TARGETED_NGS
                          } else {
                            YProfileSourceType.WGS_SHORT_READ
                          }
                        }

                        val result = processor.analyze(
                          bamPath = bamPath,
                          libraryStats = libraryStats,
                          treeType = treeType,
                          treeProviderType = treeProviderType,
                          onProgress = (message, current, total) => {
                            val pct = if (total > 0) current / total else 0.0
                            updateProgress(message, pct)
                          },
                          artifactContext = Some(artifactCtx),
                          yProfileService = yProfileService,
                          biosampleId = biosampleId,
                          yProfileSourceType = Some(yProfileSourceType)
                        )

                        result match {
                          case Right(results) if results.nonEmpty =>
                            val topResult = results.head
                            Platform.runLater {
                              // Determine technology based on test type
                              val technology = seqRun.testType match {
                                case t if t.startsWith("BIGY") || t.contains("Y_ELITE") || t.contains("Y_PRIME") =>
                                  HaplogroupTechnology.BIG_Y
                                case _ => HaplogroupTechnology.WGS
                              }

                              // Create a RunHaplogroupCall for the reconciliation system
                              val runCall = RunHaplogroupCall(
                                sourceRef = seqRun.atUri.getOrElse(s"local:sequencerun:unknown"),
                                haplogroup = topResult.name,
                                confidence = topResult.score,
                                callMethod = CallMethod.SNP_PHYLOGENETIC,
                                score = Some(topResult.score),
                                supportingSnps = Some(topResult.matchingSnps),
                                conflictingSnps = Some(topResult.mismatchingSnps),
                                noCalls = None,
                                technology = Some(technology),
                                meanCoverage = None,
                                treeProvider = Some(treeProviderType.toString.toLowerCase),
                                treeVersion = None
                              )

                              // Convert TreeType to DnaType
                              val dnaType = treeType match {
                                case TreeType.YDNA => DnaType.Y_DNA
                                case TreeType.MTDNA => DnaType.MT_DNA
                              }

                              // Add to reconciliation - this automatically updates biosample haplogroups with consensus
                              val currentState = WorkspaceState(_workspace.value)
                              workspaceOps.addHaplogroupCall(currentState, subject.sampleAccession, dnaType, runCall) match {
                                case Right((newState, _)) =>
                                  // Persist updated biosample to H2
                                  newState.workspace.main.samples.find(_.sampleAccession == subject.sampleAccession).foreach { updatedBiosample =>
                                    h2Service.updateBiosample(updatedBiosample) match {
                                      case Right(persistedBs) =>
                                        log.info(s" Biosample haplogroup updated in H2: ${persistedBs.sampleAccession}")
                                      case Left(err) =>
                                        log.error(s"Failed to update Biosample haplogroup in H2: $err")
                                    }
                                  }
                                  _workspace.value = newState.workspace
                                case Left(hapError) =>
                                  log.error(s"Error adding haplogroup call: $hapError")
                              }

                              lastHaplogroupResult.value = Some(topResult)
                              analysisInProgress.value = false
                              analysisProgress.value = "Haplogroup analysis complete"
                              analysisProgressPercent.value = 1.0
                              onComplete(Right(topResult))
                            }

                          case Right(_) =>
                            Platform.runLater {
                              analysisInProgress.value = false
                              analysisError.value = "No haplogroup matches found"
                              analysisProgress.value = "Analysis complete - no matches"
                              onComplete(Left("No haplogroup matches found"))
                            }

                          case Left(error) =>
                            Platform.runLater {
                              analysisInProgress.value = false
                              analysisError.value = error
                              analysisProgress.value = s"Analysis failed: $error"
                              onComplete(Left(error))
                            }
                        }
                      } catch {
                        case e: Exception =>
                          Platform.runLater {
                            analysisInProgress.value = false
                            analysisError.value = e.getMessage
                            analysisProgress.value = s"Analysis failed: ${e.getMessage}"
                            onComplete(Left(e.getMessage))
                          }
                      }
                    }

                  case None =>
                    onComplete(Left("No alignment file associated with this sequence run"))
                }

              case None =>
                onComplete(Left("Please run initial analysis first to detect reference build"))
            }

          case None =>
            onComplete(Left(s"Sequence run not found at index $sequenceRunIndex for subject $sampleAccession"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /**
   * Gets the current haplogroup assignments for a subject.
   */
  def getHaplogroupAssignments(sampleAccession: String): Option[HaplogroupAssignments] = {
    findSubject(sampleAccession).flatMap(_.haplogroups)
  }

  /**
   * Gets the haplogroup artifact directory for a subject/run/alignment combination.
   * This can be used to display cached haplogroup reports.
   */
  def getHaplogroupArtifactDir(sampleAccession: String, sequenceRunIndex: Int): Option[java.nio.file.Path] = {
    for {
      subject <- findSubject(sampleAccession)
      sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
      seqRun <- sequenceRuns.lift(sequenceRunIndex)
      alignments = _workspace.value.main.getAlignmentsForSequenceRun(seqRun)
      alignment <- alignments.headOption
    } yield {
      val artifactCtx = ArtifactContext(
        sampleAccession = sampleAccession,
        sequenceRunUri = seqRun.atUri,
        alignmentUri = alignment.atUri
      )
      artifactCtx.getSubdir("haplogroup")
    }
  }

  /**
   * Gets the haplogroup artifact directory for a specific alignment.
   */
  def getHaplogroupArtifactDirForAlignment(sampleAccession: String, sequenceRunIndex: Int, alignmentIndex: Int): Option[java.nio.file.Path] = {
    for {
      subject <- findSubject(sampleAccession)
      sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
      seqRun <- sequenceRuns.lift(sequenceRunIndex)
      alignments = _workspace.value.main.getAlignmentsForSequenceRun(seqRun)
      alignment <- alignments.lift(alignmentIndex)
    } yield {
      val artifactCtx = ArtifactContext(
        sampleAccession = sampleAccession,
        sequenceRunUri = seqRun.atUri,
        alignmentUri = alignment.atUri
      )
      artifactCtx.getSubdir("haplogroup")
    }
  }

  /**
   * Runs haplogroup analysis for a specific alignment.
   */
  def runHaplogroupAnalysisForAlignment(
                                         sampleAccession: String,
                                         sequenceRunIndex: Int,
                                         alignmentIndex: Int,
                                         treeType: TreeType,
                                         onComplete: Either[String, AnalysisHaplogroupResult] => Unit
                                       ): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case Some(seqRun) =>
            val alignments = _workspace.value.main.getAlignmentsForSequenceRun(seqRun)
            alignments.lift(alignmentIndex) match {
              case Some(alignment) =>
                // Use alignment's file, fall back to seqRun file
                val fileInfo = alignment.files.headOption.orElse(seqRun.files.headOption)
                fileInfo match {
                  case Some(file) =>
                    val bamPath = file.location.getOrElse("")
                    log.info(s" Starting ${treeType} haplogroup analysis for ${file.fileName} (${alignment.referenceBuild})")

                    analysisInProgress.value = true
                    analysisError.value = ""
                    analysisProgress.value = "Starting haplogroup analysis..."
                    analysisProgressPercent.value = 0.0

                    Future {
                      try {
                        val processor = new HaplogroupProcessor()

                        // Select tree provider based on user preferences
                        val treeProviderType = treeType match {
                          case TreeType.YDNA =>
                            if (UserPreferencesService.getYdnaTreeProvider.equalsIgnoreCase("decodingus"))
                              TreeProviderType.DECODINGUS
                            else TreeProviderType.FTDNA
                          case TreeType.MTDNA =>
                            if (UserPreferencesService.getMtdnaTreeProvider.equalsIgnoreCase("decodingus"))
                              TreeProviderType.DECODINGUS
                            else TreeProviderType.FTDNA
                        }

                        // Check for existing whole-genome VCF
                        val runId = seqRun.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
                        val alignId = alignment.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
                        val vcfStatus = VcfCache.getStatus(sampleAccession, runId, alignId)

                        // Extract biosampleId for YProfile population
                        val biosampleId = subject.atUri.flatMap { uri =>
                          scala.util.Try(UUID.fromString(uri.split("/").last)).toOption
                        }

                        // Infer source type from platform/test type
                        val yProfileSourceType = {
                          val testType = seqRun.testType.toLowerCase
                          val platform = seqRun.platformName.toLowerCase
                          if (testType.contains("hifi") || testType.contains("pacbio") || platform.contains("pacbio")) {
                            YProfileSourceType.WGS_LONG_READ
                          } else if (testType.contains("nanopore") || platform.contains("nanopore")) {
                            YProfileSourceType.WGS_LONG_READ
                          } else if (testType.contains("targeted") || testType.contains("panel")) {
                            YProfileSourceType.TARGETED_NGS
                          } else {
                            YProfileSourceType.WGS_SHORT_READ
                          }
                        }

                        val result = if (vcfStatus.isAvailable) {
                          log.info(s" Found cached whole-genome VCF for $sampleAccession, using it for haplogroup analysis.")
                          processor.analyzeFromCachedVcf(
                            sampleAccession = sampleAccession,
                            runId = runId,
                            alignmentId = alignId,
                            referenceBuild = alignment.referenceBuild,
                            treeType = treeType,
                            treeProviderType = treeProviderType,
                            onProgress = (message, current, total) => {
                              Platform.runLater {
                                analysisProgress.value = message
                                analysisProgressPercent.value = if (total > 0) current / total else 0.0
                              }
                            },
                            yProfileService = yProfileService,
                            biosampleId = biosampleId,
                            yProfileSourceType = Some(yProfileSourceType)
                          )
                        } else {
                          // Check for existing haplogroup artifacts to avoid redundant analysis
                          val haplogroupDir = SubjectArtifactCache.getArtifactSubdir(sampleAccession, runId, alignId, "haplogroup")
                          val prefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"
                          val callsVcf = haplogroupDir.resolve(s"${prefix}_calls.vcf")

                          val artifactsExist = if (treeType == TreeType.YDNA) {
                            val privateVcf = haplogroupDir.resolve(s"${prefix}_private_variants.vcf")
                            Files.exists(callsVcf) && Files.exists(privateVcf)
                          } else {
                            Files.exists(callsVcf)
                          }

                          if (artifactsExist) {
                            log.info(s" Found existing haplogroup artifacts for $prefix at $haplogroupDir. Analysis will reuse them.")
                          }

                          val libraryStats = LibraryStats(
                            readCount = seqRun.totalReads.map(_.toInt).getOrElse(0),
                            pairedReads = 0,
                            lengthDistribution = Map.empty,
                            insertSizeDistribution = Map.empty,
                            aligner = alignment.aligner,
                            referenceBuild = alignment.referenceBuild,
                            sampleName = subject.donorIdentifier,
                            flowCells = Map.empty,
                            instruments = Map.empty,
                            mostFrequentInstrument = seqRun.instrumentModel.getOrElse("Unknown"),
                            inferredPlatform = seqRun.platformName,
                            platformCounts = Map.empty
                          )

                          val artifactCtx = ArtifactContext(
                            sampleAccession = sampleAccession,
                            sequenceRunUri = seqRun.atUri,
                            alignmentUri = alignment.atUri
                          )

                          processor.analyze(
                            bamPath = bamPath,
                            libraryStats = libraryStats,
                            treeType = treeType,
                            treeProviderType = treeProviderType,
                            onProgress = (message, current, total) => {
                              Platform.runLater {
                                analysisProgress.value = message
                                analysisProgressPercent.value = if (total > 0) current / total else 0.0
                              }
                            },
                            artifactContext = Some(artifactCtx),
                            yProfileService = yProfileService,
                            biosampleId = biosampleId,
                            yProfileSourceType = Some(yProfileSourceType)
                          )
                        }

                        Platform.runLater {
                          analysisInProgress.value = false
                          result match {
                            case Right(results) if results.nonEmpty =>
                              val topResult = results.head

                              // Update biosample haplogroups via reconciliation system
                              val technology = seqRun.testType match {
                                case t if t.startsWith("BIGY") || t.contains("Y_ELITE") || t.contains("Y_PRIME") =>
                                  HaplogroupTechnology.BIG_Y
                                case _ => HaplogroupTechnology.WGS
                              }

                              val treeProviderStr = treeProviderType.toString.toLowerCase
                              val runCall = RunHaplogroupCall(
                                sourceRef = seqRun.atUri.getOrElse(s"local:sequencerun:unknown"),
                                haplogroup = topResult.name,
                                confidence = topResult.score,
                                callMethod = CallMethod.SNP_PHYLOGENETIC,
                                score = Some(topResult.score),
                                supportingSnps = Some(topResult.matchingSnps),
                                conflictingSnps = Some(topResult.mismatchingSnps),
                                noCalls = None,
                                technology = Some(technology),
                                meanCoverage = None,
                                treeProvider = Some(treeProviderStr),
                                treeVersion = None
                              )

                              val dnaType = treeType match {
                                case TreeType.YDNA => DnaType.Y_DNA
                                case TreeType.MTDNA => DnaType.MT_DNA
                              }

                              // Add to reconciliation - this updates biosample haplogroups with consensus
                              val currentState = WorkspaceState(_workspace.value)
                              workspaceOps.addHaplogroupCall(currentState, subject.sampleAccession, dnaType, runCall) match {
                                case Right((newState, _)) =>
                                  // Persist updated biosample to H2
                                  newState.workspace.main.samples.find(_.sampleAccession == subject.sampleAccession).foreach { updatedBiosample =>
                                    h2Service.updateBiosample(updatedBiosample) match {
                                      case Right(persistedBs) =>
                                        log.info(s" Biosample haplogroup updated in H2: ${persistedBs.sampleAccession}")
                                      case Left(err) =>
                                        log.error(s"Failed to update Biosample haplogroup in H2: $err")
                                    }
                                  }
                                  _workspace.value = newState.workspace
                                case Left(hapError) =>
                                  log.error(s"Error adding haplogroup call: $hapError")
                              }

                              lastHaplogroupResult.value = Some(topResult)
                              onComplete(Right(topResult))
                            case Right(_) =>
                              onComplete(Left("No haplogroup results"))
                            case Left(error) =>
                              onComplete(Left(error))
                          }
                        }
                      } catch {
                        case e: Exception =>
                          Platform.runLater {
                            analysisInProgress.value = false
                            onComplete(Left(e.getMessage))
                          }
                      }
                    }

                  case None =>
                    onComplete(Left("No file associated with this alignment"))
                }
              case None =>
                onComplete(Left(s"Alignment not found at index $alignmentIndex"))
            }
          case None =>
            onComplete(Left(s"Sequence run not found at index $sequenceRunIndex"))
        }
      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /**
   * Runs callable loci analysis for a specific alignment.
   */
  def runCallableLociAnalysisForAlignment(
                                           sampleAccession: String,
                                           sequenceRunIndex: Int,
                                           alignmentIndex: Int,
                                           onComplete: Either[String, (CallableLociResult, java.nio.file.Path)] => Unit
                                         ): Unit = {
    // Delegate to existing method for now - alignment index support to be added later
    runCallableLociAnalysis(sampleAccession, sequenceRunIndex, onComplete)
  }

  /**
   * Runs WGS metrics analysis for a specific alignment.
   */
  def runWgsMetricsAnalysisForAlignment(
                                         sampleAccession: String,
                                         sequenceRunIndex: Int,
                                         alignmentIndex: Int,
                                         onComplete: Either[String, WgsMetrics] => Unit
                                       ): Unit = {
    // Delegate to existing method for now - alignment index support to be added later
    runWgsMetricsAnalysis(sampleAccession, sequenceRunIndex, (result: Either[String, WgsMetrics]) => {
      onComplete(result)
    })
  }

  /**
   * Runs multiple metrics analysis for a specific alignment.
   */
  def runMultipleMetricsAnalysisForAlignment(
                                              sampleAccession: String,
                                              sequenceRunIndex: Int,
                                              alignmentIndex: Int,
                                              onComplete: Either[String, com.decodingus.analysis.ReadMetrics] => Unit
                                            ): Unit = {
    // Delegate to existing method for now - alignment index support to be added later
    runMultipleMetricsAnalysis(sampleAccession, sequenceRunIndex, onComplete)
  }

  /**
   * Runs whole-genome variant calling for a specific alignment.
   * This is a long-running operation that generates a cached VCF.
   */
  def runWholeGenomeVariantCallingForAlignment(
                                                sampleAccession: String,
                                                sequenceRunIndex: Int,
                                                alignmentIndex: Int,
                                                onComplete: Either[String, CachedVcfInfo] => Unit
                                              ): Unit = {
    analysisInProgress.value = true
    analysisProgress.value = "Starting whole-genome variant calling..."
    analysisProgressPercent.value = 0.0

    val progressHandler: AnalysisProgress => Unit = { progress =>
      Platform.runLater {
        analysisProgress.value = progress.message
        analysisProgressPercent.value = progress.percent
      }
    }

    analysisCoordinator.runWholeGenomeVariantCalling(
      currentState,
      sampleAccession,
      sequenceRunIndex,
      alignmentIndex,
      progressHandler
    ).onComplete {
      case Success(Right((newState, vcfInfo))) =>
        Platform.runLater {
          applyState(newState)
          analysisInProgress.value = false
          analysisProgress.value = "Whole-genome VCF generation complete"
          analysisProgressPercent.value = 1.0
          onComplete(Right(vcfInfo))
        }
      case Success(Left(error)) =>
        Platform.runLater {
          analysisInProgress.value = false
          analysisProgress.value = s"Error: $error"
          onComplete(Left(error))
        }
      case Failure(ex) =>
        Platform.runLater {
          analysisInProgress.value = false
          analysisProgress.value = s"Error: ${ex.getMessage}"
          onComplete(Left(ex.getMessage))
        }
    }
  }

  /**
   * Runs comprehensive analysis pipeline for a specific alignment.
   * Executes: Read Metrics → WGS Metrics → Callable Loci → Sex Inference → mtDNA → Y-DNA → Ancestry stub
   */
  def runComprehensiveAnalysisForAlignment(
                                            sampleAccession: String,
                                            sequenceRunIndex: Int,
                                            alignmentIndex: Int,
                                            onComplete: Either[String, BatchAnalysisResult] => Unit
                                          ): Unit = {
    analysisInProgress.value = true
    analysisProgress.value = "Starting comprehensive analysis..."
    analysisProgressPercent.value = 0.0

    val progressHandler: AnalysisProgress => Unit = { progress =>
      Platform.runLater {
        analysisProgress.value = progress.message
        analysisProgressPercent.value = progress.percent
      }
    }

    analysisCoordinator.runComprehensiveAnalysis(
      currentState,
      sampleAccession,
      sequenceRunIndex,
      alignmentIndex,
      progressHandler
    ).onComplete {
      case Success(Right((newState, batchResult))) =>
        Platform.runLater {
          applyState(newState)
          analysisInProgress.value = false
          analysisProgress.value = "Comprehensive analysis complete"
          analysisProgressPercent.value = 1.0
          onComplete(Right(batchResult))
        }
      case Success(Left(error)) =>
        Platform.runLater {
          analysisInProgress.value = false
          analysisProgress.value = s"Error: $error"
          onComplete(Left(error))
        }
      case Failure(ex) =>
        Platform.runLater {
          analysisInProgress.value = false
          analysisProgress.value = s"Error: ${ex.getMessage}"
          onComplete(Left(ex.getMessage))
        }
    }
  }

  // --- STR Profile CRUD Operations ---

  /**
   * Adds a new STR profile for a biosample.
   * Supports multiple profiles per subject (e.g., from different vendors like FTDNA and YSEQ).
   * Returns the URI of the new profile, or an error message.
   */
  def addStrProfile(sampleAccession: String, profile: StrProfile): Either[String, String] = {
    // Get biosample ID for H2 persistence
    getBiosampleIdByAccession(sampleAccession) match {
      case None =>
        Left(s"Biosample not found: $sampleAccession")
      case Some(biosampleId) =>
        // Persist to H2 first
        h2Service.createStrProfile(profile, biosampleId) match {
          case Right(savedProfile) =>
            // Update in-memory state with the saved profile (which has the atUri)
            workspaceOps.addStrProfile(currentState, sampleAccession, savedProfile) match {
              case Right((newState, profileUri)) =>
                applyState(newState)
                log.info(s"Added STR profile with ${profile.markers.size} markers for $sampleAccession (provider: ${profile.importedFrom.getOrElse("unknown")})")
                Right(profileUri)
              case Left(error) =>
                // H2 succeeded but in-memory failed - this shouldn't happen but log it
                log.warn(s"STR profile saved to H2 but in-memory update failed: $error")
                savedProfile.atUri.toRight(error)
            }
          case Left(error) =>
            log.error(s"Failed to persist STR profile to H2: $error")
            Left(error)
        }
    }
  }

  /**
   * Gets all STR profiles for a biosample.
   * Returns profiles from all vendors (FTDNA, YSEQ, WGS-derived, etc.)
   */
  def getStrProfilesForBiosample(sampleAccession: String): List[StrProfile] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Get profiles by refs list (preferred) or by biosampleRef matching
        val byRefs = subject.strProfileRefs.flatMap { ref =>
          _workspace.value.main.strProfiles.find(_.atUri.contains(ref))
        }
        if (byRefs.nonEmpty) byRefs
        else {
          // Fallback: find by biosampleRef for legacy data
          val biosampleUri = subject.atUri.getOrElse(s"local:biosample:$sampleAccession")
          _workspace.value.main.strProfiles.filter(_.biosampleRef == biosampleUri)
        }
      case None =>
        List.empty
    }
  }

  /**
   * Gets all STR profiles in the workspace.
   */
  def getAllStrProfiles: List[StrProfile] = {
    _workspace.value.main.strProfiles
  }

  /**
   * Updates an existing STR profile.
   */
  def updateStrProfile(profileUri: String, updatedProfile: StrProfile): Either[String, Unit] = {
    workspaceOps.updateStrProfile(currentState, profileUri, updatedProfile) match {
      case Right(newState) =>
        applyState(newState)
        Right(())
      case Left(error) =>
        Left(error)
    }
  }

  /**
   * Deletes an STR profile.
   */
  def deleteStrProfile(sampleAccession: String, profileUri: String): Either[String, Unit] = {
    workspaceOps.deleteStrProfile(currentState, sampleAccession, profileUri) match {
      case Right(newState) =>
        applyState(newState)
        log.info(s" Deleted STR profile $profileUri for $sampleAccession")
        Right(())
      case Left(error) =>
        Left(error)
    }
  }

  // --- Chip Profile CRUD Operations ---

  /**
   * Imports chip data from a file for a biosample.
   * Parses the file, computes statistics, and creates a ChipProfile record.
   *
   * @param sampleAccession The subject's accession ID
   * @param file            The chip data file to import
   * @param onComplete      Callback when import completes
   */
  def importChipData(
                      sampleAccession: String,
                      file: File,
                      onComplete: Either[String, ChipProfile] => Unit
                    ): Unit = {
    import com.decodingus.genotype.model.GenotypingTestSummary
    import com.decodingus.genotype.parser.ChipDataParser

    import java.security.MessageDigest

    findSubject(sampleAccession) match {
      case Some(subject) =>
        analysisInProgress.value = true
        analysisError.value = ""
        analysisProgress.value = "Detecting chip format..."
        analysisProgressPercent.value = 0.1

        Future {
          try {
            // Step 1: Detect and parse the file
            updateProgress("Parsing chip data file...", 0.2)

            ChipDataParser.detectParser(file) match {
              case Right((parser, detection)) =>
                updateProgress(s"Detected ${parser.vendor} format...", 0.3)

                parser.parse(file, (current, total) => {
                  val pct = 0.3 + (current.toDouble / total) * 0.4
                  updateProgress(s"Reading genotypes: ${current}/${total}", pct)
                }) match {
                  case Right(callsIterator) =>
                    updateProgress("Computing statistics...", 0.75)

                    val calls = callsIterator.toList
                    val summary = GenotypingTestSummary.fromCalls(
                      calls,
                      detection.testType.getOrElse(com.decodingus.genotype.model.TestTypes.ARRAY_23ANDME_V5),
                      detection.chipVersion,
                      Some(computeFileHash(file))
                    )

                    // Create ChipProfile
                    updateProgress("Creating chip profile...", 0.9)

                    Platform.runLater {
                      val fileInfo = FileInfo(
                        fileName = file.getName,
                        fileSizeBytes = Some(file.length()),
                        fileFormat = if (file.getName.endsWith(".csv")) "CSV" else "TXT",
                        checksum = Some(computeFileHash(file)),
                        checksumAlgorithm = Some("SHA-256"),
                        location = Some(file.getAbsolutePath)
                      )

                      val chipProfile = ChipProfile(
                        atUri = None, // Will be assigned by H2
                        meta = RecordMeta.initial,
                        biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}"),
                        vendor = parser.vendor,
                        testTypeCode = detection.testType.map(_.code).getOrElse("UNKNOWN"),
                        chipVersion = detection.chipVersion,
                        totalMarkersCalled = summary.totalMarkersCalled,
                        totalMarkersPossible = summary.totalMarkersPossible,
                        noCallRate = summary.noCallRate,
                        yMarkersCalled = summary.yMarkersCalled,
                        mtMarkersCalled = summary.mtMarkersCalled,
                        autosomalMarkersCalled = summary.autosomalMarkersCalled,
                        hetRate = summary.hetRate,
                        importDate = LocalDateTime.now(),
                        sourceFileHash = Some(computeFileHash(file)),
                        sourceFileName = Some(file.getName),
                        files = List(fileInfo)
                      )

                      // Get biosample ID for H2 persistence
                      getBiosampleIdByAccession(sampleAccession) match {
                        case None =>
                          analysisInProgress.value = false
                          analysisError.value = s"Biosample not found: $sampleAccession"
                          onComplete(Left(s"Biosample not found: $sampleAccession"))

                        case Some(biosampleId) =>
                          // Persist chip profile to H2
                          h2Service.createChipProfile(chipProfile, biosampleId) match {
                            case Right(savedProfile) =>
                              val chipProfileUri = savedProfile.atUri.getOrElse(s"local:chipprofile:${java.util.UUID.randomUUID()}")

                              // Update workspace with saved profile
                              val updatedChipProfiles = _workspace.value.main.chipProfiles :+ savedProfile
                              val updatedSubject = subject.copy(
                                genotypeRefs = subject.genotypeRefs :+ chipProfileUri,
                                meta = subject.meta.updated("genotypeRefs")
                              )

                              // Persist biosample update to H2
                              h2Service.updateBiosample(updatedSubject) match {
                                case Right(persisted) =>
                                  log.info(s"Biosample genotypeRefs updated in H2: ${persisted.sampleAccession}")
                                case Left(error) =>
                                  log.error(s"Failed to update Biosample genotypeRefs in H2: $error")
                              }

                              val updatedSamples = _workspace.value.main.samples.map { s =>
                                if (s.sampleAccession == sampleAccession) updatedSubject else s
                              }
                              val updatedContent = _workspace.value.main.copy(
                                samples = updatedSamples,
                                chipProfiles = updatedChipProfiles
                              )
                              _workspace.value = _workspace.value.copy(main = updatedContent)

                              analysisInProgress.value = false
                              analysisProgress.value = "Import complete"
                              analysisProgressPercent.value = 1.0

                              log.info(s"Imported chip data: ${parser.vendor}, ${summary.totalMarkersCalled} markers")
                              onComplete(Right(savedProfile))

                            case Left(error) =>
                              log.error(s"Failed to persist chip profile to H2: $error")
                              analysisInProgress.value = false
                              analysisError.value = error
                              onComplete(Left(error))
                          }
                      }
                    }

                  case Left(parseError) =>
                    Platform.runLater {
                      analysisInProgress.value = false
                      analysisError.value = parseError
                      onComplete(Left(parseError))
                    }
                }

              case Left(detectError) =>
                Platform.runLater {
                  analysisInProgress.value = false
                  analysisError.value = detectError
                  onComplete(Left(detectError))
                }
            }

          } catch {
            case e: Exception =>
              Platform.runLater {
                analysisInProgress.value = false
                analysisError.value = e.getMessage
                onComplete(Left(e.getMessage))
              }
          }
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /** Compute SHA-256 hash of a file */
  private def computeFileHash(file: File): String = {
    import java.nio.file.Files
    import java.security.MessageDigest
    val bytes = Files.readAllBytes(file.toPath)
    val digest = MessageDigest.getInstance("SHA-256").digest(bytes)
    digest.map("%02x".format(_)).mkString
  }

  /**
   * Gets all chip profiles for a biosample.
   */
  def getChipProfilesForBiosample(sampleAccession: String): List[ChipProfile] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Get profiles by genotypeRefs list or by biosampleRef matching
        val byRefs = subject.genotypeRefs.flatMap { ref =>
          _workspace.value.main.chipProfiles.find(_.atUri.contains(ref))
        }
        if (byRefs.nonEmpty) byRefs
        else {
          // Fallback: find by biosampleRef for legacy data
          val biosampleUri = subject.atUri.getOrElse(s"local:biosample:$sampleAccession")
          _workspace.value.main.chipProfiles.filter(_.biosampleRef == biosampleUri)
        }
      case None =>
        List.empty
    }
  }

  /**
   * Deletes a chip profile.
   */
  def deleteChipProfile(sampleAccession: String, profileUri: String): Either[String, Unit] = {
    workspaceOps.deleteChipProfile(currentState, sampleAccession, profileUri) match {
      case Right(newState) =>
        applyState(newState)
        log.info(s" Deleted chip profile $profileUri for $sampleAccession")
        Right(())
      case Left(error) =>
        Left(error)
    }
  }

  /**
   * Runs ancestry analysis on chip/array data for a biosample.
   *
   * This uses the ChipAncestryAdapter to project chip genotypes onto PCA space
   * and estimate population proportions. Unlike WGS ancestry analysis, chip data
   * already has genotypes called at known positions so we can directly project
   * without calling GATK HaplotypeCaller.
   *
   * @param sampleAccession The subject's accession ID
   * @param profileUri      The AT URI of the chip profile to analyze
   * @param panelType       AIMs (quick) or GenomeWide (detailed)
   * @param onComplete      Callback when analysis completes
   */
  def runChipAncestryAnalysis(
                               sampleAccession: String,
                               profileUri: String,
                               panelType: com.decodingus.ancestry.model.AncestryPanelType,
                               onComplete: Either[String, com.decodingus.ancestry.model.AncestryResult] => Unit
                             ): Unit = {
    import com.decodingus.genotype.parser.ChipDataParser
    import com.decodingus.genotype.processor.{ChipAncestryAdapter, ChipDataProcessor}

    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Find the chip profile
        _workspace.value.main.chipProfiles.find(_.atUri.contains(profileUri)) match {
          case Some(profile) =>
            // Need the source file to re-parse genotypes
            profile.files.headOption.flatMap(_.location) match {
              case Some(filePath) =>
                val file = new File(filePath)
                if (!file.exists()) {
                  onComplete(Left(s"Source file not found: $filePath"))
                } else {
                  analysisInProgress.value = true
                  analysisError.value = ""
                  analysisProgress.value = "Loading chip genotypes..."
                  analysisProgressPercent.value = 0.1

                  Future {
                    try {
                      // Re-parse the chip data to get genotypes
                      val processor = new ChipDataProcessor()

                      processor.process(file, (msg, done, total) => {
                        updateProgress(msg, done)
                      }) match {
                        case Right(chipResult) =>
                          updateProgress("Running ancestry analysis...", 0.5)

                          // Run ancestry analysis
                          val adapter = new ChipAncestryAdapter()
                          adapter.analyze(chipResult, panelType, (msg, done, total) => {
                            updateProgress(msg, 0.5 + done * 0.4)
                          }) match {
                            case Right(ancestryResult) =>
                              Platform.runLater {
                                updateProgress("Ancestry analysis complete.", 1.0)
                                analysisInProgress.value = false

                                // Store the result in the workspace
                                val updatedSubject = subject.copy(
                                  populationBreakdownRef = Some(s"local:ancestry:$sampleAccession:${java.util.UUID.randomUUID().toString.take(8)}"),
                                  meta = subject.meta.updated("populationBreakdownRef")
                                )

                                // Persist biosample update to H2
                                h2Service.updateBiosample(updatedSubject) match {
                                  case Right(persisted) =>
                                    log.info(s" Biosample ancestry ref updated in H2: ${persisted.sampleAccession}")
                                  case Left(error) =>
                                    log.error(s"Failed to update Biosample ancestry ref in H2: $error")
                                }

                                val updatedSamples = _workspace.value.main.samples.map { s =>
                                  if (s.sampleAccession == sampleAccession) updatedSubject else s
                                }
                                val updatedContent = _workspace.value.main.copy(samples = updatedSamples)
                                _workspace.value = _workspace.value.copy(main = updatedContent)

                                log.info(s" Ancestry analysis complete for $sampleAccession: " +
                                  s"${ancestryResult.percentages.take(3).map(p => s"${p.populationName}: ${f"${p.percentage}%.1f"}%").mkString(", ")}")
                                onComplete(Right(ancestryResult))
                              }

                            case Left(analysisError) =>
                              Platform.runLater {
                                this.analysisInProgress.value = false
                                this.analysisError.value = analysisError
                                onComplete(Left(analysisError))
                              }
                          }

                        case Left(parseError) =>
                          Platform.runLater {
                            analysisInProgress.value = false
                            analysisError.value = parseError
                            onComplete(Left(parseError))
                          }
                      }

                    } catch {
                      case e: Exception =>
                        Platform.runLater {
                          analysisInProgress.value = false
                          analysisError.value = e.getMessage
                          onComplete(Left(e.getMessage))
                        }
                    }
                  }
                }

              case None =>
                onComplete(Left("Chip profile has no source file location. Re-import the chip data."))
            }

          case None =>
            onComplete(Left(s"Chip profile not found: $profileUri"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /**
   * Runs haplogroup analysis on chip/array data for a biosample.
   *
   * This uses the ChipHaplogroupAdapter to score chip genotypes against the
   * haplogroup tree (Y-DNA or mtDNA). Unlike WGS-based analysis, chip data
   * has limited coverage of tree positions (typically 10-30%), so results
   * may be less precise.
   *
   * @param sampleAccession The subject's accession ID
   * @param profileUri      The AT URI of the chip profile to analyze
   * @param treeType        Y-DNA or MT-DNA tree type
   * @param onComplete      Callback when analysis completes
   */
  def runChipHaplogroupAnalysis(
                                 sampleAccession: String,
                                 profileUri: String,
                                 treeType: com.decodingus.haplogroup.tree.TreeType,
                                 onComplete: Either[String, com.decodingus.genotype.processor.ChipHaplogroupResult] => Unit
                               ): Unit = {
    import com.decodingus.genotype.processor.{ChipDataProcessor, ChipHaplogroupAdapter}

    val treeLabel = if (treeType == com.decodingus.haplogroup.tree.TreeType.YDNA) "Y-DNA" else "mtDNA"

    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Find the chip profile
        _workspace.value.main.chipProfiles.find(_.atUri.contains(profileUri)) match {
          case Some(profile) =>
            // Need the source file to re-parse genotypes
            profile.files.headOption.flatMap(_.location) match {
              case Some(filePath) =>
                val file = new File(filePath)
                if (!file.exists()) {
                  onComplete(Left(s"Source file not found: $filePath"))
                } else {
                  analysisInProgress.value = true
                  analysisError.value = ""
                  analysisProgress.value = s"Loading chip data for $treeLabel analysis..."
                  analysisProgressPercent.value = 0.1

                  Future {
                    try {
                      // Re-parse the chip data to get genotypes
                      val processor = new ChipDataProcessor()

                      processor.process(file, (msg, done, total) => {
                        updateProgress(msg, done * 0.3)
                      }) match {
                        case Right(chipResult) =>
                          updateProgress(s"Running $treeLabel haplogroup analysis...", 0.3)

                          // Run haplogroup analysis
                          val adapter = new ChipHaplogroupAdapter()
                          adapter.analyze(chipResult, treeType, (msg, done, total) => {
                            updateProgress(msg, 0.3 + done * 0.6)
                          }) match {
                            case Right(haplogroupResult) =>
                              Platform.runLater {
                                updateProgress(s"$treeLabel haplogroup analysis complete.", 1.0)
                                analysisInProgress.value = false

                                // Determine tree provider used (same logic as ChipHaplogroupAdapter)
                                val treeProviderName = treeType match {
                                  case com.decodingus.haplogroup.tree.TreeType.YDNA =>
                                    if (com.decodingus.config.UserPreferencesService.getYdnaTreeProvider.equalsIgnoreCase("decodingus")) "decodingus"
                                    else "ftdna"
                                  case com.decodingus.haplogroup.tree.TreeType.MTDNA =>
                                    if (com.decodingus.config.UserPreferencesService.getMtdnaTreeProvider.equalsIgnoreCase("decodingus")) "decodingus"
                                    else "ftdna"
                                }

                                // Create a RunHaplogroupCall for the reconciliation system
                                val runCall = RunHaplogroupCall(
                                  sourceRef = profile.atUri.getOrElse(s"local:chipprofile:unknown"),
                                  haplogroup = haplogroupResult.topHaplogroup,
                                  confidence = haplogroupResult.confidence,
                                  callMethod = CallMethod.SNP_PHYLOGENETIC,
                                  score = haplogroupResult.results.headOption.map(_.score),
                                  supportingSnps = Some(haplogroupResult.snpsMatched),
                                  conflictingSnps = None,
                                  noCalls = Some(haplogroupResult.snpsTotal - haplogroupResult.snpsMatched),
                                  technology = Some(HaplogroupTechnology.SNP_ARRAY),
                                  meanCoverage = None,
                                  treeProvider = Some(treeProviderName),
                                  treeVersion = None
                                )

                                // Convert TreeType to DnaType
                                val dnaType = treeType match {
                                  case com.decodingus.haplogroup.tree.TreeType.YDNA => DnaType.Y_DNA
                                  case com.decodingus.haplogroup.tree.TreeType.MTDNA => DnaType.MT_DNA
                                }

                                // Add to reconciliation - this automatically updates biosample haplogroups with consensus
                                val currentState = WorkspaceState(_workspace.value)
                                workspaceOps.addHaplogroupCall(currentState, sampleAccession, dnaType, runCall) match {
                                  case Right((newState, _)) =>
                                    // Persist updated biosample to H2
                                    newState.workspace.main.samples.find(_.sampleAccession == sampleAccession).foreach { updatedBiosample =>
                                      h2Service.updateBiosample(updatedBiosample) match {
                                        case Right(persistedBs) =>
                                          log.info(s" Biosample chip haplogroup updated in H2: ${persistedBs.sampleAccession}")
                                        case Left(err) =>
                                          log.error(s"Failed to update Biosample chip haplogroup in H2: $err")
                                      }
                                    }
                                    _workspace.value = newState.workspace
                                  case Left(hapError) =>
                                    log.error(s"Error adding chip haplogroup call: $hapError")
                                }

                                log.info(s" $treeLabel haplogroup analysis complete for $sampleAccession: " +
                                  s"${haplogroupResult.topHaplogroup} (${haplogroupResult.snpsMatched}/${haplogroupResult.snpsTotal} SNPs)")
                                onComplete(Right(haplogroupResult))
                              }

                            case Left(analysisErr) =>
                              Platform.runLater {
                                this.analysisInProgress.value = false
                                this.analysisError.value = analysisErr
                                onComplete(Left(analysisErr))
                              }
                          }

                        case Left(parseError) =>
                          Platform.runLater {
                            analysisInProgress.value = false
                            analysisError.value = parseError
                            onComplete(Left(parseError))
                          }
                      }

                    } catch {
                      case e: Exception =>
                        Platform.runLater {
                          analysisInProgress.value = false
                          analysisError.value = e.getMessage
                          onComplete(Left(e.getMessage))
                        }
                    }
                  }
                }

              case None =>
                onComplete(Left("Chip profile has no source file location. Re-import the chip data."))
            }

          case None =>
            onComplete(Left(s"Chip profile not found: $profileUri"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  /**
   * Runs CollectMultipleMetrics (alignment summary + insert size) on a sequencing run.
   *
   * This collects metrics not provided by CollectWgsMetrics:
   * - Total reads, PF reads, aligned reads
   * - Pair rates and proper pair percentages
   * - Insert size distribution (median, mean, std)
   *
   * Results are stored back into the SequenceRun model.
   *
   * @param sampleAccession  The subject's accession ID
   * @param sequenceRunIndex Index of the sequence run within the subject
   * @param onComplete       Callback when analysis completes
   */
  def runMultipleMetricsAnalysis(
                                  sampleAccession: String,
                                  sequenceRunIndex: Int,
                                  onComplete: Either[String, com.decodingus.analysis.ReadMetrics] => Unit
                                ): Unit = {
    import com.decodingus.analysis.{ArtifactContext, UnifiedMetricsProcessor}

    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.length) {
          onComplete(Left(s"Invalid sequence run index: $sequenceRunIndex"))
          return
        }

        val sequenceRun = sequenceRuns(sequenceRunIndex)
        val alignments = _workspace.value.main.getAlignmentsForSequenceRun(sequenceRun)

        if (alignments.isEmpty) {
          onComplete(Left("No alignments found for this sequence run"))
          return
        }

        val alignment = alignments.head
        val bamFile = alignment.files.find(f =>
          f.fileFormat == "BAM" || f.fileFormat == "CRAM"
        ).flatMap(_.location)

        bamFile match {
          case Some(bamPath) =>
            // Verify file exists
            if (!new File(bamPath).exists()) {
              onComplete(Left(s"BAM/CRAM file not found: $bamPath"))
              return
            }

            // Get reference path
            val referenceBuild = alignment.referenceBuild
            val referenceGateway = new com.decodingus.refgenome.ReferenceGateway((_, _) => {})

            referenceGateway.resolve(referenceBuild) match {
              case Left(error) =>
                onComplete(Left(s"Could not resolve reference genome: $error"))

              case Right(referencePath) =>
                analysisInProgress.value = true
                analysisError.value = ""
                analysisProgress.value = "Collecting read metrics..."
                analysisProgressPercent.value = 0.0

                // Create artifact context
                val artifactContext = ArtifactContext(
                  sampleAccession = sampleAccession,
                  sequenceRunUri = sequenceRun.atUri,
                  alignmentUri = alignment.atUri
                )

                Future {
                  try {
                    // Use UnifiedMetricsProcessor - no R dependency
                    val processor = new UnifiedMetricsProcessor()
                    processor.process(
                      bamPath = bamPath,
                      referencePath = referencePath.toString,
                      onProgress = (msg, done, total) => updateProgress(msg, done),
                      artifactContext = Some(artifactContext)
                    ) match {
                      case Right(metricsResult) =>
                        Platform.runLater {
                          analysisInProgress.value = false

                          // Update the sequence run with the new metrics
                          val updatedRun = sequenceRun.copy(
                            totalReads = Some(metricsResult.totalReads),
                            pfReads = Some(metricsResult.pfReads),
                            pfReadsAligned = Some(metricsResult.pfReadsAligned),
                            pctPfReadsAligned = Some(metricsResult.pctPfReadsAligned),
                            readsPaired = Some(metricsResult.readsAlignedInPairs),
                            pctReadsPaired = Some(metricsResult.pctReadsAlignedInPairs),
                            pctProperPairs = Some(metricsResult.pctProperPairs),
                            readLength = Some(metricsResult.meanReadLength.toInt),
                            maxReadLength = Some(metricsResult.maxReadLength),
                            meanInsertSize = Some(metricsResult.meanInsertSize),
                            medianInsertSize = Some(metricsResult.medianInsertSize),
                            stdInsertSize = Some(metricsResult.stdInsertSize),
                            pairOrientation = Some(metricsResult.pairOrientation),
                            meta = sequenceRun.meta.updated("readMetrics")
                          )

                          // Use the existing updateSequenceRun method
                          updateSequenceRun(sampleAccession, sequenceRunIndex, updatedRun)

                          log.info(s" ReadMetrics complete for $sampleAccession: " +
                            s"${metricsResult.totalReads} reads, ${f"${metricsResult.pctPfReadsAligned * 100}%.1f"}% aligned, " +
                            s"median insert ${metricsResult.medianInsertSize.toInt}bp")
                          onComplete(Right(metricsResult))
                        }

                      case Left(processorError) =>
                        Platform.runLater {
                          analysisInProgress.value = false
                          analysisError.value = processorError.getMessage
                          onComplete(Left(processorError.getMessage))
                        }
                    }
                  } catch {
                    case e: Exception =>
                      Platform.runLater {
                        analysisInProgress.value = false
                        analysisError.value = e.getMessage
                        onComplete(Left(e.getMessage))
                      }
                  }
                }
            }

          case None =>
            onComplete(Left("No BAM/CRAM file found in alignment"))
        }

      case None =>
        onComplete(Left(s"Subject not found: $sampleAccession"))
    }
  }

  // ============================================
  // Y Chromosome Profile Methods
  // ============================================

  /**
   * Summary data for Y Profile display in subject detail view.
   */
  case class YProfileSummary(
                              profileId: UUID,
                              consensusHaplogroup: Option[String],
                              haplogroupConfidence: Option[Double],
                              totalVariants: Int,
                              confirmedCount: Int,
                              novelCount: Int,
                              conflictCount: Int,
                              sourceCount: Int,
                              callableRegionPct: Option[Double]
                            )

  /**
   * Full Y Profile data for detail dialog.
   *
   * @param yRegionAnnotator Optional annotator for region visualization (ideogram)
   */
  case class YProfileLoadedData(
                                 profile: YChromosomeProfileEntity,
                                 variants: List[YProfileVariantEntity],
                                 sources: List[YProfileSourceEntity],
                                 variantCalls: Map[UUID, List[YVariantSourceCallEntity]],
                                 auditEntries: List[YVariantAuditEntity],
                                 yRegionAnnotator: Option[YRegionAnnotator] = None
                               )

  /**
   * Get Y Profile summary for a biosample (lightweight query).
   * Returns None if no profile exists or service unavailable.
   */
  def getYProfileSummary(biosampleId: UUID): Option[YProfileSummary] =
    yProfileService.flatMap { service =>
      service.getProfileByBiosample(biosampleId).toOption.flatten.map { profile =>
        YProfileSummary(
          profileId = profile.id,
          consensusHaplogroup = profile.consensusHaplogroup,
          haplogroupConfidence = profile.haplogroupConfidence,
          totalVariants = profile.totalVariants,
          confirmedCount = profile.confirmedCount,
          novelCount = profile.novelCount,
          conflictCount = profile.conflictCount,
          sourceCount = profile.sourceCount,
          callableRegionPct = profile.callableRegionPct
        )
      }
    }

  /**
   * Check if Y Profile service is available.
   */
  def isYProfileAvailable: Boolean = yProfileService.isDefined

  /**
   * Load full Y Profile data for a biosample (async, on-demand).
   *
   * @param biosampleId The biosample UUID
   * @param onComplete  Callback with Either[error, data]
   */
  def loadYProfileForBiosample(
                                biosampleId: UUID,
                                onComplete: Either[String, YProfileLoadedData] => Unit
                              ): Unit = {
    yProfileService match {
      case None =>
        onComplete(Left("Y Profile service not available"))

      case Some(service) =>
        Future {
          for {
            profileOpt <- service.getProfileByBiosample(biosampleId)
            profile <- profileOpt.toRight("No Y Profile found for this biosample")
            variants <- service.getVariants(profile.id)
            sources <- service.getSources(profile.id)
            // Load variant calls for each variant
            variantCalls = variants.flatMap { v =>
              service.getVariantCalls(v.id).toOption.map(calls => v.id -> calls)
            }.toMap
            // Load audit entries for all variants
            auditEntries = variants.flatMap { v =>
              service.getAuditHistory(v.id).toOption.getOrElse(Nil)
            }
            // Get reference build from first source (or default to GRCh38)
            referenceBuild = sources.flatMap(_.referenceBuild).headOption
            // Create basic annotator based on reference build (with hardcoded heterochromatin)
            annotator = createBasicYRegionAnnotator(referenceBuild)
          } yield YProfileLoadedData(profile, variants, sources, variantCalls, auditEntries, Some(annotator))
        }.onComplete {
          case Success(result) =>
            Platform.runLater {
              onComplete(result)
            }
          case Failure(ex) =>
            Platform.runLater {
              onComplete(Left(ex.getMessage))
            }
        }
    }
  }

  /**
   * Data loaded for Y Profile management dialog.
   */
  case class YProfileManagementData(
                                     yProfileService: YProfileService,
                                     profile: Option[YChromosomeProfileEntity],
                                     sources: List[YProfileSourceEntity],
                                     variants: List[YProfileVariantEntity],
                                     snpPanels: List[YSnpPanelEntity]
                                   )

  /**
   * Load Y Profile data for the management dialog.
   * Unlike loadYProfileForBiosample, this works even when no profile exists yet.
   */
  def loadYProfileManagementData(
                                  biosampleId: UUID,
                                  onComplete: Either[String, YProfileManagementData] => Unit
                                ): Unit = {
    yProfileService match {
      case None =>
        onComplete(Left("Y Profile service not available"))

      case Some(service) =>
        Future {
          val snpPanelRepo = new YSnpPanelRepository()
          val transactor = databaseContext.transactor

          // Load profile (may not exist)
          val profileOpt = service.getProfileByBiosample(biosampleId).toOption.flatten

          // Load sources and variants (if profile exists)
          val (sources, variants) = profileOpt match {
            case Some(profile) =>
              val s = service.getSources(profile.id).toOption.getOrElse(Nil)
              val v = service.getVariants(profile.id).toOption.getOrElse(Nil)
              (s, v)
            case None =>
              (Nil, Nil)
          }

          // Load SNP panels for this biosample
          val snpPanels = transactor.readOnly {
            snpPanelRepo.findByBiosample(biosampleId)
          }.toOption.getOrElse(Nil)

          Right(YProfileManagementData(service, profileOpt, sources, variants, snpPanels))
        }.onComplete {
          case Success(result) =>
            Platform.runLater {
              onComplete(result)
            }
          case Failure(ex) =>
            Platform.runLater {
              onComplete(Left(ex.getMessage))
            }
        }
    }
  }

  /**
   * Create a basic Y region annotator with hardcoded heterochromatin boundaries.
   * Uses reference build to select correct coordinates.
   */
  private def createBasicYRegionAnnotator(referenceBuild: Option[String]): YRegionAnnotator = {
    val heterochromatin = referenceBuild.map(_.toUpperCase) match {
      case Some(b) if b.contains("37") || b.contains("HG19") =>
        YRegionAnnotator.grch37Heterochromatin
      case Some(b) if b.contains("CHM13") || b.contains("T2T") || b.contains("HS1") =>
        YRegionAnnotator.chm13v2Heterochromatin
      case _ =>
        YRegionAnnotator.grch38Heterochromatin // Default to GRCh38
    }
    YRegionAnnotator.fromRegions(heterochromatin = heterochromatin)
  }

  /**
   * Look up biosample UUID by sample accession.
   * Searches the workspace samples for matching accession and extracts UUID from atUri.
   */
  def getBiosampleIdByAccession(sampleAccession: String): Option[UUID] =
    samples.find(_.sampleAccession == sampleAccession).flatMap { sample =>
      // Try to extract UUID from atUri (format: at://did:plc:xxx/collection/uuid)
      sample.atUri.flatMap { uri =>
        Try(UUID.fromString(uri.split("/").last)).toOption
      }
    }

  // --- Vendor VCF Import ---

  /**
   * Import a vendor-provided VCF (and optional target regions BED) for a sequence run.
   * This is for vendor deliverables like FTDNA Big Y that don't include BAM files.
   *
   * The VCF will be stored in the cache at the SequenceRun level and automatically
   * used for haplogroup analysis when available.
   *
   * @param sampleAccession  Sample accession
   * @param sequenceRunIndex Index of the sequence run within the subject
   * @param vcfPath          Path to the VCF file
   * @param bedPath          Optional path to the target regions BED file
   * @param vendor           The vendor that provided the files (e.g., FTDNA_BIGY)
   * @param referenceBuild   Reference genome build (e.g., "GRCh38")
   * @param notes            Optional notes about this import
   * @return Either error message or success message
   */
  def importVendorVcf(
                       sampleAccession: String,
                       sequenceRunIndex: Int,
                       vcfPath: java.nio.file.Path,
                       bedPath: Option[java.nio.file.Path],
                       vendor: VcfVendor,
                       referenceBuild: String,
                       notes: Option[String] = None
                     ): Either[String, String] = {
    import com.decodingus.analysis.{SubjectArtifactCache, VcfCache}

    findSubject(sampleAccession) match {
      case None =>
        Left(s"Subject not found: $sampleAccession")

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.size) {
          return Left(s"Invalid sequence run index: $sequenceRunIndex")
        }

        val seqRun = sequenceRuns(sequenceRunIndex)
        val runId = SubjectArtifactCache.extractIdFromUri(seqRun.atUri.getOrElse("unknown"))

        // Validate VCF file exists
        if (!java.nio.file.Files.exists(vcfPath)) {
          return Left(s"VCF file not found: $vcfPath")
        }

        // Validate BED file exists if provided
        bedPath.foreach { path =>
          if (!java.nio.file.Files.exists(path)) {
            return Left(s"BED file not found: $path")
          }
        }

        // Import the vendor VCF
        VcfCache.importRunVendorVcf(
          sampleAccession = sampleAccession,
          runId = runId,
          vcfSourcePath = vcfPath,
          bedSourcePath = bedPath,
          vendor = vendor,
          referenceBuild = referenceBuild,
          notes = notes
        ) match {
          case Right(info) =>
            log.info(s"Imported ${vendor.displayName} VCF for $sampleAccession run $runId: ${info.variantCount} variants")
            Right(s"Successfully imported ${vendor.displayName} VCF with ${info.variantCount} variants")

          case Left(error) =>
            Left(s"Failed to import vendor VCF: $error")
        }
    }
  }

  /**
   * List vendor VCFs imported for a sequence run.
   */
  def listVendorVcfsForRun(
                            sampleAccession: String,
                            sequenceRunIndex: Int
                          ): List[VendorVcfInfo] = {
    import com.decodingus.analysis.{SubjectArtifactCache, VcfCache}

    findSubject(sampleAccession) match {
      case None => List.empty

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.size) {
          return List.empty
        }

        val seqRun = sequenceRuns(sequenceRunIndex)
        val runId = SubjectArtifactCache.extractIdFromUri(seqRun.atUri.getOrElse("unknown"))

        VcfCache.listRunVendorVcfs(sampleAccession, runId)
    }
  }

  /**
   * Delete a vendor VCF from a sequence run.
   */
  def deleteVendorVcf(
                       sampleAccession: String,
                       sequenceRunIndex: Int,
                       vendor: VcfVendor
                     ): Boolean = {
    import com.decodingus.analysis.{SubjectArtifactCache, VcfCache}

    findSubject(sampleAccession) match {
      case None => false

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.size) {
          return false
        }

        val seqRun = sequenceRuns(sequenceRunIndex)
        val runId = SubjectArtifactCache.extractIdFromUri(seqRun.atUri.getOrElse("unknown"))

        VcfCache.deleteRunVendorVcf(sampleAccession, runId, vendor)
    }
  }

  // --- Vendor FASTA Management (mtDNA) ---

  /**
   * Import a vendor-provided mtDNA FASTA file for a sequence run.
   *
   * @param sampleAccession  The sample accession identifier
   * @param sequenceRunIndex Index of the sequence run
   * @param fastaPath        Path to the FASTA file
   * @param vendor           Vendor type (e.g., FTDNA_MTFULL, YSEQ)
   * @param notes            Optional notes about the import
   * @return Either error message or success message
   */
  def importVendorFasta(
                         sampleAccession: String,
                         sequenceRunIndex: Int,
                         fastaPath: java.nio.file.Path,
                         vendor: VcfVendor,
                         notes: Option[String] = None
                       ): Either[String, String] = {
    import com.decodingus.analysis.{MtDnaFastaProcessor, SubjectArtifactCache, VcfCache}

    findSubject(sampleAccession) match {
      case None =>
        Left(s"Subject not found: $sampleAccession")

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.size) {
          return Left(s"Invalid sequence run index: $sequenceRunIndex")
        }

        val seqRun = sequenceRuns(sequenceRunIndex)
        val runId = SubjectArtifactCache.extractIdFromUri(seqRun.atUri.getOrElse("unknown"))

        // Validate the FASTA file
        MtDnaFastaProcessor.readFasta(fastaPath) match {
          case Left(error) =>
            Left(s"Invalid FASTA file: $error")

          case Right(sequence) =>
            if (sequence.length < 16000 || sequence.length > 17000) {
              Left(s"Unexpected mtDNA sequence length: ${sequence.length} bp (expected ~16569)")
            } else {
              // Import the FASTA
              VcfCache.importRunFasta(sampleAccession, runId, fastaPath, vendor, notes) match {
                case Right(info) =>
                  Right(s"Imported ${vendor.displayName} mtDNA FASTA\nSequence length: ${info.sequenceLength} bp")
                case Left(error) =>
                  Left(error)
              }
            }
        }
    }
  }

  /**
   * List all imported vendor FASTAs for a sequence run.
   */
  def listVendorFastasForRun(
                              sampleAccession: String,
                              sequenceRunIndex: Int
                            ): List[VendorFastaInfo] = {
    import com.decodingus.analysis.{SubjectArtifactCache, VcfCache}

    findSubject(sampleAccession) match {
      case None => List.empty

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.size) {
          return List.empty
        }

        val seqRun = sequenceRuns(sequenceRunIndex)
        val runId = SubjectArtifactCache.extractIdFromUri(seqRun.atUri.getOrElse("unknown"))

        VcfCache.listRunFastas(sampleAccession, runId)
    }
  }

  /**
   * Delete a vendor FASTA from a sequence run.
   */
  def deleteVendorFasta(
                         sampleAccession: String,
                         sequenceRunIndex: Int,
                         vendor: VcfVendor
                       ): Boolean = {
    import com.decodingus.analysis.{SubjectArtifactCache, VcfCache}

    findSubject(sampleAccession) match {
      case None => false

      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (sequenceRunIndex < 0 || sequenceRunIndex >= sequenceRuns.size) {
          return false
        }

        val seqRun = sequenceRuns(sequenceRunIndex)
        val runId = SubjectArtifactCache.extractIdFromUri(seqRun.atUri.getOrElse("unknown"))

        VcfCache.deleteRunFasta(sampleAccession, runId, vendor)
    }
  }
}
