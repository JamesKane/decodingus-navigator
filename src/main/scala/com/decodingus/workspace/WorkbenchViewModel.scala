package com.decodingus.workspace

import com.decodingus.analysis.{HaplogroupProcessor, LibraryStatsProcessor, WgsMetricsProcessor}
import com.decodingus.haplogroup.tree.{TreeType, TreeProviderType}
import com.decodingus.auth.User
import com.decodingus.config.FeatureToggles
import com.decodingus.model.{LibraryStats, WgsMetrics}
import com.decodingus.pds.PdsClient
import com.decodingus.config.ReferenceConfigService
import com.decodingus.refgenome.{ReferenceGateway, ReferenceResolveResult}
import com.decodingus.workspace.model.{Workspace, Project, Biosample, WorkspaceContent, SyncStatus, SequenceData, AlignmentData, AlignmentMetrics, FileInfo, HaplogroupAssignments, HaplogroupResult => WorkspaceHaplogroupResult}
import com.decodingus.haplogroup.model.{HaplogroupResult => AnalysisHaplogroupResult}
import htsjdk.samtools.SamReaderFactory
import scalafx.beans.property.{ObjectProperty, ReadOnlyObjectProperty, StringProperty, BooleanProperty, DoubleProperty}
import scalafx.collections.ObservableBuffer
import scalafx.application.Platform

import java.io.File
import java.time.LocalDateTime
import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.{Success, Failure, Try}

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
    // Also refresh filtered lists
    applyFilters()
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

  // --- Subject CRUD Operations ---

  /** Creates a new subject and adds it to the workspace */
  def addSubject(newBiosample: Biosample): Unit = {
    val updatedSamples = _workspace.value.main.samples :+ newBiosample
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))

    // Select the newly added subject
    selectedSubject.value = Some(newBiosample)

    saveWorkspace()
  }

  /** Updates an existing subject identified by sampleAccession */
  def updateSubject(updatedBiosample: Biosample): Unit = {
    val updatedSamples = _workspace.value.main.samples.map { sample =>
      if (sample.sampleAccession == updatedBiosample.sampleAccession) updatedBiosample
      else sample
    }
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))

    // Update selection to reflect changes
    selectedSubject.value = Some(updatedBiosample)

    saveWorkspace()
  }

  /** Deletes a subject identified by sampleAccession */
  def deleteSubject(sampleAccession: String): Unit = {
    val updatedSamples = _workspace.value.main.samples.filterNot(_.sampleAccession == sampleAccession)
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))

    // Clear selection if the deleted subject was selected
    selectedSubject.value match {
      case Some(selected) if selected.sampleAccession == sampleAccession =>
        selectedSubject.value = None
      case _ => // Keep current selection
    }

    saveWorkspace()
  }

  /** Finds a subject by sampleAccession */
  def findSubject(sampleAccession: String): Option[Biosample] = {
    _workspace.value.main.samples.find(_.sampleAccession == sampleAccession)
  }

  // --- Project CRUD Operations ---

  def addProject(newProject: Project): Unit = {
    val updatedProjects = _workspace.value.main.projects :+ newProject
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))

    selectedProject.value = Some(newProject)

    saveWorkspace()
  }

  /** Updates an existing project identified by projectName */
  def updateProject(updatedProject: Project): Unit = {
    val updatedProjects = _workspace.value.main.projects.map { project =>
      if (project.projectName == updatedProject.projectName) updatedProject
      else project
    }
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))

    selectedProject.value = Some(updatedProject)

    saveWorkspace()
  }

  /** Deletes a project by projectName */
  def deleteProject(projectName: String): Unit = {
    val updatedProjects = _workspace.value.main.projects.filterNot(_.projectName == projectName)
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))

    selectedProject.value match {
      case Some(selected) if selected.projectName == projectName =>
        selectedProject.value = None
      case _ =>
    }

    saveWorkspace()
  }

  /** Finds a project by projectName */
  def findProject(projectName: String): Option[Project] = {
    _workspace.value.main.projects.find(_.projectName == projectName)
  }

  /** Adds a subject (by accession) to a project's members list */
  def addSubjectToProject(projectName: String, sampleAccession: String): Boolean = {
    findProject(projectName) match {
      case Some(project) =>
        if (project.members.contains(sampleAccession)) {
          println(s"[ViewModel] Subject $sampleAccession already in project $projectName")
          false
        } else {
          val updatedProject = project.copy(members = project.members :+ sampleAccession)
          updateProject(updatedProject)
          true
        }
      case None =>
        println(s"[ViewModel] Project $projectName not found")
        false
    }
  }

  /** Removes a subject (by accession) from a project's members list */
  def removeSubjectFromProject(projectName: String, sampleAccession: String): Boolean = {
    findProject(projectName) match {
      case Some(project) =>
        if (!project.members.contains(sampleAccession)) {
          println(s"[ViewModel] Subject $sampleAccession not in project $projectName")
          false
        } else {
          val updatedProject = project.copy(members = project.members.filterNot(_ == sampleAccession))
          updateProject(updatedProject)
          true
        }
      case None =>
        println(s"[ViewModel] Project $projectName not found")
        false
    }
  }

  /** Gets subjects that are members of a project */
  def getProjectMembers(projectName: String): List[Biosample] = {
    findProject(projectName) match {
      case Some(project) =>
        project.members.flatMap(accession => findSubject(accession))
      case None =>
        List.empty
    }
  }

  /** Gets subjects that are NOT members of a project (for adding) */
  def getNonProjectMembers(projectName: String): List[Biosample] = {
    findProject(projectName) match {
      case Some(project) =>
        _workspace.value.main.samples.filterNot(s => project.members.contains(s.sampleAccession))
      case None =>
        _workspace.value.main.samples
    }
  }

  // --- SequenceData CRUD Operations (nested within a Biosample) ---

  /**
   * Creates a new SequenceData entry from just a FileInfo.
   * All metadata (platform, reads, etc.) will be populated during analysis.
   * Returns the index of the new entry, or -1 if duplicate/error.
   */
  def addSequenceDataFromFile(sampleAccession: String, fileInfo: FileInfo): Int = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Check for duplicate by checksum
        val existingChecksums = subject.sequenceData.flatMap(_.files.flatMap(_.checksum)).toSet
        if (fileInfo.checksum.exists(existingChecksums.contains)) {
          println(s"[ViewModel] Duplicate file detected: ${fileInfo.fileName}")
          return -1
        }

        val newSequenceData = SequenceData(
          platformName = "Unknown", // Will be inferred during analysis
          instrumentModel = None,
          testType = "Unknown",
          libraryLayout = None,
          totalReads = None,
          readLength = None,
          meanInsertSize = None,
          files = List(fileInfo),
          alignments = List.empty
        )

        val newIndex = subject.sequenceData.size
        val updatedSubject = subject.copy(
          sequenceData = subject.sequenceData :+ newSequenceData
        )
        updateSubject(updatedSubject)
        newIndex

      case None =>
        println(s"[ViewModel] Cannot add sequence data: subject $sampleAccession not found")
        -1
    }
  }

  /**
   * Adds a file and immediately runs library stats analysis.
   * This is the primary flow for adding new sequencing data.
   *
   * Pipeline:
   * 1. Add file to subject (creates SequenceData entry)
   * 2. Detect reference build from BAM/CRAM header
   * 3. Resolve reference genome
   * 4. Run library stats analysis
   * 5. Update SequenceData with inferred metadata
   *
   * @param sampleAccession The subject's accession ID
   * @param fileInfo The file information
   * @param onProgress Progress callback (message, percent)
   * @param onComplete Completion callback with result
   */
  def addFileAndAnalyze(
    sampleAccession: String,
    fileInfo: FileInfo,
    onProgress: (String, Double) => Unit,
    onComplete: Either[String, (Int, LibraryStats)] => Unit
  ): Unit = {
    // Step 1: Add the file
    onProgress("Adding file...", 0.05)
    val index = addSequenceDataFromFile(sampleAccession, fileInfo)

    if (index < 0) {
      onComplete(Left("Duplicate file - this file has already been added"))
      return
    }

    val bamPath = fileInfo.location
    println(s"[ViewModel] Starting add+analyze pipeline for ${fileInfo.fileName}")

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

        // Step 5: Update the SequenceData with results
        onProgress("Saving results...", 0.95)
        updateProgress("Saving results...", 0.95)

        Platform.runLater {
          // Get the current sequence data (it may have been updated)
          findSubject(sampleAccession).flatMap(_.sequenceData.lift(index)) match {
            case Some(seqData) =>
              val updatedSeqData = seqData.copy(
                platformName = libraryStats.inferredPlatform,
                instrumentModel = Some(libraryStats.mostFrequentInstrument),
                testType = inferTestType(libraryStats),
                libraryLayout = Some(if (libraryStats.pairedReads > libraryStats.readCount / 2) "Paired-End" else "Single-End"),
                totalReads = Some(libraryStats.readCount.toLong),
                readLength = libraryStats.lengthDistribution.keys.maxOption,
                meanInsertSize = calculateMeanInsertSize(libraryStats.insertSizeDistribution),
                alignments = List(AlignmentData(
                  referenceBuild = libraryStats.referenceBuild,
                  aligner = libraryStats.aligner,
                  files = seqData.files,
                  metrics = None
                ))
              )
              updateSequenceData(sampleAccession, index, updatedSeqData)

              lastLibraryStats.value = Some(libraryStats)
              analysisInProgress.value = false
              analysisProgress.value = "Analysis complete"
              analysisProgressPercent.value = 1.0
              onProgress("Complete", 1.0)
              onComplete(Right((index, libraryStats)))

            case None =>
              analysisInProgress.value = false
              onComplete(Left("Sequence data entry was removed during analysis"))
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

    if (stats.inferredPlatform == "PacBio" && avgReadLength > 10000) "HiFi"
    else if (stats.inferredPlatform == "PacBio") "CLR"
    else if (stats.inferredPlatform == "Oxford Nanopore") "Nanopore"
    else "WGS" // Default assumption for Illumina/MGI
  }

  /** Calculate mean insert size from distribution */
  private def calculateMeanInsertSize(distribution: Map[Long, Int]): Option[Double] = {
    if (distribution.isEmpty) None
    else {
      val total = distribution.map { case (size, count) => size * count }.sum
      val count = distribution.values.sum
      if (count > 0) Some(total.toDouble / count) else None
    }
  }

  /** Gets all checksums for a subject's sequence data files */
  def getExistingChecksums(sampleAccession: String): Set[String] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        subject.sequenceData.flatMap(_.files.flatMap(_.checksum)).toSet
      case None =>
        Set.empty
    }
  }

  /** Adds sequencing data to the currently selected subject */
  def addSequenceData(sequenceData: SequenceData): Unit = {
    selectedSubject.value match {
      case Some(subject) =>
        val updatedSubject = subject.copy(
          sequenceData = subject.sequenceData :+ sequenceData
        )
        updateSubject(updatedSubject)
      case None =>
        println("[ViewModel] Cannot add sequence data: no subject selected")
    }
  }

  /** Adds sequencing data to a specific subject by accession */
  def addSequenceDataToSubject(sampleAccession: String, sequenceData: SequenceData): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val updatedSubject = subject.copy(
          sequenceData = subject.sequenceData :+ sequenceData
        )
        updateSubject(updatedSubject)
      case None =>
        println(s"[ViewModel] Cannot add sequence data: subject $sampleAccession not found")
    }
  }

  /** Removes sequencing data from a subject by index */
  def removeSequenceData(sampleAccession: String, index: Int): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) if index >= 0 && index < subject.sequenceData.size =>
        val updatedSeqData = subject.sequenceData.zipWithIndex.filterNot(_._2 == index).map(_._1)
        val updatedSubject = subject.copy(sequenceData = updatedSeqData)
        updateSubject(updatedSubject)
      case Some(_) =>
        println(s"[ViewModel] Cannot remove sequence data: index $index out of bounds")
      case None =>
        println(s"[ViewModel] Cannot remove sequence data: subject $sampleAccession not found")
    }
  }

  /** Updates sequencing data at a specific index for a subject */
  def updateSequenceData(sampleAccession: String, index: Int, updatedData: SequenceData): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) if index >= 0 && index < subject.sequenceData.size =>
        val updatedSeqData = subject.sequenceData.updated(index, updatedData)
        val updatedSubject = subject.copy(sequenceData = updatedSeqData)
        updateSubject(updatedSubject)
      case Some(_) =>
        println(s"[ViewModel] Cannot update sequence data: index $index out of bounds")
      case None =>
        println(s"[ViewModel] Cannot update sequence data: subject $sampleAccession not found")
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
   * @param onProgress Progress callback for download
   * @param onResolved Called with the resolved path when available
   * @param onError Called if resolution fails or is cancelled
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
                      Platform.runLater { onResolved(path.toString) }
                    case Left(error) =>
                      Platform.runLater { onError(error) }
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
                Platform.runLater { onResolved(path.toString) }
              case Left(error) =>
                Platform.runLater { onError(error) }
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
   * 4. Updates the SequenceData with results
   *
   * @param sampleAccession The subject's accession ID
   * @param sequenceDataIndex The index of the sequence data entry to analyze
   * @param onComplete Callback when analysis completes (success or failure)
   */
  def runInitialAnalysis(
    sampleAccession: String,
    sequenceDataIndex: Int,
    onComplete: Either[String, LibraryStats] => Unit
  ): Unit = {
    findSubject(sampleAccession).flatMap { subject =>
      subject.sequenceData.lift(sequenceDataIndex).map((subject, _))
    } match {
      case Some((subject, seqData)) =>
        seqData.files.headOption match {
          case Some(fileInfo) =>
            val bamPath = fileInfo.location
            println(s"[ViewModel] Starting initial analysis for ${fileInfo.fileName}")

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

                // Step 4: Update the SequenceData with results
                updateProgress("Saving results...", 0.95)
                Platform.runLater {
                  val updatedSeqData = seqData.copy(
                    platformName = if (seqData.platformName == "Other") libraryStats.inferredPlatform else seqData.platformName,
                    instrumentModel = seqData.instrumentModel.orElse(Some(libraryStats.mostFrequentInstrument)),
                    totalReads = Some(libraryStats.readCount.toLong),
                    alignments = List(AlignmentData(
                      referenceBuild = libraryStats.referenceBuild,
                      aligner = libraryStats.aligner,
                      files = seqData.files,
                      metrics = None // Will be filled by WGS metrics analysis
                    ))
                  )
                  updateSequenceData(sampleAccession, sequenceDataIndex, updatedSeqData)
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
            onComplete(Left("No alignment file associated with this sequence data"))
        }

      case None =>
        onComplete(Left(s"Sequence data not found at index $sequenceDataIndex for subject $sampleAccession"))
    }
  }

  /**
   * Runs WGS metrics analysis (deep coverage) for a sequencing run.
   * Requires that initial analysis has already been run.
   *
   * @param sampleAccession The subject's accession ID
   * @param sequenceDataIndex The index of the sequence data entry to analyze
   * @param onComplete Callback when analysis completes
   */
  def runWgsMetricsAnalysis(
    sampleAccession: String,
    sequenceDataIndex: Int,
    onComplete: Either[String, WgsMetrics] => Unit
  ): Unit = {
    findSubject(sampleAccession).flatMap { subject =>
      subject.sequenceData.lift(sequenceDataIndex).map((subject, _))
    } match {
      case Some((subject, seqData)) =>
        // Need alignment data with reference build
        seqData.alignments.headOption match {
          case Some(alignmentData) =>
            seqData.files.headOption match {
              case Some(fileInfo) =>
                val bamPath = fileInfo.location
                println(s"[ViewModel] Starting WGS metrics analysis for ${fileInfo.fileName}")

                analysisInProgress.value = true
                analysisError.value = ""
                analysisProgress.value = "Starting WGS metrics analysis..."
                analysisProgressPercent.value = 0.0

                Future {
                  try {
                    // Resolve reference
                    updateProgress("Resolving reference genome...", 0.1)
                    val referenceGateway = new ReferenceGateway((_, _) => {})
                    val referencePath = referenceGateway.resolve(alignmentData.referenceBuild) match {
                      case Right(path) => path.toString
                      case Left(error) => throw new Exception(s"Failed to resolve reference: $error")
                    }

                    // Run GATK CollectWgsMetrics
                    updateProgress("Running GATK CollectWgsMetrics (this may take a while)...", 0.2)
                    val wgsProcessor = new WgsMetricsProcessor()
                    val wgsMetrics = wgsProcessor.process(bamPath, referencePath, (message, current, total) => {
                      val pct = 0.2 + (current / total) * 0.7
                      updateProgress(message, pct)
                    }) match {
                      case Right(metrics) => metrics
                      case Left(error) => throw error
                    }

                    // Update alignment metrics
                    updateProgress("Saving results...", 0.95)
                    Platform.runLater {
                      val updatedMetrics = AlignmentMetrics(
                        genomeTerritory = Some(wgsMetrics.genomeTerritory),
                        meanCoverage = Some(wgsMetrics.meanCoverage),
                        medianCoverage = Some(wgsMetrics.medianCoverage),
                        sdCoverage = Some(wgsMetrics.sdCoverage),
                        pctExcDupe = Some(wgsMetrics.pctExcDupe),
                        pctExcMapq = Some(wgsMetrics.pctExcMapq),
                        pct10x = Some(wgsMetrics.pct10x),
                        pct20x = Some(wgsMetrics.pct20x),
                        pct30x = Some(wgsMetrics.pct30x),
                        hetSnpSensitivity = Some(wgsMetrics.hetSnpSensitivity),
                        contigs = List.empty // TODO: Add per-contig metrics if needed
                      )
                      val updatedAlignment = alignmentData.copy(metrics = Some(updatedMetrics))
                      val updatedSeqData = seqData.copy(alignments = List(updatedAlignment))
                      updateSequenceData(sampleAccession, sequenceDataIndex, updatedSeqData)

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
                onComplete(Left("No alignment file associated with this sequence data"))
            }

          case None =>
            onComplete(Left("Please run initial analysis first to detect reference build"))
        }

      case None =>
        onComplete(Left(s"Sequence data not found at index $sequenceDataIndex for subject $sampleAccession"))
    }
  }

  /** Helper to update progress on the JavaFX thread */
  private def updateProgress(message: String, percent: Double): Unit = {
    Platform.runLater {
      analysisProgress.value = message
      analysisProgressPercent.value = percent
    }
  }

  // --- Haplogroup Analysis ---

  // Store last haplogroup analysis result
  val lastHaplogroupResult: ObjectProperty[Option[AnalysisHaplogroupResult]] = ObjectProperty(None)

  /**
   * Runs haplogroup analysis for a subject using the specified tree type.
   * Uses FTDNA tree provider for both Y-DNA and MT-DNA.
   *
   * @param sampleAccession The subject's accession ID
   * @param sequenceDataIndex The index of the sequence data entry to analyze
   * @param treeType The type of haplogroup tree (YDNA or MTDNA)
   * @param onComplete Callback when analysis completes
   */
  def runHaplogroupAnalysis(
    sampleAccession: String,
    sequenceDataIndex: Int,
    treeType: TreeType,
    onComplete: Either[String, AnalysisHaplogroupResult] => Unit
  ): Unit = {
    findSubject(sampleAccession).flatMap { subject =>
      subject.sequenceData.lift(sequenceDataIndex).map((subject, _))
    } match {
      case Some((subject, seqData)) =>
        // Need alignment data with reference build for haplogroup analysis
        seqData.alignments.headOption match {
          case Some(alignmentData) =>
            seqData.files.headOption match {
              case Some(fileInfo) =>
                val bamPath = fileInfo.location
                println(s"[ViewModel] Starting ${treeType} haplogroup analysis for ${fileInfo.fileName}")

                analysisInProgress.value = true
                analysisError.value = ""
                analysisProgress.value = "Starting haplogroup analysis..."
                analysisProgressPercent.value = 0.0

                Future {
                  try {
                    // Build LibraryStats from existing data for the processor
                    val libraryStats = LibraryStats(
                      readCount = seqData.totalReads.map(_.toInt).getOrElse(0),
                      pairedReads = 0,
                      lengthDistribution = Map.empty,
                      insertSizeDistribution = Map.empty,
                      aligner = alignmentData.aligner,
                      referenceBuild = alignmentData.referenceBuild,
                      sampleName = subject.donorIdentifier,
                      flowCells = Map.empty,
                      instruments = Map.empty,
                      mostFrequentInstrument = seqData.instrumentModel.getOrElse("Unknown"),
                      inferredPlatform = seqData.platformName,
                      platformCounts = Map.empty
                    )

                    updateProgress("Loading haplogroup tree...", 0.1)

                    val processor = new HaplogroupProcessor()
                    val result = processor.analyze(
                      bamPath,
                      libraryStats,
                      treeType,
                      TreeProviderType.FTDNA, // Using FTDNA for POC
                      (message, current, total) => {
                        val pct = if (total > 0) current / total else 0.0
                        updateProgress(message, pct)
                      }
                    )

                    result match {
                      case Right(results) if results.nonEmpty =>
                        val topResult = results.head
                        Platform.runLater {
                          // Convert analysis result to workspace model and save
                          val workspaceResult = WorkspaceHaplogroupResult(
                            haplogroupName = topResult.name,
                            score = topResult.score,
                            matchingSnps = Some(topResult.matchingSnps),
                            mismatchingSnps = Some(topResult.mismatchingSnps),
                            ancestralMatches = Some(topResult.ancestralMatches),
                            treeDepth = Some(topResult.depth),
                            lineagePath = None // Could be populated if needed
                          )

                          // Update subject's haplogroup assignments
                          val currentAssignments = subject.haplogroups.getOrElse(HaplogroupAssignments(None, None))
                          val updatedAssignments = treeType match {
                            case TreeType.YDNA => currentAssignments.copy(yDna = Some(workspaceResult))
                            case TreeType.MTDNA => currentAssignments.copy(mtDna = Some(workspaceResult))
                          }
                          val updatedSubject = subject.copy(haplogroups = Some(updatedAssignments))
                          updateSubject(updatedSubject)

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
                onComplete(Left("No alignment file associated with this sequence data"))
            }

          case None =>
            onComplete(Left("Please run initial analysis first to detect reference build"))
        }

      case None =>
        onComplete(Left(s"Sequence data not found at index $sequenceDataIndex for subject $sampleAccession"))
    }
  }

  /**
   * Gets the current haplogroup assignments for a subject.
   */
  def getHaplogroupAssignments(sampleAccession: String): Option[HaplogroupAssignments] = {
    findSubject(sampleAccession).flatMap(_.haplogroups)
  }
}
