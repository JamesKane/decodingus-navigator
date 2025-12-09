package com.decodingus.workspace

import com.decodingus.analysis.{ArtifactContext, HaplogroupProcessor, LibraryStatsProcessor, WgsMetricsProcessor}
import com.decodingus.haplogroup.tree.{TreeType, TreeProviderType}
import com.decodingus.auth.User
import com.decodingus.config.FeatureToggles
import com.decodingus.model.{LibraryStats, WgsMetrics}
import com.decodingus.pds.PdsClient
import com.decodingus.config.ReferenceConfigService
import com.decodingus.refgenome.{ReferenceGateway, ReferenceResolveResult}
import com.decodingus.workspace.model.{Workspace, Project, Biosample, WorkspaceContent, SyncStatus, SequenceRun, Alignment, AlignmentMetrics, FileInfo, HaplogroupAssignments, HaplogroupResult => WorkspaceHaplogroupResult, RecordMeta, StrProfile, ChipProfile}
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
    println(s"[DEBUG] WorkbenchViewModel: _workspace changed from (samples: ${oldWorkspace.main.samples.size}, projects: ${oldWorkspace.main.projects.size}) to (samples: ${newWorkspace.main.samples.size}, projects: ${newWorkspace.main.projects.size})")
    syncBuffers(newWorkspace)
  }

  // --- Initialization ---
  // Load workspace on creation of ViewModel (AFTER onChange listener is registered)
  loadWorkspace()

  /** Synchronizes the observable buffers with the workspace state */
  private def syncBuffers(workspace: Workspace): Unit = {
    println(s"[DEBUG] WorkbenchViewModel: Syncing buffers with workspace...")

    // Preserve current selection identifiers before clearing
    val selectedProjectName = selectedProject.value.map(_.projectName)
    val selectedSampleAccession = selectedSubject.value.map(_.sampleAccession)

    projects.clear()
    projects ++= workspace.main.projects
    samples.clear()
    samples ++= workspace.main.samples
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

  private def emptyWorkspace: Workspace = Workspace.empty

  // --- Subject CRUD Operations ---

  /** Creates a new subject and adds it to the workspace */
  def addSubject(newBiosample: Biosample): Unit = {
    // Generate atUri from current user's DID
    val enrichedBiosample = currentUser.value match {
      case Some(user) =>
        val atUri = s"at://${user.did}/com.decodingus.atmosphere.biosample/${newBiosample.sampleAccession}"
        newBiosample.copy(atUri = Some(atUri))
      case None =>
        newBiosample // Keep as-is if no user logged in
    }

    val updatedSamples = _workspace.value.main.samples :+ enrichedBiosample
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))

    // Select the newly added subject
    selectedSubject.value = Some(enrichedBiosample)

    saveWorkspace()
  }

  /** Updates an existing subject identified by sampleAccession */
  def updateSubject(updatedBiosample: Biosample): Unit = {
    val updatedSamples = _workspace.value.main.samples.map { sample =>
      if (sample.sampleAccession == updatedBiosample.sampleAccession) {
        // Update meta to track the edit
        updatedBiosample.copy(meta = sample.meta.updated("edit"))
      } else sample
    }
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))

    // Update selection to reflect changes (with updated meta)
    val finalUpdated = updatedSamples.find(_.sampleAccession == updatedBiosample.sampleAccession)
    selectedSubject.value = finalUpdated

    saveWorkspace()
  }

  /** Internal: Updates a subject without modifying meta (used when meta is already updated) */
  private def updateSubjectDirect(updatedBiosample: Biosample): Unit = {
    val updatedSamples = _workspace.value.main.samples.map { sample =>
      if (sample.sampleAccession == updatedBiosample.sampleAccession) updatedBiosample
      else sample
    }
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(samples = updatedSamples))
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
    // Generate atUri for the project using current user's DID
    val enrichedProject = currentUser.value match {
      case Some(user) =>
        val rkey = java.util.UUID.randomUUID().toString
        val atUri = s"at://${user.did}/com.decodingus.atmosphere.project/$rkey"
        newProject.copy(atUri = Some(atUri))
      case None =>
        newProject // Keep as-is if no user logged in
    }

    val updatedProjects = _workspace.value.main.projects :+ enrichedProject
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))

    selectedProject.value = Some(enrichedProject)

    saveWorkspace()
  }

  /** Updates an existing project identified by projectName */
  def updateProject(updatedProject: Project): Unit = {
    val updatedProjects = _workspace.value.main.projects.map { project =>
      if (project.projectName == updatedProject.projectName) {
        // Update meta to track the edit
        updatedProject.copy(meta = project.meta.updated("edit"))
      } else project
    }
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))

    // Update selection to reflect changes (with updated meta)
    val finalUpdated = updatedProjects.find(_.projectName == updatedProject.projectName)
    selectedProject.value = finalUpdated

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

  /** Backfills atUri for any samples/projects that were created while logged out */
  private def backfillAtUris(did: String): Unit = {
    val currentWorkspace = _workspace.value
    var updated = false

    // Backfill samples missing atUri
    val updatedSamples = currentWorkspace.main.samples.map { sample =>
      if (sample.atUri.isEmpty) {
        updated = true
        val atUri = s"at://$did/com.decodingus.atmosphere.biosample/${sample.sampleAccession}"
        sample.copy(atUri = Some(atUri))
      } else sample
    }

    // Backfill projects missing atUri
    val updatedProjects = currentWorkspace.main.projects.map { project =>
      if (project.atUri.isEmpty) {
        updated = true
        val rkey = java.util.UUID.randomUUID().toString
        val atUri = s"at://$did/com.decodingus.atmosphere.project/$rkey"
        project.copy(atUri = Some(atUri))
      } else project
    }

    if (updated) {
      _workspace.value = currentWorkspace.copy(
        main = currentWorkspace.main.copy(samples = updatedSamples, projects = updatedProjects)
      )
      saveWorkspace()
      println(s"[ViewModel] Backfilled atUri for samples/projects after login")
    }
  }

  /** Adds a subject (by accession) to a project's members list */
  def addSubjectToProject(projectName: String, sampleAccession: String): Boolean = {
    findProject(projectName) match {
      case Some(project) =>
        if (project.memberRefs.contains(sampleAccession)) {
          println(s"[ViewModel] Subject $sampleAccession already in project $projectName")
          false
        } else {
          val updatedProject = project.copy(
            memberRefs = project.memberRefs :+ sampleAccession,
            meta = project.meta.updated("memberRefs")
          )
          updateProjectDirect(updatedProject)
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
        if (!project.memberRefs.contains(sampleAccession)) {
          println(s"[ViewModel] Subject $sampleAccession not in project $projectName")
          false
        } else {
          val updatedProject = project.copy(
            memberRefs = project.memberRefs.filterNot(_ == sampleAccession),
            meta = project.meta.updated("memberRefs")
          )
          updateProjectDirect(updatedProject)
          true
        }
      case None =>
        println(s"[ViewModel] Project $projectName not found")
        false
    }
  }

  /** Internal: Updates a project without modifying meta (used when meta is already updated) */
  private def updateProjectDirect(updatedProject: Project): Unit = {
    val updatedProjects = _workspace.value.main.projects.map { project =>
      if (project.projectName == updatedProject.projectName) updatedProject
      else project
    }
    _workspace.value = _workspace.value.copy(main = _workspace.value.main.copy(projects = updatedProjects))
    selectedProject.value = Some(updatedProject)
    saveWorkspace()
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
   */
  def addSequenceRunFromFile(sampleAccession: String, fileInfo: FileInfo): Int = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Check for duplicate by checksum across all sequence runs for this biosample
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        val existingChecksums = sequenceRuns.flatMap(_.files.flatMap(_.checksum)).toSet
        if (fileInfo.checksum.exists(existingChecksums.contains)) {
          println(s"[ViewModel] Duplicate file detected: ${fileInfo.fileName}")
          return -1
        }

        // Generate a unique URI for this sequence run
        val seqRunUri = s"local:sequencerun:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"

        val newSequenceRun = SequenceRun(
          atUri = Some(seqRunUri),
          meta = RecordMeta.initial,
          biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}"),
          platformName = "Unknown", // Will be inferred during analysis
          instrumentModel = None,
          testType = "Unknown",
          libraryLayout = None,
          totalReads = None,
          readLength = None,
          meanInsertSize = None,
          files = List(fileInfo),
          alignmentRefs = List.empty
        )

        val newIndex = sequenceRuns.size

        // Update workspace with new sequence run and update biosample refs
        val updatedSequenceRuns = _workspace.value.main.sequenceRuns :+ newSequenceRun
        val updatedSubject = subject.copy(
          sequenceRunRefs = subject.sequenceRunRefs :+ seqRunUri,
          meta = subject.meta.updated("sequenceRunRefs")
        )
        val updatedSamples = _workspace.value.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }
        val updatedContent = _workspace.value.main.copy(
          samples = updatedSamples,
          sequenceRuns = updatedSequenceRuns
        )
        _workspace.value = _workspace.value.copy(main = updatedContent)
        saveWorkspace()

        newIndex

      case None =>
        println(s"[ViewModel] Cannot add sequence run: subject $sampleAccession not found")
        -1
    }
  }

  // Backward compatibility alias
  def addSequenceDataFromFile(sampleAccession: String, fileInfo: FileInfo): Int =
    addSequenceRunFromFile(sampleAccession, fileInfo)

  /**
   * Adds a file and immediately runs library stats analysis.
   * This is the primary flow for adding new sequencing data.
   *
   * Pipeline:
   * 1. Add file to subject (creates SequenceRun entry)
   * 2. Detect reference build from BAM/CRAM header
   * 3. Resolve reference genome
   * 4. Run library stats analysis
   * 5. Update SequenceRun with inferred metadata and create Alignment
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
    val index = addSequenceRunFromFile(sampleAccession, fileInfo)

    if (index < 0) {
      onComplete(Left("Duplicate file - this file has already been added"))
      return
    }

    val bamPath = fileInfo.location.getOrElse("")
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

        // Step 5: Update the SequenceRun with results and create Alignment
        onProgress("Saving results...", 0.95)
        updateProgress("Saving results...", 0.95)

        Platform.runLater {
          findSubject(sampleAccession) match {
            case Some(subject) =>
              val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
              sequenceRuns.lift(index) match {
                case Some(seqRun) =>
                  // Check for existing alignment or create new one
                  val alignUri = seqRun.alignmentRefs.headOption.getOrElse(
                    s"local:alignment:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"
                  )
                  val existingAlignment = _workspace.value.main.alignments.find(_.atUri.contains(alignUri))

                  val newAlignment = Alignment(
                    atUri = Some(alignUri),
                    meta = existingAlignment.map(_.meta.updated("analysis")).getOrElse(RecordMeta.initial),
                    sequenceRunRef = seqRun.atUri.getOrElse(""),
                    biosampleRef = Some(subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")),
                    referenceBuild = libraryStats.referenceBuild,
                    aligner = libraryStats.aligner,
                    files = seqRun.files,
                    metrics = existingAlignment.flatMap(_.metrics) // Preserve existing metrics on re-run
                  )

                  // Update the SequenceRun with inferred metadata
                  val updatedSeqRun = seqRun.copy(
                    meta = seqRun.meta.updated("analysis"),
                    platformName = libraryStats.inferredPlatform,
                    instrumentModel = Some(libraryStats.mostFrequentInstrument),
                    testType = inferTestType(libraryStats),
                    libraryLayout = Some(if (libraryStats.pairedReads > libraryStats.readCount / 2) "Paired-End" else "Single-End"),
                    totalReads = Some(libraryStats.readCount.toLong),
                    readLength = calculateMeanReadLength(libraryStats.lengthDistribution),
                    maxReadLength = libraryStats.lengthDistribution.keys.maxOption,
                    meanInsertSize = calculateMeanInsertSize(libraryStats.insertSizeDistribution),
                    alignmentRefs = if (seqRun.alignmentRefs.contains(alignUri)) seqRun.alignmentRefs else seqRun.alignmentRefs :+ alignUri
                  )

                  // Update workspace - update existing alignment or add new one
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
                  saveWorkspace()

                  lastLibraryStats.value = Some(libraryStats)
                  analysisInProgress.value = false
                  analysisProgress.value = "Analysis complete"
                  analysisProgressPercent.value = 1.0
                  onProgress("Complete", 1.0)
                  onComplete(Right((index, libraryStats)))

                case None =>
                  analysisInProgress.value = false
                  onComplete(Left("Sequence run entry was removed during analysis"))
              }
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

    if (stats.inferredPlatform == "PacBio" && avgReadLength > 10000) "HiFi"
    else if (stats.inferredPlatform == "PacBio") "CLR"
    else if (stats.inferredPlatform == "Oxford Nanopore") "Nanopore"
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

  /** Gets all checksums for a subject's sequence run files */
  def getExistingChecksums(sampleAccession: String): Set[String] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.flatMap(_.files.flatMap(_.checksum)).toSet
      case None =>
        Set.empty
    }
  }

  /** Removes a sequence run from a subject by index */
  def removeSequenceData(sampleAccession: String, index: Int): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (index >= 0 && index < sequenceRuns.size) {
          val seqRunToRemove = sequenceRuns(index)

          // Remove the sequence run from workspace
          val updatedSequenceRuns = _workspace.value.main.sequenceRuns.filterNot(_.atUri == seqRunToRemove.atUri)

          // Remove any alignments that reference this sequence run
          val updatedAlignments = _workspace.value.main.alignments.filterNot { align =>
            seqRunToRemove.atUri.exists(uri => align.sequenceRunRef == uri)
          }

          // Update subject's sequenceRunRefs
          val updatedSubject = subject.copy(
            sequenceRunRefs = subject.sequenceRunRefs.filterNot(ref => seqRunToRemove.atUri.contains(ref)),
            meta = subject.meta.updated("sequenceRunRefs")
          )
          val updatedSamples = _workspace.value.main.samples.map { s =>
            if (s.sampleAccession == sampleAccession) updatedSubject else s
          }

          val updatedContent = _workspace.value.main.copy(
            samples = updatedSamples,
            sequenceRuns = updatedSequenceRuns,
            alignments = updatedAlignments
          )
          _workspace.value = _workspace.value.copy(main = updatedContent)
          saveWorkspace()
        } else {
          println(s"[ViewModel] Cannot remove sequence run: index $index out of bounds")
        }
      case None =>
        println(s"[ViewModel] Cannot remove sequence run: subject $sampleAccession not found")
    }
  }

  /** Updates a sequence run at a specific index for a subject */
  def updateSequenceRun(sampleAccession: String, index: Int, updatedRun: SequenceRun): Unit = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        if (index >= 0 && index < sequenceRuns.size) {
          val originalRun = sequenceRuns(index)
          val runWithUpdatedMeta = updatedRun.copy(
            atUri = originalRun.atUri,
            meta = originalRun.meta.updated("edit")
          )

          val updatedSequenceRuns = _workspace.value.main.sequenceRuns.map { sr =>
            if (sr.atUri == originalRun.atUri) runWithUpdatedMeta else sr
          }

          val updatedContent = _workspace.value.main.copy(sequenceRuns = updatedSequenceRuns)
          _workspace.value = _workspace.value.copy(main = updatedContent)
          saveWorkspace()
        } else {
          println(s"[ViewModel] Cannot update sequence run: index $index out of bounds")
        }
      case None =>
        println(s"[ViewModel] Cannot update sequence run: subject $sampleAccession not found")
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
   * 4. Updates the SequenceRun with results and creates Alignment
   *
   * @param sampleAccession The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param onComplete Callback when analysis completes (success or failure)
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

                    // Step 4: Update SequenceRun and create/update Alignment
                    updateProgress("Saving results...", 0.95)
                    Platform.runLater {
                      // Create or update alignment
                      val alignUri = seqRun.alignmentRefs.headOption.getOrElse(
                        s"local:alignment:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"
                      )
                      val existingAlignment = _workspace.value.main.alignments.find(_.atUri.contains(alignUri))
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
                        totalReads = Some(libraryStats.readCount.toLong),
                        readLength = calculateMeanReadLength(libraryStats.lengthDistribution).orElse(seqRun.readLength),
                        maxReadLength = libraryStats.lengthDistribution.keys.maxOption.orElse(seqRun.maxReadLength),
                        meanInsertSize = calculateMeanInsertSize(libraryStats.insertSizeDistribution).orElse(seqRun.meanInsertSize),
                        alignmentRefs = if (seqRun.alignmentRefs.contains(alignUri)) seqRun.alignmentRefs else seqRun.alignmentRefs :+ alignUri
                      )

                      // Update workspace
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
                      saveWorkspace()

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
   * @param sampleAccession The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param onComplete Callback when analysis completes
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
                        val wgsMetrics = wgsProcessor.process(
                          bamPath,
                          referencePath,
                          (message, current, total) => {
                            val pct = 0.2 + (current / total) * 0.7
                            updateProgress(message, pct)
                          },
                          seqRun.maxReadLength, // Pass max read length to handle long reads (e.g., PacBio HiFi, NovaSeq 151bp)
                          Some(artifactCtx)
                        ) match {
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
                            contigs = List.empty
                          )
                          val updatedAlignment = alignment.copy(
                            meta = alignment.meta.updated("wgsMetrics"),
                            metrics = Some(updatedMetrics)
                          )

                          // Update workspace
                          val updatedAlignments = _workspace.value.main.alignments.map { a =>
                            if (a.atUri == alignment.atUri) updatedAlignment else a
                          }
                          val updatedContent = _workspace.value.main.copy(alignments = updatedAlignments)
                          _workspace.value = _workspace.value.copy(main = updatedContent)
                          saveWorkspace()

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

  // --- Haplogroup Analysis ---

  // Store last haplogroup analysis result
  val lastHaplogroupResult: ObjectProperty[Option[AnalysisHaplogroupResult]] = ObjectProperty(None)

  /**
   * Runs haplogroup analysis for a subject using the specified tree type.
   * Uses FTDNA tree provider for both Y-DNA and MT-DNA.
   *
   * @param sampleAccession The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param treeType The type of haplogroup tree (YDNA or MTDNA)
   * @param onComplete Callback when analysis completes
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
                    println(s"[ViewModel] Starting ${treeType} haplogroup analysis for ${fileInfo.fileName}")

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
                        // Select tree provider based on configuration
                        val treeProviderType = treeType match {
                          case TreeType.YDNA =>
                            if (FeatureToggles.treeProviders.ydna.equalsIgnoreCase("decodingus"))
                              TreeProviderType.DECODINGUS
                            else TreeProviderType.FTDNA
                          case TreeType.MTDNA =>
                            if (FeatureToggles.treeProviders.mtdna.equalsIgnoreCase("decodingus"))
                              TreeProviderType.DECODINGUS
                            else TreeProviderType.FTDNA
                        }
                        val result = processor.analyze(
                          bamPath,
                          libraryStats,
                          treeType,
                          treeProviderType,
                          (message, current, total) => {
                            val pct = if (total > 0) current / total else 0.0
                            updateProgress(message, pct)
                          },
                          Some(artifactCtx)
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
                              val updatedSubject = subject.copy(
                                haplogroups = Some(updatedAssignments),
                                meta = subject.meta.updated("haplogroups")
                              )
                              updateSubjectDirect(updatedSubject)

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

  // --- STR Profile CRUD Operations ---

  /**
   * Adds a new STR profile for a biosample.
   * Supports multiple profiles per subject (e.g., from different vendors like FTDNA and YSEQ).
   * Returns the URI of the new profile, or an error message.
   */
  def addStrProfile(sampleAccession: String, profile: StrProfile): Either[String, String] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Generate a unique URI for this STR profile
        val strProfileUri = s"local:strprofile:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"

        val enrichedProfile = profile.copy(
          atUri = Some(strProfileUri),
          biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")
        )

        // Update workspace with new STR profile
        val updatedStrProfiles = _workspace.value.main.strProfiles :+ enrichedProfile

        // Update subject's strProfileRefs list (supports multiple profiles)
        val updatedSubject = subject.copy(
          strProfileRefs = subject.strProfileRefs :+ strProfileUri,
          meta = subject.meta.updated("strProfileRefs")
        )
        val updatedSamples = _workspace.value.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }

        val updatedContent = _workspace.value.main.copy(
          samples = updatedSamples,
          strProfiles = updatedStrProfiles
        )
        _workspace.value = _workspace.value.copy(main = updatedContent)
        saveWorkspace()

        println(s"[ViewModel] Added STR profile with ${enrichedProfile.markers.size} markers for $sampleAccession (provider: ${enrichedProfile.importedFrom.getOrElse("unknown")})")
        Right(strProfileUri)

      case None =>
        Left(s"Subject not found: $sampleAccession")
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
    _workspace.value.main.strProfiles.find(_.atUri.contains(profileUri)) match {
      case Some(existing) =>
        val withUpdatedMeta = updatedProfile.copy(
          atUri = existing.atUri,
          meta = existing.meta.updated("edit")
        )

        val updatedStrProfiles = _workspace.value.main.strProfiles.map { p =>
          if (p.atUri.contains(profileUri)) withUpdatedMeta else p
        }
        val updatedContent = _workspace.value.main.copy(strProfiles = updatedStrProfiles)
        _workspace.value = _workspace.value.copy(main = updatedContent)
        saveWorkspace()
        Right(())

      case None =>
        Left(s"STR profile not found: $profileUri")
    }
  }

  /**
   * Deletes an STR profile.
   */
  def deleteStrProfile(sampleAccession: String, profileUri: String): Either[String, Unit] = {
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Remove the profile from workspace
        val updatedStrProfiles = _workspace.value.main.strProfiles.filterNot(_.atUri.contains(profileUri))

        // Remove from subject's strProfileRefs list
        val updatedSubject = subject.copy(
          strProfileRefs = subject.strProfileRefs.filterNot(_ == profileUri),
          meta = subject.meta.updated("strProfileRefs")
        )
        val updatedSamples = _workspace.value.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }

        val updatedContent = _workspace.value.main.copy(
          samples = updatedSamples,
          strProfiles = updatedStrProfiles
        )
        _workspace.value = _workspace.value.copy(main = updatedContent)
        saveWorkspace()

        println(s"[ViewModel] Deleted STR profile $profileUri for $sampleAccession")
        Right(())

      case None =>
        Left(s"Subject not found: $sampleAccession")
    }
  }

  // --- Chip Profile CRUD Operations ---

  /**
   * Imports chip data from a file for a biosample.
   * Parses the file, computes statistics, and creates a ChipProfile record.
   *
   * @param sampleAccession The subject's accession ID
   * @param file The chip data file to import
   * @param onComplete Callback when import completes
   */
  def importChipData(
    sampleAccession: String,
    file: File,
    onComplete: Either[String, ChipProfile] => Unit
  ): Unit = {
    import com.decodingus.genotype.parser.ChipDataParser
    import com.decodingus.genotype.model.GenotypingTestSummary
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
                      val chipProfileUri = s"local:chipprofile:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"

                      val fileInfo = FileInfo(
                        fileName = file.getName,
                        fileSizeBytes = Some(file.length()),
                        fileFormat = if (file.getName.endsWith(".csv")) "CSV" else "TXT",
                        checksum = Some(computeFileHash(file)),
                        checksumAlgorithm = Some("SHA-256"),
                        location = Some(file.getAbsolutePath)
                      )

                      val chipProfile = ChipProfile(
                        atUri = Some(chipProfileUri),
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

                      // Update workspace
                      val updatedChipProfiles = _workspace.value.main.chipProfiles :+ chipProfile
                      val updatedSubject = subject.copy(
                        genotypeRefs = subject.genotypeRefs :+ chipProfileUri,
                        meta = subject.meta.updated("genotypeRefs")
                      )
                      val updatedSamples = _workspace.value.main.samples.map { s =>
                        if (s.sampleAccession == sampleAccession) updatedSubject else s
                      }
                      val updatedContent = _workspace.value.main.copy(
                        samples = updatedSamples,
                        chipProfiles = updatedChipProfiles
                      )
                      _workspace.value = _workspace.value.copy(main = updatedContent)
                      saveWorkspace()

                      analysisInProgress.value = false
                      analysisProgress.value = "Import complete"
                      analysisProgressPercent.value = 1.0

                      println(s"[ViewModel] Imported chip data: ${parser.vendor}, ${summary.totalMarkersCalled} markers")
                      onComplete(Right(chipProfile))
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
    import java.security.MessageDigest
    import java.nio.file.Files
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
    findSubject(sampleAccession) match {
      case Some(subject) =>
        // Remove the profile from workspace
        val updatedChipProfiles = _workspace.value.main.chipProfiles.filterNot(_.atUri.contains(profileUri))

        // Remove from subject's genotypeRefs list
        val updatedSubject = subject.copy(
          genotypeRefs = subject.genotypeRefs.filterNot(_ == profileUri),
          meta = subject.meta.updated("genotypeRefs")
        )
        val updatedSamples = _workspace.value.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }

        val updatedContent = _workspace.value.main.copy(
          samples = updatedSamples,
          chipProfiles = updatedChipProfiles
        )
        _workspace.value = _workspace.value.copy(main = updatedContent)
        saveWorkspace()

        println(s"[ViewModel] Deleted chip profile $profileUri for $sampleAccession")
        Right(())

      case None =>
        Left(s"Subject not found: $sampleAccession")
    }
  }
}
