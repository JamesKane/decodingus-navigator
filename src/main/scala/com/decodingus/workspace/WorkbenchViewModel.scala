package com.decodingus.workspace

import com.decodingus.analysis.{ArtifactContext, CallableLociProcessor, CallableLociResult, HaplogroupProcessor, LibraryStatsProcessor, UnifiedMetricsProcessor, WgsMetricsProcessor}
import com.decodingus.client.DecodingUsClient
import com.decodingus.haplogroup.tree.{TreeType, TreeProviderType}
import com.decodingus.auth.User
import com.decodingus.config.{FeatureToggles, UserPreferencesService, ReferenceConfigService}
import com.decodingus.model.{LibraryStats, WgsMetrics}
import com.decodingus.refgenome.{ReferenceGateway, ReferenceResolveResult}
import com.decodingus.workspace.model.{Workspace, Project, Biosample, WorkspaceContent, SyncStatus, SequenceRun, Alignment, AlignmentMetrics, ContigMetrics, FileInfo, HaplogroupAssignments, HaplogroupResult => WorkspaceHaplogroupResult, RecordMeta, StrProfile, ChipProfile}
import com.decodingus.haplogroup.model.{HaplogroupResult => AnalysisHaplogroupResult}
import com.decodingus.workspace.services.{WorkspaceOperations, AnalysisCoordinator, AnalysisProgress, SyncService, SyncResult, FingerprintMatchService, FingerprintMatchResult}
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

  // --- Service Instances ---
  private val workspaceOps = new WorkspaceOperations()
  private val analysisCoordinator = new AnalysisCoordinator()
  private val syncService = new SyncService(workspaceService)
  private val fingerprintMatchService = new FingerprintMatchService()

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
    syncService.syncFromPds(
      currentUser.value,
      _workspace.value,
      status => Platform.runLater { syncStatus.value = status }
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
      status => Platform.runLater { syncStatus.value = status }
    ).foreach {
      case Left(error) =>
        Platform.runLater { lastSyncError.value = error }
      case Right(SyncResult.Error(msg)) =>
        Platform.runLater { lastSyncError.value = msg }
      case Right(_) =>
        Platform.runLater { lastSyncError.value = "" }
    }
  }

  private def emptyWorkspace: Workspace = Workspace.empty

  /** Gets the current workspace state for service calls */
  private def currentState: WorkspaceState = WorkspaceState(_workspace.value)

  /** Applies a new workspace state and saves */
  private def applyState(newState: WorkspaceState): Unit = {
    _workspace.value = newState.workspace
    saveWorkspace()
  }

  // --- Subject CRUD Operations ---

  /** Creates a new subject and adds it to the workspace */
  def addSubject(newBiosample: Biosample): Unit = {
    val userDid = currentUser.value.map(_.did)
    val (newState, enrichedBiosample) = workspaceOps.addSubject(currentState, newBiosample, userDid)
    applyState(newState)
    selectedSubject.value = Some(enrichedBiosample)
  }

  /** Updates an existing subject identified by sampleAccession */
  def updateSubject(updatedBiosample: Biosample): Unit = {
    val newState = workspaceOps.updateSubject(currentState, updatedBiosample)
    applyState(newState)
    // Update selection to reflect changes
    selectedSubject.value = newState.workspace.main.samples.find(_.sampleAccession == updatedBiosample.sampleAccession)
  }

  /** Internal: Updates a subject without modifying meta (used when meta is already updated) */
  private def updateSubjectDirect(updatedBiosample: Biosample): Unit = {
    val newState = workspaceOps.updateSubjectDirect(currentState, updatedBiosample)
    applyState(newState)
    selectedSubject.value = Some(updatedBiosample)
  }

  /** Deletes a subject identified by sampleAccession */
  def deleteSubject(sampleAccession: String): Unit = {
    val newState = workspaceOps.deleteSubject(currentState, sampleAccession)
    applyState(newState)

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

  def addProject(newProject: Project): Unit = {
    val userDid = currentUser.value.map(_.did)
    val (newState, enrichedProject) = workspaceOps.addProject(currentState, newProject, userDid)
    applyState(newState)
    selectedProject.value = Some(enrichedProject)
  }

  /** Updates an existing project identified by projectName */
  def updateProject(updatedProject: Project): Unit = {
    val newState = workspaceOps.updateProject(currentState, updatedProject)
    applyState(newState)
    selectedProject.value = newState.workspace.main.projects.find(_.projectName == updatedProject.projectName)
  }

  /** Deletes a project by projectName */
  def deleteProject(projectName: String): Unit = {
    val newState = workspaceOps.deleteProject(currentState, projectName)
    applyState(newState)

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

  /** Backfills atUri for any samples/projects that were created while logged out */
  private def backfillAtUris(did: String): Unit = {
    val newState = workspaceOps.backfillAtUris(currentState, did)
    if (newState.workspace != currentState.workspace) {
      applyState(newState)
      println(s"[ViewModel] Backfilled atUri for samples/projects after login")
    }
  }

  /** Adds a subject (by accession) to a project's members list */
  def addSubjectToProject(projectName: String, sampleAccession: String): Boolean = {
    workspaceOps.addSubjectToProject(currentState, projectName, sampleAccession) match {
      case Right(newState) =>
        applyState(newState)
        // Update selection to the modified project
        selectedProject.value = newState.workspace.main.projects.find(_.projectName == projectName)
        true
      case Left(error) =>
        println(s"[ViewModel] $error")
        false
    }
  }

  /** Removes a subject (by accession) from a project's members list */
  def removeSubjectFromProject(projectName: String, sampleAccession: String): Boolean = {
    workspaceOps.removeSubjectFromProject(currentState, projectName, sampleAccession) match {
      case Right(newState) =>
        applyState(newState)
        // Update selection to the modified project
        selectedProject.value = newState.workspace.main.projects.find(_.projectName == projectName)
        true
      case Left(error) =>
        println(s"[ViewModel] $error")
        false
    }
  }

  /** Internal: Updates a project without modifying meta (used when meta is already updated) */
  private def updateProjectDirect(updatedProject: Project): Unit = {
    val newState = workspaceOps.updateProjectDirect(currentState, updatedProject)
    applyState(newState)
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
   */
  def addSequenceRunFromFile(sampleAccession: String, fileInfo: FileInfo): Int = {
    workspaceOps.addSequenceRunFromFile(currentState, sampleAccession, fileInfo) match {
      case Right((newState, _, newIndex)) =>
        applyState(newState)
        newIndex
      case Left(error) =>
        println(s"[ViewModel] $error")
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
    // Step 1: Check for exact duplicate by checksum
    onProgress("Checking for duplicates...", 0.05)

    findSubject(sampleAccession) match {
      case None =>
        onComplete(Left(s"Subject $sampleAccession not found"))
        return
      case Some(subject) =>
        val sequenceRuns = _workspace.value.main.getSequenceRunsForBiosample(subject)
        val existingChecksums = sequenceRuns.flatMap(_.files.flatMap(_.checksum)).toSet
        if (fileInfo.checksum.exists(existingChecksums.contains)) {
          onComplete(Left("Duplicate file - this file has already been added"))
          return
        }
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
                    import com.decodingus.ui.components.{FingerprintMatchDialog, GroupTogether, KeepSeparate, FingerprintMatchDecision}
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
                        println(s"[ViewModel] User confirmed LOW confidence match - adding ${libraryStats.referenceBuild} alignment to existing run")
                        Some((idx, existingRun, false))
                      case Some(KeepSeparate) =>
                        println(s"[ViewModel] User chose to keep separate - creating new sequence run")
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
                    println(s"[ViewModel] Fingerprint match found ($confidence confidence) - adding ${libraryStats.referenceBuild} alignment to existing run")
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

              // Update workspace
              val updatedSequenceRuns = _workspace.value.main.sequenceRuns.map { sr =>
                if (sr.atUri == seqRun.atUri) updatedSeqRun else sr
              }
              val updatedAlignments = _workspace.value.main.alignments :+ newAlignment
              val updatedContent = _workspace.value.main.copy(
                sequenceRuns = updatedSequenceRuns,
                alignments = updatedAlignments
              )
              _workspace.value = _workspace.value.copy(main = updatedContent)
              saveWorkspace()

              lastLibraryStats.value = Some(libraryStats)
              analysisInProgress.value = false
              analysisProgress.value = if (isNewRun) "Analysis complete" else s"Added ${libraryStats.referenceBuild} alignment to existing run"
              analysisProgressPercent.value = 1.0
              onProgress("Complete", 1.0)
              onComplete(Right((resultIndex, libraryStats)))
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

  /** Gets a specific sequence run by index for a subject */
  def getSequenceRun(sampleAccession: String, index: Int): Option[SequenceRun] = {
    findSubject(sampleAccession).flatMap { subject =>
      _workspace.value.main.getSequenceRunsForBiosample(subject).lift(index)
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
              println(s"[ViewModel] Updated facility for instrument $instrumentId: ${labInfo.labName}")
            }
          }
        }
      case Success(None) =>
        println(s"[ViewModel] No facility found for instrument ID: $instrumentId")
      case Failure(error) =>
        println(s"[ViewModel] Failed to lookup facility for instrument $instrumentId: ${error.getMessage}")
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
   * @param fingerprint The computed fingerprint to match
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
   *
   * @param sequenceRunIndex Index of the sequence run to add alignment to
   * @param newAlignment The new alignment to add
   * @param fileInfo The file info for the new alignment
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
        val alignUri = newAlignment.atUri.getOrElse("")

        // Add file to sequence run if not already present
        val updatedFiles = if (seqRun.files.exists(_.checksum == fileInfo.checksum)) {
          seqRun.files
        } else {
          seqRun.files :+ fileInfo
        }

        // Add alignment ref if not already present
        val updatedAlignmentRefs = if (seqRun.alignmentRefs.contains(alignUri)) {
          seqRun.alignmentRefs
        } else {
          seqRun.alignmentRefs :+ alignUri
        }

        val updatedSeqRun = seqRun.copy(
          files = updatedFiles,
          alignmentRefs = updatedAlignmentRefs,
          meta = seqRun.meta.updated("alignment-added")
        )

        // Update workspace
        val updatedSequenceRuns = _workspace.value.main.sequenceRuns.map { sr =>
          if (sr.atUri == seqRun.atUri) updatedSeqRun else sr
        }
        val updatedAlignments = _workspace.value.main.alignments :+ newAlignment

        val updatedContent = _workspace.value.main.copy(
          sequenceRuns = updatedSequenceRuns,
          alignments = updatedAlignments
        )
        _workspace.value = _workspace.value.copy(main = updatedContent)
        saveWorkspace()

        println(s"[ViewModel] Added ${newAlignment.referenceBuild} alignment to existing sequence run")
      }
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
                        testType = inferTestType(libraryStats),
                        libraryLayout = Some(if (libraryStats.pairedReads > libraryStats.readCount / 2) "Paired-End" else "Single-End"),
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

  // --- Callable Loci Analysis ---

  // Store last callable loci result
  val lastCallableLociResult: ObjectProperty[Option[CallableLociResult]] = ObjectProperty(None)

  /**
   * Runs callable loci analysis for a sequencing run.
   * Analyzes each contig to determine callable vs non-callable regions.
   * Requires that initial analysis has already been run.
   *
   * @param sampleAccession The subject's accession ID
   * @param sequenceRunIndex The index of the sequence run to analyze
   * @param onComplete Callback when analysis completes, returns result and artifact directory path
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
                    println(s"[ViewModel] Starting callable loci analysis for ${fileInfo.fileName}")

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

                          // Update workspace
                          val updatedAlignments = _workspace.value.main.alignments.map { a =>
                            if (a.atUri == alignment.atUri) updatedAlignment else a
                          }
                          val updatedContent = _workspace.value.main.copy(alignments = updatedAlignments)
                          _workspace.value = _workspace.value.copy(main = updatedContent)
                          saveWorkspace()

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
                              // Determine source type based on test type
                              val sourceType = seqRun.testType match {
                                case t if t.startsWith("BIGY") || t.contains("Y_ELITE") || t.contains("Y_PRIME") => "bigy"
                                case _ => "wgs"
                              }

                              // Convert analysis result to workspace model with full provenance
                              val workspaceResult = WorkspaceHaplogroupResult(
                                haplogroupName = topResult.name,
                                score = topResult.score,
                                matchingSnps = Some(topResult.matchingSnps),
                                mismatchingSnps = Some(topResult.mismatchingSnps),
                                ancestralMatches = Some(topResult.ancestralMatches),
                                treeDepth = Some(topResult.depth),
                                lineagePath = None,
                                source = Some(sourceType),
                                sourceRef = seqRun.atUri,
                                treeProvider = Some(treeProviderType.toString.toLowerCase),
                                treeVersion = None,
                                analyzedAt = Some(java.time.Instant.now())
                              )

                              // Update the consensus result in HaplogroupAssignments
                              val currentAssignments = subject.haplogroups.getOrElse(HaplogroupAssignments())
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

  // --- STR Profile CRUD Operations ---

  /**
   * Adds a new STR profile for a biosample.
   * Supports multiple profiles per subject (e.g., from different vendors like FTDNA and YSEQ).
   * Returns the URI of the new profile, or an error message.
   */
  def addStrProfile(sampleAccession: String, profile: StrProfile): Either[String, String] = {
    workspaceOps.addStrProfile(currentState, sampleAccession, profile) match {
      case Right((newState, profileUri)) =>
        applyState(newState)
        println(s"[ViewModel] Added STR profile with ${profile.markers.size} markers for $sampleAccession (provider: ${profile.importedFrom.getOrElse("unknown")})")
        Right(profileUri)
      case Left(error) =>
        Left(error)
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
        println(s"[ViewModel] Deleted STR profile $profileUri for $sampleAccession")
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
    workspaceOps.deleteChipProfile(currentState, sampleAccession, profileUri) match {
      case Right(newState) =>
        applyState(newState)
        println(s"[ViewModel] Deleted chip profile $profileUri for $sampleAccession")
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
   * @param profileUri The AT URI of the chip profile to analyze
   * @param panelType AIMs (quick) or GenomeWide (detailed)
   * @param onComplete Callback when analysis completes
   */
  def runChipAncestryAnalysis(
    sampleAccession: String,
    profileUri: String,
    panelType: com.decodingus.ancestry.model.AncestryPanelType,
    onComplete: Either[String, com.decodingus.ancestry.model.AncestryResult] => Unit
  ): Unit = {
    import com.decodingus.genotype.parser.ChipDataParser
    import com.decodingus.genotype.processor.{ChipDataProcessor, ChipAncestryAdapter}

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
                                val updatedSamples = _workspace.value.main.samples.map { s =>
                                  if (s.sampleAccession == sampleAccession) updatedSubject else s
                                }
                                val updatedContent = _workspace.value.main.copy(samples = updatedSamples)
                                _workspace.value = _workspace.value.copy(main = updatedContent)
                                saveWorkspace()

                                println(s"[ViewModel] Ancestry analysis complete for $sampleAccession: " +
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
   * @param profileUri The AT URI of the chip profile to analyze
   * @param treeType Y-DNA or MT-DNA tree type
   * @param onComplete Callback when analysis completes
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

                                // Update subject with haplogroup result (with full provenance)
                                val currentHaplogroups = subject.haplogroups.getOrElse(
                                  com.decodingus.workspace.model.HaplogroupAssignments()
                                )

                                val newHaplogroupResult = com.decodingus.workspace.model.HaplogroupResult(
                                  haplogroupName = haplogroupResult.topHaplogroup,
                                  score = haplogroupResult.results.headOption.map(_.score).getOrElse(0.0),
                                  matchingSnps = Some(haplogroupResult.snpsMatched),
                                  treeDepth = haplogroupResult.results.headOption.map(_.depth),
                                  source = Some("chip"),
                                  sourceRef = profile.atUri,
                                  treeProvider = Some(treeProviderName),
                                  treeVersion = None,
                                  analyzedAt = Some(java.time.Instant.now())
                                )

                                // Update the consensus result in HaplogroupAssignments
                                val updatedHaplogroups = treeType match {
                                  case com.decodingus.haplogroup.tree.TreeType.YDNA =>
                                    currentHaplogroups.copy(yDna = Some(newHaplogroupResult))
                                  case com.decodingus.haplogroup.tree.TreeType.MTDNA =>
                                    currentHaplogroups.copy(mtDna = Some(newHaplogroupResult))
                                }

                                val updatedSubject = subject.copy(
                                  haplogroups = Some(updatedHaplogroups),
                                  meta = subject.meta.updated("haplogroups")
                                )
                                val updatedSamples = _workspace.value.main.samples.map { s =>
                                  if (s.sampleAccession == sampleAccession) updatedSubject else s
                                }
                                val updatedContent = _workspace.value.main.copy(samples = updatedSamples)
                                _workspace.value = _workspace.value.copy(main = updatedContent)
                                saveWorkspace()

                                println(s"[ViewModel] $treeLabel haplogroup analysis complete for $sampleAccession: " +
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
   * @param sampleAccession The subject's accession ID
   * @param sequenceRunIndex Index of the sequence run within the subject
   * @param onComplete Callback when analysis completes
   */
  def runMultipleMetricsAnalysis(
    sampleAccession: String,
    sequenceRunIndex: Int,
    onComplete: Either[String, com.decodingus.analysis.ReadMetrics] => Unit
  ): Unit = {
    import com.decodingus.analysis.{UnifiedMetricsProcessor, ArtifactContext}

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

                          println(s"[ViewModel] ReadMetrics complete for $sampleAccession: " +
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
}
