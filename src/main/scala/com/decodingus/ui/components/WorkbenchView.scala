package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.layout.{VBox, HBox, Priority, StackPane, BorderPane, Region}
import scalafx.scene.control.{Label, Button, ListView, SplitPane, Alert, ListCell, Tooltip, TextField}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import com.decodingus.workspace.model.{Workspace, Project, Biosample, SyncStatus, StrProfile}
import com.decodingus.workspace.WorkbenchViewModel
import scalafx.scene.control.Alert.AlertType
import scalafx.application.Platform
import scalafx.scene.control.ControlIncludes._
import scalafx.scene.control.ButtonType
import scalafx.scene.input.{MouseEvent, DragEvent, TransferMode, ClipboardContent}
import java.util.{Timer, TimerTask}

class WorkbenchView(val viewModel: WorkbenchViewModel) extends SplitPane {
  println(s"[DEBUG] WorkbenchView: Initializing WorkbenchView. ViewModel Projects: ${viewModel.projects.size}, ViewModel Samples: ${viewModel.samples.size}")

  // Track drag state to prevent click-on-drag from triggering navigation
  private var dragInProgress = false
  // Timer for delayed selection check
  private val selectionTimer = new Timer("SelectionTimer", true) // daemon thread
  private var pendingSelectionTask: Option[TimerTask] = None
  private val clickDelayMs = 150 // milliseconds to wait before applying selection

  // Use the shared DataFormat from ProjectDetailView companion object
  private val biosampleFormat = ProjectDetailView.biosampleFormat

  // Observable buffers for UI lists - now using filtered versions from ViewModel
  private val projectBuffer: ObservableBuffer[Project] = viewModel.filteredProjects
  private val sampleBuffer: ObservableBuffer[Biosample] = viewModel.filteredSamples

  println(s"[DEBUG] WorkbenchView: After binding buffers. projectBuffer size: ${projectBuffer.size}, sampleBuffer size: ${sampleBuffer.size}")

  // Detail view for right panel
  private val detailView = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new Label("Select an item to view details") { style = "-fx-font-size: 18px; -fx-font-weight: bold;" }
    )
  }
  VBox.setVgrow(detailView, Priority.Always)

  /** Renders the detail view for a selected subject with Edit/Delete actions */
  private def renderSubjectDetail(subject: Biosample): Unit = {
    detailView.children.clear()

    val editButton = new Button("Edit") {
      onAction = _ => handleEditSubject(subject)
    }

    val deleteButton = new Button("Delete") {
      style = "-fx-text-fill: #D32F2F;"
      onAction = _ => handleDeleteSubject(subject)
    }

    val actionButtons = new HBox(10) {
      padding = Insets(10, 0, 10, 0)
      children = Seq(editButton, deleteButton)
    }

    // Subject info section
    val infoSection = new VBox(5) {
      children = Seq(
        new Label(s"Accession: ${subject.sampleAccession}"),
        new Label(s"Sex: ${subject.sex.getOrElse("N/A")}"),
        new Label(s"Center: ${subject.centerName.getOrElse("N/A")}"),
        new Label(s"Description: ${subject.description.getOrElse("N/A")}"),
        new Label(s"Created At: ${subject.meta.createdAt.toLocalDate.toString}")
      )
    }

    // Haplogroup summary - format nicely with name and score
    def formatHaplogroup(result: Option[com.decodingus.workspace.model.HaplogroupResult]): String = {
      result match {
        case Some(h) =>
          val name = h.haplogroupName
          val derived = h.matchingSnps.map(n => s"+$n").getOrElse("")
          val ancestral = h.ancestralMatches.map(n => s"-$n").getOrElse("")
          if (derived.nonEmpty || ancestral.nonEmpty) s"$name ($derived/$ancestral)" else name
        case None => "—"
      }
    }

    val haplogroupText = subject.haplogroups match {
      case Some(h) =>
        val yDisplay = s"Y: ${formatHaplogroup(h.yDna)}"
        val mtDisplay = s"MT: ${formatHaplogroup(h.mtDna)}"
        s"$yDisplay    $mtDisplay"
      case None => "Haplogroups: Not analyzed"
    }

    val haplogroupLabel = new Label(haplogroupText) {
      style = "-fx-padding: 10 0 0 0; -fx-font-family: monospace; -fx-font-size: 14px; -fx-font-weight: bold;"
    }

    // Sequence data table with callbacks
    // Get sequence runs and alignments for this subject from the workspace
    val sequenceRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
    val allAlignments = viewModel.workspace.value.main.alignments

    val sequenceTable = new SequenceDataTable(
      viewModel = viewModel,
      subject = subject,
      sequenceRuns = sequenceRuns,
      alignments = allAlignments,
      onAnalyze = (index: Int) => handleAnalyzeSequenceData(subject.sampleAccession, index),
      onRemove = (index: Int) => handleRemoveSequenceData(subject.sampleAccession, index)
    )
    VBox.setVgrow(sequenceTable, Priority.Always)

    // Chip/Array data table
    val chipProfiles = viewModel.getChipProfilesForBiosample(subject.sampleAccession)
    val chipTable = new ChipDataTable(
      viewModel = viewModel,
      subject = subject,
      chipProfiles = chipProfiles,
      onRemove = (uri: String) => handleRemoveChipProfile(subject.sampleAccession, uri)
    )

    // STR profile table
    val strProfiles = viewModel.getStrProfilesForBiosample(subject.sampleAccession)
    val strTable = new StrProfileTable(
      viewModel = viewModel,
      subject = subject,
      strProfiles = strProfiles,
      onRemove = (uri: String) => handleRemoveStrProfile(subject.sampleAccession, uri)
    )

    detailView.children.addAll(
      new Label(s"Subject: ${subject.donorIdentifier}") { style = "-fx-font-size: 20px; -fx-font-weight: bold;" },
      actionButtons,
      infoSection,
      haplogroupLabel,
      sequenceTable,
      chipTable,
      strTable
    )
  }

  /** Handles triggering analysis for a sequence run */
  private def handleAnalyzeSequenceData(sampleAccession: String, index: Int): Unit = {
    viewModel.findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(index) match {
          case Some(sequenceRun) =>
            sequenceRun.files.headOption match {
              case Some(fileInfo) =>
                // Check if initial analysis has been run (has alignments)
                val runAlignments = viewModel.workspace.value.main.getAlignmentsForSequenceRun(sequenceRun)
                val hasAlignments = runAlignments.nonEmpty
                val hasMetrics = runAlignments.exists(_.metrics.isDefined)

                if (!hasAlignments) {
                  // Run initial analysis
                  showAnalysisChoiceDialog(sampleAccession, index, fileInfo.fileName, "initial")
                } else if (!hasMetrics) {
                  // Offer to run deep coverage analysis
                  showAnalysisChoiceDialog(sampleAccession, index, fileInfo.fileName, "wgs")
                } else {
                  // Both analyses complete - offer to re-run
                  showAnalysisChoiceDialog(sampleAccession, index, fileInfo.fileName, "both_complete")
                }
              case None =>
                new Alert(AlertType.Warning) {
                  title = "No File"
                  headerText = "No alignment file associated"
                  contentText = "Please add a BAM/CRAM file to this sequencing run."
                }.showAndWait()
            }
          case None =>
            println(s"[View] Sequence run not found at index $index")
        }
      case None =>
        println(s"[View] Subject not found: $sampleAccession")
    }
  }

  /** Shows a dialog to choose which analysis to run */
  private def showAnalysisChoiceDialog(sampleAccession: String, index: Int, fileName: String, state: String): Unit = {
    val (dialogHeader, dialogContent, options) = state match {
      case "initial" =>
        ("Run Initial Analysis",
         s"Analyze $fileName to detect platform, reference build, and collect library statistics.",
         Seq(("Run Initial Analysis", () => runInitialAnalysis(sampleAccession, index))))
      case "wgs" =>
        ("Run Deep Coverage Analysis",
         s"Initial analysis complete. Would you like to run WGS metrics analysis?\n\nThis will calculate detailed coverage statistics using GATK and may take several minutes for large genomes.",
         Seq(
           ("Run WGS Metrics", () => runWgsMetricsAnalysis(sampleAccession, index)),
           ("Re-run Initial Analysis", () => runInitialAnalysis(sampleAccession, index))
         ))
      case "both_complete" =>
        ("Analysis Complete",
         s"Both initial and WGS metrics analysis have been completed for $fileName.\n\nWould you like to re-run any analysis?",
         Seq(
           ("Re-run WGS Metrics", () => runWgsMetricsAnalysis(sampleAccession, index)),
           ("Re-run Initial Analysis", () => runInitialAnalysis(sampleAccession, index))
         ))
      case _ =>
        ("Analysis", "Choose an analysis to run.", Seq.empty)
    }

    if (options.size == 1) {
      // Single option - just confirm
      val confirm = new Alert(AlertType.Confirmation) {
        title = "Analysis"
        headerText = dialogHeader
        contentText = dialogContent
      }
      confirm.showAndWait() match {
        case Some(ButtonType.OK) => options.head._2()
        case _ =>
      }
    } else if (options.nonEmpty) {
      // Multiple options - use custom buttons
      val alert = new Alert(AlertType.Confirmation) {
        title = "Analysis Options"
        headerText = dialogHeader
        contentText = dialogContent
        buttonTypes = options.map(o => new ButtonType(o._1)) :+ ButtonType.Cancel
      }
      val result = alert.showAndWait()
      result.foreach { btn =>
        options.find(_._1 == btn.text).foreach(_._2())
      }
    }
  }

  /** Runs initial analysis with progress dialog */
  private def runInitialAnalysis(sampleAccession: String, index: Int): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      "Initial Analysis",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runInitialAnalysis(sampleAccession, index, {
      case Right(libraryStats) =>
        Platform.runLater {
          // Calculate mean read length from distribution
          val meanReadLength = if (libraryStats.lengthDistribution.nonEmpty) {
            val totalReads = libraryStats.lengthDistribution.values.sum.toDouble
            val weightedSum = libraryStats.lengthDistribution.map { case (len, count) => len.toLong * count }.sum
            f"${weightedSum / totalReads}%.1f bp"
          } else "N/A"

          // Calculate mean insert size from distribution
          val meanInsertSize = if (libraryStats.insertSizeDistribution.nonEmpty) {
            val totalPairs = libraryStats.insertSizeDistribution.values.sum.toDouble
            val weightedSum = libraryStats.insertSizeDistribution.map { case (size, count) => size * count }.sum
            f"${weightedSum / totalPairs}%.1f bp"
          } else "N/A"

          new Alert(AlertType.Information) {
            title = "Analysis Complete"
            headerText = "Initial Analysis Results"
            contentText = s"""Sample: ${libraryStats.sampleName}
                             |Platform: ${libraryStats.inferredPlatform}
                             |Instrument: ${libraryStats.mostFrequentInstrument}
                             |Reference: ${libraryStats.referenceBuild}
                             |Aligner: ${libraryStats.aligner}
                             |Mean Read Length: $meanReadLength
                             |Mean Insert Size: $meanInsertSize""".stripMargin
          }.showAndWait()
          // Refresh the detail view
          viewModel.selectedSubject.value.foreach(renderSubjectDetail)
        }
      case Left(error) =>
        Platform.runLater {
          new Alert(AlertType.Error) {
            title = "Analysis Failed"
            headerText = "Initial analysis encountered an error"
            contentText = error
          }.showAndWait()
        }
    })

    progressDialog.show()
  }

  /** Runs WGS metrics analysis with progress dialog */
  private def runWgsMetricsAnalysis(sampleAccession: String, index: Int): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      "WGS Metrics Analysis",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runWgsMetricsAnalysis(sampleAccession, index, {
      case Right(wgsMetrics) =>
        Platform.runLater {
          new Alert(AlertType.Information) {
            title = "Analysis Complete"
            headerText = "WGS Metrics Results"
            contentText = f"""Mean Coverage: ${wgsMetrics.meanCoverage}%.1fx
                             |Median Coverage: ${wgsMetrics.medianCoverage}%.1fx
                             |SD Coverage: ${wgsMetrics.sdCoverage}%.2f
                             |PCT 10x: ${wgsMetrics.pct10x * 100}%.1f%%
                             |PCT 20x: ${wgsMetrics.pct20x * 100}%.1f%%
                             |PCT 30x: ${wgsMetrics.pct30x * 100}%.1f%%
                             |Het SNP Sensitivity: ${wgsMetrics.hetSnpSensitivity}%.4f""".stripMargin
          }.showAndWait()
          // Refresh the detail view
          viewModel.selectedSubject.value.foreach(renderSubjectDetail)
        }
      case Left(error) =>
        Platform.runLater {
          new Alert(AlertType.Error) {
            title = "Analysis Failed"
            headerText = "WGS metrics analysis encountered an error"
            contentText = error
          }.showAndWait()
        }
    })

    progressDialog.show()
  }

  /** Handles removing a sequence data entry */
  private def handleRemoveSequenceData(sampleAccession: String, index: Int): Unit = {
    viewModel.removeSequenceData(sampleAccession, index)
  }

  /** Handles removing an STR profile */
  private def handleRemoveStrProfile(sampleAccession: String, profileUri: String): Unit = {
    viewModel.deleteStrProfile(sampleAccession, profileUri) match {
      case Right(_) =>
        // Refresh the detail view
        viewModel.selectedSubject.value.foreach(renderSubjectDetail)
      case Left(error) =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Could not remove STR profile"
          contentText = error
        }.showAndWait()
    }
  }

  /** Handles removing a chip profile */
  private def handleRemoveChipProfile(sampleAccession: String, profileUri: String): Unit = {
    viewModel.deleteChipProfile(sampleAccession, profileUri) match {
      case Right(_) =>
        // Refresh the detail view
        viewModel.selectedSubject.value.foreach(renderSubjectDetail)
      case Left(error) =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Could not remove chip profile"
          contentText = error
        }.showAndWait()
    }
  }

  /** Renders the detail view for a selected project */
  private def renderProjectDetail(project: Project): Unit = {
    detailView.children.clear()
    val projectDetailView = new ProjectDetailView(
      viewModel = viewModel,
      project = project,
      onEdit = handleEditProject,
      onDelete = handleDeleteProject
    )
    VBox.setVgrow(projectDetailView, Priority.Always)
    detailView.children.add(projectDetailView)
  }

  /** Handles the Edit Project action */
  private def handleEditProject(project: Project): Unit = {
    val dialog = new EditProjectDialog(project)
    val result = dialog.showAndWait().asInstanceOf[Option[Option[Project]]]

    result match {
      case Some(Some(updatedProject)) =>
        viewModel.updateProject(updatedProject)
        // Refresh the detail view with updated project
        viewModel.findProject(updatedProject.projectName).foreach(renderProjectDetail)
      case _ => // User cancelled
    }
  }

  /** Handles the Delete Project action with confirmation */
  private def handleDeleteProject(project: Project): Unit = {
    val confirmDialog = new Alert(AlertType.Confirmation) {
      title = "Delete Project"
      headerText = s"Delete ${project.projectName}?"
      contentText = "This action cannot be undone. The project will be removed but subjects will remain in the workspace."
    }

    val result = confirmDialog.showAndWait()
    result match {
      case Some(ButtonType.OK) =>
        viewModel.deleteProject(project.projectName)
        renderEmptyDetail("Select an item to view details")
      case _ => // User cancelled
    }
  }

  /** Renders the empty state when nothing is selected */
  private def renderEmptyDetail(message: String): Unit = {
    detailView.children.clear()
    detailView.children.add(
      new Label(message) { style = "-fx-font-size: 18px; -fx-font-weight: bold;" }
    )
  }

  /** Handles the Edit Subject action */
  private def handleEditSubject(subject: Biosample): Unit = {
    val dialog = new EditSubjectDialog(subject)
    val result = dialog.showAndWait().asInstanceOf[Option[Option[Biosample]]]

    result match {
      case Some(Some(updatedBiosample)) =>
        viewModel.updateSubject(updatedBiosample)
      case _ => // User cancelled
    }
  }

  /** Handles the Delete Subject action with confirmation */
  private def handleDeleteSubject(subject: Biosample): Unit = {
    val confirmDialog = new Alert(AlertType.Confirmation) {
      title = "Delete Subject"
      headerText = s"Delete ${subject.donorIdentifier}?"
      contentText = "This action cannot be undone. All associated sequence data and analysis results will be removed."
    }

    val result = confirmDialog.showAndWait()
    result match {
      case Some(ButtonType.OK) =>
        viewModel.deleteSubject(subject.sampleAccession)
      case _ => // User cancelled
    }
  }

  // Unified detail view rendering based on selection state
  // This prevents race conditions between project and subject selection
  private def updateDetailView(): Unit = {
    Platform.runLater {
      (viewModel.selectedProject.value, viewModel.selectedSubject.value) match {
        case (Some(project), _) =>
          // Project takes precedence when selected
          renderProjectDetail(project)
        case (None, Some(subject)) =>
          renderSubjectDetail(subject)
        case (None, None) =>
          renderEmptyDetail("Select an item to view details")
      }
    }
  }

  // Listen to ViewModel's selectedSubject changes to update detailView
  viewModel.selectedSubject.onChange { (_, _, _) =>
    updateDetailView()
  }

  // Listen to ViewModel's selectedProject changes to update detailView
  viewModel.selectedProject.onChange { (_, _, _) =>
    updateDetailView()
  }

  // Left Panel - Navigation
  private val projectList = new ListView[Project]() {
    items = projectBuffer // Explicitly set items
    vgrow = Priority.Always
    prefHeight = 200 // Initial height for projects
    cellFactory = { (v: ListView[Project]) =>
      new ListCell[Project] {
        item.onChange { (_, _, newProject) =>
          text = if (newProject != null) newProject.projectName else null
        }
      }
    }
  }
  // UI to ViewModel sync
  projectList.selectionModel().selectedItem.onChange { (_, _, newProject) =>
    if (newProject != null) {
      // Clear subject selection first (both UI and ViewModel) to avoid race conditions
      if (viewModel.selectedSubject.value.isDefined) {
        viewModel.selectedSubject.value = None
      }
      sampleList.selectionModel().clearSelection()
      // Then set project selection
      viewModel.selectedProject.value = Some(newProject)
    } else if (viewModel.selectedProject.value.isDefined && projectList.selectionModel().getSelectedItem == null) {
      // Clear ViewModel selection if UI selection is cleared manually
      viewModel.selectedProject.value = None
    }
  }
  // ViewModel to UI sync
  viewModel.selectedProject.onChange { (_, _, newViewModelProjectOpt) =>
    if (newViewModelProjectOpt.isDefined && projectList.selectionModel().getSelectedItem != newViewModelProjectOpt.getOrElse(null)) {
      projectList.selectionModel().select(newViewModelProjectOpt.get)
    } else if (newViewModelProjectOpt.isEmpty && projectList.selectionModel().getSelectedItem != null) {
      projectList.selectionModel().clearSelection()
    }
  }

  private val sampleList = new ListView[Biosample]() {
    items = sampleBuffer // Explicitly set items
    vgrow = Priority.Always
    cellFactory = { (v: ListView[Biosample]) =>
      new ListCell[Biosample] {
        item.onChange { (_, _, newBiosample) =>
          text = if (newBiosample != null) s"${newBiosample.donorIdentifier} (${newBiosample.sampleAccession.take(8)}...)" else null
        }

        // On mouse press, schedule delayed selection check
        onMousePressed = (_: MouseEvent) => {
          dragInProgress = false
          // Cancel any existing pending task
          pendingSelectionTask.foreach(_.cancel())

          Option(item.value).foreach { biosample =>
            val task = new TimerTask {
              override def run(): Unit = {
                // Check if drag started during the delay
                if (!dragInProgress) {
                  Platform.runLater {
                    // Clear project selection first
                    if (viewModel.selectedProject.value.isDefined) {
                      viewModel.selectedProject.value = None
                    }
                    projectList.selectionModel().clearSelection()
                    // Then set subject selection
                    viewModel.selectedSubject.value = Some(biosample)
                  }
                }
                pendingSelectionTask = None
              }
            }
            pendingSelectionTask = Some(task)
            selectionTimer.schedule(task, clickDelayMs)
          }
        }

        // Drag source - enable dragging subjects to project members lists
        onDragDetected = (event: MouseEvent) => {
          Option(item.value).foreach { biosample =>
            dragInProgress = true
            // Cancel pending selection
            pendingSelectionTask.foreach(_.cancel())
            pendingSelectionTask = None

            val db = startDragAndDrop(TransferMode.Move)
            val content = new ClipboardContent()
            content.put(biosampleFormat, biosample.sampleAccession)
            content.putString(biosample.sampleAccession)
            db.setContent(content)
            event.consume()
          }
        }

        // Reset drag flag when drag completes
        onDragDone = (_: DragEvent) => {
          dragInProgress = false
        }
      }
    }
  }
  // UI to ViewModel sync - only handle deselection
  sampleList.selectionModel().selectedItem.onChange { (_, _, newBiosample) =>
    if (newBiosample == null && viewModel.selectedSubject.value.isDefined) {
      // Clear ViewModel selection if UI selection is cleared
      viewModel.selectedSubject.value = None
    }
  }
  // ViewModel to UI sync
  viewModel.selectedSubject.onChange { (_, _, newViewModelSubjectOpt) =>
    if (newViewModelSubjectOpt.isDefined && sampleList.selectionModel().getSelectedItem != newViewModelSubjectOpt.getOrElse(null)) {
      sampleList.selectionModel().select(newViewModelSubjectOpt.get)
    } else if (newViewModelSubjectOpt.isEmpty && sampleList.selectionModel().getSelectedItem != null) {
      sampleList.selectionModel().clearSelection()
    }
  }

  private val newProjectButton = new Button("New Project") {
    onAction = _ => {
      val dialog = new AddProjectDialog()
      val result = dialog.showAndWait().asInstanceOf[Option[Option[Project]]]

      result match {
        case Some(Some(newProject)) =>
          viewModel.addProject(newProject)
        case _ => // User cancelled
      }
    }
  }

  // Filter controls
  private val projectFilterField = new TextField() {
    promptText = "Filter projects..."
    prefWidth = 150
  }
  projectFilterField.text.bindBidirectional(viewModel.projectFilter)

  private val subjectFilterField = new TextField() {
    promptText = "Filter subjects..."
    prefWidth = 150
  }
  subjectFilterField.text.bindBidirectional(viewModel.subjectFilter)

  private val addSampleButton = new Button("Add Subject") {
    onAction = _ => {
      val dialog = new AddSubjectDialog()
      val result = dialog.showAndWait().asInstanceOf[Option[Option[Biosample]]]

      result match {
        case Some(Some(newBiosample)) =>
          viewModel.addSubject(newBiosample) // Delegate to ViewModel
        case _ => // User cancelled or closed dialog
      }
    }
  }

  private val saveButton = new Button("Save Workspace") {
    styleClass.add("button-primary")
    onAction = _ => viewModel.saveWorkspace() // Delegate to ViewModel
  }

  // Sync status indicator - initialize with default values directly
  private val syncStatusLabel = new Label("○ Local Only") {
    style = "-fx-font-size: 12px; -fx-padding: 2 6 2 6; -fx-text-fill: #9E9E9E;"
    tooltip = Tooltip("Using local storage only (not logged in)")
  }

  private def updateSyncStatus(status: SyncStatus): Unit = {
    status match {
      case SyncStatus.Synced =>
        syncStatusLabel.text = "✓ Synced"
        syncStatusLabel.style = "-fx-font-size: 12px; -fx-padding: 2 6 2 6; -fx-text-fill: #4CAF50;"
        syncStatusLabel.tooltip = Tooltip("Workspace is synced")
      case SyncStatus.Syncing =>
        syncStatusLabel.text = "↻ Syncing..."
        syncStatusLabel.style = "-fx-font-size: 12px; -fx-padding: 2 6 2 6; -fx-text-fill: #2196F3;"
        syncStatusLabel.tooltip = Tooltip("Syncing with PDS...")
      case SyncStatus.Pending =>
        syncStatusLabel.text = "● Pending"
        syncStatusLabel.style = "-fx-font-size: 12px; -fx-padding: 2 6 2 6; -fx-text-fill: #FF9800;"
        syncStatusLabel.tooltip = Tooltip("Changes pending sync")
      case SyncStatus.Error =>
        syncStatusLabel.text = "⚠ Sync Error"
        syncStatusLabel.style = "-fx-font-size: 12px; -fx-padding: 2 6 2 6; -fx-text-fill: #F44336;"
        syncStatusLabel.tooltip = Tooltip(s"Sync failed: ${viewModel.lastSyncError.value}")
      case SyncStatus.Offline =>
        syncStatusLabel.text = "○ Local Only"
        syncStatusLabel.style = "-fx-font-size: 12px; -fx-padding: 2 6 2 6; -fx-text-fill: #9E9E9E;"
        syncStatusLabel.tooltip = Tooltip("Using local storage only (not logged in)")
    }
  }

  // Listen to sync status changes
  viewModel.syncStatus.onChange { (_, _, newStatus) =>
    Platform.runLater {
      updateSyncStatus(newStatus)
    }
  }

  // Update tooltip when error message changes
  viewModel.lastSyncError.onChange { (_, _, newError) =>
    Platform.runLater {
      if (viewModel.syncStatus.value == SyncStatus.Error && newError.nonEmpty) {
        syncStatusLabel.tooltip = Tooltip(s"Sync failed: $newError")
      }
    }
  }

  // Apply initial status after label is constructed
  updateSyncStatus(viewModel.syncStatus.value)

  private val leftPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new HBox(10) {
        alignment = Pos.CenterLeft
        children = Seq(
          new Label("Projects:") { style = "-fx-font-weight: bold;" },
          new Region { HBox.setHgrow(this, Priority.Always) },
          projectFilterField
        )
      },
      projectList,
      newProjectButton,
      new HBox(10) {
        alignment = Pos.CenterLeft
        children = Seq(
          new Label("Subjects:") { style = "-fx-font-weight: bold;" },
          new Region { HBox.setHgrow(this, Priority.Always) },
          subjectFilterField
        )
      },
      sampleList,
      addSampleButton,
      new HBox(10) {
        alignment = Pos.CenterLeft
        children = Seq(
          syncStatusLabel,
          new Region { HBox.setHgrow(this, Priority.Always) }, // Spacer
          saveButton
        )
      }
    )
  }
  SplitPane.setResizableWithParent(leftPanel, false) // Make left panel not resize with parent by default

  // Right Panel - Details/Content Area
  private val rightPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(detailView) // rightPanel now contains the dynamic detailView
  }
  VBox.setVgrow(rightPanel, Priority.Always) // Allow right panel to grow vertically

  // Set the items of the SplitPane
  items.addAll(leftPanel, rightPanel)
  dividerPositions = 0.25 // Initial divider position

  // Listen for reference download prompts from the ViewModel
  viewModel.pendingReferenceDownload.onChange { (_, _, request) =>
    request match {
      case viewModel.PendingDownload(build, url, sizeMB, onConfirm, onCancel) =>
        Platform.runLater {
          val dialog = new ReferenceDownloadPromptDialog(build, url, sizeMB)
          dialog.showAndWait() match {
            case Some(ReferenceDownloadPromptDialog.Result.Download) =>
              onConfirm()
            case Some(ReferenceDownloadPromptDialog.Result.Configure) =>
              // Open settings dialog
              val configDialog = new ReferenceConfigDialog()
              configDialog.showAndWait()
              onCancel() // Cancel the current operation - user can retry after configuring
            case _ =>
              onCancel()
          }
        }
      case viewModel.NoDownloadPending =>
        // Nothing to do
    }
  }
}