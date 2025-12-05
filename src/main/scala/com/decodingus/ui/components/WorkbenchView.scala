package com.decodingus.ui.components

import scalafx.scene.layout.{VBox, HBox, Priority, StackPane, BorderPane, Region}
import scalafx.scene.control.{Label, Button, ListView, SplitPane, Alert, ListCell, Tooltip}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import com.decodingus.workspace.model.{Workspace, Project, Biosample, SyncStatus}
import com.decodingus.workspace.WorkbenchViewModel
import scalafx.scene.control.Alert.AlertType
import scalafx.application.Platform
import scalafx.scene.control.ControlIncludes._
import scalafx.scene.control.ButtonType

class WorkbenchView(val viewModel: WorkbenchViewModel) extends SplitPane {
  println(s"[DEBUG] WorkbenchView: Initializing WorkbenchView. ViewModel Projects: ${viewModel.projects.size}, ViewModel Samples: ${viewModel.samples.size}")

  // Observable buffers for UI lists - now directly from ViewModel
  private val projectBuffer: ObservableBuffer[Project] = viewModel.projects
  private val sampleBuffer: ObservableBuffer[Biosample] = viewModel.samples

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
        new Label(s"Created At: ${subject.createdAt.map(_.toLocalDate.toString).getOrElse("N/A")}")
      )
    }

    // Haplogroup summary
    val haplogroupLabel = new Label(
      s"Haplogroups: ${subject.haplogroups.map(h => s"Y: ${h.yDna.getOrElse("—")}, MT: ${h.mtDna.getOrElse("—")}").getOrElse("Not analyzed")}"
    ) {
      style = "-fx-padding: 10 0 0 0;"
    }

    // Sequence data table with callbacks
    val sequenceTable = new SequenceDataTable(
      viewModel = viewModel,
      subject = subject,
      onAnalyze = (index: Int) => handleAnalyzeSequenceData(subject.sampleAccession, index),
      onRemove = (index: Int) => handleRemoveSequenceData(subject.sampleAccession, index)
    )
    VBox.setVgrow(sequenceTable, Priority.Always)

    detailView.children.addAll(
      new Label(s"Subject: ${subject.donorIdentifier}") { style = "-fx-font-size: 20px; -fx-font-weight: bold;" },
      actionButtons,
      infoSection,
      haplogroupLabel,
      sequenceTable
    )
  }

  /** Handles triggering analysis for a sequence data entry */
  private def handleAnalyzeSequenceData(sampleAccession: String, index: Int): Unit = {
    viewModel.findSubject(sampleAccession).flatMap(_.sequenceData.lift(index)) match {
      case Some(seqData) =>
        seqData.files.headOption match {
          case Some(fileInfo) =>
            // For now, show a placeholder - this will be wired to actual analysis
            new Alert(AlertType.Information) {
              title = "Analysis"
              headerText = s"Starting analysis for ${fileInfo.fileName}"
              contentText = s"Platform: ${seqData.platformName}\nTest Type: ${seqData.testType}\n\nAnalysis functionality will be implemented next."
            }.showAndWait()
            // TODO: Call viewModel.analyzeLibraryStats or viewModel.analyzeDeepCoverage
          case None =>
            new Alert(AlertType.Warning) {
              title = "No File"
              headerText = "No alignment file associated"
              contentText = "Please add a BAM/CRAM file to this sequencing run."
            }.showAndWait()
        }
      case None =>
        println(s"[View] Sequence data not found at index $index")
    }
  }

  /** Handles removing a sequence data entry */
  private def handleRemoveSequenceData(sampleAccession: String, index: Int): Unit = {
    viewModel.removeSequenceData(sampleAccession, index)
  }

  /** Renders the detail view for a selected project */
  private def renderProjectDetail(project: Project): Unit = {
    detailView.children.clear()
    detailView.children.addAll(
      new Label(s"Project: ${project.projectName}") { style = "-fx-font-size: 20px; -fx-font-weight: bold;" },
      new Label(s"Description: ${project.description.getOrElse("N/A")}"),
      new Label(s"Administrator: ${project.administrator}"),
      new Label(s"Members: ${project.members.size} subjects")
    )
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

  // Listen to ViewModel's selectedSubject changes to update detailView
  viewModel.selectedSubject.onChange { (_, _, newSubjectOpt) =>
    Platform.runLater {
      newSubjectOpt match {
        case Some(subject) => renderSubjectDetail(subject)
        case None =>
          // Only show empty state if no project is selected either
          if (viewModel.selectedProject.value.isEmpty) {
            renderEmptyDetail("Select an item to view details")
          }
      }
    }
  }

  // Listen to ViewModel's selectedProject changes to update detailView
  viewModel.selectedProject.onChange { (_, _, newProjectOpt) =>
    Platform.runLater {
      newProjectOpt match {
        case Some(project) => renderProjectDetail(project)
        case None =>
          // Only show empty state if no subject is selected either
          if (viewModel.selectedSubject.value.isEmpty) {
            renderEmptyDetail("Select an item to view details")
          }
      }
    }
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
      viewModel.selectedProject.value = Some(newProject)
      // Clear subject selection when a project is selected
      if (viewModel.selectedSubject.value.isDefined) {
        viewModel.selectedSubject.value = None
      }
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
      }
    }
  }
  // UI to ViewModel sync
  sampleList.selectionModel().selectedItem.onChange { (_, _, newBiosample) =>
    if (newBiosample != null) {
      viewModel.selectedSubject.value = Some(newBiosample)
      // Clear project selection when a subject is selected
      if (viewModel.selectedProject.value.isDefined) {
        viewModel.selectedProject.value = None
      }
    } else if (viewModel.selectedSubject.value.isDefined && sampleList.selectionModel().getSelectedItem == null) {
      // Clear ViewModel selection if UI selection is cleared manually
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
      // Placeholder for new project dialog
      new Alert(AlertType.Information) {
        title = "New Project"
        headerText = "New Project functionality not yet implemented."
        contentText = "Stay tuned for project creation!"
      }.showAndWait()
    }
  }

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
      new Label("Projects:") { style = "-fx-font-weight: bold;" },
      projectList,
      newProjectButton,
      new Label("Subjects:") { style = "-fx-font-weight: bold;" },
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
}