package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.layout.{VBox, HBox, Priority, Region}
import scalafx.scene.control.{Label, Button, ListView, ListCell, Alert, ButtonType}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.input.{DragEvent, TransferMode, ClipboardContent, DataFormat}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import scalafx.application.Platform
import com.decodingus.workspace.model.{Project, Biosample}
import com.decodingus.workspace.WorkbenchViewModel

/**
 * Detail view for a selected Project showing:
 * - Project metadata with Edit/Delete actions
 * - Member list (drag-drop target for adding subjects)
 * - Available subjects list (drag source for adding to project)
 */
class ProjectDetailView(
  viewModel: WorkbenchViewModel,
  project: Project,
  onEdit: Project => Unit,
  onDelete: Project => Unit
) extends VBox(10) {

  padding = Insets(10)

  // Data format for drag-drop
  private val biosampleFormat = DataFormat("application/x-biosample-accession")

  // Observable buffers for the lists
  private val memberBuffer: ObservableBuffer[Biosample] = ObservableBuffer.from(
    viewModel.getProjectMembers(project.projectName)
  )
  private val availableBuffer: ObservableBuffer[Biosample] = ObservableBuffer.from(
    viewModel.getNonProjectMembers(project.projectName)
  )

  // Header with project name and action buttons
  private val editButton = new Button("Edit") {
    onAction = _ => onEdit(project)
  }

  private val deleteButton = new Button("Delete") {
    style = "-fx-text-fill: #D32F2F;"
    onAction = _ => onDelete(project)
  }

  private val actionButtons = new HBox(10) {
    padding = Insets(10, 0, 10, 0)
    children = Seq(editButton, deleteButton)
  }

  // Project info section
  private val infoSection = new VBox(5) {
    children = Seq(
      new Label(s"Description: ${project.description.getOrElse("No description")}"),
      new Label(s"Administrator: ${project.administrator}"),
      new Label(s"Created: ${project.projectName}") // Could add createdAt field
    )
  }

  // Members list with drag-drop support
  private val membersLabel = new Label("Project Members:") {
    style = "-fx-font-weight: bold;"
  }

  private val membersListView = new ListView[Biosample](memberBuffer) {
    prefHeight = 150
    placeholder = new Label("Drag subjects here to add them") {
      style = "-fx-text-fill: #888888;"
    }

    cellFactory = { (v: ListView[Biosample]) =>
      new ListCell[Biosample] {
        item.onChange { (_, _, biosample) =>
          if (biosample != null) {
            text = s"${biosample.donorIdentifier} (${biosample.sampleAccession.take(8)}...)"
            graphic = null
          } else {
            text = null
            graphic = null
          }
        }

        // Drag source - drag from members list to remove
        onDragDetected = (event) => {
          Option(item.value).foreach { biosample =>
            val db = startDragAndDrop(TransferMode.Move)
            val content = new ClipboardContent()
            content.put(biosampleFormat, biosample.sampleAccession)
            content.putString(biosample.sampleAccession)
            db.setContent(content)
            event.consume()
          }
        }
      }
    }

    // Drag over - accept drops from available list
    onDragOver = (event: DragEvent) => {
      if (event.gestureSource != this && event.dragboard.hasString) {
        event.acceptTransferModes(TransferMode.Move)
      }
      event.consume()
    }

    // Drag dropped - add subject to project
    onDragDropped = (event: DragEvent) => {
      val success = if (event.dragboard.hasString) {
        val accession = event.dragboard.getString
        // Only accept if it's from the available list (not already a member)
        if (availableBuffer.exists(_.sampleAccession == accession)) {
          viewModel.addSubjectToProject(project.projectName, accession)
          refreshLists()
          true
        } else false
      } else false
      event.dropCompleted = success
      event.consume()
    }

    onDragEntered = (_: DragEvent) => {
      style = "-fx-background-color: #e8f5e9;"
    }

    onDragExited = (_: DragEvent) => {
      style = ""
    }
  }

  private val removeFromProjectButton = new Button("Remove Selected") {
    disable = true
    onAction = _ => {
      Option(membersListView.selectionModel().getSelectedItem).foreach { biosample =>
        viewModel.removeSubjectFromProject(project.projectName, biosample.sampleAccession)
        refreshLists()
      }
    }
  }

  membersListView.selectionModel().selectedItem.onChange { (_, _, selected) =>
    removeFromProjectButton.disable = selected == null
  }

  // Available subjects list (not in project)
  private val availableLabel = new Label("Available Subjects:") {
    style = "-fx-font-weight: bold; -fx-padding: 15 0 0 0;"
  }

  private val availableHint = new Label("Drag subjects to the members list above to add them") {
    style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
  }

  private val availableListView = new ListView[Biosample](availableBuffer) {
    prefHeight = 150
    placeholder = new Label("All subjects are in this project") {
      style = "-fx-text-fill: #888888;"
    }

    cellFactory = { (v: ListView[Biosample]) =>
      new ListCell[Biosample] {
        item.onChange { (_, _, biosample) =>
          if (biosample != null) {
            text = s"${biosample.donorIdentifier} (${biosample.sampleAccession.take(8)}...)"
            graphic = null
          } else {
            text = null
            graphic = null
          }
        }

        // Drag source - drag from available to add
        onDragDetected = (event) => {
          Option(item.value).foreach { biosample =>
            val db = startDragAndDrop(TransferMode.Move)
            val content = new ClipboardContent()
            content.put(biosampleFormat, biosample.sampleAccession)
            content.putString(biosample.sampleAccession)
            db.setContent(content)
            event.consume()
          }
        }
      }
    }

    // Drag over - accept drops from members list (for removal)
    onDragOver = (event: DragEvent) => {
      if (event.gestureSource != this && event.dragboard.hasString) {
        event.acceptTransferModes(TransferMode.Move)
      }
      event.consume()
    }

    // Drag dropped - remove subject from project
    onDragDropped = (event: DragEvent) => {
      val success = if (event.dragboard.hasString) {
        val accession = event.dragboard.getString
        // Only accept if it's from the members list
        if (memberBuffer.exists(_.sampleAccession == accession)) {
          viewModel.removeSubjectFromProject(project.projectName, accession)
          refreshLists()
          true
        } else false
      } else false
      event.dropCompleted = success
      event.consume()
    }

    onDragEntered = (_: DragEvent) => {
      style = "-fx-background-color: #fff3e0;"
    }

    onDragExited = (_: DragEvent) => {
      style = ""
    }
  }

  private val addToProjectButton = new Button("Add Selected") {
    disable = true
    onAction = _ => {
      Option(availableListView.selectionModel().getSelectedItem).foreach { biosample =>
        viewModel.addSubjectToProject(project.projectName, biosample.sampleAccession)
        refreshLists()
      }
    }
  }

  availableListView.selectionModel().selectedItem.onChange { (_, _, selected) =>
    addToProjectButton.disable = selected == null
  }

  // Refresh lists after add/remove operations
  private def refreshLists(): Unit = {
    // Re-fetch from ViewModel to get updated state
    viewModel.findProject(project.projectName) match {
      case Some(updatedProject) =>
        memberBuffer.clear()
        memberBuffer ++= viewModel.getProjectMembers(updatedProject.projectName)
        availableBuffer.clear()
        availableBuffer ++= viewModel.getNonProjectMembers(updatedProject.projectName)
      case None =>
        // Project was deleted
        memberBuffer.clear()
        availableBuffer.clear()
    }
  }

  // Layout
  children = Seq(
    new Label(s"Project: ${project.projectName}") { style = "-fx-font-size: 20px; -fx-font-weight: bold;" },
    actionButtons,
    infoSection,
    membersLabel,
    membersListView,
    new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(removeFromProjectButton)
    },
    availableLabel,
    availableHint,
    availableListView,
    new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(addToProjectButton)
    }
  )

  VBox.setVgrow(membersListView, Priority.Sometimes)
  VBox.setVgrow(availableListView, Priority.Sometimes)
}
