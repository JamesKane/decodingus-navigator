package com.decodingus.ui.components

import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.{Biosample, Project}
import scalafx.Includes.*
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.{Button, Label, ListCell, ListView}
import scalafx.scene.input.{ClipboardContent, DataFormat, DragEvent, TransferMode}
import scalafx.scene.layout.{HBox, Priority, VBox}

/** Companion object for shared constants */
object ProjectDetailView {
  // DataFormat must be a singleton - JavaFX throws if created multiple times
  val biosampleFormat: DataFormat = DataFormat("application/x-biosample-accession")
}

/**
 * Detail view for a selected Project showing:
 * - Project metadata with Edit/Delete actions
 * - Member list (drag-drop target for adding subjects from left panel)
 */
class ProjectDetailView(
                         viewModel: WorkbenchViewModel,
                         project: Project,
                         onEdit: Project => Unit,
                         onDelete: Project => Unit
                       ) extends VBox(10) {

  padding = Insets(10)

  // Import the shared DataFormat from companion object

  import ProjectDetailView.biosampleFormat

  // Observable buffer for project members
  private val memberBuffer: ObservableBuffer[Biosample] = ObservableBuffer.from(
    viewModel.getProjectMembers(project.projectName)
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
      }
    }

    // Drag over - accept drops from the main subjects list
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
        // Accept if subject exists and is not already a member
        val isAlreadyMember = memberBuffer.exists(_.sampleAccession == accession)
        val subjectExists = viewModel.findSubject(accession).isDefined
        if (subjectExists && !isAlreadyMember) {
          viewModel.addSubjectToProject(project.projectName, accession)
          refreshMembers()
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
        refreshMembers()
      }
    }
  }

  membersListView.selectionModel().selectedItem.onChange { (_, _, selected) =>
    removeFromProjectButton.disable = selected == null
  }

  // Refresh members list after add/remove operations
  private def refreshMembers(): Unit = {
    viewModel.findProject(project.projectName) match {
      case Some(updatedProject) =>
        memberBuffer.clear()
        memberBuffer ++= viewModel.getProjectMembers(updatedProject.projectName)
      case None =>
        memberBuffer.clear()
    }
  }

  private val dragHint = new Label("Drag subjects from the list on the left to add them") {
    style = "-fx-font-size: 11px; -fx-text-fill: #888888; -fx-padding: 5 0 0 0;"
  }

  // Layout
  children = Seq(
    new Label(s"Project: ${project.projectName}") {
      style = "-fx-font-size: 20px; -fx-font-weight: bold;"
    },
    actionButtons,
    infoSection,
    membersLabel,
    membersListView,
    new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(removeFromProjectButton)
    },
    dragHint
  )

  VBox.setVgrow(membersListView, Priority.Always)
}
