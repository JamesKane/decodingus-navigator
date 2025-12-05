package com.decodingus.ui.components

import scalafx.scene.layout.{VBox, HBox, Priority, StackPane, BorderPane}
import scalafx.scene.control.{Label, Button, ListView, SplitPane, Alert, ListCell}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import com.decodingus.workspace.model.{Workspace, Project, Biosample}
import com.decodingus.workspace.{WorkspaceService, LiveWorkspaceService}
import scalafx.scene.control.Alert.AlertType
import scalafx.application.Platform

class WorkbenchView(
  var workspace: Workspace, // Mutable var for now, will be updated by UI actions
  workspaceService: WorkspaceService // Now accepts the trait
) extends SplitPane { // Change to SplitPane

  // Observable buffers for UI lists
  private val projectBuffer: ObservableBuffer[Project] = ObservableBuffer(workspace.projects*)
  private val sampleBuffer: ObservableBuffer[Biosample] = ObservableBuffer(workspace.samples*)

  // Left Panel - Navigation
  private val projectList = new ListView[Project](projectBuffer) {
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

  private val sampleList = new ListView[Biosample](sampleBuffer) {
    vgrow = Priority.Always
    cellFactory = { (v: ListView[Biosample]) =>
      new ListCell[Biosample] {
        item.onChange { (_, _, newBiosample) =>
          text = if (newBiosample != null) s"${newBiosample.donorIdentifier} (${newBiosample.sampleAccession.take(8)}...)" else null
        }
      }
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

  private val addSampleButton = new Button("Add Sample") {
    onAction = _ => {
      // Placeholder for adding sample (e.g., via file dialog or manual entry)
      new Alert(AlertType.Information) {
        title = "Add Sample"
        headerText = "Add Sample functionality not yet implemented."
        contentText = "Stay tuned for sample addition!"
      }.showAndWait()
    }
  }

  private val saveButton = new Button("Save Workspace") {
    styleClass.add("button-primary")
    onAction = _ => saveWorkspace()
  }

  private val leftPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new Label("Projects:") { style = "-fx-font-weight: bold;" },
      projectList,
      newProjectButton,
      new Label("Samples:") { style = "-fx-font-weight: bold;" },
      sampleList,
      addSampleButton,
      new HBox(10) {
        alignment = Pos.BottomRight
        children = Seq(saveButton)
        HBox.setHgrow(saveButton, Priority.Always)
      }
    )
  }
  SplitPane.setResizableWithParent(leftPanel, false) // Make left panel not resize with parent by default

  // Right Panel - Details/Content Area
  private val rightPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new Label("Details Panel") { style = "-fx-font-size: 18px; -fx-font-weight: bold;" },
      new Label("Select a project or sample to view its details here.")
    )
  }
  VBox.setVgrow(rightPanel, Priority.Always) // Allow right panel to grow vertically

  // Set the items of the SplitPane
  items.addAll(leftPanel, rightPanel)
  dividerPositions = 0.25 // Initial divider position

  // Method to update the status (used by GenomeNavigatorApp) - now updates lists
  def updateWorkspace(newWorkspace: Workspace): Unit = {
    workspace = newWorkspace
    Platform.runLater {
      projectBuffer.clear()
      projectBuffer ++= workspace.projects
      sampleBuffer.clear()
      sampleBuffer ++= workspace.samples
    }
  }

  // Method to save the current workspace
  private def saveWorkspace(): Unit = {
    workspaceService.save(workspace).fold(
      error => {
        new Alert(AlertType.Error) {
          title = "Save Error"
          headerText = "Could not save workspace"
          contentText = s"Reason: $error"
        }.showAndWait()
      },
      _ => {
        new Alert(AlertType.Information) {
          title = "Workspace Saved"
          headerText = "Workspace saved successfully!"
          contentText = "Your projects and samples have been saved to workspace.json."
        }.showAndWait()
      }
    )
  }
}