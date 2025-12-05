package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, TextArea, ButtonBar}
import scalafx.scene.layout.{GridPane, VBox, Priority}
import scalafx.geometry.Insets
import com.decodingus.workspace.model.Project

/**
 * Dialog for editing an existing Project.
 * Note: Members are managed separately via drag-drop in the project detail view.
 */
class EditProjectDialog(existingProject: Project) extends Dialog[Option[Project]] {
  title = "Edit Project"
  headerText = s"Edit ${existingProject.projectName}"

  val saveButtonType = new ButtonType("Save", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(saveButtonType, ButtonType.Cancel)

  private val nameField = new TextField() {
    text = existingProject.projectName
    promptText = "Project name"
    prefWidth = 300
  }

  private val descriptionField = new TextArea() {
    text = existingProject.description.getOrElse("")
    promptText = "Project description (optional)"
    prefWidth = 300
    prefHeight = 80
    wrapText = true
  }

  private val adminField = new TextField() {
    text = existingProject.administrator
    promptText = "Administrator name"
  }

  private val grid = new GridPane() {
    hgap = 10
    vgap = 10
    padding = Insets(20)

    add(new Label("Name:"), 0, 0)
    add(nameField, 1, 0)
    add(new Label("Description:"), 0, 1)
    add(descriptionField, 1, 1)
    add(new Label("Administrator:"), 0, 2)
    add(adminField, 1, 2)
    add(new Label("Members:"), 0, 3)
    add(new Label(s"${existingProject.members.size} subject(s) - manage in project view") {
      style = "-fx-text-fill: #888888; -fx-font-style: italic;"
    }, 1, 3)
  }

  dialogPane().content = grid

  // Disable save until name is provided
  private val saveButton = dialogPane().lookupButton(saveButtonType)
  saveButton.disable = nameField.text.value.trim.isEmpty

  nameField.text.onChange { (_, _, newValue) =>
    saveButton.disable = newValue == null || newValue.trim.isEmpty
  }

  // Focus on name field
  javafx.application.Platform.runLater(() => nameField.requestFocus())

  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType) {
      Some(existingProject.copy(
        projectName = nameField.text.value.trim,
        description = Option(descriptionField.text.value).map(_.trim).filter(_.nonEmpty),
        administrator = adminField.text.value.trim
        // members are preserved - managed separately
      ))
    } else {
      None
    }
  }
}
