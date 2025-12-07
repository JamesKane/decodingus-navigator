package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, TextArea, ButtonBar}
import scalafx.scene.layout.{GridPane, VBox, Priority}
import scalafx.geometry.Insets
import com.decodingus.workspace.model.{Project, RecordMeta}

/**
 * Dialog for creating a new Project.
 */
class AddProjectDialog(defaultAdministrator: String = "Local User") extends Dialog[Option[Project]] {
  title = "Create New Project"
  headerText = "Enter project details"

  val createButtonType = new ButtonType("Create", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(createButtonType, ButtonType.Cancel)

  private val nameField = new TextField() {
    promptText = "Project name"
    prefWidth = 300
  }

  private val descriptionField = new TextArea() {
    promptText = "Project description (optional)"
    prefWidth = 300
    prefHeight = 80
    wrapText = true
  }

  private val adminField = new TextField() {
    text = defaultAdministrator
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
  }

  dialogPane().content = grid

  // Disable create until name is provided
  private val createButton = dialogPane().lookupButton(createButtonType)
  createButton.disable = true

  nameField.text.onChange { (_, _, newValue) =>
    createButton.disable = newValue == null || newValue.trim.isEmpty
  }

  // Focus on name field
  javafx.application.Platform.runLater(() => nameField.requestFocus())

  resultConverter = dialogButton => {
    if (dialogButton == createButtonType) {
      Some(Project(
        atUri = None,
        meta = RecordMeta.initial,
        projectName = nameField.text.value.trim,
        description = Option(descriptionField.text.value).map(_.trim).filter(_.nonEmpty),
        administrator = adminField.text.value.trim,
        memberRefs = List.empty
      ))
    } else {
      None
    }
  }
}
