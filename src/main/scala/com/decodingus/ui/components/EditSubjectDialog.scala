package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ChoiceBox, ButtonBar}
import scalafx.scene.layout.{GridPane, Priority}
import scalafx.geometry.Insets
import scalafx.collections.ObservableBuffer
import com.decodingus.workspace.model.Biosample

/**
 * Dialog for editing an existing Subject/Biosample.
 * Pre-populates fields with the existing subject data.
 */
class EditSubjectDialog(existingSubject: Biosample) extends Dialog[Option[Biosample]] {
  title = "Edit Subject"
  headerText = s"Edit details for ${existingSubject.donorIdentifier}"

  val saveButtonType = new ButtonType("Save", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(saveButtonType, ButtonType.Cancel)

  val donorIdField = new TextField() {
    text = existingSubject.donorIdentifier
    promptText = "Donor Identifier (e.g., John Doe)"
  }

  val accessionField = new TextField() {
    text = existingSubject.sampleAccession
    editable = false
    disable = true
  }

  val sexChoiceBox = new ChoiceBox[String](ObservableBuffer("Male", "Female", "Other", "Unknown")) {
    value = existingSubject.sex.getOrElse("Unknown")
  }

  val centerNameField = new TextField() {
    text = existingSubject.centerName.getOrElse("")
    promptText = "Sequencing Center (Optional)"
  }

  val descriptionField = new TextField() {
    text = existingSubject.description.getOrElse("")
    promptText = "Description (Optional)"
  }

  val grid = new GridPane() {
    hgap = 10
    vgap = 10
    padding = Insets(20, 150, 10, 10)

    add(new Label("Donor ID:"), 0, 0)
    add(donorIdField, 1, 0)
    add(new Label("Accession:"), 0, 1)
    add(accessionField, 1, 1)
    add(new Label("Sex:"), 0, 2)
    add(sexChoiceBox, 1, 2)
    add(new Label("Center Name:"), 0, 3)
    add(centerNameField, 1, 3)
    add(new Label("Description:"), 0, 4)
    add(descriptionField, 1, 4)
  }

  dialogPane().content = grid

  // Request focus on the donor ID field by default
  javafx.application.Platform.runLater(() => donorIdField.requestFocus())

  // Convert the result to an updated Biosample when the save button is clicked
  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType) {
      // Preserve existing data that shouldn't change (sequenceData, haplogroups, createdAt)
      val updatedBiosample = existingSubject.copy(
        donorIdentifier = donorIdField.text.value,
        description = Option(descriptionField.text.value).filter(_.nonEmpty),
        centerName = Option(centerNameField.text.value).filter(_.nonEmpty),
        sex = Option(sexChoiceBox.value.value)
      )
      Some(updatedBiosample)
    } else {
      None
    }
  }
}
