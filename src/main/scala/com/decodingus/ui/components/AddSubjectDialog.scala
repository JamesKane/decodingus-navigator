package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ChoiceBox, ButtonBar}
import scalafx.scene.layout.{GridPane, Priority}
import scalafx.geometry.Insets
import scalafx.collections.ObservableBuffer
import com.decodingus.workspace.model.{Biosample, HaplogroupAssignments, SequenceData}
import java.time.LocalDateTime
import java.util.UUID

class AddSubjectDialog extends Dialog[Option[Biosample]] {
  title = "Add New Subject"
  headerText = "Enter details for the new Donor/Subject."

  val saveButtonType = new ButtonType("Save", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(saveButtonType, ButtonType.Cancel)

  val donorIdField = new TextField() {
    promptText = "Donor Identifier (e.g., John Doe)"
  }
  val accessionField = new TextField() {
    text = UUID.randomUUID().toString
    editable = false
    disable = true 
  }
  val sexChoiceBox = new ChoiceBox[String](ObservableBuffer("Male", "Female", "Other", "Unknown")) {
    value = "Unknown"
  }
  val centerNameField = new TextField() {
    promptText = "Sequencing Center (Optional)"
  }
  val descriptionField = new TextField() {
    promptText = "Description (Optional)"
  }

  val grid = new GridPane() {
    hgap = 10
    vgap = 10
    padding = Insets(20, 150, 10, 10)

    add(new Label("Donor ID:"), 0, 0)
    add(donorIdField, 1, 0)
    add(new Label("Accession (Auto):"), 0, 1)
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

  // Convert the result to a Biosample when the save button is clicked
  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType) {
      val newBiosample = Biosample(
        sampleAccession = accessionField.text.value,
        donorIdentifier = donorIdField.text.value,
        atUri = None,
        description = Option(descriptionField.text.value).filter(_.nonEmpty),
        centerName = Option(centerNameField.text.value).filter(_.nonEmpty),
        sex = Option(sexChoiceBox.value.value),
        sequenceData = List.empty, // Empty sequence data initially
        haplogroups = None,
        createdAt = Some(LocalDateTime.now())
      )
      Some(newBiosample)
    } else {
      None
    }
  }
}
