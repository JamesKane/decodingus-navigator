package com.decodingus.ui.components

import com.decodingus.config.LabsConfig
import com.decodingus.workspace.model.{Biosample, HaplogroupAssignments, RecordMeta}
import scalafx.Includes.*
import scalafx.collections.ObservableBuffer
import scalafx.geometry.Insets
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, Priority}

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
  // Center name combo box (editable for custom entry)
  private val centerOptions = ObservableBuffer.from(
    "" +: LabsConfig.allLabNames :+ "Other..."
  )
  val centerNameCombo = new ComboBox[String]() {
    items = centerOptions
    editable = true
    promptText = "Sequencing Center (Optional)"
    prefWidth = 200
  }

  // Handle "Other..." selection to clear for custom entry
  centerNameCombo.selectionModel().selectedItem.onChange { (_, _, newVal) =>
    if (newVal == "Other...") {
      centerNameCombo.editor.value.clear()
      centerNameCombo.editor.value.requestFocus()
    }
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
    add(centerNameCombo, 1, 3)
    add(new Label("Description:"), 0, 4)
    add(descriptionField, 1, 4)
  }

  dialogPane().content = grid

  // Request focus on the donor ID field by default
  javafx.application.Platform.runLater(() => donorIdField.requestFocus())

  // Convert the result to a Biosample when the save button is clicked
  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType) {
      // Get value from editor (for editable combo) or selection
      val centerRaw = Option(centerNameCombo.editor.value.getText).getOrElse(centerNameCombo.value.value)
      val centerName = Option(centerRaw).map(_.trim).filter(v => v.nonEmpty && v != "Other...")
      val newBiosample = Biosample(
        atUri = None,
        meta = RecordMeta.initial,
        sampleAccession = accessionField.text.value,
        donorIdentifier = donorIdField.text.value,
        description = Option(descriptionField.text.value).filter(_.nonEmpty),
        centerName = centerName,
        sex = Option(sexChoiceBox.value.value),
        haplogroups = None,
        sequenceRunRefs = List.empty // No sequence runs initially
      )
      Some(newBiosample)
    } else {
      None
    }
  }
}
