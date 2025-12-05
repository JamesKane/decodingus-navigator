package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ChoiceBox, ButtonBar, Button}
import scalafx.scene.layout.{GridPane, HBox, Priority}
import scalafx.geometry.Insets
import scalafx.collections.ObservableBuffer
import scalafx.stage.FileChooser
import com.decodingus.workspace.model.{SequenceData, FileInfo, AlignmentData}

import java.io.File

/**
 * Dialog for adding sequencing run metadata to a subject.
 * Allows user to specify platform, test type, and select BAM/CRAM file.
 */
class AddSequenceDataDialog extends Dialog[Option[SequenceData]] {
  title = "Add Sequencing Data"
  headerText = "Enter details for the sequencing run."

  val saveButtonType = new ButtonType("Add", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(saveButtonType, ButtonType.Cancel)

  // Platform selection
  val platformChoiceBox = new ChoiceBox[String](ObservableBuffer(
    "Illumina", "PacBio", "Oxford Nanopore", "Ion Torrent", "BGI/MGI", "Other"
  )) {
    value = "Illumina"
  }

  // Instrument model (optional, text field)
  val instrumentField = new TextField() {
    promptText = "e.g., NovaSeq 6000, Sequel II"
  }

  // Test type
  val testTypeChoiceBox = new ChoiceBox[String](ObservableBuffer(
    "WGS", "WES", "Targeted Panel", "RNA-Seq", "HiFi", "Other"
  )) {
    value = "WGS"
  }

  // Library layout
  val layoutChoiceBox = new ChoiceBox[String](ObservableBuffer(
    "Paired-End", "Single-End", "Unknown"
  )) {
    value = "Paired-End"
  }

  // File selection
  val filePathField = new TextField() {
    editable = false
    promptText = "Select BAM/CRAM file..."
    prefWidth = 300
  }

  private var selectedFile: Option[File] = None

  val browseButton = new Button("Browse...") {
    onAction = _ => {
      val fileChooser = new FileChooser() {
        title = "Select Alignment File"
        extensionFilters.addAll(
          new FileChooser.ExtensionFilter("Alignment Files", Seq("*.bam", "*.cram")),
          new FileChooser.ExtensionFilter("BAM Files", "*.bam"),
          new FileChooser.ExtensionFilter("CRAM Files", "*.cram"),
          new FileChooser.ExtensionFilter("All Files", "*.*")
        )
      }
      Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach { file =>
        selectedFile = Some(file)
        filePathField.text = file.getAbsolutePath
      }
    }
  }

  val fileSelectionBox = new HBox(5) {
    children = Seq(filePathField, browseButton)
  }

  val grid = new GridPane() {
    hgap = 10
    vgap = 10
    padding = Insets(20, 20, 10, 10)

    add(new Label("Platform:"), 0, 0)
    add(platformChoiceBox, 1, 0)
    add(new Label("Instrument:"), 0, 1)
    add(instrumentField, 1, 1)
    add(new Label("Test Type:"), 0, 2)
    add(testTypeChoiceBox, 1, 2)
    add(new Label("Library Layout:"), 0, 3)
    add(layoutChoiceBox, 1, 3)
    add(new Label("Alignment File:"), 0, 4)
    add(fileSelectionBox, 1, 4)
  }

  dialogPane().content = grid

  // Disable save button until file is selected
  val saveButton = dialogPane().lookupButton(saveButtonType)
  saveButton.disable = true

  filePathField.text.onChange { (_, _, newValue) =>
    saveButton.disable = newValue == null || newValue.isEmpty
  }

  // Convert the result to SequenceData when the save button is clicked
  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType && selectedFile.isDefined) {
      val file = selectedFile.get
      val fileFormat = if (file.getName.toLowerCase.endsWith(".cram")) "CRAM" else "BAM"

      val fileInfo = FileInfo(
        fileName = file.getName,
        fileSizeBytes = Some(file.length()),
        fileFormat = fileFormat,
        checksum = None, // Will be calculated during analysis
        location = file.getAbsolutePath
      )

      val sequenceData = SequenceData(
        platformName = platformChoiceBox.value.value,
        instrumentModel = Option(instrumentField.text.value).filter(_.nonEmpty),
        testType = testTypeChoiceBox.value.value,
        libraryLayout = Some(layoutChoiceBox.value.value),
        totalReads = None, // Will be populated during analysis
        readLength = None,
        meanInsertSize = None,
        files = List(fileInfo),
        alignments = List.empty // Will be populated during analysis
      )

      Some(sequenceData)
    } else {
      None
    }
  }
}
