package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ComboBox, ButtonBar}
import scalafx.scene.layout.{GridPane, VBox}
import scalafx.geometry.Insets
import scalafx.collections.ObservableBuffer
import com.decodingus.workspace.model.{SequenceRun, Alignment}

/**
 * Dialog for editing sequencing run metadata.
 * Allows users to correct or customize auto-detected values.
 */
class EditSequenceDataDialog(
  existingRun: SequenceRun,
  alignments: List[Alignment] = List.empty
) extends Dialog[Option[SequenceRun]] {
  title = "Edit Sequencing Run"
  headerText = "Edit sequencing run metadata"

  val saveButtonType = new ButtonType("Save", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(saveButtonType, ButtonType.Cancel)

  // Lab (sequencing facility) - editable for manual override
  private val labField = new TextField() {
    text = existingRun.sequencingFacility.getOrElse("")
    promptText = "e.g., Nebula Genomics, Dante Labs"
    prefWidth = 200
  }

  // Sample name (from BAM @RG SM tag) - editable for manual override
  private val sampleNameField = new TextField() {
    text = existingRun.sampleName.getOrElse("")
    promptText = "Sample name from BAM header"
    prefWidth = 200
  }

  // Platform selection
  private val platformCombo = new ComboBox[String] {
    items = ObservableBuffer(
      "Illumina",
      "PacBio",
      "Oxford Nanopore",
      "MGI",
      "Ion Torrent",
      "Complete Genomics",
      "Other"
    )
    value = existingRun.platformName
    editable = true
    prefWidth = 200
  }

  // Test type selection
  private val testTypeCombo = new ComboBox[String] {
    items = ObservableBuffer(
      "WGS",
      "WES",
      "HiFi",
      "CLR",
      "Nanopore",
      "Targeted Panel",
      "RNA-Seq",
      "Other"
    )
    value = existingRun.testType
    editable = true
    prefWidth = 200
  }

  // Instrument model
  private val instrumentField = new TextField() {
    text = existingRun.instrumentModel.getOrElse("")
    promptText = "e.g., NovaSeq 6000, Sequel II"
    prefWidth = 200
  }

  // Library layout selection
  private val layoutCombo = new ComboBox[String] {
    items = ObservableBuffer("Paired-End", "Single-End", "Unknown")
    value = existingRun.libraryLayout.getOrElse("Unknown")
    prefWidth = 200
  }

  // Read statistics (editable)
  private val totalReadsField = new TextField() {
    text = existingRun.totalReads.map(_.toString).getOrElse("")
    promptText = "Total number of reads"
    prefWidth = 200
  }

  private val readLengthField = new TextField() {
    text = existingRun.readLength.map(_.toString).getOrElse("")
    promptText = "Read length (bp)"
    prefWidth = 200
  }

  private val insertSizeField = new TextField() {
    text = existingRun.meanInsertSize.map(d => f"$d%.1f").getOrElse("")
    promptText = "Mean insert size (bp)"
    prefWidth = 200
  }

  // File info (read-only)
  private val fileNameLabel = new Label(
    existingRun.files.headOption.map(_.fileName).getOrElse("No file")
  ) {
    style = "-fx-font-style: italic;"
  }

  // Alignment info (read-only) - from the alignments list
  private val referenceLabel = new Label(
    alignments.headOption.map(_.referenceBuild).getOrElse("Not analyzed")
  ) {
    style = "-fx-font-style: italic;"
  }

  private val alignerLabel = new Label(
    alignments.headOption.map(_.aligner).getOrElse("Unknown")
  ) {
    style = "-fx-font-style: italic;"
  }

  private val grid = new GridPane() {
    hgap = 10
    vgap = 10
    padding = Insets(20)

    // Lab and Sample (at top - most likely to need manual editing)
    add(new Label("Lab:"), 0, 0)
    add(labField, 1, 0)

    add(new Label("Sample:"), 0, 1)
    add(sampleNameField, 1, 1)

    // Editable fields
    add(new Label("Platform:"), 0, 2)
    add(platformCombo, 1, 2)

    add(new Label("Test Type:"), 0, 3)
    add(testTypeCombo, 1, 3)

    add(new Label("Instrument:"), 0, 4)
    add(instrumentField, 1, 4)

    add(new Label("Library Layout:"), 0, 5)
    add(layoutCombo, 1, 5)

    // Statistics (editable)
    add(new Label("Total Reads:"), 0, 6)
    add(totalReadsField, 1, 6)

    add(new Label("Read Length:"), 0, 7)
    add(readLengthField, 1, 7)

    add(new Label("Insert Size:"), 0, 8)
    add(insertSizeField, 1, 8)

    // Read-only info
    add(new Label("") { prefHeight = 10 }, 0, 9) // Spacer

    add(new Label("File:") { style = "-fx-text-fill: #888888;" }, 0, 10)
    add(fileNameLabel, 1, 10)

    add(new Label("Reference:") { style = "-fx-text-fill: #888888;" }, 0, 11)
    add(referenceLabel, 1, 11)

    add(new Label("Aligner:") { style = "-fx-text-fill: #888888;" }, 0, 12)
    add(alignerLabel, 1, 12)
  }

  private val content = new VBox(10) {
    children = Seq(
      grid,
      new Label("Note: File, reference, and aligner are determined by analysis and cannot be edited here.") {
        style = "-fx-text-fill: #888888; -fx-font-size: 11px; -fx-font-style: italic;"
        wrapText = true
        maxWidth = 350
      }
    )
  }

  dialogPane().content = content

  // Focus on lab field (most likely to need editing)
  javafx.application.Platform.runLater(() => labField.requestFocus())

  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType) {
      // Parse numeric fields
      val totalReads = Option(totalReadsField.text.value)
        .filter(_.nonEmpty)
        .flatMap(s => scala.util.Try(s.toLong).toOption)

      val readLength = Option(readLengthField.text.value)
        .filter(_.nonEmpty)
        .flatMap(s => scala.util.Try(s.toInt).toOption)

      val insertSize = Option(insertSizeField.text.value)
        .filter(_.nonEmpty)
        .flatMap(s => scala.util.Try(s.toDouble).toOption)

      val layoutValue = layoutCombo.value.value
      val libraryLayout = if (layoutValue == "Unknown") None else Some(layoutValue)

      Some(existingRun.copy(
        sequencingFacility = Option(labField.text.value).filter(_.nonEmpty),
        sampleName = Option(sampleNameField.text.value).filter(_.nonEmpty),
        platformName = platformCombo.value.value,
        testType = testTypeCombo.value.value,
        instrumentModel = Option(instrumentField.text.value).filter(_.nonEmpty),
        libraryLayout = libraryLayout,
        totalReads = totalReads,
        readLength = readLength,
        meanInsertSize = insertSize
        // files and alignmentRefs are preserved from existingRun
      ))
    } else {
      None
    }
  }
}
