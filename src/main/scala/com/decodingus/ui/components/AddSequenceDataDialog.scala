package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ButtonBar, Button, ProgressIndicator, RadioButton, ToggleGroup}
import scalafx.scene.layout.{VBox, HBox, Priority, StackPane}
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.geometry.{Insets, Pos}
import scalafx.stage.FileChooser
import scalafx.application.Platform
import com.decodingus.workspace.model.FileInfo
import com.decodingus.analysis.AnalysisCache

import java.io.File
import java.net.URL
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.{Success, Failure}

/**
 * Result from the dialog - contains the file info and optional SHA256
 */
case class SequenceDataInput(
  fileInfo: FileInfo,
  sha256: Option[String]
)

/**
 * Simplified dialog for adding sequencing data.
 * Just handles file selection - metadata is inferred during analysis.
 * Supports:
 * - Local file selection via file browser
 * - Drag and drop of local files
 * - Cloud storage URLs (HTTP/HTTPS/S3)
 */
class AddSequenceDataDialog(
  existingChecksums: Set[String] // Checksums of already-added files to detect duplicates
) extends Dialog[Option[SequenceDataInput]] {

  title = "Add Sequencing Data"
  headerText = "Select an alignment file (BAM/CRAM)"

  val addButtonType = new ButtonType("Add", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(addButtonType, ButtonType.Cancel)

  // Source type toggle
  private val sourceToggle = new ToggleGroup()
  private val localFileRadio = new RadioButton("Local File") {
    toggleGroup = sourceToggle
    selected = true
  }
  private val cloudUrlRadio = new RadioButton("Cloud URL (HTTP/S3)") {
    toggleGroup = sourceToggle
  }

  private val sourceSelector = new HBox(20) {
    alignment = Pos.CenterLeft
    children = Seq(localFileRadio, cloudUrlRadio)
  }

  // Local file selection with drag-drop zone
  private var selectedFile: Option[File] = None
  private var computedSha256: Option[String] = None
  private var isComputingHash = false

  private val dropZoneLabel = new Label("Drag & drop BAM/CRAM file here\nor click Browse") {
    style = "-fx-font-size: 14px; -fx-text-fill: #888888;"
    wrapText = true
  }

  private val fileNameLabel = new Label() {
    style = "-fx-font-size: 12px; -fx-font-weight: bold;"
    visible = false
    managed = false
  }

  private val checksumLabel = new Label() {
    style = "-fx-font-size: 11px; -fx-text-fill: #666666;"
    visible = false
    managed = false
  }

  private val hashProgress = new ProgressIndicator() {
    prefWidth = 24
    prefHeight = 24
    visible = false
    managed = false
  }

  private val dropZone = new VBox(10) {
    alignment = Pos.Center
    padding = Insets(30)
    prefHeight = 150
    style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
    children = Seq(dropZoneLabel, fileNameLabel, new HBox(10) {
      alignment = Pos.Center
      children = Seq(hashProgress, checksumLabel)
    })

    // Drag over - accept files
    onDragOver = (event: DragEvent) => {
      if (event.dragboard.hasFiles && localFileRadio.selected.value) {
        val files = event.dragboard.getFiles
        if (files.size == 1 && isValidAlignmentFile(files.get(0))) {
          event.acceptTransferModes(TransferMode.Copy)
          style = "-fx-border-color: #4CAF50; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e8f5e9; -fx-background-radius: 8;"
        }
      }
      event.consume()
    }

    onDragExited = (_: DragEvent) => {
      style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
    }

    onDragDropped = (event: DragEvent) => {
      if (event.dragboard.hasFiles) {
        val files = event.dragboard.getFiles
        if (files.size == 1 && isValidAlignmentFile(files.get(0))) {
          handleFileSelected(files.get(0))
          event.dropCompleted = true
        }
      }
      event.consume()
    }
  }

  private val browseButton = new Button("Browse...") {
    onAction = _ => {
      val fileChooser = new FileChooser() {
        title = "Select Alignment File"
        extensionFilters.addAll(
          new FileChooser.ExtensionFilter("Alignment Files", Seq("*.bam", "*.cram")),
          new FileChooser.ExtensionFilter("BAM Files", "*.bam"),
          new FileChooser.ExtensionFilter("CRAM Files", "*.cram")
        )
      }
      Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach(handleFileSelected)
    }
  }

  private val localFilePane = new VBox(10) {
    children = Seq(dropZone, browseButton)
    alignment = Pos.Center
  }

  // Cloud URL input
  private val urlField = new TextField() {
    promptText = "https://... or s3://bucket/path/file.bam"
    prefWidth = 400
  }

  private val urlHelpLabel = new Label("Enter HTTP(S) or S3 URL to a BAM/CRAM file") {
    style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
  }

  private val cloudUrlPane = new VBox(10) {
    children = Seq(urlField, urlHelpLabel)
    visible = false
    managed = false
  }

  // Toggle between local and cloud
  sourceToggle.selectedToggle.onChange { (_, _, _) =>
    val isLocal = localFileRadio.selected.value
    localFilePane.visible = isLocal
    localFilePane.managed = isLocal
    cloudUrlPane.visible = !isLocal
    cloudUrlPane.managed = !isLocal
    updateAddButton()
  }

  private val contentPane = new VBox(15) {
    padding = Insets(20)
    children = Seq(sourceSelector, localFilePane, cloudUrlPane)
  }

  dialogPane().content = contentPane

  // Add button state
  private val addButton = dialogPane().lookupButton(addButtonType)
  addButton.disable = true

  urlField.text.onChange { (_, _, _) => updateAddButton() }

  private def updateAddButton(): Unit = {
    val isValid = if (localFileRadio.selected.value) {
      selectedFile.isDefined && !isComputingHash
    } else {
      val url = urlField.text.value
      url != null && url.nonEmpty && isValidUrl(url)
    }
    addButton.disable = !isValid
  }

  private def isValidAlignmentFile(file: File): Boolean = {
    val name = file.getName.toLowerCase
    name.endsWith(".bam") || name.endsWith(".cram")
  }

  private def isValidUrl(url: String): Boolean = {
    url.startsWith("http://") || url.startsWith("https://") || url.startsWith("s3://")
  }

  private def handleFileSelected(file: File): Unit = {
    selectedFile = Some(file)
    computedSha256 = None
    isComputingHash = true

    dropZoneLabel.visible = false
    dropZoneLabel.managed = false
    fileNameLabel.text = file.getName
    fileNameLabel.visible = true
    fileNameLabel.managed = true
    checksumLabel.text = "Computing checksum..."
    checksumLabel.visible = true
    checksumLabel.managed = true
    hashProgress.visible = true
    hashProgress.managed = true

    dropZone.style = "-fx-border-color: #4CAF50; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e8f5e9; -fx-background-radius: 8;"

    updateAddButton()

    // Calculate SHA256 in background
    Future {
      AnalysisCache.calculateSha256(file)
    }.onComplete {
      case Success(sha256) =>
        Platform.runLater {
          computedSha256 = Some(sha256)
          isComputingHash = false
          hashProgress.visible = false
          hashProgress.managed = false

          if (existingChecksums.contains(sha256)) {
            checksumLabel.text = s"⚠ Duplicate: ${sha256.take(12)}..."
            checksumLabel.style = "-fx-font-size: 11px; -fx-text-fill: #F44336;"
            dropZone.style = "-fx-border-color: #FF9800; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #fff3e0; -fx-background-radius: 8;"
          } else {
            checksumLabel.text = s"SHA256: ${sha256.take(12)}..."
            checksumLabel.style = "-fx-font-size: 11px; -fx-text-fill: #4CAF50;"
          }
          updateAddButton()
        }
      case Failure(e) =>
        Platform.runLater {
          isComputingHash = false
          hashProgress.visible = false
          hashProgress.managed = false
          checksumLabel.text = s"Checksum failed: ${e.getMessage}"
          checksumLabel.style = "-fx-font-size: 11px; -fx-text-fill: #F44336;"
          updateAddButton()
        }
    }
  }

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == addButtonType) {
      if (localFileRadio.selected.value) {
        selectedFile.map { file =>
          val fileFormat = if (file.getName.toLowerCase.endsWith(".cram")) "CRAM" else "BAM"
          val fileInfo = FileInfo(
            fileName = file.getName,
            fileSizeBytes = Some(file.length()),
            fileFormat = fileFormat,
            checksum = computedSha256,
            checksumAlgorithm = computedSha256.map(_ => "SHA-256"),
            location = Some(file.getAbsolutePath)
          )
          SequenceDataInput(fileInfo, computedSha256)
        }
      } else {
        val url = urlField.text.value
        if (url != null && url.nonEmpty) {
          val fileName = url.split("/").lastOption.getOrElse("remote_file")
          val fileFormat = if (fileName.toLowerCase.endsWith(".cram")) "CRAM" else "BAM"
          val fileInfo = FileInfo(
            fileName = fileName,
            fileSizeBytes = None, // Unknown for remote files
            fileFormat = fileFormat,
            checksum = None, // Will be computed during analysis
            checksumAlgorithm = None,
            location = Some(url)
          )
          Some(SequenceDataInput(fileInfo, None))
        } else {
          None
        }
      }
    } else {
      None
    }
  }
}
