package com.decodingus.ui.components

import com.decodingus.analysis.AnalysisCache
import com.decodingus.i18n.I18n.t
import com.decodingus.util.{DetectedFileType, FileTypeDetector}
import com.decodingus.workspace.model.FileInfo
import scalafx.Includes.*
import scalafx.application.Platform
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.{HBox, Priority, VBox}
import scalafx.stage.FileChooser

import java.io.File
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.util.{Failure, Success}

/**
 * Supported data types for import
 */
sealed trait DataType {
  def label: String
}

object DataType {
  case object Alignment extends DataType {
    val label = "Alignment (BAM/CRAM)"
  }

  case object Variants extends DataType {
    val label = "Variants (VCF)"
  }

  case object StrProfile extends DataType {
    val label = "STR Profile"
  }

  case object MtdnaFasta extends DataType {
    val label = "mtDNA FASTA"
  }

  case class ChipData(vendor: Option[String] = None) extends DataType {
    val label = vendor.map(v => s"$v Chip Data").getOrElse("Chip/SNP Data")
  }

  case object Unknown extends DataType {
    val label = "Unknown"
  }
}

/**
 * Result from the Add Data dialog
 */
case class DataInput(
  dataType: DataType,
  fileInfo: FileInfo,
  sha256: Option[String]
)

/**
 * Unified dialog for adding data files with auto-detection.
 *
 * Automatically detects file type based on:
 * - File extension (BAM, CRAM, VCF)
 * - Content fingerprinting for CSV/TXT (STR profiles vs chip data)
 */
class AddDataDialog(
  existingChecksums: Set[String] = Set.empty
) extends Dialog[Option[DataInput]] {

  title = t("data.add")
  headerText = t("data.add_file")

  val addButtonType = new ButtonType(t("action.add"), ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(addButtonType, ButtonType.Cancel)

  // State
  private var selectedFile: Option[File] = None
  private var detectedType: DataType = DataType.Unknown
  private var computedSha256: Option[String] = None
  private var isProcessing = false

  // UI Elements - Drop zone
  private val dropZoneLabel = new Label(t("data.drag_drop_or_browse")) {
    style = "-fx-font-size: 14px; -fx-text-fill: #888888;"
    wrapText = true
  }

  private val fileNameLabel = new Label() {
    style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #333333;"
    visible = false
    managed = false
  }

  private val detectedTypeLabel = new Label() {
    style = "-fx-font-size: 12px; -fx-text-fill: #4CAF50; -fx-font-weight: bold;"
    visible = false
    managed = false
  }

  private val checksumLabel = new Label() {
    style = "-fx-font-size: 11px; -fx-text-fill: #666666;"
    visible = false
    managed = false
  }

  private val progressIndicator = new ProgressIndicator() {
    prefWidth = 24
    prefHeight = 24
    visible = false
    managed = false
  }

  private val statusLabel = new Label() {
    style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
    visible = false
    managed = false
  }

  private val dropZone = new VBox(8) {
    alignment = Pos.Center
    padding = Insets(30)
    prefHeight = 180
    style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
    children = Seq(
      dropZoneLabel,
      fileNameLabel,
      detectedTypeLabel,
      new HBox(10) {
        alignment = Pos.Center
        children = Seq(progressIndicator, statusLabel)
      },
      checksumLabel
    )

    onDragOver = (event: DragEvent) => {
      if (event.dragboard.hasFiles && !isProcessing) {
        val files = event.dragboard.getFiles
        if (files.size == 1 && isSupportedFile(files.get(0))) {
          event.acceptTransferModes(TransferMode.Copy)
          style = "-fx-border-color: #4CAF50; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e8f5e9; -fx-background-radius: 8;"
        }
      }
      event.consume()
    }

    onDragExited = (_: DragEvent) => {
      if (!selectedFile.isDefined) {
        style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
      }
    }

    onDragDropped = (event: DragEvent) => {
      if (event.dragboard.hasFiles && !isProcessing) {
        val files = event.dragboard.getFiles
        if (files.size == 1 && isSupportedFile(files.get(0))) {
          handleFileSelected(files.get(0))
          event.dropCompleted = true
        }
      }
      event.consume()
    }
  }

  private val browseButton = new Button(t("action.browse")) {
    onAction = _ => {
      if (!isProcessing) {
        val fileChooser = new FileChooser() {
          title = t("data.select_file")
          extensionFilters.addAll(
            new FileChooser.ExtensionFilter("All Supported Files", Seq("*.bam", "*.cram", "*.vcf", "*.vcf.gz", "*.csv", "*.tsv", "*.txt", "*.zip")),
            new FileChooser.ExtensionFilter("Alignment Files", Seq("*.bam", "*.cram")),
            new FileChooser.ExtensionFilter("VCF Files", Seq("*.vcf", "*.vcf.gz")),
            new FileChooser.ExtensionFilter("CSV/Text Files", Seq("*.csv", "*.tsv", "*.txt")),
            new FileChooser.ExtensionFilter("All Files", "*.*")
          )
        }
        Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach(handleFileSelected)
      }
    }
  }

  private val clearButton = new Button(t("action.clear")) {
    visible = false
    managed = false
    onAction = _ => clearFileSelection()
  }

  private val buttonBox = new HBox(10) {
    alignment = Pos.Center
    children = Seq(browseButton, clearButton)
  }

  private val helpLabel = new Label(t("data.auto_detect_help")) {
    wrapText = true
    prefWidth = 380
    style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
  }

  private val contentPane = new VBox(15) {
    padding = Insets(20)
    children = Seq(dropZone, buttonBox, helpLabel)
  }

  dialogPane().content = contentPane
  dialogPane().setPrefWidth(450)

  // Add button state
  private val addButton = dialogPane().lookupButton(addButtonType)
  addButton.disable = true

  private def updateAddButton(): Unit = {
    val isValid = selectedFile.isDefined && !isProcessing && detectedType != DataType.Unknown
    addButton.disable = !isValid
  }

  private def isSupportedFile(file: File): Boolean = {
    val name = file.getName.toLowerCase
    name.endsWith(".bam") ||
      name.endsWith(".cram") ||
      name.endsWith(".vcf") ||
      name.endsWith(".vcf.gz") ||
      name.endsWith(".csv") ||
      name.endsWith(".tsv") ||
      name.endsWith(".txt") ||
      name.endsWith(".zip") ||
      name.endsWith(".csv.gz") ||
      name.endsWith(".txt.gz")
  }

  private def clearFileSelection(): Unit = {
    selectedFile = None
    detectedType = DataType.Unknown
    computedSha256 = None
    isProcessing = false

    dropZoneLabel.visible = true
    dropZoneLabel.managed = true
    fileNameLabel.visible = false
    fileNameLabel.managed = false
    detectedTypeLabel.visible = false
    detectedTypeLabel.managed = false
    checksumLabel.visible = false
    checksumLabel.managed = false
    progressIndicator.visible = false
    progressIndicator.managed = false
    statusLabel.visible = false
    statusLabel.managed = false
    clearButton.visible = false
    clearButton.managed = false

    dropZone.style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
    updateAddButton()
  }

  private def handleFileSelected(file: File): Unit = {
    selectedFile = Some(file)
    detectedType = DataType.Unknown
    computedSha256 = None
    isProcessing = true

    // Update UI to show processing state
    dropZoneLabel.visible = false
    dropZoneLabel.managed = false
    fileNameLabel.text = file.getName
    fileNameLabel.visible = true
    fileNameLabel.managed = true
    detectedTypeLabel.text = t("data.detecting_type")
    detectedTypeLabel.style = "-fx-font-size: 12px; -fx-text-fill: #888888;"
    detectedTypeLabel.visible = true
    detectedTypeLabel.managed = true
    progressIndicator.visible = true
    progressIndicator.managed = true
    statusLabel.text = t("data.analyzing_file")
    statusLabel.visible = true
    statusLabel.managed = true
    checksumLabel.visible = false
    checksumLabel.managed = false
    clearButton.visible = true
    clearButton.managed = true

    dropZone.style = "-fx-border-color: #2196F3; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e3f2fd; -fx-background-radius: 8;"
    updateAddButton()

    // Process file in background: fingerprint + checksum
    Future {
      val fingerprint = FileTypeDetector.detect(file)
      val sha256 = AnalysisCache.calculateSha256(file)
      (fingerprint, sha256)
    }.onComplete {
      case Success((fingerprint, sha256)) =>
        Platform.runLater {
          computedSha256 = Some(sha256)
          isProcessing = false
          progressIndicator.visible = false
          progressIndicator.managed = false

          // Convert fingerprint to DataType
          detectedType = fingerprint match {
            case DetectedFileType.Alignment => DataType.Alignment
            case DetectedFileType.VcfVariants => DataType.Variants
            case DetectedFileType.StrProfile => DataType.StrProfile
            case DetectedFileType.FastaMtdna => DataType.MtdnaFasta
            case DetectedFileType.ChipData(vendor) => DataType.ChipData(vendor)
            case DetectedFileType.Unknown => DataType.Unknown
          }

          // Update UI based on detection result
          if (detectedType == DataType.Unknown) {
            detectedTypeLabel.text = s"⚠ ${t("data.unknown_type")}"
            detectedTypeLabel.style = "-fx-font-size: 12px; -fx-text-fill: #F44336; -fx-font-weight: bold;"
            statusLabel.text = t("data.unknown_type_help")
            statusLabel.visible = true
            statusLabel.managed = true
            dropZone.style = "-fx-border-color: #FF9800; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #fff3e0; -fx-background-radius: 8;"
          } else {
            detectedTypeLabel.text = s"✓ ${detectedType.label}"
            detectedTypeLabel.style = "-fx-font-size: 12px; -fx-text-fill: #4CAF50; -fx-font-weight: bold;"
            statusLabel.visible = false
            statusLabel.managed = false
            dropZone.style = "-fx-border-color: #4CAF50; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e8f5e9; -fx-background-radius: 8;"
          }

          // Show checksum
          checksumLabel.visible = true
          checksumLabel.managed = true
          if (existingChecksums.contains(sha256)) {
            checksumLabel.text = s"⚠ ${t("data.duplicate")}: ${sha256.take(12)}..."
            checksumLabel.style = "-fx-font-size: 11px; -fx-text-fill: #F44336;"
          } else {
            checksumLabel.text = s"SHA256: ${sha256.take(12)}..."
            checksumLabel.style = "-fx-font-size: 11px; -fx-text-fill: #666666;"
          }

          updateAddButton()
        }

      case Failure(e) =>
        Platform.runLater {
          isProcessing = false
          progressIndicator.visible = false
          progressIndicator.managed = false
          detectedTypeLabel.text = s"⚠ ${t("error.generic")}"
          detectedTypeLabel.style = "-fx-font-size: 12px; -fx-text-fill: #F44336; -fx-font-weight: bold;"
          statusLabel.text = e.getMessage
          statusLabel.visible = true
          statusLabel.managed = true
          dropZone.style = "-fx-border-color: #F44336; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #ffebee; -fx-background-radius: 8;"
          updateAddButton()
        }
    }
  }

  private def determineFileFormat(file: File): String = {
    val name = file.getName.toLowerCase
    if (name.endsWith(".cram")) "CRAM"
    else if (name.endsWith(".bam")) "BAM"
    else if (name.endsWith(".vcf.gz")) "VCF.GZ"
    else if (name.endsWith(".vcf")) "VCF"
    else if (name.endsWith(".tsv")) "TSV"
    else if (name.endsWith(".csv") || name.endsWith(".csv.gz")) "CSV"
    else if (name.endsWith(".zip")) "ZIP"
    else "TXT"
  }

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == addButtonType && selectedFile.isDefined && detectedType != DataType.Unknown) {
      selectedFile.map { file =>
        val fileFormat = determineFileFormat(file)
        val fileInfo = FileInfo(
          fileName = file.getName,
          fileSizeBytes = Some(file.length()),
          fileFormat = fileFormat,
          checksum = computedSha256,
          checksumAlgorithm = computedSha256.map(_ => "SHA-256"),
          location = Some(file.getAbsolutePath)
        )
        DataInput(detectedType, fileInfo, computedSha256)
      }
    } else {
      None
    }
  }
}
