package com.decodingus.ui.components

import com.decodingus.analysis.AnalysisCache
import com.decodingus.i18n.I18n.t
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
  def extensions: Seq[String]
  def extensionFilter: FileChooser.ExtensionFilter
}

object DataType {
  case object Alignment extends DataType {
    val label = "Alignment (BAM/CRAM)"
    val extensions = Seq("*.bam", "*.cram")
    val extensionFilter = new FileChooser.ExtensionFilter("Alignment Files", extensions)
  }

  case object Variants extends DataType {
    val label = "Variants (VCF)"
    val extensions = Seq("*.vcf", "*.vcf.gz")
    val extensionFilter = new FileChooser.ExtensionFilter("VCF Files", extensions)
  }

  case object StrProfile extends DataType {
    val label = "STR Profile (CSV/TSV)"
    val extensions = Seq("*.csv", "*.tsv", "*.txt")
    val extensionFilter = new FileChooser.ExtensionFilter("STR Files", extensions)
  }

  case object ChipData extends DataType {
    val label = "Chip/Array Data (23andMe, Ancestry, etc.)"
    val extensions = Seq("*.txt", "*.csv", "*.zip")
    val extensionFilter = new FileChooser.ExtensionFilter("Chip Data Files", extensions)
  }

  val all: Seq[DataType] = Seq(Alignment, Variants, StrProfile, ChipData)
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
 * Unified dialog for adding data of various types.
 * Supports:
 * - Alignment files (BAM/CRAM)
 * - Variant files (VCF)
 * - STR profiles (CSV/TSV)
 * - Chip/Array data (23andMe, AncestryDNA, etc.)
 */
class AddDataDialog(
  existingChecksums: Set[String] = Set.empty
) extends Dialog[Option[DataInput]] {

  title = t("data.add")
  headerText = t("data.select_type")

  val addButtonType = new ButtonType(t("action.add"), ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(addButtonType, ButtonType.Cancel)

  // Data type selection
  private val dataTypeCombo = new ComboBox[DataType] {
    items = scalafx.collections.ObservableBuffer.from(DataType.all)
    cellFactory = (lv: scalafx.scene.control.ListView[DataType]) => new ListCell[DataType] {
      item.onChange { (_, _, newType) =>
        text = Option(newType).map(_.label).getOrElse("")
      }
    }
    buttonCell = new ListCell[DataType] {
      item.onChange { (_, _, newType) =>
        text = Option(newType).map(_.label).getOrElse("")
      }
    }
    selectionModel.value.selectFirst()
    prefWidth = 350
  }

  private val dataTypeSection = new VBox(5) {
    children = Seq(
      new Label(t("data.type")) { style = "-fx-font-weight: bold;" },
      dataTypeCombo
    )
  }

  // File selection
  private var selectedFile: Option[File] = None
  private var computedSha256: Option[String] = None
  private var isComputingHash = false

  private val dropZoneLabel = new Label(t("data.drag_drop_or_browse")) {
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

    onDragOver = (event: DragEvent) => {
      if (event.dragboard.hasFiles) {
        val files = event.dragboard.getFiles
        if (files.size == 1 && isValidFile(files.get(0))) {
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
        if (files.size == 1 && isValidFile(files.get(0))) {
          handleFileSelected(files.get(0))
          event.dropCompleted = true
        }
      }
      event.consume()
    }
  }

  private val browseButton = new Button(t("action.browse")) {
    onAction = _ => {
      val selectedType = dataTypeCombo.selectionModel.value.getSelectedItem
      val fileChooser = new FileChooser() {
        title = t("data.select_file")
        extensionFilters.add(selectedType.extensionFilter)
        extensionFilters.add(new FileChooser.ExtensionFilter("All Files", "*.*"))
      }
      Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach(handleFileSelected)
    }
  }

  private val fileSection = new VBox(10) {
    children = Seq(
      new Label(t("data.file")) { style = "-fx-font-weight: bold;" },
      dropZone,
      browseButton
    )
  }

  // Update drop zone text when data type changes
  dataTypeCombo.selectionModel.value.selectedItemProperty.onChange { (_, _, newType) =>
    if (newType != null) {
      dropZoneLabel.text = s"${t("data.drag_drop_or_browse")}\n(${newType.extensions.mkString(", ")})"
      // Clear any previously selected file
      clearFileSelection()
    }
  }

  private val contentPane = new VBox(20) {
    padding = Insets(20)
    children = Seq(dataTypeSection, fileSection)
  }

  dialogPane().content = contentPane
  dialogPane().setPrefWidth(450)

  // Add button state
  private val addButton = dialogPane().lookupButton(addButtonType)
  addButton.disable = true

  private def updateAddButton(): Unit = {
    val isValid = selectedFile.isDefined && !isComputingHash
    addButton.disable = !isValid
  }

  private def isValidFile(file: File): Boolean = {
    val selectedType = dataTypeCombo.selectionModel.value.getSelectedItem
    if (selectedType == null) return false

    val name = file.getName.toLowerCase
    selectedType.extensions.exists { ext =>
      val pattern = ext.replace("*", "").toLowerCase
      name.endsWith(pattern)
    }
  }

  private def clearFileSelection(): Unit = {
    selectedFile = None
    computedSha256 = None
    isComputingHash = false

    dropZoneLabel.visible = true
    dropZoneLabel.managed = true
    fileNameLabel.visible = false
    fileNameLabel.managed = false
    checksumLabel.visible = false
    checksumLabel.managed = false
    hashProgress.visible = false
    hashProgress.managed = false

    dropZone.style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
    updateAddButton()
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
    checksumLabel.text = t("data.computing_checksum")
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
            checksumLabel.text = s"⚠ ${t("data.duplicate")}: ${sha256.take(12)}..."
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
          checksumLabel.text = s"${t("data.checksum_failed")}: ${e.getMessage}"
          checksumLabel.style = "-fx-font-size: 11px; -fx-text-fill: #F44336;"
          updateAddButton()
        }
    }
  }

  private def determineFileFormat(file: File, dataType: DataType): String = {
    val name = file.getName.toLowerCase
    dataType match {
      case DataType.Alignment =>
        if (name.endsWith(".cram")) "CRAM" else "BAM"
      case DataType.Variants =>
        if (name.endsWith(".vcf.gz")) "VCF.GZ" else "VCF"
      case DataType.StrProfile =>
        if (name.endsWith(".tsv")) "TSV" else "CSV"
      case DataType.ChipData =>
        if (name.endsWith(".zip")) "ZIP"
        else if (name.endsWith(".csv")) "CSV"
        else "TXT"
    }
  }

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == addButtonType) {
      selectedFile.map { file =>
        val dataType = dataTypeCombo.selectionModel.value.getSelectedItem
        val fileFormat = determineFileFormat(file, dataType)
        val fileInfo = FileInfo(
          fileName = file.getName,
          fileSizeBytes = Some(file.length()),
          fileFormat = fileFormat,
          checksum = computedSha256,
          checksumAlgorithm = computedSha256.map(_ => "SHA-256"),
          location = Some(file.getAbsolutePath)
        )
        DataInput(dataType, fileInfo, computedSha256)
      }
    } else {
      None
    }
  }
}
