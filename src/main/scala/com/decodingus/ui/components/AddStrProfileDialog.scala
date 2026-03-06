package com.decodingus.ui.components

import com.decodingus.config.LabsConfig
import com.decodingus.str.StrCsvParser
import com.decodingus.str.StrCsvParser.VendorFormat
import com.decodingus.workspace.model.{RecordMeta, StrProfile}
import scalafx.Includes.*
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.*
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.{GridPane, HBox, Priority, VBox}
import scalafx.stage.FileChooser

import java.io.File

/**
 * Result from the dialog - contains the parsed STR profile
 */
case class StrProfileInput(
                            profile: StrProfile,
                            detectedFormat: VendorFormat,
                            warnings: List[String]
                          )

/**
 * Dialog for importing Y-STR data from vendor CSV files.
 * Supports FTDNA, YSEQ, and generic two-column formats.
 */
class AddStrProfileDialog(
                           biosampleRef: String,
                           existingProfileCount: Int
                         ) extends Dialog[Option[StrProfileInput]] {

  title = "Import Y-STR Profile"
  headerText = "Import Y-STR markers from a CSV file"

  val importButtonType = new ButtonType("Import", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(importButtonType, ButtonType.Cancel)
  dialogPane().setPrefSize(550, 450)

  // State
  private var selectedFile: Option[File] = None
  private var parseResult: Option[StrCsvParser.ParseResult] = None

  // Provider selection (manual override) - populated from LabsConfig
  private val strProviders = LabsConfig.strProviderNames
  private val providerCombo = new ComboBox[String] {
    items = ObservableBuffer.from("Auto-detect" +: strProviders :+ "Other")
    selectionModel().selectFirst()
    prefWidth = 180
  }

  // Drop zone for file selection
  private val dropZoneLabel = new Label("Drag & drop CSV file here\nor click Browse") {
    style = "-fx-font-size: 14px; -fx-text-fill: #888888;"
    wrapText = true
  }

  private val fileNameLabel = new Label() {
    style = "-fx-font-size: 12px; -fx-font-weight: bold;"
    visible = false
    managed = false
  }

  private val parseStatusLabel = new Label() {
    style = "-fx-font-size: 11px; -fx-text-fill: #666666;"
    visible = false
    managed = false
    wrapText = true
    maxWidth = 400
  }

  private val dropZone = new VBox(10) {
    alignment = Pos.Center
    padding = Insets(30)
    prefHeight = 150
    style = "-fx-border-color: #cccccc; -fx-border-style: dashed; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #f8f8f8; -fx-background-radius: 8;"
    children = Seq(dropZoneLabel, fileNameLabel, parseStatusLabel)

    // Drag over - accept files
    onDragOver = (event: DragEvent) => {
      if (event.dragboard.hasFiles) {
        val files = event.dragboard.getFiles
        if (files.size == 1 && isValidCsvFile(files.get(0))) {
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
        if (files.size == 1 && isValidCsvFile(files.get(0))) {
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
        title = "Select Y-STR CSV File"
        extensionFilters.addAll(
          new FileChooser.ExtensionFilter("CSV Files", Seq("*.csv", "*.CSV")),
          new FileChooser.ExtensionFilter("Text Files", Seq("*.txt", "*.TXT")),
          new FileChooser.ExtensionFilter("All Files", "*.*")
        )
      }
      Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach(handleFileSelected)
    }
  }

  // Preview of parsed markers
  private val previewLabel = new Label("") {
    style = "-fx-font-size: 11px; -fx-text-fill: #444444;"
    wrapText = true
    maxWidth = 480
  }

  private val warningsLabel = new Label("") {
    style = "-fx-font-size: 11px; -fx-text-fill: #FF9800;"
    wrapText = true
    maxWidth = 480
    visible = false
    managed = false
  }

  // Provider row
  private val providerRow = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(
      new Label("Provider:"),
      providerCombo,
      new Label("(Auto-detected from file content)") {
        style = "-fx-font-size: 10px; -fx-text-fill: #888888;"
      }
    )
  }

  // Help text
  private val helpText = new Label(
    "Supported formats:\n" +
      "  - Vendor Y-STR exports (FTDNA, YSEQ, etc.)\n" +
      "  - Two-column CSV (Marker, Value)"
  ) {
    style = "-fx-font-size: 11px; -fx-text-fill: #666666;"
    wrapText = true
  }

  // Content pane
  private val contentPane = new VBox(15) {
    padding = Insets(20)
    children = Seq(
      helpText,
      new Label("") {
        prefHeight = 5
      }, // Spacer
      providerRow,
      dropZone,
      browseButton,
      previewLabel,
      warningsLabel
    )
  }

  dialogPane().content = contentPane

  // Import button state
  private val importButton = dialogPane().lookupButton(importButtonType)
  importButton.disable = true

  private def isValidCsvFile(file: File): Boolean = {
    val name = file.getName.toLowerCase
    name.endsWith(".csv") || name.endsWith(".txt")
  }

  private def handleFileSelected(file: File): Unit = {
    selectedFile = Some(file)
    parseResult = None

    dropZoneLabel.visible = false
    dropZoneLabel.managed = false
    fileNameLabel.text = file.getName
    fileNameLabel.visible = true
    fileNameLabel.managed = true
    parseStatusLabel.text = "Parsing..."
    parseStatusLabel.visible = true
    parseStatusLabel.managed = true
    parseStatusLabel.style = "-fx-font-size: 11px; -fx-text-fill: #666666;"

    dropZone.style = "-fx-border-color: #2196F3; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e3f2fd; -fx-background-radius: 8;"

    // Parse the file
    StrCsvParser.parse(file, biosampleRef) match {
      case Right(result) =>
        parseResult = Some(result)

        // Apply manual provider override if set
        val selectedProvider = providerCombo.selectionModel().getSelectedItem
        val finalResult = selectedProvider match {
          case "Auto-detect" | "Other" => result
          case provider if strProviders.contains(provider) =>
            val updatedProfile = result.profile.copy(importedFrom = Some(provider))
            result.copy(profile = updatedProfile)
          case _ => result
        }
        parseResult = Some(finalResult)

        parseStatusLabel.text = s"Detected format: ${formatName(result.detectedFormat)}"
        parseStatusLabel.style = "-fx-font-size: 11px; -fx-text-fill: #4CAF50;"
        dropZone.style = "-fx-border-color: #4CAF50; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #e8f5e9; -fx-background-radius: 8;"

        // Show preview
        val profile = finalResult.profile
        val panelInfo = profile.panels.headOption.map(p => p.panelName).getOrElse("Unknown")
        val sampleMarkers = profile.markers.take(5).map(m => s"${m.marker}=${formatValue(m.value)}").mkString(", ")
        val moreText = if (profile.markers.size > 5) s" ... and ${profile.markers.size - 5} more" else ""
        previewLabel.text = s"Panel: $panelInfo | ${profile.markers.size} markers\nSample: $sampleMarkers$moreText"

        // Show warnings if any
        if (finalResult.warnings.nonEmpty) {
          warningsLabel.text = s"Warnings (${finalResult.warnings.size}):\n" + finalResult.warnings.take(3).mkString("\n")
          warningsLabel.visible = true
          warningsLabel.managed = true
        } else {
          warningsLabel.visible = false
          warningsLabel.managed = false
        }

        importButton.disable = false

      case Left(error) =>
        parseStatusLabel.text = s"Error: $error"
        parseStatusLabel.style = "-fx-font-size: 11px; -fx-text-fill: #F44336;"
        dropZone.style = "-fx-border-color: #F44336; -fx-border-style: solid; -fx-border-width: 2; -fx-border-radius: 8; -fx-background-color: #ffebee; -fx-background-radius: 8;"
        previewLabel.text = ""
        warningsLabel.visible = false
        warningsLabel.managed = false
        importButton.disable = true
    }
  }

  private def formatName(format: VendorFormat): String = format match {
    case VendorFormat.FTDNA => "FTDNA"
    case VendorFormat.YSEQ => "YSEQ"
    case VendorFormat.Generic => "Generic CSV"
    case VendorFormat.Unknown => "Unknown"
  }

  private def formatValue(value: com.decodingus.workspace.model.StrValue): String = {
    import com.decodingus.workspace.model.*
    value match {
      case SimpleStrValue(repeats) => repeats.toString
      case MultiCopyStrValue(copies) => copies.mkString("-")
      case ComplexStrValue(_, Some(raw)) => raw
      case ComplexStrValue(alleles, None) =>
        alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-")
    }
  }

  // Re-parse when provider selection changes
  providerCombo.selectionModel().selectedItem.onChange { (_, _, _) =>
    selectedFile.foreach(handleFileSelected)
  }

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == importButtonType) {
      parseResult.map { result =>
        StrProfileInput(result.profile, result.detectedFormat, result.warnings)
      }
    } else {
      None
    }
  }
}
