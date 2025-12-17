package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, CheckBox, ButtonBar, Button, TableView, TableColumn, Alert}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{GridPane, VBox, HBox, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.stage.{FileChooser, DirectoryChooser}
import scalafx.collections.ObservableBuffer
import scalafx.beans.property.{StringProperty, BooleanProperty}
import com.decodingus.refgenome.config.{ReferenceConfig, ReferenceGenomeConfig, ReferenceConfigService}

import java.io.File
import java.nio.file.Files

/**
 * Row model for the reference table
 */
case class ReferenceRow(
  build: String,
  localPath: StringProperty,
  autoDownload: BooleanProperty,
  status: StringProperty
)

/**
 * Dialog for configuring reference genome paths and download settings.
 */
class ReferenceConfigDialog extends Dialog[Unit] {
  title = "Reference Genome Settings"
  headerText = "Configure local reference genome paths"
  resizable = true

  dialogPane().buttonTypes = Seq(ButtonType.OK, ButtonType.Cancel)
  dialogPane().setPrefSize(700, 700)

  // Load current config
  private val config = ReferenceConfigService.load()

  // Create row data from config
  private val rowData: ObservableBuffer[ReferenceRow] = ObservableBuffer.from(
    ReferenceConfig.knownBuilds.keys.toSeq.sorted.map { build =>
      val buildConfig = config.getOrDefault(build)
      val status = if (buildConfig.hasValidLocalPath) "Local file configured"
                   else if (ReferenceConfigService.isReferenceAvailable(build)) "In cache"
                   else "Not available"
      ReferenceRow(
        build = build,
        localPath = StringProperty(buildConfig.localPath.getOrElse("")),
        autoDownload = BooleanProperty(buildConfig.autoDownload),
        status = StringProperty(status)
      )
    }
  )

  // Prompt before download checkbox
  private val promptCheckbox = new CheckBox("Prompt before downloading references") {
    selected = config.promptBeforeDownload
    tooltip = scalafx.scene.control.Tooltip("When enabled, you'll be asked before downloading large reference files")
  }

  // Cache directory display and selector
  private val cacheDirField = new TextField {
    text = config.defaultCacheDir.getOrElse(ReferenceConfigService.getCacheDir.toString)
    editable = false
    prefWidth = 400
  }

  private val browseCacheDirButton = new Button("Browse...") {
    onAction = _ => {
      val dirChooser = new DirectoryChooser {
        title = "Select Reference Cache Directory"
        initialDirectory = new File(cacheDirField.text.value)
      }
      Option(dirChooser.showDialog(dialogPane().getScene.getWindow)).foreach { dir =>
        cacheDirField.text = dir.getAbsolutePath
      }
    }
  }

  private val resetCacheDirButton = new Button("Reset to Default") {
    onAction = _ => {
      cacheDirField.text = System.getProperty("user.home") + "/.decodingus/cache/references"
    }
  }

  // Table for reference configurations
  private val table = new TableView[ReferenceRow](rowData) {
    prefHeight = 400
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Build column
    columns += new TableColumn[ReferenceRow, String] {
      text = "Reference Build"
      cellValueFactory = { row => StringProperty(row.value.build) }
      prefWidth = 100
      editable = false
    }

    // Local Path column
    columns += new TableColumn[ReferenceRow, String] {
      text = "Local Path"
      cellValueFactory = { row => row.value.localPath }
      prefWidth = 300
    }

    // Status column
    columns += new TableColumn[ReferenceRow, String] {
      text = "Status"
      cellValueFactory = { row => row.value.status }
      prefWidth = 120
    }

    // Auto-download column
    columns += new TableColumn[ReferenceRow, java.lang.Boolean] {
      text = "Auto-Download"
      cellValueFactory = { row =>
        // Convert BooleanProperty to ObjectProperty[java.lang.Boolean]
        val objProp = new scalafx.beans.property.ObjectProperty[java.lang.Boolean]()
        objProp.value = row.value.autoDownload.value
        row.value.autoDownload.onChange { (_, _, newVal) =>
          objProp.value = newVal
        }
        objProp
      }
      prefWidth = 100
    }
  }

  // Buttons for table actions
  private val browseButton = new Button("Browse...") {
    disable = true
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        val fileChooser = new FileChooser {
          title = s"Select ${row.build} Reference FASTA"
          extensionFilters.addAll(
            new FileChooser.ExtensionFilter("FASTA files", Seq("*.fa.gz", "*.fasta.gz", "*.fa", "*.fasta")),
            new FileChooser.ExtensionFilter("All files", "*.*")
          )
          // Start in directory of current path if valid
          if (row.localPath.value.nonEmpty) {
            val f = new File(row.localPath.value)
            if (f.getParentFile != null && f.getParentFile.exists()) {
              initialDirectory = f.getParentFile
            }
          }
        }
        Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach { file =>
          row.localPath.value = file.getAbsolutePath
          updateRowStatus(row)
        }
      }
    }
  }

  private val clearPathButton = new Button("Clear Path") {
    disable = true
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        row.localPath.value = ""
        updateRowStatus(row)
      }
    }
  }

  private val toggleAutoDownloadButton = new Button("Toggle Auto-Download") {
    disable = true
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        row.autoDownload.value = !row.autoDownload.value
        table.refresh()
      }
    }
  }

  // Enable/disable buttons based on selection
  table.selectionModel().selectedItem.onChange { (_, _, selected) =>
    val hasSelection = selected != null
    browseButton.disable = !hasSelection
    clearPathButton.disable = !hasSelection
    toggleAutoDownloadButton.disable = !hasSelection
  }

  private def updateRowStatus(row: ReferenceRow): Unit = {
    val path = row.localPath.value
    if (path.nonEmpty && Files.exists(java.nio.file.Paths.get(path))) {
      row.status.value = "Local file configured"
    } else if (ReferenceConfigService.isReferenceAvailable(row.build)) {
      row.status.value = "In cache"
    } else {
      row.status.value = "Not available"
    }
    table.refresh()
  }

  // Layout
  private val tableButtonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(browseButton, clearPathButton, toggleAutoDownloadButton)
  }

  private val cacheDirRow = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(
      new Label("Cache Directory:"),
      cacheDirField,
      browseCacheDirButton,
      resetCacheDirButton
    )
  }

  private val content = new VBox(15) {
    padding = Insets(20)
    children = Seq(
      new Label("Reference Genome Paths:") { style = "-fx-font-weight: bold;" },
      new Label("Specify local paths to reference FASTA files (.fa.gz). If not specified, references will be downloaded to the cache directory when needed.") {
        wrapText = true
        maxWidth = 650
        style = "-fx-text-fill: #666666;"
      },
      table,
      tableButtonBar,
      new Label("") { prefHeight = 10 }, // Spacer
      new Label("Download Settings:") { style = "-fx-font-weight: bold;" },
      promptCheckbox,
      new Label("") { prefHeight = 10 }, // Spacer
      new Label("Cache Settings:") { style = "-fx-font-weight: bold;" },
      cacheDirRow
    )
  }

  VBox.setVgrow(table, Priority.Always)

  dialogPane().content = content

  // Result converter - save config when OK is clicked
  resultConverter = dialogButton => {
    if (dialogButton == ButtonType.OK) {
      saveConfig()
    }
  }

  private def saveConfig(): Unit = {
    // Build the new config from UI state
    val references = rowData.map { row =>
      row.build -> ReferenceGenomeConfig(
        build = row.build,
        localPath = Option(row.localPath.value).filter(_.nonEmpty),
        autoDownload = row.autoDownload.value
      )
    }.toMap

    val defaultCacheDir = {
      val defaultPath = System.getProperty("user.home") + "/.decodingus/cache/references"
      if (cacheDirField.text.value == defaultPath) None
      else Some(cacheDirField.text.value)
    }

    val newConfig = ReferenceConfig(
      references = references,
      promptBeforeDownload = promptCheckbox.selected.value,
      defaultCacheDir = defaultCacheDir
    )

    ReferenceConfigService.save(newConfig) match {
      case Right(_) =>
        println("[ReferenceConfigDialog] Config saved successfully")
      case Left(error) =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Failed to save configuration"
          contentText = error
        }.showAndWait()
    }
  }
}

/**
 * Dialog shown when a reference download is required.
 * Prompts the user to confirm the download or configure a local path.
 */
class ReferenceDownloadPromptDialog(
  build: String,
  url: String,
  estimatedSizeMB: Int
) extends Dialog[ReferenceDownloadPromptDialog.Result] {

  import ReferenceDownloadPromptDialog._

  title = "Reference Genome Required"
  headerText = s"Reference genome $build is not available locally"

  val downloadButton = new ButtonType("Download Now", ButtonBar.ButtonData.OKDone)
  val configureButton = new ButtonType("Configure Path", ButtonBar.ButtonData.Other)
  dialogPane().buttonTypes = Seq(downloadButton, configureButton, ButtonType.Cancel)

  private val content = new VBox(15) {
    padding = Insets(20)
    children = Seq(
      new Label(s"The $build reference genome is required for this analysis but is not available locally.") {
        wrapText = true
        maxWidth = 450
      },
      new Label(s"Estimated download size: ${estimatedSizeMB} MB") {
        style = "-fx-font-weight: bold;"
      },
      new Label("This may take a significant amount of time on slower connections.") {
        style = "-fx-text-fill: #666666;"
      },
      new Label("") { prefHeight = 10 },
      new Label("Options:") { style = "-fx-font-weight: bold;" },
      new Label("  - Download Now: Download the reference from the internet"),
      new Label("  - Configure Path: Specify a local file path in settings"),
      new Label("  - Cancel: Abort the current operation")
    )
  }

  dialogPane().content = content

  resultConverter = dialogButton => {
    if (dialogButton == downloadButton) Result.Download
    else if (dialogButton == configureButton) Result.Configure
    else Result.Cancel
  }
}

object ReferenceDownloadPromptDialog {
  sealed trait Result
  object Result {
    case object Download extends Result
    case object Configure extends Result
    case object Cancel extends Result
  }
}
