package com.decodingus.ui.components

import com.decodingus.config.{UserPreferences, UserPreferencesService}
import com.decodingus.refgenome.config.{ReferenceConfig, ReferenceConfigService, ReferenceGenomeConfig}
import scalafx.Includes.*
import scalafx.beans.property.{BooleanProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, Priority, VBox}
import scalafx.stage.{DirectoryChooser, FileChooser}

import java.io.File
import java.nio.file.Files

/**
 * Main Settings dialog with tabs for different configuration areas.
 */
class SettingsDialog extends Dialog[Unit] {
  title = "Settings"
  headerText = "Application Settings"
  resizable = true

  dialogPane().buttonTypes = Seq(ButtonType.OK, ButtonType.Cancel)
  dialogPane().setPrefSize(800, 700)

  // Load current preferences
  private val currentPrefs = UserPreferencesService.load()

  // Tree provider selections using display names
  private val providerDisplayNames = Map(
    "FTDNA (FamilyTreeDNA)" -> "ftdna",
    "Decoding-Us" -> "decodingus"
  )
  private val providerReverseMap = providerDisplayNames.map(_.swap)

  private val ydnaProviderCombo = new ComboBox[String] {
    items = ObservableBuffer(providerDisplayNames.keys.toSeq.sorted: _*)
    value = providerReverseMap.getOrElse(currentPrefs.ydnaTreeProvider, "FTDNA (FamilyTreeDNA)")
    prefWidth = 200
  }

  private val mtdnaProviderCombo = new ComboBox[String] {
    items = ObservableBuffer(providerDisplayNames.keys.toSeq.sorted: _*)
    value = providerReverseMap.getOrElse(currentPrefs.mtdnaTreeProvider, "FTDNA (FamilyTreeDNA)")
    prefWidth = 200
  }

  // Tree Providers tab content
  private val treeProvidersContent = new VBox(20) {
    padding = Insets(20)
    children = Seq(
      new Label("Haplogroup Tree Providers") {
        style = "-fx-font-size: 16px; -fx-font-weight: bold;"
      },
      new Label("Select which haplogroup tree provider to use for Y-DNA and MT-DNA analysis.") {
        wrapText = true
        prefWidth = 500
        style = "-fx-text-fill: #666666;"
      },
      new VBox(15) {
        padding = Insets(10, 0, 0, 0)
        children = Seq(
          // Y-DNA Provider
          new VBox(5) {
            children = Seq(
              new Label("Y-DNA Tree Provider:") {
                style = "-fx-font-weight: bold;"
              },
              ydnaProviderCombo,
              new Label("  FTDNA: Uses FamilyTreeDNA's public Y-DNA haplogroup tree") {
                style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
              },
              new Label("  Decoding-Us: Uses the Decoding-Us curated tree with additional variants") {
                style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
              }
            )
          },
          // MT-DNA Provider
          new VBox(5) {
            children = Seq(
              new Label("MT-DNA Tree Provider:") {
                style = "-fx-font-weight: bold;"
              },
              mtdnaProviderCombo,
              new Label("  FTDNA: Uses FamilyTreeDNA's public MT-DNA haplogroup tree") {
                style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
              },
              new Label("  Decoding-Us: Uses the Decoding-Us curated tree") {
                style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
              }
            )
          }
        )
      },
      new Label("") {
        prefHeight = 20
      }, // Spacer
      new Label("Note: Changes take effect on the next haplogroup analysis.") {
        style = "-fx-text-fill: #888888; -fx-font-style: italic;"
      }
    )
  }

  // ============================================================================
  // Reference Genomes tab - embedded content (no separate dialog)
  // ============================================================================

  // Load current reference config
  private val refConfig = ReferenceConfigService.load()

  // Row model for the reference table
  private case class ReferenceRow(
                                   build: String,
                                   localPath: StringProperty,
                                   autoDownload: BooleanProperty,
                                   status: StringProperty
                                 )

  // Create row data from config
  private val refRowData: ObservableBuffer[ReferenceRow] = ObservableBuffer.from(
    ReferenceConfig.knownBuilds.keys.toSeq.sorted.map { build =>
      val buildConfig = refConfig.getOrDefault(build)
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
    selected = refConfig.promptBeforeDownload
    tooltip = scalafx.scene.control.Tooltip("When enabled, you'll be asked before downloading large reference files")
  }

  // Cache directory display and selector
  private val cacheDirField = new TextField {
    text = refConfig.defaultCacheDir.getOrElse(ReferenceConfigService.getCacheDir.toString)
    editable = false
    prefWidth = 350
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

  private val resetCacheDirButton = new Button("Reset") {
    onAction = _ => {
      cacheDirField.text = System.getProperty("user.home") + "/.decodingus/cache/references"
    }
  }

  // Table for reference configurations
  private val refTable = new TableView[ReferenceRow](refRowData) {
    prefHeight = 250
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
  private val refBrowseButton = new Button("Browse...") {
    disable = true
    onAction = _ => {
      Option(refTable.selectionModel().getSelectedItem).foreach { row =>
        val fileChooser = new FileChooser {
          title = s"Select ${row.build} Reference FASTA"
          extensionFilters.addAll(
            new FileChooser.ExtensionFilter("FASTA files", Seq("*.fa.gz", "*.fasta.gz", "*.fa", "*.fasta", "*.fna")),
            new FileChooser.ExtensionFilter("All files", "*.*")
          )
          if (row.localPath.value.nonEmpty) {
            val f = new File(row.localPath.value)
            if (f.getParentFile != null && f.getParentFile.exists()) {
              initialDirectory = f.getParentFile
            }
          }
        }
        Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach { file =>
          row.localPath.value = file.getAbsolutePath
          updateRefRowStatus(row)
        }
      }
    }
  }

  private val refClearPathButton = new Button("Clear Path") {
    disable = true
    onAction = _ => {
      Option(refTable.selectionModel().getSelectedItem).foreach { row =>
        row.localPath.value = ""
        updateRefRowStatus(row)
      }
    }
  }

  private val refToggleAutoDownloadButton = new Button("Toggle Auto-Download") {
    disable = true
    onAction = _ => {
      Option(refTable.selectionModel().getSelectedItem).foreach { row =>
        row.autoDownload.value = !row.autoDownload.value
        refTable.refresh()
      }
    }
  }

  // Enable/disable buttons based on selection
  refTable.selectionModel().selectedItem.onChange { (_, _, selected) =>
    val hasSelection = selected != null
    refBrowseButton.disable = !hasSelection
    refClearPathButton.disable = !hasSelection
    refToggleAutoDownloadButton.disable = !hasSelection
  }

  private def updateRefRowStatus(row: ReferenceRow): Unit = {
    val path = row.localPath.value
    if (path.nonEmpty && Files.exists(java.nio.file.Paths.get(path))) {
      row.status.value = "Local file configured"
    } else if (ReferenceConfigService.isReferenceAvailable(row.build)) {
      row.status.value = "In cache"
    } else {
      row.status.value = "Not available"
    }
    refTable.refresh()
  }

  // Layout for reference tab
  private val refTableButtonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(refBrowseButton, refClearPathButton, refToggleAutoDownloadButton)
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

  private val referenceContent = new VBox(12) {
    padding = Insets(20)
    children = Seq(
      new Label("Reference Genome Paths:") {
        style = "-fx-font-weight: bold;"
      },
      new Label("Specify local paths to reference FASTA files. If not specified, references will be downloaded to the cache directory when needed.") {
        wrapText = true
        maxWidth = 700
        style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
      },
      refTable,
      refTableButtonBar,
      new Label("") {
        prefHeight = 5
      }, // Spacer
      new Label("Download Settings:") {
        style = "-fx-font-weight: bold;"
      },
      promptCheckbox,
      new Label("") {
        prefHeight = 5
      }, // Spacer
      new Label("Cache Settings:") {
        style = "-fx-font-weight: bold;"
      },
      cacheDirRow
    )
  }

  VBox.setVgrow(refTable, Priority.Always)

  // Create tabs
  private val tabPane = new TabPane {
    tabs = Seq(
      new Tab {
        text = "Tree Providers"
        content = treeProvidersContent
        closable = false
      },
      new Tab {
        text = "Reference Genomes"
        content = referenceContent
        closable = false
      }
    )
  }

  dialogPane().content = tabPane

  // Handle OK button - save all preferences
  resultConverter = dialogButton => {
    if (dialogButton == ButtonType.OK) {
      // Save tree provider preferences
      val ydnaCode = providerDisplayNames.getOrElse(ydnaProviderCombo.value.value, "ftdna")
      val mtdnaCode = providerDisplayNames.getOrElse(mtdnaProviderCombo.value.value, "ftdna")

      val updatedPrefs = UserPreferences(
        ydnaTreeProvider = ydnaCode,
        mtdnaTreeProvider = mtdnaCode
      )
      UserPreferencesService.save(updatedPrefs) match {
        case Right(_) =>
          println(s"[SettingsDialog] Saved preferences: Y-DNA=${updatedPrefs.ydnaTreeProvider}, MT-DNA=${updatedPrefs.mtdnaTreeProvider}")
        case Left(error) =>
          println(s"[SettingsDialog] Error saving preferences: $error")
      }

      // Save reference config
      saveReferenceConfig()
    }
  }

  private def saveReferenceConfig(): Unit = {
    val references = refRowData.map { row =>
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
        println("[SettingsDialog] Reference config saved successfully")
      case Left(error) =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Failed to save reference configuration"
          contentText = error
        }.showAndWait()
    }
  }
}
