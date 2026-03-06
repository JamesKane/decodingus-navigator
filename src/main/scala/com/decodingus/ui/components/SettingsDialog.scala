package com.decodingus.ui.components

import com.decodingus.config.{UserPreferences, UserPreferencesService}
import com.decodingus.i18n.I18n.t
import com.decodingus.refgenome.config.{ReferenceConfig, ReferenceConfigService, ReferenceGenomeConfig}
import com.decodingus.ui.theme.Theme
import scalafx.Includes.*
import scalafx.beans.property.{BooleanProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, Priority, Region, VBox}
import scalafx.stage.{DirectoryChooser, FileChooser}

import java.io.File
import java.nio.file.Files

/**
 * Main Settings dialog with tabs for different configuration areas.
 * Uses dark theme styling consistent with V2 UI.
 */
class SettingsDialog extends Dialog[Unit] {
  title = t("settings.title")
  headerText = t("settings.header")
  resizable = true

  dialogPane().buttonTypes = Seq(ButtonType.OK, ButtonType.Cancel)
  dialogPane().setPrefSize(800, 700)

  // Helper to get current theme colors
  private def colors = Theme.current

  // Apply theme to dialog
  dialogPane().style = s"-fx-background-color: ${colors.background};"
  dialogPane().lookup(".header-panel").setStyle(s"-fx-background-color: ${colors.surface};")
  dialogPane().lookup(".content").setStyle(s"-fx-background-color: ${colors.background};")

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
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
  }

  private val mtdnaProviderCombo = new ComboBox[String] {
    items = ObservableBuffer(providerDisplayNames.keys.toSeq.sorted: _*)
    value = providerReverseMap.getOrElse(currentPrefs.mtdnaTreeProvider, "FTDNA (FamilyTreeDNA)")
    prefWidth = 200
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
  }

  // Tree Providers tab content
  private val treeProvidersContent = new VBox(20) {
    padding = Insets(20)
    style = s"-fx-background-color: ${colors.background};"
    children = Seq(
      new Label(t("settings.tree_providers.title")) {
        style = s"-fx-font-size: 16px; -fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
      },
      new Label(t("settings.tree_providers.description")) {
        wrapText = true
        prefWidth = 500
        style = s"-fx-text-fill: ${colors.textMuted};"
      },
      new VBox(15) {
        padding = Insets(10, 0, 0, 0)
        children = Seq(
          // Y-DNA Provider
          new VBox(5) {
            children = Seq(
              new Label(t("settings.tree_providers.ydna")) {
                style = s"-fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
              },
              ydnaProviderCombo,
              new Label(t("settings.tree_providers.ftdna_desc")) {
                style = s"-fx-text-fill: ${colors.textDisabled}; -fx-font-size: 11px;"
              },
              new Label(t("settings.tree_providers.decodingus_desc")) {
                style = s"-fx-text-fill: ${colors.textDisabled}; -fx-font-size: 11px;"
              }
            )
          },
          // MT-DNA Provider
          new VBox(5) {
            children = Seq(
              new Label(t("settings.tree_providers.mtdna")) {
                style = s"-fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
              },
              mtdnaProviderCombo,
              new Label(t("settings.tree_providers.ftdna_mt_desc")) {
                style = s"-fx-text-fill: ${colors.textDisabled}; -fx-font-size: 11px;"
              },
              new Label(t("settings.tree_providers.decodingus_mt_desc")) {
                style = s"-fx-text-fill: ${colors.textDisabled}; -fx-font-size: 11px;"
              }
            )
          }
        )
      },
      { val spacer = new Region(); spacer.prefHeight = 20; spacer }, // Spacer
      new Label(t("settings.tree_providers.note")) {
        style = s"-fx-text-fill: ${colors.textDisabled}; -fx-font-style: italic;"
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
  private val promptCheckbox = new CheckBox(t("settings.references.prompt_download")) {
    selected = refConfig.promptBeforeDownload
    tooltip = new Tooltip(t("settings.references.prompt_download_tooltip"))
    style = s"-fx-text-fill: ${colors.textPrimary};"
  }

  // Cache directory display and selector
  private val cacheDirField = new TextField {
    text = refConfig.defaultCacheDir.getOrElse(ReferenceConfigService.getCacheDir.toString)
    editable = false
    prefWidth = 350
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary}; -fx-border-color: ${colors.borderLight};"
  }

  private val browseCacheDirButton = new Button(t("settings.browse")) {
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
    onAction = _ => {
      val dirChooser = new DirectoryChooser {
        title = t("settings.references.select_cache_dir")
        initialDirectory = new File(cacheDirField.text.value)
      }
      Option(dirChooser.showDialog(dialogPane().getScene.getWindow)).foreach { dir =>
        cacheDirField.text = dir.getAbsolutePath
      }
    }
  }

  private val resetCacheDirButton = new Button(t("settings.reset")) {
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
    onAction = _ => {
      cacheDirField.text = System.getProperty("user.home") + "/.decodingus/cache/references"
    }
  }

  // Table for reference configurations
  private val refTable = new TableView[ReferenceRow](refRowData) {
    prefHeight = 250
    columnResizePolicy = TableView.ConstrainedResizePolicy
    style = s"-fx-background-color: ${colors.surface}; -fx-control-inner-background: ${colors.surface}; -fx-table-cell-border-color: ${colors.borderLight};"

    // Build column
    columns += new TableColumn[ReferenceRow, String] {
      text = t("settings.references.col_build")
      cellValueFactory = { row => StringProperty(row.value.build) }
      prefWidth = 100
      editable = false
      style = s"-fx-text-fill: ${colors.textPrimary};"
    }

    // Local Path column
    columns += new TableColumn[ReferenceRow, String] {
      text = t("settings.references.col_path")
      cellValueFactory = { row => row.value.localPath }
      prefWidth = 300
      style = s"-fx-text-fill: ${colors.textPrimary};"
    }

    // Status column
    columns += new TableColumn[ReferenceRow, String] {
      text = t("settings.references.col_status")
      cellValueFactory = { row => row.value.status }
      prefWidth = 120
      style = s"-fx-text-fill: ${colors.textPrimary};"
    }

    // Auto-download column
    columns += new TableColumn[ReferenceRow, java.lang.Boolean] {
      text = t("settings.references.col_auto_download")
      cellValueFactory = { row =>
        val objProp = new scalafx.beans.property.ObjectProperty[java.lang.Boolean]()
        objProp.value = row.value.autoDownload.value
        row.value.autoDownload.onChange { (_, _, newVal) =>
          objProp.value = newVal
        }
        objProp
      }
      prefWidth = 100
      style = s"-fx-text-fill: ${colors.textPrimary};"
    }
  }

  // Buttons for table actions
  private val refBrowseButton = new Button(t("settings.browse")) {
    disable = true
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
    onAction = _ => {
      Option(refTable.selectionModel().getSelectedItem).foreach { row =>
        val fileChooser = new FileChooser {
          title = t("settings.references.select_fasta").replace("{build}", row.build)
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

  private val refClearPathButton = new Button(t("settings.references.clear_path")) {
    disable = true
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
    onAction = _ => {
      Option(refTable.selectionModel().getSelectedItem).foreach { row =>
        row.localPath.value = ""
        updateRefRowStatus(row)
      }
    }
  }

  private val refToggleAutoDownloadButton = new Button(t("settings.references.toggle_auto_download")) {
    disable = true
    style = s"-fx-background-color: ${colors.border}; -fx-text-fill: ${colors.textPrimary};"
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
      new Label(t("settings.references.cache_dir")) {
        style = s"-fx-text-fill: ${colors.textPrimary};"
      },
      cacheDirField,
      browseCacheDirButton,
      resetCacheDirButton
    )
  }

  private val referenceContent = new VBox(12) {
    padding = Insets(20)
    style = s"-fx-background-color: ${colors.background};"
    children = Seq(
      new Label(t("settings.references.paths_title")) {
        style = s"-fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
      },
      new Label(t("settings.references.paths_description")) {
        wrapText = true
        maxWidth = 700
        style = s"-fx-text-fill: ${colors.textMuted}; -fx-font-size: 11px;"
      },
      refTable,
      refTableButtonBar,
      { val spacer = new Region(); spacer.prefHeight = 5; spacer }, // Spacer
      new Label(t("settings.references.download_settings")) {
        style = s"-fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
      },
      promptCheckbox,
      { val spacer = new Region(); spacer.prefHeight = 5; spacer }, // Spacer
      new Label(t("settings.references.cache_settings")) {
        style = s"-fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
      },
      cacheDirRow
    )
  }

  VBox.setVgrow(refTable, Priority.Always)

  // Create tabs
  private val tabPane = new TabPane {
    style = s"-fx-background-color: ${colors.background};"
    tabs = Seq(
      new Tab {
        text = t("settings.tab.tree_providers")
        content = treeProvidersContent
        closable = false
        style = s"-fx-background-color: ${colors.surface};"
      },
      new Tab {
        text = t("settings.tab.references")
        content = referenceContent
        closable = false
        style = s"-fx-background-color: ${colors.surface};"
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
          title = t("error.title")
          headerText = t("settings.references.save_failed")
          contentText = error
        }.showAndWait()
    }
  }
}
