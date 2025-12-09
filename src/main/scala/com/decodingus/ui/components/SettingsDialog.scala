package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, RadioButton, ToggleGroup, Tab, TabPane, ComboBox, ButtonBar}
import scalafx.scene.layout.{VBox, HBox, GridPane}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import com.decodingus.config.{UserPreferences, UserPreferencesService}

/**
 * Main Settings dialog with tabs for different configuration areas.
 */
class SettingsDialog extends Dialog[Unit] {
  title = "Settings"
  headerText = "Application Settings"
  resizable = true

  dialogPane().buttonTypes = Seq(ButtonType.OK, ButtonType.Cancel)
  dialogPane().setPrefSize(750, 700)

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
              new Label("Y-DNA Tree Provider:") { style = "-fx-font-weight: bold;" },
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
              new Label("MT-DNA Tree Provider:") { style = "-fx-font-weight: bold;" },
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
      new Label("") { prefHeight = 20 }, // Spacer
      new Label("Note: Changes take effect on the next haplogroup analysis.") {
        style = "-fx-text-fill: #888888; -fx-font-style: italic;"
      }
    )
  }

  // Reference Genomes tab - embed the existing dialog content
  private val referenceContent = new VBox(10) {
    padding = Insets(20)
    children = Seq(
      new Label("Reference Genomes") {
        style = "-fx-font-size: 16px; -fx-font-weight: bold;"
      },
      new Label("Configure reference genome paths and download settings.") {
        wrapText = true
        prefWidth = 500
        style = "-fx-text-fill: #666666;"
      },
      new VBox(15) {
        padding = Insets(10, 0, 0, 0)
        children = Seq(
          new scalafx.scene.control.Button("Open Reference Genome Settings...") {
            onAction = _ => {
              val refDialog = new ReferenceConfigDialog()
              refDialog.showAndWait()
            }
          },
          new Label("Click to configure local reference genome paths, cache directory, and download settings.") {
            style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
          }
        )
      }
    )
  }

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

  // Handle OK button - save preferences
  resultConverter = dialogButton => {
    if (dialogButton == ButtonType.OK) {
      // Convert display names back to provider codes
      val ydnaCode = providerDisplayNames.getOrElse(ydnaProviderCombo.value.value, "ftdna")
      val mtdnaCode = providerDisplayNames.getOrElse(mtdnaProviderCombo.value.value, "ftdna")

      // Save tree provider preferences
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
    }
  }
}
