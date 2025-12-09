package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{TableView, TableColumn, Button, Alert, ButtonType, ContextMenu, MenuItem, Tooltip, Dialog, Label, TextArea, ScrollPane}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{VBox, HBox, Priority, GridPane}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import scalafx.beans.property.StringProperty
import scalafx.application.Platform
import com.decodingus.workspace.model.{Biosample, StrProfile, StrMarkerValue, StrValue, SimpleStrValue, MultiCopyStrValue, ComplexStrValue}
import com.decodingus.workspace.WorkbenchViewModel

/**
 * Table component displaying Y-STR profiles for a subject.
 * Supports importing profiles from FTDNA, YSEQ, and other vendors.
 */
class StrProfileTable(
  viewModel: WorkbenchViewModel,
  subject: Biosample,
  strProfiles: List[StrProfile],
  onRemove: (String) => Unit  // Callback when remove is clicked, passes profile URI
) extends VBox(10) {

  padding = Insets(10, 0, 0, 0)

  // Convert the profiles to an observable buffer
  case class StrProfileRow(profile: StrProfile)

  private val tableData: ObservableBuffer[StrProfileRow] = ObservableBuffer.from(
    strProfiles.map(StrProfileRow.apply)
  )

  private val table = new TableView[StrProfileRow](tableData) {
    prefHeight = 120
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Provider column
    columns += new TableColumn[StrProfileRow, String] {
      text = "Provider"
      cellValueFactory = { row =>
        StringProperty(row.value.profile.importedFrom.getOrElse(row.value.profile.source.getOrElse("Unknown")))
      }
      prefWidth = 80
    }

    // Panel column
    columns += new TableColumn[StrProfileRow, String] {
      text = "Panel"
      cellValueFactory = { row =>
        val panel = row.value.profile.panels.headOption.map(_.panelName).getOrElse("Custom")
        StringProperty(panel)
      }
      prefWidth = 80
    }

    // Marker count column
    columns += new TableColumn[StrProfileRow, String] {
      text = "Markers"
      cellValueFactory = { row =>
        StringProperty(row.value.profile.markers.size.toString)
      }
      prefWidth = 60
    }

    // Source column (imported, WGS-derived, etc.)
    columns += new TableColumn[StrProfileRow, String] {
      text = "Source"
      cellValueFactory = { row =>
        val source = row.value.profile.source match {
          case Some("IMPORTED") => "Imported"
          case Some("WGS_DERIVED") => "WGS"
          case Some("BIG_Y_DERIVED") => "Big Y"
          case Some("MANUAL_ENTRY") => "Manual"
          case Some(s) => s
          case None => "—"
        }
        StringProperty(source)
      }
      prefWidth = 70
    }

    // File column
    columns += new TableColumn[StrProfileRow, String] {
      text = "File"
      cellValueFactory = { row =>
        val fileName = row.value.profile.files.headOption.map(_.fileName).getOrElse("—")
        StringProperty(fileName)
      }
      prefWidth = 150
    }

    // Date column
    columns += new TableColumn[StrProfileRow, String] {
      text = "Added"
      cellValueFactory = { row =>
        val date = row.value.profile.meta.createdAt.toLocalDate.toString
        StringProperty(date)
      }
      prefWidth = 90
    }

    // Context menu for row actions
    rowFactory = { _ =>
      val row = new javafx.scene.control.TableRow[StrProfileRow]()
      val contextMenu = new ContextMenu(
        new MenuItem("View Markers") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              showMarkersDialog(item.profile)
            }
          }
        },
        new MenuItem("Export CSV") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              exportToCsv(item.profile)
            }
          }
        },
        new MenuItem("Remove") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              val confirm = new Alert(AlertType.Confirmation) {
                title = "Remove STR Profile"
                headerText = "Remove this STR profile?"
                contentText = s"Provider: ${item.profile.importedFrom.getOrElse("Unknown")}, ${item.profile.markers.size} markers"
              }
              confirm.showAndWait() match {
                case Some(ButtonType.OK) =>
                  item.profile.atUri.foreach(onRemove)
                case _ =>
              }
            }
          }
        }
      )
      row.contextMenu = contextMenu
      row
    }
  }

  /** Shows a dialog with all marker values */
  private def showMarkersDialog(profile: StrProfile): Unit = {
    val dialog = new Dialog[Unit] {
      title = "Y-STR Markers"
      headerText = s"${profile.importedFrom.getOrElse("STR")} Profile - ${profile.markers.size} markers"
      dialogPane().buttonTypes = Seq(ButtonType.Close)
      dialogPane().setPrefSize(500, 500)

      // Format markers as readable text
      val markerText = profile.markers
        .sortBy(_.marker)
        .map { m =>
          val value = formatStrValue(m.value)
          f"${m.marker}%-15s = $value"
        }
        .mkString("\n")

      val textArea = new TextArea(markerText) {
        editable = false
        wrapText = false
        style = "-fx-font-family: monospace; -fx-font-size: 12px;"
      }
      VBox.setVgrow(textArea, Priority.Always)

      dialogPane().content = new VBox(10) {
        padding = Insets(10)
        children = Seq(
          new Label(s"Panel: ${profile.panels.headOption.map(_.panelName).getOrElse("Custom")}"),
          textArea
        )
      }
    }
    dialog.showAndWait()
  }

  /** Formats an STR value for display */
  private def formatStrValue(value: StrValue): String = value match {
    case SimpleStrValue(repeats) => repeats.toString
    case MultiCopyStrValue(copies) => copies.mkString("-")
    case ComplexStrValue(_, Some(raw)) => raw
    case ComplexStrValue(alleles, None) =>
      alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-")
  }

  /** Exports profile to CSV (placeholder - could open save dialog) */
  private def exportToCsv(profile: StrProfile): Unit = {
    val csvContent = profile.markers
      .sortBy(_.marker)
      .map { m =>
        s"${m.marker},${formatStrValue(m.value)}"
      }
      .mkString("Marker,Value\n", "\n", "")

    // For now, show in a dialog - could be extended to save to file
    val dialog = new Dialog[Unit] {
      title = "Export CSV"
      headerText = "STR Profile CSV Export"
      dialogPane().buttonTypes = Seq(ButtonType.Close)
      dialogPane().setPrefSize(400, 400)

      val textArea = new TextArea(csvContent) {
        editable = false
        wrapText = false
        style = "-fx-font-family: monospace; -fx-font-size: 11px;"
      }

      dialogPane().content = new VBox(10) {
        padding = Insets(10)
        children = Seq(
          new Label("Copy and paste to save:"),
          textArea
        )
      }
    }
    dialog.showAndWait()
  }

  // Action buttons
  private val importButton = new Button("Import Y-STR") {
    tooltip = Tooltip("Import Y-STR profile from CSV file (FTDNA, YSEQ, etc.)")
    onAction = _ => {
      val biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")
      val dialog = new AddStrProfileDialog(biosampleRef, strProfiles.size)
      val result = dialog.showAndWait().asInstanceOf[Option[Option[StrProfileInput]]]

      result match {
        case Some(Some(input)) =>
          viewModel.addStrProfile(subject.sampleAccession, input.profile) match {
            case Right(uri) =>
              Platform.runLater {
                // Add to local table data
                tableData += StrProfileRow(input.profile.copy(atUri = Some(uri)))

                // Show success with any warnings
                val warningText = if (input.warnings.nonEmpty) {
                  s"\n\nWarnings:\n${input.warnings.take(5).mkString("\n")}"
                } else ""

                new Alert(AlertType.Information) {
                  title = "STR Profile Imported"
                  headerText = s"Successfully imported ${input.profile.markers.size} markers"
                  contentText = s"Provider: ${input.profile.importedFrom.getOrElse("Unknown")}\n" +
                    s"Panel: ${input.profile.panels.headOption.map(_.panelName).getOrElse("Custom")}" +
                    warningText
                }.showAndWait()
              }
            case Left(error) =>
              new Alert(AlertType.Error) {
                title = "Import Failed"
                headerText = "Could not import STR profile"
                contentText = error
              }.showAndWait()
          }
        case _ => // User cancelled
      }
    }
  }

  private val viewMarkersButton = new Button("View") {
    disable = true
    tooltip = Tooltip("View all marker values")
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        showMarkersDialog(row.profile)
      }
    }
  }

  // Enable/disable buttons based on selection
  table.selectionModel().selectedItem.onChange { (_, _, selected) =>
    val hasSelection = selected != null
    viewMarkersButton.disable = !hasSelection
  }

  private val buttonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(importButton, viewMarkersButton)
  }

  children = Seq(
    new scalafx.scene.control.Label("Y-STR Profiles:") { style = "-fx-font-weight: bold;" },
    table,
    buttonBar
  )

  VBox.setVgrow(table, Priority.Always)
}
