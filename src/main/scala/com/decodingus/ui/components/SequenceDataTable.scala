package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{TableView, TableColumn, Button, Alert, ButtonType, ContextMenu, MenuItem, Tooltip}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{VBox, HBox, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import scalafx.beans.property.StringProperty
import com.decodingus.workspace.model.{Biosample, SequenceData, AlignmentMetrics}
import com.decodingus.workspace.WorkbenchViewModel

/**
 * Table component displaying sequencing runs for a subject.
 * Supports adding new runs and triggering analysis actions.
 */
class SequenceDataTable(
  viewModel: WorkbenchViewModel,
  subject: Biosample,
  onAnalyze: (Int) => Unit,  // Callback when analyze is clicked, passes index
  onRemove: (Int) => Unit    // Callback when remove is clicked, passes index
) extends VBox(10) {

  padding = Insets(10, 0, 0, 0)

  // Convert the subject's sequence data to an observable buffer with index
  case class SequenceDataRow(index: Int, data: SequenceData)

  private val tableData: ObservableBuffer[SequenceDataRow] = ObservableBuffer.from(
    subject.sequenceData.zipWithIndex.map { case (sd, idx) => SequenceDataRow(idx, sd) }
  )

  private val table = new TableView[SequenceDataRow](tableData) {
    prefHeight = 200
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Platform column
    columns += new TableColumn[SequenceDataRow, String] {
      text = "Platform"
      cellValueFactory = { row => StringProperty(row.value.data.platformName) }
      prefWidth = 100
    }

    // Test Type column
    columns += new TableColumn[SequenceDataRow, String] {
      text = "Test"
      cellValueFactory = { row => StringProperty(row.value.data.testType) }
      prefWidth = 80
    }

    // File column
    columns += new TableColumn[SequenceDataRow, String] {
      text = "File"
      cellValueFactory = { row =>
        val fileName = row.value.data.files.headOption.map(_.fileName).getOrElse("No file")
        StringProperty(fileName)
      }
      prefWidth = 200
    }

    // Coverage column (from alignment metrics if available)
    columns += new TableColumn[SequenceDataRow, String] {
      text = "Mean Cov."
      cellValueFactory = { row =>
        val coverage = row.value.data.alignments.headOption
          .flatMap(_.metrics)
          .flatMap(_.meanCoverage)
          .map(c => f"$c%.1fx")
          .getOrElse("—")
        StringProperty(coverage)
      }
      prefWidth = 80
    }

    // Reference column
    columns += new TableColumn[SequenceDataRow, String] {
      text = "Reference"
      cellValueFactory = { row =>
        val ref = row.value.data.alignments.headOption
          .map(_.referenceBuild)
          .getOrElse("—")
        StringProperty(ref)
      }
      prefWidth = 80
    }

    // Status column
    columns += new TableColumn[SequenceDataRow, String] {
      text = "Status"
      cellValueFactory = { row =>
        val hasMetrics = row.value.data.alignments.exists(_.metrics.isDefined)
        val status = if (hasMetrics) "Analyzed" else "Pending"
        StringProperty(status)
      }
      prefWidth = 80
    }

    // Context menu for row actions
    rowFactory = { _ =>
      val row = new javafx.scene.control.TableRow[SequenceDataRow]()
      val contextMenu = new ContextMenu(
        new MenuItem("Analyze") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              onAnalyze(item.index)
            }
          }
        },
        new MenuItem("Remove") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              val confirm = new Alert(AlertType.Confirmation) {
                title = "Remove Sequencing Data"
                headerText = s"Remove this sequencing run?"
                contentText = s"Platform: ${item.data.platformName}, Test: ${item.data.testType}"
              }
              confirm.showAndWait() match {
                case Some(ButtonType.OK) => onRemove(item.index)
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

  // Action buttons
  private val addButton = new Button("Add Sequencing Run") {
    onAction = _ => {
      val dialog = new AddSequenceDataDialog()
      val result = dialog.showAndWait().asInstanceOf[Option[Option[SequenceData]]]

      result match {
        case Some(Some(newSeqData)) =>
          viewModel.addSequenceDataToSubject(subject.sampleAccession, newSeqData)
        case _ => // User cancelled
      }
    }
  }

  private val analyzeSelectedButton = new Button("Analyze Selected") {
    disable = true
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        onAnalyze(row.index)
      }
    }
  }

  // Enable/disable analyze button based on selection
  table.selectionModel().selectedItem.onChange { (_, _, selected) =>
    analyzeSelectedButton.disable = selected == null
  }

  private val buttonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(addButton, analyzeSelectedButton)
  }

  children = Seq(
    new scalafx.scene.control.Label("Sequencing Runs:") { style = "-fx-font-weight: bold;" },
    table,
    buttonBar
  )

  VBox.setVgrow(table, Priority.Always)
}
