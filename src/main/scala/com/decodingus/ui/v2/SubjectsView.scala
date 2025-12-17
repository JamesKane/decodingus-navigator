package com.decodingus.ui.v2

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.i18n.Formatters
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.Biosample
import scalafx.Includes.*
import scalafx.beans.property.{ObjectProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Orientation, Pos}
import scalafx.scene.control.*
import scalafx.scene.control.cell.CheckBoxTableCell
import scalafx.scene.input.{KeyCode, MouseEvent}
import scalafx.scene.layout.*

import java.util.UUID

/**
 * Subjects view with data grid and detail panel.
 *
 * Features:
 * - Searchable/filterable data grid
 * - Multi-select for batch operations
 * - Split pane with detail view
 * - Column customization
 */
class SubjectsView(viewModel: WorkbenchViewModel) extends SplitPane {

  orientation = Orientation.Horizontal
  dividerPositions = 0.45
  styleClass += "subjects-view"

  // ============================================================================
  // State
  // ============================================================================

  private val searchText = StringProperty("")
  private val selectedSubject: ObjectProperty[Option[Biosample]] = ObjectProperty(None)
  private val selectedSubjects: ObservableBuffer[Biosample] = ObservableBuffer.empty

  // ============================================================================
  // Left Panel: Subject Grid
  // ============================================================================

  // Search bar
  private val searchField = new TextField {
    promptText = t("subjects.search.placeholder")
    prefWidth = 300
    text.onChange { (_, _, newValue) =>
      searchText.value = newValue
      applyFilter()
    }
  }

  private val addSubjectButton = new Button {
    text = t("subjects.add")
    styleClass += "button-primary"
    onAction = _ => handleAddSubject()
  }

  private val searchBar = new HBox(10) {
    alignment = Pos.CenterLeft
    padding = Insets(10)
    children = Seq(
      searchField,
      new Region { hgrow = Priority.Always },
      addSubjectButton
    )
  }

  // Data grid
  private val subjectTable = createSubjectTable()

  // Selection actions bar
  private val selectionCountLabel = new Label {
    text = ""
    styleClass += "selection-count"
  }

  private val compareButton = new Button {
    text = t("subjects.compare")
    disable = true
    onAction = _ => handleCompare()
  }

  private val batchAnalyzeButton = new Button {
    text = t("subjects.batch_analyze")
    disable = true
    onAction = _ => handleBatchAnalyze()
  }

  private val addToProjectButton = new Button {
    text = t("subjects.add_to_project")
    disable = true
    onAction = _ => handleAddToProject()
  }

  private val selectionActionsBar = new HBox(10) {
    alignment = Pos.CenterLeft
    padding = Insets(10)
    visible = false
    managed <== visible
    styleClass += "selection-actions-bar"
    children = Seq(
      selectionCountLabel,
      new Region { hgrow = Priority.Always },
      compareButton,
      batchAnalyzeButton,
      addToProjectButton
    )
  }

  private val leftPanel = new VBox {
    vgrow = Priority.Always
    children = Seq(searchBar, subjectTable, selectionActionsBar)
  }

  // ============================================================================
  // Right Panel: Detail View
  // ============================================================================

  private val detailView = new SubjectDetailView(viewModel)

  private val emptyDetailPane = new VBox {
    alignment = Pos.Center
    spacing = 10
    children = Seq(
      new Label(t("info.no_data")) {
        styleClass += "empty-state-text"
        style = "-fx-font-size: 16px; -fx-text-fill: #666666;"
      },
      new Label(t("subjects.no_results")) {
        style = "-fx-text-fill: #888888;"
      }
    )
  }

  private val rightPanel = new StackPane {
    children = Seq(emptyDetailPane, detailView)
  }

  // ============================================================================
  // Split Pane Setup
  // ============================================================================

  items.addAll(leftPanel, rightPanel)

  // ============================================================================
  // Data Binding
  // ============================================================================

  // Populate table with ViewModel data
  viewModel.samples.onChange { (_, _) =>
    applyFilter()
  }

  // Update detail view when selection changes
  selectedSubject.onChange { (_, _, newSubject) =>
    newSubject match {
      case Some(subject) =>
        detailView.setSubject(subject)
        detailView.visible = true
        emptyDetailPane.visible = false
      case None =>
        detailView.visible = false
        emptyDetailPane.visible = true
    }
  }

  // Update multi-select state
  selectedSubjects.onChange { (_, _) =>
    val count = selectedSubjects.size
    if (count > 0) {
      selectionCountLabel.text = t("subjects.selected", count)
      selectionActionsBar.visible = true
      compareButton.disable = count < 2
      batchAnalyzeButton.disable = false
      addToProjectButton.disable = false
    } else {
      selectionActionsBar.visible = false
    }
  }

  // Initial load
  applyFilter()

  // ============================================================================
  // Table Creation
  // ============================================================================

  private def createSubjectTable(): TableView[Biosample] = {
    new TableView[Biosample] {
      vgrow = Priority.Always
      placeholder = new Label(t("subjects.no_subjects"))
      styleClass += "subject-table"

      // Enable multi-select
      selectionModel.value.selectionMode = SelectionMode.Multiple

      // Columns
      columns ++= Seq(
        createColumn[String]("column.id", 100, _.accession),
        createColumn[String]("column.name", 150, s => s.donorId.getOrElse("-")),
        createColumn[String]("column.ydna", 100, s => s.yHaplogroup.getOrElse("-")),
        createColumn[String]("column.mtdna", 80, s => s.mtHaplogroup.getOrElse("-")),
        createColumn[String]("column.sex", 60, s => formatSex(s.sex)),
        createColumn[String]("column.center", 120, s => s.center.getOrElse("-")),
        createStatusColumn()
      )

      // Handle selection
      selectionModel.value.selectedItems.onChange { (buffer, _) =>
        selectedSubjects.clear()
        selectedSubjects ++= buffer.toSeq

        // Update single selection for detail view
        buffer.headOption match {
          case Some(s) => selectedSubject.value = Some(s)
          case None => selectedSubject.value = None
        }
      }

      // Double-click to open detail
      onMouseClicked = (event: MouseEvent) => {
        if (event.clickCount == 2) {
          val selected = selectionModel.value.getSelectedItem
          if (selected != null) {
            selectedSubject.value = Some(selected)
          }
        }
      }
    }
  }

  private def createColumn[T](
    headerKey: String,
    colWidth: Double,
    valueExtractor: Biosample => T
  ): TableColumn[Biosample, T] = {
    new TableColumn[Biosample, T] {
      text <== bind(headerKey)
      prefWidth = colWidth
      cellValueFactory = { cellData =>
        ObjectProperty(valueExtractor(cellData.value))
      }
    }
  }

  private def createStatusColumn(): TableColumn[Biosample, String] = {
    val col = new TableColumn[Biosample, String] {
      text <== bind("column.status")
      prefWidth = 80
      cellValueFactory = { cellData =>
        val status = determineStatus(cellData.value)
        ObjectProperty(status)
      }
    }

    col.cellFactory = { (_: TableColumn[Biosample, String]) =>
      new TableCell[Biosample, String] {
        item.onChange { (_, _, newStatus) =>
          if (newStatus != null && !newStatus.isEmpty) {
            text = newStatus
            style = newStatus match {
              case s if s == t("status.complete") => "-fx-text-fill: #4ade80;"
              case s if s == t("status.pending") => "-fx-text-fill: #fbbf24;"
              case s if s == t("status.error") => "-fx-text-fill: #f87171;"
              case _ => ""
            }
          } else {
            text = ""
            style = ""
          }
        }
      }
    }

    col
  }

  // ============================================================================
  // Filter Logic
  // ============================================================================

  private def applyFilter(): Unit = {
    val query = searchText.value.toLowerCase.trim
    val filtered = if (query.isEmpty) {
      viewModel.samples.toSeq
    } else {
      viewModel.samples.filter { s =>
        s.accession.toLowerCase.contains(query) ||
          s.donorId.exists(_.toLowerCase.contains(query)) ||
          s.yHaplogroup.exists(_.toLowerCase.contains(query)) ||
          s.mtHaplogroup.exists(_.toLowerCase.contains(query)) ||
          s.center.exists(_.toLowerCase.contains(query)) ||
          s.description.exists(_.toLowerCase.contains(query))
      }.toSeq
    }

    subjectTable.items = ObservableBuffer.from(filtered)
  }

  // ============================================================================
  // Helper Methods
  // ============================================================================

  private def formatSex(sex: Option[String]): String = {
    sex.map(_.toUpperCase) match {
      case Some("M") | Some("MALE") => t("sex.male")
      case Some("F") | Some("FEMALE") => t("sex.female")
      case _ => t("sex.unknown")
    }
  }

  private def determineStatus(subject: Biosample): String = {
    // Simple status determination based on haplogroup presence
    val hasYdna = subject.yHaplogroup.isDefined
    val hasMtdna = subject.mtHaplogroup.isDefined

    if (hasYdna || hasMtdna) t("status.complete")
    else t("status.pending")
  }

  // ============================================================================
  // Action Handlers
  // ============================================================================

  private def handleAddSubject(): Unit = {
    // TODO: Open AddSubjectDialog
    println("[SubjectsView] Add subject")
  }

  private def handleCompare(): Unit = {
    val subjects = selectedSubjects.toSeq
    if (subjects.size >= 2) {
      // TODO: Open compare view
      println(s"[SubjectsView] Compare ${subjects.size} subjects")
    }
  }

  private def handleBatchAnalyze(): Unit = {
    val subjects = selectedSubjects.toSeq
    // TODO: Run batch analysis
    println(s"[SubjectsView] Batch analyze ${subjects.size} subjects")
  }

  private def handleAddToProject(): Unit = {
    val subjects = selectedSubjects.toSeq
    // TODO: Show project picker dialog
    println(s"[SubjectsView] Add ${subjects.size} subjects to project")
  }

  // ============================================================================
  // Public API
  // ============================================================================

  /**
   * Refresh the subjects list.
   */
  def refresh(): Unit = {
    applyFilter()
  }

  /**
   * Select a subject by ID.
   */
  def selectSubject(subjectId: String): Unit = {
    viewModel.samples.find(_.id == subjectId).foreach { subject =>
      // Set the selection directly in our state
      selectedSubject.value = Some(subject)
      // Note: TableView selection sync happens through the onChange binding
    }
  }
}
