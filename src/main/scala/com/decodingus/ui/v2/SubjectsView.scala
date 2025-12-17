package com.decodingus.ui.v2

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.i18n.Formatters
import com.decodingus.ui.components.AddSubjectDialog
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.util.Logger
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

  private val log = Logger[SubjectsView]

  orientation = Orientation.Horizontal
  dividerPositions = 0.45
  styleClass += "subjects-view"
  style = "-fx-background-color: #1e1e1e;"

  // ============================================================================
  // Types
  // ============================================================================

  /** Represents an alignment to analyze for batch processing */
  private case class SubjectAlignment(
    subject: Biosample,
    seqRunIndex: Int,
    alignIndex: Int,
    referenceBuild: String
  )

  // ============================================================================
  // State
  // ============================================================================

  private val searchText = StringProperty("")
  private val selectedSubject: ObjectProperty[Option[Biosample]] = ObjectProperty(None)
  private val selectedSubjects: ObservableBuffer[Biosample] = ObservableBuffer.empty

  /** Helper to get the window for dialog ownership */
  private def getWindow: Option[javafx.stage.Window] =
    Option(this.getScene).flatMap(s => Option(s.getWindow))

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
    style = "-fx-background-color: #1e1e1e;"
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
    style = "-fx-background-color: #252525;"
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
    style = "-fx-background-color: #1e1e1e;"
    children = Seq(searchBar, subjectTable, selectionActionsBar)
  }

  // ============================================================================
  // Right Panel: Detail View
  // ============================================================================

  private val detailView = new SubjectDetailView(viewModel)

  private val emptyDetailPane = new VBox {
    alignment = Pos.Center
    spacing = 10
    style = "-fx-background-color: #1e1e1e;"
    children = Seq(
      new Label(t("info.no_data")) {
        styleClass += "empty-state-text"
        style = "-fx-font-size: 16px; -fx-text-fill: #888888;"
      },
      new Label(t("subjects.no_results")) {
        style = "-fx-text-fill: #666666;"
      }
    )
  }

  private val rightPanel = new StackPane {
    style = "-fx-background-color: #1e1e1e;"
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
    // Preserve current selection before replacing items
    val previouslySelectedAccession = selectedSubject.value.map(_.accession)

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

    // Restore selection if the previously selected subject is still in the filtered list
    previouslySelectedAccession.foreach { accession =>
      val indexOpt = filtered.zipWithIndex.find(_._1.accession == accession).map(_._2)
      indexOpt.foreach { index =>
        val subject = filtered(index)
        // Update with fresh data from filtered list (may have been modified)
        selectedSubject.value = Some(subject)
        // Also update table selection to keep UI in sync - use delegate to avoid ScalaFX wrapper issues
        subjectTable.delegate.getSelectionModel.clearAndSelect(index)
      }
    }
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
    val dialog = new AddSubjectDialog()
    dialog.showAndWait() match {
      case Some(Some(newSubject: Biosample)) =>
        viewModel.addSubject(newSubject)
        log.info(s"Added new subject: ${newSubject.sampleAccession}")
        // Select the newly added subject
        applyFilter()
        selectedSubject.value = Some(newSubject)
      case _ =>
        log.debug("Add subject cancelled")
    }
  }

  private def handleCompare(): Unit = {
    val subjects = selectedSubjects.toSeq
    if (subjects.size >= 2) {
      // TODO: Open compare view - requires new CompareView component
      log.debug(s"Compare ${subjects.size} subjects - not yet implemented")
      showInfoDialog(
        t("compare.title"),
        t("compare.not_implemented"),
        s"${subjects.size} ${t("subjects.selected_for_compare")}"
      )
    }
  }

  private def handleBatchAnalyze(): Unit = {
    val subjects = selectedSubjects.toSeq
    if (subjects.isEmpty) return

    // Collect alignments for each subject
    val subjectAlignments = subjects.flatMap { subject =>
      val seqRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
      seqRuns.zipWithIndex.flatMap { case (seqRun, seqRunIdx) =>
        val alignments = viewModel.workspace.value.main.getAlignmentsForSequenceRun(seqRun)
        alignments.zipWithIndex.map { case (align, alignIdx) =>
          SubjectAlignment(subject, seqRunIdx, alignIdx, align.referenceBuild)
        }
      }
    }

    if (subjectAlignments.isEmpty) {
      showInfoDialog(
        t("analysis.batch.title"),
        "No Alignments Found",
        "None of the selected subjects have alignments. Please add alignment files first."
      )
      return
    }

    // Get unique reference builds
    val availableBuilds = subjectAlignments.map(_.referenceBuild).distinct.sorted

    // If multiple reference builds, let user choose preference
    val selectedBuild: Option[String] = if (availableBuilds.size > 1) {
      val choices = Seq("All available") ++ availableBuilds
      val dialog = new scalafx.scene.control.ChoiceDialog[String](
        choices.head,
        choices
      ) {
        title = t("analysis.batch.title")
        headerText = "Select Reference Build"
        contentText = s"${subjects.size} subjects selected. Choose reference build to analyze:"
      }
      getWindow.foreach(w => dialog.initOwner(w))

      dialog.showAndWait() match {
        case Some("All available") => None // Analyze all
        case Some(build) => Some(build) // Filter to this build
        case _ => return // User cancelled
      }
    } else {
      None // Only one build available, use it
    }

    // Filter alignments based on selection
    val toAnalyze = selectedBuild match {
      case Some(build) => subjectAlignments.filter(_.referenceBuild == build)
      case None => subjectAlignments
    }

    // Group by subject to avoid analyzing same subject multiple times for same reference
    val uniqueAnalyses = toAnalyze
      .groupBy(sa => (sa.subject.sampleAccession, sa.referenceBuild))
      .values
      .map(_.head) // Take first alignment per subject/build combo
      .toSeq

    if (uniqueAnalyses.isEmpty) {
      showInfoDialog(
        t("analysis.batch.title"),
        "No Matching Alignments",
        s"No alignments found for the selected reference build."
      )
      return
    }

    // Show confirmation
    val buildSummary = uniqueAnalyses.groupBy(_.referenceBuild).map { case (build, items) =>
      s"$build: ${items.size} subject(s)"
    }.mkString("\n")

    val confirmDialog = new scalafx.scene.control.Alert(scalafx.scene.control.Alert.AlertType.Confirmation) {
      title = t("analysis.batch.title")
      headerText = s"Run comprehensive analysis on ${uniqueAnalyses.size} alignment(s)?"
      contentText = s"""$buildSummary

This will run the full analysis pipeline:
1. WGS Metrics
2. Callable Loci
3. Sex Inference
4. mtDNA Haplogroup
5. Y-DNA Haplogroup

This may take a while. Continue?"""
    }
    getWindow.foreach(w => confirmDialog.initOwner(w))

    confirmDialog.showAndWait() match {
      case Some(scalafx.scene.control.ButtonType.OK) =>
        runBatchAnalysis(uniqueAnalyses)
      case _ => // Cancelled
    }
  }

  /** Run batch analysis on multiple subject alignments */
  private def runBatchAnalysis(analyses: Seq[SubjectAlignment]): Unit = {
    import scala.concurrent.ExecutionContext.Implicits.global
    import scala.concurrent.Future

    val totalCount = analyses.size
    var completedCount = 0
    var failedCount = 0
    val results = scala.collection.mutable.ListBuffer[String]()

    // Create progress dialog
    val progressLabel = new scalafx.beans.property.StringProperty(s"Starting batch analysis of $totalCount subject(s)...")
    val progressValue = scalafx.beans.property.DoubleProperty(0.0)

    val progressDialog = new scalafx.scene.control.Dialog[Unit]() {
      title = t("analysis.batch.title")
      headerText = "Batch Analysis in Progress"
      dialogPane().content = new scalafx.scene.layout.VBox(15) {
        padding = scalafx.geometry.Insets(20)
        prefWidth = 400
        children = Seq(
          new scalafx.scene.control.Label {
            text <== progressLabel
          },
          new scalafx.scene.control.ProgressBar {
            progress <== progressValue
            prefWidth = 360
          }
        )
      }
      dialogPane().buttonTypes = Seq(scalafx.scene.control.ButtonType.Cancel)
    }
    getWindow.foreach(w => progressDialog.initOwner(w))

    // Run analyses sequentially
    def runNext(remaining: List[SubjectAlignment]): Unit = {
      remaining match {
        case Nil =>
          // All done
          scalafx.application.Platform.runLater {
            progressDialog.close()
            val summary = if (failedCount == 0) {
              s"Successfully analyzed $completedCount subject(s)."
            } else {
              s"Completed: $completedCount, Failed: $failedCount"
            }
            showInfoDialog(
              t("analysis.batch.title"),
              "Batch Analysis Complete",
              summary + (if (results.nonEmpty) "\n\n" + results.mkString("\n") else "")
            )
          }

        case current :: rest =>
          scalafx.application.Platform.runLater {
            progressLabel.value = s"Analyzing ${current.subject.donorIdentifier} (${current.referenceBuild})... (${completedCount + 1}/$totalCount)"
            progressValue.value = completedCount.toDouble / totalCount
          }

          viewModel.runComprehensiveAnalysisForAlignment(
            current.subject.sampleAccession,
            current.seqRunIndex,
            current.alignIndex,
            {
              case Right(result) =>
                completedCount += 1
                val hg = result.yDnaHaplogroup.map(_.name).orElse(result.mtDnaHaplogroup.map(_.name)).getOrElse("-")
                results += s"✓ ${current.subject.donorIdentifier}: $hg"
                runNext(rest)

              case Left(error) =>
                completedCount += 1
                failedCount += 1
                results += s"✗ ${current.subject.donorIdentifier}: $error"
                runNext(rest)
            }
          )
      }
    }

    progressDialog.show()
    runNext(analyses.toList)
  }

  private def handleAddToProject(): Unit = {
    val subjects = selectedSubjects.toSeq
    if (subjects.nonEmpty && viewModel.projects.nonEmpty) {
      // Show project picker using ChoiceDialog
      val projectNames = viewModel.projects.map(_.projectName).toSeq
      val dialog = new scalafx.scene.control.ChoiceDialog[String](
        projectNames.head,
        projectNames
      ) {
        title = t("projects.add_to")
        headerText = t("projects.select_project")
        contentText = s"${subjects.size} ${t("subjects.to_add")}:"
      }

      dialog.showAndWait() match {
        case Some(projectName) =>
          var addedCount = 0
          subjects.foreach { subject =>
            if (viewModel.addSubjectToProject(projectName, subject.accession)) {
              addedCount += 1
            }
          }
          log.info(s"Added $addedCount subjects to project: $projectName")
        case None =>
          log.debug("Add to project cancelled")
      }
    } else if (viewModel.projects.isEmpty) {
      showInfoDialog(
        t("projects.none"),
        t("projects.create_first"),
        t("projects.create_first.detail")
      )
    }
  }

  private def showInfoDialog(dialogTitle: String, dialogHeader: String, dialogContent: String): Unit = {
    import scalafx.scene.control.Alert
    import scalafx.scene.control.Alert.AlertType
    new Alert(AlertType.Information) {
      title = dialogTitle
      headerText = dialogHeader
      contentText = dialogContent
    }.showAndWait()
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
