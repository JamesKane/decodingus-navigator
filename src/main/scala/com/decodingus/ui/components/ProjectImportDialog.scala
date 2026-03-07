package com.decodingus.ui.components

import com.decodingus.analysis.{DiscoveredProject, DiscoveredSample, FlagstatResult, MetricsFileLoader, ProjectDirectoryScanner}
import com.decodingus.client.{EnaClient, EnaSampleMetadata}
import com.decodingus.i18n.I18n.t
import com.decodingus.util.Logger
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.AlignmentMetrics
import scalafx.Includes.*
import scalafx.application.Platform
import scalafx.beans.property.{BooleanProperty, ObjectProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.control.cell.CheckBoxTableCell
import scalafx.scene.layout.{HBox, Priority, Region, VBox}
import scalafx.stage.DirectoryChooser

import java.io.File
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future

/**
 * Entry representing a single sample to import from a project directory.
 */
case class ProjectImportEntry(
  sampleId: String,
  discoveredSample: DiscoveredSample,
  enaMetadata: Option[EnaSampleMetadata],
  localDuplicate: Boolean,
  precomputedMetrics: Option[AlignmentMetrics],
  flagstatResult: Option[FlagstatResult],
  sampleAccession: String,
  donorIdentifier: String,
  sex: Option[String],
  population: Option[String],
  superPopulation: Option[String],
  centerName: Option[String],
  description: Option[String]
)

/**
 * Dialog for importing samples from a NAS project directory.
 *
 * Workflow:
 *   1. User selects a project directory (e.g., /Volumes/nas/Genomics/PRJEB31736)
 *   2. Scanner discovers sample subdirectories and their files
 *   3. Optional: ENA metadata resolution enriches sex, population, center
 *   4. User reviews table, deselects unwanted, confirms import
 *   5. Returns list of selected entries for the caller to import
 */
class ProjectImportDialog(viewModel: WorkbenchViewModel) extends Dialog[Option[List[ProjectImportEntry]]] {

  private val log = Logger[ProjectImportDialog]

  title = t("project.import.title")
  headerText = t("project.import.header")
  resizable = true
  dialogPane().setPrefSize(900, 620)

  private val importButtonType = new ButtonType(t("project.import.button"), ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(importButtonType, ButtonType.Cancel)

  // ============================================================================
  // State
  // ============================================================================

  private val discoveredProject = ObjectProperty[Option[DiscoveredProject]](None)
  private val rows = ObservableBuffer.empty[ImportRow]
  private val enaResolved = BooleanProperty(false)
  private val resolving = BooleanProperty(false)
  private val statusText = StringProperty("")

  // ============================================================================
  // Directory Selection
  // ============================================================================

  private val dirLabel = new Label {
    text = t("project.import.no_dir")
    style = "-fx-text-fill: #888888; -fx-font-size: 12px;"
    maxWidth = 500
  }

  private val browseButton = new Button(t("project.import.browse")) {
    onAction = _ => browseDirectory()
  }

  private val projectInfoLabel = new Label {
    style = "-fx-text-fill: #aaaaaa; -fx-font-size: 12px;"
    visible = false
    managed <== visible
  }

  private val dirBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(browseButton, dirLabel)
  }

  // ============================================================================
  // Action Bar
  // ============================================================================

  private val resolveEnaButton = new Button(t("project.import.resolve_ena")) {
    disable = true
    onAction = _ => resolveFromEna()
  }

  private val selectAllButton = new Button(t("project.import.select_all")) {
    disable = true
    onAction = _ => rows.foreach(_.selected.value = true)
  }

  private val selectNoneButton = new Button(t("project.import.select_none")) {
    disable = true
    onAction = _ => rows.foreach(_.selected.value = false)
  }

  private val selectNewButton = new Button(t("project.import.select_new")) {
    disable = true
    onAction = _ => rows.foreach(r => r.selected.value = !r.isDuplicate)
  }

  private val actionBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(
      resolveEnaButton,
      new Region { HBox.setHgrow(this, Priority.Always) },
      selectAllButton,
      selectNewButton,
      selectNoneButton
    )
  }

  // ============================================================================
  // Table
  // ============================================================================

  private val tableView = createTable()
  VBox.setVgrow(tableView, Priority.Always)

  // ============================================================================
  // Status Bar
  // ============================================================================

  private val statusLabel = new Label {
    text <== statusText
    style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
  }

  private val progressIndicator = new ProgressIndicator {
    prefWidth = 16
    prefHeight = 16
    visible <== resolving
    managed <== visible
  }

  private val statusBar = new HBox(8) {
    alignment = Pos.CenterLeft
    padding = Insets(5, 0, 0, 0)
    children = Seq(progressIndicator, statusLabel)
  }

  // ============================================================================
  // Layout
  // ============================================================================

  dialogPane().content = new VBox(10) {
    padding = Insets(10)
    children = Seq(dirBar, projectInfoLabel, actionBar, tableView, statusBar)
  }

  updateImportButton()

  // ============================================================================
  // Result Converter
  // ============================================================================

  resultConverter = dialogButton => {
    if (dialogButton == importButtonType) {
      val selected = rows.filter(_.selected.value).toList
      if (selected.nonEmpty) Some(selected.map(_.toEntry)) else None
    } else {
      None
    }
  }

  // ============================================================================
  // Table Creation
  // ============================================================================

  private def createTable(): TableView[ImportRow] = {
    new TableView[ImportRow](rows) {
      columnResizePolicy = TableView.ConstrainedResizePolicy
      placeholder = new Label(t("project.import.no_samples")) {
        style = "-fx-text-fill: #888888;"
      }
      editable = true

      columns ++= Seq(
        createCheckColumn(),
        createTextColumn(t("project.import.col.sample"), 110, _.sampleId),
        createTextColumn(t("project.import.col.accession"), 120, r =>
          r.enaAccession.getOrElse(r.sampleId)),
        createTextColumn(t("project.import.col.sex"), 50, _.sexDisplay),
        createTextColumn(t("project.import.col.population"), 60, _.populationDisplay),
        createTextColumn(t("project.import.col.files"), 80, _.filesDisplay),
        createTextColumn(t("project.import.col.metrics"), 70, _.metricsDisplay),
        createStatusColumn()
      )
    }
  }

  private def createCheckColumn(): TableColumn[ImportRow, java.lang.Boolean] = {
    val col = new TableColumn[ImportRow, java.lang.Boolean] {
      text = ""
      prefWidth = 35
      editable = true
      cellValueFactory = { cellData =>
        val prop = cellData.value.selected
        prop.asInstanceOf[ObjectProperty[java.lang.Boolean]]
      }
    }
    col.cellFactory = CheckBoxTableCell.forTableColumn(col)
    col
  }

  private def createTextColumn(header: String, colWidth: Double, extractor: ImportRow => String): TableColumn[ImportRow, String] = {
    new TableColumn[ImportRow, String] {
      text = header
      prefWidth = colWidth
      cellValueFactory = r => StringProperty(extractor(r.value))
    }
  }

  private def createStatusColumn(): TableColumn[ImportRow, String] = {
    val col = new TableColumn[ImportRow, String] {
      text = t("project.import.col.status")
      prefWidth = 80
      cellValueFactory = r => StringProperty(r.value.statusDisplay)
    }

    col.cellFactory = { (_: TableColumn[ImportRow, String]) =>
      new TableCell[ImportRow, String] {
        item.onChange { (_, _, newValue) =>
          if (newValue != null && newValue.nonEmpty) {
            text = newValue
            style = if (newValue == t("project.import.status.exists"))
              "-fx-text-fill: #fbbf24;"
            else
              "-fx-text-fill: #4ade80;"
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
  // Actions
  // ============================================================================

  private def browseDirectory(): Unit = {
    val chooser = new DirectoryChooser {
      this.title = t("project.import.choose_dir")
    }
    val selected = chooser.showDialog(dialogPane().getScene.getWindow)
    if (selected != null) {
      dirLabel.text = selected.getAbsolutePath
      scanDirectory(selected)
    }
  }

  private def scanDirectory(dir: File): Unit = {
    statusText.value = t("project.import.scanning")
    rows.clear()

    ProjectDirectoryScanner.scan(dir) match {
      case Right(project) =>
        discoveredProject.value = Some(project)

        // Check local duplicates
        val existingAccessions = viewModel.samples.map(_.sampleAccession).toSet

        // Load pre-computed metrics and build rows
        val importRows = project.samples.map { sample =>
          val isDuplicate = existingAccessions.contains(sample.sampleId)
          val precomputed = MetricsFileLoader.loadMetrics(sample)
          val flagstat = MetricsFileLoader.extractSequenceRunStats(sample)

          ImportRow(
            sampleId = sample.sampleId,
            discoveredSample = sample,
            selected = BooleanProperty(!isDuplicate),
            isDuplicate = isDuplicate,
            enaMetadata = None,
            precomputedMetrics = precomputed,
            flagstatResult = flagstat
          )
        }

        rows ++= importRows

        // Update project info
        val metricsCount = project.samplesWithMetrics
        projectInfoLabel.text = s"${project.projectId}  |  ${project.sampleCount} samples  |  " +
          s"${project.totalAlignmentFiles} alignments  |  $metricsCount with pre-computed metrics"
        projectInfoLabel.visible = true

        updateStatusSummary()
        enableButtons()

      case Left(error) =>
        statusText.value = error
        discoveredProject.value = None
    }

    updateImportButton()
  }

  private def resolveFromEna(): Unit = {
    if (discoveredProject.value.isEmpty) return ()

    resolving.value = true
    resolveEnaButton.disable = true
    statusText.value = t("project.import.resolving")

    val sampleIds = rows.map(_.sampleId).toList
    val total = sampleIds.size
    var completed = 0

    def onSampleResolved(): Unit = {
      completed += 1
      Platform.runLater {
        statusText.value = s"${t("project.import.resolving")} ($completed/$total)"
        if (completed == total) {
          val snapshot = rows.toList
          rows.clear()
          rows ++= snapshot
          resolving.value = false
          enaResolved.value = true
          resolveEnaButton.disable = false
          updateStatusSummary()
          updateImportButton()
        }
      }
    }

    // Resolve each sample from ENA
    sampleIds.foreach { sampleId =>
      EnaClient.resolveSample(sampleId).map { result =>
        Platform.runLater {
          result.foreach { meta =>
            rows.find(_.sampleId == sampleId).foreach { row =>
              row.enaMetadata = Some(meta)

              // Also check if the resolved ENA accession is a duplicate
              val existingAccessions = viewModel.samples.map(_.sampleAccession).toSet
              if (meta.sampleAccession.nonEmpty && existingAccessions.contains(meta.sampleAccession)) {
                row.isDuplicate = true
                row.selected.value = false
              }
            }
          }
        }
        onSampleResolved()
      }.recover { case _ =>
        onSampleResolved()
      }
    }
  }

  // ============================================================================
  // Helpers
  // ============================================================================

  private def enableButtons(): Unit = {
    resolveEnaButton.disable = false
    selectAllButton.disable = false
    selectNoneButton.disable = false
    selectNewButton.disable = false
  }

  private def updateImportButton(): Unit = {
    val selectedCount = rows.count(_.selected.value)
    val button = dialogPane().lookupButton(importButtonType)
    button.disable = selectedCount == 0

    // Update button text with count
    button match {
      case b: javafx.scene.control.Button =>
        b.setText(if (selectedCount > 0) s"${t("project.import.button")} ($selectedCount)" else t("project.import.button"))
      case _ =>
    }
  }

  private def updateStatusSummary(): Unit = {
    val newCount = rows.count(r => !r.isDuplicate)
    val existsCount = rows.count(_.isDuplicate)
    val metricsCount = rows.count(_.precomputedMetrics.isDefined)
    val totalAlignments = rows.map(_.discoveredSample.alignmentFiles.size).sum

    val parts = List(
      s"$newCount new",
      if (existsCount > 0) s"$existsCount existing" else "",
      s"$totalAlignments alignment files",
      if (metricsCount > 0) s"$metricsCount with metrics" else ""
    ).filter(_.nonEmpty)

    statusText.value = parts.mkString("  |  ")
  }

  // Listen for selection changes to update import button
  rows.onChange { (_, _) =>
    updateImportButton()
    rows.foreach(_.selected.onChange { (_, _, _) => updateImportButton() })
  }

  // ============================================================================
  // Row Model
  // ============================================================================

  class ImportRow(
    val sampleId: String,
    val discoveredSample: DiscoveredSample,
    val selected: BooleanProperty,
    var isDuplicate: Boolean,
    var enaMetadata: Option[EnaSampleMetadata],
    val precomputedMetrics: Option[AlignmentMetrics],
    val flagstatResult: Option[FlagstatResult]
  ) {

    def enaAccession: Option[String] = enaMetadata.map(_.sampleAccession).filter(_.nonEmpty)

    def sexDisplay: String = enaMetadata.flatMap(_.sex).map(normalizeSex).getOrElse("-")

    def populationDisplay: String = enaMetadata.flatMap(_.population).getOrElse("-")

    def filesDisplay: String = {
      val a = discoveredSample.alignmentFiles.size
      val v = discoveredSample.variantFiles.size
      val parts = List(
        if (a > 0) s"${a}A" else "",
        if (v > 0) s"${v}V" else ""
      ).filter(_.nonEmpty)
      parts.mkString(" ")
    }

    def metricsDisplay: String =
      if (precomputedMetrics.isDefined) {
        precomputedMetrics.flatMap(_.meanCoverage).map(c => f"${c}%.0fx").getOrElse("Yes")
      } else "-"

    def statusDisplay: String =
      if (isDuplicate) t("project.import.status.exists")
      else t("project.import.status.new")

    def toEntry: ProjectImportEntry = {
      val accession = enaAccession.getOrElse(sampleId)
      val donor = enaMetadata.flatMap(_.sampleAlias).getOrElse(sampleId)
      val sex = enaMetadata.flatMap(_.sex).map(normalizeSex)
      val pop = enaMetadata.flatMap(_.population)
      val superPop = enaMetadata.flatMap(_.superPopulation)
      val center = enaMetadata.flatMap(_.centerName)
      val desc = enaMetadata.flatMap(_.description).orElse {
        // Build description from available metadata
        (pop, superPop) match {
          case (Some(p), Some(sp)) =>
            val popName = enaMetadata.flatMap(_.populationName).getOrElse(p)
            Some(s"$popName ($sp)")
          case _ => None
        }
      }

      ProjectImportEntry(
        sampleId = sampleId,
        discoveredSample = discoveredSample,
        enaMetadata = enaMetadata,
        localDuplicate = isDuplicate,
        precomputedMetrics = precomputedMetrics,
        flagstatResult = flagstatResult,
        sampleAccession = accession,
        donorIdentifier = donor,
        sex = sex,
        population = pop,
        superPopulation = superPop,
        centerName = center,
        description = desc
      )
    }

    private def normalizeSex(raw: String): String = raw.toLowerCase.trim match {
      case "female" | "f" => "Female"
      case "male" | "m" => "Male"
      case _ => "Unknown"
    }
  }

  object ImportRow {
    def apply(
      sampleId: String,
      discoveredSample: DiscoveredSample,
      selected: BooleanProperty,
      isDuplicate: Boolean,
      enaMetadata: Option[EnaSampleMetadata],
      precomputedMetrics: Option[AlignmentMetrics],
      flagstatResult: Option[FlagstatResult]
    ): ImportRow = new ImportRow(sampleId, discoveredSample, selected, isDuplicate,
      enaMetadata, precomputedMetrics, flagstatResult)
  }
}
