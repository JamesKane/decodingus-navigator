package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{TableView, TableColumn, Button, Alert, ButtonType, ContextMenu, MenuItem, Tooltip}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{VBox, HBox, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import scalafx.beans.property.StringProperty
import scalafx.application.Platform
import com.decodingus.analysis.CallableLociResult
import com.decodingus.workspace.model.{Biosample, SequenceRun, Alignment, AlignmentMetrics}
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.haplogroup.tree.TreeType
import com.decodingus.genotype.model.{TestTypes, TargetType}

/**
 * Table component displaying sequencing runs for a subject.
 * Supports adding new runs and triggering analysis actions.
 */
class SequenceDataTable(
  viewModel: WorkbenchViewModel,
  subject: Biosample,
  sequenceRuns: List[SequenceRun],
  alignments: List[Alignment],
  onAnalyze: (Int) => Unit,  // Callback when analyze is clicked, passes index
  onRemove: (Int) => Unit    // Callback when remove is clicked, passes index
) extends VBox(10) {

  padding = Insets(10, 0, 0, 0)

  // Helper to get alignments for a sequence run
  private def getAlignmentsForRun(run: SequenceRun): List[Alignment] = {
    run.alignmentRefs.flatMap { ref =>
      alignments.find(_.atUri.contains(ref))
    }
  }

  // Convert the sequence runs to an observable buffer with index
  case class SequenceRunRow(index: Int, run: SequenceRun, runAlignments: List[Alignment])

  private val tableData: ObservableBuffer[SequenceRunRow] = ObservableBuffer.from(
    sequenceRuns.zipWithIndex.map { case (run, idx) =>
      SequenceRunRow(idx, run, getAlignmentsForRun(run))
    }
  )

  /**
   * Calculate appropriate coverage display based on test type.
   * - WGS: genome-wide mean coverage
   * - WES: on-target coverage (higher than genome-wide)
   * - Targeted Y-DNA: Y chromosome coverage
   * - Targeted mtDNA: mtDNA coverage
   */
  private def formatCoverage(run: SequenceRun, alignment: Option[Alignment]): String = {
    val metrics = alignment.flatMap(_.metrics)
    val testType = TestTypes.byCode(run.testType)
    val targetType = testType.map(_.targetType)

    targetType match {
      case Some(TargetType.YChromosome) =>
        // For targeted Y-DNA, show Y coverage if available
        // TODO: Add yCoverage to AlignmentMetrics when we calculate per-contig coverage
        metrics.flatMap(_.meanCoverage).map(c => f"$c%.0fx").getOrElse("—")

      case Some(TargetType.MtDna) =>
        // For targeted mtDNA, show mtDNA coverage
        // TODO: Add mtCoverage to AlignmentMetrics
        metrics.flatMap(_.meanCoverage).map(c => f"$c%.0fx").getOrElse("—")

      case Some(TargetType.Autosomal) if run.testType == "WES" =>
        // For WES, coverage is on-target (typically 50-100x)
        metrics.flatMap(_.meanCoverage).map(c => f"$c%.0fx").getOrElse("—")

      case _ =>
        // WGS and others: genome-wide mean coverage
        metrics.flatMap(_.meanCoverage).map(c => f"$c%.1fx").getOrElse("—")
    }
  }

  /**
   * Get abbreviated test type display name.
   */
  private def getTestTypeAbbrev(testType: String): String = {
    SequenceRun.testTypeDisplayName(testType) match {
      case name if name == testType => testType
      case "Whole Genome Sequencing" => "WGS"
      case "Whole Exome Sequencing" => "WES"
      case "Low-Pass WGS" => "WGS-LP"
      case "PacBio HiFi WGS" => "HiFi"
      case "Nanopore WGS" => "ONT"
      case "PacBio CLR WGS" => "CLR"
      case n if n.startsWith("FTDNA Big Y") => testType
      case n if n.contains("Y Elite") => "Y Elite"
      case n if n.contains("Y-Prime") => "Y-Prime"
      case n if n.contains("mtDNA Full") => "MT Full"
      case n if n.contains("mtDNA Plus") => "MT Plus"
      case n if n.contains("Control Region") => "MT HVR"
      case _ => testType
    }
  }

  private val table = new TableView[SequenceRunRow](tableData) {
    prefHeight = 150
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Lab column (sequencing facility)
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Lab"
      cellValueFactory = { row =>
        val lab = row.value.run.sequencingFacility.getOrElse("—")
        StringProperty(lab)
      }
      prefWidth = 100
    }

    // Sample column (from BAM @RG SM tag)
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Sample"
      cellValueFactory = { row =>
        val sample = row.value.run.sampleName.getOrElse("—")
        StringProperty(sample)
      }
      prefWidth = 100
    }

    // Platform column (e.g., "Illumina", "PacBio")
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Platform"
      cellValueFactory = { row =>
        StringProperty(row.value.run.platformName)
      }
      prefWidth = 80
    }

    // Instrument Type column (e.g., "NovaSeq", "Sequel II")
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Instrument"
      cellValueFactory = { row =>
        val instrument = row.value.run.instrumentModel.getOrElse("—")
        StringProperty(instrument)
      }
      prefWidth = 90
    }

    // Library Type column (test type abbreviation)
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Library"
      cellValueFactory = { row =>
        StringProperty(getTestTypeAbbrev(row.value.run.testType))
      }
      prefWidth = 70
    }

    // Read Length + SE/PE column (e.g., "150bp PE")
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Reads"
      cellValueFactory = { row =>
        val run = row.value.run
        val readLen = run.readLength.map(r => s"${r}bp").getOrElse("")
        val layout = run.libraryLayout.map {
          case "Paired-End" => "PE"
          case "Single-End" => "SE"
          case other => other
        }.getOrElse("")
        val display = Seq(readLen, layout).filter(_.nonEmpty).mkString(" ")
        StringProperty(if (display.nonEmpty) display else "—")
      }
      prefWidth = 70
    }

    // File column - shows primary file and count if multiple
    columns += new TableColumn[SequenceRunRow, String] {
      text = "File"
      cellValueFactory = { row =>
        val files = row.value.run.files
        val display = files match {
          case Nil => "No file"
          case single :: Nil => single.fileName
          case first :: rest => s"${first.fileName} (+${rest.size})"
        }
        StringProperty(display)
      }
      prefWidth = 140
    }

    // Coverage column (smart display based on test type)
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Coverage"
      cellValueFactory = { row =>
        StringProperty(formatCoverage(row.value.run, row.value.runAlignments.headOption))
      }
      prefWidth = 65
    }

    // Reference column - shows all aligned references
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Ref"
      cellValueFactory = { row =>
        val refs = row.value.runAlignments.map(_.referenceBuild).distinct
        val display = refs match {
          case Nil => "—"
          case single :: Nil => single
          case multiple => multiple.mkString(", ")
        }
        StringProperty(display)
      }
      prefWidth = 85
    }

    // Status column
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Status"
      cellValueFactory = { row =>
        val hasMetrics = row.value.runAlignments.exists(_.metrics.isDefined)
        val status = if (hasMetrics) "Analyzed" else "Pending"
        StringProperty(status)
      }
      prefWidth = 65
    }

    // Context menu for row actions
    rowFactory = { _ =>
      val row = new javafx.scene.control.TableRow[SequenceRunRow]()
      val contextMenu = new ContextMenu(
        new MenuItem("Edit") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              handleEditSequenceRun(item.index, item.run, item.runAlignments)
            }
          }
        },
        new MenuItem("Analyze") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              onAnalyze(item.index)
            }
          }
        },
        new MenuItem("Haplogroup Analysis") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              handleHaplogroupAnalysis(item.index, item.runAlignments)
            }
          }
        },
        new MenuItem("Callable Loci") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              handleCallableLociAnalysis(item.index, item.runAlignments)
            }
          }
        },
        new MenuItem("Read/Insert Metrics") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              handleMultipleMetricsAnalysis(item.index, item.runAlignments)
            }
          }
        },
        new MenuItem("Remove") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              val confirm = new Alert(AlertType.Confirmation) {
                title = "Remove Sequencing Data"
                headerText = s"Remove this sequencing run?"
                contentText = s"Platform: ${item.run.platformName}, Test: ${item.run.testType}"
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

  /** Handles editing sequence run metadata */
  private def handleEditSequenceRun(index: Int, run: SequenceRun, runAlignments: List[Alignment]): Unit = {
    val dialog = new EditSequenceDataDialog(run, runAlignments)
    val result = dialog.showAndWait().asInstanceOf[Option[Option[SequenceRun]]]

    result match {
      case Some(Some(updatedRun)) =>
        viewModel.updateSequenceRun(subject.sampleAccession, index, updatedRun)
        // Update local table data
        tableData.update(index, SequenceRunRow(index, updatedRun, runAlignments))
      case _ => // User cancelled
    }
  }

  /** Handles launching haplogroup analysis for a sequence run */
  private def handleHaplogroupAnalysis(index: Int, runAlignments: List[Alignment]): Unit = {
    // Check if initial analysis has been run (need reference build info)
    val hasAlignments = runAlignments.nonEmpty

    if (!hasAlignments) {
      new Alert(AlertType.Warning) {
        title = "Analysis Required"
        headerText = "Initial analysis required"
        contentText = "Please run the initial analysis first to detect the reference build before running haplogroup analysis."
      }.showAndWait()
      return
    }

    // Get the test type capabilities
    val row = table.selectionModel().getSelectedItem
    val testType = row.run.testType
    val supportsY = SequenceRun.supportsYDna(testType)
    val supportsMt = SequenceRun.supportsMtDna(testType)

    // For targeted tests, auto-select the appropriate tree type
    if (supportsY && !supportsMt) {
      // Y-DNA only test (Big Y, Y Elite, etc.) - go straight to Y analysis
      runHaplogroupAnalysisForType(index, runAlignments, TreeType.YDNA)
      return
    } else if (supportsMt && !supportsY) {
      // mtDNA only test - go straight to MT analysis
      runHaplogroupAnalysisForType(index, runAlignments, TreeType.MTDNA)
      return
    } else if (!supportsY && !supportsMt) {
      // Neither supported (WES, etc.)
      new Alert(AlertType.Information) {
        title = "Haplogroup Analysis Unavailable"
        headerText = s"${SequenceRun.testTypeDisplayName(testType)} does not support haplogroup analysis"
        contentText = "This test type does not include sufficient Y-DNA or mtDNA coverage for haplogroup determination."
      }.showAndWait()
      return
    }

    // Show dialog to select tree type (for WGS and similar)
    val dialog = new HaplogroupAnalysisDialog()
    val result = dialog.showAndWait().asInstanceOf[Option[Option[TreeType]]]

    result match {
      case Some(Some(treeType)) =>
        runHaplogroupAnalysisForType(index, runAlignments, treeType)
      case _ => // User cancelled
    }
  }

  /** Runs haplogroup analysis for a specific tree type */
  private def runHaplogroupAnalysisForType(index: Int, runAlignments: List[Alignment], treeType: TreeType): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      s"${if (treeType == TreeType.YDNA) "Y-DNA" else "MT-DNA"} Haplogroup Analysis",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runHaplogroupAnalysis(
      subject.sampleAccession,
      index,
      treeType,
      onComplete = {
        case Right(haplogroupResult) =>
          Platform.runLater {
            // Show results dialog
            new HaplogroupResultDialog(
              treeType = treeType,
              haplogroupName = haplogroupResult.name,
              score = haplogroupResult.score,
              matchingSnps = haplogroupResult.matchingSnps,
              mismatchingSnps = haplogroupResult.mismatchingSnps,
              ancestralMatches = haplogroupResult.ancestralMatches,
              depth = haplogroupResult.depth
            ).showAndWait()
          }
        case Left(error) =>
          Platform.runLater {
            new Alert(AlertType.Error) {
              title = "Haplogroup Analysis Failed"
              headerText = "Could not complete haplogroup analysis"
              contentText = error
            }.showAndWait()
          }
      }
    )

    progressDialog.show()
  }

  /** Handles launching callable loci analysis for a sequence run */
  private def handleCallableLociAnalysis(index: Int, runAlignments: List[Alignment]): Unit = {
    // Check if initial analysis has been run (need reference build info)
    val hasAlignments = runAlignments.nonEmpty

    if (!hasAlignments) {
      new Alert(AlertType.Warning) {
        title = "Analysis Required"
        headerText = "Initial analysis required"
        contentText = "Please run the initial analysis first to detect the reference build before running callable loci analysis."
      }.showAndWait()
      return
    }

    val progressDialog = new AnalysisProgressDialog(
      "Callable Loci Analysis",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runCallableLociAnalysis(
      subject.sampleAccession,
      index,
      onComplete = {
        case Right((result, artifactDir)) =>
          Platform.runLater {
            // Show results dialog with artifact path for SVG viewing
            new CallableLociResultDialog(result, Some(artifactDir)).showAndWait()
          }
        case Left(error) =>
          Platform.runLater {
            new Alert(AlertType.Error) {
              title = "Callable Loci Analysis Failed"
              headerText = "Could not complete callable loci analysis"
              contentText = error
            }.showAndWait()
          }
      }
    )

    progressDialog.show()
  }

  /** Handles running CollectMultipleMetrics (read counts, insert sizes) for a sequence run */
  private def handleMultipleMetricsAnalysis(index: Int, runAlignments: List[Alignment]): Unit = {
    // Check if initial analysis has been run (need reference build info)
    val hasAlignments = runAlignments.nonEmpty

    if (!hasAlignments) {
      new Alert(AlertType.Warning) {
        title = "Analysis Required"
        headerText = "Initial analysis required"
        contentText = "Please run the initial analysis first to detect the reference build before collecting read metrics."
      }.showAndWait()
      return
    }

    val progressDialog = new AnalysisProgressDialog(
      "Collecting Read/Insert Metrics",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runMultipleMetricsAnalysis(
      subject.sampleAccession,
      index,
      onComplete = {
        case Right(metricsResult) =>
          Platform.runLater {
            // Format numbers nicely
            def fmt(n: Long): String = f"$n%,d"
            def pct(d: Double): String = f"${d * 100}%.1f%%"

            val summaryText = s"""Total Reads: ${fmt(metricsResult.totalReads)}
PF Reads: ${fmt(metricsResult.pfReads)}
Aligned: ${fmt(metricsResult.pfReadsAligned)} (${pct(metricsResult.pctPfReadsAligned)})
Paired: ${fmt(metricsResult.readsAlignedInPairs)} (${pct(metricsResult.pctReadsAlignedInPairs)})
Proper Pairs: ${pct(metricsResult.pctProperPairs)}

Read Length:
  Median: ${metricsResult.medianReadLength.toInt} bp
  Mean: ${f"${metricsResult.meanReadLength}%.1f"} bp
  Std Dev: ${f"${metricsResult.stdReadLength}%.1f"} bp
  Range: ${metricsResult.minReadLength} - ${metricsResult.maxReadLength} bp

Insert Size:
  Median: ${metricsResult.medianInsertSize.toInt} bp
  Mean: ${f"${metricsResult.meanInsertSize}%.1f"} bp
  Std Dev: ${f"${metricsResult.stdInsertSize}%.1f"} bp
  Range: ${metricsResult.minInsertSize} - ${metricsResult.maxInsertSize} bp
  Orientation: ${metricsResult.pairOrientation}"""

            new Alert(AlertType.Information) {
              title = "Read Metrics Complete"
              headerText = "Read & Insert Size Metrics"
              contentText = summaryText
            }.showAndWait()

            // Refresh the table data with updated sequence run from workspace
            viewModel.getSequenceRun(subject.sampleAccession, index).foreach { updatedRun =>
              tableData.update(index, tableData(index).copy(run = updatedRun))
            }
          }
        case Left(metricsError) =>
          Platform.runLater {
            new Alert(AlertType.Error) {
              title = "Metrics Collection Failed"
              headerText = "Could not complete metrics collection"
              contentText = metricsError
            }.showAndWait()
          }
      }
    )

    progressDialog.show()
  }

  // Action buttons
  private val addButton = new Button("Add Sequencing Run") {
    onAction = _ => {
      val existingChecksums = viewModel.getExistingChecksums(subject.sampleAccession)
      val dialog = new AddSequenceDataDialog(existingChecksums)
      val result = dialog.showAndWait().asInstanceOf[Option[Option[SequenceDataInput]]]

      result match {
        case Some(Some(input)) =>
          // Show progress dialog and run add+analyze pipeline
          val progressDialog = new AnalysisProgressDialog(
            "Adding Sequencing Data",
            viewModel.analysisProgress,
            viewModel.analysisProgressPercent,
            viewModel.analysisInProgress
          )

          viewModel.addFileAndAnalyze(
            subject.sampleAccession,
            input.fileInfo,
            onProgress = (message, _) => {
              // Progress is already bound via observable properties
            },
            onComplete = {
              case Right((index, libraryStats)) =>
                Platform.runLater {
                  new Alert(Alert.AlertType.Information) {
                    title = "Sequencing Data Added"
                    headerText = s"Successfully analyzed ${input.fileInfo.fileName}"
                    contentText = s"""Platform: ${libraryStats.inferredPlatform}
                                     |Instrument: ${libraryStats.mostFrequentInstrument}
                                     |Reference: ${libraryStats.referenceBuild}
                                     |Sample: ${libraryStats.sampleName}
                                     """.stripMargin
                  }.showAndWait()
                }
              case Left(error) =>
                Platform.runLater {
                  if (error.contains("Duplicate")) {
                    new Alert(Alert.AlertType.Warning) {
                      title = "Duplicate File"
                      headerText = "This file has already been added"
                      contentText = error
                    }.showAndWait()
                  } else {
                    new Alert(Alert.AlertType.Error) {
                      title = "Analysis Failed"
                      headerText = "Could not analyze the file"
                      contentText = error
                    }.showAndWait()
                  }
                }
            }
          )

          progressDialog.show()
        case _ => // User cancelled
      }
    }
  }

  private val editButton = new Button("Edit") {
    disable = true
    tooltip = Tooltip("Edit sequencing run metadata")
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        handleEditSequenceRun(row.index, row.run, row.runAlignments)
      }
    }
  }

  private val analyzeSelectedButton = new Button("Analyze") {
    disable = true
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        onAnalyze(row.index)
      }
    }
  }

  private val haplogroupButton = new Button("Haplogroup") {
    disable = true
    tooltip = Tooltip("Run haplogroup analysis (Y-DNA or MT-DNA)")
    onAction = _ => {
      Option(table.selectionModel().getSelectedItem).foreach { row =>
        handleHaplogroupAnalysis(row.index, row.runAlignments)
      }
    }
  }

  // Enable/disable buttons based on selection
  table.selectionModel().selectedItem.onChange { (_, _, selected) =>
    val hasSelection = selected != null
    editButton.disable = !hasSelection
    analyzeSelectedButton.disable = !hasSelection
    haplogroupButton.disable = !hasSelection
  }

  private val buttonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(addButton, editButton, analyzeSelectedButton, haplogroupButton)
  }

  children = Seq(
    new scalafx.scene.control.Label("Sequencing Runs:") { style = "-fx-font-weight: bold;" },
    table,
    buttonBar
  )

  VBox.setVgrow(table, Priority.Always)
}
