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

  private val table = new TableView[SequenceRunRow](tableData) {
    prefHeight = 150
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Platform + Instrument column (e.g., "Illumina NovaSeq")
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Platform"
      cellValueFactory = { row =>
        val platform = row.value.run.platformName
        val instrument = row.value.run.instrumentModel.getOrElse("")
        val display = if (instrument.nonEmpty) s"$platform $instrument" else platform
        StringProperty(display)
      }
      prefWidth = 140
    }

    // Test Type + Read Length + Library Layout (e.g., "WGS - 150bp PE")
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Test"
      cellValueFactory = { row =>
        val run = row.value.run
        // Use display name from TestTypes if available, otherwise raw code
        val testTypeDisplay = SequenceRun.testTypeDisplayName(run.testType) match {
          case name if name == run.testType => run.testType // No mapping found, use code
          case name =>
            // For long names, use abbreviated form
            name match {
              case "Whole Genome Sequencing" => "WGS"
              case "Whole Exome Sequencing" => "WES"
              case "Low-Pass WGS" => "WGS-LP"
              case "PacBio HiFi WGS" => "HiFi WGS"
              case "Nanopore WGS" => "ONT WGS"
              case "PacBio CLR WGS" => "CLR WGS"
              case n if n.startsWith("FTDNA Big Y") => run.testType
              case n if n.contains("Y Elite") => "Y Elite"
              case n if n.contains("Y-Prime") => "Y-Prime"
              case n if n.contains("mtDNA Full") => "MT Full"
              case n if n.contains("mtDNA Plus") => "MT Plus"
              case n if n.contains("Control Region") => "MT HVR"
              case _ => run.testType
            }
        }
        val readLen = run.readLength.map(r => s"${r}bp").getOrElse("")
        val layout = run.libraryLayout.map {
          case "Paired-End" => "PE"
          case "Single-End" => "SE"
          case other => other
        }.getOrElse("")

        val details = Seq(readLen, layout).filter(_.nonEmpty).mkString(" ")
        val display = if (details.nonEmpty) s"$testTypeDisplay - $details" else testTypeDisplay
        StringProperty(display)
      }
      prefWidth = 130
    }

    // File column
    columns += new TableColumn[SequenceRunRow, String] {
      text = "File"
      cellValueFactory = { row =>
        val fileName = row.value.run.files.headOption.map(_.fileName).getOrElse("No file")
        StringProperty(fileName)
      }
      prefWidth = 180
    }

    // Coverage column (from alignment metrics if available)
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Coverage"
      cellValueFactory = { row =>
        val coverage = row.value.runAlignments.headOption
          .flatMap(_.metrics)
          .flatMap(_.meanCoverage)
          .map(c => f"$c%.1fx")
          .getOrElse("—")
        StringProperty(coverage)
      }
      prefWidth = 70
    }

    // Reference column
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Reference"
      cellValueFactory = { row =>
        val ref = row.value.runAlignments.headOption
          .map(_.referenceBuild)
          .getOrElse("—")
        StringProperty(ref)
      }
      prefWidth = 80
    }

    // Status column
    columns += new TableColumn[SequenceRunRow, String] {
      text = "Status"
      cellValueFactory = { row =>
        val hasMetrics = row.value.runAlignments.exists(_.metrics.isDefined)
        val status = if (hasMetrics) "Analyzed" else "Pending"
        StringProperty(status)
      }
      prefWidth = 70
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
