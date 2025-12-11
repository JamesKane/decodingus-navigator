package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{TreeTableView, TreeTableColumn, TreeItem, Button, Alert, ButtonType, ContextMenu, MenuItem, Tooltip, SelectionMode}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{VBox, HBox, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import scalafx.beans.property.StringProperty
import scalafx.application.Platform
import com.decodingus.analysis.{CallableLociResult, ReadMetrics, VcfCache, VcfStatus, SubjectArtifactCache}
import com.decodingus.model.WgsMetrics
import com.decodingus.workspace.model.{Biosample, SequenceRun, Alignment, AlignmentMetrics, FileInfo}
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.haplogroup.tree.TreeType
import com.decodingus.genotype.model.{TestTypes, TargetType}
import java.nio.file.Files
import scala.jdk.CollectionConverters._

/**
 * Represents a row in the sequence data tree table.
 * Can be either a SequenceRun (parent) or an Alignment (child).
 */
sealed trait SequenceDataRow {
  def runIndex: Int
  def testType: String
}

case class SequenceRunRow(
  runIndex: Int,
  run: SequenceRun,
  alignmentCount: Int
) extends SequenceDataRow {
  def testType: String = run.testType
}

case class AlignmentRow(
  runIndex: Int,
  alignmentIndex: Int,
  run: SequenceRun,
  alignment: Alignment
) extends SequenceDataRow {
  def testType: String = run.testType
}

/**
 * Table component displaying sequencing runs for a subject with expandable alignment rows.
 * Each SequenceRun can have multiple Alignments (e.g., GRCh38, CHM13v2).
 * Analysis actions operate on specific alignments.
 */
class SequenceDataTable(
  viewModel: WorkbenchViewModel,
  subject: Biosample,
  sequenceRuns: List[SequenceRun],
  alignments: List[Alignment],
  onAnalyze: (Int) => Unit,  // Callback when analyze is clicked, passes run index
  onRemove: (Int) => Unit    // Callback when remove is clicked, passes run index
) extends VBox(10) {

  padding = Insets(10, 0, 0, 0)

  // Helper to get alignments for a sequence run
  private def getAlignmentsForRun(run: SequenceRun): List[Alignment] = {
    run.alignmentRefs.flatMap { ref =>
      alignments.find(_.atUri.contains(ref))
    }
  }

  /**
   * Build tree items for the TreeTableView.
   * SequenceRuns are parent nodes, Alignments are children.
   */
  private def buildTreeItems(): TreeItem[SequenceDataRow] = {
    val rootItem = new TreeItem[SequenceDataRow](null: SequenceDataRow)
    rootItem.setExpanded(true)

    sequenceRuns.zipWithIndex.foreach { case (run, runIdx) =>
      val runAlignments = getAlignmentsForRun(run)
      val runRow = SequenceRunRow(runIdx, run, runAlignments.size)
      val runItem = new TreeItem[SequenceDataRow](runRow)

      // Add alignment children
      runAlignments.zipWithIndex.foreach { case (alignment, alignIdx) =>
        val alignRow = AlignmentRow(runIdx, alignIdx, run, alignment)
        runItem.children += new TreeItem[SequenceDataRow](alignRow)
      }

      // Auto-expand if there are multiple alignments
      runItem.setExpanded(runAlignments.size > 1)
      rootItem.children += runItem
    }

    rootItem
  }

  /**
   * Calculate appropriate coverage display based on test type.
   */
  private def formatCoverage(row: SequenceDataRow): String = {
    val (metrics, testTypeCode) = row match {
      case sr: SequenceRunRow =>
        // For parent row, show first alignment's coverage or aggregate
        val alignments = getAlignmentsForRun(sr.run)
        (alignments.headOption.flatMap(_.metrics), sr.run.testType)
      case ar: AlignmentRow =>
        (ar.alignment.metrics, ar.run.testType)
    }

    val testType = TestTypes.byCode(testTypeCode)
    val targetType = testType.map(_.targetType)

    targetType match {
      case Some(TargetType.YChromosome) | Some(TargetType.MtDna) =>
        metrics.flatMap(_.meanCoverage).map(c => f"$c%.0fx").getOrElse("—")
      case Some(TargetType.Autosomal) if testTypeCode == "WES" =>
        metrics.flatMap(_.meanCoverage).map(c => f"$c%.0fx").getOrElse("—")
      case _ =>
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

  /**
   * Get the primary file for a row (for display and analysis).
   */
  private def getFileForRow(row: SequenceDataRow): Option[FileInfo] = {
    row match {
      case sr: SequenceRunRow =>
        // Parent row: show first alignment's file or run's file
        val alignments = getAlignmentsForRun(sr.run)
        alignments.headOption.flatMap(_.files.headOption).orElse(sr.run.files.headOption)
      case ar: AlignmentRow =>
        // Alignment row: use alignment's own file
        ar.alignment.files.headOption.orElse(ar.run.files.headOption)
    }
  }

  /**
   * Get IDs for artifact cache queries from an alignment row.
   */
  private def getArtifactIds(ar: AlignmentRow): (String, String) = {
    val runId = ar.run.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
    val alignId = ar.alignment.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
    (runId, alignId)
  }

  /**
   * Get VCF status for an alignment.
   */
  private def getVcfStatus(ar: AlignmentRow): VcfStatus = {
    val (runId, alignId) = getArtifactIds(ar)
    VcfCache.getStatus(subject.sampleAccession, runId, alignId)
  }

  /**
   * Check if callable loci data exists for an alignment.
   */
  private def hasCallableLoci(ar: AlignmentRow): Boolean = {
    val (runId, alignId) = getArtifactIds(ar)
    val callableLociDir = SubjectArtifactCache.getArtifactSubdir(subject.sampleAccession, runId, alignId, "callable_loci")
    if (Files.exists(callableLociDir)) {
      // Check if there are any .bed files
      Files.list(callableLociDir).iterator().asScala.exists(_.toString.endsWith(".callable.bed"))
    } else {
      false
    }
  }

  /**
   * Format VCF status for display.
   */
  private def formatVcfStatus(status: VcfStatus): String = {
    status match {
      case VcfStatus.Available(_) => "✓"
      case VcfStatus.InProgress(_, progress, _) => f"◐ ${progress * 100}%.0f%%"
      case VcfStatus.NotGenerated => "○"
      case VcfStatus.Incomplete => "⚠"
      case VcfStatus.Stale => "⚠"
    }
  }

  /**
   * Get tooltip text for VCF status.
   */
  private def getVcfTooltip(ar: AlignmentRow): String = {
    getVcfStatus(ar) match {
      case VcfStatus.Available(info) =>
        val sizeGb = info.fileSizeBytes / (1024.0 * 1024.0 * 1024.0)
        s"""Whole-Genome VCF
Reference: ${info.referenceBuild}
Variants: ${f"${info.variantCount}%,d"}
Size: ${f"$sizeGb%.2f"} GB
Created: ${info.createdAt}
Sex: ${info.inferredSex.getOrElse("Unknown")}"""
      case VcfStatus.InProgress(startedAt, progress, currentContig) =>
        val contigInfo = currentContig.map(c => s"\nProcessing: $c").getOrElse("")
        s"VCF Generation in Progress\nStarted: $startedAt\nProgress: ${f"${progress * 100}%.0f"}%$contigInfo"
      case VcfStatus.NotGenerated =>
        "No VCF Generated\nRight-click to generate whole-genome VCF"
      case VcfStatus.Incomplete =>
        "VCF Incomplete\nMetadata missing or files corrupted"
      case VcfStatus.Stale =>
        "VCF May Be Stale\nAlignment modified since VCF was generated"
    }
  }

  /**
   * Get tooltip text for callable loci status.
   */
  private def getCallableLociTooltip(ar: AlignmentRow): String = {
    if (hasCallableLoci(ar)) {
      "Callable Loci Available\nBED files generated for this alignment"
    } else {
      "No Callable Loci\nRight-click to run callable loci analysis"
    }
  }

  private val treeTable = new TreeTableView[SequenceDataRow](buildTreeItems()) {
    prefHeight = 200
    showRoot = false
    columnResizePolicy = TreeTableView.ConstrainedResizePolicy

    // Lab column (sequencing facility) - only for SequenceRun rows
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Lab"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow => sr.run.sequencingFacility.getOrElse("—")
          case _: AlignmentRow => ""  // Empty for alignment rows
        }
        StringProperty(value)
      }
      prefWidth = 80
    }

    // Sample column (from BAM @RG SM tag) - only for SequenceRun rows
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Sample"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow => sr.run.sampleName.getOrElse("—")
          case _: AlignmentRow => ""
        }
        StringProperty(value)
      }
      prefWidth = 80
    }

    // Platform column - only for SequenceRun rows
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Platform"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow => sr.run.platformName
          case _: AlignmentRow => ""
        }
        StringProperty(value)
      }
      prefWidth = 70
    }

    // Library Type column - only for SequenceRun rows
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Library"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow => getTestTypeAbbrev(sr.run.testType)
          case _: AlignmentRow => ""
        }
        StringProperty(value)
      }
      prefWidth = 60
    }

    // Reference column - shows reference build
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Reference"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow =>
            val alignments = getAlignmentsForRun(sr.run)
            if (alignments.isEmpty) "—"
            else if (alignments.size == 1) alignments.head.referenceBuild
            else s"${alignments.size} builds"
          case ar: AlignmentRow => ar.alignment.referenceBuild
        }
        StringProperty(value)
      }
      prefWidth = 80
    }

    // File column - shows alignment file
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "File"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow =>
            val alignments = getAlignmentsForRun(sr.run)
            if (alignments.isEmpty) {
              sr.run.files.headOption.map(_.fileName).getOrElse("No file")
            } else if (alignments.size == 1) {
              alignments.head.files.headOption.map(_.fileName).getOrElse(
                sr.run.files.headOption.map(_.fileName).getOrElse("No file")
              )
            } else {
              s"${alignments.size} alignments"
            }
          case ar: AlignmentRow =>
            ar.alignment.files.headOption.map(_.fileName).getOrElse(
              ar.run.files.headOption.map(_.fileName).getOrElse("—")
            )
        }
        StringProperty(value)
      }
      prefWidth = 150
    }

    // Coverage column
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Coverage"
      cellValueFactory = { p =>
        StringProperty(formatCoverage(p.value.getValue))
      }
      prefWidth = 65
    }

    // Status column
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "Status"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case sr: SequenceRunRow =>
            val alignments = getAlignmentsForRun(sr.run)
            if (alignments.isEmpty) "Pending"
            else if (alignments.exists(_.metrics.isDefined)) "Analyzed"
            else "Pending"
          case ar: AlignmentRow =>
            if (ar.alignment.metrics.isDefined) "Analyzed" else "Pending"
        }
        StringProperty(value)
      }
      prefWidth = 60
    }

    // VCF status column (for alignment rows only)
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "VCF"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case _: SequenceRunRow => ""  // Empty for parent rows
          case ar: AlignmentRow => formatVcfStatus(getVcfStatus(ar))
        }
        StringProperty(value)
      }
      prefWidth = 40
      style = "-fx-alignment: CENTER;"
    }

    // Callable Loci status column (for alignment rows only)
    columns += new TreeTableColumn[SequenceDataRow, String] {
      text = "CL"
      cellValueFactory = { p =>
        val value = p.value.getValue match {
          case _: SequenceRunRow => ""  // Empty for parent rows
          case ar: AlignmentRow => if (hasCallableLoci(ar)) "✓" else "○"
        }
        StringProperty(value)
      }
      prefWidth = 30
      style = "-fx-alignment: CENTER;"
    }

    // Row factory for context menus
    rowFactory = { _ =>
      val row = new javafx.scene.control.TreeTableRow[SequenceDataRow]()

      row.itemProperty().addListener { (_, _, newItem) =>
        if (newItem != null) {
          row.setContextMenu(createContextMenu(newItem).delegate)
        } else {
          row.setContextMenu(null)
        }
      }

      row
    }
  }

  /**
   * Create context menu based on row type.
   */
  private def createContextMenu(row: SequenceDataRow): ContextMenu = {
    row match {
      case sr: SequenceRunRow => createSequenceRunContextMenu(sr)
      case ar: AlignmentRow => createAlignmentContextMenu(ar)
    }
  }

  /**
   * Context menu for SequenceRun rows (parent).
   */
  private def createSequenceRunContextMenu(sr: SequenceRunRow): ContextMenu = {
    val runAlignments = getAlignmentsForRun(sr.run)

    new ContextMenu(
      new MenuItem("Edit") {
        onAction = _ => handleEditSequenceRun(sr.runIndex, sr.run, runAlignments)
      },
      new MenuItem("Analyze") {
        onAction = _ => onAnalyze(sr.runIndex)
      },
      new MenuItem("Remove") {
        onAction = _ => {
          val confirm = new Alert(AlertType.Confirmation) {
            title = "Remove Sequencing Data"
            headerText = s"Remove this sequencing run and all alignments?"
            contentText = s"Platform: ${sr.run.platformName}, Test: ${sr.run.testType}"
          }
          confirm.showAndWait() match {
            case Some(ButtonType.OK) => onRemove(sr.runIndex)
            case _ =>
          }
        }
      }
    )
  }

  /**
   * Context menu for Alignment rows (child) - includes analysis actions.
   */
  private def createAlignmentContextMenu(ar: AlignmentRow): ContextMenu = {
    val vcfStatus = getVcfStatus(ar)
    val vcfMenuItem = new MenuItem("Generate Whole-Genome VCF") {
      onAction = _ => handleGenerateVcf(ar)
      disable = vcfStatus.isAvailable || vcfStatus.isInProgress
    }

    val viewVcfStatsMenuItem = new MenuItem("View VCF Statistics") {
      onAction = _ => handleViewVcfStats(ar)
      disable = !vcfStatus.isAvailable
    }

    new ContextMenu(
      vcfMenuItem,
      viewVcfStatsMenuItem,
      new javafx.scene.control.SeparatorMenuItem(),
      new MenuItem("Haplogroup Analysis") {
        onAction = _ => handleHaplogroupAnalysis(ar.runIndex, ar.alignmentIndex, ar.alignment)
      },
      new MenuItem("View Y-DNA Report") {
        onAction = _ => showCachedHaplogroupReport(ar.runIndex, ar.alignmentIndex, TreeType.YDNA)
      },
      new MenuItem("View mtDNA Report") {
        onAction = _ => showCachedHaplogroupReport(ar.runIndex, ar.alignmentIndex, TreeType.MTDNA)
      },
      new javafx.scene.control.SeparatorMenuItem(),
      new MenuItem("Callable Loci") {
        onAction = _ => handleCallableLociAnalysis(ar.runIndex, ar.alignmentIndex, ar.alignment)
      },
      new MenuItem("WGS Metrics") {
        onAction = _ => handleWgsMetricsAnalysis(ar.runIndex, ar.alignmentIndex, ar.alignment)
      },
      new MenuItem("Read/Insert Metrics") {
        onAction = _ => handleMultipleMetricsAnalysis(ar.runIndex, ar.alignmentIndex, ar.alignment)
      }
    )
  }

  /** Handles editing sequence run metadata */
  private def handleEditSequenceRun(index: Int, run: SequenceRun, runAlignments: List[Alignment]): Unit = {
    val dialog = new EditSequenceDataDialog(run, runAlignments)
    val result = dialog.showAndWait().asInstanceOf[Option[Option[SequenceRun]]]

    result match {
      case Some(Some(updatedRun)) =>
        viewModel.updateSequenceRun(subject.sampleAccession, index, updatedRun)
        // Rebuild tree to reflect changes
        treeTable.root = buildTreeItems()
      case _ => // User cancelled
    }
  }

  /** Shows a cached haplogroup report if it exists */
  private def showCachedHaplogroupReport(runIndex: Int, alignmentIndex: Int, treeType: TreeType): Unit = {
    val dnaType = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"
    val prefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"

    viewModel.getHaplogroupArtifactDirForAlignment(subject.sampleAccession, runIndex, alignmentIndex) match {
      case Some(dir) if Files.exists(dir.resolve(s"${prefix}_haplogroup_report.txt")) =>
        new HaplogroupReportDialog(
          treeType = treeType,
          artifactDir = Some(dir),
          sampleName = Some(subject.donorIdentifier)
        ).showAndWait()

      case _ =>
        new Alert(AlertType.Information) {
          title = "No Report Available"
          headerText = s"No $dnaType haplogroup report found"
          contentText = "Run haplogroup analysis first to generate a report."
        }.showAndWait()
    }
  }

  /** Handles launching haplogroup analysis for a specific alignment */
  private def handleHaplogroupAnalysis(runIndex: Int, alignmentIndex: Int, alignment: Alignment): Unit = {
    // Get the test type capabilities from the parent run
    viewModel.getSequenceRun(subject.sampleAccession, runIndex) match {
      case None =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Sequence run not found"
        }.showAndWait()
        return

      case Some(run) =>
        val testType = run.testType
        val supportsY = SequenceRun.supportsYDna(testType)
        val supportsMt = SequenceRun.supportsMtDna(testType)

        // For targeted tests, auto-select the appropriate tree type
        if (supportsY && !supportsMt) {
          runHaplogroupAnalysisForAlignment(runIndex, alignmentIndex, alignment, TreeType.YDNA)
          return
        } else if (supportsMt && !supportsY) {
          runHaplogroupAnalysisForAlignment(runIndex, alignmentIndex, alignment, TreeType.MTDNA)
          return
        } else if (!supportsY && !supportsMt) {
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
            runHaplogroupAnalysisForAlignment(runIndex, alignmentIndex, alignment, treeType)
          case _ => // User cancelled
        }
    }
  }

  /** Runs haplogroup analysis for a specific alignment */
  private def runHaplogroupAnalysisForAlignment(runIndex: Int, alignmentIndex: Int, alignment: Alignment, treeType: TreeType): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      s"${if (treeType == TreeType.YDNA) "Y-DNA" else "MT-DNA"} Haplogroup Analysis (${alignment.referenceBuild})",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runHaplogroupAnalysisForAlignment(
      subject.sampleAccession,
      runIndex,
      alignmentIndex,
      treeType,
      onComplete = {
        case Right(haplogroupResult) =>
          Platform.runLater {
            val artifactDir = viewModel.getHaplogroupArtifactDirForAlignment(subject.sampleAccession, runIndex, alignmentIndex)

            val hasDetailedReport = artifactDir.exists { dir =>
              val prefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"
              Files.exists(dir.resolve(s"${prefix}_haplogroup_report.txt"))
            }

            if (hasDetailedReport) {
              new HaplogroupReportDialog(
                treeType = treeType,
                artifactDir = artifactDir,
                sampleName = Some(subject.donorIdentifier)
              ).showAndWait()
            } else {
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

  /** Handles callable loci analysis for a specific alignment */
  private def handleCallableLociAnalysis(runIndex: Int, alignmentIndex: Int, alignment: Alignment): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      s"Callable Loci Analysis (${alignment.referenceBuild})",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runCallableLociAnalysisForAlignment(
      subject.sampleAccession,
      runIndex,
      alignmentIndex,
      onComplete = {
        case Right((result, artifactDir)) =>
          Platform.runLater {
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

  /** Handles WGS metrics analysis for a specific alignment */
  private def handleWgsMetricsAnalysis(runIndex: Int, alignmentIndex: Int, alignment: Alignment): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      s"WGS Metrics Analysis (${alignment.referenceBuild})",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runWgsMetricsAnalysisForAlignment(
      subject.sampleAccession,
      runIndex,
      alignmentIndex,
      onComplete = {
        case Right(wgsMetrics) =>
          Platform.runLater {
            new Alert(AlertType.Information) {
              title = "WGS Metrics Complete"
              headerText = s"Coverage Analysis (${alignment.referenceBuild})"
              contentText = s"""Mean Coverage: ${f"${wgsMetrics.meanCoverage}%.1f"}x
Median Coverage: ${f"${wgsMetrics.medianCoverage}%.1f"}x
SD Coverage: ${f"${wgsMetrics.sdCoverage}%.1f"}

% Bases at 10x: ${f"${wgsMetrics.pct10x * 100}%.1f"}%
% Bases at 20x: ${f"${wgsMetrics.pct20x * 100}%.1f"}%
% Bases at 30x: ${f"${wgsMetrics.pct30x * 100}%.1f"}%"""
            }.showAndWait()

            // Rebuild tree to show updated metrics
            treeTable.root = buildTreeItems()
          }
        case Left(error) =>
          Platform.runLater {
            new Alert(AlertType.Error) {
              title = "WGS Metrics Failed"
              headerText = "Could not complete WGS metrics analysis"
              contentText = error
            }.showAndWait()
          }
      }
    )

    progressDialog.show()
  }

  /** Handles read/insert metrics for a specific alignment */
  private def handleMultipleMetricsAnalysis(runIndex: Int, alignmentIndex: Int, alignment: Alignment): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      s"Read/Insert Metrics (${alignment.referenceBuild})",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runMultipleMetricsAnalysisForAlignment(
      subject.sampleAccession,
      runIndex,
      alignmentIndex,
      onComplete = (result: Either[String, ReadMetrics]) => result match {
        case Right(metricsResult) =>
          Platform.runLater {
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
              headerText = s"Read & Insert Size Metrics (${alignment.referenceBuild})"
              contentText = summaryText
            }.showAndWait()

            // Rebuild tree to show updated data
            treeTable.root = buildTreeItems()
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

  /** Handles generating whole-genome VCF for an alignment */
  private def handleGenerateVcf(ar: AlignmentRow): Unit = {
    // Show confirmation dialog with time estimate
    val confirmDialog = new Alert(AlertType.Confirmation) {
      title = "Generate Whole-Genome VCF"
      headerText = s"Generate VCF for ${ar.alignment.referenceBuild} alignment?"
      contentText = s"""This process is CPU-intensive and may take several hours
depending on your hardware.

Reference: ${ar.alignment.referenceBuild}
Contigs: All chromosomes (1-22, X, Y, MT)

The VCF will be cached locally and used for:
• Haplogroup analysis (improved gap-aware resolution)
• Future ancestry analysis

You can continue using other features while it runs."""
    }

    confirmDialog.showAndWait() match {
      case Some(ButtonType.OK) =>
        val progressDialog = new AnalysisProgressDialog(
          s"Generating Whole-Genome VCF (${ar.alignment.referenceBuild})",
          viewModel.analysisProgress,
          viewModel.analysisProgressPercent,
          viewModel.analysisInProgress
        )

        viewModel.runWholeGenomeVariantCallingForAlignment(
          subject.sampleAccession,
          ar.runIndex,
          ar.alignmentIndex,
          onComplete = {
            case Right(vcfInfo) =>
              Platform.runLater {
                val sizeGb = vcfInfo.fileSizeBytes / (1024.0 * 1024.0 * 1024.0)
                new Alert(AlertType.Information) {
                  title = "VCF Generation Complete"
                  headerText = s"Whole-Genome VCF Generated (${ar.alignment.referenceBuild})"
                  contentText = s"""Variants: ${f"${vcfInfo.variantCount}%,d"}
Size: ${f"$sizeGb%.2f"} GB
Contigs: ${vcfInfo.contigs.size}
Inferred Sex: ${vcfInfo.inferredSex.getOrElse("Unknown")}"""
                }.showAndWait()

                // Rebuild tree to show updated status
                treeTable.root = buildTreeItems()
              }
            case Left(error) =>
              Platform.runLater {
                new Alert(AlertType.Error) {
                  title = "VCF Generation Failed"
                  headerText = "Could not generate whole-genome VCF"
                  contentText = error
                }.showAndWait()
              }
          }
        )

        progressDialog.show()
      case _ => // User cancelled
    }
  }

  /** Shows VCF statistics for an alignment */
  private def handleViewVcfStats(ar: AlignmentRow): Unit = {
    getVcfStatus(ar) match {
      case VcfStatus.Available(info) =>
        val sizeGb = info.fileSizeBytes / (1024.0 * 1024.0 * 1024.0)
        new Alert(AlertType.Information) {
          title = "VCF Statistics"
          headerText = s"Whole-Genome VCF (${info.referenceBuild})"
          contentText = s"""Path: ${info.vcfPath}
Variants: ${f"${info.variantCount}%,d"}
Size: ${f"$sizeGb%.2f"} GB
Contigs: ${info.contigs.mkString(", ")}
Inferred Sex: ${info.inferredSex.getOrElse("Unknown")}
Created: ${info.createdAt}
GATK Version: ${info.gatkVersion}
Caller: ${info.callerVersion}"""
        }.showAndWait()

      case _ =>
        new Alert(AlertType.Information) {
          title = "No VCF Available"
          headerText = "No whole-genome VCF found"
          contentText = "Generate a whole-genome VCF first to view statistics."
        }.showAndWait()
    }
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

  private val buttonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(addButton)
  }

  children = Seq(
    new scalafx.scene.control.Label("Sequencing Runs:") { style = "-fx-font-weight: bold;" },
    treeTable,
    buttonBar
  )

  VBox.setVgrow(treeTable, Priority.Always)
}
