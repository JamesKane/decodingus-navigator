package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.i18n.Formatters
import com.decodingus.workspace.model.SequenceRun
import scalafx.Includes.*
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, Priority, Region, VBox}

/**
 * Result of merge dialog - which run to keep as primary.
 */
case class MergeDecision(
  primaryRun: SequenceRun,
  primaryIndex: Int,
  secondaryRun: SequenceRun,
  secondaryIndex: Int
)

/**
 * Dialog for manually merging two sequence runs that represent the same source data
 * (e.g., same sequencing run aligned to different reference genomes).
 *
 * The user selects which run to keep as the "primary" - the other run's alignments
 * will be moved to the primary run, and the secondary run will be deleted.
 */
class MergeSequenceRunsDialog(
  run1: SequenceRun,
  run1Index: Int,
  run2: SequenceRun,
  run2Index: Int
) extends Dialog[Option[MergeDecision]] {

  title = "Merge Sequence Runs"
  headerText = "Select which run to keep as the primary"

  val mergeButtonType = new ButtonType("Merge", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(mergeButtonType, ButtonType.Cancel)

  // Radio buttons to select primary
  private val selectionToggleGroup = new ToggleGroup()

  private val run1Radio = new RadioButton() {
    toggleGroup = selectionToggleGroup
    selected = true // Default to first run
  }

  private val run2Radio = new RadioButton() {
    toggleGroup = selectionToggleGroup
  }

  private def createRunSummary(run: SequenceRun): VBox = {
    val platform = run.platformName + run.instrumentModel.map(m => s" ($m)").getOrElse("")
    val testType = SequenceRun.testTypeDisplayName(run.testType)
    val reads = run.totalReads.map(r => Formatters.formatNumber(r) + " reads").getOrElse("Unknown reads")

    // Get alignment info
    val alignmentCount = run.alignmentRefs.size
    val alignmentInfo = if (alignmentCount > 0) {
      s"$alignmentCount alignment${if (alignmentCount > 1) "s" else ""}"
    } else {
      "No alignments"
    }

    // File info
    val fileCount = run.files.size
    val fileInfo = if (fileCount > 0) {
      s"$fileCount file${if (fileCount > 1) "s" else ""}"
    } else {
      "No files"
    }

    new VBox(5) {
      padding = Insets(10)
      style = "-fx-background-color: #f5f5f5; -fx-background-radius: 5; -fx-border-color: #cccccc; -fx-border-radius: 5;"
      children = Seq(
        new Label(platform) {
          style = "-fx-font-weight: bold; -fx-font-size: 13px;"
        },
        new Label(testType) {
          style = "-fx-text-fill: #666666;"
        },
        new HBox(15) {
          children = Seq(
            new Label(reads) { style = "-fx-text-fill: #888888; -fx-font-size: 11px;" },
            new Label(alignmentInfo) { style = "-fx-text-fill: #888888; -fx-font-size: 11px;" },
            new Label(fileInfo) { style = "-fx-text-fill: #888888; -fx-font-size: 11px;" }
          )
        },
        run.sampleName.map(sm => new Label(s"Sample: $sm") { style = "-fx-text-fill: #888888; -fx-font-size: 11px;" }).getOrElse(new Region),
        run.libraryId.map(lb => new Label(s"Library: $lb") { style = "-fx-text-fill: #888888; -fx-font-size: 11px;" }).getOrElse(new Region)
      )
    }
  }

  private val run1Panel = new HBox(10) {
    alignment = Pos.TopLeft
    children = Seq(run1Radio, createRunSummary(run1))
    hgrow = Priority.Always
  }

  private val run2Panel = new HBox(10) {
    alignment = Pos.TopLeft
    children = Seq(run2Radio, createRunSummary(run2))
    hgrow = Priority.Always
  }

  private val explanationLabel = new Label(
    """These two sequence runs appear to be from the same source data.
      |
      |Select which run to keep as the "primary". The other run's alignments
      |will be moved to the primary run, and the secondary run will be deleted.
      |
      |This is useful when the same reads have been aligned to different
      |reference genomes (e.g., GRCh38 and CHM13).""".stripMargin
  ) {
    wrapText = true
    prefWidth = 400
    style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
  }

  private val warningLabel = new Label(
    "Note: This operation cannot be undone."
  ) {
    style = "-fx-text-fill: #d97706; -fx-font-weight: bold; -fx-font-size: 11px;"
  }

  dialogPane().content = new VBox(15) {
    padding = Insets(20)
    children = Seq(
      explanationLabel,
      new Label("Run 1:") { style = "-fx-font-weight: bold;" },
      run1Panel,
      new Label("Run 2:") { style = "-fx-font-weight: bold;" },
      run2Panel,
      warningLabel
    )
  }

  dialogPane().setPrefWidth(480)

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == mergeButtonType) {
      val (primary, primaryIdx, secondary, secondaryIdx) = if (run1Radio.selected.value) {
        (run1, run1Index, run2, run2Index)
      } else {
        (run2, run2Index, run1, run1Index)
      }
      Some(MergeDecision(primary, primaryIdx, secondary, secondaryIdx))
    } else {
      None
    }
  }
}
