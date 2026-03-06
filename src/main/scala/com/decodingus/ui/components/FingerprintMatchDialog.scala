package com.decodingus.ui.components

import com.decodingus.workspace.model.SequenceRun
import scalafx.Includes.*
import scalafx.geometry.Insets
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, VBox}

/**
 * Result of the fingerprint match confirmation dialog.
 */
sealed trait FingerprintMatchDecision

case object GroupTogether extends FingerprintMatchDecision

case object KeepSeparate extends FingerprintMatchDecision

/**
 * Dialog for confirming fingerprint matches when confidence is LOW.
 * Shown when the system detects a potential match based on read statistics
 * but cannot confirm with @RG header fields.
 */
class FingerprintMatchDialog(
                              existingRun: SequenceRun,
                              newReferenceBuild: String,
                              matchConfidence: String,
                              totalReads: Long,
                              sampleName: String,
                              libraryId: String
                            ) extends Dialog[Option[FingerprintMatchDecision]] {

  title = "Potential Matching Sequence Run Detected"
  headerText = "This file may be from the same sequencing run"

  val groupButtonType = new ButtonType("Group Together", ButtonBar.ButtonData.OKDone)
  val separateButtonType = new ButtonType("Keep Separate", ButtonBar.ButtonData.CancelClose)
  dialogPane().buttonTypes = Seq(groupButtonType, separateButtonType, ButtonType.Cancel)

  private val content = new VBox(15) {
    padding = Insets(20)
    children = Seq(
      new Label("The new file appears to be a different reference alignment of an existing sequencing run.") {
        wrapText = true
        prefWidth = 450
      },
      new Label(s"Match Confidence: $matchConfidence") {
        style = matchConfidence match {
          case "HIGH" => "-fx-font-weight: bold; -fx-text-fill: #4CAF50;"
          case "MEDIUM" => "-fx-font-weight: bold; -fx-text-fill: #FF9800;"
          case _ => "-fx-font-weight: bold; -fx-text-fill: #F44336;"
        }
      },
      new Label("") {
        prefHeight = 5
      },
      new Label("Comparison:") {
        style = "-fx-font-weight: bold;"
      },
      createComparisonGrid(),
      new Label("") {
        prefHeight = 5
      },
      new Label("What would you like to do?") {
        style = "-fx-font-weight: bold;"
      },
      new VBox(8) {
        children = Seq(
          new Label("• Group Together: Add this as another reference alignment of the same run") {
            style = "-fx-text-fill: #666666;"
          },
          new Label("• Keep Separate: Create a new sequence run entry") {
            style = "-fx-text-fill: #666666;"
          }
        )
      }
    )
  }

  private def createComparisonGrid(): GridPane = {
    new GridPane {
      hgap = 15
      vgap = 5
      padding = Insets(10, 0, 10, 20)

      // Header row
      add(new Label("") {
        prefWidth = 100
      }, 0, 0)
      add(new Label("Existing Run") {
        style = "-fx-font-weight: bold;"
      }, 1, 0)
      add(new Label("New File") {
        style = "-fx-font-weight: bold;"
      }, 2, 0)

      // Reference build
      add(new Label("Reference:"), 0, 1)
      add(new Label(existingRun.alignmentRefs.headOption
        .flatMap(ref => Option(ref).map(_.split(":").lift(2).getOrElse("Unknown")))
        .getOrElse(existingRun.files.headOption.map(_ => "Analyzed").getOrElse("—"))), 1, 1)
      add(new Label(newReferenceBuild), 2, 1)

      // Sample name
      add(new Label("Sample:"), 0, 2)
      add(new Label(existingRun.sampleName.getOrElse("—")), 1, 2)
      add(new Label(if (sampleName != "Unknown") sampleName else "—"), 2, 2)

      // Library ID
      add(new Label("Library ID:"), 0, 3)
      add(new Label(existingRun.libraryId.getOrElse("—")), 1, 3)
      add(new Label(if (libraryId != "Unknown") libraryId else "—"), 2, 3)

      // Total reads
      add(new Label("Total Reads:"), 0, 4)
      add(new Label(existingRun.totalReads.map(formatReads).getOrElse("—")), 1, 4)
      add(new Label(formatReads(totalReads)), 2, 4)

      // Platform
      add(new Label("Platform:"), 0, 5)
      add(new Label(existingRun.platformName), 1, 5)
      add(new Label(existingRun.platformName), 2, 5) // Same platform assumed
    }
  }

  private def formatReads(count: Long): String = {
    if (count >= 1_000_000_000) f"${count / 1_000_000_000.0}%.1fB"
    else if (count >= 1_000_000) f"${count / 1_000_000.0}%.1fM"
    else if (count >= 1_000) f"${count / 1_000.0}%.1fK"
    else count.toString
  }

  dialogPane().content = content

  resultConverter = dialogButton => {
    if (dialogButton == groupButtonType) Some(GroupTogether)
    else if (dialogButton == separateButtonType) Some(KeepSeparate)
    else None
  }
}
