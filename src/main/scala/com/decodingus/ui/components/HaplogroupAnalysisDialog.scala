package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, RadioButton, ToggleGroup, ButtonBar}
import scalafx.scene.layout.{VBox, HBox}
import scalafx.geometry.Insets
import com.decodingus.haplogroup.tree.TreeType

/**
 * Dialog for selecting haplogroup analysis type (Y-DNA or MT-DNA).
 */
class HaplogroupAnalysisDialog extends Dialog[Option[TreeType]] {
  title = "Haplogroup Analysis"
  headerText = "Select haplogroup type to analyze"

  val analyzeButtonType = new ButtonType("Analyze", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(analyzeButtonType, ButtonType.Cancel)

  private val treeTypeToggleGroup = new ToggleGroup()

  private val yDnaRadio = new RadioButton("Y-DNA (Paternal lineage)") {
    toggleGroup = treeTypeToggleGroup
    selected = true
  }

  private val mtDnaRadio = new RadioButton("MT-DNA (Maternal lineage)") {
    toggleGroup = treeTypeToggleGroup
  }

  private val content = new VBox(15) {
    padding = Insets(20)
    children = Seq(
      new Label("Choose which haplogroup lineage to analyze:") {
        style = "-fx-font-weight: bold;"
      },
      new VBox(10) {
        children = Seq(
          yDnaRadio,
          new Label("    Analyzes the Y chromosome for paternal ancestry") {
            style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
          },
          mtDnaRadio,
          new Label("    Analyzes mitochondrial DNA for maternal ancestry") {
            style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
          }
        )
      },
      new Label("") { prefHeight = 10 }, // Spacer
      new Label("Note: Analysis requires initial analysis to be completed first.") {
        style = "-fx-text-fill: #888888; -fx-font-size: 11px; -fx-font-style: italic;"
        wrapText = true
      }
    )
  }

  dialogPane().content = content

  resultConverter = dialogButton => {
    if (dialogButton == analyzeButtonType) {
      if (yDnaRadio.selected.value) Some(TreeType.YDNA)
      else Some(TreeType.MTDNA)
    } else {
      None
    }
  }
}

/**
 * Dialog showing haplogroup analysis results.
 */
class HaplogroupResultDialog(
  treeType: TreeType,
  haplogroupName: String,
  score: Double,
  matchingSnps: Int,
  mismatchingSnps: Int,
  ancestralMatches: Int,
  depth: Int
) extends Dialog[Unit] {

  title = "Haplogroup Analysis Results"
  headerText = s"${if (treeType == TreeType.YDNA) "Y-DNA" else "MT-DNA"} Haplogroup: $haplogroupName"

  dialogPane().buttonTypes = Seq(ButtonType.OK)

  private val content = new VBox(10) {
    padding = Insets(20)
    children = Seq(
      new Label(s"Top Haplogroup: $haplogroupName") {
        style = "-fx-font-size: 18px; -fx-font-weight: bold;"
      },
      new Label("") { prefHeight = 5 },
      new Label("Analysis Details:") { style = "-fx-font-weight: bold;" },
      createStatRow("Score:", f"$score%.4f"),
      createStatRow("Matching SNPs:", matchingSnps.toString),
      createStatRow("Mismatching SNPs:", mismatchingSnps.toString),
      createStatRow("Ancestral Matches:", ancestralMatches.toString),
      createStatRow("Tree Depth:", depth.toString),
      new Label("") { prefHeight = 10 },
      new Label("This result has been saved to the subject's profile.") {
        style = "-fx-text-fill: #4CAF50; -fx-font-style: italic;"
      }
    )
  }

  private def createStatRow(label: String, value: String): HBox = {
    new HBox(10) {
      children = Seq(
        new Label(label) { prefWidth = 150 },
        new Label(value) { style = "-fx-font-weight: bold;" }
      )
    }
  }

  dialogPane().content = content
}
