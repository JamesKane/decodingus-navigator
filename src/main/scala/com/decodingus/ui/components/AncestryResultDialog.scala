package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{Dialog, ButtonType, Label, TableView, TableColumn, ScrollPane, SplitPane, ProgressBar}
import scalafx.scene.layout.{VBox, HBox, Priority, GridPane, StackPane, Region}
import scalafx.scene.chart.PieChart
import scalafx.geometry.{Insets, Orientation, Pos}
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import com.decodingus.ancestry.model.{AncestryResult, Population, PopulationPercentage, SuperPopulationPercentage}

/**
 * Dialog showing ancestry analysis results.
 * Left side: pie chart of continental (super-population) breakdown
 * Right side: detailed table of sub-population percentages with confidence intervals
 */
class AncestryResultDialog(result: AncestryResult) extends Dialog[Unit] {

  title = "Ancestry Analysis Results"
  headerText = s"Population Breakdown (${result.panelType.capitalize} Panel)"

  dialogPane().buttonTypes = Seq(ButtonType.OK)
  dialogPane().setPrefSize(900, 650)

  // Super-population colors (continental groups)
  private val superPopColors = Map(
    "European" -> "#3366CC",
    "African" -> "#FF6600",
    "East Asian" -> "#00AA00",
    "South Asian" -> "#9933FF",
    "Americas" -> "#CC0066",
    "West Asian" -> "#996633",
    "Oceanian" -> "#00AAAA",
    "Central Asian" -> "#66CCCC",
    "Native American" -> "#990033"
  )

  // Create pie chart data from super-population summary
  private val pieChartData = ObservableBuffer.from(
    result.superPopulationSummary
      .filter(_.percentage >= 0.5) // Hide trace amounts
      .map { sp =>
        val slice = PieChart.Data(f"${sp.superPopulation} (${sp.percentage}%.1f%%)", sp.percentage)
        // Note: Color styling is applied after chart is rendered
        slice
      }
  )

  private val pieChart = new PieChart(pieChartData) {
    title = "Continental Breakdown"
    legendVisible = true
    labelsVisible = true
    prefWidth = 350
    prefHeight = 350
    style = "-fx-font-size: 11px;"
  }

  // Apply colors to pie slices after data is loaded
  pieChartData.zipWithIndex.foreach { case (data, idx) =>
    val superPop = result.superPopulationSummary.filter(_.percentage >= 0.5)(idx).superPopulation
    val color = superPopColors.getOrElse(superPop, "#888888")
    data.node.value.style = s"-fx-pie-color: $color;"
  }

  // Table data for detailed population breakdown
  case class PopulationRow(
    code: String,
    name: String,
    superPop: String,
    percentage: Double,
    ciLow: Double,
    ciHigh: Double,
    rank: Int
  )

  private val tableData = ObservableBuffer.from(
    result.percentages
      .filter(_.percentage >= 0.1) // Hide trace amounts in table
      .map { p =>
        val pop = Population.byCode(p.populationCode)
        PopulationRow(
          p.populationCode,
          p.populationName,
          pop.map(_.superPopulation).getOrElse("Unknown"),
          p.percentage,
          p.confidenceLow,
          p.confidenceHigh,
          p.rank
        )
      }
  )

  private val populationTable = new TableView[PopulationRow](tableData) {
    prefHeight = 400
    columnResizePolicy = TableView.ConstrainedResizePolicy

    columns ++= Seq(
      new TableColumn[PopulationRow, String] {
        text = "#"
        cellValueFactory = r => StringProperty(r.value.rank.toString)
        prefWidth = 35
      },
      new TableColumn[PopulationRow, String] {
        text = "Population"
        cellValueFactory = r => StringProperty(r.value.name)
        prefWidth = 150
      },
      new TableColumn[PopulationRow, String] {
        text = "Region"
        cellValueFactory = r => StringProperty(r.value.superPop)
        prefWidth = 90
      },
      new TableColumn[PopulationRow, String] {
        text = "Percentage"
        cellValueFactory = r => StringProperty(f"${r.value.percentage}%.1f%%")
        prefWidth = 80
      },
      new TableColumn[PopulationRow, String] {
        text = "95% CI"
        cellValueFactory = r => StringProperty(f"${r.value.ciLow}%.1f%% - ${r.value.ciHigh}%.1f%%")
        prefWidth = 100
      }
    )
  }
  VBox.setVgrow(populationTable, Priority.Always)

  // Quality metrics section
  private val qualitySection = new VBox(5) {
    padding = Insets(10, 0, 0, 0)
    children = Seq(
      new Label("Data Quality") { style = "-fx-font-weight: bold; -fx-font-size: 13px;" },
      new GridPane {
        hgap = 15
        vgap = 5
        add(new Label("SNPs Analyzed:"), 0, 0)
        add(new Label(f"${result.snpsAnalyzed}%,d"), 1, 0)
        add(new Label("SNPs with Data:"), 0, 1)
        add(new Label(f"${result.snpsWithGenotype}%,d (${result.snpsWithGenotype.toDouble / result.snpsAnalyzed * 100}%.1f%%)"), 1, 1)
        add(new Label("Missing SNPs:"), 0, 2)
        add(new Label(f"${result.snpsMissing}%,d"), 1, 2)
        add(new Label("Confidence:"), 0, 3)
        add(new HBox(5) {
          alignment = Pos.CenterLeft
          children = Seq(
            new ProgressBar {
              progress = result.confidenceLevel
              prefWidth = 100
              prefHeight = 16
            },
            new Label(f"${result.confidenceLevel * 100}%.0f%%")
          )
        }, 1, 3)
      }
    )
  }

  // Left panel: pie chart + quality metrics
  private val leftPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(pieChart, qualitySection)
    prefWidth = 380
  }

  // Right panel: detailed table
  private val rightPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new Label("Detailed Population Breakdown") { style = "-fx-font-weight: bold; -fx-font-size: 13px;" },
      populationTable
    )
  }
  HBox.setHgrow(rightPanel, Priority.Always)

  // Main layout: horizontal split
  private val mainContent = new HBox(0) {
    children = Seq(leftPanel, rightPanel)
  }

  // Version info at bottom
  private val versionInfo = new Label(s"Analysis: ${result.analysisVersion} | Reference: ${result.referenceVersion}") {
    style = "-fx-font-size: 10px; -fx-text-fill: #888;"
  }

  private val content = new VBox(5) {
    padding = Insets(5)
    children = Seq(mainContent, versionInfo)
  }
  VBox.setVgrow(content, Priority.Always)

  dialogPane().content = content

  // Make dialog resizable
  dialogPane().getScene.getWindow match {
    case stage: javafx.stage.Stage => stage.setResizable(true)
    case _ =>
  }
}

object AncestryResultDialog {
  /**
   * Create a summary string for the ancestry result.
   */
  def summaryString(result: AncestryResult): String = {
    val top = result.superPopulationSummary
      .filter(_.percentage >= 1.0)
      .take(3)
      .map(p => f"${p.superPopulation}: ${p.percentage}%.1f%%")
      .mkString(", ")
    if (top.nonEmpty) top else "Analysis complete"
  }
}
