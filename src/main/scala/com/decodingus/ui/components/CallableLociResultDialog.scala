package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{Dialog, ButtonType, Label, TableView, TableColumn, ScrollPane, SplitPane}
import scalafx.scene.layout.{VBox, HBox, Priority, StackPane}
import scalafx.scene.web.WebView
import scalafx.geometry.{Insets, Orientation}
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import com.decodingus.analysis.CallableLociResult

import java.io.File
import java.nio.file.{Files, Path}

/**
 * Dialog showing callable loci analysis results.
 * Upper half: numeric summary table
 * Lower half: SVG histogram for selected contig
 */
class CallableLociResultDialog(
  result: CallableLociResult,
  artifactDir: Option[Path] = None
) extends Dialog[Unit] {

  title = "Callable Loci Analysis Results"
  headerText = "Genome Coverage Analysis Complete"

  dialogPane().buttonTypes = Seq(ButtonType.OK)
  dialogPane().setPrefSize(800, 700)

  // Format large numbers with commas
  private def formatNumber(n: Long): String = f"$n%,d"

  // Calculate total bases and percentages
  private val totalBases = result.contigAnalysis.map { cs =>
    cs.callable + cs.noCoverage + cs.lowCoverage + cs.excessiveCoverage + cs.poorMappingQuality
  }.sum

  private val callablePercent = if (totalBases > 0) {
    (result.callableBases.toDouble / totalBases * 100)
  } else 0.0

  // Table data model
  case class ContigRow(
    contig: String,
    callable: Long,
    noCoverage: Long,
    lowCoverage: Long,
    excessiveCoverage: Long,
    poorMappingQuality: Long,
    callablePercent: Double
  )

  private val tableData = ObservableBuffer.from(
    result.contigAnalysis.map { cs =>
      val total = cs.callable + cs.noCoverage + cs.lowCoverage + cs.excessiveCoverage + cs.poorMappingQuality
      val pct = if (total > 0) cs.callable.toDouble / total * 100 else 0.0
      ContigRow(
        cs.contigName,
        cs.callable,
        cs.noCoverage,
        cs.lowCoverage,
        cs.excessiveCoverage,
        cs.poorMappingQuality,
        pct
      )
    }
  )

  // WebView for displaying SVG
  private val svgWebView = new WebView {
    prefHeight = 280
  }

  private def loadSvgForContig(contigName: String): Unit = {
    artifactDir match {
      case Some(dir) =>
        val svgFile = dir.resolve("callable_loci").resolve(s"$contigName.callable.svg")
        if (Files.exists(svgFile)) {
          val svgContent = Files.readString(svgFile)
          // Wrap SVG in HTML for proper rendering
          val html = s"""
            |<!DOCTYPE html>
            |<html>
            |<head>
            |  <style>
            |    body { margin: 0; padding: 10px; background: #222; display: flex; justify-content: center; }
            |    svg { max-width: 100%; height: auto; }
            |  </style>
            |</head>
            |<body>
            |$svgContent
            |</body>
            |</html>
          """.stripMargin
          svgWebView.engine.loadContent(html)
        } else {
          svgWebView.engine.loadContent(s"<html><body style='background:#333;color:#ccc;padding:20px;'>SVG not found: $svgFile</body></html>")
        }
      case None =>
        svgWebView.engine.loadContent("<html><body style='background:#333;color:#ccc;padding:20px;'>No artifact directory provided</body></html>")
    }
  }

  private val contigTable = new TableView[ContigRow](tableData) {
    prefHeight = 250
    columnResizePolicy = TableView.ConstrainedResizePolicy

    columns ++= Seq(
      new TableColumn[ContigRow, String] {
        text = "Contig"
        cellValueFactory = r => StringProperty(r.value.contig)
        prefWidth = 70
      },
      new TableColumn[ContigRow, String] {
        text = "Callable"
        cellValueFactory = r => StringProperty(formatNumber(r.value.callable))
        prefWidth = 95
      },
      new TableColumn[ContigRow, String] {
        text = "No Cov"
        cellValueFactory = r => StringProperty(formatNumber(r.value.noCoverage))
        prefWidth = 85
      },
      new TableColumn[ContigRow, String] {
        text = "Low Cov"
        cellValueFactory = r => StringProperty(formatNumber(r.value.lowCoverage))
        prefWidth = 85
      },
      new TableColumn[ContigRow, String] {
        text = "Excess"
        cellValueFactory = r => StringProperty(formatNumber(r.value.excessiveCoverage))
        prefWidth = 80
      },
      new TableColumn[ContigRow, String] {
        text = "Poor MapQ"
        cellValueFactory = r => StringProperty(formatNumber(r.value.poorMappingQuality))
        prefWidth = 85
      },
      new TableColumn[ContigRow, String] {
        text = "% Callable"
        cellValueFactory = r => StringProperty(f"${r.value.callablePercent}%.1f%%")
        prefWidth = 75
      }
    )

    // When a row is selected, load its SVG
    selectionModel().selectedItem.onChange { (_, _, selected) =>
      if (selected != null) {
        loadSvgForContig(selected.contig)
      }
    }
  }
  VBox.setVgrow(contigTable, Priority.Always)

  // Summary header
  private val summarySection = new VBox(5) {
    padding = Insets(0, 0, 10, 0)
    children = Seq(
      new HBox(20) {
        children = Seq(
          new Label(s"Total Callable: ${formatNumber(result.callableBases)}") {
            style = "-fx-font-size: 16px; -fx-font-weight: bold;"
          },
          new Label(f"($callablePercent%.1f%% of genome)") {
            style = "-fx-font-size: 14px;"
          }
        )
      }
    )
  }

  // Upper section: summary + table
  private val upperSection = new VBox(5) {
    padding = Insets(10)
    children = Seq(
      summarySection,
      new Label("Per-Contig Analysis (select row to view histogram):") {
        style = "-fx-font-weight: bold;"
      },
      contigTable
    )
  }
  VBox.setVgrow(upperSection, Priority.Always)

  // Lower section: SVG histogram
  private val lowerSection = new VBox(5) {
    padding = Insets(10)
    children = Seq(
      new Label("Coverage Histogram:") { style = "-fx-font-weight: bold;" },
      svgWebView
    )
  }
  VBox.setVgrow(lowerSection, Priority.Always)

  // Split pane with upper (table) and lower (histogram) sections
  private val splitPane = new SplitPane {
    orientation = Orientation.Vertical
    items.addAll(upperSection, lowerSection)
    dividerPositions = 0.55
  }
  VBox.setVgrow(splitPane, Priority.Always)

  private val artifactPathLabel = artifactDir match {
    case Some(dir) =>
      new Label(s"Artifacts: ${dir.resolve("callable_loci")}") {
        style = "-fx-font-size: 11px; -fx-text-fill: #888;"
      }
    case None =>
      new Label("") { visible = false }
  }

  private val content = new VBox(5) {
    padding = Insets(10)
    children = Seq(splitPane, artifactPathLabel)
  }
  VBox.setVgrow(content, Priority.Always)

  dialogPane().content = content

  // Make dialog resizable
  dialogPane().getScene.getWindow match {
    case stage: javafx.stage.Stage => stage.setResizable(true)
    case _ =>
  }

  // Auto-select first contig to show its histogram
  if (tableData.nonEmpty) {
    contigTable.selectionModel().selectFirst()
  }
}
