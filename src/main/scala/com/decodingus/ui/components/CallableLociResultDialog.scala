package com.decodingus.ui.components

import com.decodingus.analysis.{CallableLociResult, ContigCoverageStats}
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Orientation}
import scalafx.scene.control.*
import scalafx.scene.layout.{HBox, Priority, StackPane, VBox}
import scalafx.scene.web.WebView

import java.io.File
import java.nio.file.{Files, Path}
import scala.io.Source
import scala.util.Using

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
  dialogPane().setPrefSize(1000, 750)

  // Format large numbers with commas
  private def formatNumber(n: Long): String = f"$n%,d"

  // Calculate total bases and percentages
  private val totalBases = result.contigAnalysis.map { cs =>
    cs.callable + cs.noCoverage + cs.lowCoverage + cs.excessiveCoverage + cs.poorMappingQuality
  }.sum

  private val callablePercent = if (totalBases > 0) {
    (result.callableBases.toDouble / totalBases * 100)
  } else 0.0

  // Table data model - enriched with samtools coverage metrics
  case class ContigRow(
                        contig: String,
                        callable: Long,
                        noCoverage: Long,
                        lowCoverage: Long,
                        excessiveCoverage: Long,
                        poorMappingQuality: Long,
                        callablePercent: Double,
                        // samtools coverage-style metrics (optional)
                        numReads: Option[Long] = None,
                        covBases: Option[Long] = None,
                        coveragePct: Option[Double] = None,
                        meanDepth: Option[Double] = None,
                        meanBaseQ: Option[Double] = None,
                        meanMapQ: Option[Double] = None
                      )

  // Load samtools-style coverage stats from coverage.txt if available
  private val coverageStatsMap: Map[String, ContigCoverageStats] = artifactDir.map { dir =>
    val coverageFile = dir.resolve("callable_loci").resolve("coverage.txt")
    if (Files.exists(coverageFile)) {
      loadCoverageStats(coverageFile)
    } else {
      Map.empty[String, ContigCoverageStats]
    }
  }.getOrElse(Map.empty)

  private def loadCoverageStats(coverageFile: Path): Map[String, ContigCoverageStats] = {
    val result = scala.collection.mutable.Map[String, ContigCoverageStats]()
    Using(Source.fromFile(coverageFile.toFile)) { source =>
      for (line <- source.getLines()) {
        if (!line.startsWith("#") && line.trim.nonEmpty) {
          val fields = line.split("\\t")
          if (fields.length >= 9) {
            val stats = ContigCoverageStats(
              contig = fields(0),
              startPos = fields(1).toLong,
              endPos = fields(2).toLong,
              numReads = fields(3).toLong,
              covBases = fields(4).toLong,
              coverage = fields(5).toDouble,
              meanDepth = fields(6).toDouble,
              meanBaseQ = fields(7).toDouble,
              meanMapQ = fields(8).toDouble
            )
            result(stats.contig) = stats
          }
        }
      }
    }
    result.toMap
  }

  private val tableData = ObservableBuffer.from(
    result.contigAnalysis.map { cs =>
      val total = cs.callable + cs.noCoverage + cs.lowCoverage + cs.excessiveCoverage + cs.poorMappingQuality
      val pct = if (total > 0) cs.callable.toDouble / total * 100 else 0.0
      val covStats = coverageStatsMap.get(cs.contigName)
      ContigRow(
        cs.contigName,
        cs.callable,
        cs.noCoverage,
        cs.lowCoverage,
        cs.excessiveCoverage,
        cs.poorMappingQuality,
        pct,
        numReads = covStats.map(_.numReads),
        covBases = covStats.map(_.covBases),
        coveragePct = covStats.map(_.coverage),
        meanDepth = covStats.map(_.meanDepth),
        meanBaseQ = covStats.map(_.meanBaseQ),
        meanMapQ = covStats.map(_.meanMapQ)
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
          val html =
            s"""
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

  // Check if we have coverage stats available (determines which columns to show)
  private val hasCoverageStats = coverageStatsMap.nonEmpty

  private val contigTable = new TableView[ContigRow](tableData) {
    prefHeight = 280
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Core callable loci columns
    val coreColumns = Seq(
      new TableColumn[ContigRow, String] {
        text = "Contig"
        cellValueFactory = r => StringProperty(r.value.contig)
        prefWidth = 60
      },
      new TableColumn[ContigRow, String] {
        text = "Callable"
        cellValueFactory = r => StringProperty(formatNumber(r.value.callable))
        prefWidth = 90
      },
      new TableColumn[ContigRow, String] {
        text = "% Callable"
        cellValueFactory = r => StringProperty(f"${r.value.callablePercent}%.1f%%")
        prefWidth = 70
      }
    )

    // Coverage metrics columns (only shown if coverage.txt exists)
    val coverageColumns = if (hasCoverageStats) Seq(
      new TableColumn[ContigRow, String] {
        text = "Reads"
        cellValueFactory = r => StringProperty(r.value.numReads.map(formatNumber).getOrElse("-"))
        prefWidth = 85
      },
      new TableColumn[ContigRow, String] {
        text = "Cov Bases"
        cellValueFactory = r => StringProperty(r.value.covBases.map(formatNumber).getOrElse("-"))
        prefWidth = 90
      },
      new TableColumn[ContigRow, String] {
        text = "% Cov"
        cellValueFactory = r => StringProperty(r.value.coveragePct.map(p => f"$p%.1f%%").getOrElse("-"))
        prefWidth = 60
      },
      new TableColumn[ContigRow, String] {
        text = "Depth"
        cellValueFactory = r => StringProperty(r.value.meanDepth.map(d => f"$d%.1f").getOrElse("-"))
        prefWidth = 55
      },
      new TableColumn[ContigRow, String] {
        text = "BaseQ"
        cellValueFactory = r => StringProperty(r.value.meanBaseQ.map(q => f"$q%.1f").getOrElse("-"))
        prefWidth = 50
      },
      new TableColumn[ContigRow, String] {
        text = "MapQ"
        cellValueFactory = r => StringProperty(r.value.meanMapQ.map(q => f"$q%.1f").getOrElse("-"))
        prefWidth = 50
      }
    ) else Seq.empty

    // Additional callable loci detail columns
    val detailColumns = Seq(
      new TableColumn[ContigRow, String] {
        text = "No Cov"
        cellValueFactory = r => StringProperty(formatNumber(r.value.noCoverage))
        prefWidth = 80
      },
      new TableColumn[ContigRow, String] {
        text = "Low Cov"
        cellValueFactory = r => StringProperty(formatNumber(r.value.lowCoverage))
        prefWidth = 80
      },
      new TableColumn[ContigRow, String] {
        text = "Excess"
        cellValueFactory = r => StringProperty(formatNumber(r.value.excessiveCoverage))
        prefWidth = 70
      },
      new TableColumn[ContigRow, String] {
        text = "Poor MapQ"
        cellValueFactory = r => StringProperty(formatNumber(r.value.poorMappingQuality))
        prefWidth = 80
      }
    )

    columns ++= (coreColumns ++ coverageColumns ++ detailColumns).map(_.delegate)

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
      new Label("Coverage Histogram:") {
        style = "-fx-font-weight: bold;"
      },
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
      new Label("") {
        visible = false
      }
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
