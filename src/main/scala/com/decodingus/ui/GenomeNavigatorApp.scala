package com.decodingus.ui

import com.decodingus.analysis._
import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.tree.{TreeProviderType, TreeType}
import com.decodingus.model._
import com.decodingus.pds.PdsClient
import javafx.concurrent as jfxc
import scalafx.Includes.*
import scalafx.application.JFXApp3.PrimaryStage
import scalafx.application.{JFXApp3, Platform}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.Scene
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.*
import scalafx.scene.control.TableColumn.sfxTableColumn2jfx
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.{GridPane, HBox, StackPane, VBox}
import scalafx.scene.text.{Text, TextAlignment}
import scalafx.scene.web.WebView

import scala.concurrent.ExecutionContext.Implicits.global

case class ContigAnalysisRow(
  contig: String,
  callableBases: Long,
  noCoverage: Long,
  lowCoverage: Long,
  poorMq: Long,
  refN: Long,
  svgFile: String
)

object GenomeNavigatorApp extends JFXApp3 {
  private val mainLayout = new StackPane()
  private var currentFilePath: String = ""

  override def start(): Unit = {
    stage = new PrimaryStage {
      title = "Decoding-Us Navigator"
      scene = new Scene(800, 800) {
        root = mainLayout
        stylesheets.add(getClass.getResource("/style.css").toExternalForm)
      }
    }

    val welcomeScreen = createWelcomeScreen()
    mainLayout.children = welcomeScreen
  }

  private def createWelcomeScreen(): VBox = {
    new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      children = Seq(
        new Label("Welcome to Decoding-Us Navigator") {
          styleClass.add("title-label")
        },
        new Label("Drag and drop your BAM/CRAM file here, or click to select.") {
          styleClass.add("info-label")
          textAlignment = TextAlignment.Center
        },
        new StackPane {
          prefWidth = 400
          prefHeight = 200
          styleClass.add("drag-drop-area")
          children = new Label("Drop File Here") {
            styleClass.add("drag-drop-text")
          }

          onDragOver = (event: DragEvent) => {
            if (event.gestureSource != this && event.dragboard.hasFiles) {
              event.acceptTransferModes(TransferMode.Copy, TransferMode.Move)
            }
            event.consume()
          }

          onDragDropped = (event: DragEvent) => {
            val db = event.dragboard
            var success = false
            if (db.hasFiles) {
              db.files.headOption.foreach { file =>
                val filePath = file.getAbsolutePath
                println(s"Dropped file: $filePath")
                startAnalysis(filePath)
                success = true
              }
            }
            event.dropCompleted = success
            event.consume()
          }
        },
        new Button("Select BAM/CRAM File") {
          styleClass.add("button-select")
          onAction = _ => selectFile()
        }
      )
    }
  }

  private def selectFile(): Unit = {
    println("File selection dialog would open here.")
    startAnalysis("path/to/mock/file.bam")
  }

  private def startAnalysis(filePath: String): Unit = {
    currentFilePath = filePath
    val progressLabel = new Label("Analysis in progress...") {
      styleClass.add("progress-label")
    }
    val progressBar = new ProgressBar {
      prefWidth = 400
    }
    val progressIndicator = new ProgressIndicator

    val intermediateSummaryBox = new VBox(10) {
      alignment = Pos.Center
      visible = false
      padding = Insets(20)
      style = "-fx-border-color: #555; -fx-border-width: 1; -fx-background-color: #333;"
    }

    val progressScreen = new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      children = Seq(
        progressLabel,
        new HBox(20) {
          alignment = Pos.Center
          children = Seq(progressBar, progressIndicator)
        },
        intermediateSummaryBox
      )
    }

    mainLayout.children = progressScreen

    val jfxTask = new jfxc.Task[(CoverageSummary, List[String])]() {
      override def call(): (CoverageSummary, List[String]) = {
        try {
          val libraryStatsProcessor = new LibraryStatsProcessor()
          val wgsMetricsProcessor = new WgsMetricsProcessor()
          val callableLociProcessor = new CallableLociProcessor()
          val referencePath = "/Library/Genomics/Reference/chm13v2.0/chm13v2.0.fa.gz"

          // Phase 1: Library Stats (quick scan)
          val libraryStats = libraryStatsProcessor.process(filePath, (message, current, total) => {
            Platform.runLater { progressLabel.text = s"Library Stats: $message" }
            updateProgress(current, total * 3) // 0-33%
          })

          Platform.runLater {
            intermediateSummaryBox.children = Seq(
              new Label(s"Sample: ${libraryStats.sampleName}") { styleClass.add("info-label") },
              new Label(s"Reference: ${libraryStats.referenceBuild}") { styleClass.add("info-label") },
              new Label(s"Platform: ${libraryStats.inferredPlatform}") { styleClass.add("info-label") }
            )
            intermediateSummaryBox.visible = true
          }

          // Phase 2: WGS Metrics
          val wgsMetrics = wgsMetricsProcessor.process(filePath, referencePath, (message, current, total) => {
            Platform.runLater { progressLabel.text = message }
            updateProgress(total + current, total * 3) // 33-66%
          })

          // Phase 3: Callable Loci Analysis
          val (callableLociResult, svgStrings) = callableLociProcessor.process(filePath, referencePath, (message, current, total) => {
            Platform.runLater { progressLabel.text = s"Callable Loci: $message" }
            updateProgress(total * 2 + current, total * 3) // 66-100%
          })

          val coverageSummary = CoverageSummary(
            pdsUserId = "60820188481374", // placeholder
            libraryStats = libraryStats,
            wgsMetrics = wgsMetrics,
            callableBases = callableLociResult.callableBases,
            contigAnalysis = callableLociResult.contigAnalysis
          )
          (coverageSummary, svgStrings)
        } catch {
          case e: Exception =>
            e.printStackTrace()
            cancel()
            throw e
        }
      }
    }

    progressBar.progress <== jfxTask.progressProperty
    progressIndicator.progress <== jfxTask.progressProperty

    jfxTask.setOnSucceeded(_ => {
      val (summary, svgStrings) = jfxTask.getValue
      showResults(summary, svgStrings)
    })

    jfxTask.setOnFailed(_ => {
      val errorScreen = new VBox(20) {
        alignment = Pos.Center
        styleClass.add("root-pane")
        children = Seq(
          new Label("Analysis Failed!") {
            styleClass.add("error-label")
          },
          new Label("Please check the console for more details and ensure the reference and input files are correct.") {
            styleClass.add("info-label")
          },
          new Button("Back to Welcome") {
            onAction = _ => mainLayout.children = createWelcomeScreen()
          }
        )
      }
      mainLayout.children = errorScreen
    })

    new Thread(jfxTask).start()
  }

  private def showResults(summary: CoverageSummary, svgStrings: List[String]): Unit = {
    val resultsTitle = new Label("Analysis Results") {
      styleClass.add("title-label")
    }

    val statsGrid = new GridPane {
      hgap = 20
      vgap = 10
      padding = Insets(20)
      styleClass.add("stats-grid")
    }

    def addStat(name: String, value: String, row: Int): Unit = {
      statsGrid.add(new Text(name) { styleClass.add("stat-name") }, 0, row)
      statsGrid.add(new Text(value) { styleClass.add("stat-value") }, 1, row)
    }

    addStat("PDS User ID:", summary.pdsUserId, 0)
    addStat("Sample Name:", summary.libraryStats.sampleName, 1)
    addStat("Aligner:", summary.libraryStats.aligner, 2)
    addStat("Platform:", summary.libraryStats.inferredPlatform, 3)
    addStat("Instrument:", summary.libraryStats.mostFrequentInstrument, 4)
    addStat("Reference:", summary.libraryStats.referenceBuild, 5)
    addStat("Genome Size:", f"${summary.wgsMetrics.genomeTerritory}%,d", 6)
    addStat("Mean Coverage:", f"${summary.wgsMetrics.meanCoverage}%.2fx", 7)
    addStat("Callable Bases:", f"${summary.callableBases}%,d", 8)
    val callablePercent = if (summary.wgsMetrics.genomeTerritory > 0) (summary.callableBases.toDouble / summary.wgsMetrics.genomeTerritory * 100) else 0.0
    addStat("Callable Percentage:", f"$callablePercent%.2f%%", 9)

    val contigBreakdownTitle = new Label("🧬 Contig Breakdown") {
      styleClass.add("title-label")
      padding = Insets(20, 0, 10, 0)
    }

    val contigData = ObservableBuffer.from(summary.contigAnalysis.map { contig =>
      ContigAnalysisRow(
        contig.contigName,
        contig.callable,
        contig.noCoverage,
        contig.lowCoverage,
        contig.poorMappingQuality,
        contig.refN,
        s"${contig.contigName}.callable.svg"
      )
    })

    val contigTable = new TableView[ContigAnalysisRow](contigData) {
      columns ++= Seq(
        new TableColumn[ContigAnalysisRow, String]("Contig") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(c.value.contig)
        },
        new TableColumn[ContigAnalysisRow, String]("Callable Bases") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.callableBases}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("No Coverage") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.noCoverage}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("Low Coverage") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.lowCoverage}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("Poor MQ") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.poorMq}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("REF N") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.refN}%,d")
        }
      )
    }

    val webView = new WebView {
      prefHeight = 400
    }

    contigTable.selectionModel.value.selectedItem.onChange { (_, _, newValue) =>
      if (newValue != null) {
        val index = contigTable.items.value.indexOf(newValue)
        if (index >= 0 && index < svgStrings.length) {
          webView.engine.loadContent(svgStrings(index))
        }
      }
    }

    val haplogroupButton = new Button("Analyze Y-DNA Haplogroup") {
      onAction = _ => startHaplogroupAnalysis()
    }

    val resultsVBox = new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      padding = Insets(20)
      children = Seq(resultsTitle, statsGrid, new Separator(), contigBreakdownTitle, contigTable, webView, new Separator(), haplogroupButton)
    }

    if (FeatureToggles.pdsSubmissionEnabled) {
      val pdsBox = new VBox(15) {
        alignment = Pos.Center
        padding = Insets(20)
        children = Seq(
          new Label("Help advance research by securely contributing your anonymized summary data to your Personal Data Store.") {
            wrapText = true
            textAlignment = TextAlignment.Center
            styleClass.add("info-label")
          },
          new CheckBox("I agree to upload my anonymized summary data.") {
            selected = true
            style = "-fx-text-fill: #E0E0E0; -fx-font-size: 14px;"
          },
          new Button("Upload to PDS") {
            styleClass.add("button-upload")
            onAction = _ => {
              PdsClient.uploadSummary(summary).foreach { _ =>
                Platform.runLater {
                  text = "Upload Complete!"
                  styleClass.remove("button-upload")
                  styleClass.add("button-success")
                  disable = true
                }
              }
            }
          }
        )
      }
      resultsVBox.children.addAll(new Separator(), pdsBox)
    }

    val resultsScreen = new ScrollPane {
      content = resultsVBox
      fitToWidth = true
    }

    mainLayout.children = resultsScreen
  }

  private def startHaplogroupAnalysis(): Unit = {
    val progressDialog = new Dialog[Unit]() {
      initOwner(stage)
      title = "Haplogroup Analysis"
      headerText = "Running Y-DNA haplogroup analysis..."
      dialogPane().content = new ProgressIndicator()
    }

    val haplogroupTask = new jfxc.Task[Either[String, List[com.decodingus.haplogroup.model.HaplogroupResult]]]() {
      override def call(): Either[String, List[com.decodingus.haplogroup.model.HaplogroupResult]] = {
        val processor = new HaplogroupProcessor()
        processor.analyze(currentFilePath, "/Library/Genomics/Reference/chm13v2.0/chm13v2.0.fa.gz", TreeType.YDNA, TreeProviderType.FTDNA, (message, current, total) => {
          // Could update a progress bar here if needed
        })
      }
    }

    haplogroupTask.setOnSucceeded(_ => {
      progressDialog.close()
      haplogroupTask.getValue match {
        case Right(results) =>
          val topResult = results.headOption.map(_.name).getOrElse("Not found")
          new Alert(AlertType.Information) {
            initOwner(stage)
            title = "Haplogroup Analysis Complete"
            headerText = "Top Y-DNA Haplogroup Result:"
            contentText = topResult
          }.showAndWait()
        case Left(error) =>
          new Alert(AlertType.Error) {
            initOwner(stage)
            title = "Haplogroup Analysis Failed"
            headerText = "An error occurred during haplogroup analysis."
            contentText = error
          }.showAndWait()
      }
    })

    haplogroupTask.setOnFailed(_ => {
      progressDialog.close()
      new Alert(AlertType.Error) {
        initOwner(stage)
        title = "Haplogroup Analysis Failed"
        headerText = "A critical error occurred during haplogroup analysis."
        contentText = haplogroupTask.getException.getMessage
      }.showAndWait()
    })

    progressDialog.show()
    new Thread(haplogroupTask).start()
  }
}