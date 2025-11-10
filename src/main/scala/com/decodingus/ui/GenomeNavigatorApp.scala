package com.decodingus.ui

import com.decodingus.analysis.CallableLociProcessor
import com.decodingus.model.CoverageSummary
import com.decodingus.pds.PdsClient
import javafx.concurrent as jfxc
import scalafx.Includes.*
import scalafx.application.JFXApp3.PrimaryStage
import scalafx.application.{JFXApp3, Platform}
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.Scene
import scalafx.scene.control.*
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.{GridPane, HBox, StackPane, VBox}
import scalafx.scene.text.{Text, TextAlignment}

import scala.concurrent.ExecutionContext.Implicits.global

object GenomeNavigatorApp extends JFXApp3 {
  private val mainLayout = new StackPane()

  override def start(): Unit = {
    stage = new PrimaryStage {
      title = "Decoding-Us Navigator"
      scene = new Scene(800, 600) {
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
    val progressLabel = new Label("Analysis in progress...") {
      styleClass.add("progress-label")
    }
    val progressBar = new ProgressBar {
      prefWidth = 400
    }
    val progressIndicator = new ProgressIndicator

    val progressScreen = new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      children = Seq(
        progressLabel,
        new HBox(20) {
          alignment = Pos.Center
          children = Seq(progressBar, progressIndicator)
        }
      )
    }

    mainLayout.children = progressScreen

    val jfxTask = new jfxc.Task[CoverageSummary]() {
      override def call(): CoverageSummary = {
        try {
          val processor = new CallableLociProcessor()
          val referencePath = "/Library/Genomics/Reference/chm13v2.0/chm13v2.0.fa.gz"

          Platform.runLater {
            progressLabel.text = "Running analysis..."
          }
          updateProgress(1, 2)

          val (summary, _) = processor.process(filePath, referencePath)

          Platform.runLater {
            progressLabel.text = "Analysis complete."
          }
          updateProgress(2, 2)
          summary
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
      val results = jfxTask.getValue
      showResults(results)
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

  private def showResults(summary: CoverageSummary): Unit = {
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
    addStat("Platform Source:", summary.platformSource, 1)
    addStat("Reference:", summary.reference, 2)
    addStat("Total Bases:", f"${summary.totalBases}%,d", 3)
    addStat("Callable Bases:", f"${summary.callableBases}%,d", 4)
    val callablePercent = if (summary.totalBases > 0) (summary.callableBases.toDouble / summary.totalBases * 100) else 0.0
    addStat("Callable Percentage:", f"$callablePercent%.2f%%", 5)
    addStat("Average Depth:", f"${summary.averageDepth}%.2fx", 6)

    val optInCheck = new CheckBox("I agree to upload my anonymized summary data.") {
      selected = true
      style = "-fx-text-fill: #E0E0E0; -fx-font-size: 14px;"
    }

    val uploadButton = new Button("Upload to PDS") {
      styleClass.add("button-upload")
      disable <== !optInCheck.selected
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

    val pdsBox = new VBox(15) {
      alignment = Pos.Center
      padding = Insets(20)
      children = Seq(
        new Label("Help advance research by securely contributing your anonymized summary data to your Personal Data Store.") {
          wrapText = true
          textAlignment = TextAlignment.Center
          styleClass.add("info-label")
        },
        optInCheck,
        uploadButton
      )
    }

    val resultsScreen = new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      padding = Insets(20)
      children = Seq(resultsTitle, statsGrid, new Separator(), pdsBox)
    }

    mainLayout.children = resultsScreen
  }
}
