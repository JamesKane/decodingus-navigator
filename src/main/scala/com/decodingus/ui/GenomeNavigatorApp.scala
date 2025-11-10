package com.decodingus.ui

import com.decodingus.model.{ContigSummary, CoverageSummary}
import scalafx.application.JFXApp3
import scalafx.application.JFXApp3.PrimaryStage
import scalafx.application.Platform
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.Scene
import scalafx.scene.control.{Button, CheckBox, Label, ProgressBar, ProgressIndicator, Separator}
import scalafx.scene.layout.{GridPane, HBox, StackPane, VBox}
import scalafx.scene.text.{Text, TextAlignment}
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.concurrent.Task
import scalafx.Includes._
import javafx.{concurrent => jfxc}

object GenomeNavigatorApp extends JFXApp3 {
  private val mainLayout = new StackPane()

  override def start(): Unit = {
    stage = new PrimaryStage {
      title = "Decoding-Us Navigator"
      scene = new Scene(800, 600) {
        root = mainLayout
      }
    }

    val welcomeScreen = createWelcomeScreen()
    mainLayout.children = welcomeScreen
  }

  private def createWelcomeScreen(): VBox = {
    val welcomeLabel = new Label("Welcome to Decoding-Us Navigator") {
      style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #E0E0E0;"
    }

    val instructionsLabel = new Label("Drag and drop your BAM/CRAM file here, or click to select.") {
      style = "-fx-font-size: 14px; -fx-text-fill: #B0B0B0;"
      textAlignment = TextAlignment.Center
    }

    val selectFileButton = new Button("Select BAM/CRAM File") {
      style = "-fx-background-color: #4CAF50; -fx-text-fill: white; -fx-font-size: 16px; -fx-padding: 10 20;"
      onAction = _ => selectFile()
    }

    val dragDropArea = new StackPane {
      prefWidth = 400
      prefHeight = 200
      style = "-fx-border-color: #606060; -fx-border-width: 2; -fx-border-style: dashed; -fx-background-color: #303030;"
      children = new Label("Drop File Here") {
        style = "-fx-font-size: 18px; -fx-text-fill: #909090;"
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
    }

    new VBox(20) {
      alignment = Pos.Center
      style = "-fx-background-color: #282828; -fx-padding: 20;"
      children = Seq(welcomeLabel, instructionsLabel, dragDropArea, selectFileButton)
    }
  }

  private def selectFile(): Unit = {
    println("File selection dialog would open here.")
    startAnalysis("path/to/mock/file.bam")
  }

  private def startAnalysis(filePath: String): Unit = {
    val progressLabel = new Label("Analysis in progress...") {
      style = "-fx-font-size: 20px; -fx-font-weight: bold; -fx-text-fill: #E0E0E0;"
    }

    val progressBar = new ProgressBar {
      prefWidth = 400
    }

    val progressIndicator = new ProgressIndicator

    val progressBox = new HBox(20) {
      alignment = Pos.Center
      children = Seq(progressBar, progressIndicator)
    }

    val progressScreen = new VBox(20) {
      alignment = Pos.Center
      style = "-fx-background-color: #282828;"
      children = Seq(progressLabel, progressBox)
    }

    mainLayout.children = progressScreen

    val jfxTask = new jfxc.Task[CoverageSummary]() {
      override def call(): CoverageSummary = {
        Platform.runLater {
          progressLabel.text = "Reading reference sequence..."
        }
        updateProgress(1, 5)
        Thread.sleep(1000)

        Platform.runLater {
          progressLabel.text = "Binning intervals..."
        }
        updateProgress(2, 5)
        Thread.sleep(1000)

        Platform.runLater {
          progressLabel.text = "Generating SVG plots..."
        }
        updateProgress(3, 5)
        Thread.sleep(1000)

        Platform.runLater {
          progressLabel.text = "Aggregating results..."
        }
        updateProgress(4, 5)
        Thread.sleep(500)

        Platform.runLater {
          progressLabel.text = "Analysis complete."
        }
        updateProgress(5, 5)

        CoverageSummary(
          pdsUserId = "60820188481374",
          platformSource = "bwa-mem2",
          reference = "T2T-CHM13v2.0",
          totalReads = 21206068,
          readLength = 147,
          totalBases = 3055787501L,
          callableBases = 2897310432L,
          averageDepth = 32.5,
          contigAnalysis = List(
            ContigSummary("chr1", 0, 248956422, 0, 0, 0, 0)
          )
        )
      }
    }

    // Convert JavaFX task to ScalaFX observable values
    progressBar.progress <== jfxTask.progressProperty
    progressIndicator.progress <== jfxTask.progressProperty

    // Handle task completion
    jfxTask.setOnSucceeded(_ => {
      val results = jfxTask.getValue
      showResults(results)
    })

    new Thread(jfxTask).start()
  }

  private def showResults(summary: CoverageSummary): Unit = {
    val resultsTitle = new Label("Analysis Results") {
      style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #E0E0E0;"
    }

    val statsGrid = new GridPane {
      hgap = 20
      vgap = 10
      padding = Insets(20)
      style = "-fx-background-color: #383838; -fx-background-radius: 10;"
    }

    def addStat(name: String, value: String, row: Int): Unit = {
      val nameText = new Text(name) { style = "-fx-fill: #B0B0B0; -fx-font-size: 14px;" }
      val valueText = new Text(value) { style = "-fx-fill: #FFFFFF; -fx-font-size: 14px; -fx-font-weight: bold;" }
      statsGrid.add(nameText, 0, row)
      statsGrid.add(valueText, 1, row)
    }

    addStat("Total Bases:", f"${summary.totalBases}%,d", 0)
    addStat("Callable Bases:", f"${summary.callableBases}%,d", 1)
    val callablePercent = (summary.callableBases.toDouble / summary.totalBases * 100)
    addStat("Callable Percentage:", f"$callablePercent%.2f%%", 2)
    addStat("Average Depth:", f"${summary.averageDepth}%.2fx", 3)

    val resultsScreen = new VBox(20) {
      alignment = Pos.Center
      style = "-fx-background-color: #282828; -fx-padding: 20;"
      children = Seq(resultsTitle, statsGrid, new Separator())
    }

    mainLayout.children = resultsScreen
  }
}