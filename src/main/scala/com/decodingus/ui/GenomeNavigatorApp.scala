package com.decodingus.ui

import scalafx.Includes.*
import scalafx.application.JFXApp3
import scalafx.geometry.Pos
import scalafx.scene.Scene
import scalafx.scene.control.{Button, Label}
import scalafx.scene.input.TransferMode
import scalafx.scene.layout.{StackPane, VBox}
import scalafx.scene.text.TextAlignment

object GenomeNavigatorApp extends JFXApp3 {
  override def start(): Unit = {
    stage = new JFXApp3.PrimaryStage {
      title = "Decoding-Us Navigator"
      scene = new Scene(800, 600) {
        root = new VBox(20) {
          alignment = Pos.Center
          style = "-fx-background-color: #282828; -fx-padding: 20;"
          children = Seq(
            new Label("Welcome to Decoding-Us Navigator") {
              style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #E0E0E0;"
            },
            new Label("Drag and drop your BAM/CRAM file here, or click to select.") {
              style = "-fx-font-size: 14px; -fx-text-fill: #B0B0B0;"
              textAlignment = TextAlignment.Center
            },
            new StackPane {
              prefWidth = 400
              prefHeight = 200
              style = "-fx-border-color: #606060; -fx-border-width: 2; -fx-border-style: dashed; -fx-background-color: #303030;"
              children = new Label("Drop File Here") {
                style = "-fx-font-size: 18px; -fx-text-fill: #909090;"
              }

              onDragOver = event => {
                if (event.gestureSource != this && event.dragboard.hasFiles) {
                  event.acceptTransferModes(TransferMode.Copy, TransferMode.Move)
                }
                event.consume()
              }

              onDragDropped = event => {
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
              style = "-fx-background-color: #4CAF50; -fx-text-fill: white; -fx-font-size: 16px; -fx-padding: 10 20;"
              onAction = _ => selectFile()
            }
          )
        }
      }
    }
  }

  private def selectFile(): Unit = {
    println("File selection dialog would open here.")
  }

  private def startAnalysis(filePath: String): Unit = {
    println(s"Starting analysis for: $filePath")
  }
}