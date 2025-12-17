package com.decodingus.ui.components

import scalafx.Includes.*
import scalafx.geometry.Pos
import scalafx.scene.control.{Button, Label}
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.{StackPane, VBox}
import scalafx.scene.text.TextAlignment

class WelcomeScreen(onFileSelected: String => Unit, onSelectFileClicked: () => Unit) extends VBox(20) {
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
            onFileSelected(file.getAbsolutePath)
            success = true
          }
        }
        event.dropCompleted = success
        event.consume()
      }
    },
    new Button("Select BAM/CRAM File") {
      styleClass.add("button-select")
      onAction = _ => onSelectFileClicked()
    }
  )
}
