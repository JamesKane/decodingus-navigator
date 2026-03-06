package com.decodingus.ui.components

import scalafx.Includes.*
import scalafx.application.Platform
import scalafx.beans.property.{BooleanProperty, DoubleProperty, StringProperty}
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{HBox, VBox}

/**
 * A dialog that displays analysis progress with a progress bar and status message.
 * Binds to observable properties from the ViewModel.
 */
class AnalysisProgressDialog(
                              titleText: String,
                              progressMessage: StringProperty,
                              progressPercent: DoubleProperty,
                              inProgress: BooleanProperty
                            ) extends Dialog[Unit] {

  title = titleText
  headerText = "Analysis in progress..."
  resizable = false

  private val statusLabel = new Label() {
    text <== progressMessage
    wrapText = true
    prefWidth = 400
  }

  private val progressBar = new ProgressBar() {
    prefWidth = 400
    progress <== progressPercent
  }

  private val progressIndicator = new ProgressIndicator() {
    progress <== progressPercent
    prefWidth = 50
    prefHeight = 50
  }

  private val contentBox = new VBox(15) {
    padding = Insets(20)
    alignment = Pos.Center
    children = Seq(
      statusLabel,
      new HBox(20) {
        alignment = Pos.Center
        children = Seq(progressBar, progressIndicator)
      }
    )
  }

  dialogPane().content = contentBox

  // Only show cancel button
  dialogPane().buttonTypes = Seq(ButtonType.Cancel)

  // Auto-close when analysis completes
  inProgress.onChange { (_, _, isRunning) =>
    if (!isRunning) {
      Platform.runLater {
        result = ()
      }
    }
  }
}
