package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, TextArea}
import scalafx.scene.layout.{VBox, Priority}
import scalafx.geometry.Insets

/**
 * Utility for displaying read-only information dialogs.
 *
 * Provides a simple way to show text content in a scrollable, read-only text area.
 * Useful for displaying logs, reports, JSON data, or other read-only content.
 */
object InfoDialog {

  /**
   * Display an information dialog with the given content.
   *
   * @param dialogTitle Window title
   * @param header Header text describing the content
   * @param content The text content to display
   * @param dialogWidth Preferred dialog width (default 500)
   * @param dialogHeight Preferred dialog height (default 400)
   * @param monospaced Use monospace font for the content (default false)
   * @param enableWrap Enable text wrapping (default true)
   */
  def show(
    dialogTitle: String,
    header: String,
    content: String,
    dialogWidth: Int = 500,
    dialogHeight: Int = 400,
    monospaced: Boolean = false,
    enableWrap: Boolean = true
  ): Unit = {
    val dialog = new Dialog[Unit] {
      title = dialogTitle
      headerText = header
      dialogPane().buttonTypes = Seq(ButtonType.Close)
      dialogPane().setPrefSize(dialogWidth, dialogHeight)

      val textArea = new TextArea(content) {
        editable = false
        wrapText = enableWrap
        if (monospaced) {
          style = "-fx-font-family: monospace; -fx-font-size: 12px;"
        }
      }
      VBox.setVgrow(textArea, Priority.Always)

      dialogPane().content = new VBox(10) {
        padding = Insets(10)
        children = Seq(textArea)
      }
    }
    dialog.showAndWait()
  }

  /**
   * Display a monospaced information dialog, suitable for code or JSON.
   */
  def showCode(
    dialogTitle: String,
    header: String,
    content: String,
    dialogWidth: Int = 600,
    dialogHeight: Int = 500
  ): Unit = show(dialogTitle, header, content, dialogWidth, dialogHeight, monospaced = true, enableWrap = false)

  /**
   * Display a JSON information dialog.
   */
  def showJson(
    dialogTitle: String,
    header: String,
    json: String,
    dialogWidth: Int = 600,
    dialogHeight: Int = 500
  ): Unit = showCode(dialogTitle, header, json, dialogWidth, dialogHeight)
}
