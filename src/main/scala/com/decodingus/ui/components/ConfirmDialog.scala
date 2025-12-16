package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{Alert, ButtonType}
import scalafx.scene.control.Alert.AlertType

/**
 * Utility for displaying confirmation dialogs.
 *
 * Provides convenient methods for common confirmation patterns
 * like delete confirmations, action confirmations, etc.
 */
object ConfirmDialog {

  /**
   * Display a generic confirmation dialog.
   *
   * @param dialogTitle Window title
   * @param header Header text (the main question)
   * @param content Detailed content explaining the action
   * @return true if the user confirmed, false otherwise
   */
  def confirm(
    dialogTitle: String,
    header: String,
    content: String
  ): Boolean = {
    val alert = new Alert(AlertType.Confirmation) {
      title = dialogTitle
      headerText = header
      contentText = content
    }
    alert.showAndWait() match {
      case Some(ButtonType.OK) => true
      case _ => false
    }
  }

  /**
   * Confirm removal/deletion of an item.
   *
   * @param itemType The type of item being removed (e.g., "Subject", "Project")
   * @param details Additional details about what will be removed
   * @return true if the user confirmed, false otherwise
   */
  def confirmRemoval(itemType: String, details: String): Boolean =
    confirm(
      s"Remove $itemType",
      s"Remove this $itemType?",
      details
    )

  /**
   * Confirm a destructive action.
   *
   * @param action Description of the action (e.g., "delete all files")
   * @param warning Warning text about consequences
   * @return true if the user confirmed, false otherwise
   */
  def confirmDestructiveAction(action: String, warning: String): Boolean =
    confirm(
      "Confirm Action",
      s"Are you sure you want to $action?",
      warning
    )

  /**
   * Confirm overwriting existing data.
   *
   * @param itemDescription What will be overwritten
   * @return true if the user confirmed, false otherwise
   */
  def confirmOverwrite(itemDescription: String): Boolean =
    confirm(
      "Confirm Overwrite",
      s"Overwrite existing $itemDescription?",
      "This action cannot be undone."
    )

  /**
   * Confirm cancellation of an in-progress operation.
   *
   * @param operationName Name of the operation being cancelled
   * @return true if the user confirmed cancellation, false otherwise
   */
  def confirmCancel(operationName: String): Boolean =
    confirm(
      "Confirm Cancel",
      s"Cancel $operationName?",
      "Any progress will be lost."
    )
}
