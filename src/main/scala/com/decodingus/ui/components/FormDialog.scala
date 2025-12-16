package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ButtonBar}
import scalafx.scene.layout.{GridPane, Priority}
import scalafx.scene.Node
import scalafx.geometry.Insets

/**
 * Base trait for form-based dialogs with common setup patterns.
 *
 * Provides standard button setup, grid layout helpers, and result conversion.
 * Subclasses should override the abstract members and call `initializeDialog()`
 * at the end of their constructor.
 *
 * @tparam T The result type returned when the form is submitted successfully
 */
trait FormDialog[T] extends Dialog[Option[T]] {

  /**
   * The dialog window title.
   */
  protected def dialogTitle: String

  /**
   * The dialog header text describing the form's purpose.
   */
  protected def dialogHeader: String

  /**
   * Text shown on the primary action button.
   * Default is "Save", but can be overridden (e.g., "Create", "Add").
   */
  protected def primaryButtonText: String = "Save"

  /**
   * Build the form content to display in the dialog.
   * Typically returns a GridPane with form fields.
   */
  protected def buildFormContent(): Node

  /**
   * Construct the result object from form field values.
   * Called when the primary button is clicked and validation passes.
   */
  protected def buildResult(): T

  /**
   * Validate form fields before allowing submission.
   * Return true if the form is valid, false otherwise.
   * Default implementation always returns true.
   */
  protected def isValid: Boolean = true

  /**
   * Optional field to focus when the dialog opens.
   * Override to specify which field should receive initial focus.
   */
  protected def initialFocusField: Option[TextField] = None

  // Standard button type for the primary action
  protected val primaryButtonType: ButtonType =
    new ButtonType(primaryButtonText, ButtonBar.ButtonData.OKDone)

  /**
   * Lookup the primary button node for enabling/disabling based on validation.
   */
  protected lazy val primaryButton: Node =
    dialogPane().lookupButton(primaryButtonType)

  /**
   * Initialize the dialog with standard configuration.
   * Subclasses must call this at the end of their constructor.
   */
  protected def initializeDialog(): Unit = {
    title = dialogTitle
    headerText = dialogHeader

    dialogPane().buttonTypes = Seq(primaryButtonType, ButtonType.Cancel)
    dialogPane().content = buildFormContent()

    // Set up result converter
    resultConverter = btn => {
      if (btn == primaryButtonType && isValid) Some(buildResult())
      else None
    }

    // Set initial focus if specified
    initialFocusField.foreach { field =>
      javafx.application.Platform.runLater(() => field.requestFocus())
    }
  }

  /**
   * Helper to build a standard form grid with labeled rows.
   *
   * @param rows Pairs of (label text, form node) for each row
   * @param gridHgap Horizontal gap between columns (default 10)
   * @param gridVgap Vertical gap between rows (default 10)
   * @param gridPadding Insets around the grid
   * @return A configured GridPane with the rows added
   */
  protected def buildGrid(
    rows: Seq[(String, Node)],
    gridHgap: Int = 10,
    gridVgap: Int = 10,
    gridPadding: Insets = Insets(20, 150, 10, 10)
  ): GridPane = {
    new GridPane() {
      hgap = gridHgap
      vgap = gridVgap
      padding = gridPadding

      rows.zipWithIndex.foreach { case ((labelText, node), idx) =>
        add(new Label(s"$labelText:"), 0, idx)
        add(node, 1, idx)
      }
    }
  }

  /**
   * Helper to bind validation to enable/disable the primary button.
   * Call this after initializeDialog() with your validation condition.
   *
   * @param validProperty A property that is true when the form is valid
   */
  protected def bindPrimaryButtonEnabled(validProperty: => Boolean): Unit = {
    primaryButton.disable = !validProperty
  }
}
