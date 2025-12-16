package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{TableView, TableColumn, Label}
import scalafx.scene.layout.{VBox, HBox, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.beans.property.StringProperty
import scalafx.scene.Node

/**
 * Base trait for data table components with common layout patterns.
 *
 * Provides standard table setup with title, action buttons, and column helpers.
 * Subclasses should define their specific columns and actions.
 *
 * @tparam T The row data type displayed in the table
 */
trait DataTableBase[T] extends VBox {

  /**
   * The title displayed above the table.
   */
  protected def tableTitle: String

  /**
   * The table view instance.
   */
  protected def table: TableView[T]

  /**
   * The button bar with actions (import, export, etc.).
   * Can be empty if no actions are needed.
   */
  protected def buttonBar: HBox

  // Standard VBox spacing and padding
  spacing = 10
  padding = Insets(10, 0, 0, 0)

  /**
   * Initialize the table layout.
   * Subclasses must call this at the end of their constructor.
   */
  protected def initializeLayout(): Unit = {
    children = Seq(
      new Label(tableTitle) { style = "-fx-font-weight: bold;" },
      table,
      buttonBar
    )
    VBox.setVgrow(table, Priority.Always)
  }

  /**
   * Create a string column with a value extractor.
   *
   * @param header Column header text
   * @param colWidth Preferred column width
   * @param extractor Function to extract the string value from a row
   * @return Configured TableColumn
   */
  protected def stringColumn(
    header: String,
    colWidth: Double,
    extractor: T => String
  ): TableColumn[T, String] = new TableColumn[T, String] {
    text = header
    prefWidth = colWidth
    cellValueFactory = row => StringProperty(extractor(row.value))
  }

  /**
   * Create a string column with an optional value extractor.
   *
   * @param header Column header text
   * @param colWidth Preferred column width
   * @param extractor Function to extract the optional value from a row
   * @param default Value to display when the extracted value is None
   * @return Configured TableColumn
   */
  protected def optStringColumn(
    header: String,
    colWidth: Double,
    extractor: T => Option[String],
    default: String = "—"
  ): TableColumn[T, String] = stringColumn(header, colWidth, t => extractor(t).getOrElse(default))

  /**
   * Create a numeric column formatted as a string.
   *
   * @param header Column header text
   * @param colWidth Preferred column width
   * @param extractor Function to extract the numeric value
   * @return Configured TableColumn
   */
  protected def intColumn(
    header: String,
    colWidth: Double,
    extractor: T => Int
  ): TableColumn[T, String] = stringColumn(header, colWidth, t => extractor(t).toString)

  /**
   * Create an optional numeric column formatted as a string.
   *
   * @param header Column header text
   * @param colWidth Preferred column width
   * @param extractor Function to extract the optional numeric value
   * @param default Value to display when None
   * @return Configured TableColumn
   */
  protected def optIntColumn(
    header: String,
    colWidth: Double,
    extractor: T => Option[Int],
    default: String = "—"
  ): TableColumn[T, String] = optStringColumn(header, colWidth, t => extractor(t).map(_.toString), default)

  /**
   * Create a button bar with centered left alignment.
   *
   * @param buttons Buttons to include in the bar
   * @return Configured HBox
   */
  protected def createButtonBar(buttons: Node*): HBox = new HBox(10) {
    alignment = Pos.CenterLeft
    children = buttons
  }
}
