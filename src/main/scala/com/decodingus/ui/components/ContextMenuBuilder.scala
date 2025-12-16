package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ContextMenu, MenuItem, SeparatorMenuItem, TableView, TableRow}

/**
 * Represents a context menu action for a table row.
 *
 * @param label The menu item label text
 * @param action The action to perform when clicked
 * @param enabled Function to determine if the action is enabled for a given item
 * @param separator Whether to add a separator before this item
 */
case class MenuAction[T](
  label: String,
  action: T => Unit,
  enabled: T => Boolean = (_: T) => true,
  separator: Boolean = false
)

/**
 * Builder utility for creating context menus on table rows.
 *
 * Simplifies the creation of context menus with consistent patterns
 * for handling row selection and action execution.
 */
object ContextMenuBuilder {

  /**
   * Create a row factory with a context menu for the given actions.
   *
   * Usage:
   * {{{
   * table.rowFactory = ContextMenuBuilder.createRowFactory(Seq(
   *   MenuAction("View Details", item => showDetails(item)),
   *   MenuAction("Export", item => exportItem(item)),
   *   MenuAction("Delete", item => deleteItem(item), separator = true)
   * ))
   * }}}
   *
   * @param actions The menu actions to include
   * @tparam T The table row data type
   * @return A row factory function suitable for TableView.rowFactory
   */
  def createRowFactory[T](actions: Seq[MenuAction[T]]): TableView[T] => TableRow[T] = { _ =>
    val row = new javafx.scene.control.TableRow[T]()

    val menuItems = actions.flatMap { menuAction =>
      val items = scala.collection.mutable.ListBuffer[javafx.scene.control.MenuItem]()

      // Add separator if requested
      if (menuAction.separator) {
        items += new javafx.scene.control.SeparatorMenuItem()
      }

      // Create the menu item
      val item = new MenuItem(menuAction.label) {
        onAction = _ => {
          Option(row.getItem).filter(menuAction.enabled).foreach(menuAction.action)
        }
      }
      items += item.delegate

      items.toSeq
    }

    row.contextMenu = new ContextMenu(menuItems.map(mi => new MenuItem(mi)): _*)
    new TableRow[T](row)
  }

  /**
   * Create a simple context menu for a table view (not as row factory).
   *
   * @param items Menu items to include
   * @return A ContextMenu instance
   */
  def createMenu(items: (String, () => Unit)*): ContextMenu = {
    new ContextMenu(
      items.map { case (label, action) =>
        new MenuItem(label) {
          onAction = _ => action()
        }
      }: _*
    )
  }

  /**
   * Create a context menu with actions that receive the selected item.
   *
   * @param getSelected Function to get the currently selected item (may return null)
   * @param actions Menu actions to include
   * @tparam T The item type
   * @return A ContextMenu instance
   */
  def createMenuForSelection[T](
    getSelected: () => T,
    actions: Seq[MenuAction[T]]
  ): ContextMenu = {
    new ContextMenu(
      actions.flatMap { menuAction =>
        val items = scala.collection.mutable.ListBuffer[MenuItem]()

        if (menuAction.separator) {
          items += new MenuItem(new javafx.scene.control.SeparatorMenuItem())
        }

        items += new MenuItem(menuAction.label) {
          onAction = _ => {
            Option(getSelected()).filter(menuAction.enabled).foreach(menuAction.action)
          }
        }

        items.toSeq
      }: _*
    )
  }
}
