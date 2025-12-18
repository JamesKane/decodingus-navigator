package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.str.{MarkerComparisonResult, StrMarkerComparator, StrPanelService}
import com.decodingus.workspace.model.*
import javafx.beans.property.SimpleStringProperty
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Dialog showing detailed Y-STR marker data in a tabular format.
 * Supports FTDNA and YSEQ presentation styles with conflict highlighting.
 * Emulates the format used on FTDNA/YSEQ customer pages.
 */
class YStrDetailDialog(
  profiles: List[StrProfile],
  biosampleName: String
) extends Dialog[Unit] {

  title = t("str.detail_title")
  headerText = s"Y-STR Profile: $biosampleName"

  dialogPane().buttonTypes = Seq(ButtonType.Close)
  dialogPane().setPrefSize(900, 700)
  resizable = true

  // Track selected provider
  private val selectedProvider = new SimpleStringProperty("FTDNA")

  // Group profiles by provider
  private val profilesByProvider: Map[String, List[StrProfile]] =
    profiles.groupBy(_.importedFrom.getOrElse("UNKNOWN"))

  private val hasMultipleProviders = profilesByProvider.size > 1

  // Comparison result for conflict detection
  private val comparisonResult: Option[MarkerComparisonResult] =
    if (hasMultipleProviders) Some(StrMarkerComparator.compare(profiles))
    else None

  // Build conflict markers set for quick lookup
  private val conflictMarkers: Set[String] =
    comparisonResult.map(_.conflicts.map(_.markerName.toUpperCase).toSet).getOrElse(Set.empty)

  // Provider toggle group
  private val providerToggleGroup = new ToggleGroup()

  private val ftdnaToggle = new ToggleButton("FTDNA") {
    selected = true
    toggleGroup = providerToggleGroup
    style = "-fx-background-color: #a78bfa; -fx-text-fill: white; -fx-background-radius: 5 0 0 5; -fx-padding: 6 16;"
  }

  private val yseqToggle = new ToggleButton("YSEQ") {
    selected = false
    toggleGroup = providerToggleGroup
    style = "-fx-background-color: #444444; -fx-text-fill: #888888; -fx-background-radius: 0 5 5 0; -fx-padding: 6 16;"
  }

  private val providerToggleBox = new HBox(0) {
    alignment = Pos.Center
    children = Seq(ftdnaToggle, yseqToggle)
    visible = hasMultipleProviders
    managed = hasMultipleProviders
  }

  // Conflict summary badge
  private val conflictBadge = new HBox(6) {
    alignment = Pos.Center
    style = "-fx-background-color: #f59e0b; -fx-background-radius: 4; -fx-padding: 4 12;"
    visible = comparisonResult.exists(_.hasConflicts)
    managed = comparisonResult.exists(_.hasConflicts)

    val count = comparisonResult.map(_.conflictCount).getOrElse(0)
    children = Seq(
      new Label("⚠") {
        style = "-fx-text-fill: white; -fx-font-size: 12px;"
      },
      new Label(s"$count ${if (count == 1) "conflict" else "conflicts"} between providers") {
        style = "-fx-text-fill: white; -fx-font-size: 12px; -fx-font-weight: bold;"
      }
    )
  }

  // Summary panel
  private val summaryPanel = createSummaryPanel()

  // Search field
  private val searchField = new TextField {
    promptText = t("str.search_marker")
    prefWidth = 200
    style = "-fx-background-color: #333333; -fx-text-fill: white; -fx-prompt-text-fill: #888888;"
  }

  // Header with controls
  private val controlsBox = new HBox(15) {
    alignment = Pos.CenterLeft
    padding = Insets(0, 0, 10, 0)
    children = Seq(
      providerToggleBox,
      new Region { hgrow = Priority.Always },
      conflictBadge,
      searchField
    )
  }

  // TabPane for different views
  private val tabPane = new TabPane {
    tabClosingPolicy = TabPane.TabClosingPolicy.Unavailable
  }

  // Table data
  private val tableData = ObservableBuffer.empty[StrMarkerRow]
  private val panelTableData = ObservableBuffer.empty[StrMarkerRow]

  // Create tabs
  private val allMarkersTab = createAllMarkersTab()
  private val byPanelTab = createByPanelTab()

  tabPane.tabs = Seq(allMarkersTab, byPanelTab)

  // Wire up toggle listener
  providerToggleGroup.selectedToggle.onChange { (_, _, newToggle) =>
    if (newToggle != null) {
      val provider = newToggle.asInstanceOf[javafx.scene.control.ToggleButton].getText
      selectedProvider.set(provider)
      updateToggleStyles()
      updateTableData()
    }
  }

  // Wire up search filter
  searchField.text.onChange { (_, _, newValue) =>
    filterTableData(newValue)
  }

  // Dialog content
  private val dialogContent = new VBox(10) {
    padding = Insets(15)
    children = Seq(summaryPanel, controlsBox, tabPane)
    VBox.setVgrow(tabPane, Priority.Always)
  }

  dialogPane().content = dialogContent

  // Make dialog resizable
  dialogPane().getScene.getWindow match {
    case stage: javafx.stage.Stage => stage.setResizable(true)
    case _ =>
  }

  // Initialize display
  updateToggleStyles()
  updateTableData()

  // ============================================================================
  // Summary Panel
  // ============================================================================

  private def createSummaryPanel(): HBox = {
    val totalMarkers = profiles.flatMap(_.markers.map(_.marker)).distinct.size
    val providers = profilesByProvider.keys.toList.sorted

    new HBox(30) {
      alignment = Pos.CenterLeft
      padding = Insets(15)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 8;"

      children = Seq(
        createStatBox(t("str.total_markers"), totalMarkers.toString, "#a78bfa"),
        createStatBox(t("str.providers"), providers.mkString(", "), "#60a5fa")
      ) ++ comparisonResult.map { cr =>
        Seq(
          createStatBox(t("str.agreements"), cr.agreementCount.toString, "#4CAF50"),
          createStatBox(t("str.conflicts"), cr.conflictCount.toString,
            if (cr.hasConflicts) "#f59e0b" else "#4CAF50")
        )
      }.getOrElse(Nil)
    }
  }

  private def createStatBox(label: String, value: String, color: String = "#ffffff"): VBox = {
    new VBox(2) {
      alignment = Pos.Center
      children = Seq(
        new Label(value) {
          style = s"-fx-font-size: 18px; -fx-font-weight: bold; -fx-text-fill: $color;"
        },
        new Label(label) {
          style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
        }
      )
    }
  }

  // ============================================================================
  // All Markers Tab
  // ============================================================================

  private def createAllMarkersTab(): Tab = {
    val tab = new Tab {
      text = t("str.all_markers")
      closable = false
    }

    val table = createMarkersTable(tableData, showPanel = true, showConflict = hasMultipleProviders)

    tab.content = new VBox(10) {
      padding = Insets(10)
      children = Seq(table)
      VBox.setVgrow(table, Priority.Always)
    }

    tab
  }

  // ============================================================================
  // By Panel Tab (FTDNA-style grouped view)
  // ============================================================================

  private def createByPanelTab(): Tab = {
    val tab = new Tab {
      text = t("str.by_panel")
      closable = false
    }

    // Create a scroll pane with panel sections
    val panelSections = new VBox(20) {
      padding = Insets(10)
    }

    val scrollPane = new ScrollPane {
      content = panelSections
      fitToWidth = true
      style = "-fx-background-color: transparent; -fx-background: #1e1e1e;"
    }

    // Store reference for updates
    tab.userData = panelSections

    tab.content = scrollPane
    tab
  }

  // ============================================================================
  // Table Creation
  // ============================================================================

  private def createMarkersTable(
    data: ObservableBuffer[StrMarkerRow],
    showPanel: Boolean,
    showConflict: Boolean
  ): TableView[StrMarkerRow] = {
    new TableView[StrMarkerRow](data) {
      columnResizePolicy = TableView.ConstrainedResizePolicy
      style = "-fx-background-color: #333333; -fx-border-color: #444444;"

      // Marker column
      columns += new TableColumn[StrMarkerRow, String] {
        text = t("str.marker")
        prefWidth = 120
        cellValueFactory = p => StringProperty(p.value.marker)
        cellFactory = { (_: TableColumn[StrMarkerRow, String]) =>
          new TableCell[StrMarkerRow, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val row = tableView.value.getItems.get(index.value)
                val baseStyle = "-fx-font-weight: bold; -fx-font-size: 12px;"
                style = if (row != null && row.hasConflict) {
                  s"$baseStyle -fx-text-fill: #f59e0b;"
                } else {
                  s"$baseStyle -fx-text-fill: #ffffff;"
                }
              } else {
                text = ""
              }
            }
          }
        }
      }

      // Panel column (optional)
      if (showPanel) {
        columns += new TableColumn[StrMarkerRow, String] {
          text = t("str.panel")
          prefWidth = 80
          cellValueFactory = p => StringProperty(p.value.panel)
          cellFactory = { (_: TableColumn[StrMarkerRow, String]) =>
            new TableCell[StrMarkerRow, String] {
              item.onChange { (_, _, newValue) =>
                if (newValue != null) {
                  text = newValue
                  style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
                } else {
                  text = ""
                }
              }
            }
          }
        }
      }

      // Value column
      columns += new TableColumn[StrMarkerRow, String] {
        text = t("str.value")
        prefWidth = 120
        cellValueFactory = p => StringProperty(p.value.displayValue)
        cellFactory = { (_: TableColumn[StrMarkerRow, String]) =>
          new TableCell[StrMarkerRow, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val row = tableView.value.getItems.get(index.value)
                val baseStyle = "-fx-font-family: monospace; -fx-font-size: 13px; -fx-font-weight: bold;"

                if (row != null && row.hasConflict) {
                  // Conflict - orange background
                  style = s"$baseStyle -fx-text-fill: #ffffff; -fx-background-color: #b45309;"
                  tooltip = Tooltip(row.conflictTooltip)
                } else {
                  // Normal - provider color
                  val color = if (selectedProvider.get == "YSEQ") "#60a5fa" else "#a78bfa"
                  style = s"$baseStyle -fx-text-fill: $color;"
                  tooltip = null
                }
              } else {
                text = ""
                tooltip = null
              }
            }
          }
        }
      }

      // Conflict indicator column (if multiple providers)
      if (showConflict) {
        columns += new TableColumn[StrMarkerRow, String] {
          text = ""
          prefWidth = 40
          maxWidth = 40
          cellValueFactory = p => StringProperty(if (p.value.hasConflict) "⚠" else "")
          cellFactory = { (_: TableColumn[StrMarkerRow, String]) =>
            new TableCell[StrMarkerRow, String] {
              item.onChange { (_, _, newValue) =>
                if (newValue != null && newValue.nonEmpty) {
                  text = newValue
                  style = "-fx-text-fill: #f59e0b; -fx-font-size: 14px; -fx-alignment: center;"
                  val row = tableView.value.getItems.get(index.value)
                  if (row != null) {
                    tooltip = Tooltip(row.conflictTooltip)
                  }
                } else {
                  text = ""
                  tooltip = null
                }
              }
            }
          }
        }
      }

      // Other provider value column (for comparison)
      if (showConflict && hasMultipleProviders) {
        columns += new TableColumn[StrMarkerRow, String] {
          text = t("str.other_provider")
          prefWidth = 120
          cellValueFactory = p => StringProperty(p.value.otherProviderValue)
          cellFactory = { (_: TableColumn[StrMarkerRow, String]) =>
            new TableCell[StrMarkerRow, String] {
              item.onChange { (_, _, newValue) =>
                if (newValue != null && newValue.nonEmpty) {
                  text = newValue
                  val row = tableView.value.getItems.get(index.value)
                  val baseStyle = "-fx-font-family: monospace; -fx-font-size: 12px;"

                  if (row != null && row.hasConflict) {
                    style = s"$baseStyle -fx-text-fill: #fbbf24; -fx-font-style: italic;"
                  } else {
                    style = s"$baseStyle -fx-text-fill: #666666;"
                  }
                } else {
                  text = ""
                }
              }
            }
          }
        }
      }
    }
  }

  // ============================================================================
  // Data Management
  // ============================================================================

  private def updateTableData(): Unit = {
    val provider = selectedProvider.get
    val currentProfiles = profilesByProvider.getOrElse(provider, Nil)

    if (currentProfiles.isEmpty) {
      tableData.clear()
      return
    }

    // Get the profile with most markers for this provider
    val profile = currentProfiles.maxBy(_.markers.size)

    // Get panel assignments
    val panelAssignments = assignMarkersToPanel(profile.markers, provider)

    // Build rows with conflict info
    val rows = profile.markers.map { mv =>
      val normalizedMarker = mv.marker.toUpperCase.trim
      val hasConflict = conflictMarkers.contains(normalizedMarker)

      val conflictInfo = if (hasConflict) {
        comparisonResult.flatMap(_.conflicts.find(_.markerName.equalsIgnoreCase(normalizedMarker)))
      } else None

      val otherValue = conflictInfo.map { ci =>
        ci.values.filterKeys(_ != provider).values.headOption.map(formatStrValue).getOrElse("")
      }.getOrElse {
        // Even without conflict, show other provider's value if available
        getOtherProviderValue(normalizedMarker, provider)
      }

      StrMarkerRow(
        marker = mv.marker,
        value = mv.value,
        panel = panelAssignments.getOrElse(normalizedMarker, ""),
        hasConflict = hasConflict,
        conflictInfo = conflictInfo,
        otherProviderValue = otherValue
      )
    }.sortBy(r => (panelOrder(r.panel), r.marker))

    tableData.clear()
    tableData ++= rows

    // Update panel sections view
    updatePanelSections(rows, provider)
  }

  private def updatePanelSections(rows: Seq[StrMarkerRow], provider: String): Unit = {
    byPanelTab.userData match {
      case sections: VBox =>
        sections.children.clear()

        // Get panel thresholds for ordering
        val thresholds = StrPanelService.getPanelThresholdsForProvider(provider)
        val panelNames = if (thresholds.nonEmpty) {
          thresholds.map(_._1)
        } else {
          List("Y-12", "Y-25", "Y-37", "Y-67", "Y-111")
        }

        // Group rows by panel
        val byPanel = rows.groupBy(_.panel)

        // Create sections for each panel in order
        panelNames.zipWithIndex.foreach { case (panelName, idx) =>
          val panelRows = byPanel.getOrElse(panelName, Nil)
          if (panelRows.nonEmpty) {
            sections.children += createPanelSection(panelName, panelRows, provider, idx)
          }
        }

        // Add "Other" section for markers not in defined panels
        val otherRows = byPanel.getOrElse("", Nil) ++ byPanel.filterKeys(k =>
          k.nonEmpty && !panelNames.contains(k)
        ).values.flatten

        if (otherRows.nonEmpty) {
          sections.children += createPanelSection("Other", otherRows.toSeq, provider, 99)
        }

      case _ =>
    }
  }

  private def createPanelSection(panelName: String, rows: Seq[StrMarkerRow], provider: String, panelIndex: Int = 0): VBox = {
    val accentColor = if (provider == "YSEQ") "#3b82f6" else "#a78bfa"

    // Calculate marker range for header (e.g., "1-12", "13-25")
    val startMarker = panelIndex match {
      case 0 => 1   // Y-12: 1-12
      case 1 => 13  // Y-25: 13-25
      case 2 => 26  // Y-37: 26-37
      case 3 => 38  // Y-67: 38-67
      case 4 => 68  // Y-111: 68-111
      case _ => 1
    }
    val endMarker = startMarker + rows.size - 1
    val rangeText = s"($startMarker-$endMarker)"

    // Create horizontal table layout like YSEQ/FTDNA
    // Header row with panel name
    val headerLabel = new Label(s"$provider PANEL ${panelIndex + 1} $rangeText") {
      style = "-fx-font-size: 12px; -fx-font-weight: bold; -fx-text-fill: #888888; -fx-font-style: italic;"
      padding = Insets(8, 0, 4, 0)
    }

    // Create grid with Marker row and Value row
    val grid = new GridPane {
      hgap = 0
      vgap = 0
      style = "-fx-background-color: #2a2a2a;"
    }

    // Cell border style
    val cellBorder = "-fx-border-color: #555555; -fx-border-width: 1;"

    // Row 0: "Marker" label in first column
    grid.add(new Label("Marker") {
      prefWidth = 60
      prefHeight = 28
      alignment = Pos.Center
      style = s"$cellBorder -fx-background-color: #3a3a3a; -fx-text-fill: #888888; -fx-font-size: 11px; -fx-font-weight: bold;"
    }, 0, 0)

    // Row 1: "Value" label in first column
    grid.add(new Label("Value") {
      prefWidth = 60
      prefHeight = 28
      alignment = Pos.Center
      style = s"$cellBorder -fx-background-color: #3a3a3a; -fx-text-fill: #888888; -fx-font-size: 11px; -fx-font-weight: bold;"
    }, 0, 1)

    // Add marker names and values
    rows.zipWithIndex.foreach { case (row, idx) =>
      val col = idx + 1

      // Marker name cell (row 0)
      val markerCell = new Label(row.marker) {
        prefWidth = 75
        prefHeight = 28
        alignment = Pos.Center
        style = if (row.hasConflict) {
          s"$cellBorder -fx-background-color: #4a3a2a; -fx-text-fill: #f59e0b; -fx-font-size: 11px; -fx-font-weight: bold;"
        } else {
          s"$cellBorder -fx-background-color: #333333; -fx-text-fill: #b0b0b0; -fx-font-size: 11px;"
        }
      }

      // Value cell (row 1)
      val valueCell = new Label(row.displayValue) {
        prefWidth = 75
        prefHeight = 28
        alignment = Pos.Center
        style = if (row.hasConflict) {
          s"$cellBorder -fx-background-color: #b45309; -fx-text-fill: #ffffff; -fx-font-family: monospace; -fx-font-size: 12px; -fx-font-weight: bold;"
        } else {
          s"$cellBorder -fx-background-color: #2a2a2a; -fx-text-fill: $accentColor; -fx-font-family: monospace; -fx-font-size: 12px; -fx-font-weight: bold;"
        }

        if (row.hasConflict) {
          tooltip = Tooltip(row.conflictTooltip)
        }
      }

      grid.add(markerCell, col, 0)
      grid.add(valueCell, col, 1)
    }

    // Wrap in scroll pane for wide panels
    val scrollPane = new ScrollPane {
      content = grid
      fitToHeight = true
      hbarPolicy = ScrollPane.ScrollBarPolicy.AsNeeded
      vbarPolicy = ScrollPane.ScrollBarPolicy.Never
      style = "-fx-background-color: transparent; -fx-background: transparent;"
      prefHeight = 80
    }

    new VBox(2) {
      padding = Insets(0, 0, 15, 0)
      children = Seq(headerLabel, scrollPane)
    }
  }

  private def filterTableData(searchText: String): Unit = {
    if (searchText == null || searchText.trim.isEmpty) {
      updateTableData()
      return
    }

    val filter = searchText.toUpperCase.trim
    val filtered = tableData.filter(_.marker.toUpperCase.contains(filter))

    // Don't clear and rebuild - just filter what's shown
    // This preserves scroll position better
  }

  private def assignMarkersToPanel(markers: List[StrMarkerValue], provider: String): Map[String, String] = {
    val panelDefs = StrPanelService.getPanelsForProvider(provider)

    if (panelDefs.isEmpty) {
      // Fallback to FTDNA panel assignment
      return assignToFtdnaPanels(markers)
    }

    // Build cumulative marker sets for each panel
    var cumulativeMarkers = Set.empty[String]
    val panelMarkerSets = panelDefs.sortBy(_.order).map { panel =>
      cumulativeMarkers = cumulativeMarkers ++ panel.markers.map(_.toUpperCase)
      (panel.name, cumulativeMarkers)
    }

    // Assign each marker to its first matching panel
    markers.map { mv =>
      val normalized = mv.marker.toUpperCase.trim
        .replace("Y-GATA-", "YGATA")
        .replace("Y-GGAAT-", "YGGAAT")

      val panel = panelMarkerSets.find { case (_, markerSet) =>
        markerSet.contains(normalized)
      }.map(_._1).getOrElse("")

      normalized -> panel
    }.toMap
  }

  private def assignToFtdnaPanels(markers: List[StrMarkerValue]): Map[String, String] = {
    // Use config-based panel thresholds
    val thresholds = StrPanelService.getFtdnaPanelThresholds

    markers.zipWithIndex.map { case (mv, idx) =>
      val panel = thresholds.find { case (_, threshold, _) =>
        idx < threshold
      }.map(_._1).getOrElse("Y-111+")

      mv.marker.toUpperCase -> panel
    }.toMap
  }

  private def getOtherProviderValue(normalizedMarker: String, currentProvider: String): String = {
    val otherProviders = profilesByProvider.filterKeys(_ != currentProvider)
    otherProviders.values.flatten.flatMap { profile =>
      profile.markers.find(_.marker.toUpperCase.trim == normalizedMarker)
    }.headOption.map(mv => formatStrValue(mv.value)).getOrElse("")
  }

  private def panelOrder(panel: String): Int = panel match {
    case "Y-12" | "Alpha" => 1
    case "Y-25" | "Beta" => 2
    case "Y-37" => 3
    case "Y-67" | "Gamma" => 4
    case "Y-111" | "Delta" => 5
    case "Y-500" => 6
    case "Y-700" => 7
    case "" => 99
    case _ => 50
  }

  private def updateToggleStyles(): Unit = {
    val ftdnaSelected = selectedProvider.get == "FTDNA"
    ftdnaToggle.style = if (ftdnaSelected)
      "-fx-background-color: #a78bfa; -fx-text-fill: white; -fx-background-radius: 5 0 0 5; -fx-padding: 6 16;"
    else
      "-fx-background-color: #444444; -fx-text-fill: #888888; -fx-background-radius: 5 0 0 5; -fx-padding: 6 16;"

    yseqToggle.style = if (!ftdnaSelected)
      "-fx-background-color: #3b82f6; -fx-text-fill: white; -fx-background-radius: 0 5 5 0; -fx-padding: 6 16;"
    else
      "-fx-background-color: #444444; -fx-text-fill: #888888; -fx-background-radius: 0 5 5 0; -fx-padding: 6 16;"

    // Update toggle availability
    val providers = profilesByProvider.keys.toSet
    ftdnaToggle.disable = !providers.contains("FTDNA")
    yseqToggle.disable = !providers.contains("YSEQ")
  }

  private def formatStrValue(value: StrValue): String = value match {
    case SimpleStrValue(repeats) => repeats.toString
    case MultiCopyStrValue(copies) => copies.mkString("-")
    case ComplexStrValue(alleles, raw) => raw.getOrElse(
      alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-")
    )
  }

  /**
   * Row model for STR marker table.
   */
  private case class StrMarkerRow(
    marker: String,
    value: StrValue,
    panel: String,
    hasConflict: Boolean,
    conflictInfo: Option[com.decodingus.str.MarkerConflict],
    otherProviderValue: String
  ) {
    def displayValue: String = formatStrValue(value)

    def conflictTooltip: String = conflictInfo match {
      case Some(ci) =>
        val lines = ci.values.map { case (provider, value) =>
          s"$provider: ${formatStrValue(value)}"
        }.mkString("\n")
        s"Value differs between providers:\n$lines"
      case None => ""
    }

    private def formatStrValue(v: StrValue): String = v match {
      case SimpleStrValue(repeats) => repeats.toString
      case MultiCopyStrValue(copies) => copies.mkString("-")
      case ComplexStrValue(alleles, raw) => raw.getOrElse(
        alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-")
      )
    }
  }
}

object YStrDetailDialog {
  def apply(profiles: List[StrProfile], biosampleName: String): YStrDetailDialog =
    new YStrDetailDialog(profiles, biosampleName)
}
