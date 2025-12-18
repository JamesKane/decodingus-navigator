package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.str.{MarkerComparisonResult, StrMarkerComparator, StrPanelService}
import com.decodingus.workspace.model.StrProfile
import javafx.beans.property.SimpleStringProperty
import scalafx.Includes.*
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Panel showing Y-STR summary with provider toggle and panel indicators.
 * Supports FTDNA and YSEQ panel views with conflict detection.
 */
class YStrSummaryPanel extends VBox {

  spacing = 10
  padding = Insets(15)
  style = "-fx-background-color: #2d2a3a; -fx-background-radius: 10;"
  prefWidth = 220

  // Track selected provider
  private val selectedProvider = new SimpleStringProperty("FTDNA")

  // Store profiles grouped by provider
  private var profilesByProvider: Map[String, List[StrProfile]] = Map.empty
  private var comparisonResult: Option[MarkerComparisonResult] = None
  private var viewProfileCallback: () => Unit = () => ()

  private val titleLabel = new Label(t("str.panel_summary")) {
    style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }

  // Provider toggle buttons
  private val ftdnaToggle = new ToggleButton("FTDNA") {
    selected = true
    style = "-fx-background-color: #a78bfa; -fx-text-fill: white; -fx-background-radius: 5 0 0 5; -fx-padding: 4 10;"
    toggleGroup = providerToggleGroup
  }

  private val yseqToggle = new ToggleButton("YSEQ") {
    selected = false
    style = "-fx-background-color: #444444; -fx-text-fill: #888888; -fx-background-radius: 0 5 5 0; -fx-padding: 4 10;"
    toggleGroup = providerToggleGroup
  }

  private val providerToggleGroup = new ToggleGroup()

  private val providerToggleBox = new HBox(0) {
    alignment = Pos.Center
    children = Seq(ftdnaToggle, yseqToggle)
    visible = false
    managed = false
  }

  private val markerCountLabel = new Label("-") {
    style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #a78bfa;"
  }

  private val markersLabel = new Label(t("str.markers")) {
    style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
  }

  private val panelIndicatorsBox = new HBox(3) {
    alignment = Pos.Center
  }

  // Conflict badge
  private val conflictBadge = new HBox(4) {
    alignment = Pos.Center
    style = "-fx-background-color: #f59e0b; -fx-background-radius: 4; -fx-padding: 3 8;"
    visible = false
    managed = false
    children = Seq(
      new Label("⚠") {
        style = "-fx-text-fill: white; -fx-font-size: 10px;"
      },
      conflictCountLabel
    )
    cursor = scalafx.scene.Cursor.Hand
  }

  private val conflictCountLabel = new Label("0 conflicts") {
    style = "-fx-text-fill: white; -fx-font-size: 10px; -fx-font-weight: bold;"
  }

  private val sourceLabel = new Label {
    style = "-fx-font-size: 10px; -fx-text-fill: #666666;"
  }

  private val viewProfileButton = new Button(t("str.view_profile")) {
    styleClass += "button-link"
    visible = false
    managed = false
  }

  // Wire up toggle group listener
  providerToggleGroup.selectedToggle.onChange { (_, _, newToggle) =>
    if (newToggle != null) {
      val provider = newToggle.asInstanceOf[javafx.scene.control.ToggleButton].getText
      selectedProvider.set(provider)
      updateDisplayForProvider(provider)
    }
  }

  children = Seq(
    titleLabel,
    providerToggleBox,
    new VBox(2) {
      alignment = Pos.Center
      children = Seq(markerCountLabel, markersLabel)
    },
    panelIndicatorsBox,
    conflictBadge,
    sourceLabel,
    viewProfileButton
  )

  // Start hidden
  visible = false
  managed = false

  /**
   * Update the panel with STR profile data.
   */
  def setStrProfile(profile: Option[StrProfile], onViewProfile: () => Unit): Unit = {
    profile match {
      case Some(p) if p.markers.nonEmpty =>
        setStrProfiles(List(p), onViewProfile)
      case _ =>
        this.visible = false
        this.managed = false
    }
  }

  /**
   * Update with multiple STR profiles from potentially different providers.
   * Shows provider toggle if multiple providers present, with conflict detection.
   */
  def setStrProfiles(profiles: List[StrProfile], onViewProfile: () => Unit): Unit = {
    if (profiles.isEmpty) {
      this.visible = false
      this.managed = false
      return
    }

    viewProfileCallback = onViewProfile

    // Group profiles by provider
    profilesByProvider = profiles.groupBy(_.importedFrom.getOrElse("UNKNOWN"))

    // Detect conflicts if multiple providers
    if (profilesByProvider.size > 1) {
      val comparison = StrMarkerComparator.compare(profiles)
      comparisonResult = Some(comparison)

      // Show toggle
      providerToggleBox.visible = true
      providerToggleBox.managed = true

      // Update toggle button availability based on which providers have data
      val providers = profilesByProvider.keys.toSet
      ftdnaToggle.disable = !providers.contains("FTDNA")
      yseqToggle.disable = !providers.contains("YSEQ")

      // Update toggle button styles
      updateToggleStyles()

      // Show conflict badge if conflicts exist
      if (comparison.hasConflicts) {
        val count = comparison.conflictCount
        conflictCountLabel.text = s"$count ${if (count == 1) "conflict" else "conflicts"}"
        Tooltip.install(conflictBadge, buildConflictTooltip(comparison))
        conflictBadge.visible = true
        conflictBadge.managed = true
      } else {
        conflictBadge.visible = false
        conflictBadge.managed = false
      }
    } else {
      // Single provider - no toggle needed
      providerToggleBox.visible = false
      providerToggleBox.managed = false
      conflictBadge.visible = false
      conflictBadge.managed = false
      comparisonResult = None

      // Set selected provider to the only one present
      val provider = profilesByProvider.keys.headOption.getOrElse("FTDNA")
      selectedProvider.set(provider)
    }

    // Display for current selected provider
    updateDisplayForProvider(selectedProvider.get)

    // Wire up view profile button
    viewProfileButton.onAction = _ => onViewProfile()
    viewProfileButton.visible = true
    viewProfileButton.managed = true

    this.visible = true
    this.managed = true
  }

  /**
   * Update display for the selected provider.
   */
  private def updateDisplayForProvider(provider: String): Unit = {
    val profiles = profilesByProvider.getOrElse(provider, Nil)
    if (profiles.isEmpty) return

    // Use profile with most markers for this provider
    val profile = profiles.maxBy(_.markers.size)
    val markerCount = profile.markers.size

    markerCountLabel.text = markerCount.toString

    // Get panel thresholds for this provider
    val thresholds = StrPanelService.getPanelThresholdsForProvider(provider)

    // Create panel indicators
    panelIndicatorsBox.children.clear()
    if (thresholds.nonEmpty) {
      thresholds.foreach { case (name, actualThreshold, _) =>
        val isFilled = markerCount >= actualThreshold
        val indicator = createPanelIndicator(name, isFilled, provider)
        panelIndicatorsBox.children += indicator
      }
    } else {
      // No panel definitions - just show marker count
      panelIndicatorsBox.children += new Label(s"$markerCount markers") {
        style = "-fx-text-fill: #888888; -fx-font-size: 10px;"
      }
    }

    // Show source info
    val source = profile.importedFrom.orElse(profile.source).getOrElse("")
    sourceLabel.text = if (source.nonEmpty) s"${t("data.source")}: $source" else ""
    sourceLabel.visible = source.nonEmpty
    sourceLabel.managed = source.nonEmpty

    // Update toggle styles to reflect selection
    updateToggleStyles()
  }

  /**
   * Update toggle button styles based on selection.
   */
  private def updateToggleStyles(): Unit = {
    val ftdnaSelected = selectedProvider.get == "FTDNA"
    ftdnaToggle.style = if (ftdnaSelected)
      "-fx-background-color: #a78bfa; -fx-text-fill: white; -fx-background-radius: 5 0 0 5; -fx-padding: 4 10;"
    else
      "-fx-background-color: #444444; -fx-text-fill: #888888; -fx-background-radius: 5 0 0 5; -fx-padding: 4 10;"

    yseqToggle.style = if (!ftdnaSelected)
      "-fx-background-color: #3b82f6; -fx-text-fill: white; -fx-background-radius: 0 5 5 0; -fx-padding: 4 10;"
    else
      "-fx-background-color: #444444; -fx-text-fill: #888888; -fx-background-radius: 0 5 5 0; -fx-padding: 4 10;"

    // Also update marker count color based on provider
    markerCountLabel.style = if (ftdnaSelected)
      "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #a78bfa;"
    else
      "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #3b82f6;"
  }

  /**
   * Build tooltip content for conflict badge.
   */
  private def buildConflictTooltip(comparison: MarkerComparisonResult): Tooltip = {
    val tooltipText = new StringBuilder("Marker Value Conflicts:\n\n")

    comparison.conflicts.take(10).foreach { conflict =>
      tooltipText.append(conflict.formatForDisplay)
      tooltipText.append("\n\n")
    }

    if (comparison.conflicts.size > 10) {
      tooltipText.append(s"... and ${comparison.conflicts.size - 10} more")
    }

    new Tooltip(tooltipText.toString.trim) {
      style = "-fx-font-family: monospace; -fx-font-size: 11px;"
      showDelay = javafx.util.Duration.millis(200)
    }
  }

  private def createPanelIndicator(name: String, isFilled: Boolean, provider: String): VBox = {
    val filledColor = if (provider == "YSEQ") "#3b82f6" else "#a78bfa"
    val color = if (isFilled) filledColor else "#444444"
    val textColor = if (isFilled) "#ffffff" else "#666666"
    val tooltipText = s"$name: ${if (isFilled) t("str.panel_filled") else t("str.panel_not_filled")}"

    new VBox(1) {
      alignment = Pos.Center
      children = Seq(
        new Label {
          text = if (isFilled) "●" else "○"
          style = s"-fx-text-fill: $color; -fx-font-size: 10px;"
          tooltip = Tooltip(tooltipText)
        },
        new Label(name.replace("Y-", "")) {
          style = s"-fx-text-fill: $textColor; -fx-font-size: 8px;"
          tooltip = Tooltip(tooltipText)
        }
      )
    }
  }
}

object YStrSummaryPanel {
  def apply(): YStrSummaryPanel = new YStrSummaryPanel()
}
