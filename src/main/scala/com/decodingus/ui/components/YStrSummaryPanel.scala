package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.workspace.model.StrProfile
import scalafx.Includes.*
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Panel showing Y-STR summary with FTDNA panel indicators.
 * Displays which FTDNA panel tiers are filled based on marker count.
 */
class YStrSummaryPanel extends VBox {

  spacing = 10
  padding = Insets(15)
  style = "-fx-background-color: #2d2a3a; -fx-background-radius: 10;"
  prefWidth = 220

  // FTDNA panel definitions (marker thresholds)
  private val panels = Seq(
    ("Y-12", 12),
    ("Y-25", 25),
    ("Y-37", 37),
    ("Y-67", 67),
    ("Y-111", 111),
    ("Y-500", 500),
    ("Y-700", 700)
  )

  private val titleLabel = new Label(t("str.panel_summary")) {
    style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;"
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

  private val sourceLabel = new Label {
    style = "-fx-font-size: 10px; -fx-text-fill: #666666;"
  }

  private val viewProfileButton = new Button(t("str.view_profile")) {
    styleClass += "button-link"
    visible = false
    managed = false
  }

  children = Seq(
    titleLabel,
    new VBox(2) {
      alignment = Pos.Center
      children = Seq(markerCountLabel, markersLabel)
    },
    panelIndicatorsBox,
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
        val markerCount = p.markers.size
        markerCountLabel.text = markerCount.toString

        // Determine which panels are filled
        val filledPanels = panels.filter { case (_, threshold) => markerCount >= threshold }
        val highestPanel = filledPanels.lastOption.map(_._1).getOrElse("")

        // Create panel indicators
        panelIndicatorsBox.children.clear()
        panels.foreach { case (name, threshold) =>
          val isFilled = markerCount >= threshold
          val indicator = createPanelIndicator(name, isFilled, markerCount >= threshold)
          panelIndicatorsBox.children += indicator
        }

        // Show source info
        val source = p.importedFrom.orElse(p.source).getOrElse("")
        sourceLabel.text = if (source.nonEmpty) s"${t("data.source")}: $source" else ""
        sourceLabel.visible = source.nonEmpty
        sourceLabel.managed = source.nonEmpty

        // Wire up view profile button
        viewProfileButton.onAction = _ => onViewProfile()
        viewProfileButton.visible = true
        viewProfileButton.managed = true

        this.visible = true
        this.managed = true

      case _ =>
        this.visible = false
        this.managed = false
    }
  }

  /**
   * Update with multiple STR profiles (uses the one with most markers).
   */
  def setStrProfiles(profiles: List[StrProfile], onViewProfile: () => Unit): Unit = {
    if (profiles.isEmpty) {
      this.visible = false
      this.managed = false
    } else {
      // Use the profile with the most markers
      val bestProfile = profiles.maxBy(_.markers.size)
      setStrProfile(Some(bestProfile), onViewProfile)
    }
  }

  private def createPanelIndicator(name: String, isFilled: Boolean, isActive: Boolean): VBox = {
    val color = if (isFilled) "#a78bfa" else "#444444"
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
