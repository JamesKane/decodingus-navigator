package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.refgenome.YRegionAnnotator
import com.decodingus.yprofile.model.{YConsensusState, YProfileVariantEntity, YVariantStatus}
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Label
import scalafx.scene.layout.{Priority, VBox}
import scalafx.scene.web.WebView

/**
 * Panel showing Y chromosome ideogram with region bands and variant markers.
 * Displays the chromosome with color-coded regions (PAR, X-degenerate, Ampliconic, etc.)
 * and triangle markers for derived/novel/conflict variants.
 */
class YChromosomeIdeogramPanel extends VBox {

  spacing = 10
  padding = Insets(15)
  style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"

  private val titleLabel = new Label(t("haplogroup.ychromosome_ideogram")) {
    style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }

  private val webView = new WebView {
    prefHeight = 280
    minHeight = 280
  }
  VBox.setVgrow(webView, Priority.Always)

  private val noDataLabel = new Label(t("info.no_data")) {
    style = "-fx-font-size: 14px; -fx-text-fill: #888888;"
  }

  private val contentBox = new VBox(10) {
    alignment = Pos.Center
    children = Seq(noDataLabel)
  }
  VBox.setVgrow(contentBox, Priority.Always)

  children = Seq(titleLabel, contentBox)

  // Start hidden
  visible = false
  managed = false

  /**
   * Update the ideogram with Y profile variant data.
   *
   * @param annotator Optional Y region annotator for chromosome structure
   * @param variants  List of Y profile variants to display as markers
   */
  def setData(annotator: Option[YRegionAnnotator], variants: List[YProfileVariantEntity]): Unit = {
    annotator match {
      case Some(ann) =>
        // Convert variants to markers
        val variantMarkers = variants
          .filter(v => v.consensusState == YConsensusState.DERIVED ||
            v.status == YVariantStatus.NOVEL ||
            v.status == YVariantStatus.CONFLICT)
          .map(YChromosomeIdeogramRenderer.VariantMarker.fromVariantEntity)

        // Generate SVG
        val svgContent = YChromosomeIdeogramRenderer.render(ann, variantMarkers)
        val statsHtml = YChromosomeIdeogramRenderer.renderStatsHtml(variants, ann)

        // Wrap in HTML for WebView
        val html =
          s"""<!DOCTYPE html>
             |<html>
             |<head>
             |  <style>
             |    body { margin: 0; padding: 10px; background: #2a2a2a; font-family: system-ui, sans-serif; }
             |    svg { max-width: 100%; height: auto; }
             |  </style>
             |</head>
             |<body>
             |$svgContent
             |$statsHtml
             |</body>
             |</html>""".stripMargin

        webView.engine.loadContent(html)
        contentBox.children = Seq(webView)

        this.visible = true
        this.managed = true

      case None =>
        contentBox.children = Seq(noDataLabel)
        this.visible = false
        this.managed = false
    }
  }

  /**
   * Clear the ideogram display.
   */
  def clear(): Unit = {
    contentBox.children = Seq(noDataLabel)
    this.visible = false
    this.managed = false
  }
}

object YChromosomeIdeogramPanel {
  def apply(): YChromosomeIdeogramPanel = new YChromosomeIdeogramPanel()
}
