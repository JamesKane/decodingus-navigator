package com.decodingus.ui.components

import com.decodingus.refgenome.{GenomicRegion, RegionType, YRegionAnnotator}
import com.decodingus.yprofile.model.{YConsensusState, YProfileVariantEntity, YVariantStatus}

/**
 * Renders Y chromosome ideogram as SVG with region bands and variant markers.
 *
 * Visual layout:
 * - Triangle markers above chromosome for derived variants
 * - Rounded rectangle chromosome with color-coded region bands
 * - Scale axis below with Mbp tick marks
 * - Legend showing region colors and variant marker colors
 */
object YChromosomeIdeogramRenderer {

  // SVG dimensions
  val SVG_WIDTH = 900
  val SVG_HEIGHT = 240
  val IDEOGRAM_HEIGHT = 40
  val IDEOGRAM_Y = 70
  val MARKER_Y = IDEOGRAM_Y - 15
  val SCALE_Y = IDEOGRAM_Y + IDEOGRAM_HEIGHT + 20
  val LEGEND_Y = SCALE_Y + 35
  val MARGIN_X = 40

  // Region colors (dark theme compatible)
  private val regionColors: Map[RegionType, String] = Map(
    RegionType.PAR -> "#6B8E23",           // Olive green
    RegionType.XDegenerate -> "#228B22",   // Forest green
    RegionType.XTR -> "#CD853F",           // Tan
    RegionType.Ampliconic -> "#DAA520",    // Goldenrod
    RegionType.Palindrome -> "#FF8C00",    // Dark orange
    RegionType.Heterochromatin -> "#4A4A4A", // Dark gray
    RegionType.Centromere -> "#696969",    // Dim gray
    RegionType.STR -> "#9370DB"            // Medium purple (for STR markers if shown)
  )

  // Variant marker colors by status
  private val variantColors: Map[YVariantStatus, String] = Map(
    YVariantStatus.CONFIRMED -> "#4CAF50", // Green
    YVariantStatus.NOVEL -> "#2196F3",     // Blue
    YVariantStatus.CONFLICT -> "#F44336",  // Red
    YVariantStatus.PENDING -> "#FF9800"    // Orange
  )

  /**
   * Variant marker data for rendering.
   *
   * @param position Genomic position
   * @param status   Variant status (determines color)
   * @param label    Optional label for tooltip (variant name)
   */
  case class VariantMarker(
    position: Long,
    status: YVariantStatus,
    label: Option[String] = None
  )

  object VariantMarker {
    def fromVariantEntity(v: YProfileVariantEntity): VariantMarker =
      VariantMarker(v.position, v.status, v.variantName.orElse(v.canonicalName))
  }

  /**
   * Render the ideogram as SVG.
   *
   * @param annotator     Y region annotator with region data
   * @param variants      Optional list of variant markers to display
   * @param showAllRegions Whether to show all region types or just major ones
   * @return SVG string
   */
  def render(
    annotator: YRegionAnnotator,
    variants: List[VariantMarker] = Nil,
    showAllRegions: Boolean = true
  ): String = {
    val chromLength = annotator.getChromosomeLength
    val regions = annotator.getAllRegions
    val drawWidth = SVG_WIDTH - (2 * MARGIN_X)

    def posToX(pos: Long): Double = MARGIN_X + (pos.toDouble / chromLength * drawWidth)

    val sb = new StringBuilder
    sb.append(s"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $SVG_WIDTH $SVG_HEIGHT" style="background:#1a1a1a;">""")

    // Title
    sb.append(s"""<text x="${SVG_WIDTH / 2}" y="25" text-anchor="middle" fill="#e0e0e0" font-size="16" font-weight="bold">Y Chromosome Regions</text>""")

    // Draw chromosome base (dark background)
    val chromRadius = IDEOGRAM_HEIGHT / 2
    sb.append(s"""<rect x="$MARGIN_X" y="$IDEOGRAM_Y" width="$drawWidth" height="$IDEOGRAM_HEIGHT" rx="$chromRadius" ry="$chromRadius" fill="#333333"/>""")

    // Draw regions (with clipping to chromosome shape)
    sb.append(s"""<defs><clipPath id="chromClip"><rect x="$MARGIN_X" y="$IDEOGRAM_Y" width="$drawWidth" height="$IDEOGRAM_HEIGHT" rx="$chromRadius" ry="$chromRadius"/></clipPath></defs>""")
    sb.append("""<g clip-path="url(#chromClip)">""")

    // Layer order: X-degenerate first (background), then other regions on top
    val layerOrder = List(
      RegionType.XDegenerate,
      RegionType.Heterochromatin,
      RegionType.Centromere,
      RegionType.PAR,
      RegionType.XTR,
      RegionType.Ampliconic,
      RegionType.Palindrome
    )

    for {
      regionType <- layerOrder
      regionList <- regions.get(regionType)
      region <- regionList
    } {
      val x = posToX(region.start)
      val width = math.max(1, posToX(region.end) - x)
      val color = regionColors.getOrElse(regionType, "#666666")
      val tooltip = region.name.map(n => s"$n (${regionType.displayName})").getOrElse(regionType.displayName)
      sb.append(s"""<rect x="$x" y="$IDEOGRAM_Y" width="$width" height="$IDEOGRAM_HEIGHT" fill="$color"><title>$tooltip: ${formatPosition(region.start)} - ${formatPosition(region.end)}</title></rect>""")
    }

    sb.append("</g>")

    // Chromosome outline
    sb.append(s"""<rect x="$MARGIN_X" y="$IDEOGRAM_Y" width="$drawWidth" height="$IDEOGRAM_HEIGHT" rx="$chromRadius" ry="$chromRadius" fill="none" stroke="#666666" stroke-width="1"/>""")

    // Draw variant markers (triangles above chromosome)
    val derivedVariants = variants.filter(v =>
      v.status == YVariantStatus.CONFIRMED ||
      v.status == YVariantStatus.NOVEL ||
      v.status == YVariantStatus.CONFLICT
    )

    // Group nearby variants to avoid overcrowding
    val markerSize = 6
    val minSpacing = 3
    var lastX = Double.MinValue

    for (v <- derivedVariants.sortBy(_.position)) {
      val x = posToX(v.position)
      if (x - lastX >= minSpacing) {
        val color = variantColors.getOrElse(v.status, "#888888")
        val tooltip = v.label.getOrElse(s"pos:${v.position}") + s" (${v.status})"
        // Triangle pointing down
        val points = s"${x},${MARKER_Y + markerSize} ${x - markerSize / 2},$MARKER_Y ${x + markerSize / 2},$MARKER_Y"
        sb.append(s"""<polygon points="$points" fill="$color" stroke="none"><title>$tooltip</title></polygon>""")
        lastX = x
      }
    }

    // Draw scale axis
    renderScaleAxis(sb, chromLength, posToX)

    // Draw legend
    renderLegend(sb, regions.keySet)

    sb.append("</svg>")
    sb.toString()
  }

  private def renderScaleAxis(sb: StringBuilder, chromLength: Long, posToX: Long => Double): Unit = {
    // Axis line
    sb.append(s"""<line x1="$MARGIN_X" y1="$SCALE_Y" x2="${SVG_WIDTH - MARGIN_X}" y2="$SCALE_Y" stroke="#888888" stroke-width="1"/>""")

    // Tick marks every 10 Mbp
    val tickInterval = 10_000_000L
    var pos = 0L
    while (pos <= chromLength) {
      val x = posToX(pos)
      sb.append(s"""<line x1="$x" y1="$SCALE_Y" x2="$x" y2="${SCALE_Y + 5}" stroke="#888888" stroke-width="1"/>""")
      val label = s"${pos / 1_000_000}"
      sb.append(s"""<text x="$x" y="${SCALE_Y + 18}" text-anchor="middle" fill="#888888" font-size="10">${label}M</text>""")
      pos += tickInterval
    }
  }

  private def renderLegend(sb: StringBuilder, presentRegions: Set[RegionType]): Unit = {
    val legendItems = List(
      (RegionType.PAR, "PAR"),
      (RegionType.XDegenerate, "X-deg"),
      (RegionType.XTR, "XTR"),
      (RegionType.Ampliconic, "Ampliconic"),
      (RegionType.Palindrome, "Palindrome"),
      (RegionType.Heterochromatin, "Het")
    ).filter(item => presentRegions.contains(item._1))

    var x = MARGIN_X
    val boxSize = 12
    val spacing = 15

    // Region legend
    for ((regionType, label) <- legendItems) {
      val color = regionColors.getOrElse(regionType, "#666666")
      sb.append(s"""<rect x="$x" y="${LEGEND_Y}" width="$boxSize" height="$boxSize" fill="$color" rx="2"/>""")
      sb.append(s"""<text x="${x + boxSize + 4}" y="${LEGEND_Y + 10}" fill="#cccccc" font-size="10">$label</text>""")
      x += boxSize + label.length * 6 + spacing
    }

    // Variant marker legend (right side)
    x = SVG_WIDTH - 280
    val variantItems = List(
      (YVariantStatus.CONFIRMED, "Confirmed"),
      (YVariantStatus.NOVEL, "Novel"),
      (YVariantStatus.CONFLICT, "Conflict")
    )

    for ((status, label) <- variantItems) {
      val color = variantColors.getOrElse(status, "#888888")
      // Small triangle
      val points = s"${x + 3},${LEGEND_Y + boxSize} ${x},$LEGEND_Y ${x + 6},$LEGEND_Y"
      sb.append(s"""<polygon points="$points" fill="$color"/>""")
      sb.append(s"""<text x="${x + 10}" y="${LEGEND_Y + 10}" fill="#cccccc" font-size="10">$label</text>""")
      x += label.length * 6 + 25
    }
  }

  private def formatPosition(pos: Long): String = {
    if (pos >= 1_000_000) f"${pos / 1_000_000.0}%.2fM"
    else if (pos >= 1_000) f"${pos / 1_000.0}%.1fK"
    else pos.toString
  }

  /**
   * Render a simple stats panel as HTML to show below the ideogram.
   */
  def renderStatsHtml(
    variants: List[YProfileVariantEntity],
    annotator: YRegionAnnotator
  ): String = {
    val derived = variants.count(_.consensusState == YConsensusState.DERIVED)
    val confirmed = variants.count(_.status == YVariantStatus.CONFIRMED)
    val novel = variants.count(_.status == YVariantStatus.NOVEL)
    val conflict = variants.count(_.status == YVariantStatus.CONFLICT)
    val regions = annotator.getAllRegions

    val regionCounts = regions.map { case (rt, rs) => s"${rt.displayName}: ${rs.size}" }.mkString(" | ")

    s"""
       |<div style="margin-top: 15px; padding: 10px; background: #2a2a2a; border-radius: 5px; font-size: 12px; color: #cccccc;">
       |  <span style="color: #4CAF50; font-weight: bold;">Derived: $derived</span> &nbsp;|&nbsp;
       |  <span style="color: #4CAF50;">Confirmed: $confirmed</span> &nbsp;|&nbsp;
       |  <span style="color: #2196F3;">Novel: $novel</span> &nbsp;|&nbsp;
       |  <span style="color: #F44336;">Conflict: $conflict</span>
       |  <div style="margin-top: 8px; font-size: 11px; color: #888888;">Regions: $regionCounts</div>
       |</div>
    """.stripMargin
  }
}
