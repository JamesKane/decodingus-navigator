package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.workspace.model.{IbdSegment, MatchResult, RelationshipEstimate}
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Label
import scalafx.scene.layout.{Priority, VBox}
import scalafx.scene.web.WebView

/**
 * Chromosome browser panel showing autosomal ideograms with highlighted IBD segments.
 *
 * Renders 22 autosomal chromosomes + X as horizontal bars with colored overlays
 * for shared IBD segments. Segment color saturation indicates length (longer = brighter).
 */
class ChromosomeBrowserPanel extends VBox {

  spacing = 10
  padding = Insets(15)
  style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"

  private val titleLabel = new Label(t("ibd.chromosome_browser")) {
    style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }

  private val webView = new WebView {
    prefHeight = 420
    minHeight = 380
  }
  VBox.setVgrow(webView, Priority.Always)

  private val placeholderLabel = new Label(t("ibd.select_match_to_view")) {
    style = "-fx-font-size: 14px; -fx-text-fill: #888888;"
  }

  private val contentBox = new VBox(10) {
    alignment = Pos.Center
    children = Seq(placeholderLabel)
  }
  VBox.setVgrow(contentBox, Priority.Always)

  children = Seq(titleLabel, contentBox)

  /**
   * Display segments for a match result on the chromosome browser.
   */
  def setMatch(matchResult: MatchResult): Unit =
    if matchResult.sharedSegments.isEmpty then
      contentBox.children = Seq(placeholderLabel)
    else
      val html = ChromosomeBrowserRenderer.renderHtml(matchResult.sharedSegments, Some(matchResult))
      webView.engine.loadContent(html)
      contentBox.children = Seq(webView)

  /**
   * Display raw segments on the chromosome browser.
   */
  def setSegments(segments: List[IbdSegment]): Unit =
    if segments.isEmpty then
      contentBox.children = Seq(placeholderLabel)
    else
      val html = ChromosomeBrowserRenderer.renderHtml(segments)
      webView.engine.loadContent(html)
      contentBox.children = Seq(webView)

  /**
   * Clear the browser display.
   */
  def clear(): Unit =
    contentBox.children = Seq(placeholderLabel)
}

/**
 * Renders SVG chromosome ideograms with IBD segment overlays.
 */
object ChromosomeBrowserRenderer:

  // GRCh38 chromosome lengths in base pairs
  private val chromosomeLengths: Map[String, Long] = Map(
    "1" -> 248956422L, "2" -> 242193529L, "3" -> 198295559L,
    "4" -> 190214555L, "5" -> 181538259L, "6" -> 170805979L,
    "7" -> 159345973L, "8" -> 145138636L, "9" -> 138394717L,
    "10" -> 133797422L, "11" -> 135086622L, "12" -> 133275309L,
    "13" -> 114364328L, "14" -> 107043718L, "15" -> 101991189L,
    "16" -> 90338345L, "17" -> 83257441L, "18" -> 80373285L,
    "19" -> 58617616L, "20" -> 64444167L, "21" -> 46709983L,
    "22" -> 50818468L, "X" -> 156040895L
  )

  private val chromosomeOrder: List[String] =
    (1 to 22).map(_.toString).toList :+ "X"

  private val maxChrLength: Long = chromosomeLengths.values.max

  def renderHtml(segments: List[IbdSegment], matchResult: Option[MatchResult] = None): String =
    val segmentsByChr = segments.groupBy(_.chromosome)
    val maxSegmentCm = if segments.nonEmpty then segments.map(_.lengthCm).max else 1.0

    val svgWidth = 700
    val chrHeight = 14
    val chrSpacing = 18
    val labelWidth = 30
    val barWidth = svgWidth - labelWidth - 10
    val svgHeight = chromosomeOrder.size * chrSpacing + 30

    val svgLines = new StringBuilder
    svgLines.append(s"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $svgWidth $svgHeight" width="100%">""")
    svgLines.append("""<style>""")
    svgLines.append(""".chr-label { font: 11px system-ui, sans-serif; fill: #999; }""")
    svgLines.append(""".chr-bar { fill: #3a3a3a; rx: 4; ry: 4; }""")
    svgLines.append(""".segment { rx: 2; ry: 2; cursor: pointer; }""")
    svgLines.append(""".segment:hover { stroke: #ffffff; stroke-width: 1.5; }""")
    svgLines.append(""".tooltip-text { font: 10px monospace; fill: #ffffff; }""")
    svgLines.append("""</style>""")

    for (chr, idx) <- chromosomeOrder.zipWithIndex do
      val y = idx * chrSpacing + 10
      val chrLen = chromosomeLengths.getOrElse(chr, maxChrLength)
      val barScale = barWidth.toDouble / maxChrLength
      val scaledLen = (chrLen * barScale).toInt

      // Chromosome label
      svgLines.append(s"""<text x="${labelWidth - 4}" y="${y + chrHeight - 3}" class="chr-label" text-anchor="end">$chr</text>""")

      // Chromosome bar background
      svgLines.append(s"""<rect x="$labelWidth" y="$y" width="$scaledLen" height="$chrHeight" class="chr-bar"/>""")

      // IBD segment overlays
      segmentsByChr.getOrElse(chr, Nil).foreach { seg =>
        val xStart = labelWidth + (seg.startPosition.toDouble / maxChrLength * barWidth).toInt
        val segWidth = math.max(2, ((seg.endPosition - seg.startPosition).toDouble / maxChrLength * barWidth).toInt)
        val saturation = math.min(1.0, seg.lengthCm / maxSegmentCm)
        val color = segmentColor(saturation)
        val tooltip = segmentTooltip(seg)
        svgLines.append(s"""<rect x="$xStart" y="$y" width="$segWidth" height="$chrHeight" class="segment" fill="$color">""")
        svgLines.append(s"""<title>$tooltip</title>""")
        svgLines.append("""</rect>""")
      }

    svgLines.append("</svg>")

    // Summary stats
    val statsHtml = matchResult.map { mr =>
      val relLabel = mr.relationshipEstimate.map(_.label).getOrElse("Unknown")
      s"""<div class="stats">
         |  <span class="stat">Shared: <b>${f"${mr.totalSharedCm}%.1f"} cM</b></span>
         |  <span class="stat">Segments: <b>${mr.segmentCount}</b></span>
         |  <span class="stat">Longest: <b>${mr.longestSegmentCm.map(cm => f"$cm%.1f").getOrElse("-")} cM</b></span>
         |  <span class="stat">Estimate: <b>$relLabel</b></span>
         |</div>""".stripMargin
    }.getOrElse("")

    s"""<!DOCTYPE html>
       |<html>
       |<head>
       |  <style>
       |    body { margin: 0; padding: 10px; background: #2a2a2a; font-family: system-ui, sans-serif; color: #cccccc; }
       |    svg { max-width: 100%; height: auto; }
       |    .stats { display: flex; gap: 20px; padding: 10px 0; font-size: 12px; color: #aaaaaa; flex-wrap: wrap; }
       |    .stat b { color: #ffffff; }
       |  </style>
       |</head>
       |<body>
       |$statsHtml
       |${svgLines.toString}
       |</body>
       |</html>""".stripMargin

  private def segmentColor(saturation: Double): String =
    // Blue-green gradient from dim (#2a6090) to bright (#30d0a0) based on cM proportion
    val r = (0x2a + (0x30 - 0x2a) * saturation).toInt
    val g = (0x60 + (0xd0 - 0x60) * saturation).toInt
    val b = (0x90 + (0xa0 - 0x90) * saturation).toInt
    f"#$r%02x$g%02x$b%02x"

  private def segmentTooltip(seg: IbdSegment): String =
    val snpInfo = seg.snpCount.map(n => s", $n SNPs").getOrElse("")
    s"Chr ${seg.chromosome}: ${formatBp(seg.startPosition)}–${formatBp(seg.endPosition)} (${f"${seg.lengthCm}%.2f"} cM$snpInfo)"

  def formatBp(bp: Long): String =
    if bp >= 1_000_000 then f"${bp / 1_000_000.0}%.1f Mb"
    else if bp >= 1_000 then f"${bp / 1_000.0}%.1f kb"
    else s"${bp} bp"
