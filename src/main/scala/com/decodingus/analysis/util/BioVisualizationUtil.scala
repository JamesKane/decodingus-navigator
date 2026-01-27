package com.decodingus.analysis.util

import java.nio.file.Path
import scala.io.Source
import scala.util.Using

/**
 * Shared utility for generating bioinformatics visualizations, such as
 * callable loci coverage maps.
 */
object BioVisualizationUtil {

  // Default Constants
  val DEFAULT_MAX_SVG_WIDTH = 2000
  val DEFAULT_SVG_HEIGHT_PER_CONTIG = 50
  val DEFAULT_BAR_HEIGHT = 200
  val DEFAULT_MARGIN_TOP = 20
  val DEFAULT_MARGIN_BOTTOM = 40
  val DEFAULT_MARGIN_LEFT = 30
  val DEFAULT_MARGIN_RIGHT = 30
  val DEFAULT_STRIDE_LEN = 10000

  // Default Colors
  val COLOR_GREEN = "#007700"
  val COLOR_RED = "#770000"
  val COLOR_GREY = "#AAAAAA"
  val TEXT_COLOR = "#CCCCCC"
  val BG_COLOR = "#222222"
  val AXIS_COLOR = "#CCCCCC"
  val TICK_COLOR = "#FFFF00"

  case class SvgConfig(
                        maxSvgWidth: Int = DEFAULT_MAX_SVG_WIDTH,
                        svgHeightPerContig: Int = DEFAULT_SVG_HEIGHT_PER_CONTIG,
                        barHeight: Int = DEFAULT_BAR_HEIGHT,
                        marginTop: Int = DEFAULT_MARGIN_TOP,
                        marginBottom: Int = DEFAULT_MARGIN_BOTTOM,
                        marginLeft: Int = DEFAULT_MARGIN_LEFT,
                        marginRight: Int = DEFAULT_MARGIN_RIGHT,
                        strideLen: Int = DEFAULT_STRIDE_LEN,
                        colorGreen: String = COLOR_GREEN,
                        colorRed: String = COLOR_RED,
                        colorGrey: String = COLOR_GREY,
                        textColor: String = TEXT_COLOR,
                        bgColor: String = BG_COLOR,
                        axisColor: String = AXIS_COLOR,
                        tickColor: String = TICK_COLOR
                      ) {
    def totalFixedHeight: Int = barHeight + svgHeightPerContig + marginBottom
  }

  val defaultSvgConfig: SvgConfig = SvgConfig()

  /**
   * Bin intervals from a BED file for visualization.
   *
   * @param bedPath Path to the BED file.
   * @param contigName The specific contig to process.
   * @param contigLength Length of the contig.
   * @param strideLen Length of each bin (default 10kb).
   * @return Array of bins, where each bin is [CallableCount, PoorQualCount, OtherCount].
   */
  def binIntervalsFromBed(bedPath: Path, contigName: String, contigLength: Int, strideLen: Int = DEFAULT_STRIDE_LEN): Array[Array[Int]] = {
    val maxBin = (contigLength.toDouble / strideLen).ceil.toInt
    val binData = Array.fill(maxBin)(Array.fill(3)(0)) // 0: CALLABLE, 1: POOR_MAPPING_QUALITY, 2: Other

    Using(Source.fromFile(bedPath.toFile)) { source =>
      for (line <- source.getLines()) {
        if (!line.startsWith("#") && line.trim.nonEmpty) {
          val fields = line.split("\\s+")
          if (fields.length >= 4 && fields(0) == contigName) {
            val start = fields(1).toInt
            val stop = fields(2).toInt
            val status = fields(3)

            (start until stop).foreach {
              basePos =>
                val binIndex = basePos / strideLen
                if (binIndex < maxBin) {
                  status match {
                    case "CALLABLE" => binData(binIndex)(0) += 1
                    case "POOR_MAPPING_QUALITY" => binData(binIndex)(1) += 1
                    case _ => binData(binIndex)(2) += 1
                  }
                }
            }
          }
        }
      }
    }
    binData
  }

  /**
   * Generate an SVG string visualization for a contig's callable loci.
   *
   * @param contigName Name of the contig.
   * @param contigLength Length of the contig.
   * @param maxGenomeLength Length of the largest contig (for scaling).
   * @param binData Binned data from `binIntervalsFromBed`.
   * @param config SVG configuration (dimensions, colors).
   * @return SVG content as a String.
   */
  def generateSvgForContig(
                            contigName: String,
                            contigLength: Int,
                            maxGenomeLength: Int,
                            binData: Array[Array[Int]],
                            config: SvgConfig = defaultSvgConfig
                          ): String = {
    val scalingFactor = contigLength.toDouble / maxGenomeLength
    val currentSvgWidth = (config.maxSvgWidth * scalingFactor).max(50) + config.marginLeft + config.marginRight
    val drawableWidth = currentSvgWidth - config.marginLeft - config.marginRight
    val maxBin = (contigLength.toDouble / config.strideLen).ceil.toInt
    val pixelsPerBin = drawableWidth / maxBin

    val svg = new StringBuilder
    svg.append(
      s"""<svg width=\"${currentSvgWidth.round}\" height=\"${config.totalFixedHeight}\" viewBox=\"0 0 ${currentSvgWidth.round} ${config.totalFixedHeight}\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Arial, sans-serif\">
    <rect x=\"0\" y=\"0\" width=\"${currentSvgWidth.round}\" height=\"${config.totalFixedHeight}\" fill=\"${config.bgColor}\" />
    <text x=\"${currentSvgWidth / 2}\" y=\"${config.marginTop + 15}\" text-anchor=\"middle\" font-size=\"20\" fill=\"${config.textColor}\">$contigName (Stride: ${config.strideLen / 1000}kb)</text>
  """
    )

    val drawYOffset = config.marginTop + config.svgHeightPerContig

    binData.zipWithIndex.foreach { case (counts, index) =>
      val binXStart = config.marginLeft + (index * pixelsPerBin)
      val callableDepth = counts(0).toDouble / config.strideLen
      val poorQualDepth = counts(1).toDouble / config.strideLen
      val otherDepth = counts(2).toDouble / config.strideLen

      var yPos = drawYOffset + config.barHeight
      if (callableDepth > 0) {
        val heightPx = (callableDepth * config.barHeight).round.toInt
        yPos -= heightPx
        svg.append(s"  <rect x=\"$binXStart\" y=\"$yPos\" width=\"$pixelsPerBin\" height=\"$heightPx\" fill=\"${config.colorGreen}\" />")
      }
      if (poorQualDepth > 0) {
        val heightPx = (poorQualDepth * config.barHeight).round.toInt
        yPos -= heightPx
        svg.append(s"  <rect x=\"$binXStart\" y=\"$yPos\" width=\"$pixelsPerBin\" height=\"$heightPx\" fill=\"${config.colorRed}\" />")
      }
      if (otherDepth > 0) {
        val heightPx = (otherDepth * config.barHeight).round.toInt
        yPos -= heightPx
        svg.append(s"  <rect x=\"$binXStart\" y=\"$yPos\" width=\"$pixelsPerBin\" height=\"$heightPx\" fill=\"${config.colorGrey}\" />")
      }
    }

    svg.append(s"  <line x1=\"${config.marginLeft}\" y1=\"${drawYOffset + config.barHeight}\" x2=\"${config.marginLeft + drawableWidth}\" y2=\"${drawYOffset + config.barHeight}\" stroke=\"${config.axisColor}\" stroke-width=\"1\" />")

    val textY = drawYOffset + config.barHeight + 15
    val tickYTop = drawYOffset + config.barHeight - 2
    val tickYBottom = drawYOffset + config.barHeight + 3

    (10000000 to contigLength by 10000000).foreach {
      mbMark =>
        val markX = config.marginLeft + (mbMark.toDouble / contigLength * drawableWidth)
        svg.append(s"  <line x1=\"$markX\" y1=\"$tickYTop\" x2=\"$markX\" y2=\"$tickYBottom\" stroke=\"${config.tickColor}\" stroke-width=\"2\" />")
        svg.append(s"  <text x=\"$markX\" y=\"$textY\" text-anchor=\"middle\" font-size=\"12\" fill=\"${config.tickColor}\">${mbMark / 1000000}Mb</text>")
    }

    svg.append("</svg>")
    svg.toString()
  }
}
