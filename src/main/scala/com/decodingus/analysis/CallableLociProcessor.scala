package com.decodingus.analysis

import com.decodingus.model.{ContigSummary, CoverageSummary}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import org.broadinstitute.hellbender.Main

import java.io.File
import java.nio.file.Files
import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.util.Using

class CallableLociProcessor {

  // Constants
  private val MAX_SVG_WIDTH = 2000
  private val SVG_HEIGHT_PER_CONTIG = 50
  private val BAR_HEIGHT = 200
  private val MARGIN_TOP = 20
  private val MARGIN_BOTTOM = 40
  private val MARGIN_LEFT = 30
  private val MARGIN_RIGHT = 30
  private val TOTAL_FIXED_HEIGHT = BAR_HEIGHT + SVG_HEIGHT_PER_CONTIG + MARGIN_BOTTOM
  private val STRIDE_LEN = 10000
  private val OUTPUT_DIR_NAME = "callable_loci"

  // Colors
  private val COLOR_GREEN = "#007700"
  private val COLOR_RED = "#770000"
  private val TEXT_COLOR = "#CCCCCC"
  private val BG_COLOR = "#222222"
  private val AXIS_COLOR = "#CCCCCC"
  private val TICK_COLOR = "#FFFF00"

  def process(bamPath: String, referencePath: String, onProgress: (String, Int, Int) => Unit): (CoverageSummary, List[String]) = {
    val referenceFile = new File(referencePath)
    val dictionary = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile).getSequenceDictionary
    val contigs = dictionary.getSequences.toArray.map(_.asInstanceOf[htsjdk.samtools.SAMSequenceRecord])
    val totalContigs = contigs.length
    val contigLengths = contigs.map(s => s.getSequenceName -> s.getSequenceLength).toMap
    val maxGenomeLength = if (contigLengths.values.isEmpty) 0 else contigLengths.values.max

    val outputDir = new File(OUTPUT_DIR_NAME)
    if (!outputDir.exists()) {
      Files.createDirectories(outputDir.toPath)
    }

    val allSvgStrings = ListBuffer[String]()
    val allContigSummaries = ListBuffer[ContigSummary]()

    for ((contig, index) <- contigs.zipWithIndex) {
      val contigName = contig.getSequenceName
      val contigLength = contig.getSequenceLength

      onProgress(s"Analyzing contig: $contigName (${index + 1} of $totalContigs)", index + 1, totalContigs)

      val bedFile = new File(outputDir, s"$contigName.callable.bed")
      val summaryFile = new File(outputDir, s"$contigName.table.txt")

      val args = Array(
        "CallableLoci",
        "-I", bamPath,
        "-R", referencePath,
        "-O", bedFile.getAbsolutePath,
        "--summary", summaryFile.getAbsolutePath,
        "-L", contigName
      )
      Main.main(args)

      val binData = binIntervals(bedFile.getAbsolutePath, contigName, contigLength)
      val svgString = generateSvg(contigName, contigLength, maxGenomeLength, binData)
      allSvgStrings += svgString

      val contigSummary = parseSummary(summaryFile.getAbsolutePath, contigName)
      allContigSummaries += contigSummary
    }

    val totalBases = allContigSummaries.map(s => s.callable + s.noCoverage + s.lowCoverage + s.excessiveCoverage + s.poorMappingQuality + s.refN).sum
    val callableBases = allContigSummaries.map(_.callable).sum

    val coverageSummary = CoverageSummary(
      pdsUserId = "60820188481374",
      platformSource = "bwa-mem2",
      reference = "T2T-CHM13v2.0",
      totalReads = 21206068,
      readLength = 147,
      totalBases = totalBases,
      callableBases = callableBases,
      averageDepth = if (totalBases > 0) (21206068L * 147) / totalBases.toDouble else 0.0,
      contigAnalysis = allContigSummaries.toList
    )

    (coverageSummary, allSvgStrings.toList)
  }

  private def binIntervals(bedPath: String, contigName: String, contigLength: Int): Array[Array[Int]] = {
    val maxBin = (contigLength.toDouble / STRIDE_LEN).ceil.toInt
    val binData = Array.fill(maxBin)(Array.fill(3)(0)) // 0: CALLABLE, 1: POOR_MAPPING_QUALITY, 2: Other

    Using(Source.fromFile(bedPath)) { source =>
      for (line <- source.getLines()) {
        if (!line.startsWith("#") && line.trim.nonEmpty) {
          val fields = line.split("\\s+")
          if (fields.length >= 4 && fields(0) == contigName) {
            val start = fields(1).toInt
            val stop = fields(2).toInt
            val status = fields(3)

            (start until stop).foreach {
              basePos =>
                val binIndex = basePos / STRIDE_LEN
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

  private def generateSvg(contigName: String, contigLength: Int, maxGenomeLength: Int, binData: Array[Array[Int]]): String = {
    val scalingFactor = contigLength.toDouble / maxGenomeLength
    val currentSvgWidth = (MAX_SVG_WIDTH * scalingFactor).max(50) + MARGIN_LEFT + MARGIN_RIGHT
    val drawableWidth = currentSvgWidth - MARGIN_LEFT - MARGIN_RIGHT
    val maxBin = (contigLength.toDouble / STRIDE_LEN).ceil.toInt
    val pixelsPerBin = drawableWidth / maxBin

    val svg = new StringBuilder
    svg.append(s"""<svg width="${currentSvgWidth.round}" height="$TOTAL_FIXED_HEIGHT" viewBox="0 0 ${currentSvgWidth.round} $TOTAL_FIXED_HEIGHT" xmlns="http://www.w3.org/2000/svg" font-family="Arial, sans-serif">
    <rect x="0" y="0" width="${currentSvgWidth.round}" height="$TOTAL_FIXED_HEIGHT" fill="$BG_COLOR" />
    <text x="${currentSvgWidth / 2}" y="${MARGIN_TOP + 15}" text-anchor="middle" font-size="20" fill="$TEXT_COLOR">$contigName (Stride: ${STRIDE_LEN / 1000}kb)</text>
  """)

    val drawYOffset = MARGIN_TOP + SVG_HEIGHT_PER_CONTIG

    binData.zipWithIndex.foreach { case (counts, index) =>
      val binXStart = MARGIN_LEFT + (index * pixelsPerBin)
      val callableDepth = counts(0).toDouble / STRIDE_LEN
      val poorQualDepth = counts(1).toDouble / STRIDE_LEN
      val otherDepth = counts(2).toDouble / STRIDE_LEN

      var yPos = drawYOffset + BAR_HEIGHT
      if (callableDepth > 0) {
        val heightPx = (callableDepth * BAR_HEIGHT).round.toInt
        yPos -= heightPx
        svg.append(s"  <rect x=\"$binXStart\" y=\"$yPos\" width=\"$pixelsPerBin\" height=\"$heightPx\" fill=\"$COLOR_GREEN\" />")
      }
      if (poorQualDepth > 0) {
        val heightPx = (poorQualDepth * BAR_HEIGHT).round.toInt
        yPos -= heightPx
        svg.append(s"  <rect x=\"$binXStart\" y=\"$yPos\" width=\"$pixelsPerBin\" height=\"$heightPx\" fill=\"$COLOR_RED\" />")
      }
      if (otherDepth > 0) {
        val heightPx = (otherDepth * BAR_HEIGHT).round.toInt
        yPos -= heightPx
        svg.append(s"  <rect x=\"$binXStart\" y=\"$yPos\" width=\"$pixelsPerBin\" height=\"$heightPx\" fill=\"#AAAAAA\" />")
      }
    }

    svg.append(s"  <line x1=\"$MARGIN_LEFT\" y1=\"${drawYOffset + BAR_HEIGHT}\" x2=\"${MARGIN_LEFT + drawableWidth}\" y2=\"${drawYOffset + BAR_HEIGHT}\" stroke=\"$AXIS_COLOR\" stroke-width=\"1\" />")

    val textY = drawYOffset + BAR_HEIGHT + 15
    val tickYTop = drawYOffset + BAR_HEIGHT - 2
    val tickYBottom = drawYOffset + BAR_HEIGHT + 3

    (10000000 to contigLength by 10000000).foreach {
      mbMark =>
        val markX = MARGIN_LEFT + (mbMark.toDouble / contigLength * drawableWidth)
        svg.append(s"  <line x1=\"$markX\" y1=\"$tickYTop\" x2=\"$markX\" y2=\"$tickYBottom\" stroke=\"$TICK_COLOR\" stroke-width=\"2\" />")
        svg.append(s"  <text x=\"$markX\" y=\"$textY\" text-anchor=\"middle\" font-size=\"12\" fill=\"$TICK_COLOR\">${mbMark / 1000000}Mb</text>")
    }

    svg.append("</svg>")
    svg.toString()
  }

  private def parseSummary(summaryPath: String, contigName: String): ContigSummary = {
    val summaryMap = scala.collection.mutable.Map[String, Long]()
    Using(Source.fromFile(summaryPath)) { source =>
      for (line <- source.getLines()) {
        if (!line.strip.startsWith("state nBases") && line.strip.nonEmpty) {
          val fields = line.strip.split("\\s+")
          if (fields.length == 2) {
            summaryMap(fields(0)) = fields(1).toLong
          }
        }
      }
    }

    ContigSummary(
      contigName = contigName,
      refN = summaryMap.getOrElse("REF_N", 0L),
      callable = summaryMap.getOrElse("CALLABLE", 0L),
      noCoverage = summaryMap.getOrElse("NO_COVERAGE", 0L),
      lowCoverage = summaryMap.getOrElse("LOW_COVERAGE", 0L),
      excessiveCoverage = summaryMap.getOrElse("EXCESSIVE_COVERAGE", 0L),
      poorMappingQuality = summaryMap.getOrElse("POOR_MAPPING_QUALITY", 0L)
    )
  }
}