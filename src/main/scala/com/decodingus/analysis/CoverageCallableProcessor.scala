package com.decodingus.analysis

import com.decodingus.model.ContigSummary

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Path}
import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.util.Using

/**
 * Processor for collecting unified coverage metrics and callable loci from BAM/CRAM files.
 *
 * Replaces both WgsMetricsProcessor and CallableLociProcessor with a single-pass
 * position-level analysis using HTSJDK SamLocusIterator.
 *
 * Features:
 * - Single pass through BAM/CRAM for both coverage and callable loci
 * - Per-contig BED files for callable regions (required for Y-chromosome analysis)
 * - Per-contig summary files in GATK format
 * - Compatible with existing ContigSummary visualization
 */
class CoverageCallableProcessor {

  private val ARTIFACT_SUBDIR_NAME = "callable_loci"
  private val walker = new CoverageCallableWalker()

  /**
   * Process a BAM/CRAM file to collect coverage and callable loci metrics.
   *
   * @param bamPath         Path to BAM/CRAM file
   * @param referencePath   Path to reference genome
   * @param onProgress      Progress callback (message, current, total)
   * @param artifactContext Optional context for organizing output artifacts
   * @param minDepth        Minimum depth to consider callable (default 4, use 2 for HiFi)
   * @return Either error or tuple of (CoverageCallableResult, List of SVG strings)
   */
  def process(
    bamPath: String,
    referencePath: String,
    onProgress: (String, Double, Double) => Unit,
    artifactContext: Option[ArtifactContext] = None,
    minDepth: Int = 4
  ): Either[Throwable, (CoverageCallableResult, List[String])] = {

    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0.0, 1.0)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(new RuntimeException(error))
      case Right(_) => // index exists or was created
    }

    // Verify input file exists
    if (!new File(bamPath).exists()) {
      return Left(new RuntimeException(s"BAM/CRAM file not found: $bamPath"))
    }

    // Use artifact cache directory if context provided, otherwise use local directory
    val outputDir: Path = artifactContext match {
      case Some(ctx) => ctx.getSubdir(ARTIFACT_SUBDIR_NAME)
      case None =>
        val dir = new File(ARTIFACT_SUBDIR_NAME).toPath
        if (!Files.exists(dir)) Files.createDirectories(dir)
        dir
    }

    // Progress adapter
    val progressAdapter: (String, Long, Long) => Unit = (msg, current, total) => {
      val fraction = if (total > 0) current.toDouble / total else 0.0
      onProgress(msg, fraction * 0.9, 1.0) // Reserve last 10% for SVG generation
    }

    // Configure callable params with minDepth
    val callableParams = CallableLociParams(minDepth = minDepth)

    walker.collectCoverageAndCallable(bamPath, referencePath, outputDir, callableParams, progressAdapter) match {
      case Right(result) =>
        onProgress("Generating visualizations...", 0.9, 1.0)

        // Generate SVG visualizations for each contig (reusing existing logic pattern)
        val svgStrings = generateSvgVisualizations(outputDir, result.contigSummaries)

        onProgress("Coverage and callable analysis complete.", 1.0, 1.0)
        Right((result, svgStrings))

      case Left(error) =>
        Left(new RuntimeException(error))
    }
  }

  /**
   * Generate SVG visualizations for callable loci by reading the BED files.
   * Reuses the visualization approach from CallableLociProcessor.
   */
  private def generateSvgVisualizations(outputDir: Path, contigSummaries: List[ContigSummary]): List[String] = {
    val STRIDE_LEN = 10000
    val MAX_SVG_WIDTH = 2000
    val SVG_HEIGHT_PER_CONTIG = 50
    val BAR_HEIGHT = 200
    val MARGIN_TOP = 20
    val MARGIN_BOTTOM = 40
    val MARGIN_LEFT = 30
    val MARGIN_RIGHT = 30
    val TOTAL_FIXED_HEIGHT = BAR_HEIGHT + SVG_HEIGHT_PER_CONTIG + MARGIN_BOTTOM

    val COLOR_GREEN = "#007700"
    val COLOR_RED = "#770000"
    val TEXT_COLOR = "#CCCCCC"
    val BG_COLOR = "#222222"
    val AXIS_COLOR = "#CCCCCC"
    val TICK_COLOR = "#FFFF00"

    // Get max contig length for scaling
    val maxContigLength = contigSummaries.map { s =>
      s.refN + s.callable + s.noCoverage + s.lowCoverage + s.excessiveCoverage + s.poorMappingQuality
    }.maxOption.getOrElse(1L)

    contigSummaries.flatMap { summary =>
      val contigName = summary.contigName
      val bedFile = outputDir.resolve(s"$contigName.callable.bed")

      if (!Files.exists(bedFile)) {
        None
      } else {
        val contigLength = (summary.refN + summary.callable + summary.noCoverage +
          summary.lowCoverage + summary.excessiveCoverage + summary.poorMappingQuality).toInt

        val binData = binIntervalsFromBed(bedFile, contigName, contigLength, STRIDE_LEN)
        val svg = generateSvgForContig(
          contigName, contigLength, maxContigLength.toInt, binData,
          STRIDE_LEN, MAX_SVG_WIDTH, BAR_HEIGHT, MARGIN_TOP, MARGIN_BOTTOM,
          MARGIN_LEFT, MARGIN_RIGHT, SVG_HEIGHT_PER_CONTIG, TOTAL_FIXED_HEIGHT,
          COLOR_GREEN, COLOR_RED, TEXT_COLOR, BG_COLOR, AXIS_COLOR, TICK_COLOR
        )

        // Write SVG to file
        val svgFile = outputDir.resolve(s"$contigName.callable.svg")
        Files.writeString(svgFile, svg)

        Some(svg)
      }
    }
  }

  /**
   * Bin intervals from BED file for visualization.
   */
  private def binIntervalsFromBed(bedPath: Path, contigName: String, contigLength: Int, strideLen: Int): Array[Array[Int]] = {
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

            (start until stop).foreach { basePos =>
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
   * Generate SVG for a single contig.
   */
  private def generateSvgForContig(
    contigName: String,
    contigLength: Int,
    maxGenomeLength: Int,
    binData: Array[Array[Int]],
    strideLen: Int,
    maxSvgWidth: Int,
    barHeight: Int,
    marginTop: Int,
    marginBottom: Int,
    marginLeft: Int,
    marginRight: Int,
    svgHeightPerContig: Int,
    totalFixedHeight: Int,
    colorGreen: String,
    colorRed: String,
    textColor: String,
    bgColor: String,
    axisColor: String,
    tickColor: String
  ): String = {
    val scalingFactor = contigLength.toDouble / maxGenomeLength
    val currentSvgWidth = (maxSvgWidth * scalingFactor).max(50) + marginLeft + marginRight
    val drawableWidth = currentSvgWidth - marginLeft - marginRight
    val maxBin = (contigLength.toDouble / strideLen).ceil.toInt
    val pixelsPerBin = drawableWidth / maxBin

    val svg = new StringBuilder
    svg.append(
      s"""<svg width="${currentSvgWidth.round}" height="$totalFixedHeight" viewBox="0 0 ${currentSvgWidth.round} $totalFixedHeight" xmlns="http://www.w3.org/2000/svg" font-family="Arial, sans-serif">
    <rect x="0" y="0" width="${currentSvgWidth.round}" height="$totalFixedHeight" fill="$bgColor" />
    <text x="${currentSvgWidth / 2}" y="${marginTop + 15}" text-anchor="middle" font-size="20" fill="$textColor">$contigName (Stride: ${strideLen / 1000}kb)</text>
  """)

    val drawYOffset = marginTop + svgHeightPerContig

    binData.zipWithIndex.foreach { case (counts, index) =>
      val binXStart = marginLeft + (index * pixelsPerBin)
      val callableDepth = counts(0).toDouble / strideLen
      val poorQualDepth = counts(1).toDouble / strideLen
      val otherDepth = counts(2).toDouble / strideLen

      var yPos = drawYOffset + barHeight
      if (callableDepth > 0) {
        val heightPx = (callableDepth * barHeight).round.toInt
        yPos -= heightPx
        svg.append(s"""  <rect x="$binXStart" y="$yPos" width="$pixelsPerBin" height="$heightPx" fill="$colorGreen" />""")
      }
      if (poorQualDepth > 0) {
        val heightPx = (poorQualDepth * barHeight).round.toInt
        yPos -= heightPx
        svg.append(s"""  <rect x="$binXStart" y="$yPos" width="$pixelsPerBin" height="$heightPx" fill="$colorRed" />""")
      }
      if (otherDepth > 0) {
        val heightPx = (otherDepth * barHeight).round.toInt
        yPos -= heightPx
        svg.append(s"""  <rect x="$binXStart" y="$yPos" width="$pixelsPerBin" height="$heightPx" fill="#AAAAAA" />""")
      }
    }

    svg.append(s"""  <line x1="$marginLeft" y1="${drawYOffset + barHeight}" x2="${marginLeft + drawableWidth}" y2="${drawYOffset + barHeight}" stroke="$axisColor" stroke-width="1" />""")

    val textY = drawYOffset + barHeight + 15
    val tickYTop = drawYOffset + barHeight - 2
    val tickYBottom = drawYOffset + barHeight + 3

    (10000000 to contigLength by 10000000).foreach { mbMark =>
      val markX = marginLeft + (mbMark.toDouble / contigLength * drawableWidth)
      svg.append(s"""  <line x1="$markX" y1="$tickYTop" x2="$markX" y2="$tickYBottom" stroke="$tickColor" stroke-width="2" />""")
      svg.append(s"""  <text x="$markX" y="$textY" text-anchor="middle" font-size="12" fill="$tickColor">${mbMark / 1000000}Mb</text>""")
    }

    svg.append("</svg>")
    svg.toString()
  }
}

object CoverageCallableProcessor {

  /**
   * Load CoverageCallableResult from cached artifacts.
   * Reads the .table.txt summary files and coverage data from the callable_loci directory.
   *
   * @param callableLociDir Path to the callable_loci artifact directory
   * @return CoverageCallableResult if successful, None if not found or invalid
   */
  def loadFromCache(callableLociDir: Path): Option[CoverageCallableResult] = {
    if (!Files.exists(callableLociDir)) return None

    import scala.jdk.CollectionConverters.*

    val tableFiles = Files.list(callableLociDir).iterator().asScala
      .filter(_.toString.endsWith(".table.txt"))
      .toList

    if (tableFiles.isEmpty) return None

    val contigSummaries = ListBuffer[ContigSummary]()

    for (tableFile <- tableFiles) {
      val fileName = tableFile.getFileName.toString
      val contigName = fileName.stripSuffix(".table.txt")

      Using(Source.fromFile(tableFile.toFile)) { source =>
        val summaryMap = scala.collection.mutable.Map[String, Long]()
        for (line <- source.getLines()) {
          if (!line.strip.startsWith("state nBases") && line.strip.nonEmpty) {
            val fields = line.strip.split("\\s+")
            if (fields.length == 2) {
              summaryMap(fields(0)) = fields(1).toLong
            }
          }
        }

        contigSummaries += ContigSummary(
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

    if (contigSummaries.isEmpty) return None

    // Sort contigs by standard order
    val sortedSummaries = contigSummaries.toList.sortBy { cs =>
      val name = cs.contigName.replaceFirst("^chr", "")
      name match {
        case "X" => 23
        case "Y" => 24
        case "M" | "MT" => 25
        case n if n.forall(_.isDigit) => n.toInt
        case _ => 100
      }
    }

    val callableBases = sortedSummaries.map(_.callable).sum

    // Note: When loading from cache, we don't have the full coverage metrics
    // This is a partial result suitable for callable loci visualization
    Some(CoverageCallableResult(
      genomeTerritory = sortedSummaries.map(s => s.refN + s.callable + s.noCoverage + s.lowCoverage + s.excessiveCoverage + s.poorMappingQuality).sum,
      meanCoverage = 0.0,  // Not available from cache
      medianCoverage = 0.0,
      sdCoverage = 0.0,
      coverageHistogram = new Array[Long](256),
      pct1x = 0.0,
      pct5x = 0.0,
      pct10x = 0.0,
      pct15x = 0.0,
      pct20x = 0.0,
      pct25x = 0.0,
      pct30x = 0.0,
      pct40x = 0.0,
      pct50x = 0.0,
      callableBases = callableBases,
      contigSummaries = sortedSummaries,
      contigCoverage = Map.empty
    ))
  }

  /**
   * Check if callable loci data exists in cache.
   */
  def existsInCache(callableLociDir: Path): Boolean = {
    if (!Files.exists(callableLociDir)) return false
    import scala.jdk.CollectionConverters.*
    Files.list(callableLociDir).iterator().asScala.exists(_.toString.endsWith(".table.txt"))
  }
}
