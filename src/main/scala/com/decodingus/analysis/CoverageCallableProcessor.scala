package com.decodingus.analysis

import com.decodingus.analysis.util.BioVisualizationUtil
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

        val binData = BioVisualizationUtil.binIntervalsFromBed(bedFile, contigName, contigLength)
        val svg = BioVisualizationUtil.generateSvgForContig(
          contigName, contigLength, maxContigLength.toInt, binData
        )

        // Write SVG to file
        val svgFile = outputDir.resolve(s"$contigName.callable.svg")
        Files.writeString(svgFile, svg)

        Some(svg)
      }
    }
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
    // Try to load coverage.txt if available for samtools-style stats
    val coverageStats = loadCoverageStats(callableLociDir)

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
      contigCoverage = Map.empty,
      contigCoverageStats = coverageStats
    ))
  }

  /**
   * Load samtools-style coverage stats from coverage.txt if available.
   */
  private def loadCoverageStats(callableLociDir: Path): List[ContigCoverageStats] = {
    val coverageFile = callableLociDir.resolve("coverage.txt")
    if (!Files.exists(coverageFile)) return List.empty

    val stats = ListBuffer[ContigCoverageStats]()
    Using(Source.fromFile(coverageFile.toFile)) { source =>
      for (line <- source.getLines()) {
        if (!line.startsWith("#") && line.trim.nonEmpty) {
          val fields = line.split("\\t")
          if (fields.length >= 9) {
            stats += ContigCoverageStats(
              contig = fields(0),
              startPos = fields(1).toLong,
              endPos = fields(2).toLong,
              numReads = fields(3).toLong,
              covBases = fields(4).toLong,
              coverage = fields(5).toDouble,
              meanDepth = fields(6).toDouble,
              meanBaseQ = fields(7).toDouble,
              meanMapQ = fields(8).toDouble
            )
          }
        }
      }
    }
    stats.toList
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
