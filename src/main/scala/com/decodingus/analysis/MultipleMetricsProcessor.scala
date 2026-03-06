package com.decodingus.analysis

import java.io.File
import java.nio.file.Path
import scala.io.Source
import scala.util.Using

/**
 * Results from CollectMultipleMetrics containing alignment summary and insert size metrics.
 *
 * @param totalReads              Total number of reads (TOTAL_READS from alignment_summary_metrics)
 * @param pfReads                 Pass-filter reads (PF_READS)
 * @param pfReadsAligned          PF reads that aligned (PF_READS_ALIGNED)
 * @param pctPfReadsAligned       Percentage of PF reads aligned
 * @param pfHqAlignedReads        PF reads with high-quality alignments (PF_HQ_ALIGNED_READS)
 * @param readsPaired             Reads from paired-end libraries (READS_ALIGNED_IN_PAIRS)
 * @param pctReadsPaired          Percentage of aligned reads that are paired
 * @param pctProperPairs          Percentage of reads aligned as proper pairs (PF_READS_IMPROPER_PAIRS inverse)
 * @param meanReadLength          Mean read length (MEAN_READ_LENGTH)
 * @param strandBalance           Strand balance (PCT_PF_READS_ALIGNED / strand ratio)
 * @param pctChimeras             Percentage of reads that are chimeric (PCT_CHIMERAS)
 * @param pctAdapter              Percentage of bases that are adapter (PCT_ADAPTER)
 * @param medianInsertSize        Median insert size (MEDIAN_INSERT_SIZE from insert_size_metrics)
 * @param meanInsertSize          Mean insert size (MEAN_INSERT_SIZE)
 * @param stdInsertSize           Standard deviation of insert size (STANDARD_DEVIATION)
 * @param medianAbsoluteDeviation Median absolute deviation of insert size (MEDIAN_ABSOLUTE_DEVIATION)
 * @param minInsertSize           Minimum insert size
 * @param maxInsertSize           Maximum insert size
 * @param pairOrientation         Read pair orientation (FR, RF, TANDEM)
 */
case class MultipleMetricsResult(
                                  // Alignment Summary Metrics
                                  totalReads: Long = 0,
                                  pfReads: Long = 0,
                                  pfReadsAligned: Long = 0,
                                  pctPfReadsAligned: Double = 0.0,
                                  pfHqAlignedReads: Long = 0,
                                  readsPaired: Long = 0,
                                  pctReadsPaired: Double = 0.0,
                                  pctProperPairs: Double = 0.0,
                                  meanReadLength: Double = 0.0,
                                  strandBalance: Double = 0.0,
                                  pctChimeras: Double = 0.0,
                                  pctAdapter: Double = 0.0,

                                  // Insert Size Metrics
                                  medianInsertSize: Double = 0.0,
                                  meanInsertSize: Double = 0.0,
                                  stdInsertSize: Double = 0.0,
                                  medianAbsoluteDeviation: Double = 0.0,
                                  minInsertSize: Int = 0,
                                  maxInsertSize: Int = 0,
                                  pairOrientation: String = "FR"
                                )

/**
 * Processor for running GATK CollectMultipleMetrics to gather alignment summary
 * and insert size metrics.
 *
 * CollectMultipleMetrics can run multiple collectors in a single pass through the BAM:
 * - CollectAlignmentSummaryMetrics: read counts, alignment rates, pair rates
 * - CollectInsertSizeMetrics: insert size distribution for paired-end libraries
 *
 * These complement CollectWgsMetrics (coverage) and CallableLoci (callable bases).
 */
class MultipleMetricsProcessor extends GatkToolProcessor[MultipleMetricsResult] {

  private val ARTIFACT_SUBDIR_NAME = "multiple_metrics"

  override protected def getToolName: String = "CollectMultipleMetrics"

  /**
   * Process a BAM/CRAM file to collect multiple metrics using GATK.
   *
   * @param bamPath         Path to the BAM/CRAM file
   * @param referencePath   Path to the reference genome
   * @param onProgress      Progress callback
   * @param artifactContext Optional context for organizing output artifacts by subject/run/alignment
   * @param totalReads      Optional total read count for progress estimation
   */
  def process(
               bamPath: String,
               referencePath: String,
               onProgress: (String, Double, Double) => Unit,
               artifactContext: Option[ArtifactContext] = None,
               totalReads: Option[Long] = None
             ): Either[Throwable, MultipleMetricsResult] = {

    executeGatkTool(
      bamPath,
      referencePath,
      onProgress,
      artifactContext,
      totalReads,
      buildArgs = (bam, ref, outPrefix) => Array(
        "CollectMultipleMetrics",
        "-I", bam,
        "-R", ref,
        "-O", outPrefix,
        "--PROGRAM", "CollectAlignmentSummaryMetrics",
        "--PROGRAM", "CollectInsertSizeMetrics",
        "--VALIDATION_STRINGENCY", "SILENT"
      ),
      parseOutput = parse,
      resolveOutputPath = resolveOutput
    )
  }

  private def resolveOutput(artifactContext: Option[ArtifactContext]): (Option[Path], String) = {
    artifactContext match {
      case Some(ctx) =>
        val dir = ctx.getSubdir(ARTIFACT_SUBDIR_NAME)
        (Some(dir), dir.resolve("metrics").toString)
      case None =>
        val tempDir = java.nio.file.Files.createTempDirectory("multiple_metrics")
        tempDir.toFile.deleteOnExit()
        (None, tempDir.resolve("metrics").toString)
    }
  }

  private def parse(outputPrefix: String): MultipleMetricsResult = {
    // Parse both output files
    val alignmentSummaryFile = new File(outputPrefix + ".alignment_summary_metrics")
    val insertSizeFile = new File(outputPrefix + ".insert_size_metrics")

    val alignmentMetrics = if (alignmentSummaryFile.exists()) {
      parseAlignmentSummaryMetrics(alignmentSummaryFile)
    } else {
      println(s"[MultipleMetricsProcessor] Warning: alignment_summary_metrics not found")
      Map.empty[String, String]
    }

    val insertSizeMetrics = if (insertSizeFile.exists()) {
      parseInsertSizeMetrics(insertSizeFile)
    } else {
      println(s"[MultipleMetricsProcessor] Warning: insert_size_metrics not found (may be single-end library)")
      Map.empty[String, String]
    }

    buildResult(alignmentMetrics, insertSizeMetrics)
  }

  /**
   * Parse CollectAlignmentSummaryMetrics output.
   * The file has multiple categories (FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR, UNPAIRED).
   * We want the PAIR row for paired-end libraries, or UNPAIRED for single-end.
   */
  private def parseAlignmentSummaryMetrics(file: File): Map[String, String] = {
    Using(Source.fromFile(file)) { source =>
      val lines = source.getLines().toList
      val headerIndex = lines.indexWhere(_.startsWith("CATEGORY"))

      if (headerIndex != -1) {
        val header = lines(headerIndex).split("\\t")
        val dataLines = lines.drop(headerIndex + 1).takeWhile(_.nonEmpty)

        // Find the PAIR row (preferred) or UNPAIRED row
        val categoryIndex = header.indexOf("CATEGORY")
        val pairRow = dataLines.find { line =>
          val values = line.split("\\t")
          categoryIndex >= 0 && values.length > categoryIndex && values(categoryIndex) == "PAIR"
        }
        val unpairedRow = dataLines.find { line =>
          val values = line.split("\\t")
          categoryIndex >= 0 && values.length > categoryIndex && values(categoryIndex) == "UNPAIRED"
        }

        (pairRow orElse unpairedRow) match {
          case Some(row) =>
            val values = row.split("\\t")
            header.zip(values).toMap
          case None =>
            // Fall back to first data row
            dataLines.headOption.map { row =>
              val values = row.split("\\t")
              header.zip(values).toMap
            }.getOrElse(Map.empty)
        }
      } else {
        Map.empty
      }
    }.getOrElse(Map.empty)
  }

  /**
   * Parse CollectInsertSizeMetrics output.
   * Takes the first data row (there's usually only one unless multiple libraries).
   */
  private def parseInsertSizeMetrics(file: File): Map[String, String] = {
    Using(Source.fromFile(file)) { source =>
      val lines = source.getLines().toList
      val headerIndex = lines.indexWhere(_.startsWith("MEDIAN_INSERT_SIZE"))

      if (headerIndex != -1 && lines.length > headerIndex + 1) {
        val header = lines(headerIndex).split("\\t")
        val values = lines(headerIndex + 1).split("\\t")
        header.zip(values).toMap
      } else {
        Map.empty
      }
    }.getOrElse(Map.empty)
  }

  /**
   * Build the result object from parsed metrics maps.
   */
  private def buildResult(
                           alignment: Map[String, String],
                           insertSize: Map[String, String]
                         ): MultipleMetricsResult = {

    def getLong(map: Map[String, String], key: String): Long =
      map.get(key).flatMap(s => scala.util.Try(s.toLong).toOption).getOrElse(0L)

    def getDouble(map: Map[String, String], key: String): Double =
      map.get(key).flatMap(s => scala.util.Try(s.toDouble).toOption).getOrElse(0.0)

    def getInt(map: Map[String, String], key: String): Int =
      map.get(key).flatMap(s => scala.util.Try(s.toInt).toOption).getOrElse(0)

    def getString(map: Map[String, String], key: String, default: String): String =
      map.getOrElse(key, default)

    // Calculate proper pair percentage from improper pair percentage
    val pctImproperPairs = getDouble(alignment, "PCT_PF_READS_IMPROPER_PAIRS")
    val pctProperPairs = if (pctImproperPairs > 0) 1.0 - pctImproperPairs else 0.0

    // Calculate paired read percentage
    val pfReadsAligned = getLong(alignment, "PF_READS_ALIGNED")
    val readsPaired = getLong(alignment, "READS_ALIGNED_IN_PAIRS")
    val pctReadsPaired = if (pfReadsAligned > 0) readsPaired.toDouble / pfReadsAligned else 0.0

    MultipleMetricsResult(
      // Alignment summary metrics
      totalReads = getLong(alignment, "TOTAL_READS"),
      pfReads = getLong(alignment, "PF_READS"),
      pfReadsAligned = pfReadsAligned,
      pctPfReadsAligned = getDouble(alignment, "PCT_PF_READS_ALIGNED"),
      pfHqAlignedReads = getLong(alignment, "PF_HQ_ALIGNED_READS"),
      readsPaired = readsPaired,
      pctReadsPaired = pctReadsPaired,
      pctProperPairs = pctProperPairs,
      meanReadLength = getDouble(alignment, "MEAN_READ_LENGTH"),
      strandBalance = getDouble(alignment, "STRAND_BALANCE"),
      pctChimeras = getDouble(alignment, "PCT_CHIMERAS"),
      pctAdapter = getDouble(alignment, "PCT_ADAPTER"),

      // Insert size metrics
      medianInsertSize = getDouble(insertSize, "MEDIAN_INSERT_SIZE"),
      meanInsertSize = getDouble(insertSize, "MEAN_INSERT_SIZE"),
      stdInsertSize = getDouble(insertSize, "STANDARD_DEVIATION"),
      medianAbsoluteDeviation = getDouble(insertSize, "MEDIAN_ABSOLUTE_DEVIATION"),
      minInsertSize = getInt(insertSize, "MIN_INSERT_SIZE"),
      maxInsertSize = getInt(insertSize, "MAX_INSERT_SIZE"),
      pairOrientation = getString(insertSize, "PAIR_ORIENTATION", "FR")
    )
  }
}
