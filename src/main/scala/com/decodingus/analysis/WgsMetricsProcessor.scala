package com.decodingus.analysis

import com.decodingus.model.WgsMetrics

import java.io.File
import java.nio.file.Path
import scala.io.Source
import scala.util.{Either, Left, Right, Using}

class WgsMetricsProcessor extends GatkToolProcessor[WgsMetrics] {

  override protected def getToolName: String = "CollectWgsMetrics"

  /**
   * Process a BAM/CRAM file to collect WGS metrics using GATK.
   *
   * @param bamPath         Path to the BAM/CRAM file
   * @param referencePath   Path to the reference genome
   * @param onProgress      Progress callback
   * @param readLength      Optional read length - if > 150bp, passed to GATK to avoid crashes with long reads (e.g., PacBio HiFi)
   * @param artifactContext Optional context for organizing output artifacts by subject/run/alignment
   * @param totalReads      Optional total read count for progress estimation
   * @param countUnpaired   If true, count unpaired reads (needed for single-end long-read data like PacBio HiFi)
   */
  def process(
               bamPath: String,
               referencePath: String,
               onProgress: (String, Double, Double) => Unit,
               readLength: Option[Int] = None,
               artifactContext: Option[ArtifactContext] = None,
               totalReads: Option[Long] = None,
               countUnpaired: Boolean = false
             ): Either[Throwable, WgsMetrics] = {

    executeGatkTool(
      bamPath,
      referencePath,
      onProgress,
      artifactContext,
      totalReads,
      buildArgs = (bam, ref, out) => {
        val baseArgs = Array(
          "CollectWgsMetrics",
          "-I", bam,
          "-R", ref,
          "-O", out,
          "--USE_FAST_ALGORITHM", "true",
          "--VALIDATION_STRINGENCY", "SILENT"
        )

        val argsWithReadLength = readLength match {
          case Some(len) if len > 150 => baseArgs ++ Array("--READ_LENGTH", len.toString)
          case _ => baseArgs
        }

        if (countUnpaired) {
          argsWithReadLength ++ Array("--COUNT_UNPAIRED", "true")
        } else {
          argsWithReadLength
        }
      },
      parseOutput = parse,
      resolveOutputPath = resolveOutput
    )
  }

  private def resolveOutput(artifactContext: Option[ArtifactContext]): (Option[Path], String) = {
    artifactContext match {
      case Some(ctx) =>
        val dir = ctx.getArtifactDir
        (Some(dir), dir.resolve("wgs_metrics.txt").toFile.getAbsolutePath)
      case None =>
        val tempFile = File.createTempFile("wgs_metrics", ".txt")
        tempFile.deleteOnExit()
        (None, tempFile.getAbsolutePath)
    }
  }

  private def parse(filePath: String): WgsMetrics = {
    Using(Source.fromFile(filePath)) {
      source =>
        val lines = source.getLines().toList
        val headerIndex = lines.indexWhere(_.startsWith("GENOME_TERRITORY"))
        if (headerIndex != -1 && lines.length > headerIndex + 1) {
          val header = lines(headerIndex).split("\\t")
          val values = lines(headerIndex + 1).split("\\t")
          val metricsMap = header.zip(values).toMap

          WgsMetrics(
            genomeTerritory = metricsMap.getOrElse("GENOME_TERRITORY", "0").toLong,
            meanCoverage = metricsMap.getOrElse("MEAN_COVERAGE", "0.0").toDouble,
            sdCoverage = metricsMap.getOrElse("SD_COVERAGE", "0.0").toDouble,
            medianCoverage = metricsMap.getOrElse("MEDIAN_COVERAGE", "0.0").toDouble,
            madCoverage = metricsMap.getOrElse("MAD_COVERAGE", "0.0").toDouble,
            pctExcMapq = metricsMap.getOrElse("PCT_EXC_MAPQ", "0.0").toDouble,
            pctExcDupe = metricsMap.getOrElse("PCT_EXC_DUPE", "0.0").toDouble,
            pctExcUnpaired = metricsMap.getOrElse("PCT_EXC_UNPAIRED", "0.0").toDouble,
            pctExcBaseq = metricsMap.getOrElse("PCT_EXC_BASEQ", "0.0").toDouble,
            pctExcOverlap = metricsMap.getOrElse("PCT_EXC_OVERLAP", "0.0").toDouble,
            pctExcCapped = metricsMap.getOrElse("PCT_EXC_CAPPED", "0.0").toDouble,
            pctExcTotal = metricsMap.getOrElse("PCT_EXC_TOTAL", "0.0").toDouble,
            pct1x = metricsMap.getOrElse("PCT_1X", "0.0").toDouble,
            pct5x = metricsMap.getOrElse("PCT_5X", "0.0").toDouble,
            pct10x = metricsMap.getOrElse("PCT_10X", "0.0").toDouble,
            pct15x = metricsMap.getOrElse("PCT_15X", "0.0").toDouble,
            pct20x = metricsMap.getOrElse("PCT_20X", "0.0").toDouble,
            pct25x = metricsMap.getOrElse("PCT_25X", "0.0").toDouble,
            pct30x = metricsMap.getOrElse("PCT_30X", "0.0").toDouble,
            pct40x = metricsMap.getOrElse("PCT_40X", "0.0").toDouble,
            pct50x = metricsMap.getOrElse("PCT_50X", "0.0").toDouble,
            pct60x = metricsMap.getOrElse("PCT_60X", "0.0").toDouble,
            pct70x = metricsMap.getOrElse("PCT_70X", "0.0").toDouble,
            pct80x = metricsMap.getOrElse("PCT_80X", "0.0").toDouble,
            pct90x = metricsMap.getOrElse("PCT_90X", "0.0").toDouble,
            pct100x = metricsMap.getOrElse("PCT_100X", "0.0").toDouble,
            hetSnpSensitivity = metricsMap.getOrElse("HET_SNP_SENSITIVITY", "0.0").toDouble,
            hetSnpQ = metricsMap.getOrElse("HET_SNP_Q", "0.0").toDouble
          )
        } else {
          WgsMetrics() // Return default if parsing fails
        }
    }.getOrElse(WgsMetrics())
  }
}
