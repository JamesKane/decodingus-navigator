package com.decodingus.analysis

import com.decodingus.model.WgsMetrics
import org.broadinstitute.hellbender.Main

import java.io.File
import scala.io.Source
import scala.util.{Either, Left, Right, Using, Try}

class WgsMetricsProcessor {

  def process(bamPath: String, referencePath: String, onProgress: (String, Double, Double) => Unit): Either[Throwable, WgsMetrics] = {
    onProgress("Running GATK CollectWgsMetrics...", 0.0, 1.0)

    val outputFile = File.createTempFile("wgs_metrics", ".txt")
    outputFile.deleteOnExit()

    val args = Array(
      "CollectWgsMetrics",
      "-I", bamPath,
      "-R", referencePath,
      "-O", outputFile.getAbsolutePath,
      "--USE_FAST_ALGORITHM", "true",
      "--READ_LENGTH", "4000000", // Support ultra-long reads up to 4Mb
      // Relax reference validation - allows GRCh38 with/without alts, etc.
      "--VALIDATION_STRINGENCY", "LENIENT",
      "--disable-sequence-dictionary-validation", "true"
    )

    // Execute GATK Main and capture any exceptions
    val gatkResult = Try {
      Main.main(args)
    }

    gatkResult match {
      case scala.util.Success(_) =>
        onProgress("Parsing GATK CollectWgsMetrics output...", 0.9, 1.0)
        val metrics = parse(outputFile.getAbsolutePath)
        onProgress("GATK CollectWgsMetrics complete.", 1.0, 1.0)
        Right(metrics)
      case scala.util.Failure(exception) =>
        Left(new RuntimeException(s"GATK CollectWgsMetrics failed: ${exception.getMessage}", exception))
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
