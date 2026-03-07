package com.decodingus.analysis

import com.decodingus.model.WgsMetrics
import com.decodingus.workspace.model.AlignmentMetrics

import java.io.File
import scala.io.Source
import scala.util.Using

/**
 * Loads pre-computed GATK/Picard/samtools metrics from files in a sample directory.
 * Used during project import to populate AlignmentMetrics without re-running analysis.
 */
object MetricsFileLoader {

  /**
   * Load all available pre-computed metrics from a discovered sample's metrics files
   * and merge them into an AlignmentMetrics.
   */
  def loadMetrics(sample: DiscoveredSample): Option[AlignmentMetrics] = {
    val flagstat = sample.flagstatFiles.headOption.flatMap(f => FlagstatParser.parse(f).toOption)
    val wgsMetrics = sample.wgsMetricsFiles.headOption.flatMap(f => parseWgsMetricsFile(f).toOption)

    if (flagstat.isEmpty && wgsMetrics.isEmpty) return None

    var metrics = AlignmentMetrics()

    // Apply WGS metrics
    wgsMetrics.foreach { wgs =>
      metrics = metrics.copy(
        genomeTerritory = Some(wgs.genomeTerritory),
        meanCoverage = Some(wgs.meanCoverage),
        medianCoverage = Some(wgs.medianCoverage),
        sdCoverage = Some(wgs.sdCoverage),
        pctExcDupe = Some(wgs.pctExcDupe),
        pctExcMapq = Some(wgs.pctExcMapq),
        pct10x = Some(wgs.pct10x),
        pct20x = Some(wgs.pct20x),
        pct30x = Some(wgs.pct30x),
        hetSnpSensitivity = Some(wgs.hetSnpSensitivity)
      )
    }

    Some(metrics)
  }

  /**
   * Extract SequenceRun-level stats from flagstat results.
   * Returns (totalReads, pfReadsAligned, pctPfReadsAligned, properlyPairedPct, isPairedEnd).
   */
  def extractSequenceRunStats(sample: DiscoveredSample): Option[FlagstatResult] =
    sample.flagstatFiles.headOption.flatMap(f => FlagstatParser.parse(f).toOption)

  /**
   * Parse a GATK/Picard WGS metrics file.
   * Reuses the same TSV format as WgsMetricsProcessor output.
   */
  def parseWgsMetricsFile(file: File): Either[String, WgsMetrics] = {
    Using(Source.fromFile(file)) { source =>
      val lines = source.getLines().toList
      val headerIndex = lines.indexWhere(_.startsWith("GENOME_TERRITORY"))

      if (headerIndex != -1 && lines.length > headerIndex + 1) {
        val header = lines(headerIndex).split("\\t")
        val values = lines(headerIndex + 1).split("\\t")
        val m = header.zip(values).toMap

        WgsMetrics(
          genomeTerritory = safeLong(m, "GENOME_TERRITORY"),
          meanCoverage = safeDouble(m, "MEAN_COVERAGE"),
          sdCoverage = safeDouble(m, "SD_COVERAGE"),
          medianCoverage = safeDouble(m, "MEDIAN_COVERAGE"),
          madCoverage = safeDouble(m, "MAD_COVERAGE"),
          pctExcMapq = safeDouble(m, "PCT_EXC_MAPQ"),
          pctExcDupe = safeDouble(m, "PCT_EXC_DUPE"),
          pctExcUnpaired = safeDouble(m, "PCT_EXC_UNPAIRED"),
          pctExcBaseq = safeDouble(m, "PCT_EXC_BASEQ"),
          pctExcOverlap = safeDouble(m, "PCT_EXC_OVERLAP"),
          pctExcCapped = safeDouble(m, "PCT_EXC_CAPPED"),
          pctExcTotal = safeDouble(m, "PCT_EXC_TOTAL"),
          pct1x = safeDouble(m, "PCT_1X"),
          pct5x = safeDouble(m, "PCT_5X"),
          pct10x = safeDouble(m, "PCT_10X"),
          pct15x = safeDouble(m, "PCT_15X"),
          pct20x = safeDouble(m, "PCT_20X"),
          pct25x = safeDouble(m, "PCT_25X"),
          pct30x = safeDouble(m, "PCT_30X"),
          pct40x = safeDouble(m, "PCT_40X"),
          pct50x = safeDouble(m, "PCT_50X"),
          pct60x = safeDouble(m, "PCT_60X"),
          pct70x = safeDouble(m, "PCT_70X"),
          pct80x = safeDouble(m, "PCT_80X"),
          pct90x = safeDouble(m, "PCT_90X"),
          pct100x = safeDouble(m, "PCT_100X"),
          hetSnpSensitivity = safeDouble(m, "HET_SNP_SENSITIVITY"),
          hetSnpQ = safeDouble(m, "HET_SNP_Q")
        )
      } else {
        WgsMetrics()
      }
    }.toEither.left.map(e => s"Failed to parse WGS metrics file ${file.getName}: ${e.getMessage}")
  }

  private def safeDouble(m: Map[String, String], key: String): Double =
    m.get(key).flatMap(v => scala.util.Try(v.toDouble).toOption).getOrElse(0.0)

  private def safeLong(m: Map[String, String], key: String): Long =
    m.get(key).flatMap(v => scala.util.Try(v.toLong).toOption).getOrElse(0L)
}
