package com.decodingus.analysis

import scala.io.Source
import scala.util.Using

/**
 * Parsed result from a samtools flagstat output file.
 * All counts are QC-passed reads only (left side of the + separator).
 */
case class FlagstatResult(
  totalReads: Long = 0,
  secondary: Long = 0,
  supplementary: Long = 0,
  duplicates: Long = 0,
  mapped: Long = 0,
  mappedPercent: Option[Double] = None,
  paired: Long = 0,
  read1: Long = 0,
  read2: Long = 0,
  properlyPaired: Long = 0,
  properlyPairedPercent: Option[Double] = None,
  withItselfAndMateMapped: Long = 0,
  singletons: Long = 0,
  singletonsPercent: Option[Double] = None,
  mateMappedToDiffChr: Long = 0,
  mateMappedToDiffChrMapQ5: Long = 0
) {
  /** Primary (non-secondary, non-supplementary) reads */
  def primaryReads: Long = totalReads - secondary - supplementary

  /** Duplication rate as a fraction */
  def duplicationRate: Option[Double] =
    if (totalReads > 0) Some(duplicates.toDouble / totalReads) else None

  /** Mapping rate as a fraction */
  def mappingRate: Option[Double] =
    if (totalReads > 0) Some(mapped.toDouble / totalReads) else None

  /** Whether data appears to be paired-end */
  def isPairedEnd: Boolean = paired > 0
}

/**
 * Parser for samtools flagstat output files.
 *
 * Expected format (samtools 1.x+):
 * {{{
 * 12345678 + 0 in total (QC-passed reads + QC-failed reads)
 * 0 + 0 secondary
 * 0 + 0 supplementary
 * 123456 + 0 duplicates
 * 12000000 + 0 mapped (97.20% : N/A)
 * 12345678 + 0 paired in sequencing
 * ...
 * }}}
 */
object FlagstatParser {

  // Pattern: "12345678 + 0 <category>" with optional "(xx.xx% : N/A)"
  private val LinePattern = """^(\d+)\s*\+\s*\d+\s+(.+)$""".r
  private val PercentPattern = """\((\d+\.?\d*)%""".r

  def parse(file: java.io.File): Either[String, FlagstatResult] = {
    Using(Source.fromFile(file)) { source =>
      val lines = source.getLines().toList
      if (lines.isEmpty) then FlagstatResult()
      else {

      var result = FlagstatResult()

      lines.foreach {
        case LinePattern(countStr, rest) =>
          val count = countStr.toLong
          val pct = PercentPattern.findFirstMatchIn(rest).map(_.group(1).toDouble)
          val category = rest.toLowerCase.trim

          if (category.startsWith("in total"))
            result = result.copy(totalReads = count)
          else if (category.startsWith("secondary"))
            result = result.copy(secondary = count)
          else if (category.startsWith("supplementary"))
            result = result.copy(supplementary = count)
          else if (category.startsWith("duplicates"))
            result = result.copy(duplicates = count)
          else if (category.startsWith("mapped"))
            result = result.copy(mapped = count, mappedPercent = pct)
          else if (category.startsWith("paired in sequencing"))
            result = result.copy(paired = count)
          else if (category.startsWith("read1"))
            result = result.copy(read1 = count)
          else if (category.startsWith("read2"))
            result = result.copy(read2 = count)
          else if (category.startsWith("properly paired"))
            result = result.copy(properlyPaired = count, properlyPairedPercent = pct)
          else if (category.startsWith("with itself and mate mapped"))
            result = result.copy(withItselfAndMateMapped = count)
          else if (category.startsWith("singletons"))
            result = result.copy(singletons = count, singletonsPercent = pct)
          else if (category.contains("mate mapped to a different chr") && category.contains("mapq"))
            result = result.copy(mateMappedToDiffChrMapQ5 = count)
          else if (category.contains("mate mapped to a different chr"))
            result = result.copy(mateMappedToDiffChr = count)

        case _ => // skip non-matching lines
      }

      result
      } // end else
    }.toEither.left.map(e => s"Failed to parse flagstat file ${file.getName}: ${e.getMessage}")
  }
}
