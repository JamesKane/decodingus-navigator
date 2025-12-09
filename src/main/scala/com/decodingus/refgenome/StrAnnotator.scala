package com.decodingus.refgenome

import java.nio.file.Path
import scala.io.Source
import scala.util.Using
import scala.collection.mutable

/**
 * STR (Short Tandem Repeat) region from HipSTR reference.
 *
 * @param chrom Chromosome
 * @param start Start position (0-based)
 * @param end End position (0-based, exclusive)
 * @param period Repeat unit length (e.g., 2 for dinucleotide, 3 for trinucleotide)
 * @param numRepeats Number of repeat copies in reference
 * @param name STR locus name (if available)
 */
case class StrRegion(
  chrom: String,
  start: Long,
  end: Long,
  period: Int,
  numRepeats: Double,
  name: Option[String]
)

/**
 * Annotates genomic positions with STR (Short Tandem Repeat) information.
 * Uses HipSTR reference BED files to identify known STR regions.
 *
 * BED format from HipSTR:
 * chrom  start  end  period  num_repeats  [name]
 */
class StrAnnotator(bedPath: Path) {
  // Index STR regions by chromosome for fast lookup
  private val regionsByChrom: Map[String, IndexedSeq[StrRegion]] = loadRegions()

  private def loadRegions(): Map[String, IndexedSeq[StrRegion]] = {
    val byChrom = mutable.Map[String, mutable.ArrayBuffer[StrRegion]]()

    Using.resource(Source.fromFile(bedPath.toFile)) { source =>
      for (line <- source.getLines() if !line.startsWith("#") && line.nonEmpty) {
        val fields = line.split("\t")
        if (fields.length >= 5) {
          val chrom = fields(0)
          val start = fields(1).toLong
          val end = fields(2).toLong
          val period = fields(3).toInt
          val numRepeats = fields(4).toDouble
          val name = if (fields.length > 5) Some(fields(5)) else None

          val region = StrRegion(chrom, start, end, period, numRepeats, name)
          byChrom.getOrElseUpdate(chrom, mutable.ArrayBuffer()) += region
        }
      }
    }

    // Sort each chromosome's regions by start position for binary search
    byChrom.view.mapValues(_.sortBy(_.start).toIndexedSeq).toMap
  }

  /**
   * Find STR region overlapping a given position.
   *
   * @param chrom Chromosome (e.g., "chrY", "Y")
   * @param position 1-based genomic position
   * @return Some(StrRegion) if position falls within a known STR, None otherwise
   */
  def findOverlapping(chrom: String, position: Long): Option[StrRegion] = {
    // Try both with and without "chr" prefix
    val chromVariants = List(chrom, s"chr$chrom", chrom.stripPrefix("chr"))

    chromVariants.flatMap { c =>
      regionsByChrom.get(c).flatMap { regions =>
        // Binary search for potential overlapping region
        // Position is 1-based, BED is 0-based
        val pos0 = position - 1
        binarySearchOverlap(regions, pos0)
      }
    }.headOption
  }

  private def binarySearchOverlap(regions: IndexedSeq[StrRegion], pos: Long): Option[StrRegion] = {
    // Find the rightmost region with start <= pos
    var lo = 0
    var hi = regions.length - 1
    var result: Option[StrRegion] = None

    while (lo <= hi) {
      val mid = lo + (hi - lo) / 2
      val region = regions(mid)

      if (region.start <= pos) {
        // Check if this region contains the position
        if (pos < region.end) {
          result = Some(region)
        }
        lo = mid + 1
      } else {
        hi = mid - 1
      }
    }

    // Also check the region just before our result in case of overlapping regions
    if (result.isEmpty && hi >= 0) {
      val region = regions(hi)
      if (region.start <= pos && pos < region.end) {
        result = Some(region)
      }
    }

    result
  }

  /**
   * Get repeat type description based on period.
   */
  def repeatTypeDescription(period: Int): String = period match {
    case 1 => "homopolymer"
    case 2 => "dinucleotide"
    case 3 => "trinucleotide"
    case 4 => "tetranucleotide"
    case 5 => "pentanucleotide"
    case 6 => "hexanucleotide"
    case n => s"${n}mer"
  }

  /**
   * Format STR annotation for display.
   */
  def formatAnnotation(region: StrRegion): String = {
    val repeatType = repeatTypeDescription(region.period)
    val nameStr = region.name.map(n => s" ($n)").getOrElse("")
    f"$repeatType, ${region.numRepeats}%.1f copies$nameStr"
  }

  /** Total number of STR regions loaded */
  def regionCount: Int = regionsByChrom.values.map(_.size).sum

  /** Number of STR regions for a specific chromosome */
  def regionCount(chrom: String): Int = {
    val chromVariants = List(chrom, s"chr$chrom", chrom.stripPrefix("chr"))
    chromVariants.flatMap(c => regionsByChrom.get(c).map(_.size)).headOption.getOrElse(0)
  }
}

object StrAnnotator {
  /**
   * Create an STR annotator for the given reference build.
   * Downloads HipSTR reference if not cached.
   */
  def forBuild(referenceBuild: String): Either[String, StrAnnotator] = {
    val gateway = new StrReferenceGateway()
    gateway.resolve(referenceBuild).map(path => new StrAnnotator(path))
  }
}
