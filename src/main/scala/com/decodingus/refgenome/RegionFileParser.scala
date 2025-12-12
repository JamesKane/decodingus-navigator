package com.decodingus.refgenome

import java.nio.file.Path
import scala.io.Source
import scala.util.Using

/**
 * Record from a GFF3 file.
 *
 * GFF3 format has 9 tab-separated columns:
 * seqid, source, type, start (1-based), end (1-based, inclusive), score, strand, phase, attributes
 *
 * @param seqId Chromosome or contig (e.g., "chrY")
 * @param source Annotation source (e.g., "ybrowse")
 * @param featureType Feature type (e.g., "cytoband", "palindrome", "str")
 * @param start Start position (1-based, inclusive)
 * @param end End position (1-based, inclusive)
 * @param score Score field (often "." for not applicable)
 * @param strand Strand ("+" or "-" or ".")
 * @param phase Phase for CDS features (0, 1, 2, or ".")
 * @param attributes Key-value attribute pairs (e.g., Name=P8; Note=...)
 */
case class Gff3Record(
  seqId: String,
  source: String,
  featureType: String,
  start: Long,
  end: Long,
  score: Option[Double],
  strand: String,
  phase: Option[Int],
  attributes: Map[String, String]
) {
  def getAttribute(key: String): Option[String] = attributes.get(key)
  def name: Option[String] = getAttribute("Name").orElse(getAttribute("name"))
  def note: Option[String] = getAttribute("Note").orElse(getAttribute("note"))
}

/**
 * Record from a BED file.
 *
 * BED format uses 0-based, half-open coordinates [start, end).
 *
 * @param chrom Chromosome (e.g., "chrY")
 * @param start Start position (0-based, inclusive)
 * @param end End position (0-based, exclusive)
 * @param name Feature name (optional, column 4)
 * @param score Score (optional, column 5)
 * @param strand Strand (optional, column 6)
 */
case class BedRecord(
  chrom: String,
  start: Long,
  end: Long,
  name: Option[String] = None,
  score: Option[Double] = None,
  strand: Option[String] = None
)

/**
 * Parser for GFF3 and BED format genomic region files.
 */
object RegionFileParser {

  /**
   * Parse a GFF3 file.
   *
   * @param path Path to the GFF3 file
   * @return List of GFF3 records
   */
  def parseGff3(path: Path): Either[String, List[Gff3Record]] = {
    try {
      Using.resource(Source.fromFile(path.toFile)) { source =>
        val records = source.getLines()
          .filterNot(line => line.startsWith("#") || line.trim.isEmpty)
          .flatMap(parseGff3Line)
          .toList
        Right(records)
      }
    } catch {
      case e: Exception => Left(s"Failed to parse GFF3 file: ${e.getMessage}")
    }
  }

  /**
   * Parse a BED file.
   *
   * @param path Path to the BED file
   * @return List of BED records
   */
  def parseBed(path: Path): Either[String, List[BedRecord]] = {
    try {
      Using.resource(Source.fromFile(path.toFile)) { source =>
        val records = source.getLines()
          .filterNot(line => line.startsWith("#") || line.startsWith("track") || line.startsWith("browser") || line.trim.isEmpty)
          .flatMap(parseBedLine)
          .toList
        Right(records)
      }
    } catch {
      case e: Exception => Left(s"Failed to parse BED file: ${e.getMessage}")
    }
  }

  /**
   * Parse a single GFF3 line.
   */
  private def parseGff3Line(line: String): Option[Gff3Record] = {
    val fields = line.split("\t")
    if (fields.length < 9) {
      None
    } else {
      try {
        val seqId = fields(0)
        val source = fields(1)
        val featureType = fields(2)
        val start = fields(3).toLong
        val end = fields(4).toLong
        val score = if (fields(5) == ".") None else Some(fields(5).toDouble)
        val strand = fields(6)
        val phase = if (fields(7) == ".") None else Some(fields(7).toInt)
        val attributes = parseGff3Attributes(fields(8))

        Some(Gff3Record(seqId, source, featureType, start, end, score, strand, phase, attributes))
      } catch {
        case _: NumberFormatException => None
      }
    }
  }

  /**
   * Parse GFF3 attributes column.
   * Format: key1=value1;key2=value2;...
   */
  private def parseGff3Attributes(attrStr: String): Map[String, String] = {
    if (attrStr == "." || attrStr.trim.isEmpty) {
      Map.empty
    } else {
      attrStr.split(";")
        .flatMap { pair =>
          val eqIdx = pair.indexOf('=')
          if (eqIdx > 0) {
            val key = pair.substring(0, eqIdx).trim
            val value = java.net.URLDecoder.decode(pair.substring(eqIdx + 1).trim, "UTF-8")
            Some(key -> value)
          } else {
            None
          }
        }
        .toMap
    }
  }

  /**
   * Parse a single BED line.
   */
  private def parseBedLine(line: String): Option[BedRecord] = {
    val fields = line.split("\t")
    if (fields.length < 3) {
      None
    } else {
      try {
        val chrom = fields(0)
        val start = fields(1).toLong
        val end = fields(2).toLong
        val name = if (fields.length > 3 && fields(3).nonEmpty && fields(3) != ".") Some(fields(3)) else None
        val score = if (fields.length > 4 && fields(4).nonEmpty && fields(4) != ".") {
          try { Some(fields(4).toDouble) } catch { case _: NumberFormatException => None }
        } else None
        val strand = if (fields.length > 5 && fields(5).nonEmpty && fields(5) != ".") Some(fields(5)) else None

        Some(BedRecord(chrom, start, end, name, score, strand))
      } catch {
        case _: NumberFormatException => None
      }
    }
  }

  /**
   * Convert BED coordinates to 1-based inclusive coordinates (like GFF3/VCF).
   * BED uses 0-based, half-open [start, end).
   * GFF3/VCF uses 1-based, closed [start, end].
   *
   * @param bedStart 0-based start
   * @param bedEnd 0-based exclusive end
   * @return (1-based start, 1-based end)
   */
  def bedToOneBased(bedStart: Long, bedEnd: Long): (Long, Long) = {
    (bedStart + 1, bedEnd)  // BED end is already at the correct position when converting
  }

  /**
   * Convert 1-based inclusive coordinates to BED format.
   *
   * @param start 1-based start
   * @param end 1-based end
   * @return (0-based start, 0-based exclusive end)
   */
  def oneBasedToBed(start: Long, end: Long): (Long, Long) = {
    (start - 1, end)
  }

  /**
   * Filter records to only include Y chromosome entries.
   */
  def filterYChromosome[T](records: List[T], getChrom: T => String): List[T] = {
    val yChromNames = Set("chrY", "Y", "chrY_", "NC_000024")
    records.filter { record =>
      val chrom = getChrom(record)
      yChromNames.exists(y => chrom.startsWith(y))
    }
  }
}
