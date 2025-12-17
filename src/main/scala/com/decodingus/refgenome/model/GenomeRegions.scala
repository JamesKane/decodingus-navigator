package com.decodingus.refgenome.model

import io.circe.*
import io.circe.generic.semiauto.*

import java.time.Instant

/**
 * Centralized genome region metadata from the Decoding-Us API.
 * Provides curated centromeres, telomeres, cytobands, and Y-specific regions
 * from a single authoritative source.
 *
 * @param build       Reference genome build (GRCh38, GRCh37, CHM13v2)
 * @param version     Semantic version for cache invalidation
 * @param generatedAt Timestamp when this data was generated
 * @param chromosomes Map of chromosome name to region data
 */
case class GenomeRegions(
                          build: String,
                          version: String,
                          generatedAt: Instant,
                          chromosomes: Map[String, ChromosomeRegions]
                        )

/**
 * Region data for a single chromosome.
 *
 * @param length     Total length of the chromosome in base pairs
 * @param centromere Centromere region coordinates
 * @param telomeres  P-arm and Q-arm telomere regions
 * @param cytobands  Cytoband annotations for ideogram display
 * @param regions    Y chromosome-specific regions (only for chrY)
 * @param strMarkers Named STR markers (only for chrY)
 */
case class ChromosomeRegions(
                              length: Long,
                              centromere: Option[Region],
                              telomeres: Option[Telomeres],
                              cytobands: List[Cytoband],
                              regions: Option[YChromosomeRegions],
                              strMarkers: Option[List[StrMarker]]
                            )

/**
 * A genomic region with optional type and quality modifier.
 *
 * @param start      Start position (1-based, inclusive)
 * @param end        End position (1-based, inclusive)
 * @param regionType Optional type classification
 * @param modifier   Quality modifier for concordance weighting (1.0 = reliable, <1.0 = reduced confidence)
 */
case class Region(
                   start: Long,
                   end: Long,
                   regionType: Option[String] = None,
                   modifier: Option[Double] = None
                 ) {
  def length: Long = end - start + 1

  def contains(position: Long): Boolean = position >= start && position <= end

  def overlaps(other: Region): Boolean = start <= other.end && end >= other.start
}

/**
 * Telomere regions for both chromosome arms.
 *
 * @param p P-arm (short arm) telomere
 * @param q Q-arm (long arm) telomere
 */
case class Telomeres(
                      p: Region,
                      q: Region
                    )

/**
 * Cytoband annotation for chromosome ideogram display.
 *
 * @param name  Band name (e.g., "p36.33", "q11.21")
 * @param start Start position
 * @param end   End position
 * @param stain Giemsa stain pattern (gneg, gpos25, gpos50, gpos75, gpos100, acen, gvar, stalk)
 */
case class Cytoband(
                     name: String,
                     start: Long,
                     end: Long,
                     stain: String
                   ) {
  /** Whether this is a centromeric band */
  def isCentromeric: Boolean = stain == "acen"

  /** Whether this is a positive (dark) Giemsa band */
  def isPositive: Boolean = stain.startsWith("gpos")
}

/**
 * Named STR (Short Tandem Repeat) marker.
 *
 * @param name     Marker name (e.g., "DYS389I", "DYS456")
 * @param start    Start position
 * @param end      End position
 * @param period   Repeat unit length in base pairs
 * @param verified Whether this position has been manually verified for this build
 * @param note     Optional annotation (e.g., "Position estimated via liftover from GRCh38")
 */
case class StrMarker(
                      name: String,
                      start: Long,
                      end: Long,
                      period: Int,
                      verified: Boolean,
                      note: Option[String] = None
                    ) {
  def length: Long = end - start + 1
}

/**
 * Y chromosome-specific region annotations.
 * Includes PAR, XTR, ampliconic regions, palindromes, heterochromatin, and X-degenerate regions.
 *
 * @param par1            Pseudoautosomal region 1 (Yp)
 * @param par2            Pseudoautosomal region 2 (Yq)
 * @param xtr             X-transposed region
 * @param ampliconic      Ampliconic (high-copy) regions
 * @param palindromes     Palindromic regions (P1-P8)
 * @param heterochromatin Yq12 heterochromatin region
 * @param xDegenerate     X-degenerate (stable single-copy) regions
 */
case class YChromosomeRegions(
                               par1: Region,
                               par2: Region,
                               xtr: List[Region],
                               ampliconic: List[Region],
                               palindromes: List[NamedRegion],
                               heterochromatin: Region,
                               xDegenerate: List[Region]
                             )

/**
 * A named region with type and quality modifier.
 *
 * @param name       Region name (e.g., "P1", "P2" for palindromes)
 * @param start      Start position
 * @param end        End position
 * @param regionType Type classification (e.g., "Palindrome", "PAR")
 * @param modifier   Quality modifier for concordance weighting
 */
case class NamedRegion(
                        name: String,
                        start: Long,
                        end: Long,
                        regionType: String,
                        modifier: Double
                      ) {
  def length: Long = end - start + 1

  def toRegion: Region = Region(start, end, Some(regionType), Some(modifier))
}

/**
 * Circe JSON codecs for genome region types.
 */
object GenomeRegionsCodecs:

  // Instant codec (ISO-8601 format)
  given Encoder[Instant] = Encoder.encodeString.contramap(_.toString)

  given Decoder[Instant] = Decoder.decodeString.emap { s =>
    try Right(Instant.parse(s))
    catch case e: Exception => Left(s"Invalid timestamp: $s")
  }

  // Region codec
  given Encoder[Region] = Encoder.instance { r =>
    Json.obj(
      "start" -> Json.fromLong(r.start),
      "end" -> Json.fromLong(r.end),
      "type" -> r.regionType.fold(Json.Null)(Json.fromString),
      "modifier" -> r.modifier.fold(Json.Null)(Json.fromDoubleOrNull)
    ).dropNullValues
  }

  given Decoder[Region] = Decoder.instance { c =>
    for
      start <- c.get[Long]("start")
      end <- c.get[Long]("end")
      regionType <- c.get[Option[String]]("type")
      modifier <- c.get[Option[Double]]("modifier")
    yield Region(start, end, regionType, modifier)
  }

  // Telomeres codec
  given Encoder[Telomeres] = deriveEncoder[Telomeres]

  given Decoder[Telomeres] = deriveDecoder[Telomeres]

  // Cytoband codec
  given Encoder[Cytoband] = deriveEncoder[Cytoband]

  given Decoder[Cytoband] = deriveDecoder[Cytoband]

  // StrMarker codec
  given Encoder[StrMarker] = deriveEncoder[StrMarker]

  given Decoder[StrMarker] = deriveDecoder[StrMarker]

  // NamedRegion codec
  given Encoder[NamedRegion] = Encoder.instance { nr =>
    Json.obj(
      "name" -> Json.fromString(nr.name),
      "start" -> Json.fromLong(nr.start),
      "end" -> Json.fromLong(nr.end),
      "type" -> Json.fromString(nr.regionType),
      "modifier" -> Json.fromDoubleOrNull(nr.modifier)
    )
  }

  given Decoder[NamedRegion] = Decoder.instance { c =>
    for
      name <- c.get[String]("name")
      start <- c.get[Long]("start")
      end <- c.get[Long]("end")
      regionType <- c.get[String]("type")
      modifier <- c.get[Double]("modifier")
    yield NamedRegion(name, start, end, regionType, modifier)
  }

  // YChromosomeRegions codec
  given Encoder[YChromosomeRegions] = deriveEncoder[YChromosomeRegions]

  given Decoder[YChromosomeRegions] = deriveDecoder[YChromosomeRegions]

  // ChromosomeRegions codec
  given Encoder[ChromosomeRegions] = deriveEncoder[ChromosomeRegions]

  given Decoder[ChromosomeRegions] = deriveDecoder[ChromosomeRegions]

  // GenomeRegions codec
  given Encoder[GenomeRegions] = deriveEncoder[GenomeRegions]

  given Decoder[GenomeRegions] = deriveDecoder[GenomeRegions]
