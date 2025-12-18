package com.decodingus.workspace.model

import com.decodingus.str.StrPanelService
import java.time.LocalDateTime

/**
 * STR value types - either a simple repeat count or complex multi-allelic structure.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#strValue).
 */
sealed trait StrValue

/**
 * Simple single-value STR (e.g., DYS393 = 13).
 */
case class SimpleStrValue(repeats: Int) extends StrValue

/**
 * Multi-copy STR with ordered values (e.g., DYS385a/b = 11-14, DYS459a/b = 9-10).
 * Convention: copies are in ascending order.
 */
case class MultiCopyStrValue(copies: List[Int]) extends StrValue

/**
 * A single allele in a complex STR marker.
 *
 * @param repeats     Repeat count (can be fractional like 26.1 for partial repeats)
 * @param count       Number of copies of this allele (e.g., 2 for 'c' = cis/both copies)
 * @param designation Allele designation letter if applicable (t, c, q)
 */
case class StrAllele(
                      repeats: Double,
                      count: Int,
                      designation: Option[String] = None
                    )

/**
 * Complex multi-allelic STR with allele counts (e.g., DYF399X = 22t-25c-26.1t).
 * Used for palindromic markers.
 *
 * @param alleles     List of alleles with their repeat values and counts
 * @param rawNotation Original notation string for reference (e.g., '22t-25c-26.1t')
 */
case class ComplexStrValue(
                            alleles: List[StrAllele],
                            rawNotation: Option[String] = None
                          ) extends StrValue

/**
 * A single STR marker value. Handles simple, multi-copy, and complex multi-allelic markers.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#strMarkerValue).
 *
 * @param marker        Standard marker name (e.g., 'DYS393', 'DYS385a', 'DYF399X')
 * @param value         The marker value - simple integer or complex allele structure
 * @param startPosition GRCh38 start position on chrY (for WGS-derived or positioned markers)
 * @param endPosition   GRCh38 end position on chrY (repeat region boundary)
 * @param orderedDate   Date this specific marker was ordered/tested (for incremental panels)
 * @param panel         Which panel this marker belongs to (Y12, Y25, Y37, Y67, Y111, Y500, Y700, etc.)
 * @param quality       Call quality if available (HIGH, MEDIUM, LOW, UNCERTAIN)
 * @param readDepth     Read depth for WGS-derived STR calls
 */
case class StrMarkerValue(
                           marker: String,
                           value: StrValue,
                           startPosition: Option[Long] = None,
                           endPosition: Option[Long] = None,
                           orderedDate: Option[LocalDateTime] = None,
                           panel: Option[String] = None,
                           quality: Option[String] = None,
                           readDepth: Option[Int] = None
                         ) {
  /** Get the repeat region span in base pairs, if positions are known */
  def regionSpan: Option[Long] =
    for
      start <- startPosition
      end <- endPosition
    yield end - start + 1
}

/**
 * Metadata about an STR panel/test.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#strPanel).
 *
 * @param panelName   Panel name (e.g., 'Y-37', 'Big Y-700 STRs')
 * @param markerCount Number of markers in this panel
 * @param provider    Testing company or source (FTDNA, YSEQ, NEBULA, DANTE, WGS_DERIVED, OTHER)
 * @param testDate    When the test was performed
 */
case class StrPanel(
                     panelName: String,
                     markerCount: Int,
                     provider: Option[String] = None,
                     testDate: Option[LocalDateTime] = None
                   )

/**
 * Y-STR profile for a biosample. Can contain multiple panels from different sources.
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.strProfile).
 *
 * @param atUri            The AT URI of this STR profile record
 * @param meta             Record metadata for tracking changes and sync
 * @param biosampleRef     AT URI of the parent biosample
 * @param sequenceRunRef   AT URI of the sequence run if STRs were derived from WGS (optional)
 * @param panels           Panels/tests that contributed to this profile
 * @param markers          The STR marker values
 * @param totalMarkers     Total number of markers in this profile
 * @param source           How these STRs were obtained (DIRECT_TEST, WGS_DERIVED, BIG_Y_DERIVED, IMPORTED, MANUAL_ENTRY)
 * @param importedFrom     If imported, the original source (e.g., 'FTDNA', 'YSEQ', 'YFULL')
 * @param derivationMethod For WGS-derived STRs, the tool/method used (HIPSTR, GANGSTR, EXPANSION_HUNTER, LOBSTR, CUSTOM)
 * @param files            Source CSV/TSV files if available
 */
case class StrProfile(
                       atUri: Option[String],
                       meta: RecordMeta,
                       biosampleRef: String,
                       sequenceRunRef: Option[String] = None,
                       panels: List[StrPanel] = List.empty,
                       markers: List[StrMarkerValue] = List.empty,
                       totalMarkers: Option[Int] = None,
                       source: Option[String] = None,
                       importedFrom: Option[String] = None,
                       derivationMethod: Option[String] = None,
                       files: List[FileInfo] = List.empty
                     )

object StrProfile {
  /** Known panel values - loaded from config with fallback */
  lazy val KnownPanels: Set[String] = {
    val configPanels = StrPanelService.getKnownPanelNames
    if (configPanels.nonEmpty) configPanels + "CUSTOM"
    else Set("Y-12", "Y-25", "Y-37", "Y-67", "Y-111", "Y-500", "Y-700", "YSEQ_ALPHA", "CUSTOM")
  }

  /** Known provider values */
  val KnownProviders: Set[String] = Set(
    "FTDNA", "YSEQ", "NEBULA", "DANTE", "WGS_DERIVED", "OTHER"
  )

  /** Known source values */
  val KnownSources: Set[String] = Set(
    "DIRECT_TEST", "WGS_DERIVED", "BIG_Y_DERIVED", "IMPORTED", "MANUAL_ENTRY"
  )

  /** Known derivation method values */
  val KnownDerivationMethods: Set[String] = Set(
    "HIPSTR", "GANGSTR", "EXPANSION_HUNTER", "LOBSTR", "CUSTOM"
  )

  /** Known quality values */
  val KnownQualityValues: Set[String] = Set("HIGH", "MEDIUM", "LOW", "UNCERTAIN")
}
