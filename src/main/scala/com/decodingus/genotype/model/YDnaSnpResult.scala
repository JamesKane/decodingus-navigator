package com.decodingus.genotype.model

import io.circe.Codec
import java.time.LocalDateTime

/**
 * Result of a Y-DNA SNP test (positive, negative, or no-call).
 */
enum YSnpResult derives Codec.AsObject:
  case Positive   // Derived/mutated state
  case Negative   // Ancestral state
  case NoCall     // Unable to determine

/**
 * A single Y-DNA SNP call from a panel or pack test.
 *
 * @param snpName       Primary SNP name (e.g., "M269", "U106", "Z381")
 * @param alternateNames Other names for this SNP (e.g., "rs9786184", "S21")
 * @param position      GRCh38 position if known
 * @param result        Test result
 * @param testedDate    When this SNP was tested (for incremental results)
 * @param source        Source of this result (e.g., "YSEQ Panel", "FTDNA SNP Pack R1b")
 */
case class YDnaSnpCall(
  snpName: String,
  alternateNames: List[String] = List.empty,
  position: Option[Int] = None,
  result: YSnpResult,
  testedDate: Option[LocalDateTime] = None,
  source: Option[String] = None
) derives Codec.AsObject {

  /**
   * Check if this is a positive/derived result.
   */
  def isPositive: Boolean = result == YSnpResult.Positive

  /**
   * Check if this is a negative/ancestral result.
   */
  def isNegative: Boolean = result == YSnpResult.Negative

  /**
   * Get all names for this SNP (primary + alternates).
   */
  def allNames: List[String] = snpName :: alternateNames
}

/**
 * Summary of Y-DNA SNP panel results for a biosample.
 *
 * This model supports:
 * - Fixed panels (FTDNA SNP Packs, BISDNA)
 * - Continuous delivery (YSEQ) with incremental updates
 * - Merging results from multiple sources
 *
 * @param atUri             AT URI of this panel result record
 * @param meta              Record metadata
 * @param biosampleRef      AT URI of the parent biosample
 * @param testTypeCode      Test type code (YDNA_SNP_PACK_FTDNA, YDNA_PANEL_YSEQ, etc.)
 * @param vendor            Panel vendor
 * @param panelName         Specific panel name if applicable (e.g., "R1b Pack", "Haplogroup J Pack")
 * @param snpCalls          Individual SNP results
 * @param totalSnpsTested   Total number of SNPs tested
 * @param positiveCount     Number of positive/derived results
 * @param negativeCount     Number of negative/ancestral results
 * @param noCallCount       Number of no-calls
 * @param terminalSnp       Most downstream positive SNP (terminal haplogroup marker)
 * @param inferredHaplogroup Haplogroup inferred from panel results
 * @param firstTestedDate   Date of first results
 * @param lastUpdatedDate   Date of most recent results (for incremental panels)
 * @param sourceFiles       Source files that contributed to these results
 * @param notes             Optional notes about the results
 */
case class YDnaSnpPanelResult(
  atUri: Option[String],
  meta: com.decodingus.workspace.model.RecordMeta,
  biosampleRef: String,
  testTypeCode: String,
  vendor: String,
  panelName: Option[String] = None,
  snpCalls: List[YDnaSnpCall] = List.empty,
  totalSnpsTested: Int,
  positiveCount: Int,
  negativeCount: Int,
  noCallCount: Int,
  terminalSnp: Option[String] = None,
  inferredHaplogroup: Option[String] = None,
  firstTestedDate: LocalDateTime,
  lastUpdatedDate: LocalDateTime,
  sourceFiles: List[com.decodingus.workspace.model.FileInfo] = List.empty,
  notes: Option[String] = None
) {

  /**
   * Call rate for this panel.
   */
  def callRate: Double =
    if (totalSnpsTested > 0) (positiveCount + negativeCount).toDouble / totalSnpsTested
    else 0.0

  /**
   * Check if this panel has been updated since initial testing.
   * Note: YSEQ and other vendors provide full dumps each time,
   * so "updates" means re-importing a newer complete file.
   */
  def hasBeenReimported: Boolean = lastUpdatedDate.isAfter(firstTestedDate)

  /**
   * Get the number of new SNPs compared to a previous version.
   * Useful for showing what changed after a re-import.
   */
  def newSnpsComparedTo(previous: YDnaSnpPanelResult): Int =
    totalSnpsTested - previous.totalSnpsTested

  /**
   * Get positive SNPs sorted by likely phylogenetic depth (most upstream first).
   * This is a heuristic based on SNP naming conventions.
   */
  def positiveSnpsByDepth: List[YDnaSnpCall] =
    snpCalls.filter(_.isPositive).sortBy { call =>
      // Rough heuristic: shorter names and M/P prefixes tend to be older markers
      val name = call.snpName.toUpperCase
      val prefixWeight = name.head match {
        case 'M' => 0  // M markers are often older (M269, M170, etc.)
        case 'P' => 1  // P markers are often intermediate
        case 'L' => 2
        case 'U' => 3
        case 'S' => 4
        case 'Z' => 5
        case 'Y' => 6  // Y markers are often more recent discoveries
        case 'F' => 7  // FGC markers
        case 'A' => 8  // A markers from 1000 Genomes
        case 'B' => 9  // BY markers from Big Y
        case _ => 10
      }
      (prefixWeight, name.length, name)
    }
}

object YDnaSnpPanelResult {
  /**
   * Create from a list of SNP calls.
   */
  def fromCalls(
    biosampleRef: String,
    testTypeCode: String,
    vendor: String,
    calls: List[YDnaSnpCall],
    panelName: Option[String] = None
  ): YDnaSnpPanelResult = {
    val now = LocalDateTime.now()
    val positive = calls.count(_.isPositive)
    val negative = calls.count(_.isNegative)
    val noCall = calls.count(_.result == YSnpResult.NoCall)

    // Find terminal SNP (most downstream positive)
    val terminal = calls.filter(_.isPositive).lastOption.map(_.snpName)

    YDnaSnpPanelResult(
      atUri = None,
      meta = com.decodingus.workspace.model.RecordMeta.initial,
      biosampleRef = biosampleRef,
      testTypeCode = testTypeCode,
      vendor = vendor,
      panelName = panelName,
      snpCalls = calls,
      totalSnpsTested = calls.size,
      positiveCount = positive,
      negativeCount = negative,
      noCallCount = noCall,
      terminalSnp = terminal,
      inferredHaplogroup = None, // Set by haplogroup analysis
      firstTestedDate = now,
      lastUpdatedDate = now
    )
  }
}
