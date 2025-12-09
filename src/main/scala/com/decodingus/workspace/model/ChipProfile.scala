package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * Chip/Array genotype profile for a biosample.
 * Stores metadata and summary statistics from chip data imports (23andMe, AncestryDNA, etc.)
 *
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.genotype).
 *
 * Note: Raw genotype calls are NOT stored in the workspace - they are processed locally
 * and only summary statistics and derived data (haplogroups, ancestry) are persisted.
 *
 * @param atUri                 The AT URI of this chip profile record
 * @param meta                  Record metadata for tracking changes and sync
 * @param biosampleRef          AT URI of the parent biosample
 * @param vendor                Chip vendor (23andMe, AncestryDNA, FamilyTreeDNA, MyHeritage, LivingDNA)
 * @param testTypeCode          Test type code (e.g., ARRAY_23ANDME_V5)
 * @param chipVersion           Chip version if detected from file
 * @param totalMarkersCalled    Total number of markers with valid calls
 * @param totalMarkersPossible  Total markers in the chip file
 * @param noCallRate            Percentage of no-calls (0.0-1.0)
 * @param yMarkersCalled        Number of Y-DNA markers called (if applicable)
 * @param mtMarkersCalled       Number of mtDNA markers called (if applicable)
 * @param autosomalMarkersCalled Number of autosomal markers called
 * @param hetRate               Heterozygosity rate (autosomal only)
 * @param importDate            When this chip data was imported
 * @param sourceFileHash        SHA-256 hash of the source file for deduplication
 * @param sourceFileName        Original filename
 * @param files                 Source files info
 */
case class ChipProfile(
  atUri: Option[String],
  meta: RecordMeta,
  biosampleRef: String,
  vendor: String,
  testTypeCode: String,
  chipVersion: Option[String] = None,
  totalMarkersCalled: Int,
  totalMarkersPossible: Int,
  noCallRate: Double,
  yMarkersCalled: Option[Int] = None,
  mtMarkersCalled: Option[Int] = None,
  autosomalMarkersCalled: Int,
  hetRate: Option[Double] = None,
  importDate: LocalDateTime,
  sourceFileHash: Option[String] = None,
  sourceFileName: Option[String] = None,
  files: List[FileInfo] = List.empty
) {

  /**
   * Overall call rate.
   */
  def callRate: Double = if (totalMarkersPossible > 0) {
    totalMarkersCalled.toDouble / totalMarkersPossible
  } else 0.0

  /**
   * Check if quality is acceptable for ancestry analysis.
   */
  def isAcceptableForAncestry: Boolean =
    noCallRate < 0.05 && autosomalMarkersCalled >= 100000

  /**
   * Check if Y-DNA coverage is sufficient for haplogroup analysis.
   */
  def hasSufficientYCoverage: Boolean =
    yMarkersCalled.exists(_ >= 50)

  /**
   * Check if mtDNA coverage is sufficient for haplogroup analysis.
   */
  def hasSufficientMtCoverage: Boolean =
    mtMarkersCalled.exists(_ >= 20)

  /**
   * Get status display string.
   */
  def status: String = {
    if (noCallRate > 0.10) "Poor Quality"
    else if (noCallRate > 0.05) "Acceptable"
    else "Good"
  }
}

object ChipProfile {
  /** Known vendor values */
  val KnownVendors: Set[String] = Set(
    "23andMe", "AncestryDNA", "FamilyTreeDNA", "MyHeritage", "LivingDNA"
  )

  /** Known test type codes for chips */
  val KnownTestTypes: Set[String] = Set(
    "ARRAY_23ANDME_V5", "ARRAY_23ANDME_V4", "ARRAY_ANCESTRY_V2",
    "ARRAY_FTDNA_FF", "ARRAY_MYHERITAGE", "ARRAY_LIVINGDNA"
  )
}
