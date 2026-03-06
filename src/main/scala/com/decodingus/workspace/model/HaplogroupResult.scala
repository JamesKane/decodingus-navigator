package com.decodingus.workspace.model

import java.time.Instant

/**
 * Detailed scoring and classification result for a haplogroup.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#haplogroupResult).
 *
 * @param haplogroupName   The assigned haplogroup nomenclature (e.g., R-M269, H1a)
 * @param score            Confidence score of the assignment
 * @param matchingSnps     Count of SNPs matching the defining mutations for this haplogroup
 * @param mismatchingSnps  Count of SNPs that contradict the assignment (potential private variants)
 * @param ancestralMatches Count of ancestral state matches
 * @param treeDepth        The depth of the assigned node in the phylogenetic tree
 * @param lineagePath      The path from root to the assigned haplogroup (e.g., A -> ... -> R -> ... -> R-M269)
 * @param privateVariants  Detailed private variant calls for haplogroup discovery (optional)
 * @param source           Data source type: "wgs", "bigy", "chip"
 * @param sourceRef        AT URI of the source record (SequenceRun or ChipProfile)
 * @param treeProvider     Tree provider used: "ftdna", "decodingus"
 * @param treeVersion      Version of the tree used (e.g., "2024-12-01")
 * @param analyzedAt       Timestamp when this analysis was performed
 */
case class HaplogroupResult(
                             haplogroupName: String,
                             score: Double,
                             matchingSnps: Option[Int] = None,
                             mismatchingSnps: Option[Int] = None,
                             ancestralMatches: Option[Int] = None,
                             treeDepth: Option[Int] = None,
                             lineagePath: Option[List[String]] = None,
                             privateVariants: Option[PrivateVariantData] = None,
                             source: Option[String] = None,
                             sourceRef: Option[String] = None,
                             treeProvider: Option[String] = None,
                             treeVersion: Option[String] = None,
                             analyzedAt: Option[Instant] = None
                           ) {

  /**
   * Quality tier for reconciliation purposes.
   * Higher tier = more trusted/detailed result.
   */
  def qualityTier: Int = source match {
    case Some("wgs") => 3 // WGS - highest quality
    case Some("bigy") => 2 // Targeted Y-DNA (Big Y, Y Elite)
    case Some("chip") => 1 // SNP array/chip
    case _ => 0 // Unknown
  }

  /**
   * Check if this result is from the same analysis run as another.
   * Same source ref + same tree provider = same analysis.
   */
  def isSameAnalysis(other: HaplogroupResult): Boolean =
    sourceRef == other.sourceRef && treeProvider == other.treeProvider
}
