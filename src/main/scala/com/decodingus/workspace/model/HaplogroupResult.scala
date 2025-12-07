package com.decodingus.workspace.model

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
 */
case class HaplogroupResult(
  haplogroupName: String,
  score: Double,
  matchingSnps: Option[Int] = None,
  mismatchingSnps: Option[Int] = None,
  ancestralMatches: Option[Int] = None,
  treeDepth: Option[Int] = None,
  lineagePath: Option[List[String]] = None,
  privateVariants: Option[PrivateVariantData] = None
)
