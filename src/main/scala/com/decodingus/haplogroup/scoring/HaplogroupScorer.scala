package com.decodingus.haplogroup.scoring

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, HaplogroupScore, Locus}

import scala.collection.mutable

/**
 * Scores haplogroups by comparing sample SNP calls against the haplogroup tree.
 *
 * Scoring rules:
 * - Each branch is scored based on the proportion of derived SNPs among callable positions
 * - Branches with many phylo-equivalent SNPs (starbursts from population bottlenecks)
 *   contribute heavily when most/all SNPs are derived
 * - Branch score = (derived / callable) * callable, where callable = derived + ancestral
 *   This rewards both high match rate AND high SNP density
 * - No calls are neutral (don't affect the ratio)
 * - Stop descent if two consecutive branches lack positive evidence (all ancestral OR all no-calls)
 *   This prevents traversing deep into unconfirmed paths when branch-defining SNPs have no coverage
 */
class HaplogroupScorer {

  def score(tree: List[Haplogroup], snpCalls: Map[Long, String]): List[HaplogroupResult] = {
    val scores = mutable.ListBuffer[HaplogroupResult]()
    tree.foreach(rootNode => calculateHaplogroupScore(
      rootNode,
      snpCalls,
      scores,
      cumulativeScore = 0.0,
      cumulativeDerived = 0,
      cumulativeAncestral = 0,
      cumulativeNoCalls = 0,
      cumulativeSnps = 0,
      depth = 0,
      consecutiveAncestralBranches = 0
    ))
    // Sort by descending score, then ASCENDING depth (prefer parent when tied)
    // When a child has the same score as its parent, the child added 0 derived SNPs
    // at its branch level - choosing the parent is more conservative/accurate
    scores.toList.sortBy(r => (-r.score, r.depth))
  }

  /**
   * Recursively calculate haplogroup scores, descending through the tree.
   *
   * @param consecutiveAncestralBranches Number of consecutive branches where all calls were ancestral.
   *                                     Stop descent when this reaches 2.
   */
  private def calculateHaplogroupScore(
                                        haplogroup: Haplogroup,
                                        snpCalls: Map[Long, String],
                                        scores: mutable.ListBuffer[HaplogroupResult],
                                        cumulativeScore: Double,
                                        cumulativeDerived: Int,
                                        cumulativeAncestral: Int,
                                        cumulativeNoCalls: Int,
                                        cumulativeSnps: Int,
                                        depth: Int,
                                        consecutiveAncestralBranches: Int
                                      ): Unit = {

    var branchDerived = 0
    var branchAncestral = 0
    var branchNoCalls = 0

    for (locus <- haplogroup.loci) {
      snpCalls.get(locus.position) match {
        case Some(calledBase) =>
          // Use case-insensitive comparison - FTDNA tree and VCF can have mixed case bases
          if (calledBase.equalsIgnoreCase(locus.alt)) {
            branchDerived += 1
          } else if (calledBase.equalsIgnoreCase(locus.ref)) {
            branchAncestral += 1
          }
          // If calledBase doesn't match either ref or alt, treat as no-call
        case None =>
          branchNoCalls += 1
      }
    }

    // Calculate branch score that rewards starburst patterns
    // A branch with 50 SNPs where 48 are derived should score much higher than
    // a branch with 2 SNPs where 2 are derived
    val branchCallable = branchDerived + branchAncestral
    val branchScore = if (branchCallable > 0) {
      val matchRate = branchDerived.toDouble / branchCallable.toDouble
      // Score = matchRate * callable - (1 - matchRate) * callable
      //       = (2 * matchRate - 1) * callable
      // This gives: 100% match rate on 50 SNPs = +50
      //             50% match rate on 50 SNPs = 0
      //             0% match rate on 50 SNPs = -50
      (2.0 * matchRate - 1.0) * branchCallable
    } else {
      0.0 // No callable SNPs = neutral
    }

    // Update cumulative values
    val newCumulativeScore = cumulativeScore + branchScore
    val newCumulativeDerived = cumulativeDerived + branchDerived
    val newCumulativeAncestral = cumulativeAncestral + branchAncestral
    val newCumulativeNoCalls = cumulativeNoCalls + branchNoCalls
    val newCumulativeSnps = cumulativeSnps + haplogroup.loci.length

    scores += HaplogroupResult(
      name = haplogroup.name,
      score = newCumulativeScore,
      matchingSnps = newCumulativeDerived,
      mismatchingSnps = 0,
      ancestralMatches = newCumulativeAncestral,
      noCalls = newCumulativeNoCalls,
      totalSnps = haplogroup.loci.length,
      cumulativeSnps = newCumulativeSnps,
      depth = depth
    )

    // Determine if this branch has evidence that says to STOP descent
    // - All ancestral (with calls): we have data showing this is the wrong path
    // - No-calls alone should NOT stop descent - absence of data isn't evidence against
    // The "2 consecutive" rule provides slack for a single branch with poor data hygiene
    // (e.g., one reversed allele due to tree error) not stopping valid descent
    val branchHasCalls = branchDerived > 0 || branchAncestral > 0
    val branchAllAncestral = branchHasCalls && branchDerived == 0 && branchAncestral > 0

    val newConsecutiveAncestral = if (branchDerived > 0) {
      // Reset counter if we found any derived calls - this path is confirmed
      0
    } else if (branchAllAncestral) {
      // Increment counter: we have calls but ALL are ancestral (evidence to stop)
      consecutiveAncestralBranches + 1
    } else {
      // No-calls or empty branch - keep counter as is (no evidence to stop)
      consecutiveAncestralBranches
    }

    // Stop descent if we have two consecutive branches with evidence AGAINST this path
    // (all-ancestral calls). This allows one "bad" branch due to tree hygiene issues.
    if (newConsecutiveAncestral >= 2) {
      return
    }

    // Continue descent to children
    for (child <- haplogroup.children) {
      calculateHaplogroupScore(
        child,
        snpCalls,
        scores,
        newCumulativeScore,
        newCumulativeDerived,
        newCumulativeAncestral,
        newCumulativeNoCalls,
        newCumulativeSnps,
        depth + 1,
        newConsecutiveAncestral
      )
    }
  }
}