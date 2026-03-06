package com.decodingus.haplogroup.scoring

import com.decodingus.haplogroup.model.HaplogroupResult

/**
 * Shared confidence calculation for haplogroup assignments.
 *
 * Confidence is based on:
 * - Match quality: proportion of callable SNPs that are derived (matching the path)
 * - Ambiguity penalty: reduction if sibling/cousin haplogroups have similar scores
 *
 * We don't penalize for missing tree coverage - those are positions we simply
 * can't evaluate, not evidence against the call.
 */
object ConfidenceCalculator {

  /**
   * Calculate confidence score for a haplogroup assignment.
   *
   * @param topResult  The top-scoring haplogroup result
   * @param allResults All scored results (sorted by score descending)
   * @param maxCap     Maximum confidence cap (e.g., 0.85 for chip data, 1.0 for WGS)
   * @return Confidence score between 0.0 and maxCap
   */
  def calculate(topResult: HaplogroupResult, allResults: List[HaplogroupResult], maxCap: Double = 1.0): Double = {
    // Match quality: proportion of callable SNPs that are derived (matching the path)
    val callableSnps = topResult.matchingSnps + topResult.ancestralMatches
    val matchQuality = if (callableSnps > 0) {
      topResult.matchingSnps.toDouble / callableSnps
    } else {
      0.0
    }

    // Ambiguity penalty: if a sibling/cousin haplogroup is close to top, we're less certain
    // Don't penalize for ancestors being close - that's expected (parent score + branch = child score)
    val ambiguityPenalty = calculateAmbiguityPenalty(topResult, allResults)

    val confidence = matchQuality * (1.0 - ambiguityPenalty)
    math.min(maxCap, math.max(0.0, confidence))
  }

  /**
   * Calculate ambiguity penalty based on competing (non-ancestor) haplogroups.
   */
  private def calculateAmbiguityPenalty(topResult: HaplogroupResult, allResults: List[HaplogroupResult]): Double = {
    if (allResults.size <= 1 || topResult.score <= 0) {
      return 0.0
    }

    val topPath = topResult.lineagePath

    // Find closest non-ancestor competitor
    val closestCompetitor = allResults.tail.find { candidate =>
      val candidatePath = candidate.lineagePath
      // Not an ancestor if candidate path is NOT a prefix of top path
      !topPath.startsWith(candidatePath)
    }

    closestCompetitor match {
      case Some(competitor) if competitor.score > 0 =>
        val scoreDiff = (topResult.score - competitor.score) / topResult.score
        // If competitor is within 20% of top score, apply penalty proportional to closeness
        if (scoreDiff < 0.2) {
          (0.2 - scoreDiff) * 0.5 // Up to 10% penalty when scores are nearly identical
        } else {
          0.0
        }
      case _ =>
        0.0 // No non-ancestor competitors
    }
  }
}
