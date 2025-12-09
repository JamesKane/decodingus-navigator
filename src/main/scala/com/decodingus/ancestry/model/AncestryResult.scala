package com.decodingus.ancestry.model

import io.circe.Codec

/**
 * Percentage assignment for a single population.
 */
case class PopulationPercentage(
  populationCode: String,
  populationName: String,
  percentage: Double,                    // 0.0 to 100.0
  confidenceLow: Double,                 // Lower bound of 95% CI
  confidenceHigh: Double,                // Upper bound of 95% CI
  rank: Int                              // Sorted by percentage (1 = highest)
) derives Codec.AsObject

/**
 * Summary for a super-population (continental grouping).
 */
case class SuperPopulationPercentage(
  superPopulation: String,               // "European", "African", etc.
  percentage: Double,                    // Sum of constituent populations
  populations: List[String]              // Contributing population codes
) derives Codec.AsObject

/**
 * Result of ancestry analysis for a single sample.
 */
case class AncestryResult(
  panelType: String,                     // "aims" or "genome-wide"
  snpsAnalyzed: Int,                     // Total SNPs in panel
  snpsWithGenotype: Int,                 // SNPs with valid genotype calls
  snpsMissing: Int,                      // SNPs with no call
  percentages: List[PopulationPercentage],
  superPopulationSummary: List[SuperPopulationPercentage],
  confidenceLevel: Double,               // Overall confidence (0-1) based on data quality
  analysisVersion: String,
  referenceVersion: String,
  pcaCoordinates: Option[List[Double]]   // Optional: first N PCA coordinates for visualization
) derives Codec.AsObject

object AncestryResult {

  /**
   * Create an AncestryResult from raw population probabilities.
   *
   * @param panelType "aims" or "genome-wide"
   * @param snpsAnalyzed Total SNPs in the panel
   * @param snpsWithGenotype SNPs with valid calls
   * @param populationProbs Map of population code -> raw probability
   * @param confidenceLevel Overall confidence score
   * @param analysisVersion Version of the analysis algorithm
   * @param referenceVersion Version of the reference panel
   * @param pcaCoords Optional PCA coordinates
   */
  def fromProbabilities(
    panelType: String,
    snpsAnalyzed: Int,
    snpsWithGenotype: Int,
    populationProbs: Map[String, Double],
    confidenceLevel: Double,
    analysisVersion: String,
    referenceVersion: String,
    pcaCoords: Option[List[Double]] = None
  ): AncestryResult = {
    // Normalize probabilities to percentages
    val totalProb = populationProbs.values.sum
    val normalizedPcts = if (totalProb > 0) {
      populationProbs.map { case (code, prob) => code -> (prob / totalProb * 100.0) }
    } else {
      populationProbs.map { case (code, _) => code -> 0.0 }
    }

    // Sort by percentage and assign ranks
    val sortedPops = normalizedPcts.toList
      .sortBy(-_._2)
      .zipWithIndex
      .flatMap { case ((code, pct), idx) =>
        Population.byCode(code).map { pop =>
          // Simple confidence interval based on sample size and data completeness
          val ciWidth = calculateCiWidth(pct, snpsWithGenotype, snpsAnalyzed)
          PopulationPercentage(
            populationCode = code,
            populationName = pop.name,
            percentage = pct,
            confidenceLow = math.max(0.0, pct - ciWidth),
            confidenceHigh = math.min(100.0, pct + ciWidth),
            rank = idx + 1
          )
        }
      }

    // Calculate super-population summaries
    val superPopSummary = Population.SuperPopulations.map { case (superPop, pops) =>
      val popCodes = pops.map(_.code)
      val total = popCodes.flatMap(normalizedPcts.get).sum
      SuperPopulationPercentage(superPop, total, popCodes)
    }.toList.sortBy(-_.percentage)

    AncestryResult(
      panelType = panelType,
      snpsAnalyzed = snpsAnalyzed,
      snpsWithGenotype = snpsWithGenotype,
      snpsMissing = snpsAnalyzed - snpsWithGenotype,
      percentages = sortedPops,
      superPopulationSummary = superPopSummary,
      confidenceLevel = confidenceLevel,
      analysisVersion = analysisVersion,
      referenceVersion = referenceVersion,
      pcaCoordinates = pcaCoords
    )
  }

  /**
   * Calculate confidence interval width based on percentage and data quality.
   * Uses a simplified approximation based on binomial proportion CI.
   */
  private def calculateCiWidth(pct: Double, snpsWithData: Int, totalSnps: Int): Double = {
    val completeness = snpsWithData.toDouble / totalSnps
    // Base width scales with sqrt(p*(1-p)/n) approximation
    val p = pct / 100.0
    val baseWidth = if (snpsWithData > 0) {
      1.96 * math.sqrt(p * (1 - p) / snpsWithData) * 100.0
    } else {
      50.0 // Maximum uncertainty
    }
    // Increase width for low completeness
    baseWidth / math.max(0.5, completeness)
  }

  /**
   * Filter results to only show populations above a threshold.
   */
  def filterByThreshold(result: AncestryResult, minPercentage: Double): AncestryResult = {
    result.copy(
      percentages = result.percentages.filter(_.percentage >= minPercentage)
    )
  }

  /**
   * Get the primary ancestry (highest super-population).
   */
  def primaryAncestry(result: AncestryResult): Option[String] = {
    result.superPopulationSummary.headOption.map(_.superPopulation)
  }

  /**
   * Check if the sample appears significantly admixed.
   * Returns true if more than one super-population is above threshold.
   */
  def isAdmixed(result: AncestryResult, threshold: Double = 10.0): Boolean = {
    result.superPopulationSummary.count(_.percentage >= threshold) > 1
  }
}
