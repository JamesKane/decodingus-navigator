package com.decodingus.ancestry.processor

import com.decodingus.ancestry.model.*
import com.decodingus.config.FeatureToggles

/**
 * Core ancestry estimation algorithm using PCA projection and Gaussian mixture model.
 *
 * Algorithm overview:
 * 1. Project sample genotypes onto pre-computed PCA space from reference populations
 * 2. Calculate Mahalanobis distance to each population centroid
 * 3. Convert distances to probabilities using multivariate Gaussian PDF
 * 4. Normalize to get percentage assignments
 */
class AncestryEstimator {

  /**
   * Estimate ancestry proportions for a sample.
   *
   * @param genotypes Map of SNP ID (chr:pos) to genotype (0, 1, 2, or -1 for missing)
   * @param alleleFreqs Reference allele frequency matrix
   * @param pcaLoadings PCA loadings for projection
   * @param panelType Panel type identifier ("aims" or "genome-wide")
   * @return AncestryResult with population percentages
   */
  def estimate(
    genotypes: Map[String, Int],
    alleleFreqs: AlleleFrequencyMatrix,
    pcaLoadings: PCALoadings,
    panelType: String
  ): AncestryResult = {

    // Step 1: Project sample onto PCA space
    val pcaCoords = projectToPca(genotypes, pcaLoadings)

    // Step 2: Calculate probability for each population
    val populationProbs = calculatePopulationProbabilities(pcaCoords, pcaLoadings)

    // Step 3: Calculate confidence based on data quality
    val snpsWithData = genotypes.count(_._2 >= 0)
    val totalSnps = pcaLoadings.numSnps
    val confidence = calculateConfidence(snpsWithData, totalSnps, pcaCoords, pcaLoadings)

    // Step 4: Build result
    AncestryResult.fromProbabilities(
      panelType = panelType,
      snpsAnalyzed = totalSnps,
      snpsWithGenotype = snpsWithData,
      populationProbs = populationProbs,
      confidenceLevel = confidence,
      analysisVersion = "1.0.0",
      referenceVersion = FeatureToggles.ancestryAnalysis.referenceVersion,
      pcaCoords = Some(pcaCoords.take(3).toList) // First 3 PCs for visualization
    )
  }

  /**
   * Project sample genotypes onto reference PCA space.
   *
   * For each SNP with data:
   * 1. Center the genotype by subtracting the population mean
   * 2. Multiply by the PCA loading for each component
   * 3. Sum contributions across all SNPs
   */
  private def projectToPca(
    genotypes: Map[String, Int],
    loadings: PCALoadings
  ): Array[Double] = {
    val numComponents = loadings.numComponents
    val coords = new Array[Double](numComponents)
    var snpsUsed = 0

    // Create SNP ID to index mapping for efficient lookup
    val snpIdToIndex = loadings.snpIds.zipWithIndex.toMap

    genotypes.foreach { case (snpId, geno) =>
      if (geno >= 0) { // Only use non-missing genotypes
        snpIdToIndex.get(snpId).foreach { snpIdx =>
          val mean = loadings.snpMeans(snpIdx)
          val centered = geno.toDouble - mean

          // Add contribution to each PC
          for (pc <- 0 until numComponents) {
            coords(pc) += centered * loadings.getLoading(snpIdx, pc)
          }
          snpsUsed += 1
        }
      }
    }

    // Normalize by number of SNPs used (optional, depends on how loadings were computed)
    // For now, keep raw projection values
    coords
  }

  /**
   * Calculate probability of belonging to each population.
   *
   * Uses multivariate Gaussian probability density function:
   * P(x|pop) = exp(-0.5 * mahalanobis_distance^2)
   *
   * Returns unnormalized probabilities (caller should normalize to sum=1).
   */
  private def calculatePopulationProbabilities(
    sampleCoords: Array[Double],
    loadings: PCALoadings
  ): Map[String, Double] = {
    loadings.populations.zipWithIndex.map { case (popCode, popIdx) =>
      val centroid = loadings.getCentroid(popIdx)
      val variance = loadings.getVariance(popIdx)

      val mahalanobis = computeMahalanobisDistance(sampleCoords, centroid, variance)
      val prob = math.exp(-0.5 * mahalanobis)

      popCode -> prob
    }.toMap
  }

  /**
   * Compute squared Mahalanobis distance using diagonal covariance.
   *
   * d^2 = sum((x_i - mu_i)^2 / sigma_i^2)
   */
  private def computeMahalanobisDistance(
    sample: Array[Double],
    centroid: Array[Float],
    variance: Array[Float]
  ): Double = {
    sample.indices.map { i =>
      val diff = sample(i) - centroid(i)
      val v = variance(i).toDouble
      if (v > 1e-10) {
        (diff * diff) / v
      } else {
        0.0 // Avoid division by zero
      }
    }.sum
  }

  /**
   * Calculate overall confidence score based on data quality.
   *
   * Factors:
   * - Data completeness (fraction of SNPs with genotypes)
   * - Distinctiveness (how clearly sample clusters with one population)
   */
  private def calculateConfidence(
    snpsWithData: Int,
    totalSnps: Int,
    pcaCoords: Array[Double],
    loadings: PCALoadings
  ): Double = {
    // Base confidence from data completeness
    val completeness = snpsWithData.toDouble / totalSnps

    // Adjust for very low completeness
    val adjustedCompleteness = if (completeness < 0.5) {
      completeness * 0.5 // Penalize heavily
    } else {
      0.25 + completeness * 0.75 // Scale 0.5-1.0 to 0.625-1.0
    }

    // Could add distinctiveness factor based on how close sample is to any centroid
    // For now, just use completeness
    math.min(1.0, adjustedCompleteness)
  }

  /**
   * Estimate ancestry using a simpler allele frequency matching approach.
   * Fallback for when PCA projection is not suitable.
   *
   * For each population, calculates likelihood as:
   * Product of P(genotype | allele_freq) across all SNPs
   */
  def estimateByAlleleFrequency(
    genotypes: Map[String, Int],
    alleleFreqs: AlleleFrequencyMatrix,
    panelType: String
  ): AncestryResult = {
    val snpIdToIndex = alleleFreqs.snpIds.zipWithIndex.toMap
    val numPops = alleleFreqs.numPopulations

    // Calculate log-likelihood for each population
    val logLikelihoods = new Array[Double](numPops)

    genotypes.foreach { case (snpId, geno) =>
      if (geno >= 0) {
        snpIdToIndex.get(snpId).foreach { snpIdx =>
          for (popIdx <- 0 until numPops) {
            val p = alleleFreqs.getFrequency(popIdx, snpIdx).toDouble
            // Clamp to avoid log(0)
            val freq = math.max(0.001, math.min(0.999, p))

            // Binomial probability for diploid genotype
            val prob = geno match {
              case 0 => (1 - freq) * (1 - freq)           // AA
              case 1 => 2 * freq * (1 - freq)             // Aa
              case 2 => freq * freq                       // aa
              case _ => 1.0
            }
            logLikelihoods(popIdx) += math.log(prob)
          }
        }
      }
    }

    // Convert log-likelihoods to probabilities
    val maxLL = logLikelihoods.max
    val probs = logLikelihoods.map(ll => math.exp(ll - maxLL))
    val popProbs = alleleFreqs.populations.zip(probs).toMap

    val snpsWithData = genotypes.count(_._2 >= 0)
    val confidence = snpsWithData.toDouble / alleleFreqs.numSnps

    AncestryResult.fromProbabilities(
      panelType = panelType,
      snpsAnalyzed = alleleFreqs.numSnps,
      snpsWithGenotype = snpsWithData,
      populationProbs = popProbs,
      confidenceLevel = math.min(1.0, confidence),
      analysisVersion = "1.0.0-af",
      referenceVersion = FeatureToggles.ancestryAnalysis.referenceVersion,
      pcaCoords = None
    )
  }
}
