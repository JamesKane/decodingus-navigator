package com.decodingus.genotype.processor

import com.decodingus.ancestry.model.*
import com.decodingus.ancestry.processor.AncestryEstimator
import com.decodingus.ancestry.reference.{AncestryReferenceGateway, AncestryReferenceResult}
import com.decodingus.genotype.model.GenotypeCall

/**
 * Adapter to run ancestry analysis on chip/array genotype data.
 *
 * This bridges the chip data processing system with the ancestry
 * analysis system, allowing population percentage estimation from
 * SNP array data (23andMe, AncestryDNA, etc.)
 *
 * Unlike WGS-based ancestry analysis which uses GATK HaplotypeCaller,
 * chip data already has genotypes called at known positions, so we
 * can directly project onto PCA space.
 */
class ChipAncestryAdapter {

  private val estimator = new AncestryEstimator()

  /**
   * Run ancestry analysis on chip genotypes.
   *
   * @param chipResult The processed chip data
   * @param panelType The ancestry panel to use (AIMs or GenomeWide)
   * @param onProgress Progress callback
   * @return Either error or ancestry result
   */
  def analyze(
    chipResult: ChipProcessingResult,
    panelType: AncestryPanelType,
    onProgress: (String, Double, Double) => Unit
  ): Either[String, AncestryResult] = {

    val panelName = panelType match {
      case AncestryPanelType.Aims => "aims"
      case AncestryPanelType.GenomeWide => "genome-wide"
    }

    onProgress(s"Checking ancestry reference data for $panelName panel...", 0.0, 1.0)

    // Check if reference data is available
    // Note: For chip data, we use GRCh37 coordinates as most chips are on that build
    val referenceBuild = "GRCh37" // Most chip data uses GRCh37

    val ancestryGateway = new AncestryReferenceGateway((_, _) => {})

    ancestryGateway.checkAvailability(panelType, referenceBuild) match {
      case AncestryReferenceResult.Available(_, alleleFreqs, pcaLoadings) =>
        analyzeWithData(chipResult, alleleFreqs, pcaLoadings, panelName, onProgress)

      case AncestryReferenceResult.DownloadRequired(_, _, sizeMB, _) =>
        Left(s"Ancestry reference panel not available. Download required (~${sizeMB}MB). " +
          "Please run the download first.")

      case AncestryReferenceResult.Error(msg) =>
        Left(s"Failed to access ancestry reference data: $msg")
    }
  }

  /**
   * Analyze chip data using pre-loaded reference data.
   */
  private def analyzeWithData(
    chipResult: ChipProcessingResult,
    alleleFreqs: AlleleFrequencyMatrix,
    pcaLoadings: PCALoadings,
    panelName: String,
    onProgress: (String, Double, Double) => Unit
  ): Either[String, AncestryResult] = {

    onProgress("Converting chip genotypes to ancestry format...", 0.2, 1.0)

    // Build a map of SNP positions covered by the ancestry panel
    val panelSnpIds = pcaLoadings.snpIds.toSet

    // Convert chip calls to genotype map (chr:pos -> 0/1/2/-1)
    // We need to match chip positions to panel positions
    val genotypeMap = buildGenotypeMap(chipResult.autosomalCalls, panelSnpIds, alleleFreqs)

    val snpsMatched = genotypeMap.count(_._2 >= 0)

    onProgress(s"Matched $snpsMatched SNPs with ancestry panel.", 0.4, 1.0)

    // Check minimum SNP requirements
    val minSnps = panelName match {
      case "aims" => 1000 // Lower threshold for chip data on AIMs
      case "genome-wide" => 50000
      case _ => 1000
    }

    if (snpsMatched < minSnps) {
      Left(s"Insufficient overlap: $snpsMatched SNPs matched with ancestry panel, " +
        s"minimum required is $minSnps. The chip may not cover enough ancestry-informative markers.")
    } else {
      onProgress("Estimating ancestry proportions...", 0.6, 1.0)

      // Run the PCA projection + GMM estimation
      val result = estimator.estimate(genotypeMap, alleleFreqs, pcaLoadings, panelName)

      onProgress("Ancestry analysis complete.", 1.0, 1.0)

      // Adjust confidence based on chip data limitations
      val adjustedResult = adjustConfidenceForChipData(result, snpsMatched, pcaLoadings.numSnps)

      Right(adjustedResult)
    }
  }

  /**
   * Build genotype map from chip calls matched to panel SNP positions.
   *
   * The challenge is that chip data uses rsIDs while our panel uses chr:pos format.
   * We'll match by position where possible.
   */
  private def buildGenotypeMap(
    autosomalCalls: List[GenotypeCall],
    panelSnpIds: Set[String],
    alleleFreqs: AlleleFrequencyMatrix
  ): Map[String, Int] = {

    // Create position-based lookup
    val callsByPosition: Map[String, GenotypeCall] = autosomalCalls.map { call =>
      val snpId = s"${normalizeChromosome(call.chromosome)}:${call.position}"
      snpId -> call
    }.toMap

    // For each panel SNP, try to find matching chip call
    panelSnpIds.flatMap { snpId =>
      callsByPosition.get(snpId).map { call =>
        // Determine reference allele from allele frequency matrix
        // The panel stores minor allele frequencies, so we need to determine ref/alt
        val snpIdx = alleleFreqs.snpIds.indexOf(snpId)
        if (snpIdx >= 0) {
          // Use mean frequency across populations to infer reference
          // Allele with higher frequency is typically reference
          val meanFreq = (0 until alleleFreqs.numPopulations)
            .map(p => alleleFreqs.getFrequency(p, snpIdx))
            .sum / alleleFreqs.numPopulations

          // If MAF is stored, we use the more common allele as reference
          val refAllele = if (meanFreq < 0.5) call.allele1 else call.allele2
          val genotype = call.numericGenotype(refAllele)
          Some(snpId -> genotype)
        } else {
          None
        }
      }.flatten
    }.toMap
  }

  /**
   * Adjust confidence score to account for chip data limitations.
   *
   * Chip data typically has:
   * - Fewer SNPs than WGS
   * - No depth/quality information
   * - Fixed marker set that may not perfectly overlap ancestry panel
   */
  private def adjustConfidenceForChipData(
    result: AncestryResult,
    matchedSnps: Int,
    panelSnps: Int
  ): AncestryResult = {
    val coverageRatio = matchedSnps.toDouble / panelSnps

    // Apply coverage-based penalty
    val adjustedConfidence = result.confidenceLevel * math.min(1.0, coverageRatio * 1.2)

    // Also widen confidence intervals for chip data
    val adjustedPercentages = result.percentages.map { p =>
      val intervalWidth = p.confidenceHigh - p.confidenceLow
      val adjustedWidth = intervalWidth * (1.0 + (1.0 - coverageRatio) * 0.5)
      val midpoint = (p.confidenceHigh + p.confidenceLow) / 2
      p.copy(
        confidenceLow = math.max(0.0, midpoint - adjustedWidth / 2),
        confidenceHigh = math.min(100.0, midpoint + adjustedWidth / 2)
      )
    }

    result.copy(
      confidenceLevel = adjustedConfidence,
      percentages = adjustedPercentages
    )
  }

  /**
   * Normalize chromosome name.
   */
  private def normalizeChromosome(chr: String): String = {
    chr.toLowerCase.stripPrefix("chr") match {
      case n if n.toIntOption.isDefined => n
      case "x" => "X"
      case "y" => "Y"
      case "m" | "mt" => "MT"
      case other => other
    }
  }
}

object ChipAncestryAdapter {
  /**
   * Check if chip data is suitable for ancestry analysis.
   */
  def isSuitableForAncestry(result: ChipProcessingResult): Boolean = {
    result.summary.isAcceptableForAncestry
  }

  /**
   * Get recommended panel type based on chip data quality.
   */
  def recommendedPanelType(result: ChipProcessingResult): AncestryPanelType = {
    // For most chip data, AIMs panel is recommended as it has better
    // overlap with typical chip marker sets
    if (result.summary.autosomalMarkersCalled >= 500000) {
      AncestryPanelType.GenomeWide
    } else {
      AncestryPanelType.Aims
    }
  }
}
