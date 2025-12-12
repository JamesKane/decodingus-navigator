package com.decodingus.yprofile.concordance

import com.decodingus.yprofile.model.*

/**
 * Concordance calculation for Y chromosome variant calls.
 *
 * Uses quality-weighted voting to determine consensus from multiple test sources.
 * Different quality tiers are used for SNPs/INDELs vs STRs based on method reliability.
 *
 * Method weights are defined in YProfileSourceType as data-carrying enum values.
 */
object YVariantConcordance:

  // Callable state weights are now defined in YCallableState enum values.

  /**
   * Threshold for discordance to be classified as CONFLICT.
   * If more than 30% of weighted votes disagree, it's a conflict.
   */
  val ConflictThreshold = 0.30

  /**
   * Minimum confidence score for CONFIRMED status.
   */
  val ConfirmationThreshold = 0.70

  /**
   * Calculate concordance weight for a source call.
   *
   * Formula: methodWeight × (1 + depthBonus) × mapQFactor × callableFactor
   *
   * Method weights are sourced from YProfileSourceType enum values.
   *
   * @param sourceType     Testing method used
   * @param variantType    Type of variant (affects method weight selection)
   * @param readDepth      Read depth at position (optional, for sequencing sources)
   * @param mappingQuality Mapping quality (optional, for sequencing sources)
   * @param callableState  Callable state at position (optional)
   * @return Calculated concordance weight
   */
  def calculateWeight(
    sourceType: YProfileSourceType,
    variantType: YVariantType,
    readDepth: Option[Int] = None,
    mappingQuality: Option[Double] = None,
    callableState: Option[YCallableState] = None
  ): Double =
    // Select method weight based on variant type (from enum)
    val methodWeight = variantType match
      case YVariantType.STR => sourceType.strWeight
      case _                => sourceType.snpWeight

    // Depth bonus: min(sqrt(depth)/10, 1.0) - rewards higher coverage
    val depthBonus = readDepth match
      case Some(depth) if depth > 0 => math.min(math.sqrt(depth.toDouble) / 10.0, 1.0)
      case _                        => 0.0

    // Mapping quality factor: min(MQ/60, 1.0)
    val mapQFactor = mappingQuality match
      case Some(mq) if mq > 0 => math.min(mq / 60.0, 1.0)
      case _                  => 1.0 // Default to 1.0 for non-sequencing sources

    // Callable state factor (from enum)
    val callableFactor = callableState match
      case Some(state) => state.weight
      case None        => 1.0 // Assume callable if unknown

    // Final weight calculation
    methodWeight * (1.0 + depthBonus) * mapQFactor * callableFactor

  /**
   * Input data for a single source call.
   */
  case class SourceCallInput(
    sourceId: java.util.UUID,
    sourceType: YProfileSourceType,
    calledAllele: String,
    callState: YConsensusState,
    readDepth: Option[Int] = None,
    mappingQuality: Option[Double] = None,
    callableState: Option[YCallableState] = None,
    calledRepeatCount: Option[Int] = None
  )

  /**
   * Result of concordance calculation.
   */
  case class ConcordanceResult(
    consensusAllele: Option[String],
    consensusState: YConsensusState,
    status: YVariantStatus,
    confidenceScore: Double,
    sourceCount: Int,
    concordantCount: Int,
    discordantCount: Int,
    weightedCalls: List[(SourceCallInput, Double)]  // Calls with their calculated weights
  )

  /**
   * Calculate consensus from multiple source calls.
   *
   * @param calls       List of source calls
   * @param variantType Type of variant (affects weight calculation)
   * @param isInTree    Whether the variant is in the reference haplogroup tree
   * @return Concordance result with consensus and status
   */
  def calculateConsensus(
    calls: List[SourceCallInput],
    variantType: YVariantType,
    isInTree: Boolean = true
  ): ConcordanceResult =
    if calls.isEmpty then
      return ConcordanceResult(
        consensusAllele = None,
        consensusState = YConsensusState.NO_CALL,
        status = YVariantStatus.NO_COVERAGE,
        confidenceScore = 0.0,
        sourceCount = 0,
        concordantCount = 0,
        discordantCount = 0,
        weightedCalls = List.empty
      )

    // Calculate weights for each call
    val weightedCalls = calls.map { call =>
      val weight = calculateWeight(
        call.sourceType,
        variantType,
        call.readDepth,
        call.mappingQuality,
        call.callableState
      )
      (call, weight)
    }

    // Filter out NO_CALL states for voting (but keep them in weighted calls)
    val votingCalls = weightedCalls.filter(_._1.callState != YConsensusState.NO_CALL)

    if votingCalls.isEmpty then
      return ConcordanceResult(
        consensusAllele = None,
        consensusState = YConsensusState.NO_CALL,
        status = YVariantStatus.NO_COVERAGE,
        confidenceScore = 0.0,
        sourceCount = calls.size,
        concordantCount = 0,
        discordantCount = 0,
        weightedCalls = weightedCalls
      )

    // Group by allele and sum weights
    val alleleWeights = votingCalls
      .groupBy(_._1.calledAllele)
      .view
      .mapValues(_.map(_._2).sum)
      .toMap

    val totalWeight = alleleWeights.values.sum
    val (consensusAllele, consensusWeight) = alleleWeights.maxBy(_._2)

    // Determine consensus state from the winning allele's calls
    val winningCalls = votingCalls.filter(_._1.calledAllele == consensusAllele)
    val consensusState = winningCalls.map(_._1.callState).groupBy(identity).maxBy(_._2.size)._1

    // Calculate concordance metrics
    val concordantWeight = consensusWeight
    val discordantWeight = totalWeight - consensusWeight

    val concordantCount = votingCalls.count(_._1.calledAllele == consensusAllele)
    val discordantCount = votingCalls.size - concordantCount

    // Confidence score: concordant weight / total weight
    val confidenceScore = if totalWeight > 0 then concordantWeight / totalWeight else 0.0

    // Determine status
    val discordanceRatio = if totalWeight > 0 then discordantWeight / totalWeight else 0.0
    val status = determineStatus(
      discordanceRatio = discordanceRatio,
      confidenceScore = confidenceScore,
      isInTree = isInTree,
      hasData = votingCalls.nonEmpty
    )

    ConcordanceResult(
      consensusAllele = Some(consensusAllele),
      consensusState = consensusState,
      status = status,
      confidenceScore = confidenceScore,
      sourceCount = calls.size,
      concordantCount = concordantCount,
      discordantCount = discordantCount,
      weightedCalls = weightedCalls
    )

  /**
   * Calculate consensus for STR markers with repeat count comparison.
   *
   * For STRs, differences of more than 1 repeat are considered significant conflicts.
   *
   * @param calls       List of source calls with repeat counts
   * @param isInTree    Whether the marker is in the reference panel
   * @return Concordance result
   */
  def calculateStrConsensus(
    calls: List[SourceCallInput],
    isInTree: Boolean = true
  ): ConcordanceResult =
    if calls.isEmpty then
      return ConcordanceResult(
        consensusAllele = None,
        consensusState = YConsensusState.NO_CALL,
        status = YVariantStatus.NO_COVERAGE,
        confidenceScore = 0.0,
        sourceCount = 0,
        concordantCount = 0,
        discordantCount = 0,
        weightedCalls = List.empty
      )

    // Calculate weights for each call
    val weightedCalls = calls.map { call =>
      val weight = calculateWeight(
        call.sourceType,
        YVariantType.STR,
        call.readDepth,
        call.mappingQuality,
        call.callableState
      )
      (call, weight)
    }

    // Filter out NO_CALL states and calls without repeat counts
    val votingCalls = weightedCalls.filter { case (call, _) =>
      call.callState != YConsensusState.NO_CALL && call.calledRepeatCount.isDefined
    }

    if votingCalls.isEmpty then
      return ConcordanceResult(
        consensusAllele = None,
        consensusState = YConsensusState.NO_CALL,
        status = YVariantStatus.NO_COVERAGE,
        confidenceScore = 0.0,
        sourceCount = calls.size,
        concordantCount = 0,
        discordantCount = 0,
        weightedCalls = weightedCalls
      )

    // Group by repeat count and sum weights
    val repeatCountWeights = votingCalls
      .groupBy(_._1.calledRepeatCount.get)
      .view
      .mapValues(_.map(_._2).sum)
      .toMap

    val totalWeight = repeatCountWeights.values.sum
    val (consensusRepeatCount, consensusWeight) = repeatCountWeights.maxBy(_._2)

    // Calculate concordance with tolerance for off-by-one differences
    // Calls within 1 repeat of consensus are considered concordant
    val (concordantWeight, discordantWeight) = votingCalls.foldLeft((0.0, 0.0)) { case ((conc, disc), (call, weight)) =>
      val repeatDiff = math.abs(call.calledRepeatCount.get - consensusRepeatCount)
      if repeatDiff <= 1 then (conc + weight, disc)
      else (conc, disc + weight)
    }

    val concordantCount = votingCalls.count { case (call, _) =>
      math.abs(call.calledRepeatCount.get - consensusRepeatCount) <= 1
    }
    val discordantCount = votingCalls.size - concordantCount

    val confidenceScore = if totalWeight > 0 then concordantWeight / totalWeight else 0.0
    val discordanceRatio = if totalWeight > 0 then discordantWeight / totalWeight else 0.0

    // For STRs, use the allele string from the winning repeat count
    val winningCall = votingCalls.filter(_._1.calledRepeatCount.contains(consensusRepeatCount)).maxBy(_._2)
    val consensusAllele = winningCall._1.calledAllele

    val status = determineStatus(
      discordanceRatio = discordanceRatio,
      confidenceScore = confidenceScore,
      isInTree = isInTree,
      hasData = votingCalls.nonEmpty
    )

    ConcordanceResult(
      consensusAllele = Some(consensusAllele),
      consensusState = YConsensusState.DERIVED,  // STRs are typically considered derived
      status = status,
      confidenceScore = confidenceScore,
      sourceCount = calls.size,
      concordantCount = concordantCount,
      discordantCount = discordantCount,
      weightedCalls = weightedCalls
    )

  /**
   * Determine variant status based on concordance metrics.
   */
  private def determineStatus(
    discordanceRatio: Double,
    confidenceScore: Double,
    isInTree: Boolean,
    hasData: Boolean
  ): YVariantStatus =
    if !hasData then
      YVariantStatus.NO_COVERAGE
    else if discordanceRatio > ConflictThreshold then
      YVariantStatus.CONFLICT
    else if confidenceScore >= ConfirmationThreshold then
      if isInTree then YVariantStatus.CONFIRMED
      else YVariantStatus.NOVEL
    else
      YVariantStatus.PENDING

  /**
   * Check if a variant should be flagged as heteroplasmy.
   *
   * Heteroplasmy on the Y chromosome is rare but can occur in STRs
   * and certain repetitive regions.
   *
   * @param variantAlleleFrequency VAF from sequencing
   * @return true if VAF suggests heteroplasmy
   */
  def isLikelyHeteroplasmy(variantAlleleFrequency: Option[Double]): Boolean =
    variantAlleleFrequency.exists(vaf => vaf >= 0.15 && vaf <= 0.85)

  /**
   * Calculate overall profile confidence from variant statistics.
   *
   * @param confirmedCount Number of confirmed variants
   * @param novelCount     Number of novel variants
   * @param conflictCount  Number of conflicting variants
   * @param totalCount     Total variant count
   * @return Overall confidence score (0.0 to 1.0)
   */
  def calculateProfileConfidence(
    confirmedCount: Int,
    novelCount: Int,
    conflictCount: Int,
    totalCount: Int
  ): Double =
    if totalCount == 0 then 0.0
    else
      // Weight: confirmed contributes fully, novel partially, conflicts reduce confidence
      val positiveContribution = confirmedCount + (0.7 * novelCount)
      val negativeContribution = conflictCount * 0.5
      val score = (positiveContribution - negativeContribution) / totalCount
      math.max(0.0, math.min(1.0, score))
