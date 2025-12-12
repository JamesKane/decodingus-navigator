package com.decodingus.yprofile.service

import com.decodingus.db.Transactor
import com.decodingus.haplogroup.model.HaplogroupResult
import com.decodingus.haplogroup.scoring.HaplogroupScorer
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.repository.{
  YChromosomeProfileRepository, YProfileSourceRepository, YProfileRegionRepository,
  YProfileVariantRepository, YVariantSourceCallRepository, YVariantAuditRepository
}
import com.decodingus.yprofile.concordance.YVariantConcordance
import com.decodingus.yprofile.concordance.YVariantConcordance.{SourceCallInput, ConcordanceResult}
import com.decodingus.yprofile.model.*
import java.time.LocalDateTime
import java.util.UUID

/**
 * Service for managing Y chromosome profiles.
 *
 * Handles profile creation, source import, variant management,
 * and concordance-based reconciliation.
 */
class YProfileService(
  transactor: Transactor,
  profileRepo: YChromosomeProfileRepository,
  sourceRepo: YProfileSourceRepository,
  regionRepo: YProfileRegionRepository,
  variantRepo: YProfileVariantRepository,
  sourceCallRepo: YVariantSourceCallRepository,
  auditRepo: YVariantAuditRepository
):

  // ============================================
  // Profile Management
  // ============================================

  /**
   * Get or create a Y chromosome profile for a biosample.
   */
  def getOrCreateProfile(biosampleId: UUID): Either[String, YChromosomeProfileEntity] =
    transactor.readWrite {
      profileRepo.findByBiosample(biosampleId) match
        case Some(profile) => profile
        case None =>
          profileRepo.insert(YChromosomeProfileEntity.create(biosampleId))
    }

  /**
   * Get a profile by ID.
   */
  def getProfile(profileId: UUID): Either[String, Option[YChromosomeProfileEntity]] =
    transactor.readOnly {
      profileRepo.findById(profileId)
    }

  /**
   * Get a profile by biosample ID.
   */
  def getProfileByBiosample(biosampleId: UUID): Either[String, Option[YChromosomeProfileEntity]] =
    transactor.readOnly {
      profileRepo.findByBiosample(biosampleId)
    }

  /**
   * Delete a profile and all related data.
   */
  def deleteProfile(profileId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      profileRepo.delete(profileId)
    }

  // ============================================
  // Source Management
  // ============================================

  /**
   * Add a new test source to a profile.
   */
  def addSource(
    profileId: UUID,
    sourceType: YProfileSourceType,
    sourceRef: Option[String] = None,
    vendor: Option[String] = None,
    testName: Option[String] = None,
    testDate: Option[LocalDateTime] = None,
    alignmentId: Option[UUID] = None,
    referenceBuild: Option[String] = None,
    meanReadDepth: Option[Double] = None,
    meanMappingQuality: Option[Double] = None,
    coveragePct: Option[Double] = None
  ): Either[String, YProfileSourceEntity] =
    transactor.readWrite {
      val methodTier = YProfileSourceType.snpMethodTier(sourceType)
      val source = YProfileSourceEntity.create(
        yProfileId = profileId,
        sourceType = sourceType,
        sourceRef = sourceRef,
        vendor = vendor,
        testName = testName,
        testDate = testDate,
        methodTier = methodTier,
        meanReadDepth = meanReadDepth,
        meanMappingQuality = meanMappingQuality,
        coveragePct = coveragePct,
        alignmentId = alignmentId,
        referenceBuild = referenceBuild
      )
      val saved = sourceRepo.insert(source)

      // Update profile source count
      updateProfileSourceCount(profileId)

      saved
    }

  /**
   * Get all sources for a profile.
   */
  def getSources(profileId: UUID): Either[String, List[YProfileSourceEntity]] =
    transactor.readOnly {
      sourceRepo.findByProfile(profileId)
    }

  /**
   * Remove a source and its variant calls.
   */
  def removeSource(sourceId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      val source = sourceRepo.findById(sourceId)
      source match
        case Some(s) =>
          // Delete source calls for all variants from this source
          sourceCallRepo.deleteBySource(sourceId)
          // Delete regions
          regionRepo.deleteBySource(sourceId)
          // Delete source
          sourceRepo.delete(sourceId)
          // Update profile
          updateProfileSourceCount(s.yProfileId)
          true
        case None => false
    }

  // ============================================
  // Variant Management
  // ============================================

  /**
   * Add or update a variant with a source call.
   *
   * If the variant already exists (by position+alleles), adds a source call.
   * Otherwise creates a new variant.
   */
  def addVariantCall(
    profileId: UUID,
    sourceId: UUID,
    position: Long,
    refAllele: String,
    altAllele: String,
    calledAllele: String,
    callState: YConsensusState,
    variantType: YVariantType = YVariantType.SNP,
    variantName: Option[String] = None,
    rsId: Option[String] = None,
    markerName: Option[String] = None,
    repeatCount: Option[Int] = None,
    strMetadata: Option[StrMetadata] = None,
    readDepth: Option[Int] = None,
    qualityScore: Option[Double] = None,
    mappingQuality: Option[Double] = None,
    variantAlleleFrequency: Option[Double] = None,
    callableState: Option[YCallableState] = None,
    definingHaplogroup: Option[String] = None,
    haplogroupBranchDepth: Option[Int] = None
  ): Either[String, (YProfileVariantEntity, YVariantSourceCallEntity)] =
    transactor.readWrite {
      // Find or create variant
      val variant = variantRepo.findByPosition(profileId, position) match
        case Some(existing) =>
          // Update variant if we have new info
          val updated = existing.copy(
            variantName = variantName.orElse(existing.variantName),
            rsId = rsId.orElse(existing.rsId),
            markerName = markerName.orElse(existing.markerName),
            definingHaplogroup = definingHaplogroup.orElse(existing.definingHaplogroup),
            haplogroupBranchDepth = haplogroupBranchDepth.orElse(existing.haplogroupBranchDepth)
          )
          if updated != existing then variantRepo.update(updated) else existing
        case None =>
          variantRepo.insert(YProfileVariantEntity.create(
            yProfileId = profileId,
            position = position,
            refAllele = refAllele,
            altAllele = altAllele,
            variantType = variantType,
            variantName = variantName,
            rsId = rsId,
            markerName = markerName,
            repeatCount = repeatCount,
            strMetadata = strMetadata,
            definingHaplogroup = definingHaplogroup,
            haplogroupBranchDepth = haplogroupBranchDepth
          ))

      // Get source info for weight calculation
      val source = sourceRepo.findById(sourceId).getOrElse(
        throw new IllegalArgumentException(s"Source not found: $sourceId")
      )

      // Calculate concordance weight
      val weight = YVariantConcordance.calculateWeight(
        source.sourceType,
        variantType,
        readDepth,
        mappingQuality,
        callableState
      )

      // Add source call (replace if exists)
      sourceCallRepo.findByVariantAndSource(variant.id, sourceId) match
        case Some(existing) =>
          sourceCallRepo.delete(existing.id)
        case None => ()

      val sourceCall = sourceCallRepo.insert(YVariantSourceCallEntity.create(
        variantId = variant.id,
        sourceId = sourceId,
        calledAllele = calledAllele,
        callState = callState,
        calledRepeatCount = repeatCount,
        readDepth = readDepth,
        qualityScore = qualityScore,
        mappingQuality = mappingQuality,
        variantAlleleFrequency = variantAlleleFrequency,
        callableState = callableState,
        concordanceWeight = weight
      ))

      (variant, sourceCall)
    }

  /**
   * Get all variants for a profile.
   */
  def getVariants(profileId: UUID): Either[String, List[YProfileVariantEntity]] =
    transactor.readOnly {
      variantRepo.findByProfile(profileId)
    }

  /**
   * Get source calls for a variant.
   */
  def getVariantCalls(variantId: UUID): Either[String, List[YVariantSourceCallEntity]] =
    transactor.readOnly {
      sourceCallRepo.findByVariant(variantId)
    }

  // ============================================
  // Reconciliation
  // ============================================

  /**
   * Reconcile a single variant using concordance algorithm.
   *
   * @param variantId The variant to reconcile
   * @param isInTree Whether the variant is in the haplogroup reference tree
   * @return Updated variant entity
   */
  def reconcileVariant(variantId: UUID, isInTree: Boolean = true): Either[String, YProfileVariantEntity] =
    transactor.readWrite {
      val variant = variantRepo.findById(variantId).getOrElse(
        throw new IllegalArgumentException(s"Variant not found: $variantId")
      )

      val calls = sourceCallRepo.findByVariant(variantId)
      val sources = calls.flatMap(c => sourceRepo.findById(c.sourceId))
      val sourceMap = sources.map(s => s.id -> s).toMap

      // Convert to concordance input
      val callInputs = calls.flatMap { call =>
        sourceMap.get(call.sourceId).map { source =>
          SourceCallInput(
            sourceId = call.sourceId,
            sourceType = source.sourceType,
            calledAllele = call.calledAllele,
            callState = call.callState,
            readDepth = call.readDepth,
            mappingQuality = call.mappingQuality,
            callableState = call.callableState,
            calledRepeatCount = call.calledRepeatCount
          )
        }
      }

      // Calculate consensus based on variant type
      val result = variant.variantType match
        case YVariantType.STR =>
          YVariantConcordance.calculateStrConsensus(callInputs, isInTree)
        case _ =>
          YVariantConcordance.calculateConsensus(callInputs, variant.variantType, isInTree)

      // Update variant with consensus
      val updated = variant.copy(
        consensusAllele = result.consensusAllele,
        consensusState = result.consensusState,
        status = result.status,
        sourceCount = result.sourceCount,
        concordantCount = result.concordantCount,
        discordantCount = result.discordantCount,
        confidenceScore = result.confidenceScore,
        maxReadDepth = calls.flatMap(_.readDepth).maxOption,
        maxQualityScore = calls.flatMap(_.qualityScore).maxOption,
        lastUpdatedAt = LocalDateTime.now()
      )

      variantRepo.update(updated)
    }

  /**
   * Reconcile all variants in a profile.
   *
   * @param profileId The profile to reconcile
   * @param treeVariantNames Set of variant names that are in the haplogroup tree
   * @return Count of variants reconciled
   */
  def reconcileProfile(
    profileId: UUID,
    treeVariantNames: Set[String] = Set.empty
  ): Either[String, Int] =
    transactor.readWrite {
      val variants = variantRepo.findByProfile(profileId)

      variants.foreach { variant =>
        val isInTree = variant.variantName.exists(treeVariantNames.contains) ||
                       variant.definingHaplogroup.isDefined
        reconcileVariantInternal(variant, isInTree)
      }

      // Update profile statistics
      updateProfileStatistics(profileId)

      // Mark as reconciled
      profileRepo.markReconciled(profileId)

      variants.size
    }

  /**
   * Internal variant reconciliation (within transaction).
   */
  private def reconcileVariantInternal(variant: YProfileVariantEntity, isInTree: Boolean)(using java.sql.Connection): YProfileVariantEntity =
    val calls = sourceCallRepo.findByVariant(variant.id)
    val sources = calls.flatMap(c => sourceRepo.findById(c.sourceId))
    val sourceMap = sources.map(s => s.id -> s).toMap

    val callInputs = calls.flatMap { call =>
      sourceMap.get(call.sourceId).map { source =>
        SourceCallInput(
          sourceId = call.sourceId,
          sourceType = source.sourceType,
          calledAllele = call.calledAllele,
          callState = call.callState,
          readDepth = call.readDepth,
          mappingQuality = call.mappingQuality,
          callableState = call.callableState,
          calledRepeatCount = call.calledRepeatCount
        )
      }
    }

    val result = variant.variantType match
      case YVariantType.STR =>
        YVariantConcordance.calculateStrConsensus(callInputs, isInTree)
      case _ =>
        YVariantConcordance.calculateConsensus(callInputs, variant.variantType, isInTree)

    val updated = variant.copy(
      consensusAllele = result.consensusAllele,
      consensusState = result.consensusState,
      status = result.status,
      sourceCount = result.sourceCount,
      concordantCount = result.concordantCount,
      discordantCount = result.discordantCount,
      confidenceScore = result.confidenceScore,
      maxReadDepth = calls.flatMap(_.readDepth).maxOption,
      maxQualityScore = calls.flatMap(_.qualityScore).maxOption,
      lastUpdatedAt = LocalDateTime.now()
    )

    variantRepo.update(updated)

  // ============================================
  // Manual Curation
  // ============================================

  /**
   * Override a variant's consensus with manual curation.
   *
   * Creates an audit trail entry.
   */
  def overrideVariant(
    variantId: UUID,
    newConsensusAllele: String,
    newConsensusState: YConsensusState,
    newStatus: YVariantStatus,
    reason: String,
    userId: Option[String] = None,
    supportingEvidence: Option[String] = None
  ): Either[String, YProfileVariantEntity] =
    transactor.readWrite {
      val variant = variantRepo.findById(variantId).getOrElse(
        throw new IllegalArgumentException(s"Variant not found: $variantId")
      )

      // Create audit entry
      auditRepo.insert(YVariantAuditEntity.create(
        variantId = variantId,
        action = YAuditAction.OVERRIDE,
        reason = reason,
        previousConsensusAllele = variant.consensusAllele,
        previousConsensusState = Some(variant.consensusState),
        previousStatus = Some(variant.status),
        previousConfidence = Some(variant.confidenceScore),
        newConsensusAllele = Some(newConsensusAllele),
        newConsensusState = Some(newConsensusState),
        newStatus = Some(newStatus),
        newConfidence = Some(1.0), // Manual override = full confidence
        userId = userId,
        supportingEvidence = supportingEvidence
      ))

      // Update variant
      val updated = variant.copy(
        consensusAllele = Some(newConsensusAllele),
        consensusState = newConsensusState,
        status = newStatus,
        confidenceScore = 1.0,
        lastUpdatedAt = LocalDateTime.now()
      )

      variantRepo.update(updated)
    }

  /**
   * Revert a variant to its calculated consensus.
   */
  def revertOverride(
    variantId: UUID,
    reason: String,
    userId: Option[String] = None,
    isInTree: Boolean = true
  ): Either[String, YProfileVariantEntity] =
    transactor.readWrite {
      val variant = variantRepo.findById(variantId).getOrElse(
        throw new IllegalArgumentException(s"Variant not found: $variantId")
      )

      // Create audit entry
      auditRepo.insert(YVariantAuditEntity.create(
        variantId = variantId,
        action = YAuditAction.REVERT,
        reason = reason,
        previousConsensusAllele = variant.consensusAllele,
        previousConsensusState = Some(variant.consensusState),
        previousStatus = Some(variant.status),
        previousConfidence = Some(variant.confidenceScore),
        userId = userId
      ))

      // Re-reconcile the variant
      reconcileVariantInternal(variant, isInTree)
    }

  /**
   * Get audit history for a variant.
   */
  def getAuditHistory(variantId: UUID): Either[String, List[YVariantAuditEntity]] =
    transactor.readOnly {
      auditRepo.findByVariant(variantId)
    }

  // ============================================
  // Profile Statistics
  // ============================================

  /**
   * Update profile statistics from variants.
   */
  def updateProfileStatistics(profileId: UUID): Either[String, YChromosomeProfileEntity] =
    transactor.readWrite {
      val profile = profileRepo.findById(profileId).getOrElse(
        throw new IllegalArgumentException(s"Profile not found: $profileId")
      )

      val statusCounts = variantRepo.countByStatus(profileId)
      val variants = variantRepo.findByProfile(profileId)

      val totalVariants = variants.count(_.variantType != YVariantType.STR)
      val strMarkerCount = variants.count(_.variantType == YVariantType.STR)
      val strConfirmedCount = variants.count(v =>
        v.variantType == YVariantType.STR && v.status == YVariantStatus.CONFIRMED
      )

      val confirmedCount = statusCounts.getOrElse(YVariantStatus.CONFIRMED, 0L).toInt
      val novelCount = statusCounts.getOrElse(YVariantStatus.NOVEL, 0L).toInt
      val conflictCount = statusCounts.getOrElse(YVariantStatus.CONFLICT, 0L).toInt
      val noCoverageCount = statusCounts.getOrElse(YVariantStatus.NO_COVERAGE, 0L).toInt

      val overallConfidence = YVariantConcordance.calculateProfileConfidence(
        confirmedCount, novelCount, conflictCount, totalVariants + strMarkerCount
      )

      val updated = profile.copy(
        totalVariants = totalVariants,
        confirmedCount = confirmedCount,
        novelCount = novelCount,
        conflictCount = conflictCount,
        noCoverageCount = noCoverageCount,
        strMarkerCount = strMarkerCount,
        strConfirmedCount = strConfirmedCount,
        overallConfidence = Some(overallConfidence)
      )

      profileRepo.update(updated)
    }

  /**
   * Update haplogroup assignment on a profile.
   */
  def updateHaplogroup(
    profileId: UUID,
    haplogroup: String,
    confidence: Double,
    treeProvider: String,
    treeVersion: String
  ): Either[String, YChromosomeProfileEntity] =
    transactor.readWrite {
      val profile = profileRepo.findById(profileId).getOrElse(
        throw new IllegalArgumentException(s"Profile not found: $profileId")
      )

      val updated = profile.copy(
        consensusHaplogroup = Some(haplogroup),
        haplogroupConfidence = Some(confidence),
        haplogroupTreeProvider = Some(treeProvider),
        haplogroupTreeVersion = Some(treeVersion)
      )

      profileRepo.update(updated)
    }

  // ============================================
  // Internal Helpers
  // ============================================

  private def updateProfileSourceCount(profileId: UUID)(using java.sql.Connection): Unit =
    val sourceCount = sourceRepo.countByProfile(profileId).toInt
    val sources = sourceRepo.findByProfile(profileId)
    val primarySourceType = sources.headOption.map(_.sourceType) // Highest tier first

    profileRepo.findById(profileId).foreach { profile =>
      profileRepo.update(profile.copy(
        sourceCount = sourceCount,
        primarySourceType = primarySourceType
      ))
    }

  // ============================================
  // Haplogroup Determination
  // ============================================

  /**
   * Determine haplogroup from unified profile variants using HaplogroupScorer.
   *
   * Collects DERIVED SNP calls from the profile and scores them against
   * the haplogroup tree to determine the terminal haplogroup.
   *
   * @param profileId The profile to score
   * @param treeProviderType Tree provider (FTDNA or DecodingUs)
   * @param referenceBuild Reference build for tree coordinates (e.g., "hg38")
   * @return Either an error or the haplogroup scoring results
   */
  def determineHaplogroup(
    profileId: UUID,
    treeProviderType: TreeProviderType,
    referenceBuild: String
  ): Either[String, HaplogroupDeterminationResult] =
    // Load the haplogroup tree (outside transaction to avoid long DB locks)
    val treeProvider: TreeProvider = treeProviderType match
      case TreeProviderType.FTDNA => new FtdnaTreeProvider(TreeType.YDNA)
      case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(TreeType.YDNA)

    treeProvider.loadTree(referenceBuild).flatMap { tree =>
      // Get SNP calls from profile
      transactor.readOnly {
        val variants = variantRepo.findByProfile(profileId)

        // Collect DERIVED SNPs with CONFIRMED or NOVEL status
        // These are our high-confidence calls for haplogroup determination
        val snpCalls: Map[Long, String] = variants
          .filter { v =>
            v.variantType == YVariantType.SNP &&
            v.consensusState == YConsensusState.DERIVED &&
            (v.status == YVariantStatus.CONFIRMED || v.status == YVariantStatus.NOVEL)
          }
          .flatMap { v =>
            v.consensusAllele.map(allele => v.position -> allele)
          }
          .toMap

        // Also include ANCESTRAL calls for proper scoring
        val ancestralCalls: Map[Long, String] = variants
          .filter { v =>
            v.variantType == YVariantType.SNP &&
            v.consensusState == YConsensusState.ANCESTRAL &&
            (v.status == YVariantStatus.CONFIRMED || v.status == YVariantStatus.NOVEL)
          }
          .flatMap { v =>
            v.consensusAllele.map(allele => v.position -> allele)
          }
          .toMap

        // Combine all calls (derived takes precedence if same position)
        ancestralCalls ++ snpCalls
      }.flatMap { allCalls =>
        if allCalls.isEmpty then
          Left("No SNP calls available for haplogroup determination")
        else
          // Score haplogroups
          val scorer = new HaplogroupScorer()
          val results = scorer.score(tree, allCalls)

          if results.isEmpty then
            Left("No haplogroup results from scoring")
          else
            // Get terminal haplogroup (highest score, shallowest depth)
            val terminal = results.head
            val confidence = calculateHaplogroupConfidence(results)

            val result = HaplogroupDeterminationResult(
              terminalHaplogroup = terminal.name,
              score = terminal.score,
              confidence = confidence,
              matchingSnps = terminal.matchingSnps,
              ancestralSnps = terminal.ancestralMatches,
              noCalls = terminal.noCalls,
              totalTreeSnps = terminal.cumulativeSnps,
              depth = terminal.depth,
              topResults = results.take(10),
              treeProvider = treeProviderType.toString,
              treeVersion = treeProvider.cachePrefix,
              snpCallCount = allCalls.size
            )

            // Update the profile with the determined haplogroup
            updateHaplogroup(
              profileId,
              result.terminalHaplogroup,
              result.confidence,
              result.treeProvider,
              result.treeVersion
            ).map(_ => result)
      }
    }

  /**
   * Calculate haplogroup confidence based on scoring results.
   *
   * Higher confidence when:
   * - Terminal haplogroup score is well above next-best alternatives
   * - High proportion of tree SNPs have callable data
   * - Low number of ancestral mismatches on path to terminal
   */
  private def calculateHaplogroupConfidence(results: List[HaplogroupResult]): Double =
    if results.isEmpty then return 0.0

    val terminal = results.head
    if terminal.cumulativeSnps == 0 then return 0.0

    // Coverage factor: proportion of tree positions with calls
    val callableProportion = 1.0 - (terminal.noCalls.toDouble / terminal.cumulativeSnps.toDouble)

    // Score factor: how much better is the terminal vs alternatives
    val scoreFactor = results.lift(1) match
      case Some(nextBest) if nextBest.score > 0 =>
        math.min(terminal.score / nextBest.score, 2.0) / 2.0 // Normalize to 0-1
      case Some(_) =>
        1.0 // Next best has zero or negative score
      case None =>
        0.5 // Only one result - uncertain

    // Match factor: proportion of derived matches vs ancestral on path
    val totalPathCalls = terminal.matchingSnps + terminal.ancestralMatches
    val matchFactor = if totalPathCalls > 0 then
      terminal.matchingSnps.toDouble / totalPathCalls.toDouble
    else
      0.0

    // Combined confidence
    val confidence = callableProportion * 0.3 + scoreFactor * 0.4 + matchFactor * 0.3
    math.min(1.0, math.max(0.0, confidence))

/**
 * Result of haplogroup determination from unified profile.
 */
case class HaplogroupDeterminationResult(
  terminalHaplogroup: String,
  score: Double,
  confidence: Double,
  matchingSnps: Int,
  ancestralSnps: Int,
  noCalls: Int,
  totalTreeSnps: Int,
  depth: Int,
  topResults: List[HaplogroupResult],
  treeProvider: String,
  treeVersion: String,
  snpCallCount: Int
)
