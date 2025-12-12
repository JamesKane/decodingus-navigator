package com.decodingus.yprofile.service

import com.decodingus.analysis.{CallableLociQueryService, CallableState}
import com.decodingus.db.Transactor
import com.decodingus.haplogroup.model.HaplogroupResult
import com.decodingus.haplogroup.scoring.HaplogroupScorer
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.repository.{
  YChromosomeProfileRepository, YProfileSourceRepository, YProfileRegionRepository,
  YProfileVariantRepository, YVariantSourceCallRepository, YVariantAuditRepository,
  YSourceCallAlignmentRepository, YSnpPanelEntity, YSnpCall, YPrivateVariant
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
  auditRepo: YVariantAuditRepository,
  alignmentRepo: YSourceCallAlignmentRepository = YSourceCallAlignmentRepository()
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
   * Add or update a variant with a source call and alignment record.
   *
   * If the variant already exists (by position+alleles), adds a source call.
   * Otherwise creates a new variant. Also creates an alignment record
   * for the source call's coordinates in the specified reference build.
   *
   * @return Tuple of (variant, sourceCall, alignment)
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
    haplogroupBranchDepth: Option[Int] = None,
    referenceBuild: Option[String] = None,
    contig: String = "chrY"
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

      // Get source info for weight calculation and reference build
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
      val existingCall = sourceCallRepo.findByVariantAndSource(variant.id, sourceId)
      existingCall.foreach { existing =>
        // Delete existing alignments for this source call
        alignmentRepo.deleteBySourceCall(existing.id)
        sourceCallRepo.delete(existing.id)
      }

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

      // Create alignment record for coordinates
      // Use provided referenceBuild, or source's referenceBuild, or default to GRCh38
      val refBuild = referenceBuild
        .orElse(source.referenceBuild)
        .getOrElse("GRCh38")

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = refBuild,
        contig = contig,
        position = position,
        refAllele = refAllele,
        altAllele = altAllele,
        calledAllele = calledAllele,
        readDepth = readDepth,
        mappingQuality = mappingQuality,
        variantAlleleFrequency = variantAlleleFrequency
      ))

      (variant, sourceCall)
    }

  /**
   * Add an additional alignment to an existing source call.
   *
   * This is used when the same variant data is aligned to multiple references
   * (e.g., GRCh37, GRCh38, hs1). Each alignment represents coordinates in a
   * different reference, but they all count as ONE piece of evidence.
   */
  def addAlignmentToSourceCall(
    sourceCallId: UUID,
    referenceBuild: String,
    position: Long,
    refAllele: String,
    altAllele: String,
    calledAllele: String,
    contig: String = "chrY",
    readDepth: Option[Int] = None,
    mappingQuality: Option[Double] = None,
    variantAlleleFrequency: Option[Double] = None
  ): Either[String, YSourceCallAlignmentEntity] =
    transactor.readWrite {
      // Upsert alignment (replace if same source_call + reference_build exists)
      alignmentRepo.upsert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCallId,
        referenceBuild = referenceBuild,
        contig = contig,
        position = position,
        refAllele = refAllele,
        altAllele = altAllele,
        calledAllele = calledAllele,
        readDepth = readDepth,
        mappingQuality = mappingQuality,
        variantAlleleFrequency = variantAlleleFrequency
      ))
    }

  /**
   * Get alignments for a source call.
   */
  def getAlignments(sourceCallId: UUID): Either[String, List[YSourceCallAlignmentEntity]] =
    transactor.readOnly {
      alignmentRepo.findBySourceCall(sourceCallId)
    }

  /**
   * Get alignment for a specific reference build.
   */
  def getAlignmentForBuild(
    sourceCallId: UUID,
    referenceBuild: String
  ): Either[String, Option[YSourceCallAlignmentEntity]] =
    transactor.readOnly {
      alignmentRepo.findBySourceCallAndBuild(sourceCallId, referenceBuild)
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
  // Callable Loci Integration
  // ============================================

  /**
   * Import callable loci regions from a CallableLociQueryService for a source.
   *
   * This reads the BED file regions and stores them in y_profile_region,
   * enabling callable state queries for variant positions.
   *
   * @param profileId The Y profile ID
   * @param sourceId The source these regions belong to
   * @param callableLociService The query service with loaded BED data
   * @param contig Contig to import (default: chrY)
   * @return Count of regions imported
   */
  def importCallableRegions(
    profileId: UUID,
    sourceId: UUID,
    callableLociService: CallableLociQueryService,
    contig: String = "chrY"
  ): Either[String, Int] =
    // Check if data exists for this contig
    if !callableLociService.hasDataForContig(contig) then
      return Right(0)

    // Get the intervals by querying a range - we'll use the BED loading
    // The CallableLociQueryService loads intervals internally, but we need to extract them
    // For now, we'll store summary regions based on what we can query

    transactor.readWrite {
      // First, delete existing regions for this source
      regionRepo.deleteBySource(sourceId)

      // Get callable bases count for profile statistics
      val callableBases = callableLociService.getCallableBasesForContig(contig)

      // Store a summary region for the source
      // In a production system, we'd want to expose the intervals from CallableLociQueryService
      // For now, we store the summary with the callable bases
      val summaryRegion = YProfileRegionEntity.create(
        yProfileId = profileId,
        sourceId = sourceId,
        contig = contig,
        startPosition = 0,
        endPosition = callableBases, // Using this field to store callable base count
        callableState = YCallableState.SUMMARY,
        meanCoverage = None,
        meanMappingQuality = None
      )

      regionRepo.insert(summaryRegion)

      // Update profile's callable region percentage (inline since we're already in a transaction)
      updateCallableRegionPctInternal(profileId, contig)

      1 // Return count of regions imported
    }

  /**
   * Import detailed callable loci intervals for a source.
   *
   * This version takes pre-parsed intervals for batch import.
   *
   * @param profileId The Y profile ID
   * @param sourceId The source these regions belong to
   * @param intervals List of (start, end, state) tuples
   * @param contig Contig (default: chrY)
   * @return Count of regions imported
   */
  def importCallableIntervals(
    profileId: UUID,
    sourceId: UUID,
    intervals: List[(Long, Long, YCallableState)],
    contig: String = "chrY"
  ): Either[String, Int] =
    transactor.readWrite {
      // Delete existing regions for this source
      regionRepo.deleteBySource(sourceId)

      // Insert all intervals
      val regions = intervals.map { case (start, end, state) =>
        YProfileRegionEntity.create(
          yProfileId = profileId,
          sourceId = sourceId,
          contig = contig,
          startPosition = start,
          endPosition = end,
          callableState = state
        )
      }

      regions.foreach(regionRepo.insert)

      // Update profile statistics (inline since we're already in a transaction)
      updateCallableRegionPctInternal(profileId, contig)

      regions.size
    }

  /**
   * Query the callable state at a specific position for a profile.
   *
   * Checks stored regions first. If a live CallableLociQueryService is provided,
   * falls back to querying it directly.
   *
   * @param profileId The profile ID
   * @param position Genomic position (1-based)
   * @param contig Contig (default: chrY)
   * @param liveService Optional live query service for fallback
   * @return Callable state at position
   */
  def queryCallableState(
    profileId: UUID,
    position: Long,
    contig: String = "chrY",
    liveService: Option[CallableLociQueryService] = None
  ): Either[String, YCallableState] =
    transactor.readOnly {
      // Query stored regions first
      val overlapping = regionRepo.findOverlapping(profileId, position)
        .filter(_.contig == contig)
        .filterNot(_.callableState == YCallableState.SUMMARY) // Skip summary regions

      overlapping.headOption match
        case Some(region) => region.callableState
        case None =>
          // Fall back to live service if available
          liveService match
            case Some(service) =>
              convertCallableState(service.queryPosition(contig, position))
            case None =>
              // No stored region and no live service - unknown
              YCallableState.NO_COVERAGE
    }

  /**
   * Query callable states for multiple positions at once.
   *
   * @param profileId The profile ID
   * @param positions List of positions to query
   * @param contig Contig (default: chrY)
   * @param liveService Optional live query service for fallback
   * @return Map of position to callable state
   */
  def queryCallableStates(
    profileId: UUID,
    positions: List[Long],
    contig: String = "chrY",
    liveService: Option[CallableLociQueryService] = None
  ): Either[String, Map[Long, YCallableState]] =
    transactor.readOnly {
      // Get all relevant regions for the profile
      val allRegions = regionRepo.findByProfile(profileId)
        .filter(_.contig == contig)
        .filterNot(_.callableState == YCallableState.SUMMARY)

      positions.map { pos =>
        val state = allRegions.find(r => pos >= r.startPosition && pos <= r.endPosition) match
          case Some(region) => region.callableState
          case None =>
            liveService match
              case Some(service) =>
                convertCallableState(service.queryPosition(contig, pos))
              case None =>
                YCallableState.NO_COVERAGE
        pos -> state
      }.toMap
    }

  /**
   * Update callable region percentage on profile.
   *
   * Calculates percentage based on stored regions vs total chrY length.
   *
   * @param profileId Profile ID
   * @param contig Contig (default: chrY)
   * @param totalChrYLength Total chrY length in reference (default: GRCh38 ~57Mb)
   */
  def updateCallableRegionPct(
    profileId: UUID,
    contig: String = "chrY",
    totalChrYLength: Long = 57227415L // GRCh38 chrY length
  ): Either[String, YChromosomeProfileEntity] =
    transactor.readWrite {
      updateCallableRegionPctInternal(profileId, contig, totalChrYLength)
    }

  /**
   * Internal version that works within an existing transaction.
   */
  private def updateCallableRegionPctInternal(
    profileId: UUID,
    contig: String = "chrY",
    totalChrYLength: Long = 57227415L
  )(using java.sql.Connection): YChromosomeProfileEntity =
    val profile = profileRepo.findById(profileId).getOrElse(
      throw new IllegalArgumentException(s"Profile not found: $profileId")
    )

    // Get all CALLABLE regions for this profile
    val callableRegions = regionRepo.findByState(profileId, YCallableState.CALLABLE)
      .filter(_.contig == contig)

    // Calculate total callable bases
    val callableBases = if callableRegions.isEmpty then
      // Check for summary region
      regionRepo.findByState(profileId, YCallableState.SUMMARY)
        .filter(_.contig == contig)
        .headOption
        .map(_.endPosition) // We stored callable bases in endPosition
        .getOrElse(0L)
    else
      callableRegions.map(r => r.endPosition - r.startPosition + 1).sum

    val pct = if totalChrYLength > 0 then
      callableBases.toDouble / totalChrYLength.toDouble
    else
      0.0

    val updated = profile.copy(
      callableRegionPct = Some(pct)
    )

    profileRepo.update(updated)

  /**
   * Get callable summary for a source.
   *
   * Returns statistics about callable regions from a specific source.
   */
  def getCallableSummary(sourceId: UUID): Either[String, CallableSummary] =
    transactor.readOnly {
      val regions = regionRepo.findBySource(sourceId)

      val callableCount = regions.count(_.callableState == YCallableState.CALLABLE)
      val lowCoverageCount = regions.count(_.callableState == YCallableState.LOW_COVERAGE)
      val noCoverageCount = regions.count(_.callableState == YCallableState.NO_COVERAGE)
      val poorMappingCount = regions.count(_.callableState == YCallableState.POOR_MAPPING_QUALITY)

      val callableBases = regions
        .filter(_.callableState == YCallableState.CALLABLE)
        .map(r => r.endPosition - r.startPosition + 1)
        .sum

      val totalBases = regions
        .filterNot(_.callableState == YCallableState.SUMMARY)
        .map(r => r.endPosition - r.startPosition + 1)
        .sum

      CallableSummary(
        sourceId = sourceId,
        regionCount = regions.size,
        callableRegionCount = callableCount,
        lowCoverageRegionCount = lowCoverageCount,
        noCoverageRegionCount = noCoverageCount,
        poorMappingRegionCount = poorMappingCount,
        callableBases = callableBases,
        totalBases = totalBases,
        callablePct = if totalBases > 0 then callableBases.toDouble / totalBases.toDouble else 0.0
      )
    }

  /**
   * Convert CallableState (from CallableLociQueryService) to YCallableState.
   */
  private def convertCallableState(state: CallableState): YCallableState =
    state match
      case CallableState.Callable => YCallableState.CALLABLE
      case CallableState.NoCoverage => YCallableState.NO_COVERAGE
      case CallableState.LowCoverage => YCallableState.LOW_COVERAGE
      case CallableState.ExcessiveCoverage => YCallableState.EXCESSIVE_COVERAGE
      case CallableState.PoorMappingQuality => YCallableState.POOR_MAPPING_QUALITY
      case CallableState.RefN => YCallableState.REF_N
      case CallableState.Unknown => YCallableState.NO_COVERAGE

  // ============================================
  // Source Import
  // ============================================

  /**
   * Import variants from a Y-SNP panel into the unified profile.
   *
   * Creates a source entity and imports all SNP calls and private variants
   * as variant calls with appropriate concordance weights.
   *
   * @param profileId Profile to import into
   * @param panel The YSnpPanelEntity to import from
   * @param sourceType Source type (default: TARGETED_NGS for Big Y, CHIP for panels)
   * @return Import result with counts
   */
  def importFromSnpPanel(
    profileId: UUID,
    panel: YSnpPanelEntity,
    sourceType: YProfileSourceType = YProfileSourceType.TARGETED_NGS
  ): Either[String, SourceImportResult] =
    // First create the source outside the main transaction
    val sourceResult = addSource(
      profileId = profileId,
      sourceType = sourceType,
      vendor = panel.provider,
      testName = panel.panelName,
      testDate = panel.testDate,
      referenceBuild = Some("GRCh38") // SNP panels use GRCh38
    )

    sourceResult.flatMap { source =>
      transactor.readWrite {
        var snpCount = 0
        var privateCount = 0
        var errorCount = 0

        // Import SNP calls
        panel.snpCalls.foreach { call =>
          try {
            val callState = if call.derived then YConsensusState.DERIVED else YConsensusState.ANCESTRAL
            val variantType = call.variantType match
              case Some(com.decodingus.repository.YVariantType.INDEL) =>
                com.decodingus.yprofile.model.YVariantType.INDEL
              case _ =>
                com.decodingus.yprofile.model.YVariantType.SNP

            // For SNP panels, we typically don't have ref/alt alleles, just the called allele
            // We'll use the allele as both ref (ancestral) and alt (derived) based on state
            val (refAllele, altAllele) = if call.derived then
              ("N", call.allele) // Derived: called allele is the alt
            else
              (call.allele, "N") // Ancestral: called allele is the ref

            addVariantCallInternal(
              profileId = profileId,
              sourceId = source.id,
              position = call.startPosition,
              refAllele = refAllele,
              altAllele = altAllele,
              calledAllele = call.allele,
              callState = callState,
              variantType = variantType,
              variantName = Some(call.name),
              qualityScore = call.quality,
              referenceBuild = "GRCh38"
            )
            snpCount += 1
          } catch {
            case e: Exception =>
              errorCount += 1
          }
        }

        // Import private variants
        panel.privateVariants.foreach { pv =>
          try {
            addVariantCallInternal(
              profileId = profileId,
              sourceId = source.id,
              position = pv.position,
              refAllele = pv.refAllele,
              altAllele = pv.altAllele,
              calledAllele = pv.altAllele, // Private variants are derived
              callState = YConsensusState.DERIVED,
              variantType = com.decodingus.yprofile.model.YVariantType.SNP,
              variantName = pv.snpName,
              qualityScore = pv.quality,
              readDepth = pv.readDepth,
              referenceBuild = "GRCh38"
            )
            privateCount += 1
          } catch {
            case e: Exception =>
              errorCount += 1
          }
        }

        SourceImportResult(
          sourceId = source.id,
          sourceType = sourceType,
          snpCallsImported = snpCount,
          privateVariantsImported = privateCount,
          errorsEncountered = errorCount
        )
      }
    }

  /**
   * Import variants from a list of SNP calls (generic import).
   *
   * This is useful for importing from parsed VCF, haplogroup analysis results,
   * or other sources that provide SNP call data.
   *
   * @param profileId Profile to import into
   * @param sourceId Source that these calls belong to
   * @param calls List of (position, refAllele, altAllele, calledAllele, derived, variantName)
   * @param referenceBuild Reference build for coordinates
   * @return Import result with counts
   */
  def importVariantCalls(
    profileId: UUID,
    sourceId: UUID,
    calls: List[VariantCallInput],
    referenceBuild: String = "GRCh38"
  ): Either[String, SourceImportResult] =
    transactor.readWrite {
      var successCount = 0
      var errorCount = 0

      calls.foreach { call =>
        try {
          val callState = if call.derived then YConsensusState.DERIVED else YConsensusState.ANCESTRAL

          addVariantCallInternal(
            profileId = profileId,
            sourceId = sourceId,
            position = call.position,
            refAllele = call.refAllele,
            altAllele = call.altAllele,
            calledAllele = call.calledAllele,
            callState = callState,
            variantType = call.variantType.getOrElse(com.decodingus.yprofile.model.YVariantType.SNP),
            variantName = call.variantName,
            readDepth = call.readDepth,
            qualityScore = call.qualityScore,
            mappingQuality = call.mappingQuality,
            referenceBuild = referenceBuild
          )
          successCount += 1
        } catch {
          case e: Exception =>
            errorCount += 1
        }
      }

      SourceImportResult(
        sourceId = sourceId,
        sourceType = YProfileSourceType.WGS_SHORT_READ, // Will be overwritten
        snpCallsImported = successCount,
        privateVariantsImported = 0,
        errorsEncountered = errorCount
      )
    }

  /**
   * Internal variant call addition (within existing transaction).
   */
  private def addVariantCallInternal(
    profileId: UUID,
    sourceId: UUID,
    position: Long,
    refAllele: String,
    altAllele: String,
    calledAllele: String,
    callState: YConsensusState,
    variantType: com.decodingus.yprofile.model.YVariantType = com.decodingus.yprofile.model.YVariantType.SNP,
    variantName: Option[String] = None,
    readDepth: Option[Int] = None,
    qualityScore: Option[Double] = None,
    mappingQuality: Option[Double] = None,
    referenceBuild: String = "GRCh38"
  )(using java.sql.Connection): Unit =
    // Find or create variant
    val variant = variantRepo.findByPosition(profileId, position) match
      case Some(existing) =>
        val updated = existing.copy(
          variantName = variantName.orElse(existing.variantName)
        )
        if updated != existing then variantRepo.update(updated) else existing
      case None =>
        variantRepo.insert(YProfileVariantEntity.create(
          yProfileId = profileId,
          position = position,
          refAllele = refAllele,
          altAllele = altAllele,
          variantType = variantType,
          variantName = variantName
        ))

    // Get source for weight calculation
    val source = sourceRepo.findById(sourceId).getOrElse(
      throw new IllegalArgumentException(s"Source not found: $sourceId")
    )

    val weight = YVariantConcordance.calculateWeight(
      source.sourceType,
      variantType,
      readDepth,
      mappingQuality,
      None
    )

    // Delete existing call if present
    sourceCallRepo.findByVariantAndSource(variant.id, sourceId).foreach { existing =>
      alignmentRepo.deleteBySourceCall(existing.id)
      sourceCallRepo.delete(existing.id)
    }

    // Insert source call
    val sourceCall = sourceCallRepo.insert(YVariantSourceCallEntity.create(
      variantId = variant.id,
      sourceId = sourceId,
      calledAllele = calledAllele,
      callState = callState,
      readDepth = readDepth,
      qualityScore = qualityScore,
      mappingQuality = mappingQuality,
      concordanceWeight = weight
    ))

    // Insert alignment
    alignmentRepo.insert(YSourceCallAlignmentEntity.create(
      sourceCallId = sourceCall.id,
      referenceBuild = referenceBuild,
      position = position,
      refAllele = refAllele,
      altAllele = altAllele,
      calledAllele = calledAllele,
      readDepth = readDepth,
      mappingQuality = mappingQuality
    ))

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

/**
 * Summary of callable loci statistics for a source.
 */
case class CallableSummary(
  sourceId: UUID,
  regionCount: Int,
  callableRegionCount: Int,
  lowCoverageRegionCount: Int,
  noCoverageRegionCount: Int,
  poorMappingRegionCount: Int,
  callableBases: Long,
  totalBases: Long,
  callablePct: Double
)

/**
 * Result of importing variants from a source.
 */
case class SourceImportResult(
  sourceId: UUID,
  sourceType: YProfileSourceType,
  snpCallsImported: Int,
  privateVariantsImported: Int,
  errorsEncountered: Int
) {
  def totalImported: Int = snpCallsImported + privateVariantsImported
}

/**
 * Input for importing a variant call.
 */
case class VariantCallInput(
  position: Long,
  refAllele: String,
  altAllele: String,
  calledAllele: String,
  derived: Boolean,
  variantName: Option[String] = None,
  variantType: Option[YVariantType] = None,
  readDepth: Option[Int] = None,
  qualityScore: Option[Double] = None,
  mappingQuality: Option[Double] = None
)
