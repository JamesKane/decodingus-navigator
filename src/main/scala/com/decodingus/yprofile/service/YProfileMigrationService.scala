package com.decodingus.yprofile.service

import com.decodingus.db.Transactor
import com.decodingus.repository.{
  HaplogroupReconciliationEntity, HaplogroupReconciliationRepository,
  YChromosomeProfileRepository, YProfileSourceRepository, YProfileRegionRepository,
  YProfileVariantRepository, YVariantSourceCallRepository, YVariantAuditRepository
}
import com.decodingus.workspace.model.{
  DnaType, HaplogroupTechnology, RunHaplogroupCall, SnpConflict, SnpCallFromRun
}
import com.decodingus.yprofile.model.*
import java.time.{LocalDateTime, ZoneId}
import java.util.UUID

/**
 * Service for migrating Y-DNA data from HaplogroupReconciliation to YChromosomeProfile.
 *
 * This is a one-way migration that:
 * 1. Reads existing Y-DNA HaplogroupReconciliation records
 * 2. Creates equivalent YChromosomeProfile records
 * 3. Preserves source information and haplogroup assignments
 *
 * After migration, HaplogroupReconciliation is retained only for MT-DNA.
 */
class YProfileMigrationService(
  transactor: Transactor,
  reconciliationRepo: HaplogroupReconciliationRepository,
  profileRepo: YChromosomeProfileRepository,
  sourceRepo: YProfileSourceRepository,
  regionRepo: YProfileRegionRepository,
  variantRepo: YProfileVariantRepository,
  sourceCallRepo: YVariantSourceCallRepository,
  auditRepo: YVariantAuditRepository
):

  /**
   * Result of a single record migration.
   */
  case class MigrationResult(
    biosampleId: UUID,
    profileId: UUID,
    sourcesCreated: Int,
    variantsCreated: Int,
    success: Boolean,
    errorMessage: Option[String] = None
  )

  /**
   * Summary of migration operation.
   */
  case class MigrationSummary(
    totalRecords: Int,
    successfulMigrations: Int,
    failedMigrations: Int,
    results: List[MigrationResult]
  )

  /**
   * Migrate all Y-DNA HaplogroupReconciliation records to YChromosomeProfile.
   *
   * @param dryRun If true, report what would be migrated without making changes
   * @return Migration summary
   */
  def migrateAll(dryRun: Boolean = false): Either[String, MigrationSummary] =
    transactor.readOnly {
      reconciliationRepo.findByDnaType(DnaType.Y_DNA)
    }.flatMap { records =>
      if records.isEmpty then
        Right(MigrationSummary(0, 0, 0, List.empty))
      else if dryRun then
        // Dry run - just report what would be migrated
        val results = records.map { rec =>
          MigrationResult(
            biosampleId = rec.biosampleId,
            profileId = UUID.randomUUID(), // Would be generated
            sourcesCreated = rec.runCalls.size,
            variantsCreated = rec.snpConflicts.flatMap(_.calls).size,
            success = true
          )
        }
        Right(MigrationSummary(records.size, records.size, 0, results))
      else
        // Actual migration
        val results = records.map(migrateRecord)
        val successful = results.count(_.success)
        val failed = results.count(!_.success)
        Right(MigrationSummary(records.size, successful, failed, results))
    }

  /**
   * Migrate a single biosample's Y-DNA data.
   *
   * @param biosampleId The biosample to migrate
   * @return Migration result
   */
  def migrateBiosample(biosampleId: UUID): Either[String, MigrationResult] =
    transactor.readOnly {
      reconciliationRepo.findByBiosampleAndDnaType(biosampleId, DnaType.Y_DNA)
    }.flatMap {
      case Some(reconciliation) => Right(migrateRecord(reconciliation))
      case None => Left(s"No Y-DNA reconciliation found for biosample: $biosampleId")
    }

  /**
   * Check if a biosample has already been migrated.
   */
  def isMigrated(biosampleId: UUID): Either[String, Boolean] =
    transactor.readOnly {
      profileRepo.findByBiosample(biosampleId).isDefined
    }

  /**
   * Get migration status for all Y-DNA biosamples.
   */
  def getMigrationStatus: Either[String, List[(UUID, Boolean)]] =
    transactor.readOnly {
      val reconciliations = reconciliationRepo.findByDnaType(DnaType.Y_DNA)
      reconciliations.map { rec =>
        val migrated = profileRepo.findByBiosample(rec.biosampleId).isDefined
        (rec.biosampleId, migrated)
      }
    }

  /**
   * Migrate a single HaplogroupReconciliation record.
   */
  private def migrateRecord(reconciliation: HaplogroupReconciliationEntity): MigrationResult =
    try
      transactor.readWrite {
        // Check if already migrated
        profileRepo.findByBiosample(reconciliation.biosampleId) match
          case Some(existing) =>
            // Already migrated - update instead
            MigrationResult(
              biosampleId = reconciliation.biosampleId,
              profileId = existing.id,
              sourcesCreated = 0,
              variantsCreated = 0,
              success = true,
              errorMessage = Some("Already migrated - skipped")
            )
          case None =>
            // Create new profile
            val profile = profileRepo.insert(YChromosomeProfileEntity.create(
              biosampleId = reconciliation.biosampleId,
              consensusHaplogroup = Some(reconciliation.status.consensusHaplogroup).filter(_.nonEmpty),
              haplogroupConfidence = Some(reconciliation.status.confidence).filter(_ > 0)
            ))

            // Create sources from run calls
            val sourceEntities = reconciliation.runCalls.zipWithIndex.map { case (call, idx) =>
              createSourceFromRunCall(profile.id, call, idx)
            }

            // Create variants from SNP conflicts
            val variantEntities = reconciliation.snpConflicts.flatMap { conflict =>
              createVariantsFromConflict(profile.id, sourceEntities, conflict)
            }

            // Update profile statistics
            val updated = profile.copy(
              sourceCount = sourceEntities.size,
              totalVariants = variantEntities.size,
              primarySourceType = sourceEntities.headOption.map(_.sourceType)
            )
            profileRepo.update(updated)

            MigrationResult(
              biosampleId = reconciliation.biosampleId,
              profileId = profile.id,
              sourcesCreated = sourceEntities.size,
              variantsCreated = variantEntities.size,
              success = true
            )
      }.getOrElse(MigrationResult(
        biosampleId = reconciliation.biosampleId,
        profileId = UUID.randomUUID(),
        sourcesCreated = 0,
        variantsCreated = 0,
        success = false,
        errorMessage = Some("Transaction failed")
      ))
    catch
      case e: Exception =>
        MigrationResult(
          biosampleId = reconciliation.biosampleId,
          profileId = UUID.randomUUID(),
          sourcesCreated = 0,
          variantsCreated = 0,
          success = false,
          errorMessage = Some(e.getMessage)
        )

  /**
   * Create a YProfileSource from a RunHaplogroupCall.
   */
  private def createSourceFromRunCall(profileId: UUID, call: RunHaplogroupCall, index: Int)(using java.sql.Connection): YProfileSourceEntity =
    val sourceType = call.technology match
      case Some(HaplogroupTechnology.WGS) => YProfileSourceType.WGS_SHORT_READ
      case Some(HaplogroupTechnology.WES) => YProfileSourceType.TARGETED_NGS
      case Some(HaplogroupTechnology.BIG_Y) => YProfileSourceType.TARGETED_NGS
      case Some(HaplogroupTechnology.SNP_ARRAY) => YProfileSourceType.CHIP
      case Some(HaplogroupTechnology.AMPLICON) => YProfileSourceType.TARGETED_NGS
      case Some(HaplogroupTechnology.STR_PANEL) => YProfileSourceType.CAPILLARY_ELECTROPHORESIS
      case None => YProfileSourceType.MANUAL

    val methodTier = sourceType.snpTier

    sourceRepo.insert(YProfileSourceEntity.create(
      yProfileId = profileId,
      sourceType = sourceType,
      sourceRef = Some(call.sourceRef),
      vendor = call.treeProvider,
      methodTier = methodTier,
      meanReadDepth = call.meanCoverage
    ))

  /**
   * Create YProfileVariant and source calls from an SNP conflict.
   */
  private def createVariantsFromConflict(
    profileId: UUID,
    sources: List[YProfileSourceEntity],
    conflict: SnpConflict
  )(using java.sql.Connection): List[YProfileVariantEntity] =
    // Create the variant
    val variant = variantRepo.insert(YProfileVariantEntity.create(
      yProfileId = profileId,
      position = conflict.position,
      refAllele = conflict.calls.headOption.map(_.allele).getOrElse("N"), // Best guess
      altAllele = conflict.calls.lastOption.map(_.allele).getOrElse("N"),
      variantType = YVariantType.SNP,
      variantName = conflict.snpName,
      status = conflict.resolution match
        case Some(_) => YVariantStatus.CONFIRMED
        case None => YVariantStatus.CONFLICT
    ))

    // Create source calls for each run in the conflict
    conflict.calls.foreach { snpCall =>
      // Find the source that matches this run reference
      val matchingSource = sources.find(_.sourceRef.contains(snpCall.runRef))
      matchingSource.foreach { source =>
        val callState = if conflict.calls.exists(_.allele != snpCall.allele) then
          YConsensusState.DERIVED // Has variance = derived from consensus
        else
          YConsensusState.ANCESTRAL

        sourceCallRepo.insert(YVariantSourceCallEntity.create(
          variantId = variant.id,
          sourceId = source.id,
          calledAllele = snpCall.allele,
          callState = callState,
          readDepth = snpCall.depth,
          qualityScore = snpCall.quality,
          variantAlleleFrequency = snpCall.variantAlleleFrequency,
          concordanceWeight = 1.0 // Will be recalculated during reconciliation
        ))
      }
    }

    List(variant)
