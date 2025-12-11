package com.decodingus.service

import com.decodingus.db.Transactor
import com.decodingus.repository.{
  AnalysisArtifactRepository, AnalysisArtifactEntity, ArtifactType, ArtifactStatus,
  SourceFileRepository, SourceFileEntity, SourceFileFormat
}
import io.circe.Json
import java.nio.file.{Files, Path, Paths}
import java.util.UUID

/**
 * Service for managing analysis cache and source file tracking.
 *
 * Provides high-level operations for:
 * - Registering and tracking source files (BAM/CRAM)
 * - Managing cached analysis artifacts
 * - Cache invalidation when dependencies change
 * - Cache statistics and cleanup
 */
trait CacheService:

  // ============================================
  // Source File Management
  // ============================================

  /**
   * Register a source file for tracking.
   * If the file already exists (by checksum), updates the path.
   */
  def registerSourceFile(
    filePath: String,
    fileChecksum: String,
    fileSize: Option[Long] = None,
    fileFormat: Option[SourceFileFormat] = None
  ): Either[String, SourceFileEntity]

  /**
   * Link a source file to an alignment.
   */
  def linkSourceFileToAlignment(sourceFileId: UUID, alignmentId: UUID): Either[String, Boolean]

  /**
   * Get source file by checksum.
   */
  def getSourceFileByChecksum(checksum: String): Either[String, Option[SourceFileEntity]]

  /**
   * Mark a source file as analyzed.
   */
  def markSourceFileAnalyzed(id: UUID): Either[String, Boolean]

  /**
   * Verify source file accessibility and update status.
   */
  def verifySourceFile(id: UUID): Either[String, Boolean]

  /**
   * Get all inaccessible source files.
   */
  def getInaccessibleSourceFiles(): Either[String, List[SourceFileEntity]]

  // ============================================
  // Artifact Management
  // ============================================

  /**
   * Start tracking a new artifact (mark as in-progress).
   */
  def startArtifact(
    alignmentId: UUID,
    artifactType: ArtifactType,
    cachePath: String,
    generatorVersion: Option[String] = None,
    generationParams: Option[Json] = None,
    dependsOnSourceChecksum: Option[String] = None,
    dependsOnReferenceBuild: Option[String] = None
  ): Either[String, AnalysisArtifactEntity]

  /**
   * Complete an artifact (mark as available).
   */
  def completeArtifact(
    id: UUID,
    fileSize: Long,
    fileChecksum: String,
    fileFormat: Option[String] = None
  ): Either[String, Boolean]

  /**
   * Mark an artifact as failed.
   */
  def failArtifact(id: UUID, reason: String): Either[String, Boolean]

  /**
   * Get artifact for an alignment and type.
   */
  def getArtifact(alignmentId: UUID, artifactType: ArtifactType): Either[String, Option[AnalysisArtifactEntity]]

  /**
   * Get all artifacts for an alignment.
   */
  def getArtifactsForAlignment(alignmentId: UUID): Either[String, List[AnalysisArtifactEntity]]

  /**
   * Check if an artifact is available (not stale, not in-progress).
   */
  def isArtifactAvailable(alignmentId: UUID, artifactType: ArtifactType): Either[String, Boolean]

  /**
   * Delete an artifact and its cached file.
   */
  def deleteArtifact(id: UUID): Either[String, Boolean]

  // ============================================
  // Cache Invalidation
  // ============================================

  /**
   * Invalidate all artifacts for an alignment.
   */
  def invalidateArtifactsForAlignment(alignmentId: UUID, reason: String): Either[String, Int]

  /**
   * Invalidate artifacts when source file changes.
   */
  def invalidateBySourceChecksum(oldChecksum: String, reason: String): Either[String, Int]

  /**
   * Invalidate artifacts when reference build changes.
   */
  def invalidateByReferenceBuild(referenceBuild: String, reason: String): Either[String, Int]

  /**
   * Validate a single artifact against its dependencies.
   * Returns true if valid, marks stale and returns false if invalid.
   */
  def validateArtifact(id: UUID): Either[String, ArtifactValidationResult]

  /**
   * Validate all available artifacts.
   * Returns validation results for all checked artifacts.
   */
  def validateAllArtifacts(): Either[String, BatchValidationResult]

  /**
   * Verify all source files and invalidate artifacts for inaccessible sources.
   */
  def verifyAllSourceFiles(): Either[String, SourceFileVerificationResult]

  /**
   * Check if an artifact's cached file exists on disk.
   */
  def verifyArtifactFile(id: UUID): Either[String, Boolean]

  /**
   * Clean up artifacts whose files are missing from disk.
   */
  def cleanupMissingArtifacts(): Either[String, Int]

  /**
   * Get all stale artifacts that need regeneration.
   */
  def getStaleArtifacts(): Either[String, List[AnalysisArtifactEntity]]

  // ============================================
  // Cache Statistics
  // ============================================

  /**
   * Get cache statistics summary.
   */
  def getCacheStats(): Either[String, CacheStats]

/**
 * Cache statistics summary.
 */
case class CacheStats(
  totalArtifacts: Int,
  availableArtifacts: Int,
  staleArtifacts: Int,
  inProgressArtifacts: Int,
  errorArtifacts: Int,
  totalCacheSizeBytes: Long,
  trackedSourceFiles: Int,
  accessibleSourceFiles: Int,
  inaccessibleSourceFiles: Int,
  analyzedSourceFiles: Int
)

/**
 * Result of validating a single artifact.
 */
enum ArtifactValidationResult:
  case Valid                                    // Artifact is valid and available
  case MarkedStale(reason: String)              // Artifact was marked stale
  case AlreadyStale                             // Artifact was already stale
  case NotFound                                 // Artifact doesn't exist
  case FileNotFound                             // Cached file is missing

/**
 * Result of batch artifact validation.
 */
case class BatchValidationResult(
  checkedCount: Int,
  validCount: Int,
  markedStaleCount: Int,
  alreadyStaleCount: Int,
  missingFileCount: Int,
  staleReasons: Map[String, Int]                // Reason -> count
)

/**
 * Result of source file verification.
 */
case class SourceFileVerificationResult(
  checkedCount: Int,
  accessibleCount: Int,
  inaccessibleCount: Int,
  newlyInaccessible: Int,
  artifactsInvalidated: Int
)

/**
 * H2 database-backed implementation of CacheService.
 */
class H2CacheService(
  transactor: Transactor,
  artifactRepo: AnalysisArtifactRepository,
  sourceFileRepo: SourceFileRepository
) extends CacheService:

  private val CacheDir: Path = Paths.get(System.getProperty("user.home"), ".decodingus", "cache")

  // ============================================
  // Source File Management
  // ============================================

  override def registerSourceFile(
    filePath: String,
    fileChecksum: String,
    fileSize: Option[Long],
    fileFormat: Option[SourceFileFormat]
  ): Either[String, SourceFileEntity] =
    transactor.readWrite {
      sourceFileRepo.upsertByChecksum(filePath, fileChecksum, fileSize, fileFormat)
    }

  override def linkSourceFileToAlignment(sourceFileId: UUID, alignmentId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      sourceFileRepo.linkToAlignment(sourceFileId, alignmentId)
    }

  override def getSourceFileByChecksum(checksum: String): Either[String, Option[SourceFileEntity]] =
    transactor.readOnly {
      sourceFileRepo.findByChecksum(checksum)
    }

  override def markSourceFileAnalyzed(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      sourceFileRepo.markAnalyzed(id)
    }

  override def verifySourceFile(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      sourceFileRepo.findById(id) match
        case Some(entity) =>
          entity.filePath match
            case Some(path) =>
              val file = Paths.get(path)
              if Files.exists(file) && Files.isReadable(file) then
                sourceFileRepo.markAccessible(id)
                true
              else
                sourceFileRepo.markInaccessible(id)
                false
            case None =>
              sourceFileRepo.markInaccessible(id)
              false
        case None =>
          false
    }

  override def getInaccessibleSourceFiles(): Either[String, List[SourceFileEntity]] =
    transactor.readOnly {
      sourceFileRepo.findInaccessible()
    }

  // ============================================
  // Artifact Management
  // ============================================

  override def startArtifact(
    alignmentId: UUID,
    artifactType: ArtifactType,
    cachePath: String,
    generatorVersion: Option[String],
    generationParams: Option[Json],
    dependsOnSourceChecksum: Option[String],
    dependsOnReferenceBuild: Option[String]
  ): Either[String, AnalysisArtifactEntity] =
    transactor.readWrite {
      // Check if artifact already exists for this alignment/type
      artifactRepo.findByAlignmentAndType(alignmentId, artifactType) match
        case Some(existing) if existing.status == ArtifactStatus.InProgress =>
          // Already in progress
          existing
        case Some(existing) =>
          // Update existing to in-progress
          val updated = existing.copy(
            cachePath = cachePath,
            generatorVersion = generatorVersion,
            generationParams = generationParams,
            status = ArtifactStatus.InProgress,
            staleReason = None,
            dependsOnSourceChecksum = dependsOnSourceChecksum,
            dependsOnReferenceBuild = dependsOnReferenceBuild
          )
          artifactRepo.update(updated)
        case None =>
          // Create new
          val entity = AnalysisArtifactEntity.create(
            alignmentId = alignmentId,
            artifactType = artifactType,
            cachePath = cachePath,
            generatorVersion = generatorVersion,
            generationParams = generationParams,
            dependsOnSourceChecksum = dependsOnSourceChecksum,
            dependsOnReferenceBuild = dependsOnReferenceBuild
          )
          artifactRepo.insert(entity)
    }

  override def completeArtifact(
    id: UUID,
    fileSize: Long,
    fileChecksum: String,
    fileFormat: Option[String]
  ): Either[String, Boolean] =
    transactor.readWrite {
      artifactRepo.markAvailable(id, fileSize, fileChecksum, fileFormat)
    }

  override def failArtifact(id: UUID, reason: String): Either[String, Boolean] =
    transactor.readWrite {
      artifactRepo.markError(id, reason)
    }

  override def getArtifact(alignmentId: UUID, artifactType: ArtifactType): Either[String, Option[AnalysisArtifactEntity]] =
    transactor.readOnly {
      artifactRepo.findByAlignmentAndType(alignmentId, artifactType)
    }

  override def getArtifactsForAlignment(alignmentId: UUID): Either[String, List[AnalysisArtifactEntity]] =
    transactor.readOnly {
      artifactRepo.findByAlignment(alignmentId)
    }

  override def isArtifactAvailable(alignmentId: UUID, artifactType: ArtifactType): Either[String, Boolean] =
    transactor.readOnly {
      artifactRepo.findByAlignmentAndType(alignmentId, artifactType) match
        case Some(artifact) =>
          artifact.status == ArtifactStatus.Available &&
            artifact.cachePath.nonEmpty &&
            Files.exists(CacheDir.resolve(artifact.cachePath))
        case None =>
          false
    }

  override def deleteArtifact(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      artifactRepo.findById(id) match
        case Some(entity) =>
          // Try to delete the cached file
          val cachePath = CacheDir.resolve(entity.cachePath)
          if Files.exists(cachePath) then
            try Files.delete(cachePath)
            catch case _: Exception => () // Ignore file deletion errors

          artifactRepo.delete(id)
        case None =>
          false
    }

  // ============================================
  // Cache Invalidation
  // ============================================

  override def invalidateArtifactsForAlignment(alignmentId: UUID, reason: String): Either[String, Int] =
    transactor.readWrite {
      val artifacts = artifactRepo.findByAlignment(alignmentId)
      var count = 0
      for artifact <- artifacts if artifact.status == ArtifactStatus.Available do
        if artifactRepo.markStale(artifact.id, reason) then count += 1
      count
    }

  override def invalidateBySourceChecksum(oldChecksum: String, reason: String): Either[String, Int] =
    transactor.readWrite {
      artifactRepo.markStaleBySourceChecksum(oldChecksum, reason)
    }

  override def invalidateByReferenceBuild(referenceBuild: String, reason: String): Either[String, Int] =
    transactor.readWrite {
      val artifacts = artifactRepo.findAll().filter { artifact =>
        artifact.status == ArtifactStatus.Available &&
          artifact.dependsOnReferenceBuild.contains(referenceBuild)
      }
      var count = 0
      for artifact <- artifacts do
        if artifactRepo.markStale(artifact.id, reason) then count += 1
      count
    }

  override def validateArtifact(id: UUID): Either[String, ArtifactValidationResult] =
    transactor.readWrite {
      artifactRepo.findById(id) match
        case None =>
          ArtifactValidationResult.NotFound

        case Some(artifact) if artifact.status == ArtifactStatus.Stale =>
          ArtifactValidationResult.AlreadyStale

        case Some(artifact) if artifact.status != ArtifactStatus.Available =>
          ArtifactValidationResult.Valid // Not available yet, can't validate

        case Some(artifact) =>
          // Check if cached file exists
          val cachePath = CacheDir.resolve(artifact.cachePath)
          if !Files.exists(cachePath) then
            artifactRepo.markDeleted(artifact.id)
            ArtifactValidationResult.FileNotFound
          else
            // Check source checksum dependency
            artifact.dependsOnSourceChecksum match
              case Some(expectedChecksum) =>
                // Find source file with this checksum
                val sourceExists = sourceFileRepo.findByChecksum(expectedChecksum).exists(_.isAccessible)
                if !sourceExists then
                  val reason = s"Source file no longer accessible (checksum: ${expectedChecksum.take(12)}...)"
                  artifactRepo.markStale(artifact.id, reason)
                  ArtifactValidationResult.MarkedStale(reason)
                else
                  ArtifactValidationResult.Valid
              case None =>
                ArtifactValidationResult.Valid
    }

  override def validateAllArtifacts(): Either[String, BatchValidationResult] =
    transactor.readWrite {
      val artifacts = artifactRepo.findByStatus(ArtifactStatus.Available)
      var validCount = 0
      var markedStaleCount = 0
      var alreadyStaleCount = 0
      var missingFileCount = 0
      val staleReasons = scala.collection.mutable.Map[String, Int]()

      for artifact <- artifacts do
        val cachePath = CacheDir.resolve(artifact.cachePath)

        if !Files.exists(cachePath) then
          artifactRepo.markDeleted(artifact.id)
          missingFileCount += 1
        else
          // Check source dependency
          artifact.dependsOnSourceChecksum match
            case Some(expectedChecksum) =>
              val sourceAccessible = sourceFileRepo.findByChecksum(expectedChecksum).exists(_.isAccessible)
              if !sourceAccessible then
                val reason = "Source file no longer accessible"
                artifactRepo.markStale(artifact.id, reason)
                markedStaleCount += 1
                staleReasons.updateWith(reason)(_.map(_ + 1).orElse(Some(1)))
              else
                validCount += 1
            case None =>
              validCount += 1

      // Also count already stale
      alreadyStaleCount = artifactRepo.findStale().size

      BatchValidationResult(
        checkedCount = artifacts.size,
        validCount = validCount,
        markedStaleCount = markedStaleCount,
        alreadyStaleCount = alreadyStaleCount,
        missingFileCount = missingFileCount,
        staleReasons = staleReasons.toMap
      )
    }

  override def verifyAllSourceFiles(): Either[String, SourceFileVerificationResult] =
    transactor.readWrite {
      val sourceFiles = sourceFileRepo.findAll()
      var accessibleCount = 0
      var inaccessibleCount = 0
      var newlyInaccessible = 0
      var artifactsInvalidated = 0

      for sf <- sourceFiles do
        val isAccessible = sf.filePath.exists { path =>
          val file = Paths.get(path)
          Files.exists(file) && Files.isReadable(file)
        }

        if isAccessible then
          accessibleCount += 1
          if !sf.isAccessible then
            sourceFileRepo.markAccessible(sf.id)
        else
          inaccessibleCount += 1
          if sf.isAccessible then
            sourceFileRepo.markInaccessible(sf.id)
            newlyInaccessible += 1
            // Invalidate artifacts depending on this source
            artifactsInvalidated += artifactRepo.markStaleBySourceChecksum(
              sf.fileChecksum,
              s"Source file no longer accessible: ${sf.filePath.getOrElse("unknown")}"
            )

      SourceFileVerificationResult(
        checkedCount = sourceFiles.size,
        accessibleCount = accessibleCount,
        inaccessibleCount = inaccessibleCount,
        newlyInaccessible = newlyInaccessible,
        artifactsInvalidated = artifactsInvalidated
      )
    }

  override def verifyArtifactFile(id: UUID): Either[String, Boolean] =
    transactor.readOnly {
      artifactRepo.findById(id) match
        case Some(artifact) =>
          val cachePath = CacheDir.resolve(artifact.cachePath)
          Files.exists(cachePath) && Files.isReadable(cachePath)
        case None =>
          false
    }

  override def cleanupMissingArtifacts(): Either[String, Int] =
    transactor.readWrite {
      val available = artifactRepo.findByStatus(ArtifactStatus.Available)
      var cleanedCount = 0

      for artifact <- available do
        val cachePath = CacheDir.resolve(artifact.cachePath)
        if !Files.exists(cachePath) then
          artifactRepo.markDeleted(artifact.id)
          cleanedCount += 1

      cleanedCount
    }

  override def getStaleArtifacts(): Either[String, List[AnalysisArtifactEntity]] =
    transactor.readOnly {
      artifactRepo.findStale()
    }

  // ============================================
  // Cache Statistics
  // ============================================

  override def getCacheStats(): Either[String, CacheStats] =
    transactor.readOnly {
      val statusCounts = artifactRepo.countByStatus()
      val (accessible, inaccessible) = sourceFileRepo.countByAccessibility()
      val totalSourceFiles = accessible + inaccessible

      val analyzedCount = sourceFileRepo.findAll().count(_.hasBeenAnalyzed)

      CacheStats(
        totalArtifacts = statusCounts.values.sum.toInt,
        availableArtifacts = statusCounts.getOrElse(ArtifactStatus.Available, 0L).toInt,
        staleArtifacts = statusCounts.getOrElse(ArtifactStatus.Stale, 0L).toInt,
        inProgressArtifacts = statusCounts.getOrElse(ArtifactStatus.InProgress, 0L).toInt,
        errorArtifacts = statusCounts.getOrElse(ArtifactStatus.Error, 0L).toInt,
        totalCacheSizeBytes = artifactRepo.totalCacheSize(),
        trackedSourceFiles = totalSourceFiles.toInt,
        accessibleSourceFiles = accessible.toInt,
        inaccessibleSourceFiles = inaccessible.toInt,
        analyzedSourceFiles = analyzedCount
      )
    }
