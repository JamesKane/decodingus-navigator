package com.decodingus.service

import com.decodingus.db.Transactor
import com.decodingus.repository.{
  BiosampleRepository, SequenceRunRepository, AlignmentRepository
}
import com.decodingus.service.EntityConversions.*
import com.decodingus.util.Logger
import com.decodingus.workspace.model.{
  Biosample, SequenceRun, Alignment, FileInfo, AlignmentMetrics, RecordMeta, Workspace
}
import java.util.UUID

/**
 * Centralized manager for all SequenceRun and Alignment CRUD operations.
 *
 * This service is the SINGLE SOURCE OF TRUTH for sequence data persistence.
 * It ensures:
 * - All operations are transactional
 * - H2 is persisted first (source of truth)
 * - In-memory workspace state is updated via callback after successful persistence
 * - Consistent error handling with Either return types
 * - Proper URI management (local -> persisted)
 *
 * Usage:
 * ```
 * val manager = SequenceDataManager(transactor, biosampleRepo, sequenceRunRepo, alignmentRepo)
 * manager.setWorkspaceUpdater(workspace => _workspace.value = workspace)
 *
 * manager.createSequenceRun(biosampleId, initialData, fileInfo) match {
 *   case Right(result) => // Success, workspace already updated
 *   case Left(error) => // Handle error, workspace unchanged
 * }
 * ```
 */
class SequenceDataManager(
  transactor: Transactor,
  biosampleRepo: BiosampleRepository,
  sequenceRunRepo: SequenceRunRepository,
  alignmentRepo: AlignmentRepository
):

  private val log = Logger[SequenceDataManager]

  // Callback to update in-memory workspace after successful persistence
  private var workspaceUpdater: Workspace => Unit = _ => ()
  private var workspaceGetter: () => Workspace = () => Workspace.empty

  /**
   * Set the callback for updating in-memory workspace state.
   * This should be called once during initialization.
   */
  def setWorkspaceCallbacks(
    getter: () => Workspace,
    updater: Workspace => Unit
  ): Unit =
    workspaceGetter = getter
    workspaceUpdater = updater

  // ============================================
  // SequenceRun Operations
  // ============================================

  /**
   * Result of creating a SequenceRun - includes the persisted run and its index.
   */
  case class CreateSequenceRunResult(
    sequenceRun: SequenceRun,
    index: Int,
    biosampleId: UUID
  )

  /**
   * Creates a new SequenceRun for a biosample.
   *
   * This method:
   * 1. Validates the biosample exists
   * 2. Checks for duplicate files by checksum
   * 3. Persists to H2
   * 4. Updates in-memory workspace state
   *
   * @param sampleAccession The biosample's accession ID
   * @param initialRun Initial SequenceRun data (can have placeholder values)
   * @param fileInfo The BAM/CRAM file info
   * @return Either error message or CreateSequenceRunResult
   */
  def createSequenceRun(
    sampleAccession: String,
    initialRun: SequenceRun,
    fileInfo: FileInfo
  ): Either[String, CreateSequenceRunResult] =
    transactor.readWrite {
      // Step 1: Find biosample
      val biosampleOpt = biosampleRepo.findByAccession(sampleAccession)
      biosampleOpt match
        case None =>
          throw new IllegalArgumentException(s"Biosample not found: $sampleAccession")

        case Some(biosampleEntity) =>
          val biosampleId = biosampleEntity.id

          // Step 2: Check for duplicate files
          val existingRuns = sequenceRunRepo.findByBiosample(biosampleId)
          val existingChecksums = existingRuns.flatMap(_.files.flatMap(_.checksum)).toSet

          if fileInfo.checksum.exists(existingChecksums.contains) then
            throw new IllegalArgumentException(s"Duplicate file detected: ${fileInfo.fileName}")

          // Step 3: Prepare the SequenceRun with file
          val runWithFile = initialRun.copy(
            files = initialRun.files :+ fileInfo,
            biosampleRef = localUri("biosample", biosampleId)
          )

          // Step 4: Convert to entity and persist
          val entity = toSequenceRunEntity(runWithFile, biosampleId)
          val savedEntity = sequenceRunRepo.insert(entity)

          // Step 5: Convert back to domain model with proper URI
          val biosampleRef = localUri("biosample", biosampleId)
          val persistedRun = fromSequenceRunEntity(savedEntity, biosampleRef)

          // Step 6: Calculate index (position in biosample's runs)
          val newIndex = existingRuns.size

          log.info(s"SequenceRun created: ${persistedRun.atUri} for biosample $sampleAccession")

          CreateSequenceRunResult(persistedRun, newIndex, biosampleId)
    }.map { result =>
      // Step 7: Update in-memory workspace after successful persist
      updateWorkspaceWithNewSequenceRun(result.sequenceRun, result.biosampleId)
      result
    }

  /**
   * Updates an existing SequenceRun.
   *
   * @param updatedRun The updated SequenceRun (must have valid atUri)
   * @return Either error message or the updated SequenceRun
   */
  def updateSequenceRun(updatedRun: SequenceRun): Either[String, SequenceRun] =
    val runIdOpt = updatedRun.atUri.flatMap(parseIdFromRef)

    runIdOpt match
      case None =>
        Left("SequenceRun has no valid ID")

      case Some(runId) =>
        transactor.readWrite {
          // Find existing to get biosampleId
          sequenceRunRepo.findById(runId) match
            case None =>
              throw new IllegalArgumentException(s"SequenceRun not found: $runId")

            case Some(existing) =>
              val updatedWithMeta = updatedRun.copy(
                meta = updatedRun.meta.updated("edit")
              )
              val entity = toSequenceRunEntity(updatedWithMeta, existing.biosampleId, Some(runId))
              val savedEntity = sequenceRunRepo.update(entity)
              val biosampleRef = localUri("biosample", existing.biosampleId)
              fromSequenceRunEntity(savedEntity, biosampleRef)
        }.map { persisted =>
          updateWorkspaceWithUpdatedSequenceRun(persisted)
          persisted
        }

  /**
   * Deletes a SequenceRun and all its alignments.
   *
   * @param sampleAccession The biosample's accession ID
   * @param sequenceRunUri The URI of the SequenceRun to delete
   * @return Either error message or true if deleted
   */
  def deleteSequenceRun(sampleAccession: String, sequenceRunUri: String): Either[String, Boolean] =
    val runIdOpt = parseIdFromRef(sequenceRunUri)

    runIdOpt match
      case None =>
        Left(s"Invalid SequenceRun URI: $sequenceRunUri")

      case Some(runId) =>
        transactor.readWrite {
          // Verify it exists
          sequenceRunRepo.findById(runId) match
            case None =>
              log.warn(s"SequenceRun not found in H2: $sequenceRunUri")
              false

            case Some(existing) =>
              // CASCADE delete will remove alignments in H2
              val deleted = sequenceRunRepo.delete(runId)
              log.info(s"SequenceRun deleted from H2: $sequenceRunUri (deleted=$deleted)")
              deleted
        }.map { deleted =>
          if deleted then
            updateWorkspaceWithDeletedSequenceRun(sampleAccession, sequenceRunUri)
          deleted
        }

  /**
   * Gets a SequenceRun by its URI.
   */
  def getSequenceRun(uri: String): Either[String, Option[SequenceRun]] =
    parseIdFromRef(uri) match
      case None => Left(s"Invalid URI: $uri")
      case Some(id) =>
        transactor.readOnly {
          sequenceRunRepo.findById(id).map { entity =>
            val biosampleRef = localUri("biosample", entity.biosampleId)
            fromSequenceRunEntity(entity, biosampleRef)
          }
        }

  /**
   * Gets all SequenceRuns for a biosample.
   */
  def getSequenceRunsForBiosample(sampleAccession: String): Either[String, List[SequenceRun]] =
    transactor.readOnly {
      biosampleRepo.findByAccession(sampleAccession) match
        case None => List.empty
        case Some(biosample) =>
          val biosampleRef = localUri("biosample", biosample.id)
          sequenceRunRepo.findByBiosample(biosample.id).map { entity =>
            fromSequenceRunEntity(entity, biosampleRef)
          }
    }

  // ============================================
  // Alignment Operations
  // ============================================

  /**
   * Result of creating an Alignment.
   */
  case class CreateAlignmentResult(
    alignment: Alignment,
    sequenceRunId: UUID
  )

  /**
   * Creates a new Alignment for a SequenceRun.
   *
   * @param sequenceRunUri The parent SequenceRun's URI
   * @param alignment The alignment to create
   * @param fileInfo Optional file to add to the parent SequenceRun
   * @return Either error message or CreateAlignmentResult
   */
  def createAlignment(
    sequenceRunUri: String,
    alignment: Alignment,
    fileInfo: Option[FileInfo] = None
  ): Either[String, CreateAlignmentResult] =
    val seqRunIdOpt = parseIdFromRef(sequenceRunUri)

    seqRunIdOpt match
      case None =>
        Left(s"Invalid SequenceRun URI: $sequenceRunUri")

      case Some(seqRunId) =>
        transactor.readWrite {
          // Verify sequence run exists
          sequenceRunRepo.findById(seqRunId) match
            case None =>
              throw new IllegalArgumentException(s"SequenceRun not found: $seqRunId")

            case Some(seqRunEntity) =>
              // Prepare alignment with proper reference
              val alignmentWithRef = alignment.copy(
                sequenceRunRef = localUri("sequencerun", seqRunId)
              )

              // Persist alignment
              val entity = toAlignmentEntity(alignmentWithRef, seqRunId)
              val savedEntity = alignmentRepo.insert(entity)
              val persistedAlignment = fromAlignmentEntity(savedEntity, localUri("sequencerun", seqRunId))

              // Optionally add file to sequence run
              fileInfo.foreach { file =>
                if !seqRunEntity.files.exists(_.checksum == file.checksum) then
                  sequenceRunRepo.addFile(seqRunId, file)
              }

              log.info(s"Alignment created: ${persistedAlignment.atUri} for SequenceRun $sequenceRunUri")

              CreateAlignmentResult(persistedAlignment, seqRunId)
        }.map { result =>
          updateWorkspaceWithNewAlignment(result.alignment, sequenceRunUri, fileInfo)
          result
        }

  /**
   * Updates an existing Alignment.
   *
   * @param updatedAlignment The updated Alignment (must have valid atUri)
   * @return Either error message or the updated Alignment
   */
  def updateAlignment(updatedAlignment: Alignment): Either[String, Alignment] =
    val alignIdOpt = updatedAlignment.atUri.flatMap(parseIdFromRef)

    alignIdOpt match
      case None =>
        Left("Alignment has no valid ID")

      case Some(alignId) =>
        transactor.readWrite {
          alignmentRepo.findById(alignId) match
            case None =>
              throw new IllegalArgumentException(s"Alignment not found: $alignId")

            case Some(existing) =>
              val updatedWithMeta = updatedAlignment.copy(
                meta = updatedAlignment.meta.updated("edit")
              )
              val entity = toAlignmentEntity(updatedWithMeta, existing.sequenceRunId, Some(alignId))
              val savedEntity = alignmentRepo.update(entity)
              val seqRunRef = localUri("sequencerun", existing.sequenceRunId)
              fromAlignmentEntity(savedEntity, seqRunRef)
        }.map { persisted =>
          updateWorkspaceWithUpdatedAlignment(persisted)
          persisted
        }

  /**
   * Updates alignment metrics.
   *
   * @param alignmentUri The alignment's URI
   * @param metrics The new metrics
   * @return Either error message or true if updated
   */
  def updateAlignmentMetrics(alignmentUri: String, metrics: AlignmentMetrics): Either[String, Boolean] =
    parseIdFromRef(alignmentUri) match
      case None => Left(s"Invalid alignment URI: $alignmentUri")
      case Some(id) =>
        transactor.readWrite {
          alignmentRepo.updateMetrics(id, metrics)
        }.map { updated =>
          if updated then
            updateWorkspaceWithAlignmentMetrics(alignmentUri, metrics)
          updated
        }

  /**
   * Deletes an Alignment.
   *
   * @param alignmentUri The URI of the Alignment to delete
   * @return Either error message or true if deleted
   */
  def deleteAlignment(alignmentUri: String): Either[String, Boolean] =
    parseIdFromRef(alignmentUri) match
      case None => Left(s"Invalid alignment URI: $alignmentUri")
      case Some(id) =>
        transactor.readWrite {
          alignmentRepo.delete(id)
        }.map { deleted =>
          if deleted then
            updateWorkspaceWithDeletedAlignment(alignmentUri)
          deleted
        }

  /**
   * Gets an Alignment by its URI.
   */
  def getAlignment(uri: String): Either[String, Option[Alignment]] =
    parseIdFromRef(uri) match
      case None => Left(s"Invalid URI: $uri")
      case Some(id) =>
        transactor.readOnly {
          alignmentRepo.findById(id).map { entity =>
            val seqRunRef = localUri("sequencerun", entity.sequenceRunId)
            fromAlignmentEntity(entity, seqRunRef)
          }
        }

  /**
   * Gets all Alignments for a SequenceRun.
   */
  def getAlignmentsForSequenceRun(sequenceRunUri: String): Either[String, List[Alignment]] =
    parseIdFromRef(sequenceRunUri) match
      case None => Left(s"Invalid URI: $sequenceRunUri")
      case Some(seqRunId) =>
        transactor.readOnly {
          val seqRunRef = localUri("sequencerun", seqRunId)
          alignmentRepo.findBySequenceRun(seqRunId).map { entity =>
            fromAlignmentEntity(entity, seqRunRef)
          }
        }

  // ============================================
  // Workspace State Update Helpers
  // ============================================

  /**
   * Updates workspace after creating a new SequenceRun.
   */
  private def updateWorkspaceWithNewSequenceRun(
    sequenceRun: SequenceRun,
    biosampleId: UUID
  ): Unit =
    val workspace = workspaceGetter()
    val biosampleUri = localUri("biosample", biosampleId)
    val seqRunUri = sequenceRun.atUri.getOrElse("")

    // Add sequence run to list
    val updatedSequenceRuns = workspace.main.sequenceRuns :+ sequenceRun

    // Update biosample's sequenceRunRefs
    val updatedSamples = workspace.main.samples.map { sample =>
      if sample.atUri.contains(biosampleUri) ||
         sample.atUri.exists(uri => parseIdFromRef(uri).contains(biosampleId)) then
        sample.copy(sequenceRunRefs = sample.sequenceRunRefs :+ seqRunUri)
      else
        sample
    }

    val updatedContent = workspace.main.copy(
      samples = updatedSamples,
      sequenceRuns = updatedSequenceRuns
    )
    workspaceUpdater(workspace.copy(main = updatedContent))

  /**
   * Updates workspace after updating a SequenceRun.
   */
  private def updateWorkspaceWithUpdatedSequenceRun(sequenceRun: SequenceRun): Unit =
    val workspace = workspaceGetter()
    val updatedSequenceRuns = workspace.main.sequenceRuns.map { sr =>
      if sr.atUri == sequenceRun.atUri then sequenceRun else sr
    }
    val updatedContent = workspace.main.copy(sequenceRuns = updatedSequenceRuns)
    workspaceUpdater(workspace.copy(main = updatedContent))

  /**
   * Updates workspace after deleting a SequenceRun.
   */
  private def updateWorkspaceWithDeletedSequenceRun(
    sampleAccession: String,
    sequenceRunUri: String
  ): Unit =
    val workspace = workspaceGetter()

    // Remove the sequence run
    val updatedSequenceRuns = workspace.main.sequenceRuns.filterNot(_.atUri.contains(sequenceRunUri))

    // Remove any alignments referencing this sequence run
    val updatedAlignments = workspace.main.alignments.filterNot(_.sequenceRunRef == sequenceRunUri)

    // Remove from biosample's refs
    val updatedSamples = workspace.main.samples.map { sample =>
      if sample.sampleAccession == sampleAccession then
        sample.copy(sequenceRunRefs = sample.sequenceRunRefs.filterNot(_ == sequenceRunUri))
      else
        sample
    }

    val updatedContent = workspace.main.copy(
      samples = updatedSamples,
      sequenceRuns = updatedSequenceRuns,
      alignments = updatedAlignments
    )
    workspaceUpdater(workspace.copy(main = updatedContent))

  /**
   * Updates workspace after creating a new Alignment.
   */
  private def updateWorkspaceWithNewAlignment(
    alignment: Alignment,
    sequenceRunUri: String,
    fileInfo: Option[FileInfo]
  ): Unit =
    val workspace = workspaceGetter()
    val alignUri = alignment.atUri.getOrElse("")

    // Add alignment to list
    val updatedAlignments = workspace.main.alignments :+ alignment

    // Update sequence run's alignmentRefs and optionally add file
    val updatedSequenceRuns = workspace.main.sequenceRuns.map { sr =>
      if sr.atUri.contains(sequenceRunUri) then
        val newAlignRefs = if sr.alignmentRefs.contains(alignUri) then sr.alignmentRefs
                          else sr.alignmentRefs :+ alignUri
        val newFiles = fileInfo match
          case Some(file) if !sr.files.exists(_.checksum == file.checksum) => sr.files :+ file
          case _ => sr.files
        sr.copy(alignmentRefs = newAlignRefs, files = newFiles)
      else
        sr
    }

    val updatedContent = workspace.main.copy(
      sequenceRuns = updatedSequenceRuns,
      alignments = updatedAlignments
    )
    workspaceUpdater(workspace.copy(main = updatedContent))

  /**
   * Updates workspace after updating an Alignment.
   */
  private def updateWorkspaceWithUpdatedAlignment(alignment: Alignment): Unit =
    val workspace = workspaceGetter()
    val updatedAlignments = workspace.main.alignments.map { a =>
      if a.atUri == alignment.atUri then alignment else a
    }
    val updatedContent = workspace.main.copy(alignments = updatedAlignments)
    workspaceUpdater(workspace.copy(main = updatedContent))

  /**
   * Updates workspace after updating alignment metrics.
   */
  private def updateWorkspaceWithAlignmentMetrics(alignmentUri: String, metrics: AlignmentMetrics): Unit =
    val workspace = workspaceGetter()
    val updatedAlignments = workspace.main.alignments.map { a =>
      if a.atUri.contains(alignmentUri) then a.copy(metrics = Some(metrics)) else a
    }
    val updatedContent = workspace.main.copy(alignments = updatedAlignments)
    workspaceUpdater(workspace.copy(main = updatedContent))

  /**
   * Updates workspace after deleting an Alignment.
   */
  private def updateWorkspaceWithDeletedAlignment(alignmentUri: String): Unit =
    val workspace = workspaceGetter()

    // Remove the alignment
    val updatedAlignments = workspace.main.alignments.filterNot(_.atUri.contains(alignmentUri))

    // Remove from sequence run's refs
    val updatedSequenceRuns = workspace.main.sequenceRuns.map { sr =>
      sr.copy(alignmentRefs = sr.alignmentRefs.filterNot(_ == alignmentUri))
    }

    val updatedContent = workspace.main.copy(
      sequenceRuns = updatedSequenceRuns,
      alignments = updatedAlignments
    )
    workspaceUpdater(workspace.copy(main = updatedContent))

object SequenceDataManager:
  /**
   * Create a SequenceDataManager from a DatabaseContext.
   */
  def apply(context: DatabaseContext): SequenceDataManager =
    new SequenceDataManager(
      context.transactor,
      context.biosampleRepository,
      context.sequenceRunRepository,
      context.alignmentRepository
    )
