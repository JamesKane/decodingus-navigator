package com.decodingus.service

import com.decodingus.repository.{AlignmentEntity, BiosampleEntity, ChipProfileEntity, EntityMeta, ProjectEntity, SequenceRunEntity, StrProfileEntity, SyncStatus as RepoSyncStatus}
import com.decodingus.workspace.model.*

import java.time.LocalDateTime
import java.util.UUID

/**
 * Conversions between domain models and database entities.
 *
 * Domain models (Biosample, SequenceRun, etc.) are what the UI and business logic use.
 * Entity models (BiosampleEntity, etc.) are what's stored in the database.
 *
 * Key differences:
 * - Entities have UUID primary keys
 * - Entities have EntityMeta (sync status, version, timestamps)
 * - Domain models use AT URIs for references; entities use foreign key UUIDs
 */
object EntityConversions:

  // ============================================
  // ID Management
  // ============================================

  /**
   * Extract UUID from an AT URI or generate a new one.
   * AT URIs have format: at://did/collection/rkey
   * For local records without AT URI, we use the last segment as UUID if valid.
   */
  def extractOrGenerateId(atUri: Option[String]): UUID =
    atUri.flatMap(extractIdFromUri).getOrElse(UUID.randomUUID())

  private def extractIdFromUri(uri: String): Option[UUID] =
    try
      val parts = uri.split("/")
      if parts.length >= 3 then
        Some(UUID.fromString(parts.last))
      else
        None
    catch
      case _: IllegalArgumentException => None

  /**
   * Create an AT URI from entity type and ID.
   * Uses local:// scheme for entities not yet synced to PDS.
   */
  def localUri(entityType: String, id: UUID): String =
    s"local://$entityType/$id"

  // ============================================
  // Biosample Conversions
  // ============================================

  def toBiosampleEntity(biosample: Biosample, existingId: Option[UUID] = None): BiosampleEntity =
    val id = existingId.getOrElse(extractOrGenerateId(biosample.atUri))
    BiosampleEntity(
      id = id,
      sampleAccession = biosample.sampleAccession,
      donorIdentifier = biosample.donorIdentifier,
      description = biosample.description,
      centerName = biosample.centerName,
      sex = biosample.sex,
      citizenDid = biosample.citizenDid,
      haplogroups = biosample.haplogroups,
      meta = toEntityMeta(biosample.meta, biosample.atUri)
    )

  def fromBiosampleEntity(entity: BiosampleEntity): Biosample =
    Biosample(
      atUri = entity.meta.atUri.orElse(Some(localUri("biosample", entity.id))),
      meta = fromEntityMeta(entity.meta),
      sampleAccession = entity.sampleAccession,
      donorIdentifier = entity.donorIdentifier,
      citizenDid = entity.citizenDid,
      description = entity.description,
      centerName = entity.centerName,
      sex = entity.sex,
      haplogroups = entity.haplogroups,
      // Refs are populated separately via queries
      sequenceRunRefs = List.empty,
      genotypeRefs = List.empty,
      populationBreakdownRef = None,
      strProfileRef = None,
      strProfileRefs = List.empty,
      ySnpPanelRefs = List.empty
    )

  // ============================================
  // Project Conversions
  // ============================================

  def toProjectEntity(project: Project, existingId: Option[UUID] = None): ProjectEntity =
    val id = existingId.getOrElse(extractOrGenerateId(project.atUri))
    ProjectEntity(
      id = id,
      projectName = project.projectName,
      description = project.description,
      administratorDid = project.administrator,
      meta = toEntityMeta(project.meta, project.atUri)
    )

  def fromProjectEntity(entity: ProjectEntity, memberRefs: List[String] = List.empty): Project =
    Project(
      atUri = entity.meta.atUri.orElse(Some(localUri("project", entity.id))),
      meta = fromEntityMeta(entity.meta),
      projectName = entity.projectName,
      description = entity.description,
      administrator = entity.administratorDid,
      memberRefs = memberRefs
    )

  // ============================================
  // SequenceRun Conversions
  // ============================================

  /**
   * Normalize platform name to match DB constraint values.
   * Analysis returns values like "Illumina", "PacBio" but DB requires "ILLUMINA", "PACBIO".
   */
  def normalizePlatform(platform: String): String =
    platform.toUpperCase match
      case "ILLUMINA" | "ILLUMINA/SOLEXA" => "ILLUMINA"
      case "PACBIO" | "PACIFIC BIOSCIENCES" => "PACBIO"
      case "NANOPORE" | "OXFORD NANOPORE" => "NANOPORE"
      case "ION_TORRENT" | "ION TORRENT" | "IONTORRENT" => "ION_TORRENT"
      case "BGI" | "MGI" | "BGISEQ" | "MGISEQ" => "BGI"
      case "ELEMENT" | "ELEMENT BIOSCIENCES" => "ELEMENT"
      case "ULTIMA" | "ULTIMA GENOMICS" => "ULTIMA"
      case _ => "Unknown"

  /**
   * Normalize library layout to match DB constraint values.
   * Analysis returns "Paired-End"/"Single-End" but DB requires "PAIRED"/"SINGLE".
   */
  def normalizeLibraryLayout(layout: Option[String]): Option[String] =
    layout.map(_.toUpperCase match
      case "PAIRED-END" | "PAIRED" | "PE" => "PAIRED"
      case "SINGLE-END" | "SINGLE" | "SE" => "SINGLE"
      case other => other // Pass through unknown values
    )

  def toSequenceRunEntity(
                           sequenceRun: SequenceRun,
                           biosampleId: UUID,
                           existingId: Option[UUID] = None
                         ): SequenceRunEntity =
    val id = existingId.getOrElse(extractOrGenerateId(sequenceRun.atUri))
    SequenceRunEntity(
      id = id,
      biosampleId = biosampleId,
      platform = normalizePlatform(sequenceRun.platformName),
      instrumentModel = sequenceRun.instrumentModel,
      instrumentId = sequenceRun.instrumentId,
      testType = sequenceRun.testType,
      libraryId = sequenceRun.libraryId,
      platformUnit = sequenceRun.platformUnit,
      libraryLayout = normalizeLibraryLayout(sequenceRun.libraryLayout),
      sampleName = sequenceRun.sampleName,
      sequencingFacility = sequenceRun.sequencingFacility,
      runFingerprint = sequenceRun.runFingerprint,
      totalReads = sequenceRun.totalReads,
      pfReads = sequenceRun.pfReads,
      pfReadsAligned = sequenceRun.pfReadsAligned,
      pctPfReadsAligned = sequenceRun.pctPfReadsAligned,
      readsPaired = sequenceRun.readsPaired,
      pctReadsPaired = sequenceRun.pctReadsPaired,
      pctProperPairs = sequenceRun.pctProperPairs,
      readLength = sequenceRun.readLength,
      maxReadLength = sequenceRun.maxReadLength,
      meanInsertSize = sequenceRun.meanInsertSize,
      medianInsertSize = sequenceRun.medianInsertSize,
      stdInsertSize = sequenceRun.stdInsertSize,
      pairOrientation = sequenceRun.pairOrientation,
      flowcellId = sequenceRun.flowcellId,
      runDate = sequenceRun.runDate,
      files = sequenceRun.files,
      meta = toEntityMeta(sequenceRun.meta, sequenceRun.atUri)
    )

  def fromSequenceRunEntity(entity: SequenceRunEntity, biosampleRef: String): SequenceRun =
    SequenceRun(
      atUri = entity.meta.atUri.orElse(Some(localUri("sequencerun", entity.id))),
      meta = fromEntityMeta(entity.meta),
      biosampleRef = biosampleRef,
      platformName = entity.platform,
      instrumentModel = entity.instrumentModel,
      instrumentId = entity.instrumentId,
      sequencingFacility = entity.sequencingFacility,
      sampleName = entity.sampleName,
      libraryId = entity.libraryId,
      platformUnit = entity.platformUnit,
      runFingerprint = entity.runFingerprint,
      testType = entity.testType,
      libraryLayout = entity.libraryLayout,
      totalReads = entity.totalReads,
      pfReads = entity.pfReads,
      pfReadsAligned = entity.pfReadsAligned,
      pctPfReadsAligned = entity.pctPfReadsAligned,
      readsPaired = entity.readsPaired,
      pctReadsPaired = entity.pctReadsPaired,
      pctProperPairs = entity.pctProperPairs,
      readLength = entity.readLength,
      maxReadLength = entity.maxReadLength,
      meanInsertSize = entity.meanInsertSize,
      medianInsertSize = entity.medianInsertSize,
      stdInsertSize = entity.stdInsertSize,
      pairOrientation = entity.pairOrientation,
      flowcellId = entity.flowcellId,
      runDate = entity.runDate,
      files = entity.files,
      // Refs populated separately
      alignmentRefs = List.empty
    )

  // ============================================
  // Alignment Conversions
  // ============================================

  def toAlignmentEntity(
                         alignment: Alignment,
                         sequenceRunId: UUID,
                         existingId: Option[UUID] = None
                       ): AlignmentEntity =
    val id = existingId.getOrElse(extractOrGenerateId(alignment.atUri))
    AlignmentEntity(
      id = id,
      sequenceRunId = sequenceRunId,
      referenceBuild = alignment.referenceBuild,
      aligner = alignment.aligner,
      variantCaller = alignment.variantCaller,
      metrics = alignment.metrics,
      files = alignment.files,
      meta = toEntityMeta(alignment.meta, alignment.atUri)
    )

  def fromAlignmentEntity(entity: AlignmentEntity, sequenceRunRef: String): Alignment =
    Alignment(
      atUri = entity.meta.atUri.orElse(Some(localUri("alignment", entity.id))),
      meta = fromEntityMeta(entity.meta),
      sequenceRunRef = sequenceRunRef,
      biosampleRef = None, // Populated via query if needed
      referenceBuild = entity.referenceBuild,
      aligner = entity.aligner,
      variantCaller = entity.variantCaller,
      files = entity.files,
      metrics = entity.metrics
    )

  // ============================================
  // STR Profile Conversions
  // ============================================

  def toStrProfileEntity(
                          profile: StrProfile,
                          biosampleId: UUID,
                          existingId: Option[UUID] = None
                        ): StrProfileEntity =
    val id = existingId.getOrElse(extractOrGenerateId(profile.atUri))
    StrProfileEntity(
      id = id,
      biosampleId = biosampleId,
      sequenceRunId = profile.sequenceRunRef.flatMap(parseIdFromRef),
      panels = profile.panels,
      markers = profile.markers,
      totalMarkers = profile.totalMarkers,
      source = profile.source,
      importedFrom = profile.importedFrom,
      derivationMethod = profile.derivationMethod,
      files = profile.files,
      meta = toEntityMeta(profile.meta, profile.atUri)
    )

  def fromStrProfileEntity(entity: StrProfileEntity, biosampleRef: String): StrProfile =
    StrProfile(
      atUri = entity.meta.atUri.orElse(Some(localUri("strprofile", entity.id))),
      meta = fromEntityMeta(entity.meta),
      biosampleRef = biosampleRef,
      sequenceRunRef = entity.sequenceRunId.map(id => localUri("sequencerun", id)),
      panels = entity.panels,
      markers = entity.markers,
      totalMarkers = entity.totalMarkers,
      source = entity.source,
      importedFrom = entity.importedFrom,
      derivationMethod = entity.derivationMethod,
      files = entity.files
    )

  // ============================================
  // Chip Profile Conversions
  // ============================================

  def toChipProfileEntity(
                           profile: ChipProfile,
                           biosampleId: UUID,
                           existingId: Option[UUID] = None
                         ): ChipProfileEntity =
    val id = existingId.getOrElse(extractOrGenerateId(profile.atUri))
    ChipProfileEntity(
      id = id,
      biosampleId = biosampleId,
      vendor = profile.vendor,
      testTypeCode = profile.testTypeCode,
      chipVersion = profile.chipVersion,
      totalMarkersCalled = profile.totalMarkersCalled,
      totalMarkersPossible = profile.totalMarkersPossible,
      noCallRate = profile.noCallRate,
      yMarkersCalled = profile.yMarkersCalled,
      mtMarkersCalled = profile.mtMarkersCalled,
      autosomalMarkersCalled = profile.autosomalMarkersCalled,
      hetRate = profile.hetRate,
      importDate = profile.importDate,
      sourceFileHash = profile.sourceFileHash,
      sourceFileName = profile.sourceFileName,
      files = profile.files,
      meta = toEntityMeta(profile.meta, profile.atUri)
    )

  def fromChipProfileEntity(entity: ChipProfileEntity, biosampleRef: String): ChipProfile =
    ChipProfile(
      atUri = entity.meta.atUri.orElse(Some(localUri("chipprofile", entity.id))),
      meta = fromEntityMeta(entity.meta),
      biosampleRef = biosampleRef,
      vendor = entity.vendor,
      testTypeCode = entity.testTypeCode,
      chipVersion = entity.chipVersion,
      totalMarkersCalled = entity.totalMarkersCalled,
      totalMarkersPossible = entity.totalMarkersPossible,
      noCallRate = entity.noCallRate,
      yMarkersCalled = entity.yMarkersCalled,
      mtMarkersCalled = entity.mtMarkersCalled,
      autosomalMarkersCalled = entity.autosomalMarkersCalled,
      hetRate = entity.hetRate,
      importDate = entity.importDate,
      sourceFileHash = entity.sourceFileHash,
      sourceFileName = entity.sourceFileName,
      files = entity.files
    )

  // ============================================
  // Metadata Conversions
  // ============================================

  private def toEntityMeta(recordMeta: RecordMeta, atUri: Option[String]): EntityMeta =
    EntityMeta(
      syncStatus = if atUri.exists(_.startsWith("at://")) then RepoSyncStatus.Synced else RepoSyncStatus.Local,
      atUri = atUri.filter(_.startsWith("at://")), // Only real AT URIs
      atCid = None, // Set by sync process
      version = recordMeta.version,
      createdAt = recordMeta.createdAt,
      updatedAt = recordMeta.updatedAt.getOrElse(recordMeta.createdAt)
    )

  private def fromEntityMeta(entityMeta: EntityMeta): RecordMeta =
    RecordMeta(
      version = entityMeta.version,
      createdAt = entityMeta.createdAt,
      updatedAt = Some(entityMeta.updatedAt),
      lastModifiedField = None
    )

  // ============================================
  // ID Lookup Helpers
  // ============================================

  /**
   * Parse a UUID from a local URI or AT URI.
   */
  def parseIdFromRef(ref: String): Option[UUID] =
    if ref.startsWith("local://") then
      try Some(UUID.fromString(ref.split("/").last))
      catch case _: Exception => None
    else if ref.startsWith("at://") then
      extractIdFromUri(ref)
    else
      try Some(UUID.fromString(ref))
      catch case _: Exception => None
