package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import io.circe.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Types of analysis artifacts that can be cached.
 */
enum ArtifactType:
  case WgsMetrics
  case CallableLoci
  case HaplogroupVcf
  case WholeGenomeVcf
  case PrivateVariants
  case InsertSizeMetrics
  case AlignmentSummary
  case CoverageSummary
  case SexInference
  case DuplicateMetrics

object ArtifactType:
  def fromString(s: String): ArtifactType = s match
    case "WGS_METRICS" => WgsMetrics
    case "CALLABLE_LOCI" => CallableLoci
    case "HAPLOGROUP_VCF" => HaplogroupVcf
    case "WHOLE_GENOME_VCF" => WholeGenomeVcf
    case "PRIVATE_VARIANTS" => PrivateVariants
    case "INSERT_SIZE_METRICS" => InsertSizeMetrics
    case "ALIGNMENT_SUMMARY" => AlignmentSummary
    case "COVERAGE_SUMMARY" => CoverageSummary
    case "SEX_INFERENCE" => SexInference
    case "DUPLICATE_METRICS" => DuplicateMetrics
    case other => throw new IllegalArgumentException(s"Unknown artifact type: $other")

  def toDbString(t: ArtifactType): String = t match
    case WgsMetrics => "WGS_METRICS"
    case CallableLoci => "CALLABLE_LOCI"
    case HaplogroupVcf => "HAPLOGROUP_VCF"
    case WholeGenomeVcf => "WHOLE_GENOME_VCF"
    case PrivateVariants => "PRIVATE_VARIANTS"
    case InsertSizeMetrics => "INSERT_SIZE_METRICS"
    case AlignmentSummary => "ALIGNMENT_SUMMARY"
    case CoverageSummary => "COVERAGE_SUMMARY"
    case SexInference => "SEX_INFERENCE"
    case DuplicateMetrics => "DUPLICATE_METRICS"

/**
 * Status of a cached artifact.
 */
enum ArtifactStatus:
  case Available   // Ready to use
  case InProgress  // Currently being generated
  case Stale       // Dependencies changed, needs regeneration
  case Deleted     // File was deleted
  case Error       // Generation failed

object ArtifactStatus:
  def fromString(s: String): ArtifactStatus = s match
    case "AVAILABLE" => Available
    case "IN_PROGRESS" => InProgress
    case "STALE" => Stale
    case "DELETED" => Deleted
    case "ERROR" => Error
    case other => throw new IllegalArgumentException(s"Unknown artifact status: $other")

  def toDbString(s: ArtifactStatus): String = s match
    case Available => "AVAILABLE"
    case InProgress => "IN_PROGRESS"
    case Stale => "STALE"
    case Deleted => "DELETED"
    case Error => "ERROR"

/**
 * Analysis artifact entity for database persistence.
 *
 * Tracks cached analysis outputs (WGS metrics, VCFs, callable loci, etc.)
 * linked to their source alignment.
 */
case class AnalysisArtifactEntity(
  id: UUID,
  alignmentId: UUID,
  artifactType: ArtifactType,
  cachePath: String,
  fileSize: Option[Long],
  fileChecksum: Option[String],
  fileFormat: Option[String],
  generatedAt: LocalDateTime,
  generatorVersion: Option[String],
  generationParams: Option[Json],
  status: ArtifactStatus,
  staleReason: Option[String],
  dependsOnSourceChecksum: Option[String],
  dependsOnReferenceBuild: Option[String],
  createdAt: LocalDateTime,
  updatedAt: LocalDateTime
) extends Entity[UUID]

object AnalysisArtifactEntity:
  /**
   * Create a new artifact entity.
   */
  def create(
    alignmentId: UUID,
    artifactType: ArtifactType,
    cachePath: String,
    generatorVersion: Option[String] = None,
    generationParams: Option[Json] = None,
    dependsOnSourceChecksum: Option[String] = None,
    dependsOnReferenceBuild: Option[String] = None
  ): AnalysisArtifactEntity =
    val now = LocalDateTime.now()
    AnalysisArtifactEntity(
      id = UUID.randomUUID(),
      alignmentId = alignmentId,
      artifactType = artifactType,
      cachePath = cachePath,
      fileSize = None,
      fileChecksum = None,
      fileFormat = None,
      generatedAt = now,
      generatorVersion = generatorVersion,
      generationParams = generationParams,
      status = ArtifactStatus.InProgress,
      staleReason = None,
      dependsOnSourceChecksum = dependsOnSourceChecksum,
      dependsOnReferenceBuild = dependsOnReferenceBuild,
      createdAt = now,
      updatedAt = now
    )

/**
 * Repository for analysis artifact persistence.
 */
class AnalysisArtifactRepository:

  // ============================================
  // Core Operations
  // ============================================

  def findById(id: UUID)(using conn: Connection): Option[AnalysisArtifactEntity] =
    queryOne(
      "SELECT * FROM analysis_artifact WHERE id = ?",
      Seq(id)
    )(mapRow)

  def findAll()(using conn: Connection): List[AnalysisArtifactEntity] =
    queryList("SELECT * FROM analysis_artifact ORDER BY generated_at DESC")(mapRow)

  def insert(entity: AnalysisArtifactEntity)(using conn: Connection): AnalysisArtifactEntity =
    val paramsJson = entity.generationParams.map(j => JsonValue(j.noSpaces))

    executeUpdate(
      """INSERT INTO analysis_artifact (
        |  id, alignment_id, artifact_type, cache_path, file_size, file_checksum,
        |  file_format, generated_at, generator_version, generation_params, status,
        |  stale_reason, depends_on_source_checksum, depends_on_reference_build,
        |  created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.alignmentId,
        ArtifactType.toDbString(entity.artifactType),
        entity.cachePath,
        entity.fileSize,
        entity.fileChecksum,
        entity.fileFormat,
        entity.generatedAt,
        entity.generatorVersion,
        paramsJson,
        ArtifactStatus.toDbString(entity.status),
        entity.staleReason,
        entity.dependsOnSourceChecksum,
        entity.dependsOnReferenceBuild,
        entity.createdAt,
        entity.updatedAt
      )
    )
    entity

  def update(entity: AnalysisArtifactEntity)(using conn: Connection): AnalysisArtifactEntity =
    val paramsJson = entity.generationParams.map(j => JsonValue(j.noSpaces))
    val now = LocalDateTime.now()

    executeUpdate(
      """UPDATE analysis_artifact SET
        |  alignment_id = ?, artifact_type = ?, cache_path = ?, file_size = ?,
        |  file_checksum = ?, file_format = ?, generated_at = ?, generator_version = ?,
        |  generation_params = ?, status = ?, stale_reason = ?,
        |  depends_on_source_checksum = ?, depends_on_reference_build = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.alignmentId,
        ArtifactType.toDbString(entity.artifactType),
        entity.cachePath,
        entity.fileSize,
        entity.fileChecksum,
        entity.fileFormat,
        entity.generatedAt,
        entity.generatorVersion,
        paramsJson,
        ArtifactStatus.toDbString(entity.status),
        entity.staleReason,
        entity.dependsOnSourceChecksum,
        entity.dependsOnReferenceBuild,
        now,
        entity.id
      )
    )
    entity.copy(updatedAt = now)

  def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM analysis_artifact WHERE id = ?", Seq(id)) > 0

  // ============================================
  // Query Operations
  // ============================================

  /**
   * Find all artifacts for an alignment.
   */
  def findByAlignment(alignmentId: UUID)(using conn: Connection): List[AnalysisArtifactEntity] =
    queryList(
      "SELECT * FROM analysis_artifact WHERE alignment_id = ? ORDER BY artifact_type",
      Seq(alignmentId)
    )(mapRow)

  /**
   * Find artifact by alignment and type (unique combination).
   */
  def findByAlignmentAndType(
    alignmentId: UUID,
    artifactType: ArtifactType
  )(using conn: Connection): Option[AnalysisArtifactEntity] =
    queryOne(
      "SELECT * FROM analysis_artifact WHERE alignment_id = ? AND artifact_type = ?",
      Seq(alignmentId, ArtifactType.toDbString(artifactType))
    )(mapRow)

  /**
   * Find artifacts by status.
   */
  def findByStatus(status: ArtifactStatus)(using conn: Connection): List[AnalysisArtifactEntity] =
    queryList(
      "SELECT * FROM analysis_artifact WHERE status = ? ORDER BY updated_at DESC",
      Seq(ArtifactStatus.toDbString(status))
    )(mapRow)

  /**
   * Find artifacts by type.
   */
  def findByType(artifactType: ArtifactType)(using conn: Connection): List[AnalysisArtifactEntity] =
    queryList(
      "SELECT * FROM analysis_artifact WHERE artifact_type = ? ORDER BY generated_at DESC",
      Seq(ArtifactType.toDbString(artifactType))
    )(mapRow)

  /**
   * Find stale artifacts that need regeneration.
   */
  def findStale()(using conn: Connection): List[AnalysisArtifactEntity] =
    queryList(
      "SELECT * FROM analysis_artifact WHERE status = 'STALE' ORDER BY updated_at ASC"
    )(mapRow)

  /**
   * Find artifacts in progress.
   */
  def findInProgress()(using conn: Connection): List[AnalysisArtifactEntity] =
    queryList(
      "SELECT * FROM analysis_artifact WHERE status = 'IN_PROGRESS' ORDER BY generated_at ASC"
    )(mapRow)

  // ============================================
  // Status Update Operations
  // ============================================

  /**
   * Mark an artifact as available (generation complete).
   */
  def markAvailable(
    id: UUID,
    fileSize: Long,
    fileChecksum: String,
    fileFormat: Option[String] = None
  )(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE analysis_artifact SET
        |  status = 'AVAILABLE', file_size = ?, file_checksum = ?, file_format = ?,
        |  stale_reason = NULL, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(fileSize, fileChecksum, fileFormat, LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark an artifact as stale.
   */
  def markStale(id: UUID, reason: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE analysis_artifact SET
        |  status = 'STALE', stale_reason = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(reason, LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark an artifact as having an error.
   */
  def markError(id: UUID, reason: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE analysis_artifact SET
        |  status = 'ERROR', stale_reason = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(reason, LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark an artifact as deleted.
   */
  def markDeleted(id: UUID)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE analysis_artifact SET
        |  status = 'DELETED', updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark artifacts as stale when source checksum changes.
   */
  def markStaleBySourceChecksum(oldChecksum: String, reason: String)(using conn: Connection): Int =
    executeUpdate(
      """UPDATE analysis_artifact SET
        |  status = 'STALE', stale_reason = ?, updated_at = ?
        |WHERE depends_on_source_checksum = ? AND status = 'AVAILABLE'
      """.stripMargin,
      Seq(reason, LocalDateTime.now(), oldChecksum)
    )

  // ============================================
  // Statistics
  // ============================================

  /**
   * Count artifacts by status.
   */
  def countByStatus()(using conn: Connection): Map[ArtifactStatus, Long] =
    queryList(
      "SELECT status, COUNT(*) as cnt FROM analysis_artifact GROUP BY status"
    ) { rs =>
      val status = ArtifactStatus.fromString(rs.getString("status"))
      val count = rs.getLong("cnt")
      (status, count)
    }.toMap

  /**
   * Get total cache size in bytes.
   */
  def totalCacheSize()(using conn: Connection): Long =
    queryOne(
      "SELECT COALESCE(SUM(file_size), 0) FROM analysis_artifact WHERE status = 'AVAILABLE'"
    ) { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): AnalysisArtifactEntity =
    val paramsJson = getOptJsonString(rs, "generation_params")
    val params = paramsJson.flatMap(json => parse(json).toOption)

    AnalysisArtifactEntity(
      id = getUUID(rs, "id"),
      alignmentId = getUUID(rs, "alignment_id"),
      artifactType = ArtifactType.fromString(rs.getString("artifact_type")),
      cachePath = rs.getString("cache_path"),
      fileSize = getOptLong(rs, "file_size"),
      fileChecksum = getOptString(rs, "file_checksum"),
      fileFormat = getOptString(rs, "file_format"),
      generatedAt = getDateTime(rs, "generated_at"),
      generatorVersion = getOptString(rs, "generator_version"),
      generationParams = params,
      status = ArtifactStatus.fromString(rs.getString("status")),
      staleReason = getOptString(rs, "stale_reason"),
      dependsOnSourceChecksum = getOptString(rs, "depends_on_source_checksum"),
      dependsOnReferenceBuild = getOptString(rs, "depends_on_reference_build"),
      createdAt = getDateTime(rs, "created_at"),
      updatedAt = getDateTime(rs, "updated_at")
    )
