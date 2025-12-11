package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * File format for source files.
 */
enum SourceFileFormat:
  case Bam
  case Cram
  case Fastq
  case Vcf
  case Gvcf

object SourceFileFormat:
  def fromString(s: String): SourceFileFormat = s.toUpperCase match
    case "BAM" => Bam
    case "CRAM" => Cram
    case "FASTQ" => Fastq
    case "VCF" => Vcf
    case "GVCF" => Gvcf
    case other => throw new IllegalArgumentException(s"Unknown file format: $other")

  def toDbString(f: SourceFileFormat): String = f match
    case Bam => "BAM"
    case Cram => "CRAM"
    case Fastq => "FASTQ"
    case Vcf => "VCF"
    case Gvcf => "GVCF"

/**
 * Source file entity for tracking user's BAM/CRAM files.
 *
 * The checksum is the stable identifier (path may change as user moves files).
 * Linked to alignments when analysis is performed.
 */
case class SourceFileEntity(
  id: UUID,
  alignmentId: Option[UUID],
  filePath: Option[String],
  fileChecksum: String,
  fileSize: Option[Long],
  fileFormat: Option[SourceFileFormat],
  lastVerifiedAt: Option[LocalDateTime],
  isAccessible: Boolean,
  hasBeenAnalyzed: Boolean,
  analysisCompletedAt: Option[LocalDateTime],
  createdAt: LocalDateTime,
  updatedAt: LocalDateTime
) extends Entity[UUID]

object SourceFileEntity:
  /**
   * Create a new source file entity.
   */
  def create(
    filePath: String,
    fileChecksum: String,
    fileSize: Option[Long] = None,
    fileFormat: Option[SourceFileFormat] = None
  ): SourceFileEntity =
    val now = LocalDateTime.now()
    SourceFileEntity(
      id = UUID.randomUUID(),
      alignmentId = None,
      filePath = Some(filePath),
      fileChecksum = fileChecksum,
      fileSize = fileSize,
      fileFormat = fileFormat,
      lastVerifiedAt = Some(now),
      isAccessible = true,
      hasBeenAnalyzed = false,
      analysisCompletedAt = None,
      createdAt = now,
      updatedAt = now
    )

/**
 * Repository for source file tracking.
 */
class SourceFileRepository:

  // ============================================
  // Core Operations
  // ============================================

  def findById(id: UUID)(using conn: Connection): Option[SourceFileEntity] =
    queryOne(
      "SELECT * FROM source_file WHERE id = ?",
      Seq(id)
    )(mapRow)

  def findAll()(using conn: Connection): List[SourceFileEntity] =
    queryList("SELECT * FROM source_file ORDER BY created_at DESC")(mapRow)

  def insert(entity: SourceFileEntity)(using conn: Connection): SourceFileEntity =
    executeUpdate(
      """INSERT INTO source_file (
        |  id, alignment_id, file_path, file_checksum, file_size, file_format,
        |  last_verified_at, is_accessible, has_been_analyzed, analysis_completed_at,
        |  created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.alignmentId,
        entity.filePath,
        entity.fileChecksum,
        entity.fileSize,
        entity.fileFormat.map(SourceFileFormat.toDbString),
        entity.lastVerifiedAt,
        entity.isAccessible,
        entity.hasBeenAnalyzed,
        entity.analysisCompletedAt,
        entity.createdAt,
        entity.updatedAt
      )
    )
    entity

  def update(entity: SourceFileEntity)(using conn: Connection): SourceFileEntity =
    val now = LocalDateTime.now()

    executeUpdate(
      """UPDATE source_file SET
        |  alignment_id = ?, file_path = ?, file_checksum = ?, file_size = ?,
        |  file_format = ?, last_verified_at = ?, is_accessible = ?,
        |  has_been_analyzed = ?, analysis_completed_at = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.alignmentId,
        entity.filePath,
        entity.fileChecksum,
        entity.fileSize,
        entity.fileFormat.map(SourceFileFormat.toDbString),
        entity.lastVerifiedAt,
        entity.isAccessible,
        entity.hasBeenAnalyzed,
        entity.analysisCompletedAt,
        now,
        entity.id
      )
    )
    entity.copy(updatedAt = now)

  def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM source_file WHERE id = ?", Seq(id)) > 0

  // ============================================
  // Query Operations
  // ============================================

  /**
   * Find source file by checksum (the stable identifier).
   */
  def findByChecksum(checksum: String)(using conn: Connection): Option[SourceFileEntity] =
    queryOne(
      "SELECT * FROM source_file WHERE file_checksum = ?",
      Seq(checksum)
    )(mapRow)

  /**
   * Find source file by path.
   */
  def findByPath(path: String)(using conn: Connection): Option[SourceFileEntity] =
    queryOne(
      "SELECT * FROM source_file WHERE file_path = ?",
      Seq(path)
    )(mapRow)

  /**
   * Find source files linked to an alignment.
   */
  def findByAlignment(alignmentId: UUID)(using conn: Connection): List[SourceFileEntity] =
    queryList(
      "SELECT * FROM source_file WHERE alignment_id = ?",
      Seq(alignmentId)
    )(mapRow)

  /**
   * Find all accessible source files.
   */
  def findAccessible()(using conn: Connection): List[SourceFileEntity] =
    queryList(
      "SELECT * FROM source_file WHERE is_accessible = TRUE ORDER BY updated_at DESC"
    )(mapRow)

  /**
   * Find all inaccessible source files (moved/deleted).
   */
  def findInaccessible()(using conn: Connection): List[SourceFileEntity] =
    queryList(
      "SELECT * FROM source_file WHERE is_accessible = FALSE ORDER BY updated_at DESC"
    )(mapRow)

  /**
   * Find source files that haven't been analyzed yet.
   */
  def findNotAnalyzed()(using conn: Connection): List[SourceFileEntity] =
    queryList(
      "SELECT * FROM source_file WHERE has_been_analyzed = FALSE AND is_accessible = TRUE ORDER BY created_at ASC"
    )(mapRow)

  /**
   * Check if a file with this checksum exists.
   */
  def existsByChecksum(checksum: String)(using conn: Connection): Boolean =
    queryOne(
      "SELECT 1 FROM source_file WHERE file_checksum = ?",
      Seq(checksum)
    ) { _ => true }.isDefined

  // ============================================
  // Update Operations
  // ============================================

  /**
   * Link a source file to an alignment.
   */
  def linkToAlignment(id: UUID, alignmentId: UUID)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE source_file SET alignment_id = ?, updated_at = ? WHERE id = ?",
      Seq(alignmentId, LocalDateTime.now(), id)
    ) > 0

  /**
   * Update file path (when user moves file).
   */
  def updatePath(id: UUID, newPath: String)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE source_file SET file_path = ?, last_verified_at = ?, is_accessible = TRUE, updated_at = ? WHERE id = ?",
      Seq(newPath, LocalDateTime.now(), LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark file as accessible after verification.
   */
  def markAccessible(id: UUID)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE source_file SET is_accessible = TRUE, last_verified_at = ?, updated_at = ? WHERE id = ?",
      Seq(LocalDateTime.now(), LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark file as inaccessible.
   */
  def markInaccessible(id: UUID)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE source_file SET is_accessible = FALSE, updated_at = ? WHERE id = ?",
      Seq(LocalDateTime.now(), id)
    ) > 0

  /**
   * Mark analysis as complete.
   */
  def markAnalyzed(id: UUID)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      "UPDATE source_file SET has_been_analyzed = TRUE, analysis_completed_at = ?, updated_at = ? WHERE id = ?",
      Seq(now, now, id)
    ) > 0

  /**
   * Register or update a source file by checksum.
   * Returns the existing or newly created entity.
   */
  def upsertByChecksum(
    filePath: String,
    fileChecksum: String,
    fileSize: Option[Long],
    fileFormat: Option[SourceFileFormat]
  )(using conn: Connection): SourceFileEntity =
    findByChecksum(fileChecksum) match
      case Some(existing) =>
        // Update path if changed
        if existing.filePath != Some(filePath) then
          updatePath(existing.id, filePath)
          existing.copy(filePath = Some(filePath), lastVerifiedAt = Some(LocalDateTime.now()))
        else
          markAccessible(existing.id)
          existing.copy(isAccessible = true, lastVerifiedAt = Some(LocalDateTime.now()))
      case None =>
        val entity = SourceFileEntity.create(filePath, fileChecksum, fileSize, fileFormat)
        insert(entity)

  // ============================================
  // Statistics
  // ============================================

  /**
   * Count source files by accessibility status.
   */
  def countByAccessibility()(using conn: Connection): (Long, Long) =
    val accessible = queryOne(
      "SELECT COUNT(*) FROM source_file WHERE is_accessible = TRUE"
    )(rs => rs.getLong(1)).getOrElse(0L)
    val inaccessible = queryOne(
      "SELECT COUNT(*) FROM source_file WHERE is_accessible = FALSE"
    )(rs => rs.getLong(1)).getOrElse(0L)
    (accessible, inaccessible)

  /**
   * Get total size of tracked source files.
   */
  def totalFileSize()(using conn: Connection): Long =
    queryOne(
      "SELECT COALESCE(SUM(file_size), 0) FROM source_file WHERE is_accessible = TRUE"
    )(rs => rs.getLong(1)).getOrElse(0L)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): SourceFileEntity =
    val formatStr = getOptString(rs, "file_format")
    val format = formatStr.map(SourceFileFormat.fromString)

    SourceFileEntity(
      id = getUUID(rs, "id"),
      alignmentId = Option(rs.getObject("alignment_id", classOf[UUID])),
      filePath = getOptString(rs, "file_path"),
      fileChecksum = rs.getString("file_checksum"),
      fileSize = getOptLong(rs, "file_size"),
      fileFormat = format,
      lastVerifiedAt = getOptDateTime(rs, "last_verified_at"),
      isAccessible = rs.getBoolean("is_accessible"),
      hasBeenAnalyzed = rs.getBoolean("has_been_analyzed"),
      analysisCompletedAt = getOptDateTime(rs, "analysis_completed_at"),
      createdAt = getDateTime(rs, "created_at"),
      updatedAt = getDateTime(rs, "updated_at")
    )
