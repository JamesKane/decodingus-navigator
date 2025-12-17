package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.{AlignmentMetrics, ContigMetrics, FileInfo}
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Alignment entity for database persistence.
 *
 * Represents an alignment of sequence data to a reference genome,
 * with associated QC metrics and file references.
 */
case class AlignmentEntity(
                            id: UUID,
                            sequenceRunId: UUID,
                            referenceBuild: String,
                            aligner: String,
                            variantCaller: Option[String],
                            metrics: Option[AlignmentMetrics],
                            files: List[FileInfo],
                            meta: EntityMeta
                          ) extends Entity[UUID]

object AlignmentEntity:
  // Circe codecs for JSON serialization
  given Encoder[FileInfo] = deriveEncoder

  given Decoder[FileInfo] = deriveDecoder

  given Encoder[ContigMetrics] = deriveEncoder

  given Decoder[ContigMetrics] = deriveDecoder

  given Encoder[AlignmentMetrics] = deriveEncoder

  given Decoder[AlignmentMetrics] = deriveDecoder

  /**
   * Create a new AlignmentEntity with generated ID and initial metadata.
   */
  def create(
              sequenceRunId: UUID,
              referenceBuild: String,
              aligner: String,
              variantCaller: Option[String] = None
            ): AlignmentEntity = AlignmentEntity(
    id = UUID.randomUUID(),
    sequenceRunId = sequenceRunId,
    referenceBuild = referenceBuild,
    aligner = aligner,
    variantCaller = variantCaller,
    metrics = None,
    files = List.empty,
    meta = EntityMeta.create()
  )

/**
 * Repository for alignment persistence operations.
 */
class AlignmentRepository extends SyncableRepositoryBase[AlignmentEntity]:

  import AlignmentEntity.given

  override protected def tableName: String = "alignment"

  /**
   * Map internal reference build names to database constraint values.
   * Internal: CHM13v2 -> DB: T2T-CHM13
   */
  private def toDbReferenceBuild(build: String): String = build match {
    case "CHM13v2" => "T2T-CHM13"
    case other => other
  }

  /**
   * Map database reference build names back to internal names.
   * DB: T2T-CHM13 -> Internal: CHM13v2
   */
  private def fromDbReferenceBuild(build: String): String = build match {
    case "T2T-CHM13" => "CHM13v2"
    case other => other
  }

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[AlignmentEntity] =
    queryOne(
      "SELECT * FROM alignment WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[AlignmentEntity] =
    queryList("SELECT * FROM alignment ORDER BY created_at DESC")(mapRow)

  override def insert(entity: AlignmentEntity)(using conn: Connection): AlignmentEntity =
    val metricsJson = entity.metrics.map(m => JsonValue(m.asJson.noSpaces))
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO alignment (
        |  id, sequence_run_id, reference_build, aligner, variant_caller,
        |  metrics, files, sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.sequenceRunId,
        toDbReferenceBuild(entity.referenceBuild),
        entity.aligner,
        entity.variantCaller,
        metricsJson,
        filesJson,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: AlignmentEntity)(using conn: Connection): AlignmentEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val metricsJson = entity.metrics.map(m => JsonValue(m.asJson.noSpaces))
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """UPDATE alignment SET
        |  sequence_run_id = ?, reference_build = ?, aligner = ?, variant_caller = ?,
        |  metrics = ?, files = ?, sync_status = ?, at_uri = ?, at_cid = ?,
        |  version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.sequenceRunId,
        toDbReferenceBuild(entity.referenceBuild),
        entity.aligner,
        entity.variantCaller,
        metricsJson,
        filesJson,
        updatedMeta.syncStatus,
        updatedMeta.atUri,
        updatedMeta.atCid,
        updatedMeta.version,
        updatedMeta.updatedAt,
        entity.id
      )
    )
    entity.copy(meta = updatedMeta)

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM alignment WHERE id = ?", Seq(id)) > 0

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM alignment WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Alignment-Specific Queries
  // ============================================

  /**
   * Find all alignments for a sequence run.
   */
  def findBySequenceRun(sequenceRunId: UUID)(using conn: Connection): List[AlignmentEntity] =
    queryList(
      "SELECT * FROM alignment WHERE sequence_run_id = ? ORDER BY created_at DESC",
      Seq(sequenceRunId)
    )(mapRow)

  /**
   * Find alignments by reference build.
   */
  def findByReferenceBuild(referenceBuild: String)(using conn: Connection): List[AlignmentEntity] =
    queryList(
      "SELECT * FROM alignment WHERE reference_build = ? ORDER BY created_at DESC",
      Seq(referenceBuild)
    )(mapRow)

  /**
   * Find alignment for a specific sequence run and reference build.
   * Returns the most recent if multiple exist.
   */
  def findBySequenceRunAndReference(
                                     sequenceRunId: UUID,
                                     referenceBuild: String
                                   )(using conn: Connection): Option[AlignmentEntity] =
    queryOne(
      """SELECT * FROM alignment
        |WHERE sequence_run_id = ? AND reference_build = ?
        |ORDER BY created_at DESC
        |LIMIT 1
      """.stripMargin,
      Seq(sequenceRunId, referenceBuild)
    )(mapRow)

  /**
   * Update metrics for an alignment.
   * Called after analysis completes.
   */
  def updateMetrics(id: UUID, metrics: AlignmentMetrics)(using conn: Connection): Boolean =
    val metricsJson = JsonValue(metrics.asJson.noSpaces)
    executeUpdate(
      """UPDATE alignment SET
        |  metrics = ?,
        |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
        |  updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(metricsJson, LocalDateTime.now(), id)
    ) > 0

  /**
   * Add a file to the alignment.
   */
  def addFile(id: UUID, file: FileInfo)(using conn: Connection): Boolean =
    findById(id) match
      case Some(entity) =>
        val updatedFiles = entity.files :+ file
        val filesJson = JsonValue(updatedFiles.asJson.noSpaces)
        executeUpdate(
          """UPDATE alignment SET
            |  files = ?,
            |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
            |  updated_at = ?
            |WHERE id = ?
          """.stripMargin,
          Seq(filesJson, LocalDateTime.now(), id)
        ) > 0
      case None => false

  /**
   * Update variant caller information.
   */
  def updateVariantCaller(id: UUID, variantCaller: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE alignment SET
        |  variant_caller = ?,
        |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
        |  updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(variantCaller, LocalDateTime.now(), id)
    ) > 0

  /**
   * Find alignments with VCF files (have vcfPath in metrics).
   */
  def findWithVcf()(using conn: Connection): List[AlignmentEntity] =
    // Using JSON path query for H2 with PostgreSQL compatibility
    queryList(
      """SELECT * FROM alignment
        |WHERE metrics IS NOT NULL
        |AND JSON_VALUE(metrics, '$.vcfPath') IS NOT NULL
        |ORDER BY created_at DESC
      """.stripMargin
    )(mapRow)

  /**
   * Count alignments per sequence run.
   */
  def countBySequenceRun(sequenceRunId: UUID)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM alignment WHERE sequence_run_id = ?",
      Seq(sequenceRunId)
    ) { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  /**
   * Get all alignments for a biosample via join through sequence_run.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[AlignmentEntity] =
    queryList(
      """SELECT a.* FROM alignment a
        |JOIN sequence_run sr ON a.sequence_run_id = sr.id
        |WHERE sr.biosample_id = ?
        |ORDER BY a.created_at DESC
      """.stripMargin,
      Seq(biosampleId)
    )(mapRow)

  // ============================================
  // Result Set Mapping
  // ============================================

  override protected def mapRow(rs: ResultSet): AlignmentEntity =
    val metricsJson = getOptJsonString(rs, "metrics")
    val metrics = metricsJson.flatMap { json =>
      parse(json).flatMap(_.as[AlignmentMetrics]).toOption
    }

    val filesJson = getOptJsonString(rs, "files").getOrElse("[]")
    val files = parse(filesJson).flatMap(_.as[List[FileInfo]]).getOrElse(List.empty)

    AlignmentEntity(
      id = getUUID(rs, "id"),
      sequenceRunId = getUUID(rs, "sequence_run_id"),
      referenceBuild = fromDbReferenceBuild(rs.getString("reference_build")),
      aligner = rs.getString("aligner"),
      variantCaller = getOptString(rs, "variant_caller"),
      metrics = metrics,
      files = files,
      meta = readEntityMeta(rs)
    )
