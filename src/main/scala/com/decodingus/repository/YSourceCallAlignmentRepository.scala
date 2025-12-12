package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.yprofile.model.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Repository for Y source call alignment persistence operations.
 *
 * Alignments represent coordinate representations of a source call in different
 * reference builds. Multiple alignments of the same source call (e.g., GRCh37,
 * GRCh38, hs1) are ONE piece of evidence, not multiple.
 */
class YSourceCallAlignmentRepository extends Repository[YSourceCallAlignmentEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YSourceCallAlignmentEntity] =
    queryOne(
      "SELECT * FROM y_source_call_alignment WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YSourceCallAlignmentEntity] =
    queryList("SELECT * FROM y_source_call_alignment ORDER BY created_at")(mapRow)

  override def insert(entity: YSourceCallAlignmentEntity)(using conn: Connection): YSourceCallAlignmentEntity =
    executeUpdate(
      """INSERT INTO y_source_call_alignment (
        |  id, source_call_id, reference_build, contig, position,
        |  ref_allele, alt_allele, called_allele,
        |  read_depth, mapping_quality, base_quality, variant_allele_frequency,
        |  graph_node, graph_offset, alignment_id, created_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.sourceCallId,
        entity.referenceBuild,
        entity.contig,
        entity.position,
        entity.refAllele,
        entity.altAllele,
        entity.calledAllele,
        entity.readDepth,
        entity.mappingQuality,
        entity.baseQuality,
        entity.variantAlleleFrequency,
        entity.graphNode,
        entity.graphOffset,
        entity.alignmentId,
        entity.createdAt
      )
    )
    entity

  override def update(entity: YSourceCallAlignmentEntity)(using conn: Connection): YSourceCallAlignmentEntity =
    executeUpdate(
      """UPDATE y_source_call_alignment SET
        |  source_call_id = ?, reference_build = ?, contig = ?, position = ?,
        |  ref_allele = ?, alt_allele = ?, called_allele = ?,
        |  read_depth = ?, mapping_quality = ?, base_quality = ?, variant_allele_frequency = ?,
        |  graph_node = ?, graph_offset = ?, alignment_id = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.sourceCallId,
        entity.referenceBuild,
        entity.contig,
        entity.position,
        entity.refAllele,
        entity.altAllele,
        entity.calledAllele,
        entity.readDepth,
        entity.mappingQuality,
        entity.baseQuality,
        entity.variantAlleleFrequency,
        entity.graphNode,
        entity.graphOffset,
        entity.alignmentId,
        entity.id
      )
    )
    entity

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM y_source_call_alignment WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_source_call_alignment") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_source_call_alignment WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Alignment-Specific Queries
  // ============================================

  /**
   * Find all alignments for a source call.
   */
  def findBySourceCall(sourceCallId: UUID)(using conn: Connection): List[YSourceCallAlignmentEntity] =
    queryList(
      "SELECT * FROM y_source_call_alignment WHERE source_call_id = ? ORDER BY reference_build",
      Seq(sourceCallId)
    )(mapRow)

  /**
   * Find alignment for a specific source call and reference build.
   */
  def findBySourceCallAndBuild(sourceCallId: UUID, referenceBuild: String)(using conn: Connection): Option[YSourceCallAlignmentEntity] =
    queryOne(
      "SELECT * FROM y_source_call_alignment WHERE source_call_id = ? AND reference_build = ?",
      Seq(sourceCallId, referenceBuild)
    )(mapRow)

  /**
   * Find all alignments for a reference build.
   */
  def findByReferenceBuild(referenceBuild: String)(using conn: Connection): List[YSourceCallAlignmentEntity] =
    queryList(
      "SELECT * FROM y_source_call_alignment WHERE reference_build = ? ORDER BY position",
      Seq(referenceBuild)
    )(mapRow)

  /**
   * Find alignments by position in a specific reference build.
   */
  def findByPosition(referenceBuild: String, contig: String, position: Long)(using conn: Connection): List[YSourceCallAlignmentEntity] =
    queryList(
      """SELECT * FROM y_source_call_alignment
        |WHERE reference_build = ? AND contig = ? AND position = ?
      """.stripMargin,
      Seq(referenceBuild, contig, position)
    )(mapRow)

  /**
   * Find alignments in a position range.
   */
  def findByPositionRange(
    referenceBuild: String,
    contig: String,
    startPosition: Long,
    endPosition: Long
  )(using conn: Connection): List[YSourceCallAlignmentEntity] =
    queryList(
      """SELECT * FROM y_source_call_alignment
        |WHERE reference_build = ? AND contig = ? AND position >= ? AND position <= ?
        |ORDER BY position
      """.stripMargin,
      Seq(referenceBuild, contig, startPosition, endPosition)
    )(mapRow)

  /**
   * Delete all alignments for a source call.
   */
  def deleteBySourceCall(sourceCallId: UUID)(using conn: Connection): Int =
    executeUpdate(
      "DELETE FROM y_source_call_alignment WHERE source_call_id = ?",
      Seq(sourceCallId)
    )

  /**
   * Count alignments by reference build.
   */
  def countByReferenceBuild(referenceBuild: String)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM y_source_call_alignment WHERE reference_build = ?",
      Seq(referenceBuild)
    ) { rs => rs.getLong(1) }.getOrElse(0L)

  /**
   * Get distinct reference builds in use.
   */
  def getDistinctReferenceBuilds()(using conn: Connection): List[String] =
    queryList(
      "SELECT DISTINCT reference_build FROM y_source_call_alignment ORDER BY reference_build"
    ) { rs => rs.getString(1) }

  /**
   * Upsert alignment (insert or update based on source_call_id + reference_build).
   */
  def upsert(entity: YSourceCallAlignmentEntity)(using conn: Connection): YSourceCallAlignmentEntity =
    findBySourceCallAndBuild(entity.sourceCallId, entity.referenceBuild) match
      case Some(existing) =>
        update(entity.copy(id = existing.id, createdAt = existing.createdAt))
      case None =>
        insert(entity)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YSourceCallAlignmentEntity =
    YSourceCallAlignmentEntity(
      id = getUUID(rs, "id"),
      sourceCallId = getUUID(rs, "source_call_id"),
      referenceBuild = rs.getString("reference_build"),
      contig = rs.getString("contig"),
      position = rs.getLong("position"),
      refAllele = rs.getString("ref_allele"),
      altAllele = rs.getString("alt_allele"),
      calledAllele = rs.getString("called_allele"),
      readDepth = getOptInt(rs, "read_depth"),
      mappingQuality = getOptDouble(rs, "mapping_quality"),
      baseQuality = getOptDouble(rs, "base_quality"),
      variantAlleleFrequency = getOptDouble(rs, "variant_allele_frequency"),
      graphNode = getOptString(rs, "graph_node"),
      graphOffset = getOptInt(rs, "graph_offset"),
      alignmentId = getOptUUID(rs, "alignment_id"),
      createdAt = getDateTime(rs, "created_at")
    )
