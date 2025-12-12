package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.yprofile.model.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Repository for Y variant source call persistence operations.
 */
class YVariantSourceCallRepository extends Repository[YVariantSourceCallEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YVariantSourceCallEntity] =
    queryOne(
      "SELECT * FROM y_variant_source_call WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YVariantSourceCallEntity] =
    queryList("SELECT * FROM y_variant_source_call ORDER BY called_at DESC")(mapRow)

  override def insert(entity: YVariantSourceCallEntity)(using conn: Connection): YVariantSourceCallEntity =
    executeUpdate(
      """INSERT INTO y_variant_source_call (
        |  id, variant_id, source_id, called_allele, call_state, called_repeat_count,
        |  read_depth, quality_score, mapping_quality, variant_allele_frequency,
        |  callable_state, concordance_weight, called_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.variantId,
        entity.sourceId,
        entity.calledAllele,
        entity.callState.toString,
        entity.calledRepeatCount,
        entity.readDepth,
        entity.qualityScore,
        entity.mappingQuality,
        entity.variantAlleleFrequency,
        entity.callableState.map(_.toString),
        entity.concordanceWeight,
        entity.calledAt
      )
    )
    entity

  override def update(entity: YVariantSourceCallEntity)(using conn: Connection): YVariantSourceCallEntity =
    executeUpdate(
      """UPDATE y_variant_source_call SET
        |  variant_id = ?, source_id = ?, called_allele = ?, call_state = ?, called_repeat_count = ?,
        |  read_depth = ?, quality_score = ?, mapping_quality = ?, variant_allele_frequency = ?,
        |  callable_state = ?, concordance_weight = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.variantId,
        entity.sourceId,
        entity.calledAllele,
        entity.callState.toString,
        entity.calledRepeatCount,
        entity.readDepth,
        entity.qualityScore,
        entity.mappingQuality,
        entity.variantAlleleFrequency,
        entity.callableState.map(_.toString),
        entity.concordanceWeight,
        entity.id
      )
    )
    entity

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM y_variant_source_call WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_variant_source_call") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_variant_source_call WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Source Call-Specific Queries
  // ============================================

  /**
   * Find all source calls for a variant.
   */
  def findByVariant(variantId: UUID)(using conn: Connection): List[YVariantSourceCallEntity] =
    queryList(
      "SELECT * FROM y_variant_source_call WHERE variant_id = ? ORDER BY concordance_weight DESC",
      Seq(variantId)
    )(mapRow)

  /**
   * Find all source calls from a source.
   */
  def findBySource(sourceId: UUID)(using conn: Connection): List[YVariantSourceCallEntity] =
    queryList(
      "SELECT * FROM y_variant_source_call WHERE source_id = ? ORDER BY called_at DESC",
      Seq(sourceId)
    )(mapRow)

  /**
   * Find source call for a specific variant from a specific source.
   */
  def findByVariantAndSource(variantId: UUID, sourceId: UUID)(using conn: Connection): Option[YVariantSourceCallEntity] =
    queryOne(
      "SELECT * FROM y_variant_source_call WHERE variant_id = ? AND source_id = ?",
      Seq(variantId, sourceId)
    )(mapRow)

  /**
   * Find calls by call state.
   */
  def findByCallState(variantId: UUID, callState: YConsensusState)(using conn: Connection): List[YVariantSourceCallEntity] =
    queryList(
      "SELECT * FROM y_variant_source_call WHERE variant_id = ? AND call_state = ?",
      Seq(variantId, callState.toString)
    )(mapRow)

  /**
   * Delete all source calls for a variant.
   */
  def deleteByVariant(variantId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_variant_source_call WHERE variant_id = ?", Seq(variantId))

  /**
   * Delete all source calls from a source.
   */
  def deleteBySource(sourceId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_variant_source_call WHERE source_id = ?", Seq(sourceId))

  /**
   * Count calls for a variant.
   */
  def countByVariant(variantId: UUID)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM y_variant_source_call WHERE variant_id = ?",
      Seq(variantId)
    ) { rs => rs.getLong(1) }.getOrElse(0L)

  /**
   * Count concordant calls (matching consensus) for a variant.
   */
  def countConcordantCalls(variantId: UUID, consensusAllele: String)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM y_variant_source_call WHERE variant_id = ? AND called_allele = ?",
      Seq(variantId, consensusAllele)
    ) { rs => rs.getLong(1) }.getOrElse(0L)

  /**
   * Sum of concordance weights for calls matching an allele.
   */
  def sumWeightsForAllele(variantId: UUID, allele: String)(using conn: Connection): Double =
    queryOne(
      "SELECT COALESCE(SUM(concordance_weight), 0) FROM y_variant_source_call WHERE variant_id = ? AND called_allele = ?",
      Seq(variantId, allele)
    ) { rs => rs.getDouble(1) }.getOrElse(0.0)

  /**
   * Batch insert source calls.
   */
  def insertBatch(calls: List[YVariantSourceCallEntity])(using conn: Connection): List[YVariantSourceCallEntity] =
    calls.foreach(insert)
    calls

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YVariantSourceCallEntity =
    val callableStateStr = getOptString(rs, "callable_state")
    val callableState = callableStateStr.map(YCallableState.fromString)

    YVariantSourceCallEntity(
      id = getUUID(rs, "id"),
      variantId = getUUID(rs, "variant_id"),
      sourceId = getUUID(rs, "source_id"),
      calledAllele = rs.getString("called_allele"),
      callState = YConsensusState.fromString(rs.getString("call_state")),
      calledRepeatCount = getOptInt(rs, "called_repeat_count"),
      readDepth = getOptInt(rs, "read_depth"),
      qualityScore = getOptDouble(rs, "quality_score"),
      mappingQuality = getOptDouble(rs, "mapping_quality"),
      variantAlleleFrequency = getOptDouble(rs, "variant_allele_frequency"),
      callableState = callableState,
      concordanceWeight = rs.getDouble("concordance_weight"),
      calledAt = getDateTime(rs, "called_at")
    )
