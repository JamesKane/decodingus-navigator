package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.yprofile.model.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Repository for Y variant audit trail persistence operations.
 */
class YVariantAuditRepository extends Repository[YVariantAuditEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YVariantAuditEntity] =
    queryOne(
      "SELECT * FROM y_variant_audit WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YVariantAuditEntity] =
    queryList("SELECT * FROM y_variant_audit ORDER BY created_at DESC")(mapRow)

  override def insert(entity: YVariantAuditEntity)(using conn: Connection): YVariantAuditEntity =
    executeUpdate(
      """INSERT INTO y_variant_audit (
        |  id, variant_id, action,
        |  previous_consensus_allele, previous_consensus_state, previous_status, previous_confidence,
        |  new_consensus_allele, new_consensus_state, new_status, new_confidence,
        |  user_id, reason, supporting_evidence, created_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.variantId,
        entity.action.toString,
        entity.previousConsensusAllele,
        entity.previousConsensusState.map(_.toString),
        entity.previousStatus.map(_.toString),
        entity.previousConfidence,
        entity.newConsensusAllele,
        entity.newConsensusState.map(_.toString),
        entity.newStatus.map(_.toString),
        entity.newConfidence,
        entity.userId,
        entity.reason,
        entity.supportingEvidence,
        entity.createdAt
      )
    )
    entity

  override def update(entity: YVariantAuditEntity)(using conn: Connection): YVariantAuditEntity =
    // Audit entries are immutable by design - only update is allowed for supporting evidence
    executeUpdate(
      """UPDATE y_variant_audit SET
        |  supporting_evidence = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.supportingEvidence,
        entity.id
      )
    )
    entity

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM y_variant_audit WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_variant_audit") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_variant_audit WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Audit-Specific Queries
  // ============================================

  /**
   * Find all audit entries for a variant.
   */
  def findByVariant(variantId: UUID)(using conn: Connection): List[YVariantAuditEntity] =
    queryList(
      "SELECT * FROM y_variant_audit WHERE variant_id = ? ORDER BY created_at DESC",
      Seq(variantId)
    )(mapRow)

  /**
   * Find audit entries by action type.
   */
  def findByAction(action: YAuditAction)(using conn: Connection): List[YVariantAuditEntity] =
    queryList(
      "SELECT * FROM y_variant_audit WHERE action = ? ORDER BY created_at DESC",
      Seq(action.toString)
    )(mapRow)

  /**
   * Find audit entries by user.
   */
  def findByUser(userId: String)(using conn: Connection): List[YVariantAuditEntity] =
    queryList(
      "SELECT * FROM y_variant_audit WHERE user_id = ? ORDER BY created_at DESC",
      Seq(userId)
    )(mapRow)

  /**
   * Find recent audit entries.
   */
  def findRecent(limit: Int)(using conn: Connection): List[YVariantAuditEntity] =
    queryList(
      s"SELECT * FROM y_variant_audit ORDER BY created_at DESC LIMIT $limit"
    )(mapRow)

  /**
   * Find audit entries in a date range.
   */
  def findInDateRange(startDate: LocalDateTime, endDate: LocalDateTime)(using conn: Connection): List[YVariantAuditEntity] =
    queryList(
      "SELECT * FROM y_variant_audit WHERE created_at >= ? AND created_at <= ? ORDER BY created_at DESC",
      Seq(startDate, endDate)
    )(mapRow)

  /**
   * Find overrides for a variant.
   */
  def findOverrides(variantId: UUID)(using conn: Connection): List[YVariantAuditEntity] =
    queryList(
      "SELECT * FROM y_variant_audit WHERE variant_id = ? AND action = 'OVERRIDE' ORDER BY created_at DESC",
      Seq(variantId)
    )(mapRow)

  /**
   * Check if a variant has been manually overridden.
   */
  def hasOverride(variantId: UUID)(using conn: Connection): Boolean =
    queryOne(
      "SELECT 1 FROM y_variant_audit WHERE variant_id = ? AND action = 'OVERRIDE' LIMIT 1",
      Seq(variantId)
    ) { _ => true }.isDefined

  /**
   * Get the most recent audit entry for a variant.
   */
  def findMostRecent(variantId: UUID)(using conn: Connection): Option[YVariantAuditEntity] =
    queryOne(
      "SELECT * FROM y_variant_audit WHERE variant_id = ? ORDER BY created_at DESC LIMIT 1",
      Seq(variantId)
    )(mapRow)

  /**
   * Delete all audit entries for a variant.
   */
  def deleteByVariant(variantId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_variant_audit WHERE variant_id = ?", Seq(variantId))

  /**
   * Count audit entries for a variant.
   */
  def countByVariant(variantId: UUID)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM y_variant_audit WHERE variant_id = ?",
      Seq(variantId)
    ) { rs => rs.getLong(1) }.getOrElse(0L)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YVariantAuditEntity =
    val previousConsensusStateStr = getOptString(rs, "previous_consensus_state")
    val previousConsensusState = previousConsensusStateStr.map(YConsensusState.fromString)

    val previousStatusStr = getOptString(rs, "previous_status")
    val previousStatus = previousStatusStr.map(YVariantStatus.fromString)

    val newConsensusStateStr = getOptString(rs, "new_consensus_state")
    val newConsensusState = newConsensusStateStr.map(YConsensusState.fromString)

    val newStatusStr = getOptString(rs, "new_status")
    val newStatus = newStatusStr.map(YVariantStatus.fromString)

    YVariantAuditEntity(
      id = getUUID(rs, "id"),
      variantId = getUUID(rs, "variant_id"),
      action = YAuditAction.fromString(rs.getString("action")),
      previousConsensusAllele = getOptString(rs, "previous_consensus_allele"),
      previousConsensusState = previousConsensusState,
      previousStatus = previousStatus,
      previousConfidence = getOptDouble(rs, "previous_confidence"),
      newConsensusAllele = getOptString(rs, "new_consensus_allele"),
      newConsensusState = newConsensusState,
      newStatus = newStatus,
      newConfidence = getOptDouble(rs, "new_confidence"),
      userId = getOptString(rs, "user_id"),
      reason = rs.getString("reason"),
      supportingEvidence = getOptString(rs, "supporting_evidence"),
      createdAt = getDateTime(rs, "created_at")
    )
