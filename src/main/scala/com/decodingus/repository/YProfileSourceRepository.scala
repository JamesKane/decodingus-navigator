package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.yprofile.model.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Repository for Y profile source persistence operations.
 */
class YProfileSourceRepository extends Repository[YProfileSourceEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YProfileSourceEntity] =
    queryOne(
      "SELECT * FROM y_profile_source WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YProfileSourceEntity] =
    queryList("SELECT * FROM y_profile_source ORDER BY imported_at DESC")(mapRow)

  override def insert(entity: YProfileSourceEntity)(using conn: Connection): YProfileSourceEntity =
    executeUpdate(
      """INSERT INTO y_profile_source (
        |  id, y_profile_id, source_type, source_ref, vendor, test_name, test_date,
        |  method_tier, mean_read_depth, mean_mapping_quality, coverage_pct,
        |  variant_count, str_marker_count, novel_variant_count,
        |  alignment_id, reference_build, imported_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.yProfileId,
        entity.sourceType.toString,
        entity.sourceRef,
        entity.vendor,
        entity.testName,
        entity.testDate,
        entity.methodTier,
        entity.meanReadDepth,
        entity.meanMappingQuality,
        entity.coveragePct,
        entity.variantCount,
        entity.strMarkerCount,
        entity.novelVariantCount,
        entity.alignmentId,
        entity.referenceBuild,
        entity.importedAt
      )
    )
    entity

  override def update(entity: YProfileSourceEntity)(using conn: Connection): YProfileSourceEntity =
    executeUpdate(
      """UPDATE y_profile_source SET
        |  y_profile_id = ?, source_type = ?, source_ref = ?, vendor = ?, test_name = ?, test_date = ?,
        |  method_tier = ?, mean_read_depth = ?, mean_mapping_quality = ?, coverage_pct = ?,
        |  variant_count = ?, str_marker_count = ?, novel_variant_count = ?,
        |  alignment_id = ?, reference_build = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.yProfileId,
        entity.sourceType.toString,
        entity.sourceRef,
        entity.vendor,
        entity.testName,
        entity.testDate,
        entity.methodTier,
        entity.meanReadDepth,
        entity.meanMappingQuality,
        entity.coveragePct,
        entity.variantCount,
        entity.strMarkerCount,
        entity.novelVariantCount,
        entity.alignmentId,
        entity.referenceBuild,
        entity.id
      )
    )
    entity

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM y_profile_source WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_profile_source") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_profile_source WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Source-Specific Queries
  // ============================================

  /**
   * Find all sources for a Y chromosome profile.
   */
  def findByProfile(yProfileId: UUID)(using conn: Connection): List[YProfileSourceEntity] =
    queryList(
      "SELECT * FROM y_profile_source WHERE y_profile_id = ? ORDER BY method_tier DESC, imported_at DESC",
      Seq(yProfileId)
    )(mapRow)

  /**
   * Find sources by type.
   */
  def findByType(sourceType: YProfileSourceType)(using conn: Connection): List[YProfileSourceEntity] =
    queryList(
      "SELECT * FROM y_profile_source WHERE source_type = ? ORDER BY imported_at DESC",
      Seq(sourceType.toString)
    )(mapRow)

  /**
   * Find sources by alignment.
   */
  def findByAlignment(alignmentId: UUID)(using conn: Connection): List[YProfileSourceEntity] =
    queryList(
      "SELECT * FROM y_profile_source WHERE alignment_id = ? ORDER BY imported_at DESC",
      Seq(alignmentId)
    )(mapRow)

  /**
   * Find sources by vendor.
   */
  def findByVendor(vendor: String)(using conn: Connection): List[YProfileSourceEntity] =
    queryList(
      "SELECT * FROM y_profile_source WHERE vendor = ? ORDER BY imported_at DESC",
      Seq(vendor)
    )(mapRow)

  /**
   * Delete all sources for a profile.
   */
  def deleteByProfile(yProfileId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_profile_source WHERE y_profile_id = ?", Seq(yProfileId))

  /**
   * Count sources for a profile.
   */
  def countByProfile(yProfileId: UUID)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM y_profile_source WHERE y_profile_id = ?",
      Seq(yProfileId)
    ) { rs => rs.getLong(1) }.getOrElse(0L)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YProfileSourceEntity =
    val testDateTs = rs.getTimestamp("test_date")
    val testDate = if rs.wasNull() then None else Some(testDateTs.toLocalDateTime)

    YProfileSourceEntity(
      id = getUUID(rs, "id"),
      yProfileId = getUUID(rs, "y_profile_id"),
      sourceType = YProfileSourceType.fromString(rs.getString("source_type")),
      sourceRef = getOptString(rs, "source_ref"),
      vendor = getOptString(rs, "vendor"),
      testName = getOptString(rs, "test_name"),
      testDate = testDate,
      methodTier = rs.getInt("method_tier"),
      meanReadDepth = getOptDouble(rs, "mean_read_depth"),
      meanMappingQuality = getOptDouble(rs, "mean_mapping_quality"),
      coveragePct = getOptDouble(rs, "coverage_pct"),
      variantCount = rs.getInt("variant_count"),
      strMarkerCount = rs.getInt("str_marker_count"),
      novelVariantCount = rs.getInt("novel_variant_count"),
      alignmentId = getOptUUID(rs, "alignment_id"),
      referenceBuild = getOptString(rs, "reference_build"),
      importedAt = getDateTime(rs, "imported_at")
    )
