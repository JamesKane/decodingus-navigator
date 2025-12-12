package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.yprofile.model.*
import com.decodingus.yprofile.model.YProfileCodecs.given
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Repository for Y profile variant persistence operations.
 */
class YProfileVariantRepository extends Repository[YProfileVariantEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YProfileVariantEntity] =
    queryOne(
      "SELECT * FROM y_profile_variant WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YProfileVariantEntity] =
    queryList("SELECT * FROM y_profile_variant ORDER BY position")(mapRow)

  override def insert(entity: YProfileVariantEntity)(using conn: Connection): YProfileVariantEntity =
    val strMetadataJson = entity.strMetadata.map(m => JsonValue(m.asJson.noSpaces))

    executeUpdate(
      """INSERT INTO y_profile_variant (
        |  id, y_profile_id, contig, position, end_position, ref_allele, alt_allele,
        |  variant_type, variant_name, rs_id, marker_name, repeat_count, str_metadata,
        |  consensus_allele, consensus_state, status,
        |  source_count, concordant_count, discordant_count, confidence_score,
        |  max_read_depth, max_quality_score, defining_haplogroup, haplogroup_branch_depth,
        |  last_updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.yProfileId,
        entity.contig,
        entity.position,
        entity.endPosition,
        entity.refAllele,
        entity.altAllele,
        entity.variantType.toString,
        entity.variantName,
        entity.rsId,
        entity.markerName,
        entity.repeatCount,
        strMetadataJson,
        entity.consensusAllele,
        entity.consensusState.toString,
        entity.status.toString,
        entity.sourceCount,
        entity.concordantCount,
        entity.discordantCount,
        entity.confidenceScore,
        entity.maxReadDepth,
        entity.maxQualityScore,
        entity.definingHaplogroup,
        entity.haplogroupBranchDepth,
        entity.lastUpdatedAt
      )
    )
    entity

  override def update(entity: YProfileVariantEntity)(using conn: Connection): YProfileVariantEntity =
    val strMetadataJson = entity.strMetadata.map(m => JsonValue(m.asJson.noSpaces))

    executeUpdate(
      """UPDATE y_profile_variant SET
        |  y_profile_id = ?, contig = ?, position = ?, end_position = ?, ref_allele = ?, alt_allele = ?,
        |  variant_type = ?, variant_name = ?, rs_id = ?, marker_name = ?, repeat_count = ?, str_metadata = ?,
        |  consensus_allele = ?, consensus_state = ?, status = ?,
        |  source_count = ?, concordant_count = ?, discordant_count = ?, confidence_score = ?,
        |  max_read_depth = ?, max_quality_score = ?, defining_haplogroup = ?, haplogroup_branch_depth = ?,
        |  last_updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.yProfileId,
        entity.contig,
        entity.position,
        entity.endPosition,
        entity.refAllele,
        entity.altAllele,
        entity.variantType.toString,
        entity.variantName,
        entity.rsId,
        entity.markerName,
        entity.repeatCount,
        strMetadataJson,
        entity.consensusAllele,
        entity.consensusState.toString,
        entity.status.toString,
        entity.sourceCount,
        entity.concordantCount,
        entity.discordantCount,
        entity.confidenceScore,
        entity.maxReadDepth,
        entity.maxQualityScore,
        entity.definingHaplogroup,
        entity.haplogroupBranchDepth,
        LocalDateTime.now(),
        entity.id
      )
    )
    entity.copy(lastUpdatedAt = LocalDateTime.now())

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM y_profile_variant WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_profile_variant") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_profile_variant WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Variant-Specific Queries
  // ============================================

  /**
   * Find all variants for a Y chromosome profile.
   */
  def findByProfile(yProfileId: UUID)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE y_profile_id = ? ORDER BY position",
      Seq(yProfileId)
    )(mapRow)

  /**
   * Find variants by type.
   */
  def findByType(yProfileId: UUID, variantType: YVariantType)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE y_profile_id = ? AND variant_type = ? ORDER BY position",
      Seq(yProfileId, variantType.toString)
    )(mapRow)

  /**
   * Find variants by status.
   */
  def findByStatus(yProfileId: UUID, status: YVariantStatus)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE y_profile_id = ? AND status = ? ORDER BY position",
      Seq(yProfileId, status.toString)
    )(mapRow)

  /**
   * Find variants by consensus state.
   */
  def findByConsensusState(yProfileId: UUID, state: YConsensusState)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE y_profile_id = ? AND consensus_state = ? ORDER BY position",
      Seq(yProfileId, state.toString)
    )(mapRow)

  /**
   * Find variant at a specific position.
   */
  def findByPosition(yProfileId: UUID, position: Long)(using conn: Connection): Option[YProfileVariantEntity] =
    queryOne(
      "SELECT * FROM y_profile_variant WHERE y_profile_id = ? AND position = ?",
      Seq(yProfileId, position)
    )(mapRow)

  /**
   * Find variant by name (marker name like M269, L21).
   */
  def findByVariantName(name: String)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE variant_name = ? ORDER BY position",
      Seq(name)
    )(mapRow)

  /**
   * Find STR variant by marker name (like DYS393).
   */
  def findByMarkerName(markerName: String)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE marker_name = ? ORDER BY position",
      Seq(markerName)
    )(mapRow)

  /**
   * Find variants defining a haplogroup.
   */
  def findByDefiningHaplogroup(haplogroup: String)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      "SELECT * FROM y_profile_variant WHERE defining_haplogroup = ? ORDER BY haplogroup_branch_depth, position",
      Seq(haplogroup)
    )(mapRow)

  /**
   * Find derived variants (positive for haplogroup determination).
   */
  def findDerivedVariants(yProfileId: UUID)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      """SELECT * FROM y_profile_variant
        |WHERE y_profile_id = ? AND consensus_state = 'DERIVED'
        |ORDER BY position
      """.stripMargin,
      Seq(yProfileId)
    )(mapRow)

  /**
   * Find novel/private variants.
   */
  def findNovelVariants(yProfileId: UUID)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      """SELECT * FROM y_profile_variant
        |WHERE y_profile_id = ? AND status = 'NOVEL'
        |ORDER BY position
      """.stripMargin,
      Seq(yProfileId)
    )(mapRow)

  /**
   * Find variants with conflicts.
   */
  def findConflictingVariants(yProfileId: UUID)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      """SELECT * FROM y_profile_variant
        |WHERE y_profile_id = ? AND status = 'CONFLICT'
        |ORDER BY discordant_count DESC, position
      """.stripMargin,
      Seq(yProfileId)
    )(mapRow)

  /**
   * Find STR markers.
   */
  def findStrMarkers(yProfileId: UUID)(using conn: Connection): List[YProfileVariantEntity] =
    queryList(
      """SELECT * FROM y_profile_variant
        |WHERE y_profile_id = ? AND variant_type = 'STR'
        |ORDER BY marker_name, position
      """.stripMargin,
      Seq(yProfileId)
    )(mapRow)

  /**
   * Delete all variants for a profile.
   */
  def deleteByProfile(yProfileId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_profile_variant WHERE y_profile_id = ?", Seq(yProfileId))

  /**
   * Count variants by status for a profile.
   */
  def countByStatus(yProfileId: UUID)(using conn: Connection): Map[YVariantStatus, Long] =
    queryList(
      """SELECT status, COUNT(*) as cnt FROM y_profile_variant
        |WHERE y_profile_id = ?
        |GROUP BY status
      """.stripMargin,
      Seq(yProfileId)
    ) { rs =>
      (YVariantStatus.fromString(rs.getString("status")), rs.getLong("cnt"))
    }.toMap

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YProfileVariantEntity =
    val strMetadataJson = getOptJsonString(rs, "str_metadata")
    val strMetadata = strMetadataJson.flatMap { json =>
      parse(json).flatMap(_.as[StrMetadata]).toOption
    }

    YProfileVariantEntity(
      id = getUUID(rs, "id"),
      yProfileId = getUUID(rs, "y_profile_id"),
      contig = rs.getString("contig"),
      position = rs.getLong("position"),
      endPosition = getOptLong(rs, "end_position"),
      refAllele = rs.getString("ref_allele"),
      altAllele = rs.getString("alt_allele"),
      variantType = YVariantType.fromString(rs.getString("variant_type")),
      variantName = getOptString(rs, "variant_name"),
      rsId = getOptString(rs, "rs_id"),
      markerName = getOptString(rs, "marker_name"),
      repeatCount = getOptInt(rs, "repeat_count"),
      strMetadata = strMetadata,
      consensusAllele = getOptString(rs, "consensus_allele"),
      consensusState = YConsensusState.fromString(rs.getString("consensus_state")),
      status = YVariantStatus.fromString(rs.getString("status")),
      sourceCount = rs.getInt("source_count"),
      concordantCount = rs.getInt("concordant_count"),
      discordantCount = rs.getInt("discordant_count"),
      confidenceScore = rs.getDouble("confidence_score"),
      maxReadDepth = getOptInt(rs, "max_read_depth"),
      maxQualityScore = getOptDouble(rs, "max_quality_score"),
      definingHaplogroup = getOptString(rs, "defining_haplogroup"),
      haplogroupBranchDepth = getOptInt(rs, "haplogroup_branch_depth"),
      lastUpdatedAt = getDateTime(rs, "last_updated_at")
    )
