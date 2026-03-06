package com.decodingus.yprofile.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.repository.{Repository, SqlHelpers}
import com.decodingus.yprofile.model.*

import java.sql.{Connection, ResultSet}
import java.util.UUID

/**
 * Repository for Y profile callable region persistence operations.
 */
class YProfileRegionRepository extends Repository[YProfileRegionEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YProfileRegionEntity] =
    queryOne(
      "SELECT * FROM y_profile_region WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YProfileRegionEntity] =
    queryList("SELECT * FROM y_profile_region ORDER BY start_position")(mapRow)

  override def insert(entity: YProfileRegionEntity)(using conn: Connection): YProfileRegionEntity =
    executeUpdate(
      """INSERT INTO y_profile_region (
        |  id, y_profile_id, source_id, contig, start_position, end_position,
        |  callable_state, mean_coverage, mean_mapping_quality, callable_loci_cache_id
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.yProfileId,
        entity.sourceId,
        entity.contig,
        entity.startPosition,
        entity.endPosition,
        entity.callableState.toString,
        entity.meanCoverage,
        entity.meanMappingQuality,
        entity.callableLociCacheId
      )
    )
    entity

  override def update(entity: YProfileRegionEntity)(using conn: Connection): YProfileRegionEntity =
    executeUpdate(
      """UPDATE y_profile_region SET
        |  y_profile_id = ?, source_id = ?, contig = ?, start_position = ?, end_position = ?,
        |  callable_state = ?, mean_coverage = ?, mean_mapping_quality = ?, callable_loci_cache_id = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.yProfileId,
        entity.sourceId,
        entity.contig,
        entity.startPosition,
        entity.endPosition,
        entity.callableState.toString,
        entity.meanCoverage,
        entity.meanMappingQuality,
        entity.callableLociCacheId,
        entity.id
      )
    )
    entity

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM y_profile_region WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_profile_region") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_profile_region WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Region-Specific Queries
  // ============================================

  /**
   * Find all regions for a Y chromosome profile.
   */
  def findByProfile(yProfileId: UUID)(using conn: Connection): List[YProfileRegionEntity] =
    queryList(
      "SELECT * FROM y_profile_region WHERE y_profile_id = ? ORDER BY start_position",
      Seq(yProfileId)
    )(mapRow)

  /**
   * Find regions for a specific source.
   */
  def findBySource(sourceId: UUID)(using conn: Connection): List[YProfileRegionEntity] =
    queryList(
      "SELECT * FROM y_profile_region WHERE source_id = ? ORDER BY start_position",
      Seq(sourceId)
    )(mapRow)

  /**
   * Find regions by callable state.
   */
  def findByState(yProfileId: UUID, state: YCallableState)(using conn: Connection): List[YProfileRegionEntity] =
    queryList(
      "SELECT * FROM y_profile_region WHERE y_profile_id = ? AND callable_state = ? ORDER BY start_position",
      Seq(yProfileId, state.toString)
    )(mapRow)

  /**
   * Find regions overlapping a position.
   */
  def findOverlapping(yProfileId: UUID, position: Long)(using conn: Connection): List[YProfileRegionEntity] =
    queryList(
      """SELECT * FROM y_profile_region
        |WHERE y_profile_id = ? AND start_position <= ? AND end_position >= ?
        |ORDER BY start_position
      """.stripMargin,
      Seq(yProfileId, position, position)
    )(mapRow)

  /**
   * Find regions overlapping a range.
   */
  def findOverlappingRange(yProfileId: UUID, start: Long, end: Long)(using conn: Connection): List[YProfileRegionEntity] =
    queryList(
      """SELECT * FROM y_profile_region
        |WHERE y_profile_id = ? AND start_position <= ? AND end_position >= ?
        |ORDER BY start_position
      """.stripMargin,
      Seq(yProfileId, end, start)
    )(mapRow)

  /**
   * Delete all regions for a profile.
   */
  def deleteByProfile(yProfileId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_profile_region WHERE y_profile_id = ?", Seq(yProfileId))

  /**
   * Delete all regions for a source.
   */
  def deleteBySource(sourceId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM y_profile_region WHERE source_id = ?", Seq(sourceId))

  /**
   * Batch insert regions.
   */
  def insertBatch(regions: List[YProfileRegionEntity])(using conn: Connection): List[YProfileRegionEntity] =
    regions.foreach(insert)
    regions

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YProfileRegionEntity =
    YProfileRegionEntity(
      id = getUUID(rs, "id"),
      yProfileId = getUUID(rs, "y_profile_id"),
      sourceId = getUUID(rs, "source_id"),
      contig = rs.getString("contig"),
      startPosition = rs.getLong("start_position"),
      endPosition = rs.getLong("end_position"),
      callableState = YCallableState.fromString(rs.getString("callable_state")),
      meanCoverage = getOptDouble(rs, "mean_coverage"),
      meanMappingQuality = getOptDouble(rs, "mean_mapping_quality"),
      callableLociCacheId = getOptUUID(rs, "callable_loci_cache_id")
    )
