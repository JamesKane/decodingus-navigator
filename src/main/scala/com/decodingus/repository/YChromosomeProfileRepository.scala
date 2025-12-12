package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.yprofile.model.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Repository for Y chromosome profile persistence operations.
 */
class YChromosomeProfileRepository extends SyncableRepository[YChromosomeProfileEntity, UUID]:

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YChromosomeProfileEntity] =
    queryOne(
      "SELECT * FROM y_chromosome_profile WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YChromosomeProfileEntity] =
    queryList("SELECT * FROM y_chromosome_profile ORDER BY created_at DESC")(mapRow)

  override def insert(entity: YChromosomeProfileEntity)(using conn: Connection): YChromosomeProfileEntity =
    executeUpdate(
      """INSERT INTO y_chromosome_profile (
        |  id, biosample_id, consensus_haplogroup, haplogroup_confidence,
        |  haplogroup_tree_provider, haplogroup_tree_version,
        |  total_variants, confirmed_count, novel_count, conflict_count, no_coverage_count,
        |  str_marker_count, str_confirmed_count,
        |  overall_confidence, callable_region_pct, mean_coverage,
        |  source_count, primary_source_type, last_reconciled_at,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.consensusHaplogroup,
        entity.haplogroupConfidence,
        entity.haplogroupTreeProvider,
        entity.haplogroupTreeVersion,
        entity.totalVariants,
        entity.confirmedCount,
        entity.novelCount,
        entity.conflictCount,
        entity.noCoverageCount,
        entity.strMarkerCount,
        entity.strConfirmedCount,
        entity.overallConfidence,
        entity.callableRegionPct,
        entity.meanCoverage,
        entity.sourceCount,
        entity.primarySourceType.map(_.toString),
        entity.lastReconciledAt,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: YChromosomeProfileEntity)(using conn: Connection): YChromosomeProfileEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)

    executeUpdate(
      """UPDATE y_chromosome_profile SET
        |  biosample_id = ?, consensus_haplogroup = ?, haplogroup_confidence = ?,
        |  haplogroup_tree_provider = ?, haplogroup_tree_version = ?,
        |  total_variants = ?, confirmed_count = ?, novel_count = ?, conflict_count = ?, no_coverage_count = ?,
        |  str_marker_count = ?, str_confirmed_count = ?,
        |  overall_confidence = ?, callable_region_pct = ?, mean_coverage = ?,
        |  source_count = ?, primary_source_type = ?, last_reconciled_at = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.consensusHaplogroup,
        entity.haplogroupConfidence,
        entity.haplogroupTreeProvider,
        entity.haplogroupTreeVersion,
        entity.totalVariants,
        entity.confirmedCount,
        entity.novelCount,
        entity.conflictCount,
        entity.noCoverageCount,
        entity.strMarkerCount,
        entity.strConfirmedCount,
        entity.overallConfidence,
        entity.callableRegionPct,
        entity.meanCoverage,
        entity.sourceCount,
        entity.primarySourceType.map(_.toString),
        entity.lastReconciledAt,
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
    executeUpdate("DELETE FROM y_chromosome_profile WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_chromosome_profile") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_chromosome_profile WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[YChromosomeProfileEntity] =
    queryList(
      "SELECT * FROM y_chromosome_profile WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE y_chromosome_profile SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE y_chromosome_profile SET
        |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // Profile-Specific Queries
  // ============================================

  /**
   * Find Y chromosome profile by biosample ID.
   * Returns at most one profile per biosample (unique constraint).
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): Option[YChromosomeProfileEntity] =
    queryOne(
      "SELECT * FROM y_chromosome_profile WHERE biosample_id = ?",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Find profiles by consensus haplogroup.
   */
  def findByHaplogroup(haplogroup: String)(using conn: Connection): List[YChromosomeProfileEntity] =
    queryList(
      "SELECT * FROM y_chromosome_profile WHERE consensus_haplogroup = ? ORDER BY created_at DESC",
      Seq(haplogroup)
    )(mapRow)

  /**
   * Find profiles under a haplogroup branch (prefix match).
   */
  def findByHaplogroupBranch(branchPrefix: String)(using conn: Connection): List[YChromosomeProfileEntity] =
    queryList(
      "SELECT * FROM y_chromosome_profile WHERE consensus_haplogroup LIKE ? ORDER BY consensus_haplogroup",
      Seq(s"$branchPrefix%")
    )(mapRow)

  /**
   * Find profiles with conflicts needing resolution.
   */
  def findWithConflicts()(using conn: Connection): List[YChromosomeProfileEntity] =
    queryList(
      "SELECT * FROM y_chromosome_profile WHERE conflict_count > 0 ORDER BY conflict_count DESC"
    )(mapRow)

  /**
   * Find profiles pending sync.
   */
  def findPendingSync()(using conn: Connection): List[YChromosomeProfileEntity] =
    queryList(
      """SELECT * FROM y_chromosome_profile
        |WHERE sync_status IN ('Local', 'Modified')
        |ORDER BY updated_at ASC
      """.stripMargin
    )(mapRow)

  /**
   * Update reconciliation timestamp.
   */
  def markReconciled(id: UUID)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE y_chromosome_profile SET last_reconciled_at = ?, updated_at = ? WHERE id = ?",
      Seq(LocalDateTime.now(), LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YChromosomeProfileEntity =
    val lastReconciledTs = rs.getTimestamp("last_reconciled_at")
    val lastReconciledAt = if rs.wasNull() then None else Some(lastReconciledTs.toLocalDateTime)

    val primarySourceTypeStr = getOptString(rs, "primary_source_type")
    val primarySourceType = primarySourceTypeStr.map(YProfileSourceType.fromString)

    YChromosomeProfileEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      consensusHaplogroup = getOptString(rs, "consensus_haplogroup"),
      haplogroupConfidence = getOptDouble(rs, "haplogroup_confidence"),
      haplogroupTreeProvider = getOptString(rs, "haplogroup_tree_provider"),
      haplogroupTreeVersion = getOptString(rs, "haplogroup_tree_version"),
      totalVariants = rs.getInt("total_variants"),
      confirmedCount = rs.getInt("confirmed_count"),
      novelCount = rs.getInt("novel_count"),
      conflictCount = rs.getInt("conflict_count"),
      noCoverageCount = rs.getInt("no_coverage_count"),
      strMarkerCount = rs.getInt("str_marker_count"),
      strConfirmedCount = rs.getInt("str_confirmed_count"),
      overallConfidence = getOptDouble(rs, "overall_confidence"),
      callableRegionPct = getOptDouble(rs, "callable_region_pct"),
      meanCoverage = getOptDouble(rs, "mean_coverage"),
      sourceCount = rs.getInt("source_count"),
      primarySourceType = primarySourceType,
      lastReconciledAt = lastReconciledAt,
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
