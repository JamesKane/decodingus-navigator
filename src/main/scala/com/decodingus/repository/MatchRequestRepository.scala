package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Match Request entity for database persistence.
 * Maps to the match_request table (V012 migration).
 */
case class MatchRequestEntity(
                                id: UUID,
                                fromBiosampleRef: String,
                                toBiosampleRef: String,
                                status: String,
                                requestType: String,
                                message: Option[String],
                                sharedAncestorHint: Option[String],
                                discoveryReason: Option[String],
                                createdAt: LocalDateTime,
                                expiresAt: Option[LocalDateTime],
                                respondedAt: Option[LocalDateTime],
                                meta: EntityMeta
                              ) extends Entity[UUID]

object MatchRequestEntity:
  def create(
              fromBiosampleRef: String,
              toBiosampleRef: String,
              requestType: String = "AUTOSOMAL",
              message: Option[String] = None,
              sharedAncestorHint: Option[String] = None,
              discoveryReason: Option[String] = None,
              expiresAt: Option[LocalDateTime] = None
            ): MatchRequestEntity = MatchRequestEntity(
    id = UUID.randomUUID(),
    fromBiosampleRef = fromBiosampleRef,
    toBiosampleRef = toBiosampleRef,
    status = "PENDING",
    requestType = requestType,
    message = message,
    sharedAncestorHint = sharedAncestorHint,
    discoveryReason = discoveryReason,
    createdAt = LocalDateTime.now(),
    expiresAt = expiresAt.orElse(Some(LocalDateTime.now().plusDays(30))),
    respondedAt = None,
    meta = EntityMeta.create()
  )

/**
 * Repository for Match Request persistence operations.
 */
class MatchRequestRepository extends SyncableRepository[MatchRequestEntity, UUID]:

  override def findById(id: UUID)(using conn: Connection): Option[MatchRequestEntity] =
    queryOne("SELECT * FROM match_request WHERE id = ?", Seq(id))(mapRow)

  override def findAll()(using conn: Connection): List[MatchRequestEntity] =
    queryList("SELECT * FROM match_request ORDER BY created_at DESC")(mapRow)

  override def insert(entity: MatchRequestEntity)(using conn: Connection): MatchRequestEntity =
    executeUpdate(
      """INSERT INTO match_request (
        |  id, from_biosample_ref, to_biosample_ref, status, request_type,
        |  message, shared_ancestor_hint, discovery_reason,
        |  created_at, expires_at, responded_at,
        |  sync_status, at_uri, at_cid, version, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id, entity.fromBiosampleRef, entity.toBiosampleRef,
        entity.status, entity.requestType,
        entity.message, entity.sharedAncestorHint, entity.discoveryReason,
        entity.createdAt, entity.expiresAt, entity.respondedAt,
        entity.meta.syncStatus, entity.meta.atUri, entity.meta.atCid,
        entity.meta.version, entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: MatchRequestEntity)(using conn: Connection): MatchRequestEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    executeUpdate(
      """UPDATE match_request SET
        |  status = ?, message = ?, shared_ancestor_hint = ?,
        |  expires_at = ?, responded_at = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.status, entity.message, entity.sharedAncestorHint,
        entity.expiresAt, entity.respondedAt,
        updatedMeta.syncStatus, updatedMeta.atUri, updatedMeta.atCid,
        updatedMeta.version, updatedMeta.updatedAt,
        entity.id
      )
    )
    entity.copy(meta = updatedMeta)

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM match_request WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM match_request")(_.getLong(1)).getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM match_request WHERE id = ?", Seq(id))(_ => true).isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[MatchRequestEntity] =
    queryList("SELECT * FROM match_request WHERE sync_status = ?", Seq(status.toString))(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE match_request SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE match_request SET sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ? WHERE id = ?",
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // MatchRequest-Specific Queries
  // ============================================

  def findByFromBiosample(biosampleRef: String)(using conn: Connection): List[MatchRequestEntity] =
    queryList(
      "SELECT * FROM match_request WHERE from_biosample_ref = ? ORDER BY created_at DESC",
      Seq(biosampleRef)
    )(mapRow)

  def findByToBiosample(biosampleRef: String)(using conn: Connection): List[MatchRequestEntity] =
    queryList(
      "SELECT * FROM match_request WHERE to_biosample_ref = ? ORDER BY created_at DESC",
      Seq(biosampleRef)
    )(mapRow)

  def findPending(biosampleRef: String)(using conn: Connection): List[MatchRequestEntity] =
    queryList(
      """SELECT * FROM match_request
        |WHERE to_biosample_ref = ? AND status = 'PENDING'
        |AND (expires_at IS NULL OR expires_at > ?)
        |ORDER BY created_at DESC
      """.stripMargin,
      Seq(biosampleRef, LocalDateTime.now())
    )(mapRow)

  def updateRequestStatus(id: UUID, status: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    val respondedAt = if status != "PENDING" then Some(now) else None
    executeUpdate(
      "UPDATE match_request SET status = ?, responded_at = ?, updated_at = ? WHERE id = ?",
      Seq(status, respondedAt, now, id)
    ) > 0

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): MatchRequestEntity =
    MatchRequestEntity(
      id = getUUID(rs, "id"),
      fromBiosampleRef = rs.getString("from_biosample_ref"),
      toBiosampleRef = rs.getString("to_biosample_ref"),
      status = rs.getString("status"),
      requestType = rs.getString("request_type"),
      message = getOptString(rs, "message"),
      sharedAncestorHint = getOptString(rs, "shared_ancestor_hint"),
      discoveryReason = getOptString(rs, "discovery_reason"),
      createdAt = getDateTime(rs, "created_at"),
      expiresAt = getOptDateTime(rs, "expires_at"),
      respondedAt = getOptDateTime(rs, "responded_at"),
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
