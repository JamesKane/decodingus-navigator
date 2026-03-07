package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import io.circe.*
import io.circe.parser.*
import io.circe.syntax.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Match Consent entity for database persistence.
 * Maps to the match_consent table (V012 migration).
 */
case class MatchConsentEntity(
                               id: UUID,
                               biosampleId: UUID,
                               consentLevel: String,
                               allowedMatchTypes: List[String],
                               minimumSegmentCm: Double,
                               shareContactInfo: Boolean,
                               consentedAt: LocalDateTime,
                               expiresAt: Option[LocalDateTime],
                               meta: EntityMeta
                             ) extends Entity[UUID]

object MatchConsentEntity:
  def create(
              biosampleId: UUID,
              consentLevel: String,
              allowedMatchTypes: List[String] = List("IBD"),
              minimumSegmentCm: Double = 7.0,
              shareContactInfo: Boolean = false,
              expiresAt: Option[LocalDateTime] = None
            ): MatchConsentEntity = MatchConsentEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    consentLevel = consentLevel,
    allowedMatchTypes = allowedMatchTypes,
    minimumSegmentCm = minimumSegmentCm,
    shareContactInfo = shareContactInfo,
    consentedAt = LocalDateTime.now(),
    expiresAt = expiresAt,
    meta = EntityMeta.create()
  )

/**
 * Repository for Match Consent persistence operations.
 */
class MatchConsentRepository extends SyncableRepository[MatchConsentEntity, UUID]:

  override def findById(id: UUID)(using conn: Connection): Option[MatchConsentEntity] =
    queryOne("SELECT * FROM match_consent WHERE id = ?", Seq(id))(mapRow)

  override def findAll()(using conn: Connection): List[MatchConsentEntity] =
    queryList("SELECT * FROM match_consent ORDER BY consented_at DESC")(mapRow)

  override def insert(entity: MatchConsentEntity)(using conn: Connection): MatchConsentEntity =
    val matchTypesJson = JsonValue(entity.allowedMatchTypes.asJson.noSpaces)
    executeUpdate(
      """INSERT INTO match_consent (
        |  id, biosample_id, consent_level, allowed_match_types, minimum_segment_cm,
        |  share_contact_info, consented_at, expires_at,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id, entity.biosampleId, entity.consentLevel, matchTypesJson,
        entity.minimumSegmentCm, entity.shareContactInfo, entity.consentedAt,
        entity.expiresAt,
        entity.meta.syncStatus, entity.meta.atUri, entity.meta.atCid,
        entity.meta.version, entity.meta.createdAt, entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: MatchConsentEntity)(using conn: Connection): MatchConsentEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val matchTypesJson = JsonValue(entity.allowedMatchTypes.asJson.noSpaces)
    executeUpdate(
      """UPDATE match_consent SET
        |  consent_level = ?, allowed_match_types = ?, minimum_segment_cm = ?,
        |  share_contact_info = ?, expires_at = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.consentLevel, matchTypesJson, entity.minimumSegmentCm,
        entity.shareContactInfo, entity.expiresAt,
        updatedMeta.syncStatus, updatedMeta.atUri, updatedMeta.atCid,
        updatedMeta.version, updatedMeta.updatedAt,
        entity.id
      )
    )
    entity.copy(meta = updatedMeta)

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM match_consent WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM match_consent")(_.getLong(1)).getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM match_consent WHERE id = ?", Seq(id))(_ => true).isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[MatchConsentEntity] =
    queryList("SELECT * FROM match_consent WHERE sync_status = ?", Seq(status.toString))(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE match_consent SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE match_consent SET sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ? WHERE id = ?",
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // MatchConsent-Specific Queries
  // ============================================

  def findByBiosample(biosampleId: UUID)(using conn: Connection): Option[MatchConsentEntity] =
    queryOne("SELECT * FROM match_consent WHERE biosample_id = ?", Seq(biosampleId))(mapRow)

  def deleteByBiosample(biosampleId: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM match_consent WHERE biosample_id = ?", Seq(biosampleId)) > 0

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): MatchConsentEntity =
    val matchTypesJson = getOptJsonString(rs, "allowed_match_types")
    val matchTypes = matchTypesJson.flatMap(json =>
      parse(json).flatMap(_.as[List[String]]).toOption
    ).getOrElse(List("IBD"))

    MatchConsentEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      consentLevel = rs.getString("consent_level"),
      allowedMatchTypes = matchTypes,
      minimumSegmentCm = rs.getDouble("minimum_segment_cm"),
      shareContactInfo = rs.getBoolean("share_contact_info"),
      consentedAt = getDateTime(rs, "consented_at"),
      expiresAt = getOptDateTime(rs, "expires_at"),
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
