package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.IbdSegment
import io.circe.*
import io.circe.parser.*
import io.circe.syntax.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Circe codecs for IbdSegment JSON persistence.
 */
object IbdSegmentCodecs:
  given Encoder[IbdSegment] = Encoder.instance { s =>
    Json.obj(
      "chromosome" -> Json.fromString(s.chromosome),
      "startPosition" -> Json.fromInt(s.startPosition),
      "endPosition" -> Json.fromInt(s.endPosition),
      "lengthCm" -> Json.fromDoubleOrNull(s.lengthCm),
      "snpCount" -> s.snpCount.fold(Json.Null)(Json.fromInt),
      "isHalfIdentical" -> s.isHalfIdentical.fold(Json.Null)(Json.fromBoolean)
    )
  }

  given Decoder[IbdSegment] = Decoder.instance { c =>
    for
      chromosome <- c.get[String]("chromosome")
      startPosition <- c.get[Int]("startPosition")
      endPosition <- c.get[Int]("endPosition")
      lengthCm <- c.get[Double]("lengthCm")
      snpCount <- c.get[Option[Int]]("snpCount")
      isHalfIdentical <- c.get[Option[Boolean]]("isHalfIdentical")
    yield IbdSegment(chromosome, startPosition, endPosition, lengthCm, snpCount, isHalfIdentical)
  }

/**
 * Match Result entity for database persistence.
 * Maps to the match_result table (V012 migration).
 */
case class MatchResultEntity(
                              id: UUID,
                              biosampleId: UUID,
                              matchedBiosampleRef: String,
                              matchedCitizenDid: Option[String],
                              relationshipEstimate: Option[String],
                              totalSharedCm: Double,
                              longestSegmentCm: Option[Double],
                              segmentCount: Int,
                              sharedSegments: List[IbdSegment],
                              xMatchSharedCm: Option[Double],
                              matchedAt: LocalDateTime,
                              attestationHash: Option[String],
                              meta: EntityMeta
                            ) extends Entity[UUID]

object MatchResultEntity:
  def create(
              biosampleId: UUID,
              matchedBiosampleRef: String,
              totalSharedCm: Double,
              segmentCount: Int,
              matchedCitizenDid: Option[String] = None,
              relationshipEstimate: Option[String] = None,
              longestSegmentCm: Option[Double] = None,
              sharedSegments: List[IbdSegment] = List.empty,
              xMatchSharedCm: Option[Double] = None,
              attestationHash: Option[String] = None
            ): MatchResultEntity = MatchResultEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    matchedBiosampleRef = matchedBiosampleRef,
    matchedCitizenDid = matchedCitizenDid,
    relationshipEstimate = relationshipEstimate,
    totalSharedCm = totalSharedCm,
    longestSegmentCm = longestSegmentCm,
    segmentCount = segmentCount,
    sharedSegments = sharedSegments,
    xMatchSharedCm = xMatchSharedCm,
    matchedAt = LocalDateTime.now(),
    attestationHash = attestationHash,
    meta = EntityMeta.create()
  )

/**
 * Repository for Match Result persistence operations.
 */
class MatchResultRepository extends SyncableRepository[MatchResultEntity, UUID]:

  import IbdSegmentCodecs.given

  override def findById(id: UUID)(using conn: Connection): Option[MatchResultEntity] =
    queryOne("SELECT * FROM match_result WHERE id = ?", Seq(id))(mapRow)

  override def findAll()(using conn: Connection): List[MatchResultEntity] =
    queryList("SELECT * FROM match_result ORDER BY total_shared_cm DESC")(mapRow)

  override def insert(entity: MatchResultEntity)(using conn: Connection): MatchResultEntity =
    val segmentsJson = JsonValue(entity.sharedSegments.asJson.noSpaces)
    executeUpdate(
      """INSERT INTO match_result (
        |  id, biosample_id, matched_biosample_ref, matched_citizen_did,
        |  relationship_estimate, total_shared_cm, longest_segment_cm, segment_count,
        |  shared_segments, x_match_shared_cm, matched_at, attestation_hash,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id, entity.biosampleId, entity.matchedBiosampleRef, entity.matchedCitizenDid,
        entity.relationshipEstimate, entity.totalSharedCm, entity.longestSegmentCm, entity.segmentCount,
        segmentsJson, entity.xMatchSharedCm, entity.matchedAt, entity.attestationHash,
        entity.meta.syncStatus, entity.meta.atUri, entity.meta.atCid,
        entity.meta.version, entity.meta.createdAt, entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: MatchResultEntity)(using conn: Connection): MatchResultEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val segmentsJson = JsonValue(entity.sharedSegments.asJson.noSpaces)
    executeUpdate(
      """UPDATE match_result SET
        |  matched_citizen_did = ?, relationship_estimate = ?,
        |  total_shared_cm = ?, longest_segment_cm = ?, segment_count = ?,
        |  shared_segments = ?, x_match_shared_cm = ?, attestation_hash = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.matchedCitizenDid, entity.relationshipEstimate,
        entity.totalSharedCm, entity.longestSegmentCm, entity.segmentCount,
        segmentsJson, entity.xMatchSharedCm, entity.attestationHash,
        updatedMeta.syncStatus, updatedMeta.atUri, updatedMeta.atCid,
        updatedMeta.version, updatedMeta.updatedAt,
        entity.id
      )
    )
    entity.copy(meta = updatedMeta)

  override def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM match_result WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM match_result")(_.getLong(1)).getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM match_result WHERE id = ?", Seq(id))(_ => true).isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[MatchResultEntity] =
    queryList("SELECT * FROM match_result WHERE sync_status = ?", Seq(status.toString))(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE match_result SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE match_result SET sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ? WHERE id = ?",
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // MatchResult-Specific Queries
  // ============================================

  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[MatchResultEntity] =
    queryList(
      "SELECT * FROM match_result WHERE biosample_id = ? ORDER BY total_shared_cm DESC",
      Seq(biosampleId)
    )(mapRow)

  def findByBiosampleAboveThreshold(biosampleId: UUID, minCm: Double)(using conn: Connection): List[MatchResultEntity] =
    queryList(
      "SELECT * FROM match_result WHERE biosample_id = ? AND total_shared_cm >= ? ORDER BY total_shared_cm DESC",
      Seq(biosampleId, minCm)
    )(mapRow)

  def findByMatchedBiosample(matchedRef: String)(using conn: Connection): List[MatchResultEntity] =
    queryList(
      "SELECT * FROM match_result WHERE matched_biosample_ref = ? ORDER BY total_shared_cm DESC",
      Seq(matchedRef)
    )(mapRow)

  def deleteByBiosample(biosampleId: UUID)(using conn: Connection): Int =
    executeUpdate("DELETE FROM match_result WHERE biosample_id = ?", Seq(biosampleId))

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): MatchResultEntity =
    val segmentsJson = getOptJsonString(rs, "shared_segments")
    val segments = segmentsJson.flatMap(json =>
      parse(json).flatMap(_.as[List[IbdSegment]]).toOption
    ).getOrElse(List.empty)

    MatchResultEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      matchedBiosampleRef = rs.getString("matched_biosample_ref"),
      matchedCitizenDid = getOptString(rs, "matched_citizen_did"),
      relationshipEstimate = getOptString(rs, "relationship_estimate"),
      totalSharedCm = rs.getDouble("total_shared_cm"),
      longestSegmentCm = getOptDouble(rs, "longest_segment_cm"),
      segmentCount = rs.getInt("segment_count"),
      sharedSegments = segments,
      xMatchSharedCm = getOptDouble(rs, "x_match_shared_cm"),
      matchedAt = getDateTime(rs, "matched_at"),
      attestationHash = getOptString(rs, "attestation_hash"),
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
