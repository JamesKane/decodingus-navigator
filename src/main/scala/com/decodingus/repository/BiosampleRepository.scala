package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.{HaplogroupAssignments, HaplogroupResult, PrivateVariantData, VariantCall}
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Biosample entity for database persistence.
 *
 * This is a database-layer entity that wraps biosample data with
 * persistence metadata (sync status, version, timestamps).
 */
case class BiosampleEntity(
  id: UUID,
  sampleAccession: String,
  donorIdentifier: String,
  description: Option[String],
  centerName: Option[String],
  sex: Option[String],
  citizenDid: Option[String],
  haplogroups: Option[HaplogroupAssignments],
  meta: EntityMeta
) extends Entity[UUID]

object BiosampleEntity:
  // Circe codecs for JSON serialization of haplogroups
  // Note: These could be moved to a shared Codecs object if needed elsewhere
  given Encoder[VariantCall] = deriveEncoder
  given Decoder[VariantCall] = deriveDecoder
  given Encoder[PrivateVariantData] = deriveEncoder
  given Decoder[PrivateVariantData] = deriveDecoder
  given Encoder[HaplogroupResult] = deriveEncoder
  given Decoder[HaplogroupResult] = deriveDecoder
  given Encoder[HaplogroupAssignments] = deriveEncoder
  given Decoder[HaplogroupAssignments] = deriveDecoder

  /**
   * Create a new BiosampleEntity with generated ID and initial metadata.
   */
  def create(
    sampleAccession: String,
    donorIdentifier: String,
    description: Option[String] = None,
    centerName: Option[String] = None,
    sex: Option[String] = None,
    citizenDid: Option[String] = None,
    haplogroups: Option[HaplogroupAssignments] = None
  ): BiosampleEntity = BiosampleEntity(
    id = UUID.randomUUID(),
    sampleAccession = sampleAccession,
    donorIdentifier = donorIdentifier,
    description = description,
    centerName = centerName,
    sex = sex,
    citizenDid = citizenDid,
    haplogroups = haplogroups,
    meta = EntityMeta.create()
  )

/**
 * Repository for biosample persistence operations.
 */
class BiosampleRepository extends SyncableRepositoryBase[BiosampleEntity]:

  import BiosampleEntity.given

  override protected def tableName: String = "biosample"

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[BiosampleEntity] =
    queryOne(
      "SELECT * FROM biosample WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[BiosampleEntity] =
    queryList("SELECT * FROM biosample ORDER BY created_at DESC")(mapRow)

  override def insert(entity: BiosampleEntity)(using conn: Connection): BiosampleEntity =
    val haplogroupsJson = entity.haplogroups.map(h => JsonValue(h.asJson.noSpaces))

    executeUpdate(
      """INSERT INTO biosample (
        |  id, sample_accession, donor_identifier, description, center_name, sex,
        |  citizen_did, haplogroups, sync_status, at_uri, at_cid, version,
        |  created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.sampleAccession,
        entity.donorIdentifier,
        entity.description,
        entity.centerName,
        entity.sex,
        entity.citizenDid,
        haplogroupsJson,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: BiosampleEntity)(using conn: Connection): BiosampleEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val haplogroupsJson = entity.haplogroups.map(h => JsonValue(h.asJson.noSpaces))

    executeUpdate(
      """UPDATE biosample SET
        |  sample_accession = ?, donor_identifier = ?, description = ?, center_name = ?,
        |  sex = ?, citizen_did = ?, haplogroups = ?, sync_status = ?, at_uri = ?,
        |  at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.sampleAccession,
        entity.donorIdentifier,
        entity.description,
        entity.centerName,
        entity.sex,
        entity.citizenDid,
        haplogroupsJson,
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
    executeUpdate("DELETE FROM biosample WHERE id = ?", Seq(id)) > 0

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM biosample WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Biosample-Specific Queries
  // ============================================

  /**
   * Find a biosample by its sample accession.
   */
  def findByAccession(accession: String)(using conn: Connection): Option[BiosampleEntity] =
    queryOne(
      "SELECT * FROM biosample WHERE sample_accession = ?",
      Seq(accession)
    )(mapRow)

  /**
   * Find all biosamples for a donor.
   */
  def findByDonor(donorIdentifier: String)(using conn: Connection): List[BiosampleEntity] =
    queryList(
      "SELECT * FROM biosample WHERE donor_identifier = ? ORDER BY created_at DESC",
      Seq(donorIdentifier)
    )(mapRow)

  /**
   * Find all biosamples owned by a citizen.
   */
  def findByCitizen(citizenDid: String)(using conn: Connection): List[BiosampleEntity] =
    queryList(
      "SELECT * FROM biosample WHERE citizen_did = ? ORDER BY created_at DESC",
      Seq(citizenDid)
    )(mapRow)

  /**
   * Update haplogroup assignments for a biosample.
   */
  def updateHaplogroups(id: UUID, haplogroups: HaplogroupAssignments)(using conn: Connection): Boolean =
    val json = JsonValue(haplogroups.asJson.noSpaces)
    executeUpdate(
      """UPDATE biosample SET
        |  haplogroups = ?, sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
        |  updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(json, LocalDateTime.now(), id)
    ) > 0

  /**
   * Search biosamples by sample accession prefix.
   */
  def searchByAccession(prefix: String)(using conn: Connection): List[BiosampleEntity] =
    queryList(
      "SELECT * FROM biosample WHERE sample_accession LIKE ? ORDER BY sample_accession",
      Seq(s"$prefix%")
    )(mapRow)

  // ============================================
  // Result Set Mapping
  // ============================================

  override protected def mapRow(rs: ResultSet): BiosampleEntity =
    val haplogroupsJson = getOptJsonString(rs, "haplogroups")
    val haplogroups = haplogroupsJson.flatMap { json =>
      parse(json).flatMap(_.as[HaplogroupAssignments]).toOption
    }

    BiosampleEntity(
      id = getUUID(rs, "id"),
      sampleAccession = rs.getString("sample_accession"),
      donorIdentifier = rs.getString("donor_identifier"),
      description = getOptString(rs, "description"),
      centerName = getOptString(rs, "center_name"),
      sex = getOptString(rs, "sex"),
      citizenDid = getOptString(rs, "citizen_did"),
      haplogroups = haplogroups,
      meta = readEntityMeta(rs)
    )
