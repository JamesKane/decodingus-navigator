package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.{ChipProfile, FileInfo}
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Chip Profile entity for database persistence.
 *
 * Stores DNA microarray (chip) testing results with marker statistics.
 */
case class ChipProfileEntity(
  id: UUID,
  biosampleId: UUID,
  vendor: String,
  testTypeCode: String,
  chipVersion: Option[String],
  totalMarkersCalled: Int,
  totalMarkersPossible: Int,
  noCallRate: Double,
  yMarkersCalled: Option[Int],
  mtMarkersCalled: Option[Int],
  autosomalMarkersCalled: Int,
  hetRate: Option[Double],
  importDate: LocalDateTime,
  sourceFileHash: Option[String],
  sourceFileName: Option[String],
  files: List[FileInfo],
  meta: EntityMeta
) extends Entity[UUID]

object ChipProfileEntity:
  import ChipProfileCodecs.given

  /**
   * Create a new ChipProfileEntity with generated ID and initial metadata.
   */
  def create(
    biosampleId: UUID,
    vendor: String,
    testTypeCode: String,
    totalMarkersCalled: Int,
    totalMarkersPossible: Int,
    noCallRate: Double,
    autosomalMarkersCalled: Int,
    importDate: LocalDateTime,
    chipVersion: Option[String] = None,
    yMarkersCalled: Option[Int] = None,
    mtMarkersCalled: Option[Int] = None,
    hetRate: Option[Double] = None,
    sourceFileHash: Option[String] = None,
    sourceFileName: Option[String] = None,
    files: List[FileInfo] = List.empty
  ): ChipProfileEntity = ChipProfileEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    vendor = vendor,
    testTypeCode = testTypeCode,
    chipVersion = chipVersion,
    totalMarkersCalled = totalMarkersCalled,
    totalMarkersPossible = totalMarkersPossible,
    noCallRate = noCallRate,
    yMarkersCalled = yMarkersCalled,
    mtMarkersCalled = mtMarkersCalled,
    autosomalMarkersCalled = autosomalMarkersCalled,
    hetRate = hetRate,
    importDate = importDate,
    sourceFileHash = sourceFileHash,
    sourceFileName = sourceFileName,
    files = files,
    meta = EntityMeta.create()
  )

  /**
   * Create entity from a workspace model ChipProfile.
   */
  def fromModel(biosampleId: UUID, profile: ChipProfile): ChipProfileEntity =
    ChipProfileEntity(
      id = UUID.randomUUID(),
      biosampleId = biosampleId,
      vendor = profile.vendor,
      testTypeCode = profile.testTypeCode,
      chipVersion = profile.chipVersion,
      totalMarkersCalled = profile.totalMarkersCalled,
      totalMarkersPossible = profile.totalMarkersPossible,
      noCallRate = profile.noCallRate,
      yMarkersCalled = profile.yMarkersCalled,
      mtMarkersCalled = profile.mtMarkersCalled,
      autosomalMarkersCalled = profile.autosomalMarkersCalled,
      hetRate = profile.hetRate,
      importDate = profile.importDate,
      sourceFileHash = profile.sourceFileHash,
      sourceFileName = profile.sourceFileName,
      files = profile.files,
      meta = EntityMeta.create()
    )

/**
 * Circe codecs for Chip profile JSON fields.
 */
object ChipProfileCodecs:
  given Encoder[FileInfo] = deriveEncoder
  given Decoder[FileInfo] = deriveDecoder

/**
 * Repository for Chip profile persistence operations.
 */
class ChipProfileRepository extends SyncableRepositoryBase[ChipProfileEntity]:

  import ChipProfileCodecs.given

  override protected def tableName: String = "chip_profile"

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[ChipProfileEntity] =
    queryOne(
      "SELECT * FROM chip_profile WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[ChipProfileEntity] =
    queryList("SELECT * FROM chip_profile ORDER BY created_at DESC")(mapRow)

  override def insert(entity: ChipProfileEntity)(using conn: Connection): ChipProfileEntity =
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO chip_profile (
        |  id, biosample_id, vendor, test_type_code, chip_version,
        |  total_markers_called, total_markers_possible, no_call_rate,
        |  y_markers_called, mt_markers_called, autosomal_markers_called, het_rate,
        |  import_date, source_file_hash, source_file_name, files,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.vendor,
        entity.testTypeCode,
        entity.chipVersion,
        entity.totalMarkersCalled,
        entity.totalMarkersPossible,
        entity.noCallRate,
        entity.yMarkersCalled,
        entity.mtMarkersCalled,
        entity.autosomalMarkersCalled,
        entity.hetRate,
        entity.importDate,
        entity.sourceFileHash,
        entity.sourceFileName,
        filesJson,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: ChipProfileEntity)(using conn: Connection): ChipProfileEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """UPDATE chip_profile SET
        |  biosample_id = ?, vendor = ?, test_type_code = ?, chip_version = ?,
        |  total_markers_called = ?, total_markers_possible = ?, no_call_rate = ?,
        |  y_markers_called = ?, mt_markers_called = ?, autosomal_markers_called = ?, het_rate = ?,
        |  import_date = ?, source_file_hash = ?, source_file_name = ?, files = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.vendor,
        entity.testTypeCode,
        entity.chipVersion,
        entity.totalMarkersCalled,
        entity.totalMarkersPossible,
        entity.noCallRate,
        entity.yMarkersCalled,
        entity.mtMarkersCalled,
        entity.autosomalMarkersCalled,
        entity.hetRate,
        entity.importDate,
        entity.sourceFileHash,
        entity.sourceFileName,
        filesJson,
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
    executeUpdate("DELETE FROM chip_profile WHERE id = ?", Seq(id)) > 0

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM chip_profile WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Chip Profile-Specific Queries
  // ============================================

  /**
   * Find all chip profiles for a biosample.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[ChipProfileEntity] =
    queryList(
      "SELECT * FROM chip_profile WHERE biosample_id = ? ORDER BY import_date DESC",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Find chip profiles by vendor.
   */
  def findByVendor(vendor: String)(using conn: Connection): List[ChipProfileEntity] =
    queryList(
      "SELECT * FROM chip_profile WHERE vendor = ? ORDER BY import_date DESC",
      Seq(vendor)
    )(mapRow)

  /**
   * Find chip profiles by test type.
   */
  def findByTestType(testTypeCode: String)(using conn: Connection): List[ChipProfileEntity] =
    queryList(
      "SELECT * FROM chip_profile WHERE test_type_code = ? ORDER BY import_date DESC",
      Seq(testTypeCode)
    )(mapRow)

  /**
   * Find chip profile by source file hash (for deduplication).
   */
  def findBySourceFileHash(hash: String)(using conn: Connection): Option[ChipProfileEntity] =
    queryOne(
      "SELECT * FROM chip_profile WHERE source_file_hash = ?",
      Seq(hash)
    )(mapRow)

  // ============================================
  // Result Set Mapping
  // ============================================

  override protected def mapRow(rs: ResultSet): ChipProfileEntity =
    val filesJson = getOptJsonString(rs, "files")
    val files = filesJson.flatMap(json => parse(json).flatMap(_.as[List[FileInfo]]).toOption).getOrElse(List.empty)

    ChipProfileEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      vendor = rs.getString("vendor"),
      testTypeCode = rs.getString("test_type_code"),
      chipVersion = getOptString(rs, "chip_version"),
      totalMarkersCalled = rs.getInt("total_markers_called"),
      totalMarkersPossible = rs.getInt("total_markers_possible"),
      noCallRate = rs.getDouble("no_call_rate"),
      yMarkersCalled = getOptInt(rs, "y_markers_called"),
      mtMarkersCalled = getOptInt(rs, "mt_markers_called"),
      autosomalMarkersCalled = rs.getInt("autosomal_markers_called"),
      hetRate = getOptDouble(rs, "het_rate"),
      importDate = getDateTime(rs, "import_date"),
      sourceFileHash = getOptString(rs, "source_file_hash"),
      sourceFileName = getOptString(rs, "source_file_name"),
      files = files,
      meta = readEntityMeta(rs)
    )
