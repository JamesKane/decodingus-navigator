package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.{FileInfo, StrMarkerValue, StrPanel, StrProfile}
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * STR Profile entity for database persistence.
 *
 * Stores Y-STR marker profiles with panel information and values.
 */
case class StrProfileEntity(
  id: UUID,
  biosampleId: UUID,
  sequenceRunId: Option[UUID],
  panels: List[StrPanel],
  markers: List[StrMarkerValue],
  totalMarkers: Option[Int],
  source: Option[String],
  importedFrom: Option[String],
  derivationMethod: Option[String],
  files: List[FileInfo],
  meta: EntityMeta
) extends Entity[UUID]

object StrProfileEntity:
  import StrProfileCodecs.given

  /**
   * Create a new StrProfileEntity with generated ID and initial metadata.
   */
  def create(
    biosampleId: UUID,
    sequenceRunId: Option[UUID] = None,
    panels: List[StrPanel] = List.empty,
    markers: List[StrMarkerValue] = List.empty,
    totalMarkers: Option[Int] = None,
    source: Option[String] = None,
    importedFrom: Option[String] = None,
    derivationMethod: Option[String] = None,
    files: List[FileInfo] = List.empty
  ): StrProfileEntity = StrProfileEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    sequenceRunId = sequenceRunId,
    panels = panels,
    markers = markers,
    totalMarkers = totalMarkers,
    source = source,
    importedFrom = importedFrom,
    derivationMethod = derivationMethod,
    files = files,
    meta = EntityMeta.create()
  )

  /**
   * Create entity from a workspace model StrProfile.
   */
  def fromModel(biosampleId: UUID, profile: StrProfile): StrProfileEntity =
    StrProfileEntity(
      id = UUID.randomUUID(),
      biosampleId = biosampleId,
      sequenceRunId = None, // Would need to resolve from sequenceRunRef
      panels = profile.panels,
      markers = profile.markers,
      totalMarkers = profile.totalMarkers,
      source = profile.source,
      importedFrom = profile.importedFrom,
      derivationMethod = profile.derivationMethod,
      files = profile.files,
      meta = EntityMeta.create()
    )

/**
 * Circe codecs for STR profile JSON fields.
 */
object StrProfileCodecs:
  import com.decodingus.workspace.model.*

  // StrValue hierarchy
  given Encoder[SimpleStrValue] = deriveEncoder
  given Decoder[SimpleStrValue] = deriveDecoder
  given Encoder[MultiCopyStrValue] = deriveEncoder
  given Decoder[MultiCopyStrValue] = deriveDecoder
  given Encoder[StrAllele] = deriveEncoder
  given Decoder[StrAllele] = deriveDecoder
  given Encoder[ComplexStrValue] = deriveEncoder
  given Decoder[ComplexStrValue] = deriveDecoder

  given Encoder[StrValue] = Encoder.instance {
    case s: SimpleStrValue => s.asJson.deepMerge(Json.obj("_type" -> Json.fromString("simple")))
    case m: MultiCopyStrValue => m.asJson.deepMerge(Json.obj("_type" -> Json.fromString("multiCopy")))
    case c: ComplexStrValue => c.asJson.deepMerge(Json.obj("_type" -> Json.fromString("complex")))
  }

  given Decoder[StrValue] = Decoder.instance { cursor =>
    cursor.get[String]("_type").flatMap {
      case "simple" => cursor.as[SimpleStrValue]
      case "multiCopy" => cursor.as[MultiCopyStrValue]
      case "complex" => cursor.as[ComplexStrValue]
      case other => Left(DecodingFailure(s"Unknown StrValue type: $other", cursor.history))
    }
  }

  given Encoder[StrMarkerValue] = deriveEncoder
  given Decoder[StrMarkerValue] = deriveDecoder
  given Encoder[StrPanel] = deriveEncoder
  given Decoder[StrPanel] = deriveDecoder
  given Encoder[FileInfo] = deriveEncoder
  given Decoder[FileInfo] = deriveDecoder

/**
 * Repository for STR profile persistence operations.
 */
class StrProfileRepository extends SyncableRepository[StrProfileEntity, UUID]:

  import StrProfileCodecs.given

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[StrProfileEntity] =
    queryOne(
      "SELECT * FROM str_profile WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[StrProfileEntity] =
    queryList("SELECT * FROM str_profile ORDER BY created_at DESC")(mapRow)

  override def insert(entity: StrProfileEntity)(using conn: Connection): StrProfileEntity =
    val panelsJson = JsonValue(entity.panels.asJson.noSpaces)
    val markersJson = JsonValue(entity.markers.asJson.noSpaces)
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO str_profile (
        |  id, biosample_id, sequence_run_id, panels, markers, total_markers,
        |  source, imported_from, derivation_method, files,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.sequenceRunId,
        panelsJson,
        markersJson,
        entity.totalMarkers,
        entity.source,
        entity.importedFrom,
        entity.derivationMethod,
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

  override def update(entity: StrProfileEntity)(using conn: Connection): StrProfileEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val panelsJson = JsonValue(entity.panels.asJson.noSpaces)
    val markersJson = JsonValue(entity.markers.asJson.noSpaces)
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """UPDATE str_profile SET
        |  biosample_id = ?, sequence_run_id = ?, panels = ?, markers = ?, total_markers = ?,
        |  source = ?, imported_from = ?, derivation_method = ?, files = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.sequenceRunId,
        panelsJson,
        markersJson,
        entity.totalMarkers,
        entity.source,
        entity.importedFrom,
        entity.derivationMethod,
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
    executeUpdate("DELETE FROM str_profile WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM str_profile") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM str_profile WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[StrProfileEntity] =
    queryList(
      "SELECT * FROM str_profile WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE str_profile SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE str_profile SET
        |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // STR Profile-Specific Queries
  // ============================================

  /**
   * Find all STR profiles for a biosample.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[StrProfileEntity] =
    queryList(
      "SELECT * FROM str_profile WHERE biosample_id = ? ORDER BY created_at DESC",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Find STR profiles derived from a sequence run.
   */
  def findBySequenceRun(sequenceRunId: UUID)(using conn: Connection): List[StrProfileEntity] =
    queryList(
      "SELECT * FROM str_profile WHERE sequence_run_id = ? ORDER BY created_at DESC",
      Seq(sequenceRunId)
    )(mapRow)

  /**
   * Find STR profiles by source type.
   */
  def findBySource(source: String)(using conn: Connection): List[StrProfileEntity] =
    queryList(
      "SELECT * FROM str_profile WHERE source = ? ORDER BY created_at DESC",
      Seq(source)
    )(mapRow)

  /**
   * Find STR profiles imported from a specific provider.
   */
  def findByImportedFrom(provider: String)(using conn: Connection): List[StrProfileEntity] =
    queryList(
      "SELECT * FROM str_profile WHERE imported_from = ? ORDER BY created_at DESC",
      Seq(provider)
    )(mapRow)

  /**
   * Find profiles pending sync.
   */
  def findPendingSync()(using conn: Connection): List[StrProfileEntity] =
    queryList(
      """SELECT * FROM str_profile
        |WHERE sync_status IN ('Local', 'Modified')
        |ORDER BY updated_at ASC
      """.stripMargin
    )(mapRow)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): StrProfileEntity =
    val panelsJson = getOptJsonString(rs, "panels")
    val panels = panelsJson.flatMap(json => parse(json).flatMap(_.as[List[StrPanel]]).toOption).getOrElse(List.empty)

    val markersJson = getOptJsonString(rs, "markers")
    val markers = markersJson.flatMap(json => parse(json).flatMap(_.as[List[StrMarkerValue]]).toOption).getOrElse(List.empty)

    val filesJson = getOptJsonString(rs, "files")
    val files = filesJson.flatMap(json => parse(json).flatMap(_.as[List[FileInfo]]).toOption).getOrElse(List.empty)

    StrProfileEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      sequenceRunId = getOptUUID(rs, "sequence_run_id"),
      panels = panels,
      markers = markers,
      totalMarkers = getOptInt(rs, "total_markers"),
      source = getOptString(rs, "source"),
      importedFrom = getOptString(rs, "imported_from"),
      derivationMethod = getOptString(rs, "derivation_method"),
      files = files,
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
