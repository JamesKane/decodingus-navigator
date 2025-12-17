package com.decodingus.yprofile.repository

import com.decodingus.repository.{SyncableRepository, Entity, EntityMeta, SqlHelpers, SyncStatus}
import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.FileInfo
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.{LocalDateTime, Instant}
import java.util.UUID

/**
 * Variant type for Y-DNA markers.
 */
enum YVariantType:
  case SNP, INDEL

/**
 * A single SNP or INDEL call from a Y-SNP panel test.
 *
 * Supports both SNPs (single position) and INDELs (position range).
 *
 * @param name          Marker name (e.g., "M343", "L21", "A1133")
 * @param startPosition GRCh38 start position on chrY
 * @param endPosition   GRCh38 end position (same as start for SNPs, different for INDELs)
 * @param allele        Called allele (e.g., "A", "G", "ins", "del") without +/- suffix
 * @param derived       Whether this is the derived (mutant) state (+ = true, - = false)
 * @param variantType   Type of variant: SNP or INDEL
 * @param orderedDate   Date this specific marker was ordered/tested (for incremental panels)
 * @param quality       Quality score if available
 */
case class YSnpCall(
  name: String,
  startPosition: Long,
  endPosition: Option[Long] = None,
  allele: String,
  derived: Boolean,
  variantType: Option[YVariantType] = None,
  orderedDate: Option[LocalDateTime] = None,
  quality: Option[Double] = None
) {
  /** Convenience: get end position, defaulting to start for SNPs */
  def effectiveEndPosition: Long = endPosition.getOrElse(startPosition)

  /** Check if this is an INDEL (has a range) */
  def isIndel: Boolean = variantType.contains(YVariantType.INDEL) || endPosition.exists(_ != startPosition)
}

/**
 * A private/novel Y-DNA variant not in the reference tree.
 *
 * @param position    GRCh38 position on chrY
 * @param refAllele   Reference allele
 * @param altAllele   Alternate (called) allele
 * @param snpName     Assigned name if available
 * @param quality     Quality score if available
 * @param readDepth   Read depth at this position
 */
case class YPrivateVariant(
  position: Long,
  refAllele: String,
  altAllele: String,
  snpName: Option[String] = None,
  quality: Option[Double] = None,
  readDepth: Option[Int] = None
)

/**
 * Y-SNP Panel entity for database persistence.
 *
 * Stores Y-DNA SNP panel testing results from various providers.
 */
case class YSnpPanelEntity(
  id: UUID,
  biosampleId: UUID,
  alignmentId: Option[UUID],
  panelName: Option[String],
  provider: Option[String],
  testDate: Option[LocalDateTime],
  totalSnpsTested: Option[Int],
  derivedCount: Option[Int],
  ancestralCount: Option[Int],
  noCallCount: Option[Int],
  terminalHaplogroup: Option[String],
  confidence: Option[Double],
  snpCalls: List[YSnpCall],
  privateVariants: List[YPrivateVariant],
  files: List[FileInfo],
  meta: EntityMeta
) extends Entity[UUID]

object YSnpPanelEntity:
  import YSnpPanelCodecs.given

  /**
   * Create a new YSnpPanelEntity with generated ID and initial metadata.
   */
  def create(
    biosampleId: UUID,
    alignmentId: Option[UUID] = None,
    panelName: Option[String] = None,
    provider: Option[String] = None,
    testDate: Option[LocalDateTime] = None,
    totalSnpsTested: Option[Int] = None,
    derivedCount: Option[Int] = None,
    ancestralCount: Option[Int] = None,
    noCallCount: Option[Int] = None,
    terminalHaplogroup: Option[String] = None,
    confidence: Option[Double] = None,
    snpCalls: List[YSnpCall] = List.empty,
    privateVariants: List[YPrivateVariant] = List.empty,
    files: List[FileInfo] = List.empty
  ): YSnpPanelEntity = YSnpPanelEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    alignmentId = alignmentId,
    panelName = panelName,
    provider = provider,
    testDate = testDate,
    totalSnpsTested = totalSnpsTested,
    derivedCount = derivedCount,
    ancestralCount = ancestralCount,
    noCallCount = noCallCount,
    terminalHaplogroup = terminalHaplogroup,
    confidence = confidence,
    snpCalls = snpCalls,
    privateVariants = privateVariants,
    files = files,
    meta = EntityMeta.create()
  )

/**
 * Circe codecs for Y-SNP panel JSON fields.
 */
object YSnpPanelCodecs:
  // YVariantType enum codec
  given Encoder[YVariantType] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YVariantType] = Decoder.decodeString.emap { s =>
    try Right(YVariantType.valueOf(s))
    catch case _: IllegalArgumentException => Left(s"Invalid YVariantType: $s")
  }

  // YSnpCall with backwards compatibility for old "position" field
  given Encoder[YSnpCall] = Encoder.instance { call =>
    Json.obj(
      "name" -> Json.fromString(call.name),
      "startPosition" -> Json.fromLong(call.startPosition),
      "endPosition" -> call.endPosition.fold(Json.Null)(Json.fromLong),
      "allele" -> Json.fromString(call.allele),
      "derived" -> Json.fromBoolean(call.derived),
      "variantType" -> call.variantType.fold(Json.Null)(vt => Json.fromString(vt.toString)),
      "orderedDate" -> call.orderedDate.fold(Json.Null)(dt => Json.fromString(dt.toString)),
      "quality" -> call.quality.fold(Json.Null)(Json.fromDoubleOrNull)
    ).dropNullValues
  }

  given Decoder[YSnpCall] = Decoder.instance { cursor =>
    for
      name <- cursor.get[String]("name")
      // Support both old "position" and new "startPosition" field names
      startPosition <- cursor.get[Long]("startPosition").orElse(cursor.get[Long]("position"))
      endPosition <- cursor.get[Option[Long]]("endPosition")
      allele <- cursor.get[String]("allele")
      derived <- cursor.get[Boolean]("derived")
      variantType <- cursor.get[Option[YVariantType]]("variantType")
      orderedDate <- cursor.get[Option[LocalDateTime]]("orderedDate")
      quality <- cursor.get[Option[Double]]("quality")
    yield YSnpCall(name, startPosition, endPosition, allele, derived, variantType, orderedDate, quality)
  }

  given Encoder[YPrivateVariant] = deriveEncoder
  given Decoder[YPrivateVariant] = deriveDecoder
  given Encoder[FileInfo] = deriveEncoder
  given Decoder[FileInfo] = deriveDecoder

/**
 * Repository for Y-SNP panel persistence operations.
 */
class YSnpPanelRepository extends SyncableRepository[YSnpPanelEntity, UUID]:

  import YSnpPanelCodecs.given

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[YSnpPanelEntity] =
    queryOne(
      "SELECT * FROM y_snp_panel WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[YSnpPanelEntity] =
    queryList("SELECT * FROM y_snp_panel ORDER BY created_at DESC")(mapRow)

  override def insert(entity: YSnpPanelEntity)(using conn: Connection): YSnpPanelEntity =
    val snpCallsJson = JsonValue(entity.snpCalls.asJson.noSpaces)
    val privateVariantsJson = JsonValue(entity.privateVariants.asJson.noSpaces)
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO y_snp_panel (
        |  id, biosample_id, alignment_id, panel_name, provider, test_date,
        |  total_snps_tested, derived_count, ancestral_count, no_call_count,
        |  terminal_haplogroup, confidence, snp_calls, private_variants, files,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.alignmentId,
        entity.panelName,
        entity.provider,
        entity.testDate,
        entity.totalSnpsTested,
        entity.derivedCount,
        entity.ancestralCount,
        entity.noCallCount,
        entity.terminalHaplogroup,
        entity.confidence,
        snpCallsJson,
        privateVariantsJson,
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

  override def update(entity: YSnpPanelEntity)(using conn: Connection): YSnpPanelEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val snpCallsJson = JsonValue(entity.snpCalls.asJson.noSpaces)
    val privateVariantsJson = JsonValue(entity.privateVariants.asJson.noSpaces)
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """UPDATE y_snp_panel SET
        |  biosample_id = ?, alignment_id = ?, panel_name = ?, provider = ?, test_date = ?,
        |  total_snps_tested = ?, derived_count = ?, ancestral_count = ?, no_call_count = ?,
        |  terminal_haplogroup = ?, confidence = ?, snp_calls = ?, private_variants = ?, files = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.alignmentId,
        entity.panelName,
        entity.provider,
        entity.testDate,
        entity.totalSnpsTested,
        entity.derivedCount,
        entity.ancestralCount,
        entity.noCallCount,
        entity.terminalHaplogroup,
        entity.confidence,
        snpCallsJson,
        privateVariantsJson,
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
    executeUpdate("DELETE FROM y_snp_panel WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM y_snp_panel") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM y_snp_panel WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      "SELECT * FROM y_snp_panel WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE y_snp_panel SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE y_snp_panel SET
        |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // Y-SNP Panel-Specific Queries
  // ============================================

  /**
   * Find all Y-SNP panels for a biosample.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      "SELECT * FROM y_snp_panel WHERE biosample_id = ? ORDER BY test_date DESC NULLS LAST",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Find Y-SNP panels derived from an alignment.
   */
  def findByAlignment(alignmentId: UUID)(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      "SELECT * FROM y_snp_panel WHERE alignment_id = ? ORDER BY created_at DESC",
      Seq(alignmentId)
    )(mapRow)

  /**
   * Find Y-SNP panels by provider.
   */
  def findByProvider(provider: String)(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      "SELECT * FROM y_snp_panel WHERE provider = ? ORDER BY test_date DESC NULLS LAST",
      Seq(provider)
    )(mapRow)

  /**
   * Find Y-SNP panels by terminal haplogroup.
   */
  def findByHaplogroup(haplogroup: String)(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      "SELECT * FROM y_snp_panel WHERE terminal_haplogroup = ? ORDER BY created_at DESC",
      Seq(haplogroup)
    )(mapRow)

  /**
   * Find Y-SNP panels with haplogroups under a branch (prefix match).
   */
  def findByHaplogroupBranch(branchPrefix: String)(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      "SELECT * FROM y_snp_panel WHERE terminal_haplogroup LIKE ? ORDER BY terminal_haplogroup",
      Seq(s"$branchPrefix%")
    )(mapRow)

  /**
   * Find profiles pending sync.
   */
  def findPendingSync()(using conn: Connection): List[YSnpPanelEntity] =
    queryList(
      """SELECT * FROM y_snp_panel
        |WHERE sync_status IN ('Local', 'Modified')
        |ORDER BY updated_at ASC
      """.stripMargin
    )(mapRow)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): YSnpPanelEntity =
    val snpCallsJson = getOptJsonString(rs, "snp_calls")
    val snpCalls = snpCallsJson.flatMap(json =>
      parse(json).flatMap(_.as[List[YSnpCall]]).toOption
    ).getOrElse(List.empty)

    val privateVariantsJson = getOptJsonString(rs, "private_variants")
    val privateVariants = privateVariantsJson.flatMap(json =>
      parse(json).flatMap(_.as[List[YPrivateVariant]]).toOption
    ).getOrElse(List.empty)

    val filesJson = getOptJsonString(rs, "files")
    val files = filesJson.flatMap(json =>
      parse(json).flatMap(_.as[List[FileInfo]]).toOption
    ).getOrElse(List.empty)

    val testDateTs = rs.getTimestamp("test_date")
    val testDate = if rs.wasNull() then None else Some(testDateTs.toLocalDateTime)

    YSnpPanelEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      alignmentId = getOptUUID(rs, "alignment_id"),
      panelName = getOptString(rs, "panel_name"),
      provider = getOptString(rs, "provider"),
      testDate = testDate,
      totalSnpsTested = getOptInt(rs, "total_snps_tested"),
      derivedCount = getOptInt(rs, "derived_count"),
      ancestralCount = getOptInt(rs, "ancestral_count"),
      noCallCount = getOptInt(rs, "no_call_count"),
      terminalHaplogroup = getOptString(rs, "terminal_haplogroup"),
      confidence = getOptDouble(rs, "confidence"),
      snpCalls = snpCalls,
      privateVariants = privateVariants,
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
