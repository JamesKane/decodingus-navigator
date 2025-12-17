package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.{CallMethod, CompatibilityLevel, ConflictResolution, DnaType, HaplogroupReconciliation, HaplogroupTechnology, ReconciliationStatus, RunHaplogroupCall, SnpCallFromRun, SnpConflict}
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*

import java.sql.{Connection, ResultSet}
import java.time.{Instant, LocalDateTime}
import java.util.UUID

/**
 * Haplogroup Reconciliation entity for database persistence.
 *
 * Stores reconciliation of haplogroup calls across multiple runs for a biosample.
 */
case class HaplogroupReconciliationEntity(
                                           id: UUID,
                                           biosampleId: UUID,
                                           dnaType: DnaType,
                                           status: ReconciliationStatus,
                                           runCalls: List[RunHaplogroupCall],
                                           snpConflicts: List[SnpConflict],
                                           lastReconciliationAt: Option[Instant],
                                           meta: EntityMeta
                                         ) extends Entity[UUID]

object HaplogroupReconciliationEntity:

  import HaplogroupReconciliationCodecs.given

  /**
   * Create a new HaplogroupReconciliationEntity with generated ID and initial metadata.
   */
  def create(
              biosampleId: UUID,
              dnaType: DnaType,
              status: ReconciliationStatus,
              runCalls: List[RunHaplogroupCall] = List.empty,
              snpConflicts: List[SnpConflict] = List.empty,
              lastReconciliationAt: Option[Instant] = None
            ): HaplogroupReconciliationEntity = HaplogroupReconciliationEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    dnaType = dnaType,
    status = status,
    runCalls = runCalls,
    snpConflicts = snpConflicts,
    lastReconciliationAt = lastReconciliationAt,
    meta = EntityMeta.create()
  )

  /**
   * Create entity from a workspace model HaplogroupReconciliation.
   */
  def fromModel(biosampleId: UUID, reconciliation: HaplogroupReconciliation): HaplogroupReconciliationEntity =
    HaplogroupReconciliationEntity(
      id = UUID.randomUUID(),
      biosampleId = biosampleId,
      dnaType = reconciliation.dnaType,
      status = reconciliation.status,
      runCalls = reconciliation.runCalls,
      snpConflicts = reconciliation.snpConflicts,
      lastReconciliationAt = reconciliation.lastReconciliationAt,
      meta = EntityMeta.create()
    )

/**
 * Circe codecs for Haplogroup reconciliation JSON fields.
 */
object HaplogroupReconciliationCodecs:
  // Enum codecs
  given Encoder[CompatibilityLevel] = Encoder.encodeString.contramap(_.toString)

  given Decoder[CompatibilityLevel] = Decoder.decodeString.emap { s =>
    try Right(CompatibilityLevel.valueOf(s))
    catch case _: IllegalArgumentException => Left(s"Invalid CompatibilityLevel: $s")
  }

  given Encoder[HaplogroupTechnology] = Encoder.encodeString.contramap(_.toString)

  given Decoder[HaplogroupTechnology] = Decoder.decodeString.emap { s =>
    try Right(HaplogroupTechnology.valueOf(s))
    catch case _: IllegalArgumentException => Left(s"Invalid HaplogroupTechnology: $s")
  }

  given Encoder[CallMethod] = Encoder.encodeString.contramap(_.toString)

  given Decoder[CallMethod] = Decoder.decodeString.emap { s =>
    try Right(CallMethod.valueOf(s))
    catch case _: IllegalArgumentException => Left(s"Invalid CallMethod: $s")
  }

  given Encoder[ConflictResolution] = Encoder.encodeString.contramap(_.toString)

  given Decoder[ConflictResolution] = Decoder.decodeString.emap { s =>
    try Right(ConflictResolution.valueOf(s))
    catch case _: IllegalArgumentException => Left(s"Invalid ConflictResolution: $s")
  }

  // Complex type codecs
  given Encoder[ReconciliationStatus] = deriveEncoder

  given Decoder[ReconciliationStatus] = deriveDecoder

  given Encoder[RunHaplogroupCall] = deriveEncoder

  given Decoder[RunHaplogroupCall] = deriveDecoder

  given Encoder[SnpCallFromRun] = deriveEncoder

  given Decoder[SnpCallFromRun] = deriveDecoder

  given Encoder[SnpConflict] = deriveEncoder

  given Decoder[SnpConflict] = deriveDecoder

/**
 * Repository for Haplogroup reconciliation persistence operations.
 */
class HaplogroupReconciliationRepository extends SyncableRepository[HaplogroupReconciliationEntity, UUID]:

  import HaplogroupReconciliationCodecs.given

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[HaplogroupReconciliationEntity] =
    queryOne(
      "SELECT * FROM haplogroup_reconciliation WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[HaplogroupReconciliationEntity] =
    queryList("SELECT * FROM haplogroup_reconciliation ORDER BY created_at DESC")(mapRow)

  override def insert(entity: HaplogroupReconciliationEntity)(using conn: Connection): HaplogroupReconciliationEntity =
    val statusJson = JsonValue(entity.status.asJson.noSpaces)
    val runCallsJson = JsonValue(entity.runCalls.asJson.noSpaces)
    val snpConflictsJson = JsonValue(entity.snpConflicts.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO haplogroup_reconciliation (
        |  id, biosample_id, dna_type, status, run_calls, snp_conflicts,
        |  last_reconciliation_at, sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.dnaType.toString,
        statusJson,
        runCallsJson,
        snpConflictsJson,
        entity.lastReconciliationAt.map(i => java.sql.Timestamp.from(i)),
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: HaplogroupReconciliationEntity)(using conn: Connection): HaplogroupReconciliationEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val statusJson = JsonValue(entity.status.asJson.noSpaces)
    val runCallsJson = JsonValue(entity.runCalls.asJson.noSpaces)
    val snpConflictsJson = JsonValue(entity.snpConflicts.asJson.noSpaces)

    executeUpdate(
      """UPDATE haplogroup_reconciliation SET
        |  biosample_id = ?, dna_type = ?, status = ?, run_calls = ?, snp_conflicts = ?,
        |  last_reconciliation_at = ?, sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.dnaType.toString,
        statusJson,
        runCallsJson,
        snpConflictsJson,
        entity.lastReconciliationAt.map(i => java.sql.Timestamp.from(i)),
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
    executeUpdate("DELETE FROM haplogroup_reconciliation WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM haplogroup_reconciliation") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM haplogroup_reconciliation WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[HaplogroupReconciliationEntity] =
    queryList(
      "SELECT * FROM haplogroup_reconciliation WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE haplogroup_reconciliation SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE haplogroup_reconciliation SET
        |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // Haplogroup Reconciliation-Specific Queries
  // ============================================

  /**
   * Find reconciliation for a biosample and DNA type (Y or MT).
   */
  def findByBiosampleAndDnaType(biosampleId: UUID, dnaType: DnaType)(using conn: Connection): Option[HaplogroupReconciliationEntity] =
    queryOne(
      "SELECT * FROM haplogroup_reconciliation WHERE biosample_id = ? AND dna_type = ?",
      Seq(biosampleId, dnaType.toString)
    )(mapRow)

  /**
   * Find all reconciliations for a biosample.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[HaplogroupReconciliationEntity] =
    queryList(
      "SELECT * FROM haplogroup_reconciliation WHERE biosample_id = ? ORDER BY dna_type",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Find reconciliations by DNA type.
   */
  def findByDnaType(dnaType: DnaType)(using conn: Connection): List[HaplogroupReconciliationEntity] =
    queryList(
      "SELECT * FROM haplogroup_reconciliation WHERE dna_type = ? ORDER BY created_at DESC",
      Seq(dnaType.toString)
    )(mapRow)

  /**
   * Find reconciliations with conflicts (divergence or incompatibility).
   */
  def findWithConflicts()(using conn: Connection): List[HaplogroupReconciliationEntity] =
    // Query reconciliations where status JSON contains divergence or incompatible
    queryList(
      """SELECT * FROM haplogroup_reconciliation
        |WHERE status LIKE '%MAJOR_DIVERGENCE%' OR status LIKE '%INCOMPATIBLE%'
        |ORDER BY updated_at DESC
      """.stripMargin
    )(mapRow)

  /**
   * Find profiles pending sync.
   */
  def findPendingSync()(using conn: Connection): List[HaplogroupReconciliationEntity] =
    queryList(
      """SELECT * FROM haplogroup_reconciliation
        |WHERE sync_status IN ('Local', 'Modified')
        |ORDER BY updated_at ASC
      """.stripMargin
    )(mapRow)

  /**
   * Upsert reconciliation (insert if not exists, update if exists).
   * Uses the unique constraint on (biosample_id, dna_type).
   */
  def upsert(entity: HaplogroupReconciliationEntity)(using conn: Connection): HaplogroupReconciliationEntity =
    findByBiosampleAndDnaType(entity.biosampleId, entity.dnaType) match
      case Some(existing) =>
        update(entity.copy(id = existing.id, meta = existing.meta))
      case None =>
        insert(entity)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): HaplogroupReconciliationEntity =
    val statusJson = getOptJsonString(rs, "status")
    val status = statusJson.flatMap(json =>
      parse(json).flatMap(_.as[ReconciliationStatus]).toOption
    ).getOrElse(ReconciliationStatus(
      compatibilityLevel = CompatibilityLevel.COMPATIBLE,
      consensusHaplogroup = "",
      confidence = 0.0,
      runCount = 0
    ))

    val runCallsJson = getOptJsonString(rs, "run_calls")
    val runCalls = runCallsJson.flatMap(json =>
      parse(json).flatMap(_.as[List[RunHaplogroupCall]]).toOption
    ).getOrElse(List.empty)

    val snpConflictsJson = getOptJsonString(rs, "snp_conflicts")
    val snpConflicts = snpConflictsJson.flatMap(json =>
      parse(json).flatMap(_.as[List[SnpConflict]]).toOption
    ).getOrElse(List.empty)

    val lastReconTs = rs.getTimestamp("last_reconciliation_at")
    val lastReconciliationAt = if rs.wasNull() then None else Some(lastReconTs.toInstant)

    val dnaTypeStr = rs.getString("dna_type")
    val dnaType = try DnaType.valueOf(dnaTypeStr) catch case _: IllegalArgumentException => DnaType.Y_DNA

    HaplogroupReconciliationEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      dnaType = dnaType,
      status = status,
      runCalls = runCalls,
      snpConflicts = snpConflicts,
      lastReconciliationAt = lastReconciliationAt,
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
