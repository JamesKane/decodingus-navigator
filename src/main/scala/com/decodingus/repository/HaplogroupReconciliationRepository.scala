package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.{AuditEntry, CallMethod, CompatibilityLevel, ConflictResolution, DnaType, HaplogroupReconciliation, HaplogroupTechnology, HeteroplasmyObservation, IdentityVerification, ManualOverride, ReconciliationStatus, RunHaplogroupCall, SnpCallFromRun, SnpConflict}
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
                                           heteroplasmyObservations: List[HeteroplasmyObservation] = List.empty,
                                           identityVerification: Option[IdentityVerification] = None,
                                           manualOverride: Option[ManualOverride] = None,
                                           auditLog: List[AuditEntry] = List.empty,
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
              heteroplasmyObservations: List[HeteroplasmyObservation] = List.empty,
              identityVerification: Option[IdentityVerification] = None,
              manualOverride: Option[ManualOverride] = None,
              auditLog: List[AuditEntry] = List.empty,
              lastReconciliationAt: Option[Instant] = None
            ): HaplogroupReconciliationEntity = HaplogroupReconciliationEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    dnaType = dnaType,
    status = status,
    runCalls = runCalls,
    snpConflicts = snpConflicts,
    heteroplasmyObservations = heteroplasmyObservations,
    identityVerification = identityVerification,
    manualOverride = manualOverride,
    auditLog = auditLog,
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
      heteroplasmyObservations = reconciliation.heteroplasmyObservations,
      identityVerification = reconciliation.identityVerification,
      manualOverride = reconciliation.manualOverride,
      auditLog = reconciliation.auditLog,
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

  given Encoder[HeteroplasmyObservation] = Encoder.instance { h =>
    Json.obj(
      "position" -> Json.fromInt(h.position),
      "majorAllele" -> Json.fromString(h.majorAllele),
      "minorAllele" -> Json.fromString(h.minorAllele),
      "majorAlleleFrequency" -> Json.fromDoubleOrNull(h.majorAlleleFrequency),
      "depth" -> h.depth.fold(Json.Null)(Json.fromInt),
      "isDefiningSnp" -> h.isDefiningSnp.fold(Json.Null)(Json.fromBoolean),
      "affectedHaplogroup" -> h.affectedHaplogroup.fold(Json.Null)(Json.fromString)
    )
  }

  given Decoder[HeteroplasmyObservation] = Decoder.instance { c =>
    for
      pos <- c.get[Int]("position")
      maj <- c.get[String]("majorAllele")
      min <- c.get[String]("minorAllele")
      freq <- c.get[Double]("majorAlleleFrequency")
      depth <- c.get[Option[Int]]("depth")
      isDef <- c.get[Option[Boolean]]("isDefiningSnp")
      haplo <- c.get[Option[String]]("affectedHaplogroup")
    yield HeteroplasmyObservation(pos, maj, min, freq, depth, isDef, haplo)
  }

  given Encoder[IdentityVerification] = Encoder.instance { iv =>
    Json.obj(
      "kinshipCoefficient" -> iv.kinshipCoefficient.fold(Json.Null)(Json.fromDoubleOrNull),
      "fingerprintSnpConcordance" -> iv.fingerprintSnpConcordance.fold(Json.Null)(Json.fromDoubleOrNull),
      "yStrDistance" -> iv.yStrDistance.fold(Json.Null)(Json.fromInt),
      "verificationStatus" -> iv.verificationStatus.fold(Json.Null)(Json.fromString),
      "verificationMethod" -> iv.verificationMethod.fold(Json.Null)(Json.fromString)
    )
  }

  given Decoder[IdentityVerification] = Decoder.instance { c =>
    for
      kinship <- c.get[Option[Double]]("kinshipCoefficient")
      fingerprint <- c.get[Option[Double]]("fingerprintSnpConcordance")
      yStr <- c.get[Option[Int]]("yStrDistance")
      status <- c.get[Option[String]]("verificationStatus")
      method <- c.get[Option[String]]("verificationMethod")
    yield IdentityVerification(kinship, fingerprint, yStr, status, method)
  }

  given Encoder[ManualOverride] = Encoder.instance { mo =>
    Json.obj(
      "overriddenHaplogroup" -> Json.fromString(mo.overriddenHaplogroup),
      "reason" -> mo.reason.fold(Json.Null)(Json.fromString),
      "overriddenAt" -> mo.overriddenAt.fold(Json.Null)(t => Json.fromString(t.toString)),
      "overriddenBy" -> mo.overriddenBy.fold(Json.Null)(Json.fromString)
    )
  }

  given Decoder[ManualOverride] = Decoder.instance { c =>
    for
      haplo <- c.get[String]("overriddenHaplogroup")
      reason <- c.get[Option[String]]("reason")
      at <- c.get[Option[String]]("overriddenAt").map(_.flatMap(s =>
        try Some(LocalDateTime.parse(s)) catch case _: Exception => None
      ))
      by <- c.get[Option[String]]("overriddenBy")
    yield ManualOverride(haplo, reason, at, by)
  }

  given Encoder[AuditEntry] = Encoder.instance { ae =>
    Json.obj(
      "timestamp" -> Json.fromString(ae.timestamp.toString),
      "action" -> Json.fromString(ae.action),
      "previousConsensus" -> ae.previousConsensus.fold(Json.Null)(Json.fromString),
      "newConsensus" -> ae.newConsensus.fold(Json.Null)(Json.fromString),
      "runRef" -> ae.runRef.fold(Json.Null)(Json.fromString),
      "notes" -> ae.notes.fold(Json.Null)(Json.fromString)
    )
  }

  given Decoder[AuditEntry] = Decoder.instance { c =>
    for
      tsStr <- c.get[String]("timestamp")
      ts = try LocalDateTime.parse(tsStr) catch case _: Exception => LocalDateTime.MIN
      action <- c.get[String]("action")
      prev <- c.get[Option[String]]("previousConsensus")
      next <- c.get[Option[String]]("newConsensus")
      runRef <- c.get[Option[String]]("runRef")
      notes <- c.get[Option[String]]("notes")
    yield AuditEntry(ts, action, prev, next, runRef, notes)
  }

/**
 * Repository for Haplogroup reconciliation persistence operations.
 */
class HaplogroupReconciliationRepository extends SyncableRepositoryBase[HaplogroupReconciliationEntity]:

  import HaplogroupReconciliationCodecs.given

  override protected def tableName: String = "haplogroup_reconciliation"

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
    val heteroplasmyJson = JsonValue(entity.heteroplasmyObservations.asJson.noSpaces)
    val identityJson = entity.identityVerification.map(iv => JsonValue(iv.asJson.noSpaces))
    val overrideJson = entity.manualOverride.map(mo => JsonValue(mo.asJson.noSpaces))
    val auditLogJson = JsonValue(entity.auditLog.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO haplogroup_reconciliation (
        |  id, biosample_id, dna_type, status, run_calls, snp_conflicts,
        |  heteroplasmy_observations, identity_verification, manual_override, audit_log,
        |  last_reconciliation_at, sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.dnaType.toString,
        statusJson,
        runCallsJson,
        snpConflictsJson,
        heteroplasmyJson,
        identityJson,
        overrideJson,
        auditLogJson,
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
    val heteroplasmyJson = JsonValue(entity.heteroplasmyObservations.asJson.noSpaces)
    val identityJson = entity.identityVerification.map(iv => JsonValue(iv.asJson.noSpaces))
    val overrideJson = entity.manualOverride.map(mo => JsonValue(mo.asJson.noSpaces))
    val auditLogJson = JsonValue(entity.auditLog.asJson.noSpaces)

    executeUpdate(
      """UPDATE haplogroup_reconciliation SET
        |  biosample_id = ?, dna_type = ?, status = ?, run_calls = ?, snp_conflicts = ?,
        |  heteroplasmy_observations = ?, identity_verification = ?, manual_override = ?, audit_log = ?,
        |  last_reconciliation_at = ?, sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.dnaType.toString,
        statusJson,
        runCallsJson,
        snpConflictsJson,
        heteroplasmyJson,
        identityJson,
        overrideJson,
        auditLogJson,
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

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM haplogroup_reconciliation WHERE id = ?", Seq(id)) { _ => true }.isDefined

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

  override protected def mapRow(rs: ResultSet): HaplogroupReconciliationEntity =
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

    val heteroplasmyJson = getOptJsonString(rs, "heteroplasmy_observations")
    val heteroplasmyObs = heteroplasmyJson.flatMap(json =>
      parse(json).flatMap(_.as[List[HeteroplasmyObservation]]).toOption
    ).getOrElse(List.empty)

    val identityJson = getOptJsonString(rs, "identity_verification")
    val identityVerification = identityJson.flatMap(json =>
      parse(json).flatMap(_.as[IdentityVerification]).toOption
    )

    val overrideJson = getOptJsonString(rs, "manual_override")
    val manualOverride = overrideJson.flatMap(json =>
      parse(json).flatMap(_.as[ManualOverride]).toOption
    )

    val auditLogJson = getOptJsonString(rs, "audit_log")
    val auditLog = auditLogJson.flatMap(json =>
      parse(json).flatMap(_.as[List[AuditEntry]]).toOption
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
      heteroplasmyObservations = heteroplasmyObs,
      identityVerification = identityVerification,
      manualOverride = manualOverride,
      auditLog = auditLog,
      lastReconciliationAt = lastReconciliationAt,
      meta = readEntityMeta(rs)
    )
