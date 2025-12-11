package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import io.circe.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Suggested conflict resolution strategy.
 */
enum ConflictResolution:
  case KeepLocal
  case AcceptRemote
  case Merge
  case Manual

object ConflictResolution:
  def fromString(s: String): ConflictResolution = s match
    case "KEEP_LOCAL" => KeepLocal
    case "ACCEPT_REMOTE" => AcceptRemote
    case "MERGE" => Merge
    case "MANUAL" => Manual
    case other => throw new IllegalArgumentException(s"Unknown resolution: $other")

  def toDbString(r: ConflictResolution): String = r match
    case KeepLocal => "KEEP_LOCAL"
    case AcceptRemote => "ACCEPT_REMOTE"
    case Merge => "MERGE"
    case Manual => "MANUAL"

/**
 * Conflict status.
 */
enum ConflictStatus:
  case Unresolved
  case Resolved
  case Dismissed

object ConflictStatus:
  def fromString(s: String): ConflictStatus = s match
    case "UNRESOLVED" => Unresolved
    case "RESOLVED" => Resolved
    case "DISMISSED" => Dismissed
    case other => throw new IllegalArgumentException(s"Unknown conflict status: $other")

  def toDbString(s: ConflictStatus): String = s match
    case Unresolved => "UNRESOLVED"
    case Resolved => "RESOLVED"
    case Dismissed => "DISMISSED"

/**
 * Resolution action taken.
 */
enum ResolutionAction:
  case KeptLocal
  case AcceptedRemote
  case Merged
  case ManualEdit

object ResolutionAction:
  def fromString(s: String): ResolutionAction = s match
    case "KEPT_LOCAL" => KeptLocal
    case "ACCEPTED_REMOTE" => AcceptedRemote
    case "MERGED" => Merged
    case "MANUAL_EDIT" => ManualEdit
    case other => throw new IllegalArgumentException(s"Unknown action: $other")

  def toDbString(a: ResolutionAction): String = a match
    case KeptLocal => "KEPT_LOCAL"
    case AcceptedRemote => "ACCEPTED_REMOTE"
    case Merged => "MERGED"
    case ManualEdit => "MANUAL_EDIT"

/**
 * Sync conflict entity.
 */
case class SyncConflictEntity(
  id: UUID,
  entityType: SyncEntityType,
  entityId: UUID,
  atUri: Option[String],
  detectedAt: LocalDateTime,
  localVersion: Int,
  remoteVersion: Int,
  localChanges: Option[Json],
  remoteChanges: Option[Json],
  overlappingFields: Option[Json],
  suggestedResolution: Option[ConflictResolution],
  resolutionReason: Option[String],
  status: ConflictStatus,
  resolvedAt: Option[LocalDateTime],
  resolvedBy: Option[String],
  resolutionAction: Option[ResolutionAction],
  localSnapshot: Option[Json],
  remoteSnapshot: Option[Json],
  createdAt: LocalDateTime,
  updatedAt: LocalDateTime
) extends Entity[UUID]

object SyncConflictEntity:
  /**
   * Create a new conflict entry.
   */
  def create(
    entityType: SyncEntityType,
    entityId: UUID,
    localVersion: Int,
    remoteVersion: Int,
    atUri: Option[String] = None,
    localChanges: Option[Json] = None,
    remoteChanges: Option[Json] = None,
    overlappingFields: Option[Json] = None,
    suggestedResolution: Option[ConflictResolution] = None,
    resolutionReason: Option[String] = None,
    localSnapshot: Option[Json] = None,
    remoteSnapshot: Option[Json] = None
  ): SyncConflictEntity =
    val now = LocalDateTime.now()
    SyncConflictEntity(
      id = UUID.randomUUID(),
      entityType = entityType,
      entityId = entityId,
      atUri = atUri,
      detectedAt = now,
      localVersion = localVersion,
      remoteVersion = remoteVersion,
      localChanges = localChanges,
      remoteChanges = remoteChanges,
      overlappingFields = overlappingFields,
      suggestedResolution = suggestedResolution,
      resolutionReason = resolutionReason,
      status = ConflictStatus.Unresolved,
      resolvedAt = None,
      resolvedBy = None,
      resolutionAction = None,
      localSnapshot = localSnapshot,
      remoteSnapshot = remoteSnapshot,
      createdAt = now,
      updatedAt = now
    )

/**
 * Repository for sync conflicts.
 */
class SyncConflictRepository:

  // ============================================
  // Core Operations
  // ============================================

  def findById(id: UUID)(using conn: Connection): Option[SyncConflictEntity] =
    queryOne(
      "SELECT * FROM sync_conflict WHERE id = ?",
      Seq(id)
    )(mapRow)

  def findAll()(using conn: Connection): List[SyncConflictEntity] =
    queryList("SELECT * FROM sync_conflict ORDER BY detected_at DESC")(mapRow)

  def insert(entity: SyncConflictEntity)(using conn: Connection): SyncConflictEntity =
    executeUpdate(
      """INSERT INTO sync_conflict (
        |  id, entity_type, entity_id, at_uri, detected_at, local_version, remote_version,
        |  local_changes, remote_changes, overlapping_fields, suggested_resolution, resolution_reason,
        |  status, resolved_at, resolved_by, resolution_action, local_snapshot, remote_snapshot,
        |  created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        SyncEntityType.toDbString(entity.entityType),
        entity.entityId,
        entity.atUri,
        entity.detectedAt,
        entity.localVersion,
        entity.remoteVersion,
        entity.localChanges.map(j => JsonValue(j.noSpaces)),
        entity.remoteChanges.map(j => JsonValue(j.noSpaces)),
        entity.overlappingFields.map(j => JsonValue(j.noSpaces)),
        entity.suggestedResolution.map(ConflictResolution.toDbString).orNull,
        entity.resolutionReason,
        ConflictStatus.toDbString(entity.status),
        entity.resolvedAt,
        entity.resolvedBy,
        entity.resolutionAction.map(ResolutionAction.toDbString).orNull,
        entity.localSnapshot.map(j => JsonValue(j.noSpaces)),
        entity.remoteSnapshot.map(j => JsonValue(j.noSpaces)),
        entity.createdAt,
        entity.updatedAt
      )
    )
    entity

  def update(entity: SyncConflictEntity)(using conn: Connection): SyncConflictEntity =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_conflict SET
        |  entity_type = ?, entity_id = ?, at_uri = ?, detected_at = ?,
        |  local_version = ?, remote_version = ?, local_changes = ?, remote_changes = ?,
        |  overlapping_fields = ?, suggested_resolution = ?, resolution_reason = ?,
        |  status = ?, resolved_at = ?, resolved_by = ?, resolution_action = ?,
        |  local_snapshot = ?, remote_snapshot = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        SyncEntityType.toDbString(entity.entityType),
        entity.entityId,
        entity.atUri,
        entity.detectedAt,
        entity.localVersion,
        entity.remoteVersion,
        entity.localChanges.map(j => JsonValue(j.noSpaces)),
        entity.remoteChanges.map(j => JsonValue(j.noSpaces)),
        entity.overlappingFields.map(j => JsonValue(j.noSpaces)),
        entity.suggestedResolution.map(ConflictResolution.toDbString).orNull,
        entity.resolutionReason,
        ConflictStatus.toDbString(entity.status),
        entity.resolvedAt,
        entity.resolvedBy,
        entity.resolutionAction.map(ResolutionAction.toDbString).orNull,
        entity.localSnapshot.map(j => JsonValue(j.noSpaces)),
        entity.remoteSnapshot.map(j => JsonValue(j.noSpaces)),
        now,
        entity.id
      )
    )
    entity.copy(updatedAt = now)

  def delete(id: UUID)(using conn: Connection): Boolean =
    executeUpdate("DELETE FROM sync_conflict WHERE id = ?", Seq(id)) > 0

  // ============================================
  // Query Operations
  // ============================================

  /**
   * Find unresolved conflicts.
   */
  def findUnresolved()(using conn: Connection): List[SyncConflictEntity] =
    queryList(
      "SELECT * FROM sync_conflict WHERE status = 'UNRESOLVED' ORDER BY detected_at ASC"
    )(mapRow)

  /**
   * Find conflicts by entity.
   */
  def findByEntity(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): List[SyncConflictEntity] =
    queryList(
      "SELECT * FROM sync_conflict WHERE entity_type = ? AND entity_id = ? ORDER BY detected_at DESC",
      Seq(SyncEntityType.toDbString(entityType), entityId)
    )(mapRow)

  /**
   * Find active (unresolved) conflict for an entity.
   */
  def findActiveConflict(entityType: SyncEntityType, entityId: UUID)(using conn: Connection): Option[SyncConflictEntity] =
    queryOne(
      "SELECT * FROM sync_conflict WHERE entity_type = ? AND entity_id = ? AND status = 'UNRESOLVED' ORDER BY detected_at DESC LIMIT 1",
      Seq(SyncEntityType.toDbString(entityType), entityId)
    )(mapRow)

  /**
   * Find conflicts by status.
   */
  def findByStatus(status: ConflictStatus)(using conn: Connection): List[SyncConflictEntity] =
    queryList(
      "SELECT * FROM sync_conflict WHERE status = ? ORDER BY detected_at DESC",
      Seq(ConflictStatus.toDbString(status))
    )(mapRow)

  // ============================================
  // Resolution Operations
  // ============================================

  /**
   * Resolve a conflict by keeping local version.
   */
  def resolveKeepLocal(id: UUID, resolvedBy: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_conflict SET
        |  status = 'RESOLVED', resolved_at = ?, resolved_by = ?,
        |  resolution_action = 'KEPT_LOCAL', updated_at = ?
        |WHERE id = ? AND status = 'UNRESOLVED'
      """.stripMargin,
      Seq(now, resolvedBy, now, id)
    ) > 0

  /**
   * Resolve a conflict by accepting remote version.
   */
  def resolveAcceptRemote(id: UUID, resolvedBy: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_conflict SET
        |  status = 'RESOLVED', resolved_at = ?, resolved_by = ?,
        |  resolution_action = 'ACCEPTED_REMOTE', updated_at = ?
        |WHERE id = ? AND status = 'UNRESOLVED'
      """.stripMargin,
      Seq(now, resolvedBy, now, id)
    ) > 0

  /**
   * Resolve a conflict by merging changes.
   */
  def resolveMerge(id: UUID, resolvedBy: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_conflict SET
        |  status = 'RESOLVED', resolved_at = ?, resolved_by = ?,
        |  resolution_action = 'MERGED', updated_at = ?
        |WHERE id = ? AND status = 'UNRESOLVED'
      """.stripMargin,
      Seq(now, resolvedBy, now, id)
    ) > 0

  /**
   * Resolve a conflict with manual edit.
   */
  def resolveManual(id: UUID, resolvedBy: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_conflict SET
        |  status = 'RESOLVED', resolved_at = ?, resolved_by = ?,
        |  resolution_action = 'MANUAL_EDIT', updated_at = ?
        |WHERE id = ? AND status = 'UNRESOLVED'
      """.stripMargin,
      Seq(now, resolvedBy, now, id)
    ) > 0

  /**
   * Dismiss a conflict (ignore it).
   */
  def dismiss(id: UUID, resolvedBy: String)(using conn: Connection): Boolean =
    val now = LocalDateTime.now()
    executeUpdate(
      """UPDATE sync_conflict SET
        |  status = 'DISMISSED', resolved_at = ?, resolved_by = ?, updated_at = ?
        |WHERE id = ? AND status = 'UNRESOLVED'
      """.stripMargin,
      Seq(now, resolvedBy, now, id)
    ) > 0

  // ============================================
  // Statistics
  // ============================================

  /**
   * Count unresolved conflicts.
   */
  def countUnresolved()(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM sync_conflict WHERE status = 'UNRESOLVED'"
    )(_.getLong(1)).getOrElse(0L)

  /**
   * Count conflicts by status.
   */
  def countByStatus()(using conn: Connection): Map[ConflictStatus, Long] =
    queryList(
      "SELECT status, COUNT(*) as cnt FROM sync_conflict GROUP BY status"
    ) { rs =>
      val status = ConflictStatus.fromString(rs.getString("status"))
      val count = rs.getLong("cnt")
      (status, count)
    }.toMap

  /**
   * Cleanup resolved conflicts older than specified days.
   */
  def cleanupResolved(olderThanDays: Int)(using conn: Connection): Int =
    val cutoff = LocalDateTime.now().minusDays(olderThanDays)
    executeUpdate(
      "DELETE FROM sync_conflict WHERE status IN ('RESOLVED', 'DISMISSED') AND resolved_at < ?",
      Seq(cutoff)
    )

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): SyncConflictEntity =
    SyncConflictEntity(
      id = getUUID(rs, "id"),
      entityType = SyncEntityType.fromString(rs.getString("entity_type")),
      entityId = getUUID(rs, "entity_id"),
      atUri = getOptString(rs, "at_uri"),
      detectedAt = getDateTime(rs, "detected_at"),
      localVersion = rs.getInt("local_version"),
      remoteVersion = rs.getInt("remote_version"),
      localChanges = getOptJsonString(rs, "local_changes").flatMap(parse(_).toOption),
      remoteChanges = getOptJsonString(rs, "remote_changes").flatMap(parse(_).toOption),
      overlappingFields = getOptJsonString(rs, "overlapping_fields").flatMap(parse(_).toOption),
      suggestedResolution = getOptString(rs, "suggested_resolution").map(ConflictResolution.fromString),
      resolutionReason = getOptString(rs, "resolution_reason"),
      status = ConflictStatus.fromString(rs.getString("status")),
      resolvedAt = getOptDateTime(rs, "resolved_at"),
      resolvedBy = getOptString(rs, "resolved_by"),
      resolutionAction = getOptString(rs, "resolution_action").map(ResolutionAction.fromString),
      localSnapshot = getOptJsonString(rs, "local_snapshot").flatMap(parse(_).toOption),
      remoteSnapshot = getOptJsonString(rs, "remote_snapshot").flatMap(parse(_).toOption),
      createdAt = getDateTime(rs, "created_at"),
      updatedAt = getDateTime(rs, "updated_at")
    )
