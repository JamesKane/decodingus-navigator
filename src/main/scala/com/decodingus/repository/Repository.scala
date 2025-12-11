package com.decodingus.repository

import java.sql.{Connection, PreparedStatement, ResultSet, Timestamp}
import java.time.LocalDateTime
import java.util.UUID
import scala.util.Using

/**
 * Sync status for PDS integration.
 */
enum SyncStatus:
  case Local     // Never synced to PDS
  case Synced    // In sync with PDS
  case Modified  // Local changes pending sync
  case Conflict  // Diverged from PDS, needs resolution

object SyncStatus:
  def fromString(s: String): SyncStatus = s match
    case "Local"    => Local
    case "Synced"   => Synced
    case "Modified" => Modified
    case "Conflict" => Conflict
    case other      => throw new IllegalArgumentException(s"Unknown sync status: $other")

/**
 * Base trait for entities with a typed identifier.
 */
trait Entity[ID]:
  def id: ID

/**
 * Metadata common to all persisted entities.
 */
case class EntityMeta(
  syncStatus: SyncStatus,
  atUri: Option[String],
  atCid: Option[String],
  version: Int,
  createdAt: LocalDateTime,
  updatedAt: LocalDateTime
)

object EntityMeta:
  def create(): EntityMeta = EntityMeta(
    syncStatus = SyncStatus.Local,
    atUri = None,
    atCid = None,
    version = 1,
    createdAt = LocalDateTime.now(),
    updatedAt = LocalDateTime.now()
  )

  def forUpdate(existing: EntityMeta): EntityMeta = existing.copy(
    syncStatus = if existing.syncStatus == SyncStatus.Synced then SyncStatus.Modified else existing.syncStatus,
    version = existing.version + 1,
    updatedAt = LocalDateTime.now()
  )

/**
 * Base repository operations.
 *
 * @tparam E Entity type
 * @tparam ID Identifier type
 */
trait Repository[E <: Entity[ID], ID]:
  /**
   * Find an entity by its identifier.
   */
  def findById(id: ID)(using Connection): Option[E]

  /**
   * Find all entities.
   */
  def findAll()(using Connection): List[E]

  /**
   * Insert a new entity.
   */
  def insert(entity: E)(using Connection): E

  /**
   * Update an existing entity.
   */
  def update(entity: E)(using Connection): E

  /**
   * Delete an entity by its identifier.
   */
  def delete(id: ID)(using Connection): Boolean

  /**
   * Count all entities.
   */
  def count()(using Connection): Long

  /**
   * Check if an entity exists by its identifier.
   */
  def exists(id: ID)(using Connection): Boolean

/**
 * Repository with sync status tracking.
 */
trait SyncableRepository[E <: Entity[ID], ID] extends Repository[E, ID]:
  /**
   * Find entities by sync status.
   */
  def findByStatus(status: SyncStatus)(using Connection): List[E]

  /**
   * Update sync status for an entity.
   */
  def updateStatus(id: ID, status: SyncStatus)(using Connection): Boolean

  /**
   * Mark an entity as synced with PDS.
   */
  def markSynced(id: ID, atUri: String, atCid: String)(using Connection): Boolean

/**
 * SQL helper utilities for repositories.
 */
object SqlHelpers:

  /**
   * Set a parameter on a PreparedStatement, handling Option types.
   */
  def setParam(ps: PreparedStatement, index: Int, value: Any): Unit =
    value match
      case null          => ps.setNull(index, java.sql.Types.NULL)
      case None          => ps.setNull(index, java.sql.Types.NULL)
      case Some(v)       => setParam(ps, index, v)
      case JsonValue(j)  => ps.setBytes(index, jsonToBytes(j))  // H2 JSON columns use byte[]
      case s: String     => ps.setString(index, s)
      case i: Int        => ps.setInt(index, i)
      case l: Long       => ps.setLong(index, l)
      case d: Double     => ps.setDouble(index, d)
      case b: Boolean    => ps.setBoolean(index, b)
      case u: UUID       => ps.setObject(index, u)
      case t: Timestamp  => ps.setTimestamp(index, t)
      case ldt: LocalDateTime => ps.setTimestamp(index, Timestamp.valueOf(ldt))
      case ss: SyncStatus => ps.setString(index, ss.toString)
      case other         => ps.setObject(index, other)

  /**
   * Get an optional string from a ResultSet.
   */
  def getOptString(rs: ResultSet, column: String): Option[String] =
    Option(rs.getString(column))

  /**
   * Get an optional JSON string from a ResultSet.
   * H2's JSON columns return byte[] when read with getBytes().
   * This helper properly extracts the JSON content as a UTF-8 string.
   */
  def getOptJsonString(rs: ResultSet, column: String): Option[String] =
    val bytes = rs.getBytes(column)
    if bytes == null then None
    else Some(new String(bytes, java.nio.charset.StandardCharsets.UTF_8))

  /**
   * Wrapper type for JSON string values to distinguish from regular strings.
   * Use this when passing JSON data to setParam to ensure proper H2 JSON column handling.
   */
  case class JsonValue(json: String)

  /**
   * Convert a JSON string to bytes for H2 JSON column storage.
   * H2 JSON type maps to byte[] in Java, so we convert the string to UTF-8 bytes.
   */
  def jsonToBytes(json: String): Array[Byte] =
    if json == null then null
    else json.getBytes(java.nio.charset.StandardCharsets.UTF_8)

  /**
   * Get an optional int from a ResultSet.
   */
  def getOptInt(rs: ResultSet, column: String): Option[Int] =
    val value = rs.getInt(column)
    if rs.wasNull() then None else Some(value)

  /**
   * Get an optional long from a ResultSet.
   */
  def getOptLong(rs: ResultSet, column: String): Option[Long] =
    val value = rs.getLong(column)
    if rs.wasNull() then None else Some(value)

  /**
   * Get an optional double from a ResultSet.
   */
  def getOptDouble(rs: ResultSet, column: String): Option[Double] =
    val value = rs.getDouble(column)
    if rs.wasNull() then None else Some(value)

  /**
   * Get an optional timestamp as LocalDateTime from a ResultSet.
   */
  def getOptDateTime(rs: ResultSet, column: String): Option[LocalDateTime] =
    Option(rs.getTimestamp(column)).map(_.toLocalDateTime)

  /**
   * Get a required timestamp as LocalDateTime from a ResultSet.
   */
  def getDateTime(rs: ResultSet, column: String): LocalDateTime =
    rs.getTimestamp(column).toLocalDateTime

  /**
   * Get a UUID from a ResultSet.
   */
  def getUUID(rs: ResultSet, column: String): UUID =
    rs.getObject(column, classOf[UUID])

  /**
   * Get an optional UUID from a ResultSet.
   */
  def getOptUUID(rs: ResultSet, column: String): Option[UUID] =
    val value = rs.getObject(column, classOf[UUID])
    if rs.wasNull() then None else Option(value)

  /**
   * Execute a query and map results to a list.
   */
  def queryList[A](sql: String, params: Seq[Any] = Seq.empty)(mapper: ResultSet => A)(using conn: Connection): List[A] =
    Using.resource(conn.prepareStatement(sql)) { ps =>
      params.zipWithIndex.foreach { case (param, idx) =>
        setParam(ps, idx + 1, param)
      }
      Using.resource(ps.executeQuery()) { rs =>
        val results = scala.collection.mutable.ListBuffer.empty[A]
        while rs.next() do
          results += mapper(rs)
        results.toList
      }
    }

  /**
   * Execute a query and return an optional single result.
   */
  def queryOne[A](sql: String, params: Seq[Any] = Seq.empty)(mapper: ResultSet => A)(using conn: Connection): Option[A] =
    Using.resource(conn.prepareStatement(sql)) { ps =>
      params.zipWithIndex.foreach { case (param, idx) =>
        setParam(ps, idx + 1, param)
      }
      Using.resource(ps.executeQuery()) { rs =>
        if rs.next() then Some(mapper(rs)) else None
      }
    }

  /**
   * Execute an update and return the number of affected rows.
   */
  def executeUpdate(sql: String, params: Seq[Any])(using conn: Connection): Int =
    Using.resource(conn.prepareStatement(sql)) { ps =>
      params.zipWithIndex.foreach { case (param, idx) =>
        setParam(ps, idx + 1, param)
      }
      ps.executeUpdate()
    }
