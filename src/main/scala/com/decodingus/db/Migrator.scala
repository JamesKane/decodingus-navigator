package com.decodingus.db

import com.decodingus.util.Logger

import java.io.{BufferedReader, InputStreamReader}
import java.nio.charset.StandardCharsets
import java.security.MessageDigest
import java.sql.{Connection, Timestamp}
import java.time.LocalDateTime
import scala.collection.mutable.ListBuffer
import scala.util.Using

/**
 * Schema migration manager for the H2 database.
 *
 * Migrations are versioned SQL scripts stored in resources/db/migration/
 * Named: V{version}__description.sql (e.g., V001__initial_schema.sql)
 *
 * Features:
 * - Forward-only migrations (no rollback)
 * - Checksum validation to detect modified scripts
 * - Idempotent execution (tracks applied versions)
 */
object Migrator:

  private val log = Logger("Migrator")
  private val MigrationPath = "db/migration"
  private val SchemaVersionTable = "schema_version"

  /**
   * A migration script to apply.
   */
  case class Migration(
                        version: Int,
                        description: String,
                        sql: String,
                        checksum: String
                      )

  /**
   * A migration that has been applied.
   */
  case class AppliedMigration(
                               version: Int,
                               description: String,
                               appliedAt: LocalDateTime,
                               checksum: String
                             )

  /**
   * Run all pending migrations.
   *
   * @param database The database to migrate
   * @return Either an error or the count of migrations applied
   */
  def migrate(database: Database): Either[String, Int] =
    database.connection { conn =>
      try
        conn.setAutoCommit(false)

        // Ensure tracking table exists
        ensureSchemaVersionTable(conn)

        // Load applied and available migrations
        val applied = getAppliedMigrations(conn)
        val appliedVersions = applied.map(_.version).toSet
        val available = loadMigrations()

        // Validate checksums of applied migrations
        validateChecksums(applied, available)

        // Apply pending migrations in order
        val pending = available
          .filterNot(m => appliedVersions.contains(m.version))
          .sortBy(_.version)

        for migration <- pending do
          applyMigration(conn, migration)

        conn.commit()
        Right(pending.size)
      catch
        case e: Exception =>
          conn.rollback()
          Left(s"Migration failed: ${e.getMessage}")
      finally
        conn.setAutoCommit(true)
    }

  /**
   * Get the current schema version.
   */
  def getCurrentVersion(database: Database): Either[String, Int] =
    try
      Right(database.connection { conn =>
        if !tableExists(conn, SchemaVersionTable) then 0
        else
          Using.resource(conn.createStatement()) { stmt =>
            Using.resource(stmt.executeQuery(s"SELECT MAX(version) FROM $SchemaVersionTable")) { rs =>
              if rs.next() then rs.getInt(1) else 0
            }
          }
      })
    catch
      case e: Exception => Left(e.getMessage)

  private def ensureSchemaVersionTable(conn: Connection): Unit =
    if !tableExists(conn, SchemaVersionTable) then
      Using.resource(conn.createStatement()) { stmt =>
        stmt.execute(
          """
          CREATE TABLE schema_version (
            version INT PRIMARY KEY,
            description VARCHAR(255) NOT NULL,
            applied_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            checksum VARCHAR(64) NOT NULL
          )
        """)
      }

  private def tableExists(conn: Connection, tableName: String): Boolean =
    // H2 in PostgreSQL mode stores table names lowercase, so check both cases
    val rsUpper = conn.getMetaData.getTables(null, null, tableName.toUpperCase, Array("TABLE"))
    val existsUpper = Using.resource(rsUpper)(_.next())
    if existsUpper then true
    else
      val rsLower = conn.getMetaData.getTables(null, null, tableName.toLowerCase, Array("TABLE"))
      Using.resource(rsLower)(_.next())

  private def getAppliedMigrations(conn: Connection): List[AppliedMigration] =
    if !tableExists(conn, SchemaVersionTable) then List.empty
    else
      Using.resource(conn.createStatement()) { stmt =>
        Using.resource(stmt.executeQuery(
          s"SELECT version, description, applied_at, checksum FROM $SchemaVersionTable ORDER BY version"
        )) { rs =>
          val migrations = ListBuffer.empty[AppliedMigration]
          while rs.next() do
            migrations += AppliedMigration(
              version = rs.getInt("version"),
              description = rs.getString("description"),
              appliedAt = rs.getTimestamp("applied_at").toLocalDateTime,
              checksum = rs.getString("checksum")
            )
          migrations.toList
        }
      }

  private def loadMigrations(): List[Migration] =
    // Known migrations - in production, could scan classpath
    val migrationFiles = List(
      "V001__initial_schema.sql",
      "V002__file_cache_tables.sql",
      "V003__sync_queue_tables.sql",
      "V004__add_artifact_unique_constraint.sql",
      "V005__phase2_entity_tables.sql",
      "V006__y_chromosome_unified_profile.sql",
      "V007__multi_reference_variant_model.sql",
      "V008__add_sequence_run_metrics_columns.sql"
    )
    val pattern = """V(\d+)__(.+)\.sql""".r

    migrationFiles.flatMap { fileName =>
      fileName match
        case pattern(versionStr, desc) =>
          loadResourceAsString(s"$MigrationPath/$fileName").map { sql =>
            Migration(
              version = versionStr.toInt,
              description = desc.replace("_", " "),
              sql = sql,
              checksum = calculateChecksum(sql)
            )
          }
        case _ => None
    }

  private def loadResourceAsString(path: String): Option[String] =
    val stream = getClass.getClassLoader.getResourceAsStream(path)
    if stream == null then None
    else
      Some(Using.resource(new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) { reader =>
        val sb = new StringBuilder()
        var line = reader.readLine()
        while line != null do
          sb.append(line).append("\n")
          line = reader.readLine()
        sb.toString()
      })

  private def calculateChecksum(sql: String): String =
    val digest = MessageDigest.getInstance("SHA-256")
    val bytes = digest.digest(sql.getBytes(StandardCharsets.UTF_8))
    bytes.map("%02x".format(_)).mkString.take(16)

  private def validateChecksums(applied: List[AppliedMigration], available: List[Migration]): Unit =
    val availableMap = available.map(m => m.version -> m.checksum).toMap
    for am <- applied do
      availableMap.get(am.version).foreach { expectedChecksum =>
        if expectedChecksum != am.checksum then
          throw new RuntimeException(
            s"Checksum mismatch for V${am.version}: script was modified after being applied"
          )
      }

  /**
   * Parse SQL into individual statements.
   * Handles semicolons inside comments correctly by first stripping comments.
   */
  private def parseStatements(sql: String): List[String] =
    // First, strip all comment lines from the entire SQL
    val noComments = sql.linesIterator
      .filterNot(_.trim.startsWith("--"))
      .mkString("\n")

    // Now split by semicolons (safe because comments are gone)
    noComments
      .split(";")
      .map(_.trim)
      .filter(_.nonEmpty)
      .toList

  private def applyMigration(conn: Connection, migration: Migration): Unit =
    log.info(s"Applying V${migration.version} - ${migration.description}")

    // Execute migration SQL statements
    Using.resource(conn.createStatement()) { stmt =>
      val statements = parseStatements(migration.sql)

      for sql <- statements do
        stmt.execute(sql)
    }

    // Record the migration
    Using.resource(conn.prepareStatement(
      s"INSERT INTO $SchemaVersionTable (version, description, checksum) VALUES (?, ?, ?)"
    )) { ps =>
      ps.setInt(1, migration.version)
      ps.setString(2, migration.description)
      ps.setString(3, migration.checksum)
      ps.executeUpdate()
    }

    log.info(s"Applied V${migration.version}")
