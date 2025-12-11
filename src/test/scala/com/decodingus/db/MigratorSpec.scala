package com.decodingus.db

import munit.FunSuite

class MigratorSpec extends FunSuite:

  test("migrate applies all migrations to fresh database") {
    Database.initializeInMemory() match
      case Right(db) =>
        try
          Migrator.migrate(db) match
            case Right(count) =>
              assert(count >= 3, s"Expected at least 3 migrations, got $count")
            case Left(error) =>
              fail(s"Migration failed: $error")

          // Verify tables were created (case-insensitive for H2 PostgreSQL mode)
          db.connection { conn =>
            val meta = conn.getMetaData
            val tables = meta.getTables(null, null, "%", Array("TABLE"))
            val tableNames = collection.mutable.Set.empty[String]
            while tables.next() do
              tableNames += tables.getString("TABLE_NAME").toUpperCase
            tables.close()

            // Core tables from V001
            assert(tableNames.contains("BIOSAMPLE"), "BIOSAMPLE table should exist")
            assert(tableNames.contains("PROJECT"), "PROJECT table should exist")
            assert(tableNames.contains("SEQUENCE_RUN"), "SEQUENCE_RUN table should exist")
            assert(tableNames.contains("ALIGNMENT"), "ALIGNMENT table should exist")

            // Cache tables from V002
            assert(tableNames.contains("SOURCE_FILE"), "SOURCE_FILE table should exist")
            assert(tableNames.contains("ANALYSIS_ARTIFACT"), "ANALYSIS_ARTIFACT table should exist")

            // Sync tables from V003
            assert(tableNames.contains("SYNC_QUEUE"), "SYNC_QUEUE table should exist")
            assert(tableNames.contains("SYNC_HISTORY"), "SYNC_HISTORY table should exist")
            assert(tableNames.contains("SYNC_CONFLICT"), "SYNC_CONFLICT table should exist")
          }
        finally
          db.shutdown()
      case Left(error) =>
        fail(s"Database init failed: $error")
  }

  test("migrate is idempotent") {
    Database.initializeInMemory() match
      case Right(db) =>
        try
          // First migration
          val first = Migrator.migrate(db)
          assert(first.isRight, s"First migration failed: ${first.left.getOrElse("")}")

          // Second migration should apply 0 new migrations
          Migrator.migrate(db) match
            case Right(count) =>
              assertEquals(count, 0, "Second migration should apply 0 new migrations")
            case Left(error) =>
              fail(s"Second migration failed: $error")
        finally
          db.shutdown()
      case Left(error) =>
        fail(s"Database init failed: $error")
  }

  test("getCurrentVersion returns correct version after migrations") {
    Database.initializeInMemory() match
      case Right(db) =>
        try
          // Before migration
          Migrator.getCurrentVersion(db) match
            case Right(version) =>
              assertEquals(version, 0, "Initial version should be 0")
            case Left(error) =>
              fail(s"Failed to get version: $error")

          // Apply migrations
          Migrator.migrate(db)

          // After migration
          Migrator.getCurrentVersion(db) match
            case Right(version) =>
              assert(version >= 3, s"Version should be at least 3, got $version")
            case Left(error) =>
              fail(s"Failed to get version: $error")
        finally
          db.shutdown()
      case Left(error) =>
        fail(s"Database init failed: $error")
  }

  test("schema_version table tracks applied migrations") {
    Database.initializeInMemory() match
      case Right(db) =>
        try
          Migrator.migrate(db)

          db.connection { conn =>
            val stmt = conn.createStatement()
            val rs = stmt.executeQuery("SELECT version, description FROM schema_version ORDER BY version")

            val migrations = collection.mutable.ListBuffer.empty[(Int, String)]
            while rs.next() do
              migrations += ((rs.getInt("version"), rs.getString("description")))
            rs.close()
            stmt.close()

            assert(migrations.nonEmpty, "Should have recorded migrations")
            assertEquals(migrations.head._1, 1, "First migration should be V1")
            assert(migrations.exists(_._2.contains("initial")), "Should have initial schema migration")
          }
        finally
          db.shutdown()
      case Left(error) =>
        fail(s"Database init failed: $error")
  }
