package com.decodingus.db

import munit.FunSuite

class DatabaseSpec extends FunSuite:

  test("initializeInMemory creates a working database") {
    Database.initializeInMemory() match
      case Right(db) =>
        try
          // Verify we can get a connection and execute a query
          db.connection { conn =>
            val stmt = conn.createStatement()
            val rs = stmt.executeQuery("SELECT 1")
            assert(rs.next())
            assertEquals(rs.getInt(1), 1)
            rs.close()
            stmt.close()
          }
        finally
          db.shutdown()
      case Left(error) =>
        fail(s"Failed to initialize in-memory database: $error")
  }

  test("database connection executes queries correctly") {
    Database.initializeInMemory() match
      case Right(db) =>
        try
          val result = db.connection { conn =>
            val stmt = conn.createStatement()
            stmt.execute("CREATE TABLE test_table (id INT PRIMARY KEY, name VARCHAR(100))")
            stmt.execute("INSERT INTO test_table VALUES (1, 'test')")
            val rs = stmt.executeQuery("SELECT name FROM test_table WHERE id = 1")
            rs.next()
            val name = rs.getString("name")
            rs.close()
            stmt.close()
            name
          }
          assertEquals(result, "test")
        finally
          db.shutdown()
      case Left(error) =>
        fail(s"Failed: $error")
  }

  test("shutdown closes connections properly") {
    Database.initializeInMemory() match
      case Right(db) =>
        // Get a connection to verify it works
        db.connection { conn =>
          assert(!conn.isClosed)
        }
        // Shutdown
        db.shutdown()
        // Attempting to get a new connection after shutdown should fail
        // (HikariCP will throw an exception)
        intercept[Exception] {
          db.connection { _ => () }
        }
      case Left(error) =>
        fail(s"Failed: $error")
  }
