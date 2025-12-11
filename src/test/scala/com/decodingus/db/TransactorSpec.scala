package com.decodingus.db

import munit.FunSuite
import java.sql.Connection

class TransactorSpec extends FunSuite with DatabaseTestSupport:

  testDatabase.test("readOnly executes queries successfully") { db =>
    val transactor = Transactor(db)

    val result = transactor.readOnly {
      val conn = summon[Connection]
      val stmt = conn.createStatement()
      val rs = stmt.executeQuery("SELECT COUNT(*) FROM biosample")
      rs.next()
      val count = rs.getLong(1)
      rs.close()
      stmt.close()
      count
    }

    assertEquals(result, Right(0L))
  }

  testDatabase.test("readWrite commits successful transactions") { db =>
    val transactor = Transactor(db)

    // Insert in a transaction
    val insertResult = transactor.readWrite {
      val conn = summon[Connection]
      val stmt = conn.prepareStatement(
        """INSERT INTO biosample (id, sample_accession, donor_identifier, sync_status, version, created_at, updated_at)
          |VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)""".stripMargin
      )
      stmt.setObject(1, java.util.UUID.randomUUID())
      stmt.setString(2, "TEST001")
      stmt.setString(3, "DONOR001")
      stmt.setString(4, "Local")
      stmt.setInt(5, 1)
      val rows = stmt.executeUpdate()
      stmt.close()
      rows
    }

    assertEquals(insertResult, Right(1))

    // Verify it was committed
    val count = transactor.readOnly {
      val conn = summon[Connection]
      val stmt = conn.createStatement()
      val rs = stmt.executeQuery("SELECT COUNT(*) FROM biosample WHERE sample_accession = 'TEST001'")
      rs.next()
      val c = rs.getLong(1)
      rs.close()
      stmt.close()
      c
    }

    assertEquals(count, Right(1L))
  }

  testDatabase.test("readWrite rolls back failed transactions") { db =>
    val transactor = Transactor(db)

    // First, insert a valid record
    transactor.readWrite {
      val conn = summon[Connection]
      val stmt = conn.prepareStatement(
        """INSERT INTO biosample (id, sample_accession, donor_identifier, sync_status, version, created_at, updated_at)
          |VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)""".stripMargin
      )
      stmt.setObject(1, java.util.UUID.randomUUID())
      stmt.setString(2, "EXISTING001")
      stmt.setString(3, "DONOR001")
      stmt.setString(4, "Local")
      stmt.setInt(5, 1)
      stmt.executeUpdate()
      stmt.close()
    }

    // Try to insert a duplicate (should fail and rollback)
    val failedResult = transactor.readWrite {
      val conn = summon[Connection]
      // First insert something new
      val stmt1 = conn.prepareStatement(
        """INSERT INTO biosample (id, sample_accession, donor_identifier, sync_status, version, created_at, updated_at)
          |VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)""".stripMargin
      )
      stmt1.setObject(1, java.util.UUID.randomUUID())
      stmt1.setString(2, "NEWRECORD")
      stmt1.setString(3, "DONOR002")
      stmt1.setString(4, "Local")
      stmt1.setInt(5, 1)
      stmt1.executeUpdate()
      stmt1.close()

      // Now try to insert a duplicate accession (violates unique constraint)
      val stmt2 = conn.prepareStatement(
        """INSERT INTO biosample (id, sample_accession, donor_identifier, sync_status, version, created_at, updated_at)
          |VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)""".stripMargin
      )
      stmt2.setObject(1, java.util.UUID.randomUUID())
      stmt2.setString(2, "EXISTING001") // Duplicate!
      stmt2.setString(3, "DONOR003")
      stmt2.setString(4, "Local")
      stmt2.setInt(5, 1)
      stmt2.executeUpdate() // This will throw
      stmt2.close()
    }

    assert(failedResult.isLeft, "Transaction should have failed")

    // Verify NEWRECORD was NOT committed (rollback worked)
    val count = transactor.readOnly {
      val conn = summon[Connection]
      val stmt = conn.createStatement()
      val rs = stmt.executeQuery("SELECT COUNT(*) FROM biosample WHERE sample_accession = 'NEWRECORD'")
      rs.next()
      val c = rs.getLong(1)
      rs.close()
      stmt.close()
      c
    }

    assertEquals(count, Right(0L), "NEWRECORD should have been rolled back")
  }

  testDatabase.test("readOnly prevents write operations") { db =>
    val transactor = Transactor(db)

    // readOnly should work for SELECT queries
    val readResult = transactor.readOnly {
      val conn = summon[Connection]
      val stmt = conn.createStatement()
      val rs = stmt.executeQuery("SELECT 1")
      rs.next()
      val v = rs.getInt(1)
      rs.close()
      stmt.close()
      v
    }

    assertEquals(readResult, Right(1))
  }

  testDatabase.test("nested operations work correctly") { db =>
    val transactor = Transactor(db)
    val id = java.util.UUID.randomUUID()

    // Create in one transaction
    transactor.readWrite {
      val conn = summon[Connection]
      val stmt = conn.prepareStatement(
        """INSERT INTO biosample (id, sample_accession, donor_identifier, sync_status, version, created_at, updated_at)
          |VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)""".stripMargin
      )
      stmt.setObject(1, id)
      stmt.setString(2, "NESTED001")
      stmt.setString(3, "DONOR001")
      stmt.setString(4, "Local")
      stmt.setInt(5, 1)
      stmt.executeUpdate()
      stmt.close()
    }

    // Read in another transaction
    val accession = transactor.readOnly {
      val conn = summon[Connection]
      val stmt = conn.prepareStatement("SELECT sample_accession FROM biosample WHERE id = ?")
      stmt.setObject(1, id)
      val rs = stmt.executeQuery()
      rs.next()
      val result = rs.getString(1)
      rs.close()
      stmt.close()
      result
    }

    assertEquals(accession, Right("NESTED001"))
  }
