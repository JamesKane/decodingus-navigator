package com.decodingus.db

import munit.FunSuite

/**
 * Base trait for database tests providing in-memory H2 database fixture.
 *
 * Each test gets a fresh database with migrations applied.
 */
trait DatabaseTestSupport extends FunSuite:

  /**
   * Fixture that provides a fresh Database instance for each test.
   */
  val testDatabase: FunFixture[Database] = FunFixture[Database](
    setup = { _ =>
      Database.initializeInMemory() match
        case Right(db) =>
          Migrator.migrate(db) match
            case Right(_) => db
            case Left(error) =>
              db.shutdown()
              throw new RuntimeException(s"Migration failed: $error")
        case Left(error) =>
          throw new RuntimeException(s"Database init failed: $error")
    },
    teardown = { db =>
      db.shutdown()
    }
  )

  /**
   * Fixture that provides a Database and Transactor pair.
   */
  val testTransactor: FunFixture[(Database, Transactor)] = FunFixture[(Database, Transactor)](
    setup = { _ =>
      Database.initializeInMemory() match
        case Right(db) =>
          Migrator.migrate(db) match
            case Right(_) => (db, Transactor(db))
            case Left(error) =>
              db.shutdown()
              throw new RuntimeException(s"Migration failed: $error")
        case Left(error) =>
          throw new RuntimeException(s"Database init failed: $error")
    },
    teardown = { case (db, _) =>
      db.shutdown()
    }
  )
