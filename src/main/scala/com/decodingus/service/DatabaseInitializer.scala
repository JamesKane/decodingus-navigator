package com.decodingus.service

import com.decodingus.db.{Database, Migrator, Transactor}
import com.decodingus.repository.*

/**
 * Database initialization and service wiring.
 *
 * Provides a clean entry point for application startup:
 * 1. Initialize the database connection pool
 * 2. Run schema migrations
 * 3. Create and wire up repositories and services
 *
 * Usage:
 * ```scala
 * DatabaseInitializer.initialize() match
 *   case Right(context) =>
 *     val workspaceService = context.workspaceService
 *     // Use the service...
 *     // On shutdown:
 *     context.shutdown()
 *   case Left(error) =>
 *     // Handle initialization error
 * ```
 */
object DatabaseInitializer:

  /**
   * Initialize the database and create the application context.
   *
   * @param inMemory If true, use an in-memory database (for testing)
   * @return Either an error message or the initialized context
   */
  def initialize(inMemory: Boolean = false): Either[String, DatabaseContext] =
    for
      database <- if inMemory then Database.initializeInMemory() else Database.initialize()
      _ <- runMigrations(database)
      context <- createContext(database)
    yield context

  /**
   * Initialize with an existing database instance (for testing).
   */
  def initializeWithDatabase(database: Database): Either[String, DatabaseContext] =
    for
      _ <- runMigrations(database)
      context <- createContext(database)
    yield context

  private def runMigrations(database: Database): Either[String, Int] =
    Migrator.migrate(database)

  private def createContext(database: Database): Either[String, DatabaseContext] =
    try
      val transactor = Transactor(database)

      // Create repositories
      val biosampleRepo = BiosampleRepository()
      val projectRepo = ProjectRepository()
      val sequenceRunRepo = SequenceRunRepository()
      val alignmentRepo = AlignmentRepository()

      // Create the workspace service
      val workspaceService = H2WorkspaceService(
        transactor = transactor,
        biosampleRepo = biosampleRepo,
        projectRepo = projectRepo,
        sequenceRunRepo = sequenceRunRepo,
        alignmentRepo = alignmentRepo
      )

      Right(DatabaseContext(
        database = database,
        transactor = transactor,
        workspaceService = workspaceService,
        biosampleRepository = biosampleRepo,
        projectRepository = projectRepo,
        sequenceRunRepository = sequenceRunRepo,
        alignmentRepository = alignmentRepo
      ))
    catch
      case e: Exception =>
        database.shutdown()
        Left(s"Failed to create application context: ${e.getMessage}")

/**
 * Application context holding all database-related components.
 *
 * Provides access to:
 * - The high-level WorkspaceService (recommended for most use cases)
 * - Individual repositories (for advanced/direct database access)
 * - The transactor (for custom transactions)
 */
case class DatabaseContext(
  database: Database,
  transactor: Transactor,
  workspaceService: WorkspaceService,
  biosampleRepository: BiosampleRepository,
  projectRepository: ProjectRepository,
  sequenceRunRepository: SequenceRunRepository,
  alignmentRepository: AlignmentRepository
):
  /**
   * Shutdown the database connection pool.
   * Call this on application exit.
   */
  def shutdown(): Unit =
    database.shutdown()

  /**
   * Check if the database is initialized and has the schema.
   */
  def isInitialized: Boolean =
    Migrator.getCurrentVersion(database) match
      case Right(version) => version > 0
      case Left(_) => false

  /**
   * Get the current schema version.
   */
  def schemaVersion: Either[String, Int] =
    Migrator.getCurrentVersion(database)
