package com.decodingus.db

import com.zaxxer.hikari.{HikariConfig, HikariDataSource}
import java.nio.file.{Files, Path}
import java.sql.Connection
import scala.util.Using

/**
 * Database connection management with clean lifecycle.
 *
 * Manages an H2 embedded database with HikariCP connection pooling.
 * The database file is stored at ~/.decodingus/data/workspace.mv.db
 */
final class Database private (dataSource: HikariDataSource):

  /**
   * Execute a function with a connection from the pool.
   * The connection is automatically returned to the pool after use.
   */
  def connection[A](f: Connection => A): A =
    Using.resource(dataSource.getConnection)(f)

  /**
   * Get the underlying data source for advanced operations.
   */
  def getDataSource: HikariDataSource = dataSource

  /**
   * Shutdown the connection pool gracefully.
   */
  def shutdown(): Unit =
    if !dataSource.isClosed then
      dataSource.close()

object Database:
  private val DbDir: Path = Path.of(System.getProperty("user.home"), ".decodingus", "data")

  /**
   * Initialize the database with file-based storage.
   * Creates the data directory if it doesn't exist.
   */
  def initialize(): Either[String, Database] =
    for
      _ <- createDirectories()
      ds <- createDataSource(fileUrl)
    yield Database(ds)

  /**
   * Initialize an in-memory database for testing.
   */
  def initializeInMemory(): Either[String, Database] =
    createDataSource(inMemoryUrl).map(Database(_))

  /**
   * Check if the database file exists.
   */
  def databaseExists: Boolean =
    Files.exists(DbDir.resolve("workspace.mv.db"))

  private def fileUrl: String =
    s"jdbc:h2:file:${DbDir.resolve("workspace")};MODE=PostgreSQL;DATABASE_TO_LOWER=TRUE"

  private def inMemoryUrl: String =
    "jdbc:h2:mem:test;MODE=PostgreSQL;DATABASE_TO_LOWER=TRUE;DB_CLOSE_DELAY=-1"

  private def createDirectories(): Either[String, Unit] =
    try
      Files.createDirectories(DbDir)
      Right(())
    catch
      case e: Exception =>
        Left(s"Failed to create data directory: ${e.getMessage}")

  private def createDataSource(jdbcUrl: String): Either[String, HikariDataSource] =
    try
      val config = HikariConfig()
      config.setJdbcUrl(jdbcUrl)
      config.setUsername("sa")
      config.setPassword("")
      config.setDriverClassName("org.h2.Driver")

      // Pool configuration tuned for desktop application
      config.setMaximumPoolSize(5)
      config.setMinimumIdle(1)
      config.setIdleTimeout(300000)      // 5 minutes
      config.setMaxLifetime(600000)      // 10 minutes
      config.setConnectionTimeout(10000) // 10 seconds
      config.setPoolName("DUNavigator-H2")

      Right(HikariDataSource(config))
    catch
      case e: Exception =>
        Left(s"Failed to create connection pool: ${e.getMessage}")
