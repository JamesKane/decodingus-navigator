package com.decodingus.db

import com.decodingus.util.Logger
import java.sql.{Connection, SQLException}

/**
 * Transaction management with clean Scala 3 context parameters.
 *
 * Provides read-only and read-write transaction boundaries.
 * Uses `Connection ?=>` context functions for implicit connection passing.
 */
final class Transactor(database: Database):
  private val log = Logger[Transactor]

  /**
   * Execute a read-only operation.
   * The connection is marked read-only for potential optimizations.
   *
   * @param f Function that receives a Connection via context parameter
   * @return Either an error message or the result
   */
  def readOnly[A](f: Connection ?=> A): Either[String, A] =
    try
      Right(database.connection { conn =>
        conn.setReadOnly(true)
        try
          f(using conn)
        finally
          conn.setReadOnly(false)
      })
    catch
      case e: SQLException =>
        log.error(s"Read operation failed - SQL State: ${e.getSQLState}, Error Code: ${e.getErrorCode}", e)
        Left(s"Read operation failed: ${e.getMessage} [SQLState: ${e.getSQLState}]")
      case e: Exception =>
        log.error(s"Read operation failed: ${e.getMessage}", e)
        Left(s"Read operation failed: ${e.getMessage}")

  /**
   * Execute a read-write operation within a transaction.
   * Automatically commits on success, rolls back on failure.
   *
   * @param f Function that receives a Connection via context parameter
   * @return Either an error message or the result
   */
  def readWrite[A](f: Connection ?=> A): Either[String, A] =
    try
      Right(database.connection { conn =>
        conn.setAutoCommit(false)
        try
          val result = f(using conn)
          conn.commit()
          result
        catch
          case e: Exception =>
            conn.rollback()
            throw e
        finally
          conn.setAutoCommit(true)
      })
    catch
      case e: SQLException =>
        log.error(s"Transaction failed - SQL State: ${e.getSQLState}, Error Code: ${e.getErrorCode}", e)
        Left(s"Transaction failed: ${e.getMessage} [SQLState: ${e.getSQLState}]")
      case e: Exception =>
        log.error(s"Transaction failed: ${e.getMessage}", e)
        Left(s"Transaction failed: ${e.getMessage}")

  /**
   * Execute multiple operations in a single transaction.
   * All operations succeed or all are rolled back.
   *
   * @param operations Sequence of operations to execute
   * @return Either an error message or the list of results
   */
  def batch[A](operations: Seq[Connection ?=> Either[String, A]]): Either[String, List[A]] =
    readWrite {
      val results = operations.foldLeft(Right(List.empty[A]): Either[String, List[A]]) {
        case (Right(acc), op) =>
          op match
            case Right(a) => Right(acc :+ a)
            case Left(err) => Left(err)
        case (left, _) => left
      }
      results match
        case Right(list) => list
        case Left(err) => throw new RuntimeException(err)
    }
