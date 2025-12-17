package com.decodingus.util

import org.apache.logging.log4j.{LogManager, Logger as Log4jLogger}

/**
 * Simple logging utility wrapping Log4j2 for consistent logging across the application.
 *
 * Usage:
 * {{{
 * // Option 1: Create a logger for a class
 * class MyClass {
 *   private val log = Logger[MyClass]
 *
 *   def doSomething(): Unit = {
 *     log.debug("Starting operation")
 *     log.info("Processing complete")
 *     log.warn("Something unexpected happened")
 *     log.error("Operation failed", exception)
 *   }
 * }
 *
 * // Option 2: Use the companion object for one-off logging
 * Logger.info("Application started")
 *
 * // Option 3: Create named logger
 * val log = Logger("MyComponent")
 * }}}
 */
class Logger private(private val underlying: Log4jLogger) {

  def trace(message: => String): Unit =
    if (underlying.isTraceEnabled) underlying.trace(message)

  def trace(message: => String, throwable: Throwable): Unit =
    if (underlying.isTraceEnabled) underlying.trace(message, throwable)

  def debug(message: => String): Unit =
    if (underlying.isDebugEnabled) underlying.debug(message)

  def debug(message: => String, throwable: Throwable): Unit =
    if (underlying.isDebugEnabled) underlying.debug(message, throwable)

  def info(message: => String): Unit =
    if (underlying.isInfoEnabled) underlying.info(message)

  def info(message: => String, throwable: Throwable): Unit =
    if (underlying.isInfoEnabled) underlying.info(message, throwable)

  def warn(message: => String): Unit =
    if (underlying.isWarnEnabled) underlying.warn(message)

  def warn(message: => String, throwable: Throwable): Unit =
    if (underlying.isWarnEnabled) underlying.warn(message, throwable)

  def error(message: => String): Unit =
    underlying.error(message)

  def error(message: => String, throwable: Throwable): Unit =
    underlying.error(message, throwable)

  def isDebugEnabled: Boolean = underlying.isDebugEnabled

  def isTraceEnabled: Boolean = underlying.isTraceEnabled

  def isInfoEnabled: Boolean = underlying.isInfoEnabled
}

object Logger {
  /** Create a logger for a specific class */
  def apply[T](implicit ct: scala.reflect.ClassTag[T]): Logger =
    new Logger(LogManager.getLogger(ct.runtimeClass))

  /** Create a logger with a specific name */
  def apply(name: String): Logger =
    new Logger(LogManager.getLogger(name))

  /** Create a logger for an object's class */
  def forObject(obj: AnyRef): Logger =
    new Logger(LogManager.getLogger(obj.getClass))

  // Root logger for static/one-off logging
  private lazy val root = new Logger(LogManager.getRootLogger)

  /** Log at trace level using root logger */
  def trace(message: => String): Unit = root.trace(message)

  /** Log at debug level using root logger */
  def debug(message: => String): Unit = root.debug(message)

  /** Log at info level using root logger */
  def info(message: => String): Unit = root.info(message)

  /** Log at warn level using root logger */
  def warn(message: => String): Unit = root.warn(message)

  /** Log at error level using root logger */
  def error(message: => String): Unit = root.error(message)

  /** Log at error level with exception using root logger */
  def error(message: => String, throwable: Throwable): Unit = root.error(message, throwable)
}
