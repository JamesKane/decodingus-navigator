package com.decodingus.analysis

/**
 * Unified progress reporting utilities for analysis processors.
 *
 * Provides a consistent interface for progress callbacks across different
 * analysis stages, with adapters for different callback signatures.
 */
object ProgressReporter {

  /**
   * Standard progress callback type using fractional progress.
   *
   * Parameters:
   *   - message: Human-readable status message
   *   - current: Current progress value (0.0 to 1.0 when total is 1.0)
   *   - total: Total value (typically 1.0 for fractional, or step count for discrete)
   */
  type ProgressCallback = (String, Double, Double) => Unit

  /**
   * No-op progress callback for when progress reporting is not needed.
   */
  val NoOp: ProgressCallback = (_, _, _) => ()

  /**
   * Convert a step-based callback to a fraction-based callback.
   *
   * Use this when you have discrete steps but want to use the standard callback type.
   *
   * @param callback The standard progress callback
   * @return A step-based progress function
   */
  def fromSteps(callback: ProgressCallback): (String, Int, Int) => Unit =
    (msg, current, total) => {
      val fraction = if (total > 0) current.toDouble / total else 0.0
      callback(msg, fraction, 1.0)
    }

  /**
   * Convert a standard callback to a step-based callback.
   *
   * Use this adapter when passing progress to code that expects step counts.
   *
   * @param callback The step-based progress callback
   * @return A standard progress callback
   */
  def toSteps(callback: (String, Int, Int) => Unit): ProgressCallback =
    (msg, current, total) => {
      val currentStep = if (total > 0) (current / total * 100).toInt else 0
      callback(msg, currentStep, 100)
    }

  /**
   * Create a scoped progress reporter for a sub-task.
   *
   * Maps progress from 0.0-1.0 to a specific range within the parent's progress.
   * Useful when a task has multiple sub-tasks that each contribute a portion
   * of the overall progress.
   *
   * @param callback The parent progress callback
   * @param start    Starting fraction within parent (0.0 to 1.0)
   * @param end      Ending fraction within parent (0.0 to 1.0)
   * @return A scoped progress callback for the sub-task
   */
  def scoped(callback: ProgressCallback, start: Double, end: Double): ProgressCallback =
    (msg, current, total) => {
      val fraction = if (total > 0) current / total else 0.0
      val mappedProgress = start + fraction * (end - start)
      callback(msg, mappedProgress, 1.0)
    }

  /**
   * Create a named progress reporter that prefixes messages.
   *
   * Useful when multiple sub-systems report to the same callback.
   *
   * @param callback The underlying progress callback
   * @param prefix   Prefix to add to all messages
   * @return A prefixed progress callback
   */
  def prefixed(callback: ProgressCallback, prefix: String): ProgressCallback =
    (msg, current, total) => callback(s"$prefix: $msg", current, total)

  /**
   * Create a throttled progress reporter that only reports at intervals.
   *
   * Reduces UI update frequency when processing many items.
   *
   * @param callback      The underlying progress callback
   * @param minIntervalMs Minimum milliseconds between updates
   * @return A throttled progress callback
   */
  def throttled(callback: ProgressCallback, minIntervalMs: Long = 100): ProgressCallback = {
    var lastUpdateTime = 0L

    (msg, current, total) => {
      val now = System.currentTimeMillis()
      if (now - lastUpdateTime >= minIntervalMs || current >= total) {
        lastUpdateTime = now
        callback(msg, current, total)
      }
    }
  }

  /**
   * Combine multiple progress callbacks into one.
   *
   * @param callbacks The callbacks to combine
   * @return A callback that invokes all provided callbacks
   */
  def combine(callbacks: ProgressCallback*): ProgressCallback =
    (msg, current, total) => callbacks.foreach(cb => cb(msg, current, total))

  /**
   * Create a progress reporter that logs to the console.
   *
   * Useful for debugging progress reporting.
   *
   * @param prefix Optional prefix for log messages
   * @return A console-logging progress callback
   */
  def console(prefix: String = ""): ProgressCallback = {
    val pfx = if (prefix.nonEmpty) s"[$prefix] " else ""
    (msg, current, total) => {
      val pct = if (total > 0) (current / total * 100).toInt else 0
      println(s"${pfx}Progress: $pct% - $msg")
    }
  }
}
