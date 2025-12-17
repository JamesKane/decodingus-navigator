package com.decodingus.i18n

import java.text.NumberFormat
import java.time.{Instant, LocalDate, LocalDateTime, ZoneId}
import java.time.format.{DateTimeFormatter, FormatStyle}
import java.time.temporal.ChronoUnit

/**
 * Locale-aware formatters for numbers, dates, and domain-specific values.
 *
 * All formatters use the current locale from [[I18n.currentLocale]] to ensure
 * consistent formatting throughout the application.
 *
 * Usage:
 * {{{
 * import com.decodingus.i18n.Formatters._
 *
 * formatNumber(1234567)      // "1,234,567" (en) or "1.234.567" (de)
 * formatPercent(0.982)       // "98.2%" (en) or "98,2 %" (fr)
 * formatCoverage(32.4)       // "32.4x"
 * formatCentimorgans(1847)   // "1,847 cM"
 * }}}
 */
object Formatters {

  /**
   * Format an integer or long with locale-appropriate thousands separators.
   * Example: 1234567 -> "1,234,567" (en) or "1.234.567" (de)
   */
  def formatNumber(n: Number): String = {
    NumberFormat.getNumberInstance(I18n.currentLocale.value).format(n)
  }

  /**
   * Format a decimal number with specified precision.
   *
   * @param d        the decimal value
   * @param decimals number of decimal places (default: 1)
   */
  def formatDecimal(d: Double, decimals: Int = 1): String = {
    val formatter = NumberFormat.getNumberInstance(I18n.currentLocale.value)
    formatter.setMinimumFractionDigits(decimals)
    formatter.setMaximumFractionDigits(decimals)
    formatter.format(d)
  }

  /**
   * Format a value as a percentage.
   * Example: 0.982 -> "98.2%" (en) or "98,2 %" (fr)
   *
   * @param d value between 0.0 and 1.0
   */
  def formatPercent(d: Double): String = {
    NumberFormat.getPercentInstance(I18n.currentLocale.value).format(d)
  }

  /**
   * Format a percentage from a value already in percent form (0-100).
   * Example: 98.2 -> "98.2%"
   */
  def formatPercentValue(d: Double, decimals: Int = 1): String = {
    s"${formatDecimal(d, decimals)}%"
  }

  /**
   * Format a date in medium style.
   * Example: "Dec 15, 2024" (en) or "15. Dez. 2024" (de)
   */
  def formatDate(date: LocalDate): String = {
    DateTimeFormatter
      .ofLocalizedDate(FormatStyle.MEDIUM)
      .withLocale(I18n.currentLocale.value)
      .format(date)
  }

  /**
   * Format a date in short style.
   * Example: "12/15/24" (en-US) or "15.12.24" (de)
   */
  def formatDateShort(date: LocalDate): String = {
    DateTimeFormatter
      .ofLocalizedDate(FormatStyle.SHORT)
      .withLocale(I18n.currentLocale.value)
      .format(date)
  }

  /**
   * Format a relative time string (e.g., "2h ago", "Yesterday").
   * Returns a localized string using i18n keys.
   */
  def formatRelativeTime(date: LocalDate): String = {
    val today = LocalDate.now()
    val daysBetween = ChronoUnit.DAYS.between(date, today)

    daysBetween match {
      case 0 => I18n.t("time.today")
      case 1 => I18n.t("time.yesterday")
      case n if n < 7 => I18n.t("time.ago.days", n)
      case _ => formatDate(date)
    }
  }

  /**
   * Format a date and time in medium style.
   * Example: "Dec 15, 2024, 2:30 PM" (en) or "15. Dez. 2024, 14:30" (de)
   */
  def formatDateTime(dateTime: LocalDateTime): String = {
    DateTimeFormatter
      .ofLocalizedDateTime(FormatStyle.MEDIUM, FormatStyle.SHORT)
      .withLocale(I18n.currentLocale.value)
      .format(dateTime)
  }

  /**
   * Format an Instant as local date/time.
   * Example: "Dec 15, 2024, 2:30 PM" (en)
   */
  def formatInstant(instant: Instant): String = {
    val localDateTime = LocalDateTime.ofInstant(instant, ZoneId.systemDefault())
    formatDateTime(localDateTime)
  }

  // ===========================================================================
  // Domain-Specific Formatters (Genomics)
  // ===========================================================================

  /**
   * Format sequencing coverage with 'x' suffix.
   * The number formatting is locale-specific, but 'x' is universal.
   * Example: 32.4 -> "32.4x" (en) or "32,4x" (de)
   */
  def formatCoverage(coverage: Double): String = {
    s"${formatDecimal(coverage)}x"
  }

  /**
   * Format a centimorgan value with unit.
   * Example: 1847 -> "1,847 cM"
   */
  def formatCentimorgans(cm: Double): String = {
    if (cm >= 1.0) {
      s"${formatNumber(cm.toLong)} cM"
    } else {
      s"${formatDecimal(cm, 2)} cM"
    }
  }

  /**
   * Format a centimorgan value (integer).
   */
  def formatCentimorgans(cm: Int): String = {
    s"${formatNumber(cm)} cM"
  }

  /**
   * Format megabases with unit.
   * Example: 2.5 -> "2.5 Mb"
   */
  def formatMegabases(mb: Double): String = {
    s"${formatDecimal(mb)} Mb"
  }

  /**
   * Format base pairs with appropriate unit (bp, Kb, Mb, Gb).
   */
  def formatBasePairs(bp: Long): String = {
    if (bp >= 1_000_000_000L) {
      s"${formatDecimal(bp / 1_000_000_000.0)} Gb"
    } else if (bp >= 1_000_000L) {
      s"${formatDecimal(bp / 1_000_000.0)} Mb"
    } else if (bp >= 1_000L) {
      s"${formatDecimal(bp / 1_000.0)} Kb"
    } else {
      s"${formatNumber(bp)} bp"
    }
  }

  /**
   * Format SNP counts (derived/ancestral).
   * Example: (847, 12) -> "847 / 12"
   */
  def formatSnpCounts(derived: Int, ancestral: Int): String = {
    s"${formatNumber(derived)} / ${formatNumber(ancestral)}"
  }

  /**
   * Format a chromosome position.
   * Example: ("chrY", 2781479) -> "chrY:2,781,479"
   */
  def formatPosition(contig: String, position: Long): String = {
    s"$contig:${formatNumber(position)}"
  }

  /**
   * Format a genomic interval.
   * Example: ("chrY", 2781479, 2781580) -> "chrY:2,781,479-2,781,580"
   */
  def formatInterval(contig: String, start: Long, end: Long): String = {
    s"$contig:${formatNumber(start)}-${formatNumber(end)}"
  }
}
