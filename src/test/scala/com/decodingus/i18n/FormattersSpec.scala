package com.decodingus.i18n

import munit.FunSuite
import java.time.LocalDate
import java.util.Locale

class FormattersSpec extends FunSuite:

  // Reset to English before each test for consistent results
  override def beforeEach(context: BeforeEach): Unit =
    I18n.setLocale(Locale.US)

  // ==========================================================================
  // Number Formatting Tests
  // ==========================================================================

  test("formatNumber formats integers with thousands separators") {
    assertEquals(Formatters.formatNumber(1234567), "1,234,567")
  }

  test("formatNumber formats small numbers without separators") {
    assertEquals(Formatters.formatNumber(123), "123")
  }

  test("formatDecimal formats with specified precision") {
    assertEquals(Formatters.formatDecimal(32.456, 1), "32.5")
    assertEquals(Formatters.formatDecimal(32.456, 2), "32.46")
  }

  test("formatPercent formats ratio as percentage") {
    val result = Formatters.formatPercent(0.982)
    assert(result.contains("98"), s"Expected percentage containing '98', got: $result")
  }

  test("formatPercentValue formats percentage value directly") {
    assertEquals(Formatters.formatPercentValue(98.2, 1), "98.2%")
  }

  // ==========================================================================
  // Locale-Specific Number Formatting
  // ==========================================================================

  test("formatNumber uses locale-specific separators for German") {
    I18n.setLocale(Locale.GERMAN)
    assertEquals(Formatters.formatNumber(1234567), "1.234.567")
  }

  test("formatDecimal uses locale-specific decimal separator for German") {
    I18n.setLocale(Locale.GERMAN)
    assertEquals(Formatters.formatDecimal(32.4, 1), "32,4")
  }

  // ==========================================================================
  // Date Formatting Tests
  // ==========================================================================

  test("formatDate returns medium-style date") {
    val date = LocalDate.of(2024, 12, 15)
    val result = Formatters.formatDate(date)
    // US locale should produce something like "Dec 15, 2024"
    assert(result.contains("2024"), s"Expected date containing '2024', got: $result")
    assert(result.contains("15"), s"Expected date containing '15', got: $result")
  }

  test("formatDateShort returns short-style date") {
    val date = LocalDate.of(2024, 12, 15)
    val result = Formatters.formatDateShort(date)
    // Should be shorter than medium format
    assert(result.length < 15, s"Expected short date, got: $result")
  }

  // ==========================================================================
  // Relative Time Formatting Tests
  // ==========================================================================

  test("formatRelativeTime returns 'Today' for today's date") {
    val today = LocalDate.now()
    assertEquals(Formatters.formatRelativeTime(today), "Today")
  }

  test("formatRelativeTime returns 'Yesterday' for yesterday") {
    val yesterday = LocalDate.now().minusDays(1)
    assertEquals(Formatters.formatRelativeTime(yesterday), "Yesterday")
  }

  test("formatRelativeTime returns days ago for recent dates") {
    val threeDaysAgo = LocalDate.now().minusDays(3)
    val result = Formatters.formatRelativeTime(threeDaysAgo)
    assert(result.contains("3"), s"Expected '3' in result, got: $result")
  }

  // ==========================================================================
  // Domain-Specific Formatters (Genomics)
  // ==========================================================================

  test("formatCoverage appends 'x' suffix") {
    assertEquals(Formatters.formatCoverage(32.4), "32.4x")
  }

  test("formatCentimorgans formats with cM unit") {
    assertEquals(Formatters.formatCentimorgans(1847), "1,847 cM")
  }

  test("formatCentimorgans handles decimal values") {
    assertEquals(Formatters.formatCentimorgans(0.75), "0.75 cM")
  }

  test("formatMegabases formats with Mb unit") {
    assertEquals(Formatters.formatMegabases(2.5), "2.5 Mb")
  }

  test("formatBasePairs uses appropriate unit for large values") {
    assertEquals(Formatters.formatBasePairs(500), "500 bp")
    assert(Formatters.formatBasePairs(5000).contains("Kb"))
    assert(Formatters.formatBasePairs(5000000).contains("Mb"))
    assert(Formatters.formatBasePairs(5000000000L).contains("Gb"))
  }

  test("formatSnpCounts formats derived and ancestral counts") {
    assertEquals(Formatters.formatSnpCounts(847, 12), "847 / 12")
  }

  test("formatPosition formats chromosome position") {
    assertEquals(Formatters.formatPosition("chrY", 2781479), "chrY:2,781,479")
  }

  test("formatInterval formats genomic interval") {
    assertEquals(
      Formatters.formatInterval("chrY", 2781479, 2781580),
      "chrY:2,781,479-2,781,580"
    )
  }
