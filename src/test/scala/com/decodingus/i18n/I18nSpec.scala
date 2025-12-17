package com.decodingus.i18n

import munit.FunSuite
import java.util.Locale

class I18nSpec extends FunSuite:

  // Reset to English before each test
  override def beforeEach(context: BeforeEach): Unit =
    I18n.setLocale(Locale.ENGLISH)

  // ==========================================================================
  // Basic Translation Tests
  // ==========================================================================

  test("t returns translated string for valid key") {
    val result = I18n.t("nav.dashboard")
    assertEquals(result, "Dashboard")
  }

  test("t returns !key! for missing key") {
    val result = I18n.t("nonexistent.key")
    assertEquals(result, "!nonexistent.key!")
  }

  test("t with parameters substitutes values") {
    val result = I18n.t("dashboard.subjects.count", 47)
    assertEquals(result, "47 Subjects")
  }

  test("t with multiple parameters substitutes all values") {
    val result = I18n.t("ibd.closest_match", "Bob Smith", "1,847 cM")
    assertEquals(result, "Closest: Bob Smith (1,847 cM)")
  }

  // ==========================================================================
  // Locale Switching Tests
  // ==========================================================================

  test("setLocale changes current locale") {
    I18n.setLocale(Locale.GERMAN)
    assertEquals(I18n.currentLocale.value, Locale.GERMAN)
  }

  test("isRTL returns false for LTR locales") {
    I18n.setLocale(Locale.ENGLISH)
    assert(!I18n.isRTL)

    I18n.setLocale(Locale.GERMAN)
    assert(!I18n.isRTL)
  }

  test("isRTL returns true for RTL locales") {
    I18n.setLocale(new Locale("ar"))
    assert(I18n.isRTL)

    I18n.setLocale(new Locale("he"))
    assert(I18n.isRTL)
  }

  // ==========================================================================
  // Supported Locales Tests
  // ==========================================================================

  test("supportedLocales includes English") {
    assert(I18n.supportedLocales.contains(Locale.ENGLISH))
  }

  test("supportedLocales includes German") {
    assert(I18n.supportedLocales.contains(Locale.GERMAN))
  }

  test("getLocaleDisplayName returns localized name") {
    assertEquals(I18n.getLocaleDisplayName(Locale.ENGLISH), "English")
    assertEquals(I18n.getLocaleDisplayName(Locale.GERMAN), "Deutsch")
  }

  // ==========================================================================
  // Reactive Binding Tests
  // ==========================================================================

  test("bind creates StringProperty with translated value") {
    val prop = I18n.bind("nav.dashboard")
    assertEquals(prop.value, "Dashboard")
  }

  test("bind with args creates parameterized StringProperty") {
    val count = 42
    val prop = I18n.bind("dashboard.subjects.count", () => Seq(count))
    assertEquals(prop.value, "42 Subjects")
  }

  // ==========================================================================
  // Message Coverage Tests
  // ==========================================================================

  test("navigation keys exist") {
    assertNotEquals(I18n.t("nav.dashboard"), "!nav.dashboard!")
    assertNotEquals(I18n.t("nav.subjects"), "!nav.subjects!")
    assertNotEquals(I18n.t("nav.projects"), "!nav.projects!")
  }

  test("action keys exist") {
    assertNotEquals(I18n.t("action.edit"), "!action.edit!")
    assertNotEquals(I18n.t("action.delete"), "!action.delete!")
    assertNotEquals(I18n.t("action.save"), "!action.save!")
    assertNotEquals(I18n.t("action.cancel"), "!action.cancel!")
  }

  test("haplogroup keys exist") {
    assertNotEquals(I18n.t("haplogroup.ydna.title"), "!haplogroup.ydna.title!")
    assertNotEquals(I18n.t("haplogroup.mtdna.title"), "!haplogroup.mtdna.title!")
    assertNotEquals(I18n.t("haplogroup.derived"), "!haplogroup.derived!")
    assertNotEquals(I18n.t("haplogroup.ancestral"), "!haplogroup.ancestral!")
  }

  test("status keys exist") {
    assertNotEquals(I18n.t("status.complete"), "!status.complete!")
    assertNotEquals(I18n.t("status.pending"), "!status.pending!")
    assertNotEquals(I18n.t("status.error"), "!status.error!")
  }
