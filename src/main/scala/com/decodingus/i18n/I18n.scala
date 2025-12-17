package com.decodingus.i18n

import java.text.MessageFormat
import java.util.{Locale, MissingResourceException, ResourceBundle}
import scalafx.beans.property.{ObjectProperty, StringProperty, ReadOnlyStringProperty}

/**
 * Internationalization support for the application.
 *
 * Provides locale-aware message retrieval with support for:
 * - Simple string lookup via `t(key)`
 * - Parameterized messages via `t(key, args...)`
 * - Reactive bindings that update when locale changes via `bind(key)`
 * - RTL (right-to-left) locale detection
 *
 * Usage:
 * {{{
 * import com.decodingus.i18n.I18n.{t, bind}
 *
 * // One-time lookup
 * val message = t("nav.dashboard")
 *
 * // With parameters
 * val count = t("dashboard.subjects.count", 47)
 *
 * // Reactive binding for UI
 * label.text <== bind("nav.dashboard")
 * }}}
 */
object I18n {

  /** Observable locale property - UI components can bind to this for reactive updates */
  val currentLocale: ObjectProperty[Locale] = ObjectProperty(Locale.getDefault)

  // Cached bundle reference, reloaded on locale change
  private var cachedBundle: ResourceBundle = loadBundle()

  private def loadBundle(): ResourceBundle = {
    ResourceBundle.getBundle("i18n.messages", currentLocale.value)
  }

  // Reload bundle when locale changes
  currentLocale.onChange { (_, _, _) =>
    cachedBundle = loadBundle()
  }

  /**
   * Get a translated string for the given key.
   * Returns "!key!" if the key is missing (makes missing translations visible in UI).
   */
  def t(key: String): String = {
    try {
      cachedBundle.getString(key)
    } catch {
      case _: MissingResourceException =>
        System.err.println(s"[I18n] Missing key: $key")
        s"!$key!"
    }
  }

  /**
   * Get a translated string with ICU MessageFormat parameters.
   *
   * Example:
   * {{{
   * // messages.properties: dashboard.subjects.count={0} Subjects
   * t("dashboard.subjects.count", 47) // -> "47 Subjects"
   * }}}
   */
  def t(key: String, args: Any*): String = {
    val pattern = t(key)
    if (args.isEmpty) pattern
    else MessageFormat.format(pattern, args.map(_.asInstanceOf[AnyRef]): _*)
  }

  /**
   * Create a reactive StringProperty that updates when the locale changes.
   * Use for static labels without parameters.
   */
  def bind(key: String): StringProperty = {
    val prop = StringProperty(t(key))
    currentLocale.onChange { (_, _, _) =>
      prop.value = t(key)
    }
    prop
  }

  /**
   * Create a reactive binding with dynamic parameters.
   * The args function is re-evaluated on each locale change.
   *
   * Example:
   * {{{
   * label.text <== bind("subjects.selected", () => Seq(selectedCount))
   * }}}
   */
  def bind(key: String, args: () => Seq[Any]): StringProperty = {
    val prop = StringProperty(t(key, args(): _*))
    currentLocale.onChange { (_, _, _) =>
      prop.value = t(key, args(): _*)
    }
    prop
  }

  /**
   * Create a read-only binding (useful for exposing to other components).
   */
  def bindReadOnly(key: String): ReadOnlyStringProperty = bind(key)

  /**
   * Check if the current locale uses right-to-left text direction.
   */
  def isRTL: Boolean = {
    val rtlLanguages = Set("ar", "he", "fa", "ur")
    rtlLanguages.contains(currentLocale.value.getLanguage)
  }

  /**
   * Switch to a new locale at runtime.
   * All reactive bindings will automatically update.
   */
  def setLocale(locale: Locale): Unit = {
    currentLocale.value = locale
    Locale.setDefault(locale)
  }

  /**
   * Get the list of supported locales for the language picker UI.
   */
  def supportedLocales: Seq[Locale] = Seq(
    Locale.ENGLISH,
    Locale.GERMAN,
    new Locale("es"),
    Locale.FRENCH,
    new Locale("pt", "BR"),
    Locale.JAPANESE,
    Locale.SIMPLIFIED_CHINESE
    // RTL locales can be added when support is implemented:
    // new Locale("ar"),
    // new Locale("he")
  )

  /**
   * Get the display name of a locale in its own language.
   * Example: Locale.GERMAN -> "Deutsch"
   */
  def getLocaleDisplayName(locale: Locale): String = {
    locale.getDisplayLanguage(locale).capitalize
  }

  /**
   * Initialize from saved preferences.
   * Call this at application startup.
   */
  def initializeFromPreferences(): Unit = {
    import com.decodingus.config.UserPreferencesService
    UserPreferencesService.getLocale match {
      case Some(savedLocale) =>
        setLocale(savedLocale)
        println(s"[I18n] Initialized locale from preferences: ${savedLocale.toLanguageTag}")
      case None =>
        println(s"[I18n] Using system default locale: ${Locale.getDefault.toLanguageTag}")
    }
  }

  /**
   * Save the current locale to preferences.
   */
  def saveLocalePreference(): Unit = {
    import com.decodingus.config.UserPreferencesService
    UserPreferencesService.setLocale(currentLocale.value) match {
      case Right(_) =>
        println(s"[I18n] Saved locale preference: ${currentLocale.value.toLanguageTag}")
      case Left(error) =>
        System.err.println(s"[I18n] Failed to save locale preference: $error")
    }
  }
}
