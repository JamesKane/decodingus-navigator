package com.decodingus.ui.theme

import com.decodingus.config.UserPreferencesService
import scalafx.beans.property.ObjectProperty

/**
 * Application theming system supporting dark and light modes.
 *
 * Usage:
 * {{{
 * // Get current color
 * val bgColor = Theme.current.background
 *
 * // Apply to style
 * style = s"-fx-background-color: ${Theme.current.background};"
 *
 * // Listen for theme changes
 * Theme.currentTheme.onChange { (_, _, _) => updateStyles() }
 * }}}
 */
object Theme {

  /** Available theme modes */
  sealed trait ThemeMode
  case object Dark extends ThemeMode
  case object Light extends ThemeMode

  /** Current theme mode (observable for reactive UI updates) */
  val currentMode: ObjectProperty[ThemeMode] = ObjectProperty(Dark)

  /** Get the current color scheme */
  def current: ColorScheme = currentMode.value match {
    case Dark => DarkColors
    case Light => LightColors
  }

  /** Initialize theme from user preferences */
  def initializeFromPreferences(): Unit = {
    val prefs = UserPreferencesService.load()
    val mode = prefs.theme.getOrElse("dark") match {
      case "light" => Light
      case _ => Dark
    }
    currentMode.value = mode
  }

  /** Set the theme mode and save to preferences */
  def setTheme(mode: ThemeMode): Unit = {
    currentMode.value = mode
    val themeName = mode match {
      case Dark => "dark"
      case Light => "light"
    }
    UserPreferencesService.setTheme(themeName)
  }

  /** Toggle between dark and light modes */
  def toggle(): Unit = {
    currentMode.value match {
      case Dark => setTheme(Light)
      case Light => setTheme(Dark)
    }
  }

  /** Check if currently in dark mode */
  def isDark: Boolean = currentMode.value == Dark
}

/**
 * Color scheme trait defining all theme colors.
 */
trait ColorScheme {
  // Base backgrounds
  def background: String
  def backgroundAlt: String
  def surface: String
  def surfaceAlt: String

  // Text colors
  def textPrimary: String
  def textSecondary: String
  def textMuted: String
  def textDisabled: String

  // Border colors
  def border: String
  def borderLight: String
  def divider: String

  // Interactive states
  def hover: String
  def selected: String
  def selectedHover: String

  // Accent colors (same for both themes)
  def accent: String = "#4a9eff"
  def accentHover: String = "#5aafff"
  def accentPressed: String = "#3a8eef"

  // Semantic colors
  def success: String = "#4ade80"
  def warning: String = "#fbbf24"
  def error: String = "#f87171"
  def danger: String = "#e53935"
  def dangerHover: String = "#f54945"

  // Haplogroup-specific colors
  def ydnaAccent: String = "#4ade80"
  def mtdnaAccent: String = "#60a5fa"

  // Card backgrounds with subtle color tints
  def ydnaCardBg: String
  def mtdnaCardBg: String
}

/**
 * Dark theme color scheme.
 */
object DarkColors extends ColorScheme {
  // Base backgrounds
  val background = "#1e1e1e"
  val backgroundAlt = "#252525"
  val surface = "#2a2a2a"
  val surfaceAlt = "#333333"

  // Text colors
  val textPrimary = "#ffffff"
  val textSecondary = "#e0e0e0"
  val textMuted = "#b0b0b0"
  val textDisabled = "#666666"

  // Border colors
  val border = "#3a3a3a"
  val borderLight = "#4a4a4a"
  val divider = "#3a3a3a"

  // Interactive states
  val hover = "#333333"
  val selected = "#3a5a8a"
  val selectedHover = "#4a6a9a"

  // Haplogroup card backgrounds
  val ydnaCardBg = "#2d3a2d"
  val mtdnaCardBg = "#2d2d3a"
}

/**
 * Light theme color scheme.
 */
object LightColors extends ColorScheme {
  // Base backgrounds
  val background = "#f5f5f5"
  val backgroundAlt = "#eeeeee"
  val surface = "#ffffff"
  val surfaceAlt = "#fafafa"

  // Text colors
  val textPrimary = "#1a1a1a"
  val textSecondary = "#333333"
  val textMuted = "#666666"
  val textDisabled = "#999999"

  // Border colors
  val border = "#d0d0d0"
  val borderLight = "#e0e0e0"
  val divider = "#e0e0e0"

  // Interactive states
  val hover = "#e8e8e8"
  val selected = "#cce5ff"
  val selectedHover = "#b3d9ff"

  // Haplogroup card backgrounds
  val ydnaCardBg = "#e8f5e9"
  val mtdnaCardBg = "#e3f2fd"
}
