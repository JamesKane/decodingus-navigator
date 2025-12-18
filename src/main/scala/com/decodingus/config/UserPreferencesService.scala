package com.decodingus.config

import io.circe.parser.*
import io.circe.syntax.*

import java.nio.file.{Files, Path, Paths, StandardOpenOption}
import java.util.Locale
import scala.util.{Failure, Success, Try}

/**
 * Service for loading and saving user preferences.
 * Stores preferences in ~/.decodingus/config/user_preferences.json
 */
object UserPreferencesService {

  private val CONFIG_DIR = Paths.get(System.getProperty("user.home"), ".decodingus", "config")
  private val CONFIG_FILE = "user_preferences.json"

  private def configFilePath: Path = CONFIG_DIR.resolve(CONFIG_FILE)

  // In-memory cache of the current preferences
  @volatile private var cachedPreferences: Option[UserPreferences] = None

  /**
   * Loads the user preferences from disk.
   * Returns default preferences if file doesn't exist.
   */
  def load(): UserPreferences = {
    cachedPreferences.getOrElse {
      val prefs = loadFromDisk()
      cachedPreferences = Some(prefs)
      prefs
    }
  }

  /**
   * Forces a reload from disk, bypassing the cache.
   */
  def reload(): UserPreferences = {
    val prefs = loadFromDisk()
    cachedPreferences = Some(prefs)
    prefs
  }

  /**
   * Saves the user preferences to disk.
   */
  def save(prefs: UserPreferences): Either[String, Unit] = {
    ensureConfigDir()
    val jsonString = prefs.asJson.spaces2

    Try {
      Files.writeString(
        configFilePath,
        jsonString,
        StandardOpenOption.CREATE,
        StandardOpenOption.TRUNCATE_EXISTING
      )
    } match {
      case Success(_) =>
        cachedPreferences = Some(prefs)
        println(s"[UserPreferencesService] Saved preferences to $configFilePath")
        Right(())
      case Failure(e) =>
        Left(s"Failed to save user preferences: ${e.getMessage}")
    }
  }

  /**
   * Gets the current Y-DNA tree provider.
   */
  def getYdnaTreeProvider: String = load().ydnaTreeProvider

  /**
   * Gets the current MT-DNA tree provider.
   */
  def getMtdnaTreeProvider: String = load().mtdnaTreeProvider

  /**
   * Sets the Y-DNA tree provider.
   */
  def setYdnaTreeProvider(provider: String): Either[String, Unit] = {
    if (!UserPreferences.ValidTreeProviders.contains(provider.toLowerCase)) {
      Left(s"Invalid tree provider: $provider")
    } else {
      val currentPrefs = load()
      save(currentPrefs.copy(ydnaTreeProvider = provider.toLowerCase))
    }
  }

  /**
   * Sets the MT-DNA tree provider.
   */
  def setMtdnaTreeProvider(provider: String): Either[String, Unit] = {
    if (!UserPreferences.ValidTreeProviders.contains(provider.toLowerCase)) {
      Left(s"Invalid tree provider: $provider")
    } else {
      val currentPrefs = load()
      save(currentPrefs.copy(mtdnaTreeProvider = provider.toLowerCase))
    }
  }

  /**
   * Gets the saved locale, or None if using system default.
   */
  def getLocale: Option[Locale] = {
    load().locale.map(Locale.forLanguageTag)
  }

  /**
   * Sets the UI locale.
   * Pass None to use system default.
   */
  def setLocale(locale: Option[Locale]): Either[String, Unit] = {
    val currentPrefs = load()
    save(currentPrefs.copy(locale = locale.map(_.toLanguageTag)))
  }

  /**
   * Sets the UI locale by Locale object.
   */
  def setLocale(locale: Locale): Either[String, Unit] = {
    setLocale(Some(locale))
  }

  /**
   * Gets the saved theme, or None if using default (dark).
   */
  def getTheme: Option[String] = {
    load().theme
  }

  /**
   * Sets the UI theme ("dark" or "light").
   */
  def setTheme(theme: String): Either[String, Unit] = {
    val validThemes = Set("dark", "light")
    if (!validThemes.contains(theme.toLowerCase)) {
      Left(s"Invalid theme: $theme. Valid options: dark, light")
    } else {
      val currentPrefs = load()
      save(currentPrefs.copy(theme = Some(theme.toLowerCase)))
    }
  }

  /**
   * Gets saved dialog size for a given dialog ID.
   */
  def getDialogSize(dialogId: String): Option[DialogSize] = {
    load().dialogSizes.get(dialogId)
  }

  /**
   * Saves dialog size for a given dialog ID.
   */
  def setDialogSize(dialogId: String, width: Double, height: Double): Either[String, Unit] = {
    val currentPrefs = load()
    val updatedSizes = currentPrefs.dialogSizes + (dialogId -> DialogSize(width, height))
    save(currentPrefs.copy(dialogSizes = updatedSizes))
  }

  private def loadFromDisk(): UserPreferences = {
    if (!Files.exists(configFilePath)) {
      println(s"[UserPreferencesService] Preferences file not found, using defaults")
      return UserPreferences.default
    }

    Try(Files.readString(configFilePath)) match {
      case Success(jsonString) =>
        parse(jsonString).flatMap(_.as[UserPreferences]) match {
          case Right(prefs) =>
            println(s"[UserPreferencesService] Loaded user preferences")
            prefs
          case Left(error) =>
            println(s"[UserPreferencesService] Failed to parse preferences: ${error.getMessage}, using defaults")
            UserPreferences.default
        }
      case Failure(e) =>
        println(s"[UserPreferencesService] Failed to read preferences file: ${e.getMessage}, using defaults")
        UserPreferences.default
    }
  }

  private def ensureConfigDir(): Unit = {
    if (!Files.exists(CONFIG_DIR)) {
      Try(Files.createDirectories(CONFIG_DIR)) match {
        case Failure(e) =>
          println(s"[UserPreferencesService] Warning: Could not create config directory: ${e.getMessage}")
        case _ =>
      }
    }
  }
}
