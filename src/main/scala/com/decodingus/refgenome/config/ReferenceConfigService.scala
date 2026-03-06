package com.decodingus.refgenome.config

import io.circe.parser.*
import io.circe.syntax.*

import java.nio.file.{Files, Path, Paths, StandardOpenOption}
import scala.util.{Failure, Success, Try}

/**
 * Service for loading and saving reference genome configuration.
 * Stores config in ~/.decodingus/config/reference_config.json
 */
object ReferenceConfigService {

  private val CONFIG_DIR = Paths.get(System.getProperty("user.home"), ".decodingus", "config")
  private val CONFIG_FILE = "reference_config.json"

  private def configFilePath: Path = CONFIG_DIR.resolve(CONFIG_FILE)

  // In-memory cache of the current config
  @volatile private var cachedConfig: Option[ReferenceConfig] = None

  /**
   * Loads the reference configuration from disk.
   * Returns default config if file doesn't exist.
   */
  def load(): ReferenceConfig = {
    cachedConfig.getOrElse {
      val config = loadFromDisk()
      cachedConfig = Some(config)
      config
    }
  }

  /**
   * Forces a reload from disk, bypassing the cache.
   */
  def reload(): ReferenceConfig = {
    val config = loadFromDisk()
    cachedConfig = Some(config)
    config
  }

  /**
   * Saves the reference configuration to disk.
   */
  def save(config: ReferenceConfig): Either[String, Unit] = {
    ensureConfigDir()
    val jsonString = config.asJson.spaces2

    Try {
      Files.writeString(
        configFilePath,
        jsonString,
        StandardOpenOption.CREATE,
        StandardOpenOption.TRUNCATE_EXISTING
      )
    } match {
      case Success(_) =>
        cachedConfig = Some(config)
        println(s"[ReferenceConfigService] Saved config to $configFilePath")
        Right(())
      case Failure(e) =>
        Left(s"Failed to save reference config: ${e.getMessage}")
    }
  }

  /**
   * Updates the config for a specific reference build.
   */
  def updateReference(buildConfig: ReferenceGenomeConfig): Either[String, Unit] = {
    val currentConfig = load()
    val updatedConfig = currentConfig.withReference(buildConfig)
    save(updatedConfig)
  }

  /**
   * Sets whether to prompt before downloading references.
   */
  def setPromptBeforeDownload(prompt: Boolean): Either[String, Unit] = {
    val currentConfig = load()
    save(currentConfig.copy(promptBeforeDownload = prompt))
  }

  /**
   * Gets the resolved path for a reference build.
   * Checks in order:
   * 1. User-specified local path (if valid)
   * 2. Default cache directory
   *
   * Returns None if reference is not available locally.
   */
  def getReferencePath(build: String): Option[Path] = {
    val config = load()
    val buildConfig = config.getOrDefault(build)

    // First check user-specified local path
    buildConfig.getValidLocalPath.orElse {
      // Then check default cache
      val cachePath = config.getCacheDir.resolve(s"$build.fa.gz")
      if (Files.exists(cachePath) && Files.isRegularFile(cachePath)) {
        Some(cachePath)
      } else {
        None
      }
    }
  }

  /**
   * Checks if a reference is available locally (either user path or cache).
   */
  def isReferenceAvailable(build: String): Boolean = {
    getReferencePath(build).isDefined
  }

  /**
   * Gets the cache directory path for storing downloaded references.
   */
  def getCacheDir: Path = load().getCacheDir

  private def loadFromDisk(): ReferenceConfig = {
    if (!Files.exists(configFilePath)) {
      println(s"[ReferenceConfigService] Config file not found, using defaults")
      return ReferenceConfig.default
    }

    Try(Files.readString(configFilePath)) match {
      case Success(jsonString) =>
        parse(jsonString).flatMap(_.as[ReferenceConfig]) match {
          case Right(config) =>
            println(s"[ReferenceConfigService] Loaded config with ${config.references.size} references")
            config
          case Left(error) =>
            println(s"[ReferenceConfigService] Failed to parse config: ${error.getMessage}, using defaults")
            ReferenceConfig.default
        }
      case Failure(e) =>
        println(s"[ReferenceConfigService] Failed to read config file: ${e.getMessage}, using defaults")
        ReferenceConfig.default
    }
  }

  private def ensureConfigDir(): Unit = {
    if (!Files.exists(CONFIG_DIR)) {
      Try(Files.createDirectories(CONFIG_DIR)) match {
        case Failure(e) =>
          println(s"[ReferenceConfigService] Warning: Could not create config directory: ${e.getMessage}")
        case _ =>
      }
    }
  }
}
