package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.model.{DefiningHaplogroup, NamedVariant, VariantAliases, VariantCoordinate}
import io.circe.generic.auto.*
import io.circe.parser.decode
import sttp.client3.{HttpURLConnectionBackend, asByteArray, basicRequest}

import java.io.*
import java.nio.file.{Files, Path, Paths}
import java.time.{Duration, Instant}
import java.util.zip.GZIPInputStream
import scala.collection.mutable
import scala.io.Source
import scala.util.{Try, Using}

/**
 * Cache for named variants from the Decoding Us variant database.
 * Downloads and caches the daily export file (gzipped JSONL) and provides
 * efficient lookups by position, name, or rsId.
 */
class NamedVariantCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "variants")
    try {
      Files.createDirectories(dir)
    } catch {
      case _: Exception =>
        Files.createTempDirectory("variant-cache")
    }
    dir
  }

  private val exportUrl = "https://decoding-us.com/api/v1/variants/export"
  private val cacheFileName = "variants-export.jsonl.gz"
  private val maxCacheAge = Duration.ofHours(24)

  // In-memory index for fast lookups (populated on first access)
  @volatile private var variantsByPosition: Map[String, Map[Long, NamedVariant]] = Map.empty
  @volatile private var variantsByName: Map[String, NamedVariant] = Map.empty
  @volatile private var variantsByRsId: Map[String, NamedVariant] = Map.empty
  @volatile private var loaded: Boolean = false

  private val cacheFile: File = cacheDir.resolve(cacheFileName).toFile

  /**
   * Check if the cache needs to be refreshed (older than 24 hours or missing).
   */
  def needsRefresh: Boolean = {
    if (!cacheFile.exists()) {
      true
    } else {
      val lastModified = Instant.ofEpochMilli(cacheFile.lastModified())
      Duration.between(lastModified, Instant.now()).compareTo(maxCacheAge) > 0
    }
  }

  /**
   * Download and cache the variant export file.
   * Returns Right on success, Left with error message on failure.
   */
  def refresh(progressCallback: String => Unit = _ => ()): Either[String, Unit] = {
    progressCallback("Downloading named variants database...")

    val backend = HttpURLConnectionBackend()
    val response = basicRequest
      .get(sttp.model.Uri.unsafeParse(exportUrl))
      .response(asByteArray)
      .send(backend)

    response.body match {
      case Right(bytes) =>
        try {
          progressCallback(s"Saving ${bytes.length / 1024 / 1024}MB to cache...")
          Using.resource(new FileOutputStream(cacheFile)) { fos =>
            fos.write(bytes)
          }
          // Clear in-memory index to force reload
          loaded = false
          variantsByPosition = Map.empty
          variantsByName = Map.empty
          variantsByRsId = Map.empty
          progressCallback("Named variants database updated successfully")
          Right(())
        } catch {
          case e: Exception =>
            Left(s"Failed to save variant cache: ${e.getMessage}")
        }
      case Left(error) =>
        Left(s"Failed to download variant export: $error")
    }
  }

  /**
   * Ensure the cache is loaded into memory.
   * Downloads if missing or stale, then parses into lookup indices.
   */
  def ensureLoaded(progressCallback: String => Unit = _ => ()): Either[String, Unit] = {
    if (loaded) {
      Right(())
    } else {
      // Download if needed
      if (needsRefresh) {
        refresh(progressCallback) match {
          case Left(error) => return Left(error)
          case Right(_) => // continue to load
        }
      }

      // Parse the gzipped JSONL file
      progressCallback("Loading named variants into memory...")
      loadFromCache(progressCallback)
    }
  }

  /**
   * Parse the cached gzipped JSONL file and build indices.
   */
  private def loadFromCache(progressCallback: String => Unit): Either[String, Unit] = {
    if (!cacheFile.exists()) {
      return Left("Variant cache file not found. Call refresh() first.")
    }

    try {
      val positionIndex = mutable.Map[String, mutable.Map[Long, NamedVariant]]()
      val nameIndex = mutable.Map[String, NamedVariant]()
      val rsIdIndex = mutable.Map[String, NamedVariant]()
      var count = 0

      Using.resource(new BufferedReader(new InputStreamReader(
        new GZIPInputStream(new FileInputStream(cacheFile))
      ))) { reader =>
        var line = reader.readLine()
        while (line != null) {
          decode[NamedVariant](line) match {
            case Right(variant) =>
              // Index by position for each build
              variant.coordinates.foreach { case (build, coord) =>
                val buildMap = positionIndex.getOrElseUpdate(build, mutable.Map.empty)
                buildMap(coord.position.toLong) = variant
              }

              // Index by canonical name
              variant.canonicalName.foreach { name =>
                nameIndex(name.toUpperCase) = variant
              }

              // Index by all common names
              variant.aliases.commonNames.foreach { name =>
                nameIndex(name.toUpperCase) = variant
              }

              // Index by rsIds
              variant.aliases.rsIds.foreach { rsId =>
                rsIdIndex(rsId.toLowerCase) = variant
              }

              count += 1
              if (count % 50000 == 0) {
                progressCallback(s"Loaded $count variants...")
              }
            case Left(error) =>
              // Log but continue - some lines might be malformed
              System.err.println(s"Failed to parse variant: ${error.getMessage}")
          }
          line = reader.readLine()
        }
      }

      variantsByPosition = positionIndex.map { case (build, map) => build -> map.toMap }.toMap
      variantsByName = nameIndex.toMap
      variantsByRsId = rsIdIndex.toMap
      loaded = true

      progressCallback(s"Loaded $count named variants")
      Right(())
    } catch {
      case e: Exception =>
        Left(s"Failed to load variant cache: ${e.getMessage}")
    }
  }

  /**
   * Look up a variant by genomic position.
   *
   * @param build    Reference build (e.g., "GRCh38")
   * @param position Genomic position
   * @return The variant if found
   */
  def getByPosition(build: String, position: Long): Option[NamedVariant] = {
    ensureLoaded()
    variantsByPosition.get(build).flatMap(_.get(position))
  }

  /**
   * Look up a variant by name (canonical or alias).
   *
   * @param name Variant name (case-insensitive)
   * @return The variant if found
   */
  def getByName(name: String): Option[NamedVariant] = {
    ensureLoaded()
    variantsByName.get(name.toUpperCase)
  }

  /**
   * Look up a variant by rsId.
   *
   * @param rsId The rsId (e.g., "rs12345")
   * @return The variant if found
   */
  def getByRsId(rsId: String): Option[NamedVariant] = {
    ensureLoaded()
    variantsByRsId.get(rsId.toLowerCase)
  }

  /**
   * Get all variants defining a specific haplogroup.
   *
   * @param haplogroupName The haplogroup name
   * @return List of variants that define this haplogroup
   */
  def getVariantsForHaplogroup(haplogroupName: String): List[NamedVariant] = {
    ensureLoaded()
    variantsByName.values.filter { variant =>
      variant.definingHaplogroup.exists(_.haplogroupName == haplogroupName)
    }.toList
  }

  /**
   * Check if the cache is loaded and ready for lookups.
   */
  def isLoaded: Boolean = loaded

  /**
   * Get the number of cached variants.
   */
  def variantCount: Int = variantsByName.size

  /**
   * Get the cache file location.
   */
  def getCacheFile: File = cacheFile

  /**
   * Clear the in-memory cache (does not delete the file).
   */
  def clearMemoryCache(): Unit = {
    loaded = false
    variantsByPosition = Map.empty
    variantsByName = Map.empty
    variantsByRsId = Map.empty
  }
}

object NamedVariantCache {
  // Singleton instance for convenience
  private lazy val instance = new NamedVariantCache()

  def apply(): NamedVariantCache = instance
}
