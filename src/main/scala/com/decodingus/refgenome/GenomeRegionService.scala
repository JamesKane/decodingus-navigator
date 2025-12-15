package com.decodingus.refgenome

import com.decodingus.client.DecodingUsClient
import com.decodingus.config.FeatureToggles
import com.decodingus.refgenome.model.{GenomeRegions, GenomeRegionsCodecs}
import io.circe.parser.decode
import io.circe.syntax.*

import java.io.IOException
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path, Paths, StandardOpenOption}
import java.time.{Duration, Instant}
import scala.concurrent.{ExecutionContext, Future}
import scala.io.Source
import scala.util.{Try, Using}

/**
 * Service for fetching and caching genome region metadata.
 *
 * This service provides centralized access to genomic region annotations
 * (centromeres, telomeres, cytobands, Y-specific regions) with:
 *
 * - Feature toggle gating (disabled by default until API is deployed)
 * - Local caching with configurable expiration (default 7 days)
 * - Bundled GRCh38 fallback for offline use
 * - Graceful degradation to existing file-based downloads
 *
 * Cache structure: ~/.decodingus/cache/genome-regions/{build}.json
 */
object GenomeRegionService {
  import GenomeRegionsCodecs.given

  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "genome-regions")
    try {
      Files.createDirectories(dir)
    } catch {
      case _: IOException =>
        Paths.get(System.getProperty("java.io.tmpdir"), "genome-regions-cache")
    }
    dir
  }

  // In-memory cache with timestamp
  @volatile private var memoryCache: Map[String, (GenomeRegions, Instant)] = Map.empty

  /**
   * Get genome regions for a reference build.
   *
   * Resolution order:
   * 1. In-memory cache (if not expired)
   * 2. Disk cache (if not expired)
   * 3. API fetch (if feature toggle enabled)
   * 4. Bundled resource fallback (for GRCh38 only)
   *
   * @param build Reference genome build (GRCh38, GRCh37, CHM13v2)
   * @param ec    Execution context
   * @return Future containing Either error message or GenomeRegions
   */
  def getRegions(build: String)(implicit ec: ExecutionContext): Future[Either[String, GenomeRegions]] = {
    val normalizedBuild = normalizeBuild(build)

    // Check in-memory cache first
    memoryCache.get(normalizedBuild) match {
      case Some((regions, cachedAt)) if !isExpired(cachedAt) =>
        return Future.successful(Right(regions))
      case _ => // Continue to disk cache or fetch
    }

    // Try disk cache
    loadFromDiskCache(normalizedBuild) match {
      case Some(regions) =>
        // Update memory cache and return
        memoryCache = memoryCache + (normalizedBuild -> (regions, Instant.now()))
        return Future.successful(Right(regions))
      case None => // Continue to API or fallback
    }

    // Feature toggle check for API access
    if (FeatureToggles.genomeRegionsApi.enabled) {
      fetchFromApiWithFallback(normalizedBuild)
    } else {
      // Feature disabled - try bundled resource only
      loadBundledResource(normalizedBuild) match {
        case Some(regions) =>
          memoryCache = memoryCache + (normalizedBuild -> (regions, Instant.now()))
          Future.successful(Right(regions))
        case None =>
          Future.successful(Left(s"Genome regions API disabled and no bundled resource for $normalizedBuild"))
      }
    }
  }

  /**
   * Fetches from API, falling back to bundled resource on error.
   */
  private def fetchFromApiWithFallback(build: String)(implicit ec: ExecutionContext): Future[Either[String, GenomeRegions]] = {
    DecodingUsClient.getGenomeRegions(build).map {
      case Right(regions) =>
        // Cache to disk and memory
        saveToDiskCache(build, regions)
        memoryCache = memoryCache + (build -> (regions, Instant.now()))
        Right(regions)

      case Left(apiError) if FeatureToggles.genomeRegionsApi.fallbackEnabled =>
        // API failed - try bundled resource
        println(s"[GenomeRegionService] API fetch failed ($apiError), trying bundled fallback")
        loadBundledResource(build) match {
          case Some(regions) =>
            memoryCache = memoryCache + (build -> (regions, Instant.now()))
            Right(regions)
          case None =>
            Left(s"API failed and no bundled resource available for $build: $apiError")
        }

      case Left(apiError) =>
        Left(apiError)
    }.recover { case e: Exception =>
      if (FeatureToggles.genomeRegionsApi.fallbackEnabled) {
        loadBundledResource(build) match {
          case Some(regions) =>
            memoryCache = memoryCache + (build -> (regions, Instant.now()))
            Right(regions)
          case None =>
            Left(s"Network error and no bundled resource: ${e.getMessage}")
        }
      } else {
        Left(s"Network error: ${e.getMessage}")
      }
    }
  }

  /**
   * Load from disk cache if not expired.
   */
  private def loadFromDiskCache(build: String): Option[GenomeRegions] = {
    val cachePath = cacheDir.resolve(s"$build.json")
    if (!Files.exists(cachePath)) return None

    // Check file modification time for expiration
    val modifiedTime = Files.getLastModifiedTime(cachePath).toInstant
    if (isExpired(modifiedTime)) {
      // Cache expired - delete and return None
      Files.deleteIfExists(cachePath)
      return None
    }

    // Parse cached JSON
    Using(Source.fromFile(cachePath.toFile, "UTF-8")) { source =>
      decode[GenomeRegions](source.mkString)
    }.toOption.flatMap {
      case Right(regions) => Some(regions)
      case Left(error) =>
        println(s"[GenomeRegionService] Failed to parse cached file: $error")
        Files.deleteIfExists(cachePath)
        None
    }
  }

  /**
   * Save regions to disk cache.
   */
  private def saveToDiskCache(build: String, regions: GenomeRegions): Unit = {
    try {
      val cachePath = cacheDir.resolve(s"$build.json")
      Files.createDirectories(cacheDir)
      val json = regions.asJson.noSpaces
      Files.writeString(cachePath, json, StandardCharsets.UTF_8,
        StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
    } catch {
      case e: IOException =>
        println(s"[GenomeRegionService] Failed to write cache: ${e.getMessage}")
    }
  }

  /**
   * Load bundled resource fallback.
   * Only GRCh38 is bundled for offline use.
   */
  private def loadBundledResource(build: String): Option[GenomeRegions] = {
    val resourcePath = s"/genome-regions/${build.toLowerCase}.json"
    Option(getClass.getResourceAsStream(resourcePath)).flatMap { stream =>
      Using(Source.fromInputStream(stream, "UTF-8")) { source =>
        decode[GenomeRegions](source.mkString)
      }.toOption.flatMap {
        case Right(regions) =>
          println(s"[GenomeRegionService] Loaded bundled resource for $build")
          Some(regions)
        case Left(error) =>
          println(s"[GenomeRegionService] Failed to parse bundled resource: $error")
          None
      }
    }
  }

  /**
   * Check if cached data has expired.
   */
  private def isExpired(cachedAt: Instant): Boolean = {
    val cacheDays = FeatureToggles.genomeRegionsApi.cacheDays
    val expiration = cachedAt.plus(Duration.ofDays(cacheDays))
    Instant.now().isAfter(expiration)
  }

  /**
   * Normalize reference build name to canonical form.
   */
  private def normalizeBuild(build: String): String = {
    build.toLowerCase match {
      case "grch38" | "hg38"                      => "GRCh38"
      case "grch37" | "hg19"                      => "GRCh37"
      case "chm13v2" | "chm13" | "t2t-chm13" | "hs1" => "CHM13v2"
      case _                                       => build
    }
  }

  /**
   * Clear all caches (memory and disk).
   */
  def clearCache(): Unit = {
    memoryCache = Map.empty
    if (Files.exists(cacheDir)) {
      Files.list(cacheDir).forEach(Files.deleteIfExists)
    }
  }

  /**
   * Clear cache for a specific build.
   */
  def clearCache(build: String): Unit = {
    val normalizedBuild = normalizeBuild(build)
    memoryCache = memoryCache - normalizedBuild
    val cachePath = cacheDir.resolve(s"$normalizedBuild.json")
    Files.deleteIfExists(cachePath)
  }

  /**
   * Check if regions are available for a build (cached or bundled).
   */
  def isAvailable(build: String): Boolean = {
    val normalizedBuild = normalizeBuild(build)
    memoryCache.contains(normalizedBuild) ||
      Files.exists(cacheDir.resolve(s"$normalizedBuild.json")) ||
      getClass.getResourceAsStream(s"/genome-regions/${normalizedBuild.toLowerCase}.json") != null
  }
}
