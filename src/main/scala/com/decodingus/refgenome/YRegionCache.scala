package com.decodingus.refgenome

import java.io.IOException
import java.nio.file.{Files, Path, Paths, StandardCopyOption}

/**
 * Cache for Y chromosome region annotation files.
 *
 * Manages cached GFF3 and BED files for:
 * - Cytobands (from ybrowse.org)
 * - Palindromes (from ybrowse.org)
 * - STR regions (from ybrowse.org)
 * - PAR regions (from GIAB stratifications)
 * - XTR regions (from GIAB stratifications)
 * - Ampliconic regions (from GIAB stratifications)
 * - Centromeres (from UCSC)
 * - Heterochromatin (hardcoded Yq12 boundaries)
 *
 * Files are stored at: ~/.decodingus/cache/yregions/{build}_{type}.{ext}
 */
class YRegionCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "yregions")
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"[YRegionCache] Failed to create cache directory: $dir")
        Files.createTempDirectory("yregions-cache")
    }
    dir
  }

  /**
   * Get cached file path for a region type and reference build.
   *
   * @param regionType     Type of region (e.g., "cytobands", "palindromes", "par")
   * @param referenceBuild Reference genome build (e.g., "GRCh38", "GRCh37")
   * @return Some(path) if cached file exists, None otherwise
   */
  def getPath(regionType: String, referenceBuild: String): Option[Path] = {
    val fileName = buildFileName(regionType, referenceBuild)
    val path = cacheDir.resolve(fileName)
    if (Files.exists(path) && Files.size(path) > 0) Some(path) else None
  }

  /**
   * Store a region file in the cache.
   *
   * @param regionType     Type of region
   * @param referenceBuild Reference genome build
   * @param sourceFile     Path to the file to cache (will be moved)
   * @return Path to the cached file
   */
  def put(regionType: String, referenceBuild: String, sourceFile: Path): Path = {
    val fileName = buildFileName(regionType, referenceBuild)
    val targetPath = cacheDir.resolve(fileName)
    Files.move(sourceFile, targetPath, StandardCopyOption.REPLACE_EXISTING)
    targetPath
  }

  /**
   * Copy a file to the cache (preserving original).
   *
   * @param regionType     Type of region
   * @param referenceBuild Reference genome build
   * @param sourceFile     Path to the file to cache (will be copied)
   * @return Path to the cached file
   */
  def putCopy(regionType: String, referenceBuild: String, sourceFile: Path): Path = {
    val fileName = buildFileName(regionType, referenceBuild)
    val targetPath = cacheDir.resolve(fileName)
    Files.copy(sourceFile, targetPath, StandardCopyOption.REPLACE_EXISTING)
    targetPath
  }

  /**
   * Check if all required region files are cached for a build.
   *
   * @param referenceBuild Reference genome build
   * @return True if all region files are cached
   */
  def isComplete(referenceBuild: String): Boolean = {
    YRegionCache.requiredRegionTypes.forall(getPath(_, referenceBuild).isDefined)
  }

  /**
   * Get list of missing region types for a build.
   *
   * @param referenceBuild Reference genome build
   * @return List of region types that need to be downloaded
   */
  def getMissing(referenceBuild: String): List[String] = {
    YRegionCache.requiredRegionTypes.filterNot(rt => getPath(rt, referenceBuild).isDefined)
  }

  /**
   * Clear cache for a specific build.
   */
  def clearBuild(referenceBuild: String): Unit = {
    YRegionCache.allRegionTypes.foreach { regionType =>
      val fileName = buildFileName(regionType, referenceBuild)
      val path = cacheDir.resolve(fileName)
      Files.deleteIfExists(path)
    }
  }

  /**
   * Clear entire cache.
   */
  def clearAll(): Unit = {
    if (Files.exists(cacheDir)) {
      Files.list(cacheDir).forEach(Files.deleteIfExists)
    }
  }

  def getCacheDir: Path = cacheDir

  private def buildFileName(regionType: String, referenceBuild: String): String = {
    val ext = if (YRegionCache.gff3Types.contains(regionType)) "gff3" else "bed"
    s"${referenceBuild}_$regionType.$ext"
  }
}

object YRegionCache {
  // GFF3 format region types (from ybrowse.org)
  val gff3Types: Set[String] = Set("cytobands", "palindromes", "strs")

  // BED format region types (from GIAB, UCSC)
  val bedTypes: Set[String] = Set("par", "xtr", "ampliconic", "centromeres", "xdegenerate")

  // Required region types for full annotation
  val requiredRegionTypes: List[String] = List(
    "cytobands",
    "palindromes",
    "strs",
    "par",
    "xtr",
    "ampliconic"
  )

  // Optional region types (may not be available for all builds)
  val optionalRegionTypes: List[String] = List(
    "centromeres",
    "xdegenerate"
  )

  // All region types
  val allRegionTypes: List[String] = requiredRegionTypes ++ optionalRegionTypes
}
