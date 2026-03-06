package com.decodingus.refgenome

import java.nio.file.{Files, Path}

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
class YRegionCache extends FileCache {
  protected def cacheSubdir = "yregions"

  def getPath(regionType: String, referenceBuild: String): Option[Path] =
    cachedPath(buildFileName(regionType, referenceBuild))

  def put(regionType: String, referenceBuild: String, sourceFile: Path): Path =
    moveToCache(sourceFile, buildFileName(regionType, referenceBuild))

  def putCopy(regionType: String, referenceBuild: String, sourceFile: Path): Path =
    copyToCache(sourceFile, buildFileName(regionType, referenceBuild))

  def isComplete(referenceBuild: String): Boolean =
    YRegionCache.requiredRegionTypes.forall(getPath(_, referenceBuild).isDefined)

  def getMissing(referenceBuild: String): List[String] =
    YRegionCache.requiredRegionTypes.filterNot(rt => getPath(rt, referenceBuild).isDefined)

  def clearBuild(referenceBuild: String): Unit =
    YRegionCache.allRegionTypes.foreach { regionType =>
      Files.deleteIfExists(resolve(buildFileName(regionType, referenceBuild)))
    }

  def clearAll(): Unit =
    if (Files.exists(cacheDir)) {
      Files.list(cacheDir).forEach(Files.deleteIfExists)
    }

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
