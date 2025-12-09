package com.decodingus.refgenome

import java.io.IOException
import java.nio.file.{Files, Path, Paths, StandardCopyOption}

/**
 * Cache for STR (Short Tandem Repeat) reference BED files.
 * Uses HipSTR reference files which contain known STR regions.
 */
class StrReferenceCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "str")
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"Failed to create STR cache directory: $dir")
        Files.createTempDirectory("str-cache")
    }
    dir
  }

  def getPath(referenceBuild: String): Option[Path] = {
    val bedFileName = s"$referenceBuild.hipstr_reference.bed"
    val bedPath = cacheDir.resolve(bedFileName)
    if (Files.exists(bedPath) && Files.size(bedPath) > 0) Some(bedPath) else None
  }

  def put(referenceBuild: String, file: Path): Path = {
    val bedFileName = s"$referenceBuild.hipstr_reference.bed"
    val targetPath = cacheDir.resolve(bedFileName)
    Files.move(file, targetPath, StandardCopyOption.REPLACE_EXISTING)
    targetPath
  }

  def getCacheDir: Path = cacheDir
}
