package com.decodingus.refgenome

import com.decodingus.config.ReferenceConfigService

import java.io.IOException
import java.nio.file.{Files, Path, Paths}

class ReferenceCache {
  private def cacheDir: Path = {
    val dir = ReferenceConfigService.getCacheDir
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"Failed to create cache directory: $dir")
        // Fallback to a temporary directory
        Files.createTempDirectory("ref-cache")
    }
    dir
  }

  /**
   * Gets the path for a reference build.
   * Checks in order:
   * 1. User-specified local path from config
   * 2. Default cache directory
   */
  def getPath(referenceBuild: String): Option[Path] = {
    // Use the config service which checks user paths first, then cache
    ReferenceConfigService.getReferencePath(referenceBuild)
  }

  /**
   * Stores a reference file in the cache directory.
   */
  def put(referenceBuild: String, file: Path): Path = {
    val targetPath = cacheDir.resolve(s"$referenceBuild.fa.gz")
    Files.move(file, targetPath)
    targetPath
  }

  /**
   * Gets the cache directory path (for UI display).
   */
  def getCacheDirectory: Path = cacheDir
}
