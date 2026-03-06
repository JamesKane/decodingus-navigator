package com.decodingus.refgenome

import java.io.IOException
import java.nio.file.{Files, Path, Paths, StandardCopyOption}

/**
 * Base trait for file-based caches in ~/.decodingus/cache/{subdir}.
 *
 * Provides common directory initialization and file management operations.
 */
trait FileCache {
  protected def cacheSubdir: String

  protected lazy val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", cacheSubdir)
    try {
      Files.createDirectories(dir)
    } catch {
      case _: IOException =>
        println(s"Failed to create cache directory: $dir")
        Files.createTempDirectory(s"$cacheSubdir-cache")
    }
    dir
  }

  def getCacheDir: Path = cacheDir

  protected def resolve(fileName: String): Path = cacheDir.resolve(fileName)

  protected def cachedPath(fileName: String): Option[Path] = {
    val path = resolve(fileName)
    if (Files.exists(path) && Files.size(path) > 0) Some(path) else None
  }

  protected def moveToCache(source: Path, fileName: String): Path = {
    val target = resolve(fileName)
    Files.move(source, target, StandardCopyOption.REPLACE_EXISTING)
    target
  }

  protected def copyToCache(source: Path, fileName: String): Path = {
    val target = resolve(fileName)
    Files.copy(source, target, StandardCopyOption.REPLACE_EXISTING)
    target
  }
}
