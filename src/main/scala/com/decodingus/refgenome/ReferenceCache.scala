package com.decodingus.refgenome

import java.io.IOException
import java.nio.file.{Files, Path, Paths}

class ReferenceCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".cache", "decodingus-tools", "references")
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

  def getPath(referenceBuild: String): Option[Path] = {
    val refPath = cacheDir.resolve(s"$referenceBuild.fa.gz")
    if (Files.exists(refPath)) Some(refPath) else None
  }

  def put(referenceBuild: String, file: Path): Path = {
    val targetPath = cacheDir.resolve(s"$referenceBuild.fa.gz")
    Files.move(file, targetPath)
    targetPath
  }
}
