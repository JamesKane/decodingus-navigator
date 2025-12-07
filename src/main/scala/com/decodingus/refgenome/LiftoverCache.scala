package com.decodingus.refgenome

import java.io.IOException
import java.nio.file.{Files, Path, Paths, StandardCopyOption}

class LiftoverCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "liftover")
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"Failed to create cache directory: $dir")
        Files.createTempDirectory("liftover-cache")
    }
    dir
  }

  def getPath(from: String, to: String): Option[Path] = {
    val chainFileName = s"${from}To${to.capitalize}.over.chain.gz"
    val chainPath = cacheDir.resolve(chainFileName)
    if (Files.exists(chainPath)) Some(chainPath) else None
  }

  def put(from: String, to: String, file: Path): Path = {
    val chainFileName = s"${from}To${to.capitalize}.over.chain.gz"
    val targetPath = cacheDir.resolve(chainFileName)
    Files.move(file, targetPath, StandardCopyOption.REPLACE_EXISTING)
    targetPath
  }
}
