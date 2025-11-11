package com.decodingus.haplogroup.tree

import java.io.IOException
import java.nio.file.{Files, Path, Paths, StandardCopyOption}
import scala.io.Source

class TreeCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".cache", "decodingus-tools", "trees")
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"Failed to create cache directory: $dir")
        Files.createTempDirectory("tree-cache")
    }
    dir
  }

  def get(prefix: String): Option[String] = {
    val path = cacheDir.resolve(s"$prefix.json")
    if (Files.exists(path)) {
      val source = Source.fromFile(path.toFile)
      try Some(source.mkString) finally source.close()
    } else {
      None
    }
  }

  def put(prefix: String, data: String): Unit = {
    val path = cacheDir.resolve(s"$prefix.json")
    Files.write(path, data.getBytes)
  }
}
