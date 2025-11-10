package com.decodingus.haplogroup.tree

import com.decodingus.haplogroup.model.HaplogroupTree
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import sttp.client3._
import sttp.client3.circe._

import java.nio.file.{Files, Path, Paths}
import java.time.Instant
import java.time.temporal.ChronoUnit

object TreeProviderType extends Enumeration {
  val FTDNA, DecodingUs = Value
}

class TreeCache(treeType: TreeType, providerType: TreeProviderType.Value) {
  val provider: TreeProvider = providerType match {
    case TreeProviderType.FTDNA => new FtdnaTreeProvider
    case TreeProviderType.DecodingUs => new DecodingUsTreeProvider
  }

  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".cache", "decodingus-tools", provider.cachePrefix(treeType))
    Files.createDirectories(dir)
    dir
  }

  private def getCachePath: Path = {
    val now = Instant.now()
    cacheDir.resolve(s"${provider.cachePrefix(treeType)}_${now.getEpochSecond}.json")
  }

  private def isCacheValid(path: Path): Boolean = {
    if (Files.exists(path)) {
      val modifiedTime = Files.getLastModifiedTime(path).toInstant
      val sevenDaysAgo = Instant.now().minus(7, ChronoUnit.DAYS)
      modifiedTime.isAfter(sevenDaysAgo)
    } else {
      false
    }
  }

  def getTree: Either[String, HaplogroupTree] = {
    val cachePath = getCachePath
    if (isCacheValid(cachePath)) {
      println(s"Using cached tree from $cachePath")
      val contents = Files.readString(cachePath)
      provider.parseTree(contents).left.map(_.getMessage)
    } else {
      println(s"Downloading tree from ${provider.url(treeType)}")
      val request = basicRequest.get(uri"${provider.url(treeType)}").response(asString)
      val backend = HttpURLConnectionBackend()
      val response = request.send(backend)

      response.body match {
        case Right(treeJson) =>
          println(s"Caching tree to $cachePath")
          Files.writeString(cachePath, treeJson)
          provider.parseTree(treeJson).left.map(_.getMessage)
        case Left(error) =>
          Left(s"Failed to download tree: $error")
      }
    }
  }
}