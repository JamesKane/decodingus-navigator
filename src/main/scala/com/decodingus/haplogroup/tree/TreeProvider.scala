package com.decodingus.haplogroup.tree

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupTree}
import sttp.client3.{HttpURLConnectionBackend, basicRequest}

import scala.collection.mutable

/**
 * Global in-memory cache for parsed haplogroup trees.
 * The FTDNA Y-DNA tree is 113MB JSON and takes significant time to parse.
 * Caching the parsed tree in memory avoids repeated parsing.
 */
object ParsedTreeCache {
  // Key: (cachePrefix, targetBuild) -> parsed tree
  private val cache: mutable.Map[(String, String), List[Haplogroup]] = mutable.Map.empty

  def get(cachePrefix: String, targetBuild: String): Option[List[Haplogroup]] = {
    cache.get((cachePrefix, targetBuild))
  }

  def put(cachePrefix: String, targetBuild: String, tree: List[Haplogroup]): Unit = {
    cache.put((cachePrefix, targetBuild), tree)
  }

  def clear(): Unit = {
    cache.clear()
  }
}

abstract class TreeProvider(val treeType: TreeType) {
  def url: String

  def cachePrefix: String

  def progressMessage: String

  def parseTree(data: String, targetBuild: String): Either[String, HaplogroupTree]

  def buildTree(tree: HaplogroupTree): List[Haplogroup]

  def supportedBuilds: List[String]

  def sourceBuild: String // This is the native build of the tree data

  def loadTree(targetBuild: String): Either[String, List[Haplogroup]] = {
    // Check in-memory cache first (avoids parsing 113MB JSON repeatedly)
    ParsedTreeCache.get(cachePrefix, targetBuild) match {
      case Some(tree) =>
        println(s"Using in-memory cached tree for ${cachePrefix}.")
        Right(tree)
      case None =>
        // Not in memory, check disk cache or download
        val cache = new TreeCache()
        val result = cache.get(cachePrefix) match {
          case Some(data) =>
            println(s"Found ${cachePrefix} in disk cache. Parsing...")
            parseTree(data, targetBuild).map(buildTree)
          case None =>
            println(s"Downloading ${cachePrefix}...")
            val backend = HttpURLConnectionBackend()
            val response = basicRequest.get(sttp.model.Uri.unsafeParse(url)).send(backend)
            response.body.flatMap { data =>
              println("Download complete. Caching tree.")
              cache.put(cachePrefix, data)
              parseTree(data, targetBuild).map(buildTree)
            }
        }

        // Store in memory cache for subsequent calls
        result.foreach { tree =>
          println(s"Caching parsed tree in memory for ${cachePrefix}.")
          ParsedTreeCache.put(cachePrefix, targetBuild, tree)
        }

        result
    }
  }
}
