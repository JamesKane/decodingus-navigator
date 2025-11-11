package com.decodingus.haplogroup.tree

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupTree}
import sttp.client3.{HttpURLConnectionBackend, basicRequest}

abstract class TreeProvider(val treeType: TreeType) {
  def url: String

  def cachePrefix: String

  def progressMessage: String

  def parseTree(data: String, targetBuild: String): Either[String, HaplogroupTree]

  def buildTree(tree: HaplogroupTree): List[Haplogroup]

  def supportedBuilds: List[String]

  def sourceBuild: String // This is the native build of the tree data

  def loadTree(targetBuild: String): Either[String, List[Haplogroup]] = {
    val cache = new TreeCache()
    cache.get(cachePrefix) match {
      case Some(data) =>
        println(s"Found ${cachePrefix} in cache.")
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
  }
}
