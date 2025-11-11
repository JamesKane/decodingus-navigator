package com.decodingus.haplogroup.tree

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupTree}
import sttp.client3.{HttpURLConnectionBackend, basicRequest}

trait TreeProvider {
  def url(treeType: TreeType): String

  def cachePrefix(treeType: TreeType): String

  def progressMessage(treeType: TreeType): String

  def parseTree(data: String, targetBuild: String): Either[String, HaplogroupTree]

  def buildTree(tree: HaplogroupTree): List[Haplogroup]

  def supportedBuilds: List[String]

  def sourceBuild: String // This is the native build of the tree data

  def loadTree(treeType: TreeType, targetBuild: String): Either[String, List[Haplogroup]] = {
    val cache = new TreeCache()
    cache.get(cachePrefix(treeType)) match {
      case Some(data) =>
        println(s"Found ${cachePrefix(treeType)} in cache.")
        parseTree(data, targetBuild).map(buildTree)
      case None =>
        println(s"Downloading ${cachePrefix(treeType)}...")
        val backend = HttpURLConnectionBackend()
        val response = basicRequest.get(sttp.model.Uri.unsafeParse(url(treeType))).send(backend)
        response.body.flatMap { data =>
          println("Download complete. Caching tree.")
          cache.put(cachePrefix(treeType), data)
          parseTree(data, targetBuild).map(buildTree)
        }
    }
  }
}
