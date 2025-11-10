package com.decodingus.haplogroup.tree

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupTree}

enum TreeType {
  case YDNA, MTDNA
}

trait TreeProvider {
  def url(treeType: TreeType): String
  def cachePrefix(treeType: TreeType): String
  def progressMessage(treeType: TreeType): String
  def parseTree(data: String): Either[io.circe.Error, HaplogroupTree]
  def buildTree(tree: HaplogroupTree, nodeId: Long, treeType: TreeType): Option[Haplogroup]
  def supportedBuilds: List[String]
}
