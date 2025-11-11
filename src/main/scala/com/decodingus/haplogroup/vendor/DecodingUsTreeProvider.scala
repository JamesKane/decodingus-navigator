package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.model.*
import com.decodingus.haplogroup.tree.{TreeProvider, TreeType}
import io.circe.generic.auto.*
import io.circe.parser.decode

import scala.collection.mutable

case class ApiCoordinate(
                          start: Long,
                          stop: Long,
                          anc: String,
                          der: String
                        )

case class ApiVariant(
                       name: String,
                       coordinates: Map[String, ApiCoordinate],
                       variantType: String
                     )

case class ApiNode(
                    name: String,
                    parentName: Option[String],
                    variants: List[ApiVariant],
                    lastUpdated: String,
                    isBackbone: Boolean
                  )

class DecodingUsTreeProvider extends TreeProvider {
  override def url(treeType: TreeType): String = treeType match {
    case TreeType.YDNA => "https://decoding-us.com/api/v1/y-tree"
    case TreeType.MTDNA => throw new UnsupportedOperationException("MT-DNA tree not yet supported by DecodingUs")
  }

  override def cachePrefix(treeType: TreeType): String = treeType match {
    case TreeType.YDNA => "decodingus-ytree"
    case TreeType.MTDNA => throw new UnsupportedOperationException("MT-DNA tree not yet supported by DecodingUs")
  }

  override def progressMessage(treeType: TreeType): String = treeType match {
    case TreeType.YDNA => "Downloading DecodingUs Y-DNA tree..."
    case TreeType.MTDNA => throw new UnsupportedOperationException("MT-DNA tree not yet supported by DecodingUs")
  }

  override def parseTree(data: String, targetBuild: String): Either[String, HaplogroupTree] = {
    val buildMap = Map(
      "CM000686.2" -> "GRCh38",
      "CM000686.1" -> "GRCh37",
      "CP086569.2" -> "CHM13v2"
    )

    decode[List[ApiNode]](data).left.map(_.toString).map { apiNodes =>
      val nameToId = apiNodes.zipWithIndex.map { case (node, i) => node.name -> i.toLong }.toMap
      val rootId = apiNodes.zipWithIndex.find(_._1.parentName.isEmpty).map(_._2.toLong).getOrElse(0L)

      val allNodes = apiNodes.zipWithIndex.map { case (node, i) =>
        val haplogroupId = i.toLong
        val parentId = node.parentName.flatMap(nameToId.get).getOrElse(if (haplogroupId == rootId) 0L else rootId)

        val loci = node.variants.flatMap { v =>
          v.coordinates.headOption.flatMap { case (apiBuild, coord) =>
            buildMap.get(apiBuild).flatMap { internalBuild =>
              if (internalBuild == targetBuild) {
                Some(Locus(v.name, coord.start, coord.anc, coord.der))
              } else {
                None
              }
            }
          }
        }

        haplogroupId.toString -> HaplogroupNode(haplogroupId, parentId, node.name, haplogroupId == rootId, loci, List())
      }.toMap

      val childrenMap = mutable.Map[Long, List[Long]]()
      allNodes.values.foreach { node =>
        if (node.parent_id != 0) {
          childrenMap(node.parent_id) = node.haplogroup_id :: childrenMap.getOrElse(node.parent_id, List())
        }
      }

      val finalNodes = allNodes.map { case (id, node) =>
        id -> node.copy(children = childrenMap.getOrElse(node.haplogroup_id, List()))
      }

      HaplogroupTree(finalNodes)
    }
  }

  override def buildTree(tree: HaplogroupTree): List[Haplogroup] = {
    val rootNodes = tree.allNodes.values.filter(_.is_root).toList
    rootNodes.map(root => buildSubTree(root.haplogroup_id, tree, None))
  }

  private def buildSubTree(nodeId: Long, tree: HaplogroupTree, parentName: Option[String]): Haplogroup = {
    val node = tree.allNodes(nodeId.toString)
    val children = node.children.map(childId => buildSubTree(childId, tree, Some(node.name)))
    Haplogroup(node.name, parentName, node.loci, children)
  }

  override def supportedBuilds: List[String] = List("GRCh38", "GRCh37", "CHM13v2")

  override def sourceBuild: String = "GRCh38"
}
