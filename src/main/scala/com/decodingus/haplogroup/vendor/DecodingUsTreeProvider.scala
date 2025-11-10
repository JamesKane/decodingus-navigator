package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.model._
import com.decodingus.haplogroup.tree.{TreeProvider, TreeType}
import io.circe.generic.auto._
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

  override def parseTree(data: String): Either[io.circe.Error, HaplogroupTree] = {
    decode[List[ApiNode]](data).map { apiNodes =>
      val nameToId = apiNodes.zipWithIndex.map { case (node, i) => node.name -> i.toLong }.toMap
      val rootId = apiNodes.zipWithIndex.find(_._1.parentName.isEmpty).map(_._2.toLong).getOrElse(0L)

      val allNodes = apiNodes.zipWithIndex.map { case (node, i) =>
        val haplogroupId = i.toLong
        val parentId = node.parentName.flatMap(nameToId.get).getOrElse(if (haplogroupId == rootId) 0L else rootId)
        val loci = node.variants.map { v =>
          val lociType = if (v.variantType == "SNP") LociType.SNP else LociType.INDEL
          val coordinates = v.coordinates.map { case (build, coord) =>
            val buildId = build match {
              case "CM000686.2" | "NC_000024.10" => "GRCh38"
              case "NC_060948.1" | "CP086569.2" => "T2T-CHM13v2.0"
              case "CM000686.1" => "GRCh37"
              case _ => build
            }
            val chromosome = if (buildId == "GRCh37") "Y" else "chrY"
            buildId -> LociCoordinate(coord.start, chromosome, coord.anc, coord.der)
          }
          Locus(v.name, lociType, coordinates)
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

  override def buildTree(tree: HaplogroupTree, nodeId: Long, treeType: TreeType): Option[Haplogroup] = {
    val nodeStr = nodeId.toString
    tree.allNodes.get(nodeStr).map { node =>
      val children = node.children.flatMap(childId => buildTree(tree, childId, treeType))
      val parentName = if (node.parent_id == 0) None else tree.allNodes.get(node.parent_id.toString).map(_.name)
      Haplogroup(node.name, parentName, node.loci, children)
    }
  }

  override def supportedBuilds: List[String] = List("GRCh38", "GRCh37", "T2T-CHM13v2.0")
}