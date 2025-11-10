package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.model._
import com.decodingus.haplogroup.tree.{TreeProvider, TreeType}
import io.circe.generic.auto._
import io.circe.parser.decode

import scala.collection.mutable

case class FtdnaVariant(
  variant: String,
  position: Option[Int],
  ancestral: String,
  derived: String,
  region: String,
  id: Option[Long]
)

case class FtdnaNode(
  haplogroupId: Long,
  parentId: Long,
  name: String,
  isRoot: Boolean,
  root: String,
  kitsCount: Int,
  subBranches: Int,
  bigYCount: Int,
  variants: List[FtdnaVariant],
  children: List[Long]
)

case class FtdnaTreeJson(
  allNodes: Map[String, FtdnaNode]
)

class FtdnaTreeProvider extends TreeProvider {
  override def url(treeType: TreeType): String = treeType match {
    case TreeType.YDNA => "https://www.familytreedna.com/public/y-dna-haplotree/get"
    case TreeType.MTDNA => "https://www.familytreedna.com/public/mt-dna-haplotree/get"
  }

  override def cachePrefix(treeType: TreeType): String = treeType match {
    case TreeType.YDNA => "ftdna-ytree"
    case TreeType.MTDNA => "ftdna-mttree"
  }

  override def progressMessage(treeType: TreeType): String = treeType match {
    case TreeType.YDNA => "Downloading FTDNA Y-DNA tree..."
    case TreeType.MTDNA => "Downloading FTDNA MT-DNA tree..."
  }

  override def parseTree(data: String): Either[io.circe.Error, HaplogroupTree] = {
    decode[FtdnaTreeJson](data).map { ftdnaTree =>
      val allNodes = ftdnaTree.allNodes.map { case (id, node) =>
        val loci = node.variants.map { v =>
          val coordinates = v.position.map { pos =>
            "GRCh38" -> LociCoordinate(pos.toLong, "chrY", v.ancestral, v.derived)
          }.toMap
          Locus(v.variant, LociType.SNP, coordinates)
        }
        id -> HaplogroupNode(node.haplogroupId, node.parentId, node.name, node.isRoot, loci, node.children)
      }
      HaplogroupTree(allNodes)
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

  override def supportedBuilds: List[String] = List("GRCh38", "rCRS")
}