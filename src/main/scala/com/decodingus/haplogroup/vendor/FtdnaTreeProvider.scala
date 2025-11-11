package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.model.*
import com.decodingus.haplogroup.tree.{TreeProvider, TreeType}
import io.circe.generic.auto.*
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

  override def parseTree(data: String, targetBuild: String): Either[String, HaplogroupTree] = {
    decode[FtdnaTreeJson](data).left.map(_.toString).map { ftdnaTree =>
      val allNodes = ftdnaTree.allNodes.map { case (id, node) =>
        val loci = node.variants.flatMap { v =>
          v.position.map { pos =>
            Locus(v.variant, pos.toLong, v.ancestral, v.derived)
          }
        }
        id -> HaplogroupNode(node.haplogroupId, node.parentId, node.name, node.isRoot, loci, node.children)
      }
      HaplogroupTree(allNodes)
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

  override def supportedBuilds: List[String] = List("GRCh38", "rCRS")

  override def sourceBuild: String = "GRCh38"
}
