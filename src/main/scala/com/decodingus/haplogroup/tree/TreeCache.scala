package com.decodingus.haplogroup.tree

import com.decodingus.refgenome.FileCache

import java.io.File
import java.nio.file.{Files, Path}
import scala.io.Source

class TreeCache extends FileCache {
  protected def cacheSubdir = "trees"

  def get(prefix: String): Option[String] = {
    val path = resolve(s"$prefix.json")
    if (Files.exists(path)) {
      val source = Source.fromFile(path.toFile)
      try Some(source.mkString) finally source.close()
    } else {
      None
    }
  }

  def put(prefix: String, data: String): Unit = {
    Files.write(resolve(s"$prefix.json"), data.getBytes)
  }

  /**
   * Get the path for a sites VCF file associated with a tree.
   * Sites VCFs are cached alongside their source tree and can be reused.
   */
  def getSitesVcfPath(treePrefix: String, referenceBuild: String): File = {
    resolve(s"$treePrefix-$referenceBuild-sites.vcf").toFile
  }

  /** Check if a sites VCF exists and is newer than the tree JSON. */
  def isSitesVcfValid(treePrefix: String, referenceBuild: String): Boolean = {
    val sitesVcf = getSitesVcfPath(treePrefix, referenceBuild)
    val treeJson = resolve(s"$treePrefix.json").toFile
    sitesVcf.exists() && treeJson.exists() && sitesVcf.lastModified() >= treeJson.lastModified()
  }

  /**
   * Get the path for a lifted sites VCF file (tree sites lifted to a different reference build).
   */
  def getLiftedSitesVcfPath(treePrefix: String, sourceBuild: String, targetBuild: String): File = {
    resolve(s"$treePrefix-$sourceBuild-to-$targetBuild-sites.vcf").toFile
  }

  /** Check if a lifted sites VCF exists and is newer than the source sites VCF. */
  def isLiftedSitesVcfValid(treePrefix: String, sourceBuild: String, targetBuild: String): Boolean = {
    val liftedVcf = getLiftedSitesVcfPath(treePrefix, sourceBuild, targetBuild)
    val sourceVcf = getSitesVcfPath(treePrefix, sourceBuild)
    liftedVcf.exists() && sourceVcf.exists() && liftedVcf.lastModified() >= sourceVcf.lastModified()
  }
}
