package com.decodingus.haplogroup.tree

import java.io.{File, IOException}
import java.nio.file.{Files, Path, Paths, StandardCopyOption}
import scala.io.Source

class TreeCache {
  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "trees")
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"Failed to create cache directory: $dir")
        Files.createTempDirectory("tree-cache")
    }
    dir
  }

  /** Get the cache directory path */
  def getCacheDir: Path = cacheDir

  def get(prefix: String): Option[String] = {
    val path = cacheDir.resolve(s"$prefix.json")
    if (Files.exists(path)) {
      val source = Source.fromFile(path.toFile)
      try Some(source.mkString) finally source.close()
    } else {
      None
    }
  }

  def put(prefix: String, data: String): Unit = {
    val path = cacheDir.resolve(s"$prefix.json")
    Files.write(path, data.getBytes)
  }

  /**
   * Get the path for a sites VCF file associated with a tree.
   * Sites VCFs are cached alongside their source tree and can be reused.
   *
   * @param treePrefix     The tree cache prefix (e.g., "ftdna-ytree", "ftdna-mttree")
   * @param referenceBuild The reference build (e.g., "GRCh38")
   * @return Path to the sites VCF file
   */
  def getSitesVcfPath(treePrefix: String, referenceBuild: String): File = {
    cacheDir.resolve(s"$treePrefix-$referenceBuild-sites.vcf").toFile
  }

  /**
   * Check if a sites VCF exists and is newer than the tree JSON.
   */
  def isSitesVcfValid(treePrefix: String, referenceBuild: String): Boolean = {
    val sitesVcf = getSitesVcfPath(treePrefix, referenceBuild)
    val treeJson = cacheDir.resolve(s"$treePrefix.json").toFile
    sitesVcf.exists() && treeJson.exists() && sitesVcf.lastModified() >= treeJson.lastModified()
  }

  /**
   * Get the path for a lifted sites VCF file (tree sites lifted to a different reference build).
   * These are cached to avoid re-lifting for every sample using the same tree/target build combination.
   *
   * @param treePrefix  The tree cache prefix (e.g., "ftdna-ytree", "ftdna-mttree")
   * @param sourceBuild The source reference build (tree's native build, e.g., "GRCh38")
   * @param targetBuild The target reference build (BAM's build, e.g., "CHM13v2")
   * @return Path to the lifted sites VCF file
   */
  def getLiftedSitesVcfPath(treePrefix: String, sourceBuild: String, targetBuild: String): File = {
    cacheDir.resolve(s"$treePrefix-$sourceBuild-to-$targetBuild-sites.vcf").toFile
  }

  /**
   * Check if a lifted sites VCF exists and is newer than the source sites VCF.
   */
  def isLiftedSitesVcfValid(treePrefix: String, sourceBuild: String, targetBuild: String): Boolean = {
    val liftedVcf = getLiftedSitesVcfPath(treePrefix, sourceBuild, targetBuild)
    val sourceVcf = getSitesVcfPath(treePrefix, sourceBuild)
    liftedVcf.exists() && sourceVcf.exists() && liftedVcf.lastModified() >= sourceVcf.lastModified()
  }
}
