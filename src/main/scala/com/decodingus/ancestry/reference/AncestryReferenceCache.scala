package com.decodingus.ancestry.reference

import com.decodingus.ancestry.model.AncestryPanelType
import com.decodingus.config.FeatureToggles

import java.io.IOException
import java.nio.file.{Files, Path, Paths}

/**
 * Manages cached ancestry reference data.
 *
 * Cache structure:
 * ~/.decodingus/cache/ancestry/{version}/
 * ├── populations.json
 * ├── aims/
 * │   ├── {build}_sites.vcf.gz
 * │   ├── {build}_sites.vcf.gz.tbi
 * │   ├── allele_freqs.bin
 * │   └── pca_loadings.bin
 * └── genome-wide/
 * └── ...
 */
class AncestryReferenceCache {

  private val cacheDir: Path = {
    val dir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "ancestry")
    try {
      Files.createDirectories(dir)
    } catch {
      case e: IOException =>
        println(s"Failed to create cache directory: $dir")
        Files.createTempDirectory("ancestry-cache")
    }
    dir
  }

  /** Get the base cache directory */
  def getCacheDir: Path = cacheDir

  /** Get version-specific directory */
  def getVersionDir(version: String = FeatureToggles.ancestryAnalysis.referenceVersion): Path = {
    val dir = cacheDir.resolve(version)
    Files.createDirectories(dir)
    dir
  }

  /** Get panel-specific directory (aims or genome-wide) */
  def getPanelDir(panelType: AncestryPanelType, version: String = FeatureToggles.ancestryAnalysis.referenceVersion): Path = {
    val panelName = panelType match {
      case AncestryPanelType.Aims => "aims"
      case AncestryPanelType.GenomeWide => "genome-wide"
    }
    val dir = getVersionDir(version).resolve(panelName)
    Files.createDirectories(dir)
    dir
  }

  /** Get path to sites VCF for a panel and reference build */
  def getSitesVcfPath(panelType: AncestryPanelType, referenceBuild: String): Path = {
    getPanelDir(panelType).resolve(s"${referenceBuild}_sites.vcf.gz")
  }

  /** Get path to sites VCF index */
  def getSitesVcfIndexPath(panelType: AncestryPanelType, referenceBuild: String): Path = {
    getPanelDir(panelType).resolve(s"${referenceBuild}_sites.vcf.gz.tbi")
  }

  /** Get path to allele frequency matrix */
  def getAlleleFreqPath(panelType: AncestryPanelType): Path = {
    getPanelDir(panelType).resolve("allele_freqs.bin")
  }

  /** Get path to PCA loadings */
  def getPcaLoadingsPath(panelType: AncestryPanelType): Path = {
    getPanelDir(panelType).resolve("pca_loadings.bin")
  }

  /** Get path to populations definition */
  def getPopulationsPath(version: String = FeatureToggles.ancestryAnalysis.referenceVersion): Path = {
    getVersionDir(version).resolve("populations.json")
  }

  /**
   * Check if a panel is fully available for a reference build.
   */
  def isPanelAvailable(panelType: AncestryPanelType, referenceBuild: String): Boolean = {
    val sitesVcf = getSitesVcfPath(panelType, referenceBuild)
    val sitesIndex = getSitesVcfIndexPath(panelType, referenceBuild)
    val alleleFreq = getAlleleFreqPath(panelType)
    val pcaLoadings = getPcaLoadingsPath(panelType)

    Files.exists(sitesVcf) &&
      Files.exists(sitesIndex) &&
      Files.exists(alleleFreq) &&
      Files.exists(pcaLoadings)
  }

  /**
   * List available panels with their status.
   */
  def listAvailablePanels(referenceBuild: String): Map[AncestryPanelType, Boolean] = {
    Map(
      AncestryPanelType.Aims -> isPanelAvailable(AncestryPanelType.Aims, referenceBuild),
      AncestryPanelType.GenomeWide -> isPanelAvailable(AncestryPanelType.GenomeWide, referenceBuild)
    )
  }

  /**
   * Get cache size for a panel (in bytes).
   */
  def getPanelCacheSize(panelType: AncestryPanelType, referenceBuild: String): Long = {
    val files = List(
      getSitesVcfPath(panelType, referenceBuild),
      getSitesVcfIndexPath(panelType, referenceBuild),
      getAlleleFreqPath(panelType),
      getPcaLoadingsPath(panelType)
    )
    files.filter(Files.exists(_)).map(Files.size(_)).sum
  }

  /**
   * Delete cached panel data.
   */
  def deletePanel(panelType: AncestryPanelType, referenceBuild: String): Unit = {
    val files = List(
      getSitesVcfPath(panelType, referenceBuild),
      getSitesVcfIndexPath(panelType, referenceBuild)
    )
    files.filter(Files.exists(_)).foreach(Files.delete(_))
  }
}
