package com.decodingus.ancestry.reference

import com.decodingus.ancestry.model.AncestryPanelType
import com.decodingus.config.FeatureToggles
import com.decodingus.refgenome.FileCache

import java.nio.file.{Files, Path}

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
class AncestryReferenceCache extends FileCache {
  protected def cacheSubdir = "ancestry"

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

  def getSitesVcfPath(panelType: AncestryPanelType, referenceBuild: String): Path =
    getPanelDir(panelType).resolve(s"${referenceBuild}_sites.vcf.gz")

  def getSitesVcfIndexPath(panelType: AncestryPanelType, referenceBuild: String): Path =
    getPanelDir(panelType).resolve(s"${referenceBuild}_sites.vcf.gz.tbi")

  def getAlleleFreqPath(panelType: AncestryPanelType): Path =
    getPanelDir(panelType).resolve("allele_freqs.bin")

  def getPcaLoadingsPath(panelType: AncestryPanelType): Path =
    getPanelDir(panelType).resolve("pca_loadings.bin")

  def getPopulationsPath(version: String = FeatureToggles.ancestryAnalysis.referenceVersion): Path =
    getVersionDir(version).resolve("populations.json")

  def isPanelAvailable(panelType: AncestryPanelType, referenceBuild: String): Boolean = {
    Files.exists(getSitesVcfPath(panelType, referenceBuild)) &&
      Files.exists(getSitesVcfIndexPath(panelType, referenceBuild)) &&
      Files.exists(getAlleleFreqPath(panelType)) &&
      Files.exists(getPcaLoadingsPath(panelType))
  }

  def listAvailablePanels(referenceBuild: String): Map[AncestryPanelType, Boolean] = {
    Map(
      AncestryPanelType.Aims -> isPanelAvailable(AncestryPanelType.Aims, referenceBuild),
      AncestryPanelType.GenomeWide -> isPanelAvailable(AncestryPanelType.GenomeWide, referenceBuild)
    )
  }

  def getPanelCacheSize(panelType: AncestryPanelType, referenceBuild: String): Long = {
    val files = List(
      getSitesVcfPath(panelType, referenceBuild),
      getSitesVcfIndexPath(panelType, referenceBuild),
      getAlleleFreqPath(panelType),
      getPcaLoadingsPath(panelType)
    )
    files.filter(Files.exists(_)).map(Files.size(_)).sum
  }

  def deletePanel(panelType: AncestryPanelType, referenceBuild: String): Unit = {
    val files = List(
      getSitesVcfPath(panelType, referenceBuild),
      getSitesVcfIndexPath(panelType, referenceBuild)
    )
    files.filter(Files.exists(_)).foreach(Files.delete(_))
  }
}
