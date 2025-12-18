package com.decodingus.config

import com.typesafe.config.{Config, ConfigFactory}

import scala.jdk.CollectionConverters.*
import scala.util.Try

/**
 * Configuration for labs, sequencing centers, and genotyping vendors.
 * Provides configurable display names and abbreviations (up to 6 chars) for the Data Sources tab.
 *
 * Labs are loaded from labs.conf and can be looked up by:
 * - ID (the HOCON key, e.g., "familytreedna")
 * - Display name (e.g., "FamilyTreeDNA")
 * - Alias (e.g., "FTDNA", "Family Tree DNA")
 *
 * Labs can be filtered by capability for use in context-specific dialogs.
 */
object LabsConfig {

  /** Known capability types for labs */
  object Capability {
    val Wgs = "wgs"
    val YDna = "y-dna"
    val MtDna = "mt-dna"
    val Str = "str"
    val Chip = "chip"
    val VcfDownload = "vcf-download"
    val BamDownload = "bam-download"
  }

  /**
   * Represents a lab, sequencing center, or vendor configuration.
   *
   * @param id           Unique identifier (the HOCON key)
   * @param displayName  Full display name
   * @param abbreviation Short code (up to 6 chars) for UI indicators
   * @param category     Category: commercial-lab, consumer-vendor, sequencing-platform, academic
   * @param capabilities List of capabilities: wgs, y-dna, mt-dna, str, chip, vcf-download
   * @param website      Optional website URL
   * @param aliases      Alternative names for matching
   */
  case class Lab(
    id: String,
    displayName: String,
    abbreviation: String,
    category: String,
    capabilities: Set[String],
    website: Option[String],
    aliases: List[String]
  ) {
    /** Returns abbreviation truncated to maxLength chars (default 6) */
    def abbreviationTruncated(maxLength: Int = 6): String = abbreviation.take(maxLength)

    /** Check if this lab has a specific capability */
    def hasCapability(capability: String): Boolean = capabilities.contains(capability)

    /** Check if this lab provides WGS services */
    def providesWgs: Boolean = hasCapability(Capability.Wgs)

    /** Check if this lab provides Y-DNA testing */
    def providesYDna: Boolean = hasCapability(Capability.YDna)

    /** Check if this lab provides mtDNA testing */
    def providesMtDna: Boolean = hasCapability(Capability.MtDna)

    /** Check if this lab provides STR testing */
    def providesStr: Boolean = hasCapability(Capability.Str)

    /** Check if this lab provides chip/array genotyping */
    def providesChip: Boolean = hasCapability(Capability.Chip)

    /** Check if this lab provides VCF downloads */
    def providesVcfDownload: Boolean = hasCapability(Capability.VcfDownload)

    /** Check if this lab provides BAM downloads */
    def providesBamDownload: Boolean = hasCapability(Capability.BamDownload)
  }

  private val config: Config = ConfigFactory.load("labs.conf")

  private val labsConfig: Config = if (config.hasPath("labs")) {
    config.getConfig("labs")
  } else {
    ConfigFactory.empty()
  }

  /** All configured labs, keyed by ID */
  val labs: Map[String, Lab] = {
    labsConfig.root().keySet().asScala.flatMap { id =>
      Try {
        val labConfig = labsConfig.getConfig(id)
        val displayName = labConfig.getString("display-name")
        val abbreviation = if (labConfig.hasPath("abbreviation")) {
          labConfig.getString("abbreviation")
        } else {
          // Default: first 6 chars of display name, uppercase
          displayName.take(6).toUpperCase
        }
        val category = if (labConfig.hasPath("category")) {
          labConfig.getString("category")
        } else {
          "unknown"
        }
        val capabilities = if (labConfig.hasPath("capabilities")) {
          labConfig.getStringList("capabilities").asScala.toSet
        } else {
          Set.empty[String]
        }
        val website = if (labConfig.hasPath("website")) {
          Some(labConfig.getString("website"))
        } else {
          None
        }
        val aliases = if (labConfig.hasPath("aliases")) {
          labConfig.getStringList("aliases").asScala.toList
        } else {
          List.empty
        }

        id -> Lab(id, displayName, abbreviation, category, capabilities, website, aliases)
      }.toOption
    }.toMap
  }

  /** Lookup index by display name (case-insensitive) */
  private val byDisplayName: Map[String, Lab] =
    labs.values.map(lab => lab.displayName.toLowerCase -> lab).toMap

  /** Lookup index by alias (case-insensitive) */
  private val byAlias: Map[String, Lab] = labs.values.flatMap { lab =>
    lab.aliases.map(alias => alias.toLowerCase -> lab)
  }.toMap

  /**
   * Find a lab by any identifier: ID, display name, or alias.
   * Case-insensitive matching.
   *
   * @param identifier The lab identifier to search for
   * @return Some(Lab) if found, None otherwise
   */
  def findLab(identifier: String): Option[Lab] = {
    val lower = identifier.toLowerCase
    labs.get(lower)
      .orElse(byDisplayName.get(lower))
      .orElse(byAlias.get(lower))
      .orElse {
        // Fuzzy match: check if any display name or alias contains the identifier
        labs.values.find { lab =>
          lab.displayName.toLowerCase.contains(lower) ||
            lab.aliases.exists(_.toLowerCase.contains(lower))
        }
      }
  }

  /**
   * Get the abbreviation for a lab/vendor name.
   * If the lab is not found in config, returns a default abbreviation
   * (first 6 chars, uppercase).
   *
   * @param name     The lab name to look up
   * @param maxChars Maximum characters for the abbreviation (default 6)
   * @return The configured or default abbreviation
   */
  def getAbbreviation(name: String, maxChars: Int = 6): String = {
    findLab(name) match {
      case Some(lab) => lab.abbreviationTruncated(maxChars)
      case None => name.take(maxChars).toUpperCase
    }
  }

  /**
   * Get the full display name for a lab identifier.
   * Returns the original identifier if not found.
   *
   * @param identifier The lab identifier
   * @return The display name or the original identifier
   */
  def getDisplayName(identifier: String): String = {
    findLab(identifier).map(_.displayName).getOrElse(identifier)
  }

  // ============================================================================
  // Category-based filtering
  // ============================================================================

  /**
   * Get all labs in a specific category.
   *
   * @param category The category to filter by
   * @return List of labs in that category
   */
  def byCategory(category: String): List[Lab] =
    labs.values.filter(_.category == category).toList.sortBy(_.displayName)

  /** All commercial testing labs */
  def commercialLabs: List[Lab] = byCategory("commercial-lab")

  /** All consumer genotyping vendors */
  def consumerVendors: List[Lab] = byCategory("consumer-vendor")

  /** All sequencing platform vendors */
  def sequencingPlatforms: List[Lab] = byCategory("sequencing-platform")

  /** All academic institutions */
  def academicInstitutions: List[Lab] = byCategory("academic")

  // ============================================================================
  // Capability-based filtering
  // ============================================================================

  /**
   * Get all labs with a specific capability.
   *
   * @param capability The capability to filter by
   * @return List of labs with that capability
   */
  def withCapability(capability: String): List[Lab] =
    labs.values.filter(_.hasCapability(capability)).toList.sortBy(_.displayName)

  /** Labs that provide WGS services */
  def wgsProviders: List[Lab] = withCapability(Capability.Wgs)

  /** Labs that provide Y-DNA testing */
  def yDnaProviders: List[Lab] = withCapability(Capability.YDna)

  /** Labs that provide mtDNA testing */
  def mtDnaProviders: List[Lab] = withCapability(Capability.MtDna)

  /** Labs that provide STR testing */
  def strProviders: List[Lab] = withCapability(Capability.Str)

  /** Labs that provide chip/array genotyping */
  def chipProviders: List[Lab] = withCapability(Capability.Chip)

  /** Labs that provide VCF downloads */
  def vcfDownloadProviders: List[Lab] = withCapability(Capability.VcfDownload)

  // ============================================================================
  // Dialog-specific lab lists (display names for ComboBox items)
  // ============================================================================

  /**
   * Get lab display names for sequence run editing.
   * Includes commercial labs, sequencing platforms, and academic institutions.
   */
  def sequenceRunLabNames: List[String] =
    (commercialLabs ++ sequencingPlatforms ++ academicInstitutions)
      .distinctBy(_.displayName)
      .sortBy(_.displayName)
      .map(_.displayName)

  /**
   * Get lab display names for subject/biosample center selection.
   * Includes all labs.
   */
  def allLabNames: List[String] =
    labs.values.toList.sortBy(_.displayName).map(_.displayName)

  /**
   * Get lab display names for STR profile imports.
   * Returns labs that provide STR testing.
   */
  def strProviderNames: List[String] =
    strProviders.map(_.displayName)

  /**
   * Get lab display names for Y-DNA profile sources.
   * Returns labs that provide Y-DNA testing (includes chip-based).
   */
  def yDnaProviderNames: List[String] =
    (yDnaProviders ++ chipProviders.filter(_.providesYDna))
      .distinctBy(_.displayName)
      .sortBy(_.displayName)
      .map(_.displayName)

  /**
   * Get lab display names for VCF imports.
   * Returns labs that provide VCF downloads.
   */
  def vcfImportProviderNames: List[String] =
    vcfDownloadProviders.map(_.displayName)

  /**
   * Get lab display names for chip/array data.
   * Returns consumer vendors and labs that provide chip genotyping.
   */
  def chipProviderNames: List[String] =
    (consumerVendors ++ chipProviders)
      .distinctBy(_.displayName)
      .sortBy(_.displayName)
      .map(_.displayName)
}
