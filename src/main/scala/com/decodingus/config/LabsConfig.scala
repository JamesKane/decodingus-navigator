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
 */
object LabsConfig {

  /**
   * Represents a lab, sequencing center, or vendor configuration.
   *
   * @param id           Unique identifier (the HOCON key)
   * @param displayName  Full display name
   * @param abbreviation Short code (up to 6 chars) for UI indicators
   * @param category     Category: commercial-lab, consumer-vendor, sequencing-platform, academic
   * @param website      Optional website URL
   * @param aliases      Alternative names for matching
   */
  case class Lab(
    id: String,
    displayName: String,
    abbreviation: String,
    category: String,
    website: Option[String],
    aliases: List[String]
  ) {
    /** Returns abbreviation truncated to maxLength chars (default 6) */
    def abbreviationTruncated(maxLength: Int = 6): String = abbreviation.take(maxLength)
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

        id -> Lab(id, displayName, abbreviation, category, website, aliases)
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

  /**
   * Get all labs in a specific category.
   *
   * @param category The category to filter by
   * @return List of labs in that category
   */
  def byCategory(category: String): List[Lab] =
    labs.values.filter(_.category == category).toList

  /** All commercial testing labs */
  def commercialLabs: List[Lab] = byCategory("commercial-lab")

  /** All consumer genotyping vendors */
  def consumerVendors: List[Lab] = byCategory("consumer-vendor")

  /** All sequencing platform vendors */
  def sequencingPlatforms: List[Lab] = byCategory("sequencing-platform")

  /** All academic institutions */
  def academicInstitutions: List[Lab] = byCategory("academic")
}
