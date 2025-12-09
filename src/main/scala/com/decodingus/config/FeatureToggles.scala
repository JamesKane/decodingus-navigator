package com.decodingus.config

import com.typesafe.config.ConfigFactory
import scala.jdk.CollectionConverters.*

object FeatureToggles {
  private val config = ConfigFactory.load("feature_toggles.conf")
  val pdsSubmissionEnabled: Boolean = config.getBoolean("pds-submission.enabled")
  val authEnabled: Boolean = config.hasPath("auth.enabled") && config.getBoolean("auth.enabled")
  val atProtocolEnabled: Boolean = config.hasPath("at-protocol.enabled") && config.getBoolean("at-protocol.enabled")

  object developerFeatures {
    private val devConfig = config.getConfig("developer-features")
    val saveJsonEnabled: Boolean = devConfig.getBoolean("save-json-enabled")
  }

  /**
   * Reference genome haplogroup mappings for Y-DNA calling optimization.
   * Maps reference build names to their known Y-DNA haplogroup name variants.
   * Multiple name variants are supported for compatibility with different tree providers.
   */
  object referenceHaplogroups {
    private val refConfig = if (config.hasPath("reference-haplogroups")) {
      config.getConfig("reference-haplogroups")
    } else {
      ConfigFactory.empty()
    }

    private val mappings: Map[String, List[String]] = {
      refConfig.entrySet().asScala.map { entry =>
        val key = entry.getKey.stripPrefix("\"").stripSuffix("\"")
        val value = entry.getValue.unwrapped().toString
        // Split comma-separated name variants
        key -> value.split(",").map(_.trim).toList
      }.toMap
    }

    /**
     * Get the known Y-DNA haplogroup name variants for a reference genome build.
     * Returns multiple name variants for compatibility with different tree providers
     * (e.g., "R1b-U152" for Decoding-Us, "R-U152" for FTDNA).
     *
     * @param referenceBuild The reference build name (e.g., "GRCh38", "T2T-CHM13")
     * @return Some(List of haplogroup name variants) if known, None otherwise
     */
    def getHaplogroups(referenceBuild: String): Option[List[String]] = {
      // Try exact match first, then case-insensitive
      mappings.get(referenceBuild)
        .orElse(mappings.find { case (k, _) => k.equalsIgnoreCase(referenceBuild) }.map(_._2))
    }

    /**
     * Check if we have haplogroup info for a reference build.
     */
    def hasHaplogroup(referenceBuild: String): Boolean = getHaplogroups(referenceBuild).isDefined
  }
}
