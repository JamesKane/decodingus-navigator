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
   * Tree provider configuration for Y-DNA and MT-DNA analysis.
   */
  object treeProviders {
    private val treeConfig = if (config.hasPath("tree-providers")) {
      config.getConfig("tree-providers")
    } else {
      ConfigFactory.empty()
    }

    /** Y-DNA tree provider: "ftdna" or "decodingus" */
    val ydna: String = if (treeConfig.hasPath("ydna")) treeConfig.getString("ydna") else "ftdna"

    /** MT-DNA tree provider: "ftdna" or "decodingus" */
    val mtdna: String = if (treeConfig.hasPath("mtdna")) treeConfig.getString("mtdna") else "ftdna"
  }

  /**
   * Ancestry analysis configuration for population percentage estimation.
   */
  object ancestryAnalysis {
    private val ancestryConfig = if (config.hasPath("ancestry-analysis")) {
      config.getConfig("ancestry-analysis")
    } else {
      ConfigFactory.empty()
    }

    /** Whether ancestry analysis is enabled */
    val enabled: Boolean =
      if (ancestryConfig.hasPath("enabled")) ancestryConfig.getBoolean("enabled") else true

    /** Default panel type: "aims" or "genome-wide" */
    val defaultPanel: String =
      if (ancestryConfig.hasPath("default-panel")) ancestryConfig.getString("default-panel") else "aims"

    /** Minimum SNPs required for AIMs panel analysis */
    val minSnpsAims: Int =
      if (ancestryConfig.hasPath("min-snps-aims")) ancestryConfig.getInt("min-snps-aims") else 3000

    /** Minimum SNPs required for genome-wide analysis */
    val minSnpsGenomeWide: Int =
      if (ancestryConfig.hasPath("min-snps-genome-wide")) ancestryConfig.getInt("min-snps-genome-wide") else 100000

    /** Minimum percentage to display in results */
    val displayThreshold: Double =
      if (ancestryConfig.hasPath("display-threshold")) ancestryConfig.getDouble("display-threshold") else 0.5

    /** Reference data version */
    val referenceVersion: String =
      if (ancestryConfig.hasPath("reference-version")) ancestryConfig.getString("reference-version") else "v1"
  }

  /**
   * Chip/Array genotype data configuration.
   */
  object chipData {
    private val chipConfig = if (config.hasPath("chip-data")) {
      config.getConfig("chip-data")
    } else {
      ConfigFactory.empty()
    }

    /** Whether chip data import is enabled */
    val enabled: Boolean =
      if (chipConfig.hasPath("enabled")) chipConfig.getBoolean("enabled") else true

    /** Minimum marker count to accept a chip file */
    val minMarkerCount: Int =
      if (chipConfig.hasPath("min-marker-count")) chipConfig.getInt("min-marker-count") else 100000

    /** Maximum acceptable no-call rate (0.05 = 5%) */
    val maxNoCallRate: Double =
      if (chipConfig.hasPath("max-no-call-rate")) chipConfig.getDouble("max-no-call-rate") else 0.05

    /** Minimum Y-DNA markers for haplogroup estimation */
    val minYMarkers: Int =
      if (chipConfig.hasPath("min-y-markers")) chipConfig.getInt("min-y-markers") else 50

    /** Minimum mtDNA markers for haplogroup estimation */
    val minMtMarkers: Int =
      if (chipConfig.hasPath("min-mt-markers")) chipConfig.getInt("min-mt-markers") else 20

    /** Supported vendors list */
    val supportedVendors: List[String] =
      if (chipConfig.hasPath("supported-vendors"))
        chipConfig.getStringList("supported-vendors").asScala.toList
      else
        List("23andMe", "AncestryDNA", "FamilyTreeDNA", "MyHeritage", "LivingDNA")
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
