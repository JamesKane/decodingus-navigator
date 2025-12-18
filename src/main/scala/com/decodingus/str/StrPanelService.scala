package com.decodingus.str

import com.typesafe.config.{Config, ConfigFactory}

import scala.jdk.CollectionConverters.*
import scala.util.Try

/**
 * Service for loading STR panel configuration and classifying profiles.
 *
 * Key responsibilities:
 * - Load panel definitions from bundled HOCON resource
 * - Detect provider from exclusive markers
 * - Classify profile into highest matching panel tier
 * - Calculate cumulative marker counts for panel matching
 */
object StrPanelService {

  // Lazy-loaded and cached configuration
  @volatile private var cachedConfig: Option[StrPanelConfig] = None

  /**
   * Get the STR panel configuration, loading from resources if needed.
   */
  def getConfig: Either[String, StrPanelConfig] = {
    cachedConfig match {
      case Some(config) => Right(config)
      case None => loadConfig()
    }
  }

  /**
   * Load configuration from bundled resource.
   */
  private def loadConfig(): Either[String, StrPanelConfig] = {
    Try {
      val config = ConfigFactory.load("str-panels.conf")
      val strPanels = config.getConfig("str-panels")
      parseConfig(strPanels)
    }.toEither.left.map(_.getMessage).flatten match {
      case Right(config) =>
        cachedConfig = Some(config)
        Right(config)
      case Left(error) =>
        Left(s"Failed to load STR panel config: $error")
    }
  }

  /**
   * Parse the HOCON config into typed case classes.
   */
  private def parseConfig(config: Config): Either[String, StrPanelConfig] = {
    Try {
      val version = config.getString("version")

      // Parse providers
      val providersConfig = config.getConfig("providers")
      val providers = providersConfig.root().keySet().asScala.map { key =>
        val pConfig = providersConfig.getConfig(key)
        val providerDef = StrProviderDef(
          displayName = pConfig.getString("display-name"),
          cumulativePanels = pConfig.getBoolean("cumulative-panels"),
          exclusiveMarkers = pConfig.getStringList("exclusive-markers").asScala.toList
        )
        key -> providerDef
      }.toMap

      // Parse panels
      val panelsList = config.getConfigList("panels").asScala.toList
      val panels = panelsList.map { pConfig =>
        StrPanelDef(
          id = pConfig.getString("id"),
          name = pConfig.getString("name"),
          provider = pConfig.getString("provider"),
          marketingCount = pConfig.getInt("marketing-count"),
          actualCount = if (pConfig.hasPath("actual-count") && !pConfig.getIsNull("actual-count"))
            Some(pConfig.getInt("actual-count"))
          else
            None,
          order = pConfig.getInt("order"),
          markers = pConfig.getStringList("markers").asScala.toList
        )
      }

      // Parse multi-value markers
      val multiValueConfig = config.getConfig("multi-value-markers")
      val multiValueMarkers = multiValueConfig.root().keySet().asScala.map { key =>
        val mConfig = multiValueConfig.getConfig(key)
        val markerDef = MultiValueMarkerDef(
          copies = mConfig.getInt("copies"),
          marketingCount = mConfig.getInt("marketing-count")
        )
        key -> markerDef
      }.toMap

      StrPanelConfig(version, providers, panels, multiValueMarkers)
    }.toEither.left.map(_.getMessage)
  }

  /**
   * Detect provider from marker set using exclusive markers.
   *
   * @param markers Set of marker names (already normalized to uppercase at import)
   * @return Detected provider key or None
   */
  def detectProvider(markers: Set[String]): Option[String] = {
    getConfig.toOption.flatMap { config =>
      // Check each provider's exclusive markers
      config.providers.collectFirst {
        case (providerKey, providerDef)
          if providerDef.exclusiveMarkers.nonEmpty &&
            providerDef.exclusiveMarkers.exists(m => markers.contains(m)) =>
          providerKey
      }
    }
  }

  /**
   * Classify a profile into the highest matching panel tier.
   *
   * For cumulative panels (FTDNA), we calculate cumulative marker counts
   * and find the highest panel where the profile has enough markers.
   *
   * @param markers  Set of marker names in the profile
   * @param provider Provider key (if known), used to filter panels
   * @return Panel classification result
   */
  def classifyPanel(
    markers: Set[String],
    provider: Option[String]
  ): PanelClassificationResult = {
    getConfig match {
      case Left(error) =>
        PanelClassificationResult(
          detectedPanel = None,
          markerCount = markers.size,
          provider = provider,
          error = Some(error)
        )
      case Right(config) =>
        // Detect provider from exclusive markers if not provided
        val effectiveProvider = provider.orElse(detectProvider(markers))

        // Get panels for this provider, sorted by order
        val providerPanels = effectiveProvider match {
          case Some(p) => config.panels.filter(_.provider == p).sortBy(_.order)
          case None => config.panels.filter(_.provider == "FTDNA").sortBy(_.order) // Default to FTDNA
        }

        // For cumulative panels, calculate running total of markers
        val providerDef = effectiveProvider.flatMap(config.providers.get)
        val isCumulative = providerDef.exists(_.cumulativePanels)

        val (matchedPanel, matchedCount) = if (isCumulative) {
          findHighestCumulativePanel(markers, providerPanels)
        } else {
          findBestMatchingPanel(markers, providerPanels)
        }

        PanelClassificationResult(
          detectedPanel = matchedPanel,
          markerCount = markers.size,
          provider = effectiveProvider,
          matchedMarkers = matchedCount
        )
    }
  }

  /**
   * Find highest panel where cumulative marker threshold is met.
   * Returns the panel and the count of matched markers.
   */
  private def findHighestCumulativePanel(
    markers: Set[String],
    panels: List[StrPanelDef]
  ): (Option[StrPanelDef], Int) = {
    var cumulativeMarkers = Set.empty[String]
    var highestMatched: Option[StrPanelDef] = None
    var matchedCount = 0

    for (panel <- panels) {
      // Add this panel's markers to cumulative set
      cumulativeMarkers = cumulativeMarkers ++ panel.markers

      // Count how many of profile's markers match cumulative set
      val matchCount = markers.count(m => cumulativeMarkers.contains(m))

      // Use actualCount for threshold (not marketing count)
      val threshold = panel.actualCount.getOrElse(panel.marketingCount)

      // Require at least 90% of threshold markers to match
      // This accounts for slight naming variations
      if (matchCount >= (threshold * 0.9).toInt && matchCount > matchedCount) {
        highestMatched = Some(panel)
        matchedCount = matchCount
      }
    }

    (highestMatched, matchedCount)
  }

  /**
   * Find best matching panel for non-cumulative providers.
   */
  private def findBestMatchingPanel(
    markers: Set[String],
    panels: List[StrPanelDef]
  ): (Option[StrPanelDef], Int) = {
    val matchedPanels = panels.map { panel =>
      val panelMarkers = panel.markers.toSet
      val overlapCount = markers.intersect(panelMarkers).size
      (panel, overlapCount)
    }.filter { case (panel, count) =>
      count >= (panel.markers.size * 0.8)
    }

    matchedPanels.maxByOption(_._1.order) match {
      case Some((panel, count)) => (Some(panel), count)
      case None => (None, 0)
    }
  }

  /**
   * Get panels for a specific provider, sorted by order.
   */
  def getPanelsForProvider(provider: String): List[StrPanelDef] = {
    getConfig.toOption.map(_.panelsForProvider(provider)).getOrElse(Nil)
  }

  /**
   * Get FTDNA panels with their actual thresholds for UI display.
   * Returns tuples of (displayName, actualThreshold, marketingCount).
   */
  def getFtdnaPanelThresholds: List[(String, Int, Int)] = {
    getConfig.toOption.map { config =>
      config.ftdnaPanels.map { p =>
        (p.name, p.threshold, p.marketingCount)
      }
    }.getOrElse {
      // Fallback hardcoded values if config fails
      List(
        ("Y-12", 11, 12),
        ("Y-25", 20, 25),
        ("Y-37", 30, 37),
        ("Y-67", 58, 67),
        ("Y-111", 102, 111),
        ("Y-500", 450, 500),
        ("Y-700", 630, 700)
      )
    }
  }

  /**
   * Get all provider definitions.
   */
  def getProviders: Map[String, StrProviderDef] = {
    getConfig.toOption.map(_.providers).getOrElse(Map.empty)
  }

  /**
   * Get list of available providers that have defined panels.
   * Filters out providers with only placeholder/empty panels.
   */
  def getAvailableProviders: List[String] = {
    getConfig.toOption.map { config =>
      config.providers.keys.toList.filter { provider =>
        // Provider is available if it has at least one panel with markers or defined thresholds
        config.panelsForProvider(provider).exists { panel =>
          panel.markers.nonEmpty || panel.actualCount.isDefined
        }
      }.sorted
    }.getOrElse(List("FTDNA"))
  }

  /**
   * Get panel thresholds for a specific provider.
   * Returns tuples of (displayName, actualThreshold, marketingCount).
   */
  def getPanelThresholdsForProvider(provider: String): List[(String, Int, Int)] = {
    getConfig.toOption.map { config =>
      config.panelsForProvider(provider)
        .filter(p => p.markers.nonEmpty || p.actualCount.isDefined) // Only panels with data
        .map { p =>
          (p.name, p.threshold, p.marketingCount)
        }
    }.getOrElse {
      if (provider == "FTDNA") getFtdnaPanelThresholds
      else Nil
    }
  }

  /**
   * Check if a provider has panel definitions with actual marker data.
   */
  def providerHasPanelData(provider: String): Boolean = {
    getConfig.toOption.exists { config =>
      config.panelsForProvider(provider).exists(_.markers.nonEmpty)
    }
  }

  /**
   * Get all known panel names for validation.
   */
  def getKnownPanelNames: Set[String] = {
    getConfig.toOption.map(_.panels.map(_.name).toSet).getOrElse(Set.empty)
  }

  /**
   * Get exclusive markers for a provider.
   */
  def getExclusiveMarkers(provider: String): Set[String] = {
    getConfig.toOption
      .flatMap(_.providers.get(provider))
      .map(_.exclusiveMarkers.toSet)
      .getOrElse(Set.empty)
  }

  /**
   * Clear cached configuration (useful for testing).
   */
  def clearCache(): Unit = {
    cachedConfig = None
  }
}
