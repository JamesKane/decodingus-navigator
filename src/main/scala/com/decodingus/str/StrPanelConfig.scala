package com.decodingus.str

/**
 * STR panel configuration data models.
 * Loaded from str-panels.conf HOCON configuration.
 */

/**
 * Root configuration for STR panels.
 *
 * @param version             Configuration version
 * @param providers           Map of provider key to provider definition
 * @param panels              List of panel definitions
 * @param multiValueMarkers   Map of marker name to multi-value info
 */
case class StrPanelConfig(
  version: String,
  providers: Map[String, StrProviderDef],
  panels: List[StrPanelDef],
  multiValueMarkers: Map[String, MultiValueMarkerDef]
) {
  /** Get panels for a specific provider, sorted by order */
  def panelsForProvider(provider: String): List[StrPanelDef] =
    panels.filter(_.provider == provider).sortBy(_.order)

  /** Get provider definition by key */
  def provider(key: String): Option[StrProviderDef] = providers.get(key)

  /** Get all FTDNA panels sorted by order */
  def ftdnaPanels: List[StrPanelDef] = panelsForProvider("FTDNA")

  /** Get all YSEQ panels sorted by order */
  def yseqPanels: List[StrPanelDef] = panelsForProvider("YSEQ")
}

/**
 * Provider definition with exclusive markers for detection.
 *
 * @param displayName       Human-readable provider name
 * @param cumulativePanels  Whether panels are cumulative (Y-25 includes Y-12 markers)
 * @param exclusiveMarkers  Markers unique to this provider (for auto-detection)
 */
case class StrProviderDef(
  displayName: String,
  cumulativePanels: Boolean,
  exclusiveMarkers: List[String]
)

/**
 * Panel definition with both marketing and actual marker counts.
 *
 * @param id               Unique identifier (e.g., "FTDNA_Y37")
 * @param name             Display name (e.g., "Y-37")
 * @param provider         Provider key (e.g., "FTDNA")
 * @param marketingCount   What vendor advertises (37)
 * @param actualCount      Distinct marker keys (30), None if unknown
 * @param order            Sort order for cumulative panel determination
 * @param markers          List of marker names in this panel section
 */
case class StrPanelDef(
  id: String,
  name: String,
  provider: String,
  marketingCount: Int,
  actualCount: Option[Int],
  order: Int,
  markers: List[String]
) {
  /** Get the threshold to use for classification (actual if known, else marketing) */
  def threshold: Int = actualCount.getOrElse(marketingCount)
}

/**
 * Multi-value marker definition.
 * Some markers (like DYS464) have multiple allele copies that count
 * as multiple markers in marketing materials.
 *
 * @param copies         Number of actual copies (e.g., DYS464 has 4)
 * @param marketingCount How many markers this counts as in marketing
 */
case class MultiValueMarkerDef(
  copies: Int,
  marketingCount: Int
)

/**
 * Result of panel classification.
 */
case class PanelClassificationResult(
  detectedPanel: Option[StrPanelDef],
  markerCount: Int,
  provider: Option[String],
  matchedMarkers: Int = 0,
  error: Option[String] = None
) {
  /** Marketing panel name (e.g., "Y-111") */
  def panelName: Option[String] = detectedPanel.map(_.name)

  /** Marketing marker count (e.g., 111) */
  def marketingCount: Option[Int] = detectedPanel.map(_.marketingCount)

  /** Actual marker count threshold (e.g., 102) */
  def actualThreshold: Option[Int] = detectedPanel.map(_.threshold)
}
