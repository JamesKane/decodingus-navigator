package com.decodingus.workspace.model

import com.decodingus.ancestry.model.{AncestryResult, PopulationComponent, SuperPopulationSummary}

import java.time.LocalDateTime

/**
 * First-class Atmosphere Lexicon record for ancestry composition analysis.
 * NSID: com.decodingus.atmosphere.populationBreakdown
 *
 * Wraps the analysis-layer AncestryResult with AT Protocol record metadata
 * (biosample reference, analysis method, record meta) for persistence and sync.
 */
case class PopulationBreakdown(
                                meta: RecordMeta,
                                biosampleRef: String, // AT URI or local ref of parent biosample
                                analysisMethod: String, // e.g. "PCA_PROJECTION_GMM"
                                panelType: String, // "aims" or "genome-wide"
                                referencePopulations: String, // e.g. "1000G_HGDP_v1"
                                snpsAnalyzed: Int,
                                snpsWithGenotype: Int,
                                snpsMissing: Int,
                                confidenceLevel: Double,
                                components: List[PopulationComponent],
                                superPopulationSummary: List[SuperPopulationSummary],
                                pcaCoordinates: Option[List[Double]],
                                analysisDate: Option[LocalDateTime],
                                pipelineVersion: String,
                                referenceVersion: String
                              )

object PopulationBreakdown:

  /**
   * Create a PopulationBreakdown from an AncestryResult.
   */
  def fromAncestryResult(
                           biosampleRef: String,
                           result: AncestryResult,
                           analysisMethod: String = "PCA_PROJECTION_GMM",
                           referencePopulations: String = "1000G_HGDP_v1"
                         ): PopulationBreakdown =
    PopulationBreakdown(
      meta = RecordMeta.initial,
      biosampleRef = biosampleRef,
      analysisMethod = analysisMethod,
      panelType = result.panelType,
      referencePopulations = referencePopulations,
      snpsAnalyzed = result.snpsAnalyzed,
      snpsWithGenotype = result.snpsWithGenotype,
      snpsMissing = result.snpsMissing,
      confidenceLevel = result.confidenceLevel,
      components = result.components,
      superPopulationSummary = result.superPopulationSummary,
      pcaCoordinates = result.pcaCoordinates,
      analysisDate = Some(LocalDateTime.now()),
      pipelineVersion = result.pipelineVersion,
      referenceVersion = result.referenceVersion
    )

  /**
   * Convert back to an AncestryResult for display/reporting.
   */
  def toAncestryResult(breakdown: PopulationBreakdown): AncestryResult =
    AncestryResult(
      panelType = breakdown.panelType,
      snpsAnalyzed = breakdown.snpsAnalyzed,
      snpsWithGenotype = breakdown.snpsWithGenotype,
      snpsMissing = breakdown.snpsMissing,
      components = breakdown.components,
      superPopulationSummary = breakdown.superPopulationSummary,
      confidenceLevel = breakdown.confidenceLevel,
      pipelineVersion = breakdown.pipelineVersion,
      referenceVersion = breakdown.referenceVersion,
      pcaCoordinates = breakdown.pcaCoordinates
    )
