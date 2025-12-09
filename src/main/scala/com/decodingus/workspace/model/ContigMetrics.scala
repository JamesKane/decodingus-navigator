package com.decodingus.workspace.model

/**
 * Metrics for a single contig from callable loci analysis.
 * Compatible with GATK CallableLoci output.
 */
case class ContigMetrics(
  contigName: String,
  callable: Long = 0,
  noCoverage: Long = 0,
  lowCoverage: Long = 0,
  excessiveCoverage: Long = 0,
  poorMappingQuality: Long = 0,
  refN: Long = 0,
  meanCoverage: Option[Double] = None
) {
  /** Alias for backward compatibility */
  def callableBases: Long = callable
}
