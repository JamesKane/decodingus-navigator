package com.decodingus.model

/**
 * Represents the analysis summary for a single contig.
 *
 * @param contigName The name of the contig (e.g., "chr1").
 * @param refN The number of bases that are 'N' in the reference.
 * @param callable The number of bases classified as CALLABLE.
 * @param noCoverage The number of bases with NO_COVERAGE.
 * @param lowCoverage The number of bases with LOW_COVERAGE.
 * @param excessiveCoverage The number of bases with EXCESSIVE_COVERAGE.
 * @param poorMappingQuality The number of bases with POOR_MAPPING_QUALITY.
 */
case class ContigSummary(
  contigName: String,
  refN: Long,
  callable: Long,
  noCoverage: Long,
  lowCoverage: Long,
  excessiveCoverage: Long,
  poorMappingQuality: Long
)

/**
 * Represents the final JSON report for the callable loci analysis.
 *
 * @param pdsUserId The user's Personal Data Store ID.
 * @param libraryStats The statistics gathered from the BAM/CRAM library.
 * @param wgsMetrics The whole genome sequencing metrics from GATK.
 * @param callableBases The total number of callable bases across all contigs.
 * @param contigAnalysis A list of summaries for each processed contig.
 */
case class CoverageSummary(
  pdsUserId: String,
  libraryStats: LibraryStats,
  wgsMetrics: WgsMetrics,
  callableBases: Long,
  contigAnalysis: List[ContigSummary]
)