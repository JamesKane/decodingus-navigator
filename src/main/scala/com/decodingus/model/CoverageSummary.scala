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
 * @param platformSource The sequencing platform source (e.g., "bwa-mem2").
 * @param reference The reference genome used (e.g., "T2T-CHM13v2.0").
 * @param totalReads The total number of reads in the input file.
 * @param readLength The average read length.
 * @param totalBases The total number of bases analyzed across all contigs.
 * @param callableBases The total number of callable bases across all contigs.
 * @param averageDepth The estimated average sequencing depth.
 * @param contigAnalysis A list of summaries for each processed contig.
 */
case class CoverageSummary(
  pdsUserId: String,
  platformSource: String,
  reference: String,
  totalReads: Long,
  readLength: Int,
  totalBases: Long,
  callableBases: Long,
  averageDepth: Double,
  contigAnalysis: List[ContigSummary]
)
