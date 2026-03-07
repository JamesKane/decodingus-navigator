package com.decodingus.model

import com.decodingus.workspace.model.ContigMetrics

/**
 * Represents the final JSON report for the callable loci analysis.
 *
 * @param pdsUserId      The user's Personal Data Store ID.
 * @param libraryStats   The statistics gathered from the BAM/CRAM library.
 * @param wgsMetrics     The whole genome sequencing metrics from GATK.
 * @param callableBases  The total number of callable bases across all contigs.
 * @param contigAnalysis A list of metrics for each processed contig.
 */
case class CoverageSummary(
                            pdsUserId: String,
                            biosampleId: String,
                            libraryStats: LibraryStats,
                            wgsMetrics: WgsMetrics,
                            callableBases: Long,
                            contigAnalysis: List[ContigMetrics]
                          )