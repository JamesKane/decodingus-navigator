package com.decodingus.workspace.model

/**
 * An identical-by-descent (IBD) segment shared between two samples.
 *
 * @param chromosome      Chromosome number (e.g., "1", "22", "X")
 * @param startPosition   Start position in base pairs
 * @param endPosition     End position in base pairs
 * @param lengthCm        Length in centiMorgans
 * @param snpCount        Number of SNPs in the segment
 * @param isHalfIdentical True if half-identical (one allele matches), false if fully identical
 */
case class IbdSegment(
                       chromosome: String,
                       startPosition: Int,
                       endPosition: Int,
                       lengthCm: Double,
                       snpCount: Option[Int] = None,
                       isHalfIdentical: Option[Boolean] = None
                     )
