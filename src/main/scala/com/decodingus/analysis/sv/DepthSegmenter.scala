package com.decodingus.analysis.sv

import scala.collection.mutable

/**
 * Segments depth bins into copy number variant calls using z-score analysis.
 *
 * Algorithm:
 * 1. Calculate expected depth per bin based on mean coverage and bin size
 * 2. Compute z-score for each bin using Poisson-like variance (sqrt(expected))
 * 3. Identify bins with significant deviation (|z| > threshold)
 * 4. Merge adjacent aberrant bins into segments
 * 5. Filter segments by minimum size
 *
 * References:
 * - Z-score based CNV detection: Xie & Tammi. "CNV-seq, a new method to detect copy
 *   number variation using high-throughput sequencing." BMC Bioinformatics 10.1 (2009): 80.
 *   https://doi.org/10.1186/1471-2105-10-80
 *
 * - Segmentation approach inspired by: Olshen et al. "Circular binary segmentation for
 *   the analysis of array-based DNA copy number data." Biostatistics 5.4 (2004): 557-572.
 *   https://doi.org/10.1093/biostatistics/kxh008
 *
 * - Read depth normalization: Abyzov et al. "CNVnator: an approach to discover, genotype,
 *   and characterize typical and atypical CNVs from family and population genome sequencing."
 *   Genome Research 21.6 (2011): 974-984.
 *   https://doi.org/10.1101/gr.114876.110
 */
class DepthSegmenter(config: SvCallerConfig = SvCallerConfig.default) {

  private val BIN_SIZE = config.binSize
  private val MIN_Z_SCORE = config.minDepthZScore
  private val MIN_CNV_SIZE = config.minCnvSize

  /**
   * Segment depth bins into CNV calls.
   *
   * @param depthBins      Map of contig to read counts per bin
   * @param contigLengths  Map of contig to total length
   * @param meanCoverage   Overall mean coverage of the sample
   * @param readLength     Mean read length (for expected calculation)
   * @return List of depth segments representing CNVs
   */
  def segment(
    depthBins: Map[String, Array[Int]],
    contigLengths: Map[String, Long],
    meanCoverage: Double,
    readLength: Double = 150.0
  ): List[DepthSegment] = {

    // Calculate expected reads per bin based on coverage
    // Expected = (coverage * binSize) / readLength
    val expectedReadsPerBin = (meanCoverage * BIN_SIZE) / readLength

    val segments = mutable.ListBuffer[DepthSegment]()

    depthBins.foreach { case (contig, bins) =>
      val contigLength = contigLengths.getOrElse(contig, bins.length.toLong * BIN_SIZE)

      // Calculate z-scores for each bin
      val zScores = bins.map { readCount =>
        if (expectedReadsPerBin > 0) {
          // Use Poisson-like variance approximation
          val variance = math.max(expectedReadsPerBin, 1.0)
          (readCount - expectedReadsPerBin) / math.sqrt(variance)
        } else {
          0.0
        }
      }

      // Find runs of aberrant bins
      var i = 0
      while (i < zScores.length) {
        val z = zScores(i)

        if (math.abs(z) >= MIN_Z_SCORE) {
          // Start of an aberrant region
          val isDelection = z < 0
          val startBin = i
          var endBin = i
          var sumZ = z
          var sumDepth = bins(i).toDouble
          var count = 1

          // Extend the segment while z-score has same sign and magnitude > threshold
          // Allow some tolerance for noise (require 2/3 of bins to be aberrant)
          while (endBin + 1 < zScores.length) {
            val nextZ = zScores(endBin + 1)
            val sameSide = (isDelection && nextZ < -MIN_Z_SCORE * 0.5) ||
                          (!isDelection && nextZ > MIN_Z_SCORE * 0.5)

            if (sameSide) {
              endBin += 1
              sumZ += nextZ
              sumDepth += bins(endBin)
              count += 1
            } else {
              // Check if we should continue through a dip
              val lookAhead = math.min(endBin + 3, zScores.length - 1)
              val futureAberrant = (endBin + 1 to lookAhead).count { j =>
                val fz = zScores(j)
                (isDelection && fz < -MIN_Z_SCORE * 0.5) ||
                (!isDelection && fz > MIN_Z_SCORE * 0.5)
              }
              if (futureAberrant >= 2) {
                // Continue through the noise
                endBin += 1
                sumZ += nextZ
                sumDepth += bins(endBin)
                count += 1
              } else {
                // End of segment
                endBin = endBin // no change, exit loop
                i = endBin // will be incremented at end
                endBin = zScores.length // exit condition
              }
            }
          }

          if (endBin >= zScores.length) endBin = count + startBin - 1

          val segmentStart = startBin.toLong * BIN_SIZE
          val segmentEnd = math.min((endBin + 1).toLong * BIN_SIZE, contigLength)
          val segmentLength = segmentEnd - segmentStart

          // Only report if segment meets minimum size
          if (segmentLength >= MIN_CNV_SIZE) {
            val meanDepthInSegment = sumDepth / count
            val meanZScore = sumZ / count

            // Calculate log2 ratio
            val log2Ratio = if (expectedReadsPerBin > 0) {
              math.log(meanDepthInSegment / expectedReadsPerBin) / math.log(2)
            } else {
              0.0
            }

            val svType = if (meanZScore < 0) SvType.DEL else SvType.DUP

            segments += DepthSegment(
              chrom = contig,
              start = segmentStart,
              end = segmentEnd,
              meanDepth = meanDepthInSegment,
              log2Ratio = log2Ratio,
              zScore = meanZScore,
              numBins = count,
              svType = svType
            )
          }

          i = startBin + count
        } else {
          i += 1
        }
      }
    }

    // Sort segments by position
    segments.toList.sortBy(s => (s.chrom, s.start))
  }

  /**
   * Merge nearby segments of the same type.
   *
   * @param segments    List of depth segments
   * @param maxGap      Maximum gap between segments to merge (in bp)
   * @return Merged segments
   */
  def mergeNearbySegments(segments: List[DepthSegment], maxGap: Int = 50000): List[DepthSegment] = {
    if (segments.isEmpty) return segments

    val sorted = segments.sortBy(s => (s.chrom, s.start))
    val merged = mutable.ListBuffer[DepthSegment]()
    var current = sorted.head

    sorted.tail.foreach { next =>
      if (current.chrom == next.chrom &&
          current.svType == next.svType &&
          next.start - current.end <= maxGap) {
        // Merge segments
        val totalBins = current.numBins + next.numBins
        val weightedDepth = (current.meanDepth * current.numBins + next.meanDepth * next.numBins) / totalBins
        val weightedZ = (current.zScore * current.numBins + next.zScore * next.numBins) / totalBins
        val weightedLog2 = (current.log2Ratio * current.numBins + next.log2Ratio * next.numBins) / totalBins

        current = DepthSegment(
          chrom = current.chrom,
          start = current.start,
          end = next.end,
          meanDepth = weightedDepth,
          log2Ratio = weightedLog2,
          zScore = weightedZ,
          numBins = totalBins,
          svType = current.svType
        )
      } else {
        merged += current
        current = next
      }
    }
    merged += current

    merged.toList
  }

  /**
   * Convert depth segments to SV calls.
   *
   * @param segments List of depth segments
   * @param config   SV caller configuration
   * @return List of SV calls
   */
  def toSvCalls(segments: List[DepthSegment]): List[SvCall] = {
    segments.zipWithIndex.map { case (seg, idx) =>
      val quality = math.min(math.abs(seg.zScore) * 10, 99.0) // Phred-like quality

      // Estimate confidence intervals based on number of bins
      val ciSize = math.max(BIN_SIZE / 2, 100)

      // Determine genotype based on log2 ratio
      val genotype = seg.svType match {
        case SvType.DEL =>
          if (seg.log2Ratio < -0.9) "1/1"  // Homozygous deletion
          else "0/1"                        // Heterozygous deletion
        case SvType.DUP =>
          if (seg.log2Ratio > 0.7) "1/1"   // High copy gain
          else "0/1"                        // Single copy gain
        case _ => "0/1"
      }

      val svLen = if (seg.svType == SvType.DEL) {
        -(seg.end - seg.start)
      } else {
        seg.end - seg.start
      }

      SvCall(
        id = s"CNV_${seg.chrom}_${seg.start}_$idx",
        chrom = seg.chrom,
        start = seg.start,
        end = seg.end,
        svType = seg.svType,
        svLen = svLen,
        ciPos = (-ciSize, ciSize),
        ciEnd = (-ciSize, ciSize),
        quality = quality,
        pairedEndSupport = 0,  // Depth-only call
        splitReadSupport = 0,
        relativeDepth = Some(math.pow(2, seg.log2Ratio)),
        mateChrom = None,
        matePos = None,
        filter = if (quality >= config.minQuality) "PASS" else "LowQual",
        genotype = genotype
      )
    }
  }
}
