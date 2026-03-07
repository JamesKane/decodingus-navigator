package com.decodingus.ibd.engine

import com.decodingus.workspace.model.IbdSegment
import scala.collection.mutable.ArrayBuffer

/**
 * Configuration for the pairwise IBD detector.
 *
 * @param minSegmentCm      Minimum segment length in cM to report
 * @param minSnpCount        Minimum SNP count in a segment to report
 * @param windowSize         Sliding window size (number of SNPs)
 * @param ibsThreshold       Minimum IBS-2 fraction in a window to call IBD
 * @param errorTolerance     Fraction of IBS-0 allowed within a segment (genotyping error tolerance)
 * @param minGapBp           Minimum gap in bp between segments to keep separate (else merge)
 */
case class IbdDetectorConfig(
                              minSegmentCm: Double = 7.0,
                              minSnpCount: Int = 100,
                              windowSize: Int = 100,
                              ibsThreshold: Double = 0.70,
                              errorTolerance: Double = 0.01,
                              minGapBp: Int = 1_000_000
                            )

/**
 * Genotype data for a single individual at a set of positions on one chromosome.
 *
 * @param chromosome Chromosome name
 * @param positions  Sorted physical positions in base pairs
 * @param genotypes  Genotype values: 0 (hom-ref), 1 (het), 2 (hom-alt), -1 (no-call)
 */
case class ChromosomeGenotypes(
                                chromosome: String,
                                positions: Array[Int],
                                genotypes: Array[Byte]
                              ):
  require(positions.length == genotypes.length,
    s"Positions (${positions.length}) and genotypes (${genotypes.length}) must be same length")

  def size: Int = positions.length

/**
 * Pairwise IBD segment detector using Identity-by-State (IBS) analysis.
 *
 * Algorithm:
 * 1. At each SNP position, classify allele sharing as IBS-0, IBS-1, or IBS-2
 * 2. Slide a window across positions, tracking IBS-2 fraction
 * 3. Mark runs of high-IBS sites as candidate IBD segments
 * 4. Apply error tolerance (allow sparse IBS-0 within segments)
 * 5. Convert physical positions to genetic distance (cM) via genetic map
 * 6. Filter segments below minimum cM threshold
 * 7. Merge nearby segments separated by small gaps
 */
class PairwiseIbdDetector(config: IbdDetectorConfig = IbdDetectorConfig()):

  /**
   * Detect IBD segments between two individuals across all shared chromosomes.
   *
   * @param sample1    Genotypes per chromosome for individual 1
   * @param sample2    Genotypes per chromosome for individual 2
   * @param geneticMap Genetic map for bp→cM conversion
   * @return List of detected IBD segments above threshold
   */
  def detectSegments(
                      sample1: Map[String, ChromosomeGenotypes],
                      sample2: Map[String, ChromosomeGenotypes],
                      geneticMap: GeneticMap
                    ): List[IbdSegment] =
    val sharedChromosomes = sample1.keySet.intersect(sample2.keySet)
    sharedChromosomes.toList.sorted.flatMap { chr =>
      detectChromosomeSegments(sample1(chr), sample2(chr), geneticMap)
    }

  /**
   * Detect IBD segments on a single chromosome.
   */
  def detectChromosomeSegments(
                                geno1: ChromosomeGenotypes,
                                geno2: ChromosomeGenotypes,
                                geneticMap: GeneticMap
                              ): List[IbdSegment] =
    // Intersect positions to find shared sites
    val (alignedPos, alignedG1, alignedG2) = intersectPositions(geno1, geno2)
    if alignedPos.length < config.minSnpCount then return Nil

    // Classify IBS at each position
    val ibsStates = classifyIbs(alignedG1, alignedG2)

    // Find candidate segments using sliding window
    val candidates = findCandidateSegments(alignedPos, ibsStates)

    // Merge nearby segments
    val merged = mergeSegments(candidates)

    // Convert to cM and filter
    merged.flatMap { case (start, end, snpCount, ibs2Count) =>
      val startBp = alignedPos(start)
      val endBp = alignedPos(end)
      geneticMap.intervalCm(geno1.chromosome, startBp, endBp).flatMap { cm =>
        if cm >= config.minSegmentCm && snpCount >= config.minSnpCount then
          Some(IbdSegment(
            chromosome = GeneticMap.normalizeChromosome(geno1.chromosome),
            startPosition = startBp,
            endPosition = endBp,
            lengthCm = math.round(cm * 100.0) / 100.0, // Round to 2 decimals
            snpCount = Some(snpCount),
            isHalfIdentical = Some(true) // IBS-based detection finds IBD1 (half-identical)
          ))
        else None
      }
    }

  /**
   * Intersect two genotype arrays to aligned positions.
   * Only includes positions present in both samples with valid genotypes.
   */
  private def intersectPositions(
                                  geno1: ChromosomeGenotypes,
                                  geno2: ChromosomeGenotypes
                                ): (Array[Int], Array[Byte], Array[Byte]) =
    val pos = ArrayBuffer.empty[Int]
    val g1 = ArrayBuffer.empty[Byte]
    val g2 = ArrayBuffer.empty[Byte]

    var i = 0
    var j = 0
    while i < geno1.size && j < geno2.size do
      if geno1.positions(i) == geno2.positions(j) then
        // Only include if both have valid calls
        if geno1.genotypes(i) >= 0 && geno2.genotypes(j) >= 0 then
          pos += geno1.positions(i)
          g1 += geno1.genotypes(i)
          g2 += geno2.genotypes(j)
        i += 1
        j += 1
      else if geno1.positions(i) < geno2.positions(j) then
        i += 1
      else
        j += 1

    (pos.toArray, g1.toArray, g2.toArray)

  /**
   * Classify IBS state at each position.
   * Returns: 0 = IBS-0 (no alleles shared), 1 = IBS-1, 2 = IBS-2 (both alleles shared)
   */
  private def classifyIbs(g1: Array[Byte], g2: Array[Byte]): Array[Byte] =
    val n = g1.length
    val ibs = new Array[Byte](n)
    var i = 0
    while i < n do
      ibs(i) = ibsState(g1(i), g2(i))
      i += 1
    ibs

  /**
   * Compute IBS state between two diploid genotypes.
   *
   * Genotypes: 0=AA, 1=AB, 2=BB
   * IBS-0: AA vs BB or BB vs AA (0 alleles shared)
   * IBS-1: AA vs AB, AB vs BB, etc. (1 allele shared)
   * IBS-2: AA vs AA, AB vs AB, BB vs BB (2 alleles shared)
   */
  private def ibsState(g1: Byte, g2: Byte): Byte =
    val diff = math.abs(g1 - g2)
    if diff == 0 then 2      // Same genotype → IBS-2
    else if diff == 1 then 1 // Adjacent → IBS-1
    else 0                   // Opposite homozygotes → IBS-0

  /**
   * Find candidate IBD segments using a sliding window approach.
   * Returns list of (startIdx, endIdx, snpCount, ibs2Count) tuples.
   */
  private def findCandidateSegments(
                                     positions: Array[Int],
                                     ibsStates: Array[Byte]
                                   ): List[(Int, Int, Int, Int)] =
    val n = positions.length
    if n < config.windowSize then return Nil

    val candidates = ArrayBuffer.empty[(Int, Int, Int, Int)]

    // Track IBS-0 count and IBS-2 count in current window
    var inSegment = false
    var segStart = 0
    var ibs0Count = 0
    var ibs2Count = 0
    var segSnpCount = 0

    val halfWindow = config.windowSize / 2

    var i = 0
    while i < n do
      // Check if current region has high IBS
      val lookBack = math.max(0, i - halfWindow)
      val lookForward = math.min(n - 1, i + halfWindow)

      // Simple IBS-2 + IBS-1 fraction in local window
      var localIbs2 = 0
      var localIbs0 = 0
      var localTotal = 0
      var j = lookBack
      while j <= lookForward do
        if ibsStates(j) == 2 then localIbs2 += 1
        else if ibsStates(j) == 0 then localIbs0 += 1
        localTotal += 1
        j += 1

      val ibsFraction = if localTotal > 0 then (localIbs2.toDouble + 0.5 * (localTotal - localIbs2 - localIbs0)) / localTotal else 0.0

      if ibsFraction >= config.ibsThreshold then
        if !inSegment then
          // Start new segment
          inSegment = true
          segStart = i
          ibs0Count = 0
          ibs2Count = 0
          segSnpCount = 0

        segSnpCount += 1
        if ibsStates(i) == 0 then ibs0Count += 1
        if ibsStates(i) == 2 then ibs2Count += 1
      else if inSegment then
        // Check if we've exceeded error tolerance
        val errorRate = if segSnpCount > 0 then ibs0Count.toDouble / segSnpCount else 1.0
        if errorRate <= config.errorTolerance || segSnpCount < config.windowSize then
          // Still within tolerance, extend
          segSnpCount += 1
          if ibsStates(i) == 0 then ibs0Count += 1
          if ibsStates(i) == 2 then ibs2Count += 1
        else
          // End segment
          if segSnpCount >= config.minSnpCount then
            candidates += ((segStart, i - 1, segSnpCount, ibs2Count))
          inSegment = false

      i += 1

    // Close any open segment
    if inSegment && segSnpCount >= config.minSnpCount then
      candidates += ((segStart, n - 1, segSnpCount, ibs2Count))

    candidates.toList

  /**
   * Merge nearby candidate segments separated by small gaps.
   */
  private def mergeSegments(
                             candidates: List[(Int, Int, Int, Int)]
                           ): List[(Int, Int, Int, Int)] =
    if candidates.size <= 1 then return candidates

    val merged = ArrayBuffer(candidates.head)
    for seg <- candidates.tail do
      val (_, prevEnd, prevSnps, prevIbs2) = merged.last
      val (curStart, curEnd, curSnps, curIbs2) = seg

      // Merge if gap between segments is small relative to segment size
      if curStart - prevEnd <= config.windowSize then
        val mergedSnps = prevSnps + curSnps + (curStart - prevEnd)
        merged(merged.length - 1) = (merged.last._1, curEnd, mergedSnps, prevIbs2 + curIbs2)
      else
        merged += seg

    merged.toList
