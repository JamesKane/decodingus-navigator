package com.decodingus.analysis.sv

import scala.collection.mutable

/**
 * Clusters discordant read pairs and split reads into structural variant calls.
 *
 * Algorithm:
 * 1. Group evidence by approximate breakpoint location
 * 2. Cluster nearby evidence points
 * 3. Infer SV type from evidence patterns
 * 4. Merge with depth-based CNV calls where applicable
 *
 * References:
 * - Evidence clustering approach: Layer et al. "LUMPY: a probabilistic framework for
 *   structural variant discovery." Genome Biology 15.6 (2014): R84.
 *   https://doi.org/10.1186/gb-2014-15-6-r84
 *
 * - SV type inference from pair orientation: Rausch et al. "DELLY: structural variant
 *   discovery by integrated paired-end and split-read analysis."
 *   Bioinformatics 28.18 (2012): i333-i339.
 *   https://doi.org/10.1093/bioinformatics/bts378
 *
 * - Multi-evidence integration: Chiang et al. "SpeedSeq: ultra-fast personal genome
 *   analysis and interpretation." Nature Methods 12.10 (2015): 966-968.
 *   https://doi.org/10.1038/nmeth.3505
 */
class SvEvidenceClusterer(config: SvCallerConfig = SvCallerConfig.default) {

  private val CLUSTER_DISTANCE = config.maxClusterDistance

  /**
   * Cluster SV evidence into structural variant calls.
   *
   * @param evidence        Collected SV evidence
   * @param depthSegments   Depth-based CNV segments (for integration)
   * @return List of SV calls
   */
  def cluster(
    evidence: SvEvidenceCollection,
    depthSegments: List[DepthSegment] = Nil
  ): List[SvCall] = {

    val calls = mutable.ListBuffer[SvCall]()
    var callIndex = 0

    // 1. Handle inter-chromosomal events (translocations)
    val translocations = clusterTranslocations(evidence.interChromosomalPairs)
    calls ++= translocations.map { cluster =>
      callIndex += 1
      breakpointClusterToSvCall(cluster, SvType.BND, callIndex)
    }

    // 2. Handle intra-chromosomal events by grouping PE and SR evidence
    val intraPairs = evidence.discordantPairs.filterNot(_.reason == DiscordantReason.InterChromosomal)

    // Group by chromosome
    val pairsByChrom = intraPairs.groupBy(_.chrom1)
    val splitsByChrom = evidence.splitReads.groupBy(_.primaryChrom)

    val allChroms = (pairsByChrom.keys ++ splitsByChrom.keys).toSet

    allChroms.foreach { chrom =>
      val pairs = pairsByChrom.getOrElse(chrom, Nil)
      val splits = splitsByChrom.getOrElse(chrom, Nil)

      // Find clusters of evidence
      val clusters = clusterIntraChromosomalEvidence(chrom, pairs, splits)

      clusters.foreach { cluster =>
        // Infer SV type from evidence pattern
        val svType = inferSvType(cluster)

        if (cluster.totalSupport >= config.minTotalSupport) {
          callIndex += 1
          calls += breakpointClusterToSvCall(cluster, svType, callIndex)
        }
      }
    }

    // 3. Integrate with depth segments
    // Mark depth-only calls that also have PE/SR support
    val integratedCalls = integratePeSrWithDepth(calls.toList, depthSegments)

    // Sort by position
    integratedCalls.sortBy(c => (c.chrom, c.start))
  }

  /**
   * Cluster inter-chromosomal discordant pairs into translocation breakpoints.
   */
  private def clusterTranslocations(pairs: List[DiscordantPair]): List[BreakpointCluster] = {
    if (pairs.isEmpty) return Nil

    // Group by chromosome pair
    val byChromPair = pairs.groupBy(p => (p.chrom1, p.chrom2))

    byChromPair.flatMap { case ((chrom1, chrom2), pairsInGroup) =>
      // Cluster by position on first chromosome
      clusterByPosition(pairsInGroup.map(p => (p.pos1, p))).map { (pos, clusteredPairs) =>
        // Find corresponding position on second chromosome
        val matePositions = clusteredPairs.map(_._2.pos2)
        val matePos = matePositions.sum / matePositions.size

        BreakpointCluster(
          chrom = chrom1,
          position = pos,
          ciLow = -CLUSTER_DISTANCE,
          ciHigh = CLUSTER_DISTANCE,
          discordantPairs = clusteredPairs.map(_._2),
          splitReads = Nil,
          mateChrom = Some(chrom2),
          matePosition = Some(matePos)
        )
      }
    }.toList
  }

  /**
   * Cluster intra-chromosomal evidence into breakpoint clusters.
   */
  private def clusterIntraChromosomalEvidence(
    chrom: String,
    pairs: List[DiscordantPair],
    splits: List[SplitRead]
  ): List[BreakpointCluster] = {

    // Create evidence points with positions
    val pairPoints = pairs.map(p => (p.pos1, Left(p): Either[DiscordantPair, SplitRead]))
    val splitPoints = splits.map(s => (s.primaryPos, Right(s): Either[DiscordantPair, SplitRead]))
    val allPoints = (pairPoints ++ splitPoints).sortBy(_._1)

    if (allPoints.isEmpty) return Nil

    // Cluster nearby evidence
    val clusters = mutable.ListBuffer[BreakpointCluster]()
    var currentCluster = mutable.ListBuffer[(Long, Either[DiscordantPair, SplitRead])]()
    var clusterStart = allPoints.head._1

    allPoints.foreach { point =>
      if (point._1 - clusterStart <= CLUSTER_DISTANCE || currentCluster.isEmpty) {
        currentCluster += point
      } else {
        // Finalize current cluster and start new one
        if (currentCluster.nonEmpty) {
          clusters += createBreakpointCluster(chrom, currentCluster.toList)
        }
        currentCluster.clear()
        currentCluster += point
        clusterStart = point._1
      }
    }

    // Don't forget the last cluster
    if (currentCluster.nonEmpty) {
      clusters += createBreakpointCluster(chrom, currentCluster.toList)
    }

    clusters.toList
  }

  /**
   * Create a breakpoint cluster from evidence points.
   */
  private def createBreakpointCluster(
    chrom: String,
    points: List[(Long, Either[DiscordantPair, SplitRead])]
  ): BreakpointCluster = {
    val positions = points.map(_._1)
    val meanPos = positions.sum / positions.size
    val minPos = positions.min
    val maxPos = positions.max

    val pairs = points.flatMap(_._2.left.toOption)
    val splits = points.flatMap(_._2.toOption)

    BreakpointCluster(
      chrom = chrom,
      position = meanPos,
      ciLow = (minPos - meanPos).toInt,
      ciHigh = (maxPos - meanPos).toInt,
      discordantPairs = pairs,
      splitReads = splits,
      mateChrom = None,
      matePosition = None
    )
  }

  /**
   * Cluster positions that are within CLUSTER_DISTANCE of each other.
   */
  private def clusterByPosition[T](items: List[(Long, T)]): List[(Long, List[(Long, T)])] = {
    if (items.isEmpty) return Nil

    val sorted = items.sortBy(_._1)
    val clusters = mutable.ListBuffer[(Long, List[(Long, T)])]()
    var current = mutable.ListBuffer[(Long, T)]()
    var clusterStart = sorted.head._1

    sorted.foreach { item =>
      if (item._1 - clusterStart <= CLUSTER_DISTANCE || current.isEmpty) {
        current += item
      } else {
        val meanPos = current.map(_._1).sum / current.size
        clusters += ((meanPos, current.toList))
        current.clear()
        current += item
        clusterStart = item._1
      }
    }

    if (current.nonEmpty) {
      val meanPos = current.map(_._1).sum / current.size
      clusters += ((meanPos, current.toList))
    }

    clusters.toList
  }

  /**
   * Infer SV type from evidence pattern.
   */
  private def inferSvType(cluster: BreakpointCluster): SvType = {
    // If we have inter-chromosomal evidence, it's a translocation
    if (cluster.mateChrom.isDefined) return SvType.BND

    // Analyze pair orientations
    val orientations = cluster.discordantPairs.map { p =>
      (p.strand1, p.strand2, p.pos1 < p.pos2)
    }

    // Count orientation patterns
    val frCount = orientations.count { case (s1, s2, upstream) =>
      (upstream && s1 == '+' && s2 == '-') || (!upstream && s1 == '-' && s2 == '+')
    }
    val rfCount = orientations.count { case (s1, s2, upstream) =>
      (upstream && s1 == '-' && s2 == '+') || (!upstream && s1 == '+' && s2 == '-')
    }
    val sameStrandCount = orientations.count { case (s1, s2, _) => s1 == s2 }

    // Check for inversions (same strand orientation)
    if (sameStrandCount > orientations.size / 2) {
      return SvType.INV
    }

    // Check for deletions vs duplications based on insert size
    val insertSizeOutliers = cluster.discordantPairs.filter(_.reason == DiscordantReason.InsertSizeOutlier)
    if (insertSizeOutliers.nonEmpty) {
      // Large insert sizes typically indicate deletions
      // Small insert sizes typically indicate tandem duplications
      val avgInsert = insertSizeOutliers.map(_.insertSize).sum.toDouble / insertSizeOutliers.size
      val expectedInsert = cluster.discordantPairs.headOption
        .map(_.insertSize.toDouble * 0.5) // Rough estimate
        .getOrElse(400.0)

      if (avgInsert > expectedInsert * 2) {
        return SvType.DEL
      } else if (avgInsert < expectedInsert * 0.5) {
        return SvType.DUP
      }
    }

    // Default to deletion for large insert size events
    SvType.DEL
  }

  /**
   * Convert a breakpoint cluster to an SV call.
   */
  private def breakpointClusterToSvCall(
    cluster: BreakpointCluster,
    svType: SvType,
    index: Int
  ): SvCall = {
    // Calculate quality based on support
    val quality = math.min(cluster.totalSupport * 5.0 + cluster.meanMapQ * 0.5, 99.0)

    // Estimate SV length and end position
    val (svLen, end, mateChrom, matePos) = svType match {
      case SvType.BND =>
        (0L, cluster.position, cluster.mateChrom, cluster.matePosition)

      case _ =>
        // For other types, estimate from mate positions in discordant pairs
        val matePositions = cluster.discordantPairs.map(_.pos2).filter(_ != cluster.position)
        if (matePositions.nonEmpty) {
          val avgMatePos = matePositions.sum / matePositions.size
          val length = math.abs(avgMatePos - cluster.position)
          val endPos = cluster.position + length
          val signedLen = if (svType == SvType.DEL) -length else length
          (signedLen, endPos, None, None)
        } else {
          // Fallback: use cluster distance as estimate
          (CLUSTER_DISTANCE.toLong, cluster.position + CLUSTER_DISTANCE, None, None)
        }
    }

    // Determine genotype based on support level
    val genotype = if (cluster.totalSupport >= 10) "1/1" else "0/1"

    // Filter based on minimum support
    val filter = if (cluster.peSupport >= config.minPairedEndSupport ||
                     cluster.srSupport >= config.minSplitReadSupport) {
      "PASS"
    } else {
      "LowSupport"
    }

    SvCall(
      id = s"${svType}_${cluster.chrom}_${cluster.position}_$index",
      chrom = cluster.chrom,
      start = cluster.position,
      end = end,
      svType = svType,
      svLen = svLen,
      ciPos = (cluster.ciLow, cluster.ciHigh),
      ciEnd = (cluster.ciLow, cluster.ciHigh),
      quality = quality,
      pairedEndSupport = cluster.peSupport,
      splitReadSupport = cluster.srSupport,
      relativeDepth = None,
      mateChrom = mateChrom,
      matePos = matePos,
      filter = filter,
      genotype = genotype
    )
  }

  /**
   * Integrate PE/SR evidence with depth-based CNV calls.
   *
   * If a PE/SR call overlaps a depth segment, add the depth evidence.
   * If a depth segment has no PE/SR support, keep it as depth-only.
   */
  private def integratePeSrWithDepth(
    peSrCalls: List[SvCall],
    depthSegments: List[DepthSegment]
  ): List[SvCall] = {

    if (depthSegments.isEmpty) return peSrCalls

    val segmenter = new DepthSegmenter(config)
    val depthCalls = segmenter.toSvCalls(depthSegments)

    // Find which depth calls overlap with PE/SR calls
    val usedDepthIndices = mutable.Set[Int]()

    val enhancedPeSrCalls = peSrCalls.map { call =>
      // Find overlapping depth segment
      depthCalls.zipWithIndex.find { case (depthCall, idx) =>
        depthCall.chrom == call.chrom &&
        depthCall.svType == call.svType &&
        overlaps(call.start, call.end, depthCall.start, depthCall.end) &&
        !usedDepthIndices.contains(idx)
      } match {
        case Some((depthCall, idx)) =>
          usedDepthIndices += idx
          // Merge depth evidence into PE/SR call
          call.copy(relativeDepth = depthCall.relativeDepth)
        case None =>
          call
      }
    }

    // Add depth-only calls that don't overlap with PE/SR calls
    val depthOnlyCalls = depthCalls.zipWithIndex
      .filterNot { case (_, idx) => usedDepthIndices.contains(idx) }
      .map(_._1)

    enhancedPeSrCalls ++ depthOnlyCalls
  }

  /**
   * Check if two intervals overlap.
   */
  private def overlaps(start1: Long, end1: Long, start2: Long, end2: Long): Boolean = {
    start1 <= end2 && start2 <= end1
  }
}
