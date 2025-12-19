package com.decodingus.analysis

import com.decodingus.model.ContigSummary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.{Interval, IntervalList, SamLocusIterator}
import htsjdk.samtools.{SamReaderFactory, ValidationStringency}

import java.io.{File, PrintWriter}
import java.nio.file.Path
import scala.collection.mutable
import scala.jdk.CollectionConverters.*
import scala.util.Using

// Note: CallableState enum is defined in CallableLociQueryService.scala
// and includes: Callable, NoCoverage, LowCoverage, ExcessiveCoverage, PoorMappingQuality, RefN, Unknown

/**
 * Parameters for callable loci determination.
 * Defaults match GATK CallableLoci.
 */
case class CallableLociParams(
  minDepth: Int = 4,                    // GATK default: 4
  maxDepth: Option[Int] = None,         // GATK default: unlimited
  minMappingQuality: Int = 10,          // GATK default: 10
  minBaseQuality: Int = 20,             // GATK default: 20
  maxLowMapQ: Int = 1,                  // MAPQ threshold for "low" reads
  maxFractionLowMapQ: Double = 0.1      // GATK default: 0.1 (10%)
)

/**
 * Per-chromosome coverage metrics for visualization.
 */
case class ContigCoverageMetrics(
  contig: String,
  length: Long,
  meanCoverage: Double,
  medianCoverage: Double,
  pct1x: Double,
  pct10x: Double,
  pct20x: Double,
  pct30x: Double,
  coverageHistogram: Array[Long]  // depth 0-255+ binned counts
)

/**
 * Combined result from coverage and callable analysis.
 * Replaces both WgsMetrics and CallableLociResult.
 */
case class CoverageCallableResult(
  // Global coverage metrics (replaces WgsMetrics)
  genomeTerritory: Long,
  meanCoverage: Double,
  medianCoverage: Double,
  sdCoverage: Double,
  coverageHistogram: Array[Long],  // depth 0-255+
  pct1x: Double,
  pct5x: Double,
  pct10x: Double,
  pct15x: Double,
  pct20x: Double,
  pct25x: Double,
  pct30x: Double,
  pct40x: Double,
  pct50x: Double,

  // Callable loci summary (replaces CallableLociResult)
  callableBases: Long,
  contigSummaries: List[ContigSummary],

  // Per-contig coverage for visualizations
  contigCoverage: Map[String, ContigCoverageMetrics]
)

/**
 * Accumulator for per-contig callable state counts.
 */
private class ContigCallableAccumulator {
  var refN: Long = 0L
  var noCoverage: Long = 0L
  var poorMappingQuality: Long = 0L
  var lowCoverage: Long = 0L
  var excessiveCoverage: Long = 0L
  var callable: Long = 0L

  def update(state: CallableState): Unit = state match {
    case CallableState.RefN => refN += 1
    case CallableState.NoCoverage => noCoverage += 1
    case CallableState.PoorMappingQuality => poorMappingQuality += 1
    case CallableState.LowCoverage => lowCoverage += 1
    case CallableState.ExcessiveCoverage => excessiveCoverage += 1
    case CallableState.Callable => callable += 1
    case CallableState.Unknown => () // Should not occur during analysis
  }

  def toContigSummary(contigName: String): ContigSummary = ContigSummary(
    contigName = contigName,
    refN = refN,
    callable = callable,
    noCoverage = noCoverage,
    lowCoverage = lowCoverage,
    excessiveCoverage = excessiveCoverage,
    poorMappingQuality = poorMappingQuality
  )
}

/**
 * Accumulator for per-contig coverage statistics using Welford's algorithm.
 */
private class ContigCoverageAccumulator(val contigName: String, val contigLength: Long) {
  private val histogram = new Array[Long](256)  // depth 0-255
  private var count = 0L
  private var mean = 0.0
  private var m2 = 0.0  // For Welford's variance

  def addDepth(depth: Int): Unit = {
    val clampedDepth = math.min(depth, 255)
    histogram(clampedDepth) += 1
    count += 1

    // Welford's online algorithm for mean and variance
    val delta = depth - mean
    mean += delta / count
    val delta2 = depth - mean
    m2 += delta * delta2
  }

  def getMean: Double = mean

  def getVariance: Double = if (count < 2) 0.0 else m2 / count

  def getStdDev: Double = math.sqrt(getVariance)

  def getMedian: Double = {
    if (count == 0) 0.0
    else {
      val halfCount = count / 2
      var cumulative = 0L
      var result = 255.0
      var found = false
      var depth = 0
      while (depth < histogram.length && !found) {
        cumulative += histogram(depth)
        if (cumulative >= halfCount) {
          result = depth.toDouble
          found = true
        }
        depth += 1
      }
      result
    }
  }

  def getHistogram: Array[Long] = histogram.clone()

  def getPctAtLeast(minDepth: Int): Double = {
    if (count == 0) return 0.0
    val atLeast = (minDepth until histogram.length).map(histogram(_)).sum
    atLeast.toDouble / count
  }

  def toContigCoverageMetrics: ContigCoverageMetrics = ContigCoverageMetrics(
    contig = contigName,
    length = contigLength,
    meanCoverage = getMean,
    medianCoverage = getMedian,
    pct1x = getPctAtLeast(1),
    pct10x = getPctAtLeast(10),
    pct20x = getPctAtLeast(20),
    pct30x = getPctAtLeast(30),
    coverageHistogram = getHistogram
  )
}

/**
 * Writer for BED file output with interval coalescing.
 */
private class BedFileWriter(outputPath: Path) extends AutoCloseable {
  private val writer = new PrintWriter(outputPath.toFile)
  private var currentContig: String = ""
  private var intervalStart: Long = 0
  private var currentState: CallableState = CallableState.NoCoverage

  def update(contig: String, position: Long, state: CallableState): Unit = {
    if (contig != currentContig) {
      // New contig - flush previous and reset
      if (currentContig.nonEmpty) {
        emitInterval()
      }
      currentContig = contig
      intervalStart = position
      currentState = state
    } else if (state != currentState) {
      // State change within contig
      emitInterval()
      intervalStart = position
      currentState = state
    }
    // Otherwise, same contig and state - extend current interval
  }

  def flushContig(contigLength: Long): Unit = {
    if (currentContig.nonEmpty) {
      // Emit final interval up to contig end
      writer.println(s"$currentContig\t${intervalStart - 1}\t$contigLength\t${stateToGatkName(currentState)}")
    }
  }

  /** Convert CallableState to GATK-style name for BED output */
  private def stateToGatkName(state: CallableState): String = state match {
    case CallableState.Callable => "CALLABLE"
    case CallableState.NoCoverage => "NO_COVERAGE"
    case CallableState.LowCoverage => "LOW_COVERAGE"
    case CallableState.ExcessiveCoverage => "EXCESSIVE_COVERAGE"
    case CallableState.PoorMappingQuality => "POOR_MAPPING_QUALITY"
    case CallableState.RefN => "REF_N"
    case CallableState.Unknown => "UNKNOWN"
  }

  private def emitInterval(): Unit = {
    // BED format is 0-based, half-open: [start, end)
    // Our positions are 1-based, so subtract 1 from start
    // The end is the position before the state change
    // This will be called when we hit position P with new state, so previous interval ended at P-1
    // Nothing to emit here - we emit on state change with the NEXT position as end
  }

  def emitIntervalTo(endPosition: Long): Unit = {
    if (currentContig.nonEmpty && endPosition > intervalStart) {
      writer.println(s"$currentContig\t${intervalStart - 1}\t${endPosition - 1}\t${stateToGatkName(currentState)}")
    }
  }

  override def close(): Unit = writer.close()
}

/**
 * Single-pass position-level walker that collects both coverage metrics
 * and callable loci state. Replaces both CollectWgsMetrics and CallableLoci.
 *
 * Uses HTSJDK SamLocusIterator for pileup without full GATK overhead.
 */
class CoverageCallableWalker {

  // Main assembly contigs only - excludes alts, decoys, HLA, etc.
  private val mainAssemblyPattern = "^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$".r

  private def isMainAssemblyContig(name: String): Boolean = {
    mainAssemblyPattern.findFirstIn(name).isDefined
  }

  /**
   * Process a BAM/CRAM file using SamLocusIterator for pileup.
   * Collects coverage histogram and callable state simultaneously.
   *
   * @param bamPath        Path to BAM/CRAM file
   * @param referencePath  Path to reference genome
   * @param outputDir      Directory for BED file output
   * @param callableParams Parameters for callable loci determination
   * @param onProgress     Progress callback (message, current, total)
   * @return Either error message or combined result
   */
  def collectCoverageAndCallable(
    bamPath: String,
    referencePath: String,
    outputDir: Path,
    callableParams: CallableLociParams = CallableLociParams(),
    onProgress: (String, Long, Long) => Unit
  ): Either[String, CoverageCallableResult] = {
    try {
      onProgress("Initializing coverage analysis...", 0, 1)

      val samReaderFactory = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)
        .referenceSequence(new File(referencePath))

      val samReader = samReaderFactory.open(new File(bamPath))
      val header = samReader.getFileHeader
      val sequenceDict = header.getSequenceDictionary

      // Build interval list for main assembly contigs only
      val intervalList = new IntervalList(header)
      val mainContigs = sequenceDict.getSequences.asScala
        .filter(seq => isMainAssemblyContig(seq.getSequenceName))
        .toList

      mainContigs.foreach { seq =>
        intervalList.add(new Interval(seq.getSequenceName, 1, seq.getSequenceLength))
      }

      // Calculate total bases for progress reporting
      val totalBases = mainContigs.map(_.getSequenceLength.toLong).sum
      val contigLengths = mainContigs.map(c => c.getSequenceName -> c.getSequenceLength.toLong).toMap

      onProgress(s"Analyzing ${mainContigs.size} contigs, $totalBases bases...", 0, totalBases)

      // Open reference for N-base detection
      val referenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath))

      // Reference base cache - load one contig at a time
      var currentRefContig: String = ""
      var currentRefBases: Array[Byte] = Array.empty

      def getRefBase(contig: String, position: Long): Byte = {
        if (contig != currentRefContig) {
          val seq = referenceFile.getSequence(contig)
          currentRefBases = seq.getBases
          currentRefContig = contig
        }
        currentRefBases((position - 1).toInt)
      }

      // Accumulators
      val globalCoverageHist = new Array[Long](256)
      var globalCount = 0L
      var globalMean = 0.0
      var globalM2 = 0.0

      val contigCallable = mutable.Map[String, ContigCallableAccumulator]()
      val contigCoverage = mutable.Map[String, ContigCoverageAccumulator]()
      val bedWriters = mutable.Map[String, BedFileWriter]()

      // Initialize per-contig accumulators
      mainContigs.foreach { seq =>
        val name = seq.getSequenceName
        contigCallable(name) = new ContigCallableAccumulator()
        contigCoverage(name) = new ContigCoverageAccumulator(name, seq.getSequenceLength)
      }

      // Create output directory
      outputDir.toFile.mkdirs()

      // Configure locus iterator
      val locusIterator = new SamLocusIterator(samReader, intervalList, true)
      locusIterator.setEmitUncoveredLoci(true)
      locusIterator.setIncludeIndels(false)
      locusIterator.setMappingQualityScoreCutoff(0)  // We handle MAPQ filtering ourselves

      var processedBases = 0L
      var lastProgressUpdate = 0L
      val progressInterval = 1000000L  // Report every 1M bases
      var lastContig = ""
      var lastPosition = 0L

      val iter = locusIterator.iterator()
      while (iter.hasNext) {
        val locus = iter.next()
        val contig = locus.getSequenceName
        val position = locus.getPosition.toLong

        // Handle contig transitions for BED output
        if (contig != lastContig) {
          // Flush previous contig's BED writer
          bedWriters.get(lastContig).foreach { writer =>
            writer.flushContig(contigLengths.getOrElse(lastContig, lastPosition))
            writer.close()
          }
          // Create new BED writer for this contig
          val bedPath = outputDir.resolve(s"$contig.callable.bed")
          bedWriters(contig) = new BedFileWriter(bedPath)
          lastContig = contig
        }

        lastPosition = position
        processedBases += 1

        // Progress reporting
        if (processedBases - lastProgressUpdate >= progressInterval) {
          onProgress(s"Processing $contig:$position...", processedBases, totalBases)
          lastProgressUpdate = processedBases
        }

        val pileup = locus.getRecordAndOffsets
        val rawDepth = pileup.size()

        // Update coverage statistics
        val clampedDepth = math.min(rawDepth, 255)
        globalCoverageHist(clampedDepth) += 1
        globalCount += 1

        // Welford's algorithm for global stats
        val delta = rawDepth - globalMean
        globalMean += delta / globalCount
        val delta2 = rawDepth - globalMean
        globalM2 += delta * delta2

        // Per-contig coverage
        contigCoverage(contig).addDepth(rawDepth)

        // Determine callable state
        val refBase = getRefBase(contig, position)
        val state = determineCallableState(refBase, pileup, callableParams)

        // Update callable counts
        contigCallable(contig).update(state)

        // Update BED writer
        bedWriters(contig).update(contig, position, state)
      }

      // Flush final contig's BED writer
      bedWriters.get(lastContig).foreach { writer =>
        writer.flushContig(contigLengths.getOrElse(lastContig, lastPosition))
        writer.close()
      }

      locusIterator.close()
      samReader.close()
      referenceFile.close()

      onProgress("Calculating final statistics...", totalBases, totalBases)

      // Calculate global percentiles from histogram
      def getPctAtLeast(hist: Array[Long], total: Long, minDepth: Int): Double = {
        if (total == 0) return 0.0
        val atLeast = (minDepth until hist.length).map(hist(_)).sum
        atLeast.toDouble / total
      }

      def getMedian(hist: Array[Long], total: Long): Double = {
        if (total == 0) 0.0
        else {
          val halfCount = total / 2
          var cumulative = 0L
          var result = 255.0
          var found = false
          var depth = 0
          while (depth < hist.length && !found) {
            cumulative += hist(depth)
            if (cumulative >= halfCount) {
              result = depth.toDouble
              found = true
            }
            depth += 1
          }
          result
        }
      }

      val globalStdDev = if (globalCount < 2) 0.0 else math.sqrt(globalM2 / globalCount)

      // Build contig summaries
      val contigSummaries = mainContigs.map { seq =>
        contigCallable(seq.getSequenceName).toContigSummary(seq.getSequenceName)
      }

      // Write per-contig summary files (GATK format)
      contigSummaries.foreach { summary =>
        val summaryPath = outputDir.resolve(s"${summary.contigName}.table.txt")
        Using.resource(new PrintWriter(summaryPath.toFile)) { writer =>
          writer.println("state nBases")
          writer.println(s"REF_N ${summary.refN}")
          writer.println(s"CALLABLE ${summary.callable}")
          writer.println(s"NO_COVERAGE ${summary.noCoverage}")
          writer.println(s"LOW_COVERAGE ${summary.lowCoverage}")
          writer.println(s"EXCESSIVE_COVERAGE ${summary.excessiveCoverage}")
          writer.println(s"POOR_MAPPING_QUALITY ${summary.poorMappingQuality}")
        }
      }

      val result = CoverageCallableResult(
        genomeTerritory = globalCount,
        meanCoverage = globalMean,
        medianCoverage = getMedian(globalCoverageHist, globalCount),
        sdCoverage = globalStdDev,
        coverageHistogram = globalCoverageHist,
        pct1x = getPctAtLeast(globalCoverageHist, globalCount, 1),
        pct5x = getPctAtLeast(globalCoverageHist, globalCount, 5),
        pct10x = getPctAtLeast(globalCoverageHist, globalCount, 10),
        pct15x = getPctAtLeast(globalCoverageHist, globalCount, 15),
        pct20x = getPctAtLeast(globalCoverageHist, globalCount, 20),
        pct25x = getPctAtLeast(globalCoverageHist, globalCount, 25),
        pct30x = getPctAtLeast(globalCoverageHist, globalCount, 30),
        pct40x = getPctAtLeast(globalCoverageHist, globalCount, 40),
        pct50x = getPctAtLeast(globalCoverageHist, globalCount, 50),
        callableBases = contigSummaries.map(_.callable).sum,
        contigSummaries = contigSummaries,
        contigCoverage = contigCoverage.view.mapValues(_.toContigCoverageMetrics).toMap
      )

      onProgress("Coverage and callable analysis complete.", totalBases, totalBases)
      Right(result)

    } catch {
      case e: Exception =>
        Left(s"Failed to collect coverage/callable metrics: ${e.getMessage}")
    }
  }

  /**
   * Determine callable state for a locus following GATK CallableLoci logic.
   * Evaluation is hierarchical - first failing condition wins.
   */
  private def determineCallableState(
    refBase: Byte,
    pileup: java.util.List[SamLocusIterator.RecordAndOffset],
    params: CallableLociParams
  ): CallableState = {
    // 1. Check reference - N bases are never callable
    if (refBase == 'N' || refBase == 'n') {
      return CallableState.RefN
    }

    val rawDepth = pileup.size()

    // 2. No coverage
    if (rawDepth == 0) {
      return CallableState.NoCoverage
    }

    // 3. Count QC-passing and low-MAPQ reads
    var qcPassCount = 0
    var lowMapQCount = 0

    val iter = pileup.iterator()
    while (iter.hasNext) {
      val rec = iter.next()
      val read = rec.getRecord
      val mapQ = read.getMappingQuality
      val baseQ = rec.getBaseQuality.toInt

      if (mapQ >= params.minMappingQuality && baseQ >= params.minBaseQuality) {
        qcPassCount += 1
      }
      if (mapQ <= params.maxLowMapQ) {
        lowMapQCount += 1
      }
    }

    // 4. Poor mapping quality (too many low-MAPQ reads)
    val lowMapQFraction = lowMapQCount.toDouble / rawDepth
    if (lowMapQFraction > params.maxFractionLowMapQ) {
      return CallableState.PoorMappingQuality
    }

    // 5. Low coverage (after QC filtering)
    if (qcPassCount < params.minDepth) {
      return CallableState.LowCoverage
    }

    // 6. Excessive coverage
    if (params.maxDepth.exists(qcPassCount > _)) {
      return CallableState.ExcessiveCoverage
    }

    // 7. Callable
    CallableState.Callable
  }
}
