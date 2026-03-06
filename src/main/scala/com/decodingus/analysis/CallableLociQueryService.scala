package com.decodingus.analysis

import com.decodingus.util.Logger

import java.io.{BufferedReader, File, FileReader}
import java.nio.file.{Files, Path}
import scala.collection.mutable
import scala.jdk.CollectionConverters.*
import scala.util.Using

/**
 * Callable loci states from GATK CallableLoci.
 */
enum CallableState:
  case Callable // Position is callable (good coverage, quality)
  case NoCoverage // No reads covering this position
  case LowCoverage // Coverage below threshold
  case ExcessiveCoverage // Coverage too high (possible mapping issue)
  case PoorMappingQuality // Reads have poor mapping quality
  case RefN // Reference is N at this position
  case Unknown // Position not in BED file

  def isCallable: Boolean = this == Callable

  def canInferReference: Boolean = this == Callable || this == RefN

/**
 * Service for querying callable loci BED files at specific genomic positions.
 *
 * BED files are expected in the format:
 * {contig}\t{start}\t{end}\t{state}
 *
 * Where state is one of: CALLABLE, NO_COVERAGE, LOW_COVERAGE, EXCESSIVE_COVERAGE, POOR_MAPPING_QUALITY, REF_N
 *
 * Performance: The service loads contig data lazily and caches it in memory.
 * For batch operations, use preloadContigs() to load all needed contigs upfront.
 */
class CallableLociQueryService(callableLociDir: Path) {

  private val log = Logger[CallableLociQueryService]

  // Cache of loaded BED intervals per contig
  // Uses Array for O(1) random access in binary search (List would be O(n) per access!)
  private val contigIntervals: mutable.Map[String, Array[CallableInterval]] = mutable.Map.empty

  case class CallableInterval(start: Long, end: Long, state: CallableState)

  /**
   * Query the callable state at a specific position.
   *
   * @param contig   Chromosome name (e.g., "chr1")
   * @param position 1-based genomic position
   * @return The callable state at that position
   */
  def queryPosition(contig: String, position: Long): CallableState = {
    ensureContigLoaded(contig)
    contigIntervals.get(contig) match {
      case None => CallableState.Unknown
      case Some(intervals) =>
        // Use binary search for efficiency since intervals are sorted by start position
        findInterval(intervals, position) match {
          case Some(interval) => interval.state
          case None => CallableState.Unknown
        }
    }
  }

  /**
   * Performs a binary search to find the CallableInterval that contains the given position.
   * Assumes the array of intervals is sorted by start position.
   * Uses Array for O(1) random access - using List would make each access O(n)!
   */
  private def findInterval(intervals: Array[CallableInterval], position: Long): Option[CallableInterval] = {
    if (intervals.isEmpty) return None

    var low = 0
    var high = intervals.length - 1

    while (low <= high) {
      val mid = low + (high - low) / 2
      val interval = intervals(mid) // O(1) with Array

      if (position >= interval.start && position <= interval.end) {
        return Some(interval)
      } else if (position < interval.start) {
        high = mid - 1
      } else {
        low = mid + 1
      }
    }
    None
  }

  /**
   * Pre-load callable loci data for specified contigs.
   * Call this before batch queries for better performance.
   *
   * @param contigs List of contig names to preload
   */
  def preloadContigs(contigs: List[String]): Unit = {
    val startTime = System.currentTimeMillis()
    var totalIntervals = 0
    contigs.foreach { contig =>
      if (!contigIntervals.contains(contig)) {
        loadContigBed(contig)
        totalIntervals += contigIntervals.getOrElse(contig, Array.empty[CallableInterval]).length
      }
    }
    val elapsed = System.currentTimeMillis() - startTime
    if (contigs.nonEmpty) {
      log.info(s"Preloaded ${contigs.size} contigs with $totalIntervals intervals in ${elapsed}ms")
    }
  }

  /**
   * Query multiple positions at once.
   * For best performance, call preloadContigs() first with the relevant contigs.
   *
   * @param positions List of (contig, position) tuples
   * @return Map of position to callable state
   */
  def queryPositions(positions: List[(String, Long)]): Map[(String, Long), CallableState] = {
    val startTime = System.currentTimeMillis()

    // Group by contig to minimize BED file loading
    val byContig = positions.groupBy(_._1)

    // Ensure all contigs are loaded first (batch loading is more efficient)
    byContig.keys.foreach(ensureContigLoaded)

    val results = byContig.flatMap { case (contig, contigPositions) =>
      val intervals = contigIntervals.getOrElse(contig, Array.empty[CallableInterval])
      contigPositions.map { case (c, pos) =>
        val state = findInterval(intervals, pos).map(_.state).getOrElse(CallableState.Unknown)
        (c, pos) -> state
      }
    }

    val elapsed = System.currentTimeMillis() - startTime
    if (positions.size > 100) {
      log.info(s"Batch queried ${positions.size} callable loci positions in ${elapsed}ms")
    }

    results
  }

  /**
   * Check if a position is callable (can be used for variant calling).
   */
  def isCallable(contig: String, position: Long): Boolean = {
    queryPosition(contig, position).isCallable
  }

  /**
   * Check if reference can be inferred at a position.
   * True if CALLABLE or REF_N (we know what's there, even if it's N).
   */
  def canInferReference(contig: String, position: Long): Boolean = {
    queryPosition(contig, position).canInferReference
  }

  /**
   * Get callable base count for a contig.
   */
  def getCallableBasesForContig(contig: String): Long = {
    ensureContigLoaded(contig)
    contigIntervals.get(contig) match {
      case None => 0L
      case Some(intervals) =>
        intervals.filter(_.state == CallableState.Callable)
          .map(i => i.end - i.start + 1)
          .sum
    }
  }

  /**
   * Get total callable bases across all loaded contigs.
   */
  def getTotalCallableBases: Long = {
    // Load all contigs
    listAvailableContigs.foreach(ensureContigLoaded)
    contigIntervals.values.flatMap(_.iterator.filter(_.state == CallableState.Callable))
      .map(i => i.end - i.start + 1)
      .sum
  }

  /**
   * List available contigs (based on BED files present).
   */
  def listAvailableContigs: List[String] = {
    if (!Files.exists(callableLociDir)) {
      List.empty
    } else {
      Files.list(callableLociDir).iterator().asScala
        .filter(_.toString.endsWith(".callable.bed"))
        .map { path =>
          val name = path.getFileName.toString
          name.stripSuffix(".callable.bed")
        }
        .toList
    }
  }

  /**
   * Check if callable loci data exists for a contig.
   */
  def hasDataForContig(contig: String): Boolean = {
    val bedFile = callableLociDir.resolve(s"$contig.callable.bed")
    Files.exists(bedFile)
  }

  /**
   * Ensure a contig's BED file is loaded into memory.
   */
  private def ensureContigLoaded(contig: String): Unit = {
    if (!contigIntervals.contains(contig)) {
      loadContigBed(contig)
    }
  }

  /**
   * Load a contig's BED file.
   */
  private def loadContigBed(contig: String): Unit = {
    val bedFile = callableLociDir.resolve(s"$contig.callable.bed")

    if (!Files.exists(bedFile)) {
      contigIntervals(contig) = Array.empty
      return
    }

    val startTime = System.currentTimeMillis()
    val intervals = mutable.ArrayBuffer[CallableInterval]()

    Using(new BufferedReader(new FileReader(bedFile.toFile))) { reader =>
      var line = reader.readLine()
      while (line != null) {
        if (!line.startsWith("#") && line.trim.nonEmpty) {
          val fields = line.split("\\t")
          if (fields.length >= 4) {
            val intervalContig = fields(0)
            if (intervalContig == contig) {
              val start = fields(1).toLong + 1 // BED is 0-based, convert to 1-based
              val end = fields(2).toLong // BED end is exclusive, but we use inclusive
              val state = parseState(fields(3))
              intervals += CallableInterval(start, end, state)
            }
          }
        }
        line = reader.readLine()
      }
    }

    // Sort and convert to Array for O(1) random access in binary search
    val sortedIntervals = intervals.sortBy(_.start).toArray
    contigIntervals(contig) = sortedIntervals
    val elapsed = System.currentTimeMillis() - startTime
    log.info(s"Loaded callable loci for $contig: ${sortedIntervals.length} intervals in ${elapsed}ms")
  }

  /**
   * Parse a callable loci state string.
   */
  private def parseState(s: String): CallableState = s.toUpperCase match {
    case "CALLABLE" => CallableState.Callable
    case "NO_COVERAGE" => CallableState.NoCoverage
    case "LOW_COVERAGE" => CallableState.LowCoverage
    case "EXCESSIVE_COVERAGE" => CallableState.ExcessiveCoverage
    case "POOR_MAPPING_QUALITY" => CallableState.PoorMappingQuality
    case "REF_N" => CallableState.RefN
    case _ => CallableState.Unknown
  }

  /**
   * Clear cached intervals to free memory.
   */
  def clearCache(): Unit = {
    contigIntervals.clear()
  }
}

object CallableLociQueryService {

  /**
   * Create a query service from an alignment's artifact directory.
   */
  def fromAlignment(
                     sampleAccession: String,
                     runId: String,
                     alignmentId: String
                   ): Option[CallableLociQueryService] = {
    val callableLociDir = SubjectArtifactCache.getArtifactSubdir(
      sampleAccession, runId, alignmentId, "callable_loci"
    )
    if (Files.exists(callableLociDir) && Files.list(callableLociDir).findFirst().isPresent) {
      Some(new CallableLociQueryService(callableLociDir))
    } else {
      None
    }
  }

  /**
   * Create a query service from AT URIs.
   */
  def fromUris(
                sampleAccession: String,
                sequenceRunUri: Option[String],
                alignmentUri: Option[String]
              ): Option[CallableLociQueryService] = {
    val runId = sequenceRunUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
    val alignId = alignmentUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
    fromAlignment(sampleAccession, runId, alignId)
  }

  /**
   * Quick query for a single position.
   */
  def quickQuery(
                  sampleAccession: String,
                  runId: String,
                  alignmentId: String,
                  contig: String,
                  position: Long
                ): CallableState = {
    fromAlignment(sampleAccession, runId, alignmentId) match {
      case Some(service) => service.queryPosition(contig, position)
      case None => CallableState.Unknown
    }
  }
}
