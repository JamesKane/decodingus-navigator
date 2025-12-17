package com.decodingus.haplogroup.processor

import com.decodingus.analysis.{CallableLociQueryService, CallableState, VcfQueryService}
import com.decodingus.refgenome.MultiContigReferenceQuerier
import io.circe.Codec

/**
 * Source of a SNP call for haplogroup analysis.
 */
enum CallSource derives Codec.AsObject:
  case FromVcf(depth: Int, quality: Double)
  case InferredReference(callableState: String) // CALLABLE, REF_N
  case NoCall(reason: String) // LOW_COVERAGE, NO_COVERAGE, etc.

/**
 * An enhanced SNP call with source information.
 */
case class EnhancedSnpCall(
                            contig: String,
                            position: Long,
                            referenceBuild: String,
                            allele: String,
                            source: CallSource,
                            quality: Option[Double]
                          ) derives Codec.AsObject {

  def isFromVcf: Boolean = source match {
    case CallSource.FromVcf(_, _) => true
    case _ => false
  }

  def isInferredReference: Boolean = source match {
    case CallSource.InferredReference(_) => true
    case _ => false
  }

  def isNoCall: Boolean = source match {
    case CallSource.NoCall(_) => true
    case _ => false
  }

  def hasCall: Boolean = !isNoCall
}

/**
 * Resolves SNP calls by combining VCF queries with callable loci data.
 *
 * When a position is not in the VCF:
 * - If callable loci shows CALLABLE: query reference genome for actual allele
 * - If callable loci shows NO_COVERAGE/LOW_COVERAGE: mark as no-call
 * - If callable loci shows REF_N: the reference is N, cannot determine
 *
 * This allows haplogroup analysis to work even when the VCF has gaps
 * (positions where no variant was called because the sample matches reference).
 *
 * IMPORTANT: The reference allele is queried from the actual reference genome,
 * NOT assumed to be the tree's ancestral allele. This is critical because the
 * reference genome may have derived alleles at haplogroup-defining positions.
 *
 * Performance: The resolver pre-loads VCF and callable loci data for the target
 * contigs into memory for fast O(1) lookups during position resolution.
 */
class GapAwareHaplogroupResolver(
                                  vcfService: VcfQueryService,
                                  callableLociService: Option[CallableLociQueryService],
                                  referenceQuerier: Option[MultiContigReferenceQuerier],
                                  referenceBuild: String
                                ) {

  /**
   * Pre-load data for specified contigs to optimize subsequent queries.
   * Call this before resolvePositions() for best performance.
   *
   * @param contigs List of contig names to preload (e.g., List("chrY") for Y-DNA)
   */
  def preloadContigs(contigs: List[String]): Unit = {
    // Pre-load callable loci data
    callableLociService.foreach(_.preloadContigs(contigs))
  }

  /**
   * Get the underlying VcfQueryService for additional queries.
   * Useful for extracting private variants after position resolution.
   */
  def getVcfService: VcfQueryService = vcfService

  /**
   * Close the underlying VCF reader, reference querier, and free resources.
   * Call this when done with the resolver.
   */
  def close(): Unit = {
    vcfService.close()
    referenceQuerier.foreach(_.close())
  }

  /**
   * Resolve SNP calls at specified haplogroup tree positions.
   *
   * Performance optimized: Batch-queries both VCF and callable loci to minimize
   * per-position overhead when resolving hundreds of thousands of positions.
   *
   * @param positions List of (contig, position, referenceAllele) for tree positions
   * @return Map of position to enhanced SNP call
   */
  def resolvePositions(
                        positions: List[(String, Long, String)]
                      ): Map[(String, Long), EnhancedSnpCall] = {
    val startTime = System.currentTimeMillis()

    // Query VCF for all positions (already optimized with in-memory cache)
    val vcfResults = vcfService.queryPositions(
      referenceBuild,
      positions.map { case (c, p, _) => (c, p) }
    ) match {
      case Left(mismatch) =>
        // Build mismatch - return empty results
        println(s"[GapAwareResolver] Build mismatch: expected ${mismatch.expected}, got ${mismatch.actual}")
        return Map.empty
      case Right(results) => results
    }

    // Identify positions NOT in VCF that need callable loci lookup
    val positionsNotInVcf = positions.filter { case (c, p, _) =>
      vcfResults.get((c, p)).flatten.isEmpty
    }

    // Batch-query callable loci for all missing positions at once
    val callableLociResults: Map[(String, Long), CallableState] = callableLociService match {
      case Some(service) if positionsNotInVcf.nonEmpty =>
        service.queryPositions(positionsNotInVcf.map { case (c, p, _) => (c, p) })
      case _ =>
        Map.empty
    }

    val vcfQueryTime = System.currentTimeMillis() - startTime

    // Build result map efficiently
    val resultBuilder = scala.collection.mutable.Map.empty[(String, Long), EnhancedSnpCall]

    positions.foreach { case (contig, position, _) =>
      val call = vcfResults.get((contig, position)).flatten match {
        case Some(vcfCall) =>
          // Position found in VCF - use the VCF allele directly
          val allele = if (vcfCall.isVariant) vcfCall.alt else vcfCall.ref
          EnhancedSnpCall(
            contig = contig,
            position = position,
            referenceBuild = referenceBuild,
            allele = allele,
            source = CallSource.FromVcf(vcfCall.depth.getOrElse(0), vcfCall.quality.getOrElse(0.0)),
            quality = vcfCall.quality
          )

        case None =>
          // Position not in VCF - check callable loci and query reference genome if callable
          val state = callableLociResults.getOrElse((contig, position), CallableState.Unknown)
          callableStateToEnhancedCall(contig, position, state)
      }
      resultBuilder((contig, position)) = call
    }

    val totalTime = System.currentTimeMillis() - startTime
    if (positions.size > 1000) {
      println(s"[GapAwareResolver] Resolved ${positions.size} positions in ${totalTime}ms (VCF: ${vcfQueryTime}ms)")
    }

    resultBuilder.toMap
  }

  /**
   * Convert a callable loci state to an EnhancedSnpCall.
   * For CALLABLE positions, queries the actual reference genome to get the allele.
   */
  private def callableStateToEnhancedCall(
                                           contig: String,
                                           position: Long,
                                           state: CallableState
                                         ): EnhancedSnpCall = {
    state match {
      case CallableState.Callable =>
        // Query the actual reference genome for the allele at this position
        val refAllele = referenceQuerier.flatMap(_.getBase(contig, position).map(_.toString.toUpperCase))
          .getOrElse(".")
        if (refAllele == ".") {
          // Couldn't get reference allele - treat as no-call
          EnhancedSnpCall(contig, position, referenceBuild, ".",
            CallSource.NoCall("Reference lookup failed"), None)
        } else {
          EnhancedSnpCall(contig, position, referenceBuild, refAllele,
            CallSource.InferredReference("CALLABLE"), None)
        }
      case CallableState.RefN =>
        EnhancedSnpCall(contig, position, referenceBuild, "N",
          CallSource.InferredReference("REF_N"), None)
      case CallableState.NoCoverage =>
        EnhancedSnpCall(contig, position, referenceBuild, ".",
          CallSource.NoCall("NO_COVERAGE"), None)
      case CallableState.LowCoverage =>
        EnhancedSnpCall(contig, position, referenceBuild, ".",
          CallSource.NoCall("LOW_COVERAGE"), None)
      case CallableState.ExcessiveCoverage =>
        EnhancedSnpCall(contig, position, referenceBuild, ".",
          CallSource.NoCall("EXCESSIVE_COVERAGE"), None)
      case CallableState.PoorMappingQuality =>
        EnhancedSnpCall(contig, position, referenceBuild, ".",
          CallSource.NoCall("POOR_MAPPING_QUALITY"), None)
      case CallableState.Unknown =>
        EnhancedSnpCall(contig, position, referenceBuild, ".",
          CallSource.NoCall("Position not in callable loci data"), None)
    }
  }

  /**
   * Resolve a single position.
   */
  def resolvePosition(contig: String, position: Long, refAllele: String): EnhancedSnpCall = {
    resolvePositions(List((contig, position, refAllele))).getOrElse(
      (contig, position),
      EnhancedSnpCall(
        contig = contig,
        position = position,
        referenceBuild = referenceBuild,
        allele = ".",
        source = CallSource.NoCall("Resolution failed"),
        quality = None
      )
    )
  }

  /**
   * Get statistics about the resolution results.
   */
  def getResolutionStats(results: Map[(String, Long), EnhancedSnpCall]): ResolutionStats = {
    val fromVcf = results.values.count(_.isFromVcf)
    val inferred = results.values.count(_.isInferredReference)
    val noCalls = results.values.count(_.isNoCall)

    ResolutionStats(
      totalPositions = results.size,
      fromVcf = fromVcf,
      inferredReference = inferred,
      noCalls = noCalls,
      callRate = if (results.nonEmpty) (fromVcf + inferred).toDouble / results.size else 0.0
    )
  }
}

case class ResolutionStats(
                            totalPositions: Int,
                            fromVcf: Int,
                            inferredReference: Int,
                            noCalls: Int,
                            callRate: Double
                          )

object GapAwareHaplogroupResolver {

  /**
   * Create a resolver from cached artifacts.
   *
   * @param sampleAccession Sample accession
   * @param runId           Sequence run ID
   * @param alignmentId     Alignment ID
   * @param referenceBuild  Reference build name
   * @param referencePath   Optional path to reference genome FASTA for inferring reference alleles at callable positions
   */
  def fromCache(
                 sampleAccession: String,
                 runId: String,
                 alignmentId: String,
                 referenceBuild: String,
                 referencePath: Option[String] = None
               ): Either[String, GapAwareHaplogroupResolver] = {
    VcfQueryService.fromCache(sampleAccession, runId, alignmentId).map { vcfService =>
      val callableService = CallableLociQueryService.fromAlignment(sampleAccession, runId, alignmentId)
      val refQuerier = referencePath.map(path => new MultiContigReferenceQuerier(path))
      new GapAwareHaplogroupResolver(vcfService, callableService, refQuerier, referenceBuild)
    }
  }

  /**
   * Create a resolver from AT URIs.
   *
   * @param sampleAccession Sample accession
   * @param sequenceRunUri  Optional sequence run URI
   * @param alignmentUri    Optional alignment URI
   * @param referenceBuild  Reference build name
   * @param referencePath   Optional path to reference genome FASTA for inferring reference alleles at callable positions
   */
  def fromUris(
                sampleAccession: String,
                sequenceRunUri: Option[String],
                alignmentUri: Option[String],
                referenceBuild: String,
                referencePath: Option[String] = None
              ): Either[String, GapAwareHaplogroupResolver] = {
    VcfQueryService.fromUris(sampleAccession, sequenceRunUri, alignmentUri).map { vcfService =>
      val callableService = CallableLociQueryService.fromUris(sampleAccession, sequenceRunUri, alignmentUri)
      val refQuerier = referencePath.map(path => new MultiContigReferenceQuerier(path))
      new GapAwareHaplogroupResolver(vcfService, callableService, refQuerier, referenceBuild)
    }
  }

  /**
   * Create a resolver from a VCF file path directly.
   * Used for vendor-provided VCFs (e.g., FTDNA Big Y).
   *
   * Note: Vendor VCFs typically don't have callable loci data, so reference inference
   * will be limited. Positions not in the VCF will be marked as no-call.
   *
   * @param vcfPath        Path to the VCF file (must be indexed)
   * @param referenceBuild Reference build of the VCF
   * @param referencePath  Optional path to reference genome FASTA for inferring reference alleles at callable positions
   * @param targetBedPath  Optional path to target regions BED (not callable loci - just capture targets)
   */
  def fromVcfPath(
                   vcfPath: String,
                   referenceBuild: String,
                   referencePath: Option[String] = None,
                   targetBedPath: Option[String] = None
                 ): Either[String, GapAwareHaplogroupResolver] = {
    VcfQueryService.fromVcfPath(vcfPath, referenceBuild).map { vcfService =>
      // Vendor VCFs don't have callable loci data - target BED is capture regions, not callable loci
      // We could potentially use target BED to infer that positions in targets but not in VCF are ref,
      // but that's not as reliable as true callable loci analysis
      val refQuerier = referencePath.map(path => new MultiContigReferenceQuerier(path))
      new GapAwareHaplogroupResolver(vcfService, None, refQuerier, referenceBuild)
    }
  }
}
