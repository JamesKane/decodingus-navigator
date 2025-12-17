package com.decodingus.haplogroup.processor

import io.circe.Codec
import com.decodingus.analysis.{CallableLociQueryService, CallableState, VcfQueryService}

/**
 * Source of a SNP call for haplogroup analysis.
 */
enum CallSource derives Codec.AsObject:
  case FromVcf(depth: Int, quality: Double)
  case InferredReference(callableState: String)  // CALLABLE, REF_N
  case NoCall(reason: String)                    // LOW_COVERAGE, NO_COVERAGE, etc.

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
 * - If callable loci shows CALLABLE: infer homozygous reference
 * - If callable loci shows NO_COVERAGE/LOW_COVERAGE: mark as no-call
 * - If callable loci shows REF_N: the reference is N, cannot determine
 *
 * This allows haplogroup analysis to work even when the VCF has gaps
 * (positions where no variant was called because the sample matches reference).
 */
class GapAwareHaplogroupResolver(
  vcfService: VcfQueryService,
  callableLociService: Option[CallableLociQueryService],
  referenceBuild: String
) {

  /**
   * Resolve SNP calls at specified haplogroup tree positions.
   *
   * @param positions List of (contig, position, referenceAllele) for tree positions
   * @return Map of position to enhanced SNP call
   */
  def resolvePositions(
    positions: List[(String, Long, String)]
  ): Map[(String, Long), EnhancedSnpCall] = {
    // Query VCF for all positions
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

    positions.map { case (contig, position, refAllele) =>
      val call = vcfResults.get((contig, position)).flatten match {
        case Some(vcfCall) =>
          // Position found in VCF
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
          // Position not in VCF - check callable loci
          resolveFromCallableLoci(contig, position, refAllele)
      }
      (contig, position) -> call
    }.toMap
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
   * When VCF has no call, check callable loci to determine if we can infer reference.
   */
  private def resolveFromCallableLoci(
    contig: String,
    position: Long,
    refAllele: String
  ): EnhancedSnpCall = {
    callableLociService match {
      case None =>
        // No callable loci data - mark as no-call
        EnhancedSnpCall(
          contig = contig,
          position = position,
          referenceBuild = referenceBuild,
          allele = ".",
          source = CallSource.NoCall("No callable loci data available"),
          quality = None
        )

      case Some(service) =>
        val state = service.queryPosition(contig, position)
        state match {
          case CallableState.Callable =>
            // Position is callable but not in VCF = homozygous reference
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = refAllele,
              source = CallSource.InferredReference("CALLABLE"),
              quality = None
            )

          case CallableState.RefN =>
            // Reference is N - cannot determine allele
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = "N",
              source = CallSource.InferredReference("REF_N"),
              quality = None
            )

          case CallableState.NoCoverage =>
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = ".",
              source = CallSource.NoCall("NO_COVERAGE"),
              quality = None
            )

          case CallableState.LowCoverage =>
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = ".",
              source = CallSource.NoCall("LOW_COVERAGE"),
              quality = None
            )

          case CallableState.ExcessiveCoverage =>
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = ".",
              source = CallSource.NoCall("EXCESSIVE_COVERAGE"),
              quality = None
            )

          case CallableState.PoorMappingQuality =>
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = ".",
              source = CallSource.NoCall("POOR_MAPPING_QUALITY"),
              quality = None
            )

          case CallableState.Unknown =>
            EnhancedSnpCall(
              contig = contig,
              position = position,
              referenceBuild = referenceBuild,
              allele = ".",
              source = CallSource.NoCall("Position not in callable loci data"),
              quality = None
            )
        }
    }
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
   */
  def fromCache(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    referenceBuild: String
  ): Either[String, GapAwareHaplogroupResolver] = {
    VcfQueryService.fromCache(sampleAccession, runId, alignmentId).map { vcfService =>
      val callableService = CallableLociQueryService.fromAlignment(sampleAccession, runId, alignmentId)
      new GapAwareHaplogroupResolver(vcfService, callableService, referenceBuild)
    }
  }

  /**
   * Create a resolver from AT URIs.
   */
  def fromUris(
    sampleAccession: String,
    sequenceRunUri: Option[String],
    alignmentUri: Option[String],
    referenceBuild: String
  ): Either[String, GapAwareHaplogroupResolver] = {
    VcfQueryService.fromUris(sampleAccession, sequenceRunUri, alignmentUri).map { vcfService =>
      val callableService = CallableLociQueryService.fromUris(sampleAccession, sequenceRunUri, alignmentUri)
      new GapAwareHaplogroupResolver(vcfService, callableService, referenceBuild)
    }
  }

  /**
   * Create a resolver from a VCF file path directly.
   * Used for vendor-provided VCFs (e.g., FTDNA Big Y).
   *
   * Note: Vendor VCFs typically don't have callable loci data, so reference inference
   * will be limited. Positions not in the VCF will be marked as no-call.
   *
   * @param vcfPath Path to the VCF file (must be indexed)
   * @param referenceBuild Reference build of the VCF
   * @param targetBedPath Optional path to target regions BED (not callable loci - just capture targets)
   */
  def fromVcfPath(
    vcfPath: String,
    referenceBuild: String,
    targetBedPath: Option[String] = None
  ): Either[String, GapAwareHaplogroupResolver] = {
    VcfQueryService.fromVcfPath(vcfPath, referenceBuild).map { vcfService =>
      // Vendor VCFs don't have callable loci data - target BED is capture regions, not callable loci
      // We could potentially use target BED to infer that positions in targets but not in VCF are ref,
      // but that's not as reliable as true callable loci analysis
      new GapAwareHaplogroupResolver(vcfService, None, referenceBuild)
    }
  }
}
