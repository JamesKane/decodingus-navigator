package com.decodingus.analysis

import com.decodingus.util.Logger
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import java.io.File
import java.nio.file.Path
import scala.collection.mutable
import scala.jdk.CollectionConverters.*
import scala.util.{Either, Left, Right, Try, Using}

/**
 * Error when VCF reference build doesn't match expected build.
 */
case class BuildMismatch(expected: String, actual: String)

/**
 * A variant call from the VCF.
 */
case class VariantCall(
                        contig: String,
                        position: Long,
                        ref: String,
                        alt: String,
                        genotype: String,
                        depth: Option[Int],
                        quality: Option[Double],
                        filter: String
                      ) {
  def isVariant: Boolean = alt != ref && alt != "."

  def isHomRef: Boolean = genotype == "0/0" || genotype == "0"

  def isHomAlt: Boolean = genotype == "1/1" || genotype == "1"

  def isHet: Boolean = genotype == "0/1" || genotype == "1/0"
}

/**
 * Service for querying a cached VCF at specific genomic positions.
 * Uses tabix index for efficient random access.
 *
 * Build-aware: validates that query build matches VCF metadata.
 *
 * Performance optimization: When querying many positions on the same contig,
 * the service loads all variants for that contig into memory and performs
 * lookups from the in-memory map. This reduces tabix queries from O(n) to O(1).
 */
class VcfQueryService(vcfInfo: CachedVcfInfo) {

  private val log = Logger[VcfQueryService]

  private lazy val reader: VCFFileReader = {
    val vcfFile = new File(vcfInfo.vcfPath)
    new VCFFileReader(vcfFile, true) // true = require index
  }

  // In-memory cache of variants by contig -> position -> VariantCall
  // This dramatically speeds up batch queries on the same contig
  private val contigCache: mutable.Map[String, Map[Long, VariantCall]] = mutable.Map.empty

  /**
   * Load all variants for a contig into memory for fast lookup.
   * This is more efficient than individual tabix queries when querying many positions.
   */
  private def ensureContigLoaded(contig: String): Map[Long, VariantCall] = {
    contigCache.getOrElseUpdate(contig, {
      val startTime = System.currentTimeMillis()
      val variants = Try {
        reader.query(contig, 1, Int.MaxValue).asScala.map { vc =>
          val call = variantContextToCall(vc)
          call.position -> call
        }.toMap
      }.getOrElse(Map.empty)
      val elapsed = System.currentTimeMillis() - startTime
      log.info(s"Loaded ${variants.size} variants for $contig in ${elapsed}ms")
      variants
    })
  }

  /**
   * Clear the contig cache to free memory.
   */
  def clearCache(): Unit = {
    contigCache.clear()
  }

  /**
   * Query a single position in the VCF.
   *
   * @param build    Expected reference build
   * @param contig   Chromosome name (e.g., "chr1")
   * @param position 1-based genomic position
   * @return Either build mismatch error or optional variant call
   */
  def queryPosition(
                     build: String,
                     contig: String,
                     position: Long
                   ): Either[BuildMismatch, Option[VariantCall]] = {
    validateBuild(build).map { _ =>
      Try {
        val results = reader.query(contig, position.toInt, position.toInt)
        if (results.hasNext) {
          Some(variantContextToCall(results.next()))
        } else {
          None
        }
      }.getOrElse(None)
    }
  }

  /**
   * Query multiple positions in the VCF (batch query).
   *
   * Performance optimized: Groups positions by contig and loads each contig's
   * variants into memory once, then performs O(1) lookups. This is much faster
   * than individual tabix queries when querying many positions.
   *
   * @param build     Expected reference build
   * @param positions List of (contig, position) tuples
   * @return Either build mismatch error or map of position to optional variant call
   */
  def queryPositions(
                      build: String,
                      positions: List[(String, Long)]
                    ): Either[BuildMismatch, Map[(String, Long), Option[VariantCall]]] = {
    validateBuild(build).map { _ =>
      val startTime = System.currentTimeMillis()

      // Group positions by contig for efficient batch loading
      val byContig = positions.groupBy(_._1)

      // Load each contig's variants into memory and perform lookups
      val results = byContig.flatMap { case (contig, contigPositions) =>
        val contigVariants = ensureContigLoaded(contig)
        contigPositions.map { case (c, pos) =>
          (c, pos) -> contigVariants.get(pos)
        }
      }

      val elapsed = System.currentTimeMillis() - startTime
      if (positions.size > 100) {
        log.info(s"Batch queried ${positions.size} positions across ${byContig.size} contigs in ${elapsed}ms")
      }

      results
    }
  }

  /**
   * Query all variants in a region.
   *
   * @param build  Expected reference build
   * @param contig Chromosome name
   * @param start  Start position (1-based, inclusive)
   * @param end    End position (1-based, inclusive)
   * @return Either build mismatch error or list of variant calls
   */
  def queryRegion(
                   build: String,
                   contig: String,
                   start: Long,
                   end: Long
                 ): Either[BuildMismatch, List[VariantCall]] = {
    validateBuild(build).map { _ =>
      Try {
        reader.query(contig, start.toInt, end.toInt).asScala.map(variantContextToCall).toList
      }.getOrElse(List.empty)
    }
  }

  /**
   * Get all variants on a contig.
   * Uses the in-memory cache if available for consistent performance.
   */
  def queryContig(build: String, contig: String): Either[BuildMismatch, Iterator[VariantCall]] = {
    validateBuild(build).map { _ =>
      // Use cache if already loaded, otherwise load from VCF
      val contigVariants = ensureContigLoaded(contig)
      contigVariants.values.iterator
    }
  }

  /**
   * Check if a position has a variant call (either ref or alt).
   */
  def hasCallAt(build: String, contig: String, position: Long): Either[BuildMismatch, Boolean] = {
    queryPosition(build, contig, position).map(_.isDefined)
  }

  /**
   * Get the allele at a position (alt if variant, ref if not).
   */
  def getAlleleAt(build: String, contig: String, position: Long): Either[BuildMismatch, Option[String]] = {
    queryPosition(build, contig, position).map { opt =>
      opt.map { call =>
        if (call.isVariant) call.alt else call.ref
      }
    }
  }

  /**
   * Validate that the query build matches the VCF's reference build.
   */
  private def validateBuild(build: String): Either[BuildMismatch, Unit] = {
    // Normalize build names for comparison
    val normalizedQuery = normalizeBuildName(build)
    val normalizedVcf = normalizeBuildName(vcfInfo.referenceBuild)

    if (normalizedQuery == normalizedVcf) {
      Right(())
    } else {
      Left(BuildMismatch(expected = build, actual = vcfInfo.referenceBuild))
    }
  }

  /**
   * Normalize reference build names for comparison.
   * e.g., "hg38" -> "GRCh38", "hg19" -> "GRCh37"
   */
  private def normalizeBuildName(build: String): String = {
    build.toLowerCase match {
      case "hg38" | "grch38" | "grch38_full_analysis_set_plus_decoy_hla" => "grch38"
      case "hg19" | "grch37" | "b37" => "grch37"
      case "chm13" | "t2t-chm13" | "chm13v2" => "chm13"
      case other => other.toLowerCase
    }
  }

  /**
   * Convert HTSJDK VariantContext to our VariantCall model.
   */
  private def variantContextToCall(vc: VariantContext): VariantCall = {
    val genotype = if (vc.getGenotypes.isEmpty) "." else {
      val gt = vc.getGenotype(0)
      if (gt.isNoCall) "."
      else gt.getGenotypeString
    }

    val depth = if (vc.getGenotypes.isEmpty) None else {
      val gt = vc.getGenotype(0)
      if (gt.hasDP) Some(gt.getDP) else None
    }

    val altAllele = if (vc.getAlternateAlleles.isEmpty) "."
    else vc.getAlternateAllele(0).getBaseString

    VariantCall(
      contig = vc.getContig,
      position = vc.getStart.toLong,
      ref = vc.getReference.getBaseString,
      alt = altAllele,
      genotype = genotype,
      depth = depth,
      quality = if (vc.hasLog10PError) Some(-10 * vc.getLog10PError) else None,
      filter = if (vc.isFiltered) vc.getFilters.asScala.mkString(",") else "PASS"
    )
  }

  /**
   * Close the VCF reader.
   */
  def close(): Unit = {
    reader.close()
  }
}

object VcfQueryService {

  /**
   * Create a VcfQueryService from cached VCF metadata.
   */
  def fromCache(
                 sampleAccession: String,
                 runId: String,
                 alignmentId: String
               ): Either[String, VcfQueryService] = {
    VcfCache.loadMetadata(sampleAccession, runId, alignmentId).map { info =>
      new VcfQueryService(info)
    }
  }

  /**
   * Create a VcfQueryService from a VCF file path directly.
   * Used for vendor-provided VCFs that aren't in the standard cache structure.
   *
   * @param vcfPath        Path to the VCF file (must be indexed)
   * @param referenceBuild Reference build of the VCF
   * @return VcfQueryService or error
   */
  def fromVcfPath(
                   vcfPath: String,
                   referenceBuild: String
                 ): Either[String, VcfQueryService] = {
    val vcfFile = new File(vcfPath)
    if (!vcfFile.exists()) {
      Left(s"VCF file not found: $vcfPath")
    } else {
      // Create a minimal CachedVcfInfo for the VCF path
      val info = CachedVcfInfo(
        vcfPath = vcfPath,
        indexPath = vcfPath + ".tbi", // Assume tabix index
        referenceBuild = referenceBuild,
        callerVersion = "vendor",
        gatkVersion = "N/A",
        createdAt = java.time.LocalDateTime.now().format(java.time.format.DateTimeFormatter.ISO_LOCAL_DATE_TIME),
        fileSizeBytes = vcfFile.length(),
        variantCount = 0L,
        contigs = List.empty,
        inferredSex = None
      )
      Right(new VcfQueryService(info))
    }
  }

  /**
   * Create a VcfQueryService from AT URIs.
   */
  def fromUris(
                sampleAccession: String,
                sequenceRunUri: Option[String],
                alignmentUri: Option[String]
              ): Either[String, VcfQueryService] = {
    VcfCache.loadMetadataFromUris(sampleAccession, sequenceRunUri, alignmentUri).map { info =>
      new VcfQueryService(info)
    }
  }

  /**
   * Query a position using a cached VCF, automatically loading and closing.
   */
  def quickQuery(
                  sampleAccession: String,
                  runId: String,
                  alignmentId: String,
                  build: String,
                  contig: String,
                  position: Long
                ): Either[String, Option[VariantCall]] = {
    fromCache(sampleAccession, runId, alignmentId).flatMap { service =>
      try {
        service.queryPosition(build, contig, position) match {
          case Left(mismatch) => Left(s"Build mismatch: expected ${mismatch.expected}, got ${mismatch.actual}")
          case Right(result) => Right(result)
        }
      } finally {
        service.close()
      }
    }
  }
}
