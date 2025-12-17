package com.decodingus.liftover

import com.decodingus.refgenome.LiftoverGateway
import htsjdk.samtools.liftover.LiftOver
import htsjdk.samtools.util.Interval

/**
 * Lifted genotype call with position in target coordinates.
 *
 * @param chromosome             Target chromosome name
 * @param position               Target position (1-based)
 * @param allele1                First allele (reverse-complemented if needed)
 * @param allele2                Second allele (reverse-complemented if needed)
 * @param wasReverseComplemented True if alleles were reverse-complemented
 */
case class LiftedGenotype(
                           chromosome: String,
                           position: Int,
                           allele1: Char,
                           allele2: Char,
                           wasReverseComplemented: Boolean
                         )

/**
 * Result of lifting a batch of genotypes.
 *
 * @param lifted      Successfully lifted genotypes
 * @param failedCount Number of positions that couldn't be lifted
 */
case class LiftoverResult(
                           lifted: List[(Long, LiftedGenotype)], // Original position -> Lifted genotype
                           failedCount: Int
                         )

/**
 * Utility for lifting genotype coordinates between reference builds using htsjdk.
 *
 * Handles:
 * - Coordinate conversion via UCSC chain files
 * - Reverse-complement of alleles when mapping to negative strand
 *
 * @see https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/liftover/LiftOver.java
 */
class GenotypeLiftover(fromBuild: String, toBuild: String) {

  private val liftoverGateway = new LiftoverGateway((_, _) => ())

  // Complement map for reverse-complement operations
  private val complementMap: Map[Char, Char] = Map(
    'A' -> 'T', 'T' -> 'A',
    'C' -> 'G', 'G' -> 'C',
    'a' -> 't', 't' -> 'a',
    'c' -> 'g', 'g' -> 'c',
    'N' -> 'N', 'n' -> 'n',
    '-' -> '-', '0' -> '0'
  )

  /**
   * Initialize the liftover chain file.
   *
   * @return Either error message or the initialized LiftOver instance
   */
  def initialize(): Either[String, LiftOver] = {
    liftoverGateway.resolve(fromBuild, toBuild).map { chainPath =>
      new LiftOver(chainPath.toFile)
    }
  }

  // Track first few failures for debugging
  private var debugFailureCount = 0
  private val maxDebugFailures = 5

  /**
   * Lift a single genotype position.
   *
   * @param liftOver   The initialized LiftOver instance
   * @param chromosome Source chromosome (e.g., "Y", "chrY", "24")
   * @param position   Source position (1-based)
   * @param allele1    First allele
   * @param allele2    Second allele
   * @return Some(LiftedGenotype) if successful, None if position couldn't be lifted
   */
  def liftGenotype(
                    liftOver: LiftOver,
                    chromosome: String,
                    position: Int,
                    allele1: Char,
                    allele2: Char
                  ): Option[LiftedGenotype] = {
    // Normalize chromosome name for liftover (htsjdk expects certain formats)
    val normalizedChr = normalizeChromosome(chromosome, fromBuild)

    // Create interval for the SNP position (1-based, inclusive)
    val sourceInterval = new Interval(normalizedChr, position, position, false, s"$normalizedChr:$position")

    val result = Option(liftOver.liftOver(sourceInterval))

    // Debug logging for first few failures
    if (result.isEmpty && debugFailureCount < maxDebugFailures) {
      println(s"[GenotypeLiftover] Failed to lift: $chromosome -> $normalizedChr:$position")
      debugFailureCount += 1
    }

    result.map { targetInterval =>
      val needsReverseComplement = targetInterval.isNegativeStrand

      val (liftedAllele1, liftedAllele2) = if (needsReverseComplement) {
        (complement(allele1), complement(allele2))
      } else {
        (allele1, allele2)
      }

      LiftedGenotype(
        chromosome = targetInterval.getContig,
        position = targetInterval.getStart,
        allele1 = liftedAllele1,
        allele2 = liftedAllele2,
        wasReverseComplemented = needsReverseComplement
      )
    }
  }

  /**
   * Lift a batch of genotypes efficiently.
   *
   * @param liftOver   The initialized LiftOver instance
   * @param genotypes  List of (chromosome, position, allele1, allele2)
   * @param onProgress Progress callback (lifted count, total count)
   * @return LiftoverResult with lifted genotypes and failure count
   */
  def liftGenotypes(
                     liftOver: LiftOver,
                     genotypes: List[(String, Int, Char, Char)],
                     onProgress: (Int, Int) => Unit = (_, _) => ()
                   ): LiftoverResult = {
    val total = genotypes.size
    var lifted = List.empty[(Long, LiftedGenotype)]
    var failedCount = 0
    var processedCount = 0

    genotypes.foreach { case (chr, pos, a1, a2) =>
      liftGenotype(liftOver, chr, pos, a1, a2) match {
        case Some(liftedGeno) =>
          lifted = (pos.toLong, liftedGeno) :: lifted
        case None =>
          failedCount += 1
      }

      processedCount += 1
      if (processedCount % 10000 == 0) {
        onProgress(processedCount, total)
      }
    }

    onProgress(total, total)
    LiftoverResult(lifted.reverse, failedCount)
  }

  /**
   * Get the complement of a nucleotide.
   */
  private def complement(allele: Char): Char = {
    complementMap.getOrElse(allele, allele)
  }

  /**
   * Normalize chromosome name for liftover.
   * UCSC chain files use "chr" prefix and "chrM" (not chrMT) for mitochondria.
   */
  private def normalizeChromosome(chr: String, build: String): String = {
    val normalized = chr.toUpperCase.stripPrefix("CHR")

    // Convert numeric/letter chromosomes to UCSC format
    // Chain files use chrM (not chrMT) for mitochondria
    val baseChr = normalized match {
      case "23" | "X" => "X"
      case "24" | "Y" => "Y"
      case "25" | "MT" | "M" | "MITO" => "M" // UCSC uses chrM, not chrMT
      case other => other
    }

    // All UCSC chain files use "chr" prefix
    s"chr$baseChr"
  }
}

object GenotypeLiftover {

  /**
   * Quick check if liftover is needed between two builds.
   */
  def needsLiftover(fromBuild: String, toBuild: String): Boolean = {
    normalizeBuildName(fromBuild) != normalizeBuildName(toBuild)
  }

  /**
   * Normalize build name for comparison.
   */
  def normalizeBuildName(build: String): String = {
    build.toLowerCase match {
      case "grch37" | "hg19" | "b37" => "GRCh37"
      case "grch38" | "hg38" | "b38" => "GRCh38"
      case "chm13v2" | "chm13" | "t2t" | "hs1" => "CHM13v2"
      case other => other
    }
  }
}
