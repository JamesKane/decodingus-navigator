package com.decodingus.analysis.sv

import io.circe.Codec
import io.circe.generic.semiauto.*

/**
 * Evidence of a discordant read pair indicating a potential SV.
 *
 * @param readName   Name of the read
 * @param chrom1     Chromosome of first read
 * @param pos1       Position of first read
 * @param strand1    Strand of first read ('+' or '-')
 * @param chrom2     Chromosome of second read (mate)
 * @param pos2       Position of second read (mate)
 * @param strand2    Strand of second read (mate)
 * @param insertSize Observed insert size
 * @param mapQ       Minimum mapping quality of the pair
 * @param reason     Why this pair is considered discordant
 */
case class DiscordantPair(
  readName: String,
  chrom1: String,
  pos1: Long,
  strand1: Char,
  chrom2: String,
  pos2: Long,
  strand2: Char,
  insertSize: Int,
  mapQ: Int,
  reason: DiscordantReason
)

/**
 * Reason a read pair is considered discordant.
 */
enum DiscordantReason:
  case InsertSizeOutlier   // Insert size significantly larger/smaller than expected
  case WrongOrientation    // Unexpected pair orientation (not FR for standard library)
  case InterChromosomal    // Reads map to different chromosomes

object DiscordantReason:
  import io.circe.{Encoder, Decoder}
  given Encoder[DiscordantReason] = Encoder.encodeString.contramap(_.toString)
  given Decoder[DiscordantReason] = Decoder.decodeString.emap { s =>
    DiscordantReason.values.find(_.toString == s).toRight(s"Unknown DiscordantReason: $s")
  }

/**
 * Evidence of a split read indicating a potential SV breakpoint.
 *
 * A split read has its alignment split across two locations,
 * often indicated by the SA (supplementary alignment) tag.
 *
 * @param readName      Name of the read
 * @param primaryChrom  Chromosome of primary alignment
 * @param primaryPos    Position of primary alignment
 * @param primaryStrand Strand of primary alignment
 * @param suppChrom     Chromosome of supplementary alignment
 * @param suppPos       Position of supplementary alignment
 * @param suppStrand    Strand of supplementary alignment
 * @param clipLength    Length of the clipped portion
 * @param mapQ          Minimum mapping quality
 */
case class SplitRead(
  readName: String,
  primaryChrom: String,
  primaryPos: Long,
  primaryStrand: Char,
  suppChrom: String,
  suppPos: Long,
  suppStrand: Char,
  clipLength: Int,
  mapQ: Int
)

/**
 * Read depth information for a genomic bin.
 *
 * @param chrom     Chromosome
 * @param binStart  Start position of the bin (0-based)
 * @param binEnd    End position of the bin (exclusive)
 * @param readCount Number of reads starting in this bin
 */
case class DepthBin(
  chrom: String,
  binStart: Long,
  binEnd: Long,
  readCount: Int
)

/**
 * A segment of the genome with abnormal copy number.
 *
 * @param chrom       Chromosome
 * @param start       Start position
 * @param end         End position
 * @param meanDepth   Mean depth in the segment
 * @param log2Ratio   Log2 ratio of observed/expected depth
 * @param zScore      Z-score of the depth deviation
 * @param numBins     Number of bins in the segment
 * @param svType      Inferred SV type (DEL or DUP)
 */
case class DepthSegment(
  chrom: String,
  start: Long,
  end: Long,
  meanDepth: Double,
  log2Ratio: Double,
  zScore: Double,
  numBins: Int,
  svType: SvType
)

/**
 * Collection of all SV evidence gathered from a BAM file.
 *
 * @param discordantPairs List of discordant read pairs
 * @param splitReads      List of split reads
 * @param depthBins       Per-chromosome depth bins (contig -> bins)
 * @param sampleName      Name of the sample
 * @param expectedInsertSize Expected insert size (from library)
 * @param insertSizeSd    Insert size standard deviation
 */
case class SvEvidenceCollection(
  discordantPairs: List[DiscordantPair],
  splitReads: List[SplitRead],
  depthBins: Map[String, Array[Int]],
  sampleName: String,
  expectedInsertSize: Double,
  insertSizeSd: Double
):
  /**
   * Total number of discordant pairs.
   */
  def totalDiscordantPairs: Long = discordantPairs.size.toLong

  /**
   * Total number of split reads.
   */
  def totalSplitReads: Long = splitReads.size.toLong

  /**
   * Get discordant pairs that suggest inter-chromosomal events.
   */
  def interChromosomalPairs: List[DiscordantPair] =
    discordantPairs.filter(_.reason == DiscordantReason.InterChromosomal)

  /**
   * Get discordant pairs grouped by approximate breakpoint location.
   */
  def groupDiscordantPairsByBreakpoint(maxDistance: Int): Map[(String, Long), List[DiscordantPair]] =
    discordantPairs.groupBy { dp =>
      // Group by chromosome and binned position
      (dp.chrom1, dp.pos1 / maxDistance * maxDistance)
    }

  /**
   * Get split reads grouped by approximate breakpoint location.
   */
  def groupSplitReadsByBreakpoint(maxDistance: Int): Map[(String, Long), List[SplitRead]] =
    splitReads.groupBy { sr =>
      // Group by chromosome and binned position
      (sr.primaryChrom, sr.primaryPos / maxDistance * maxDistance)
    }

/**
 * Breakpoint evidence cluster - grouped evidence supporting a single breakpoint.
 *
 * @param chrom           Chromosome of the breakpoint
 * @param position        Estimated position of the breakpoint
 * @param ciLow           Confidence interval - low bound
 * @param ciHigh          Confidence interval - high bound
 * @param discordantPairs Supporting discordant pairs
 * @param splitReads      Supporting split reads
 * @param mateChrom       For inter-chromosomal: mate chromosome
 * @param matePosition    For inter-chromosomal: mate position
 */
case class BreakpointCluster(
  chrom: String,
  position: Long,
  ciLow: Int,
  ciHigh: Int,
  discordantPairs: List[DiscordantPair],
  splitReads: List[SplitRead],
  mateChrom: Option[String] = None,
  matePosition: Option[Long] = None
):
  /**
   * Total evidence support.
   */
  def totalSupport: Int = discordantPairs.size + splitReads.size

  /**
   * Paired-end support count.
   */
  def peSupport: Int = discordantPairs.size

  /**
   * Split-read support count.
   */
  def srSupport: Int = splitReads.size

  /**
   * Average mapping quality of supporting evidence.
   */
  def meanMapQ: Double =
    val allMapQ = discordantPairs.map(_.mapQ) ++ splitReads.map(_.mapQ)
    if allMapQ.isEmpty then 0.0
    else allMapQ.sum.toDouble / allMapQ.size
