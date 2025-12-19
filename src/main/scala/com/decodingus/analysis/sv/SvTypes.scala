package com.decodingus.analysis.sv

import io.circe.{Codec, Decoder, Encoder}
import io.circe.generic.semiauto.*

import java.time.Instant

/**
 * Structural variant types.
 */
enum SvType:
  case DEL   // Deletion
  case DUP   // Duplication
  case INV   // Inversion
  case BND   // Breakend (translocation)
  case INS   // Insertion

object SvType:
  given Encoder[SvType] = Encoder.encodeString.contramap(_.toString)
  given Decoder[SvType] = Decoder.decodeString.emap { s =>
    SvType.values.find(_.toString == s).toRight(s"Unknown SvType: $s")
  }

/**
 * A called structural variant.
 *
 * @param id            Unique identifier for this SV
 * @param chrom         Chromosome of the SV
 * @param start         Start position (1-based)
 * @param end           End position (1-based)
 * @param svType        Type of structural variant
 * @param svLen         Length of the SV (negative for deletions)
 * @param ciPos         Confidence interval around start position (left, right)
 * @param ciEnd         Confidence interval around end position (left, right)
 * @param quality       Quality score (Phred-scaled)
 * @param pairedEndSupport Number of supporting discordant read pairs
 * @param splitReadSupport Number of supporting split reads
 * @param relativeDepth For CNVs: observed/expected depth ratio
 * @param mateChrom     For translocations: mate chromosome
 * @param matePos       For translocations: mate position
 * @param filter        Filter status (PASS or filter reason)
 * @param genotype      Genotype call (0/1, 1/1, etc.)
 */
case class SvCall(
  id: String,
  chrom: String,
  start: Long,
  end: Long,
  svType: SvType,
  svLen: Long,
  ciPos: (Int, Int),
  ciEnd: (Int, Int),
  quality: Double,
  pairedEndSupport: Int,
  splitReadSupport: Int,
  relativeDepth: Option[Double],
  mateChrom: Option[String],
  matePos: Option[Long],
  filter: String,
  genotype: String
) derives Codec.AsObject

object SvCall:
  /**
   * Calculate a confidence score for an SV call based on evidence.
   * Returns a value between 0.0 and 1.0.
   */
  def calculateConfidence(call: SvCall): Double =
    val peWeight = 0.3
    val srWeight = 0.4
    val depthWeight = 0.3

    val peScore = math.min(call.pairedEndSupport / 10.0, 1.0)
    val srScore = math.min(call.splitReadSupport / 5.0, 1.0)
    val depthScore = call.relativeDepth.map { rd =>
      val deviation = math.abs(1.0 - rd)
      math.min(deviation / 0.5, 1.0) // Full confidence at 50% deviation
    }.getOrElse(0.0)

    peScore * peWeight + srScore * srWeight + depthScore * depthWeight

/**
 * Result of structural variant analysis.
 *
 * @param svCalls              List of called SVs
 * @param totalDiscordantPairs Total discordant pairs observed
 * @param totalSplitReads      Total split reads observed
 * @param cnvSegments          Number of CNV segments detected
 * @param analysisTimestamp    When the analysis was performed
 * @param referenceBuild       Reference genome build used
 * @param meanCoverage         Mean coverage of the sample
 */
case class SvAnalysisResult(
  svCalls: List[SvCall],
  totalDiscordantPairs: Long,
  totalSplitReads: Long,
  cnvSegments: Int,
  analysisTimestamp: Instant,
  referenceBuild: String,
  meanCoverage: Double
)

object SvAnalysisResult:
  given Encoder[Instant] = Encoder.encodeString.contramap(_.toString)
  given Decoder[Instant] = Decoder.decodeString.emap { s =>
    try Right(Instant.parse(s))
    catch case e: Exception => Left(s"Invalid timestamp: $s")
  }
  given Codec[SvAnalysisResult] = deriveCodec

/**
 * Metadata about a cached SV VCF file.
 *
 * @param vcfPath          Path to the VCF file
 * @param indexPath        Path to the VCF index file
 * @param referenceBuild   Reference genome build
 * @param createdAt        When the VCF was created
 * @param svCallCount      Total number of SV calls
 * @param deletionCount    Number of deletions
 * @param duplicationCount Number of duplications
 * @param inversionCount   Number of inversions
 * @param translocationCount Number of translocations
 */
case class CachedSvInfo(
  vcfPath: String,
  indexPath: String,
  referenceBuild: String,
  createdAt: String,
  svCallCount: Int,
  deletionCount: Int,
  duplicationCount: Int,
  inversionCount: Int,
  translocationCount: Int
) derives Codec.AsObject

/**
 * Configuration for SV calling.
 */
case class SvCallerConfig(
  // Depth-based CNV detection
  binSize: Int = 1000,                    // Size of depth bins in bp
  minDepthZScore: Double = 2.5,           // Minimum z-score for CNV call
  minCnvSize: Int = 10000,                // Minimum CNV size to report (10kb)

  // Discordant pair detection
  insertSizeZThreshold: Double = 4.0,     // Insert sizes > mean + 4*SD are discordant
  minMapQ: Int = 20,                      // Minimum mapping quality

  // Evidence clustering
  maxClusterDistance: Int = 500,          // Max distance to cluster evidence
  minPairedEndSupport: Int = 2,           // Minimum PE support for call
  minSplitReadSupport: Int = 1,           // Minimum SR support for call
  minTotalSupport: Int = 3,               // Minimum total support for non-depth calls

  // Quality thresholds
  minQuality: Double = 10.0               // Minimum quality score to report
) derives Codec.AsObject

object SvCallerConfig:
  val default: SvCallerConfig = SvCallerConfig()
