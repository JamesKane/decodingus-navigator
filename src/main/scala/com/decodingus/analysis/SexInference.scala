package com.decodingus.analysis

import htsjdk.samtools.{SamReader, SamReaderFactory, ValidationStringency}

import java.io.File
import scala.jdk.CollectionConverters.*
import scala.util.{Either, Left, Right, Using}

/**
 * Infers biological sex from BAM/CRAM coverage on chrX vs autosomes.
 *
 * Rationale:
 * - Males (XY) have ~0.5x coverage on chrX relative to autosomes
 * - Females (XX) have ~1.0x coverage on chrX relative to autosomes
 *
 * This is used to determine ploidy for variant calling:
 * - chrX: diploid for females, haploid for males
 * - chrY: skip for females, haploid for males
 */
object SexInference {

  /** Inferred sex with confidence level */
  enum InferredSex:
    case Male
    case Female
    case Unknown

  case class SexInferenceResult(
                                 inferredSex: InferredSex,
                                 xAutosomeRatio: Double,
                                 autosomeMeanCoverage: Double,
                                 xCoverage: Double,
                                 confidence: String // "high", "medium", "low"
                               ) {
    def isMale: Boolean = inferredSex == InferredSex.Male

    def isFemale: Boolean = inferredSex == InferredSex.Female

    def isUnknown: Boolean = inferredSex == InferredSex.Unknown
  }

  // Thresholds for sex determination
  private val MaleRatioThreshold = 0.65 // X:autosome ratio below this suggests male
  private val FemaleRatioThreshold = 0.85 // X:autosome ratio above this suggests female
  private val MinAutosomeCoverage = 5.0 // Minimum coverage for confident inference

  // Pattern for autosomal chromosomes (chr1-22 or 1-22)
  private val autosomePattern = "^(chr)?([1-9]|1[0-9]|2[0-2])$".r

  // Pattern for chrX
  private val chrXPattern = "^(chr)?X$".r

  /**
   * Infer sex from a BAM/CRAM file by comparing chrX coverage to autosomal coverage.
   *
   * @param bamPath    Path to BAM/CRAM file
   * @param onProgress Optional progress callback
   * @return Sex inference result or error
   */
  def inferFromBam(
                    bamPath: String,
                    onProgress: (String, Double) => Unit = (_, _) => ()
                  ): Either[String, SexInferenceResult] = {
    onProgress("Opening BAM file for sex inference...", 0.0)

    val bamFile = new File(bamPath)
    if (!bamFile.exists()) {
      return Left(s"BAM file not found: $bamPath")
    }

    Using(SamReaderFactory.makeDefault()
      .validationStringency(ValidationStringency.SILENT)
      .open(bamFile)) { reader =>
      inferFromReader(reader, onProgress)
    }.fold(
      error => Left(s"Failed to read BAM file: ${error.getMessage}"),
      identity
    )
  }

  /**
   * Infer sex from an open SamReader.
   */
  private def inferFromReader(
                               reader: SamReader,
                               onProgress: (String, Double) => Unit
                             ): Either[String, SexInferenceResult] = {
    onProgress("Reading sequence dictionary...", 0.1)

    val header = reader.getFileHeader
    val sequences = header.getSequenceDictionary.getSequences.asScala.toList

    // Get autosomes and chrX
    val autosomes = sequences.filter(s => autosomePattern.findFirstIn(s.getSequenceName).isDefined)
    val chrXSeq = sequences.find(s => chrXPattern.findFirstIn(s.getSequenceName).isDefined)

    if (autosomes.isEmpty) {
      return Left("No autosomal chromosomes found in BAM header")
    }

    chrXSeq match {
      case None => Left("chrX not found in BAM header")
      case Some(chrX) =>
        onProgress("Calculating autosomal coverage...", 0.2)

        // Calculate total autosome length
        val autosomeLength = autosomes.map(_.getSequenceLength.toLong).sum

        // Get index statistics if available (much faster than scanning reads)
        val indexStats: Option[Map[String, Long]] = if (reader.hasIndex) {
          try {
            val index = reader.indexing().getIndex
            Some(header.getSequenceDictionary.getSequences.asScala.map { seq =>
              val meta = index.getMetaData(seq.getSequenceIndex)
              val count: Long = if (meta != null) meta.getAlignedRecordCount else 0L
              seq.getSequenceName -> count
            }.toMap)
          } catch {
            case _: Exception => None
          }
        } else None

        indexStats match {
          case Some(stats) =>
            // Use index statistics (fast path)
            calculateFromIndexStats(stats, autosomes.map(_.getSequenceName), chrX.getSequenceName,
              autosomeLength, chrX.getSequenceLength.toLong, onProgress)
          case None =>
            // Fall back to scanning reads (slow path)
            Left("BAM index not available - sex inference requires indexed BAM/CRAM")
        }
    }
  }

  private def calculateFromIndexStats(
                                       readCounts: Map[String, Long],
                                       autosomeNames: List[String],
                                       chrXName: String,
                                       autosomeLength: Long,
                                       chrXLength: Long,
                                       onProgress: (String, Double) => Unit
                                     ): Either[String, SexInferenceResult] = {
    onProgress("Computing coverage ratios...", 0.7)

    val autosomeReads = autosomeNames.map(name => readCounts.getOrElse(name, 0L)).sum
    val chrXReads = readCounts.getOrElse(chrXName, 0L)

    if (autosomeReads == 0) {
      return Left("No autosomal reads found - cannot infer sex")
    }

    // Normalize by length to get coverage
    val autosomeCoverage = autosomeReads.toDouble / autosomeLength * 100 // reads per 100bp
    val chrXCoverage = chrXReads.toDouble / chrXLength * 100

    val ratio = if (autosomeCoverage > 0) chrXCoverage / autosomeCoverage else 0.0

    onProgress("Sex inference complete", 1.0)

    val (sex, confidence) = determineSex(ratio, autosomeCoverage)

    Right(SexInferenceResult(
      inferredSex = sex,
      xAutosomeRatio = ratio,
      autosomeMeanCoverage = autosomeCoverage,
      xCoverage = chrXCoverage,
      confidence = confidence
    ))
  }

  private def determineSex(ratio: Double, autosomeCoverage: Double): (InferredSex, String) = {
    // Low coverage = low confidence
    if (autosomeCoverage < MinAutosomeCoverage) {
      if (ratio < MaleRatioThreshold) {
        (InferredSex.Male, "low")
      } else if (ratio > FemaleRatioThreshold) {
        (InferredSex.Female, "low")
      } else {
        (InferredSex.Unknown, "low")
      }
    } else {
      // High coverage = can determine with more confidence
      if (ratio < MaleRatioThreshold) {
        // Clear male signal
        val conf = if (ratio < 0.55) "high" else "medium"
        (InferredSex.Male, conf)
      } else if (ratio > FemaleRatioThreshold) {
        // Clear female signal
        val conf = if (ratio > 0.95) "high" else "medium"
        (InferredSex.Female, conf)
      } else {
        // Ambiguous
        (InferredSex.Unknown, "low")
      }
    }
  }

  /**
   * Get the ploidy to use for a given contig based on inferred sex.
   *
   * @param contigName The contig name (e.g., "chrX", "chrY", "chr1")
   * @param sexResult  The sex inference result
   * @return Ploidy (1 for haploid, 2 for diploid), or None if contig should be skipped
   */
  def ploidyForContig(contigName: String, sexResult: SexInferenceResult): Option[Int] = {
    val chrXPattern = "^(chr)?X$".r
    val chrYPattern = "^(chr)?Y$".r
    val chrMPattern = "^(chr)?(M|MT)$".r

    contigName match {
      case chrXPattern(_) =>
        // chrX: diploid for females, haploid for males
        sexResult.inferredSex match {
          case InferredSex.Female => Some(2)
          case InferredSex.Male => Some(1)
          case InferredSex.Unknown => Some(2) // Default to diploid if unknown
        }
      case chrYPattern(_) =>
        // chrY: skip for females, haploid for males
        sexResult.inferredSex match {
          case InferredSex.Female => None // Skip chrY for females
          case InferredSex.Male => Some(1)
          case InferredSex.Unknown => Some(1) // Include chrY but haploid if unknown
        }
      case chrMPattern(_, _) =>
        // Mitochondrial DNA is always haploid
        Some(1)
      case _ =>
        // Autosomes are always diploid
        Some(2)
    }
  }
}
