package com.decodingus.analysis

import com.decodingus.genotype.model.{TestTypeDefinition, TestTypes}
import com.decodingus.util.Logger
import htsjdk.samtools.{SamReader, SamReaderFactory, ValidationStringency}

import java.io.File
import scala.jdk.CollectionConverters.*
import scala.util.Using

/**
 * Infers test type (WGS, targeted Y, targeted MT, etc.) from BAM/CRAM coverage patterns.
 *
 * Uses BAM index statistics for fast inference without scanning reads:
 * - Targeted Y tests (Big Y, Y Elite): High Y coverage, very low autosomal
 * - Targeted MT tests: High MT coverage, very low autosomal/Y
 * - WGS: Coverage across all chromosomes
 * - WES: High autosomal exonic coverage (detected by coverage depth pattern)
 */
object TestTypeInference {

  private val log = Logger[TestTypeInference.type]

  /**
   * Coverage statistics per chromosome group.
   * Coverage values are estimated from aligned read counts divided by chromosome length.
   */
  case class ChromosomeCoverageStats(
    autosomalCoverage: Double,
    xCoverage: Double,
    yCoverage: Double,
    mtCoverage: Double,
    totalReads: Long,
    autosomalReads: Long,
    yReads: Long,
    mtReads: Long
  ) {
    /** Y:autosome coverage ratio */
    def yAutoRatio: Double = if (autosomalCoverage > 0) yCoverage / autosomalCoverage else 0.0

    /** MT:autosome coverage ratio */
    def mtAutoRatio: Double = if (autosomalCoverage > 0) mtCoverage / autosomalCoverage else 0.0

    /** Whether this looks like targeted Y (high Y, low autosomal) */
    def isLikelyTargetedY: Boolean = yCoverage > 1.0 && autosomalCoverage < 1.0

    /** Whether this looks like targeted MT (high MT, low autosomal/Y) */
    def isLikelyTargetedMt: Boolean = mtCoverage > 10.0 && autosomalCoverage < 1.0 && yCoverage < 1.0

    /** Whether this looks like WGS (coverage across all chromosomes) */
    def isLikelyWgs: Boolean = autosomalCoverage > 1.0
  }

  // Chromosome name patterns
  private val autosomePattern = "^(chr)?([1-9]|1[0-9]|2[0-2])$".r
  private val chrXPattern = "^(chr)?X$".r
  private val chrYPattern = "^(chr)?Y$".r
  private val chrMTPattern = "^(chr)?(M|MT)$".r

  // Average read length assumption for coverage calculation
  // This is a rough estimate; actual coverage will be refined in Phase 2
  private val ASSUMED_READ_LENGTH = 150

  /**
   * Infer test type from BAM index statistics.
   * This is fast (uses index metadata only) and suitable for early detection.
   *
   * @param bamPath Path to BAM/CRAM file
   * @param platform Optional platform hint (PacBio, Nanopore, Illumina, etc.)
   * @param vendor Optional vendor hint for refined test type selection
   * @param meanReadLength Optional mean read length (improves coverage estimate)
   * @return Inferred test type or error
   */
  def inferFromBamIndex(
    bamPath: String,
    platform: Option[String] = None,
    vendor: Option[String] = None,
    meanReadLength: Option[Int] = None
  ): Either[String, TestTypeDefinition] = {
    calculateChromosomeCoverage(bamPath, meanReadLength) match {
      case Right(stats) =>
        val inferred = TestTypes.inferFromCoverage(
          yCoverage = Some(stats.yCoverage),
          mtCoverage = Some(stats.mtCoverage),
          autosomalCoverage = Some(stats.autosomalCoverage),
          totalReads = stats.totalReads,
          vendor = vendor,
          platform = platform,
          meanReadLength = meanReadLength
        )
        log.info(s"Inferred test type: ${inferred.code} (Y=${f"${stats.yCoverage}%.1f"}x, auto=${f"${stats.autosomalCoverage}%.1f"}x, MT=${f"${stats.mtCoverage}%.1f"}x)")
        Right(inferred)
      case Left(error) =>
        Left(error)
    }
  }

  /**
   * Calculate per-chromosome group coverage from BAM index statistics.
   * This is a fast estimate based on aligned read counts per chromosome.
   *
   * @param bamPath Path to BAM/CRAM file
   * @param meanReadLength Optional mean read length (defaults to 150bp)
   * @return Coverage statistics or error
   */
  def calculateChromosomeCoverage(
    bamPath: String,
    meanReadLength: Option[Int] = None
  ): Either[String, ChromosomeCoverageStats] = {
    val bamFile = new File(bamPath)
    if (!bamFile.exists()) {
      return Left(s"BAM file not found: $bamPath")
    }

    Using(SamReaderFactory.makeDefault()
      .validationStringency(ValidationStringency.SILENT)
      .open(bamFile)) { reader =>
      calculateFromReader(reader, meanReadLength)
    }.fold(
      error => Left(s"Failed to read BAM file: ${error.getMessage}"),
      identity
    )
  }

  private def calculateFromReader(
    reader: SamReader,
    meanReadLength: Option[Int]
  ): Either[String, ChromosomeCoverageStats] = {
    val header = reader.getFileHeader
    val sequences = header.getSequenceDictionary.getSequences.asScala.toList

    // Categorize chromosomes
    val autosomes = sequences.filter(s => autosomePattern.findFirstIn(s.getSequenceName).isDefined)
    val chrX = sequences.find(s => chrXPattern.findFirstIn(s.getSequenceName).isDefined)
    val chrY = sequences.find(s => chrYPattern.findFirstIn(s.getSequenceName).isDefined)
    val chrMT = sequences.find(s => chrMTPattern.findFirstIn(s.getSequenceName).isDefined)

    if (autosomes.isEmpty) {
      return Left("No autosomal chromosomes found in BAM header")
    }

    // Calculate total lengths
    val autosomeLength = autosomes.map(_.getSequenceLength.toLong).sum
    val xLength = chrX.map(_.getSequenceLength.toLong).getOrElse(0L)
    val yLength = chrY.map(_.getSequenceLength.toLong).getOrElse(0L)
    val mtLength = chrMT.map(_.getSequenceLength.toLong).getOrElse(16569L) // rCRS length default

    // Get index statistics if available
    if (!reader.hasIndex) {
      return Left("BAM index not available - test type inference requires indexed BAM/CRAM")
    }

    try {
      val index = reader.indexing().getIndex
      val readCounts = header.getSequenceDictionary.getSequences.asScala.map { seq =>
        val meta = index.getMetaData(seq.getSequenceIndex)
        val count: Long = if (meta != null) meta.getAlignedRecordCount else 0L
        seq.getSequenceName -> count
      }.toMap

      val effectiveReadLength = meanReadLength.getOrElse(ASSUMED_READ_LENGTH)

      // Calculate read counts per group
      val autosomeReads = autosomes.map(s => readCounts.getOrElse(s.getSequenceName, 0L)).sum
      val xReads = chrX.map(s => readCounts.getOrElse(s.getSequenceName, 0L)).getOrElse(0L)
      val yReads = chrY.map(s => readCounts.getOrElse(s.getSequenceName, 0L)).getOrElse(0L)
      val mtReads = chrMT.map(s => readCounts.getOrElse(s.getSequenceName, 0L)).getOrElse(0L)
      val totalReads = readCounts.values.sum

      // Calculate coverage estimates: (reads * read_length) / chromosome_length
      val autosomalCoverage = if (autosomeLength > 0)
        (autosomeReads.toDouble * effectiveReadLength) / autosomeLength else 0.0
      val xCoverage = if (xLength > 0)
        (xReads.toDouble * effectiveReadLength) / xLength else 0.0
      val yCoverage = if (yLength > 0)
        (yReads.toDouble * effectiveReadLength) / yLength else 0.0
      val mtCoverage = if (mtLength > 0)
        (mtReads.toDouble * effectiveReadLength) / mtLength else 0.0

      log.debug(s"Coverage stats: autosomal=${f"$autosomalCoverage%.2f"}x ($autosomeReads reads), " +
        s"X=${f"$xCoverage%.2f"}x ($xReads reads), Y=${f"$yCoverage%.2f"}x ($yReads reads), " +
        s"MT=${f"$mtCoverage%.2f"}x ($mtReads reads)")

      Right(ChromosomeCoverageStats(
        autosomalCoverage = autosomalCoverage,
        xCoverage = xCoverage,
        yCoverage = yCoverage,
        mtCoverage = mtCoverage,
        totalReads = totalReads,
        autosomalReads = autosomeReads,
        yReads = yReads,
        mtReads = mtReads
      ))

    } catch {
      case e: Exception =>
        Left(s"Failed to read BAM index: ${e.getMessage}")
    }
  }

  /**
   * Check if the inferred test type differs significantly from an existing type.
   * Used to determine if we should update the test type after Phase 2 analysis.
   */
  def shouldUpdateTestType(current: String, inferred: TestTypeDefinition): Boolean = {
    val currentUpper = current.toUpperCase
    val inferredCode = inferred.code.toUpperCase

    // If both are WGS variants, don't update unless significantly different
    val wgsVariants = Set("WGS", "WGS_HIFI", "WGS_NANOPORE", "WGS_CLR", "WGS_LOW_PASS")
    val targetedY = Set("BIG_Y_500", "BIG_Y_700", "Y_ELITE", "Y_PRIME")
    val targetedMt = Set("MT_FULL_SEQUENCE", "MT_PLUS", "MT_CR_ONLY")

    // Update if switching between major categories
    if (wgsVariants.contains(currentUpper) && targetedY.contains(inferredCode)) true
    else if (wgsVariants.contains(currentUpper) && targetedMt.contains(inferredCode)) true
    else if (targetedY.contains(currentUpper) && wgsVariants.contains(inferredCode)) true
    else if (targetedMt.contains(currentUpper) && wgsVariants.contains(inferredCode)) true
    // Update if going from generic WGS to specific WGS type
    else if (currentUpper == "WGS" && inferredCode != "WGS" && wgsVariants.contains(inferredCode)) true
    // Update if going from "Unknown" to anything
    else if (currentUpper == "UNKNOWN") true
    else false
  }
}
