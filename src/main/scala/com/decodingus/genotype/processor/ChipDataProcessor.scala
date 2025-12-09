package com.decodingus.genotype.processor

import com.decodingus.genotype.model.*
import com.decodingus.genotype.parser.{ChipDataParser, FormatDetectionResult}

import java.io.File
import java.nio.file.{Files, Path}
import java.security.MessageDigest
import scala.util.Using

/**
 * Result of processing a chip data file.
 */
case class ChipProcessingResult(
  parser: ChipDataParser,
  detection: FormatDetectionResult,
  summary: GenotypingTestSummary,
  genotypeCalls: List[GenotypeCall],
  yDnaCalls: List[GenotypeCall],
  mtDnaCalls: List[GenotypeCall],
  autosomalCalls: List[GenotypeCall]
) {
  /**
   * Total markers parsed.
   */
  def totalMarkers: Int = genotypeCalls.size

  /**
   * Check if data is suitable for ancestry analysis.
   */
  def isSuitableForAncestry: Boolean = summary.isAcceptableForAncestry

  /**
   * Check if data is suitable for Y-DNA haplogroup determination.
   */
  def isSuitableForYHaplogroup: Boolean = summary.hasSufficientYCoverage

  /**
   * Check if data is suitable for mtDNA haplogroup determination.
   */
  def isSuitableForMtHaplogroup: Boolean = summary.hasSufficientMtCoverage
}

/**
 * Processor for chip/array genotype data files.
 *
 * This is the main entry point for processing raw chip data exports
 * from vendors like 23andMe, AncestryDNA, FTDNA, etc.
 *
 * All processing happens locally - raw genotype data never leaves
 * the user's machine. Only metadata and Y/mtDNA variants can be
 * submitted to DecodingUs for tree building.
 *
 * @see multi-test-type-roadmap.md for architecture details
 */
class ChipDataProcessor {

  /**
   * Process a chip data file.
   *
   * @param file The raw data export file
   * @param onProgress Progress callback (message, current, total)
   * @return Either error message or processing result
   */
  def process(
    file: File,
    onProgress: (String, Double, Double) => Unit = (_, _, _) => ()
  ): Either[String, ChipProcessingResult] = {

    onProgress("Detecting file format...", 0.0, 1.0)

    // Calculate file hash for deduplication
    val fileHash = calculateFileHash(file)

    // Detect format and parse
    ChipDataParser.detectParser(file).flatMap { case (parser, detection) =>
      onProgress(s"Detected ${parser.vendor} format. Parsing genotypes...", 0.1, 1.0)

      var lastProgress = 0.1
      parser.parse(file, (current, total) => {
        val progress = 0.1 + (current.toDouble / total) * 0.7
        if (progress - lastProgress > 0.05) {
          lastProgress = progress
          onProgress(s"Parsing... ${current} markers", progress, 1.0)
        }
      }).map { callsIterator =>
        val calls = callsIterator.toList

        onProgress("Computing summary statistics...", 0.85, 1.0)

        // Partition by chromosome type
        val yDnaCalls = calls.filter(_.isYChromosome)
        val mtDnaCalls = calls.filter(_.isMitochondrial)
        val autosomalCalls = calls.filter(_.isAutosomal)

        // Compute summary
        val summary = GenotypingTestSummary.fromCalls(
          calls = calls,
          testType = detection.testType.getOrElse(TestTypes.ARRAY_23ANDME_V5),
          chipVersion = detection.chipVersion,
          sourceFileHash = fileHash
        )

        onProgress("Processing complete.", 1.0, 1.0)

        ChipProcessingResult(
          parser = parser,
          detection = detection,
          summary = summary,
          genotypeCalls = calls,
          yDnaCalls = yDnaCalls,
          mtDnaCalls = mtDnaCalls,
          autosomalCalls = autosomalCalls
        )
      }
    }
  }

  /**
   * Process a file and return only the summary (lighter weight for quick checks).
   */
  def processSummaryOnly(
    file: File,
    onProgress: (String, Double, Double) => Unit = (_, _, _) => ()
  ): Either[String, GenotypingTestSummary] = {
    process(file, onProgress).map(_.summary)
  }

  /**
   * Convert chip genotypes to a format suitable for ancestry analysis.
   *
   * Returns a map of SNP ID (chr:pos format) to numeric genotype value:
   * - 0 = homozygous reference
   * - 1 = heterozygous
   * - 2 = homozygous alternate
   * - -1 = no call
   *
   * @param result The chip processing result
   * @param referenceAlleles Map of chr:pos to reference allele
   */
  def toAncestryGenotypes(
    result: ChipProcessingResult,
    referenceAlleles: Map[String, Char]
  ): Map[String, Int] = {
    result.autosomalCalls.map { call =>
      val snpId = s"${normalizeChromosome(call.chromosome)}:${call.position}"
      val refAllele = referenceAlleles.getOrElse(snpId, call.allele1) // Default to first allele if unknown
      val genotype = call.numericGenotype(refAllele)
      snpId -> genotype
    }.toMap
  }

  /**
   * Extract Y-DNA variants suitable for haplogroup determination.
   *
   * @param result The chip processing result
   * @return List of variant calls with position and allele info
   */
  def extractYDnaVariants(result: ChipProcessingResult): List[ChipVariantCall] = {
    result.yDnaCalls
      .filterNot(_.isNoCall)
      .map { call =>
        ChipVariantCall(
          chromosome = "Y",
          position = call.position,
          rsId = Some(call.markerId).filter(_.startsWith("rs")),
          allele = call.allele1.toString, // Y is haploid
          isHaploid = true
        )
      }
  }

  /**
   * Extract mtDNA variants suitable for haplogroup determination.
   *
   * @param result The chip processing result
   * @return List of variant calls with position and allele info
   */
  def extractMtDnaVariants(result: ChipProcessingResult): List[ChipVariantCall] = {
    result.mtDnaCalls
      .filterNot(_.isNoCall)
      .map { call =>
        // mtDNA can show heteroplasmy but is typically reported as haploid
        ChipVariantCall(
          chromosome = "MT",
          position = call.position,
          rsId = Some(call.markerId).filter(_.startsWith("rs")),
          allele = call.allele1.toString,
          isHaploid = true
        )
      }
  }

  /**
   * Calculate SHA-256 hash of file contents.
   */
  private def calculateFileHash(file: File): Option[String] = {
    scala.util.Try {
      val bytes = Files.readAllBytes(file.toPath)
      val digest = MessageDigest.getInstance("SHA-256")
      val hash = digest.digest(bytes)
      hash.map("%02x".format(_)).mkString
    }.toOption
  }

  /**
   * Normalize chromosome name to standard format (1-22, X, Y, MT).
   */
  private def normalizeChromosome(chr: String): String = {
    val cleaned = chr.toLowerCase.stripPrefix("chr")
    cleaned match {
      case "m" | "mt" | "mito" => "MT"
      case "x" => "X"
      case "y" => "Y"
      case n if n.toIntOption.isDefined => n
      case _ => chr
    }
  }
}

/**
 * A variant call extracted from chip data for haplogroup analysis.
 */
case class ChipVariantCall(
  chromosome: String,
  position: Int,
  rsId: Option[String],
  allele: String,
  isHaploid: Boolean
)

object ChipDataProcessor {
  /**
   * Supported file extensions for chip data.
   */
  val supportedExtensions: Set[String] = Set(".txt", ".csv")

  /**
   * Check if a file appears to be a chip data file.
   */
  def isChipDataFile(file: File): Boolean = {
    val name = file.getName.toLowerCase
    supportedExtensions.exists(name.endsWith)
  }
}
