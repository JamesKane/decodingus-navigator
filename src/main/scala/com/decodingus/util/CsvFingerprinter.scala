package com.decodingus.util

import java.io.{BufferedReader, File, FileInputStream, InputStreamReader}
import java.util.zip.GZIPInputStream
import scala.util.{Try, Using}

/**
 * Detected CSV file type based on content fingerprinting.
 */
sealed trait CsvFileType {
  def description: String
}

object CsvFileType {
  /** Y-STR profile data (FTDNA, YSEQ, etc.) */
  case object StrProfile extends CsvFileType {
    val description = "Y-STR Profile"
  }

  /** Chip/SNP array data (23andMe, AncestryDNA, etc.) */
  case class ChipData(vendor: Option[String]) extends CsvFileType {
    val description = vendor.map(v => s"$v Chip Data").getOrElse("Chip/SNP Data")
  }

  /** VCF variant file (detected by extension, not content) */
  case object VcfVariants extends CsvFileType {
    val description = "VCF Variants"
  }

  /** BAM/CRAM alignment file (detected by extension) */
  case object Alignment extends CsvFileType {
    val description = "Alignment File"
  }

  /** Unknown or unrecognized format */
  case object Unknown extends CsvFileType {
    val description = "Unknown Format"
  }
}

/**
 * Fingerprints CSV/text files to detect their type based on content patterns.
 *
 * Detects:
 * - Y-STR profiles (FTDNA, YSEQ, generic)
 * - Chip/SNP data (23andMe, AncestryDNA, MyHeritage, LivingDNA, FTDNA)
 * - VCF files (by extension)
 * - BAM/CRAM files (by extension)
 */
object CsvFingerprinter {

  private val log = Logger[CsvFingerprinter.type]

  // Number of lines to sample for fingerprinting
  private val SampleLines = 50

  /**
   * Fingerprints a file to detect its type.
   *
   * @param file The file to analyze
   * @return The detected file type
   */
  def fingerprint(file: File): CsvFileType = {
    val fileName = file.getName.toLowerCase

    // First check by extension for binary/non-CSV formats
    if (fileName.endsWith(".bam") || fileName.endsWith(".cram")) {
      return CsvFileType.Alignment
    }
    if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
      return CsvFileType.VcfVariants
    }

    // For text files, analyze content
    readSampleLines(file) match {
      case Right(lines) if lines.nonEmpty =>
        fingerprintFromContent(lines, fileName)
      case Right(_) =>
        log.warn(s"Empty file: ${file.getName}")
        CsvFileType.Unknown
      case Left(error) =>
        log.error(s"Failed to read file ${file.getName}: $error")
        CsvFileType.Unknown
    }
  }

  /**
   * Reads the first N lines from a file (supports gzip).
   */
  private def readSampleLines(file: File): Either[String, List[String]] = {
    Try {
      val fileName = file.getName.toLowerCase
      val inputStream = if (fileName.endsWith(".gz")) {
        new GZIPInputStream(new FileInputStream(file))
      } else {
        new FileInputStream(file)
      }

      Using.resource(new BufferedReader(new InputStreamReader(inputStream, "UTF-8"))) { reader =>
        Iterator.continually(reader.readLine())
          .takeWhile(_ != null)
          .take(SampleLines)
          .toList
      }
    }.toEither.left.map(_.getMessage)
  }

  /**
   * Analyzes content lines to determine file type.
   */
  private def fingerprintFromContent(lines: List[String], fileName: String): CsvFileType = {
    // Skip comment lines for analysis but keep them for vendor detection
    val commentLines = lines.filter(l => l.startsWith("#") || l.startsWith("\"#"))
    val dataLines = lines.filterNot(l => l.startsWith("#") || l.startsWith("\"#") || l.trim.isEmpty)

    if (dataLines.isEmpty) {
      return CsvFileType.Unknown
    }

    // Check for STR profile markers
    val strScore = calculateStrScore(dataLines, fileName)
    val chipScore = calculateChipScore(dataLines, commentLines, fileName)

    log.debug(s"Fingerprint scores for ${fileName}: STR=$strScore, Chip=$chipScore")

    if (strScore > chipScore && strScore >= 3) {
      CsvFileType.StrProfile
    } else if (chipScore > strScore && chipScore >= 3) {
      val vendor = detectChipVendor(dataLines, commentLines, fileName)
      CsvFileType.ChipData(vendor)
    } else if (strScore >= 2) {
      // Lower threshold for STR if it's the only match
      CsvFileType.StrProfile
    } else if (chipScore >= 2) {
      val vendor = detectChipVendor(dataLines, commentLines, fileName)
      CsvFileType.ChipData(vendor)
    } else {
      CsvFileType.Unknown
    }
  }

  /**
   * Calculates a score indicating likelihood of STR profile data.
   */
  private def calculateStrScore(lines: List[String], fileName: String): Int = {
    var score = 0

    // Check filename hints
    val lowerFileName = fileName.toLowerCase
    if (lowerFileName.contains("str") || lowerFileName.contains("ystr")) score += 2
    if (lowerFileName.contains("ftdna") || lowerFileName.contains("yseq")) score += 1

    // Analyze content
    val allContent = lines.mkString("\n").toUpperCase

    // Strong indicators: DYS marker names
    val dysMarkerPattern = """DYS\d{2,3}""".r
    val dysMatches = dysMarkerPattern.findAllIn(allContent).toList.distinct
    if (dysMatches.size >= 10) score += 4
    else if (dysMatches.size >= 5) score += 3
    else if (dysMatches.size >= 2) score += 2
    else if (dysMatches.nonEmpty) score += 1

    // Other Y-STR marker patterns
    val otherStrMarkers = Set("DYF", "GATA", "YCAII", "CDY", "Y-GATA", "Y-GGAAT")
    val hasOtherMarkers = otherStrMarkers.exists(m => allContent.contains(m))
    if (hasOtherMarkers) score += 2

    // Check for typical STR value patterns (small integers 8-40)
    val firstDataLine = lines.headOption.getOrElse("")
    val cols = splitCsvLine(firstDataLine)
    if (cols.length == 2) {
      // Two-column format is typical for STR
      val possibleValue = cols.lift(1).map(_.trim)
      if (possibleValue.exists(v => v.matches("\\d{1,2}") || v.matches("\\d{1,2}-\\d{1,2}"))) {
        score += 1
      }
    }

    // Horizontal format check - many DYS columns
    if (cols.length > 20 && cols.count(c => c.toUpperCase.startsWith("DYS")) > 10) {
      score += 3
    }

    score
  }

  /**
   * Calculates a score indicating likelihood of chip/SNP data.
   */
  private def calculateChipScore(lines: List[String], commentLines: List[String], fileName: String): Int = {
    var score = 0

    // Check filename hints
    val lowerFileName = fileName.toLowerCase
    if (lowerFileName.contains("23andme")) score += 3
    if (lowerFileName.contains("ancestry")) score += 3
    if (lowerFileName.contains("myheritage")) score += 3
    if (lowerFileName.contains("livingdna")) score += 3
    if (lowerFileName.contains("ftdna") && lowerFileName.contains("raw")) score += 2
    if (lowerFileName.contains("snp") || lowerFileName.contains("chip") || lowerFileName.contains("array")) score += 1
    if (lowerFileName.contains("genome")) score += 1

    // Check comment lines for vendor signatures
    val commentContent = commentLines.mkString("\n").toLowerCase
    if (commentContent.contains("23andme")) score += 3
    if (commentContent.contains("ancestrydna") || commentContent.contains("ancestry dna")) score += 3
    if (commentContent.contains("myheritage")) score += 3
    if (commentContent.contains("living dna")) score += 3

    // Analyze content
    val allContent = lines.take(20).mkString("\n")
    val upperContent = allContent.toUpperCase

    // Strong indicators: rsID patterns
    val rsidPattern = """rs\d{4,}""".r
    val rsidMatches = rsidPattern.findAllIn(allContent.toLowerCase).toList.distinct
    if (rsidMatches.size >= 10) score += 4
    else if (rsidMatches.size >= 5) score += 3
    else if (rsidMatches.size >= 2) score += 2
    else if (rsidMatches.nonEmpty) score += 1

    // Check for chromosome column headers
    val headerIndicators = Set("CHROMOSOME", "CHROM", "CHR", "POSITION", "POS", "GENOTYPE", "ALLELE1", "ALLELE2", "RSID")
    val firstLine = lines.headOption.getOrElse("").toUpperCase
    val headerMatches = headerIndicators.count(h => firstLine.contains(h))
    if (headerMatches >= 3) score += 3
    else if (headerMatches >= 2) score += 2
    else if (headerMatches >= 1) score += 1

    // Check for genotype patterns (AA, AG, CC, TT, etc.)
    val genotypePattern = """[ACGT]{2}""".r
    val dataLines = lines.filterNot(_.startsWith("#"))
    if (dataLines.nonEmpty) {
      val sampleLine = dataLines.head
      val cols = splitCsvLine(sampleLine)
      // Chip data usually has 4-5 columns with last being genotype
      if (cols.length >= 4 && cols.length <= 6) {
        val lastCol = cols.last.trim.toUpperCase
        if (lastCol.matches("[ACGT]{2}") || lastCol.matches("[ACGT]/[ACGT]") || lastCol == "--" || lastCol == "NC") {
          score += 2
        }
      }
    }

    // Check for chromosome values (1-22, X, Y, MT)
    val chromPattern = """^(chr)?(1[0-9]|2[0-2]|[1-9]|X|Y|MT|M)$""".r
    val hasChromValues = lines.exists { line =>
      val cols = splitCsvLine(line)
      cols.exists(c => chromPattern.findFirstIn(c.trim.toUpperCase).isDefined)
    }
    if (hasChromValues) score += 2

    score
  }

  /**
   * Detects the chip vendor from content patterns.
   */
  private def detectChipVendor(lines: List[String], commentLines: List[String], fileName: String): Option[String] = {
    val lowerFileName = fileName.toLowerCase
    val commentContent = commentLines.mkString("\n").toLowerCase
    val firstDataLine = lines.filterNot(_.startsWith("#")).headOption.getOrElse("")

    // Check filename first
    if (lowerFileName.contains("23andme")) return Some("23andMe")
    if (lowerFileName.contains("ancestry")) return Some("AncestryDNA")
    if (lowerFileName.contains("myheritage")) return Some("MyHeritage")
    if (lowerFileName.contains("livingdna") || lowerFileName.contains("living_dna")) return Some("LivingDNA")
    if (lowerFileName.contains("ftdna") || lowerFileName.contains("familytree")) return Some("FTDNA")

    // Check comments
    if (commentContent.contains("23andme")) return Some("23andMe")
    if (commentContent.contains("ancestrydna") || commentContent.contains("ancestry dna")) return Some("AncestryDNA")
    if (commentContent.contains("myheritage")) return Some("MyHeritage")
    if (commentContent.contains("living dna") || commentContent.contains("livingdna")) return Some("LivingDNA")

    // Check format patterns
    val cols = splitCsvLine(firstDataLine)

    // 23andMe: rsid, chromosome, position, genotype (4 columns, tab-separated often)
    if (cols.length == 4) {
      val lastCol = cols.last.trim.toUpperCase
      if (lastCol.matches("[ACGT]{2}") || lastCol == "--") {
        return Some("23andMe")
      }
    }

    // AncestryDNA: rsid, chromosome, position, allele1, allele2 (5 columns)
    if (cols.length == 5) {
      val allele1 = cols(3).trim.toUpperCase
      val allele2 = cols(4).trim.toUpperCase
      if ((allele1.matches("[ACGT0]") || allele1 == "-") && (allele2.matches("[ACGT0]") || allele2 == "-")) {
        return Some("AncestryDNA")
      }
    }

    None
  }

  /**
   * Splits a CSV/TSV line handling quoted values.
   */
  private def splitCsvLine(line: String): Array[String] = {
    val result = scala.collection.mutable.ArrayBuffer.empty[String]
    var current = new StringBuilder()
    var inQuotes = false

    line.foreach { c =>
      c match {
        case '"' =>
          inQuotes = !inQuotes
        case ',' | '\t' if !inQuotes =>
          result += current.toString.trim
          current = new StringBuilder()
        case _ =>
          current += c
      }
    }
    result += current.toString.trim
    result.toArray
  }

  /**
   * Quick check if a file is likely a CSV/text file (vs binary).
   */
  def isTextFile(file: File): Boolean = {
    val fileName = file.getName.toLowerCase
    fileName.endsWith(".csv") ||
      fileName.endsWith(".tsv") ||
      fileName.endsWith(".txt") ||
      fileName.endsWith(".csv.gz") ||
      fileName.endsWith(".tsv.gz") ||
      fileName.endsWith(".txt.gz")
  }
}
