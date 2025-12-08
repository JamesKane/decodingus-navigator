package com.decodingus.str

import com.decodingus.workspace.model._
import java.io.File
import java.time.LocalDateTime
import scala.io.Source
import scala.util.{Try, Using}

/**
 * Parses Y-STR CSV files from various vendors (FTDNA, YSEQ, etc.)
 * into StrProfile objects.
 */
object StrCsvParser {

  /** Detected vendor format */
  sealed trait VendorFormat
  object VendorFormat {
    case object FTDNA extends VendorFormat
    case object YSEQ extends VendorFormat
    case object Generic extends VendorFormat
    case object Unknown extends VendorFormat
  }

  /** Result of parsing a CSV file */
  case class ParseResult(
    profile: StrProfile,
    detectedFormat: VendorFormat,
    warnings: List[String]
  )

  /**
   * Parses a Y-STR CSV file and returns an StrProfile.
   *
   * @param file The CSV file to parse
   * @param biosampleRef The AT URI of the parent biosample
   * @return Either an error message or the parsed StrProfile with metadata
   */
  def parse(file: File, biosampleRef: String): Either[String, ParseResult] = {
    Using(Source.fromFile(file)) { source =>
      val lines = source.getLines().toList
      if (lines.isEmpty) {
        Left("CSV file is empty")
      } else {
        parseLines(lines, file.getName, biosampleRef)
      }
    }.toEither.left.map(_.getMessage).flatten
  }

  /**
   * Parses CSV content from a string.
   */
  def parseString(content: String, fileName: String, biosampleRef: String): Either[String, ParseResult] = {
    val lines = content.split("\n").toList
    if (lines.isEmpty) {
      Left("CSV content is empty")
    } else {
      parseLines(lines, fileName, biosampleRef)
    }
  }

  private def parseLines(lines: List[String], fileName: String, biosampleRef: String): Either[String, ParseResult] = {
    val (format, headerIdx) = detectFormat(lines)
    format match {
      case VendorFormat.Unknown =>
        Left("Could not detect CSV format. Expected FTDNA, YSEQ, or a two-column Marker,Value format.")
      case _ =>
        parseWithFormat(lines, headerIdx, format, fileName, biosampleRef)
    }
  }

  /**
   * Detects the vendor format from CSV headers.
   * Returns the format and the index of the header row.
   */
  private def detectFormat(lines: List[String]): (VendorFormat, Int) = {
    // FTDNA format typically has: "Marker Name" or "DYS#" in first column, "Allele" or value in second
    // Sometimes has a kit number header row first
    // YSEQ format: might have "Marker", "Value" or similar headers

    val detected = lines.zipWithIndex.collectFirst {
      case (line, idx) =>
        val lower = line.toLowerCase.trim
        val cols = splitCsvLine(line).map(_.toLowerCase.trim)

        // FTDNA format detection
        if (lower.contains("marker name") && lower.contains("allele")) {
          Some((VendorFormat.FTDNA, idx))
        } else if (cols.headOption.exists(c => c == "marker name" || c == "marker") &&
                   cols.lift(1).exists(c => c == "allele" || c == "value" || c == "alleles")) {
          // Check if it looks like FTDNA based on other columns
          if (cols.exists(c => c.contains("ftdna") || c.contains("ystr"))) {
            Some((VendorFormat.FTDNA, idx))
          } else {
            Some((VendorFormat.Generic, idx))
          }
        } else if (lower.contains("yseq") || cols.exists(_.contains("yseq"))) {
          // YSEQ format detection - often has specific headers
          Some((VendorFormat.YSEQ, idx))
        } else if (cols.size >= 2 && idx < 3) {
          // Check for two-column format (Marker, Value)
          val firstCol = cols.head
          if (firstCol.startsWith("dys") || firstCol.startsWith("marker") || firstCol == "name") {
            Some((VendorFormat.Generic, idx))
          } else {
            None
          }
        } else {
          None
        }
    }.flatten

    detected.getOrElse {
      // If we have at least 2 columns and first row looks like data (DYS marker)
      if (lines.nonEmpty) {
        val firstLine = splitCsvLine(lines.head)
        if (firstLine.size >= 2 && firstLine.head.toUpperCase.startsWith("DYS")) {
          (VendorFormat.Generic, -1) // No header, data starts at 0
        } else {
          (VendorFormat.Unknown, -1)
        }
      } else {
        (VendorFormat.Unknown, -1)
      }
    }
  }

  private def parseWithFormat(
    lines: List[String],
    headerIdx: Int,
    format: VendorFormat,
    fileName: String,
    biosampleRef: String
  ): Either[String, ParseResult] = {
    val dataLines = if (headerIdx >= 0) lines.drop(headerIdx + 1) else lines
    var warnings = List.empty[String]
    var markers = List.empty[StrMarkerValue]
    var panels = Set.empty[String]

    dataLines.foreach { line =>
      if (line.trim.nonEmpty) {
        parseMarkerLine(line, format) match {
          case Right(Some(marker)) =>
            markers = markers :+ marker
            marker.panel.foreach(p => panels += p)
          case Right(None) =>
            // Empty or skipped line
          case Left(warning) =>
            warnings = warnings :+ warning
        }
      }
    }

    if (markers.isEmpty) {
      Left("No valid STR markers found in the file")
    } else {
      // Detect provider from format
      val provider = format match {
        case VendorFormat.FTDNA => Some("FTDNA")
        case VendorFormat.YSEQ => Some("YSEQ")
        case VendorFormat.Generic => None
        case VendorFormat.Unknown => None
      }

      // Infer panel from marker count
      val inferredPanel = inferPanel(markers.size, provider)
      val strPanels = List(StrPanel(
        panelName = inferredPanel,
        markerCount = markers.size,
        provider = provider,
        testDate = Some(LocalDateTime.now())
      ))

      val fileInfo = FileInfo(
        fileName = fileName,
        fileSizeBytes = None,
        fileFormat = "CSV",
        checksum = None,
        checksumAlgorithm = None,
        location = None
      )

      val profile = StrProfile(
        atUri = None, // Will be set when saved
        meta = RecordMeta.initial,
        biosampleRef = biosampleRef,
        sequenceRunRef = None,
        panels = strPanels,
        markers = markers,
        totalMarkers = Some(markers.size),
        source = Some("IMPORTED"),
        importedFrom = provider,
        derivationMethod = None,
        files = List(fileInfo)
      )

      Right(ParseResult(profile, format, warnings))
    }
  }

  /**
   * Parses a single marker line from the CSV.
   * Returns Right(Some(marker)) if valid, Right(None) if skipped, Left(warning) if error.
   */
  private def parseMarkerLine(line: String, format: VendorFormat): Either[String, Option[StrMarkerValue]] = {
    val cols = splitCsvLine(line)
    if (cols.size < 2) {
      Right(None) // Skip empty or malformed lines
    } else {
      val markerName = cols.head.trim
      val valueStr = cols(1).trim

      // Skip non-marker rows (like section headers)
      if (!looksLikeMarkerName(markerName)) {
        Right(None)
      } else if (valueStr.isEmpty || valueStr == "-" || valueStr.toLowerCase == "null" || valueStr.toLowerCase == "n/a") {
        Right(None) // No value for this marker
      } else {
        parseStrValue(markerName, valueStr) match {
          case Right(value) =>
            val panel = inferMarkerPanel(markerName)
            Right(Some(StrMarkerValue(
              marker = normalizeMarkerName(markerName),
              value = value,
              panel = panel,
              quality = None,
              readDepth = None
            )))
          case Left(err) =>
            Left(s"Warning: Could not parse marker $markerName='$valueStr': $err")
        }
      }
    }
  }

  /**
   * Parses a string value into an StrValue.
   * Handles:
   * - Simple integers: "13" -> SimpleStrValue(13)
   * - Multi-copy: "11-14" or "11,14" -> MultiCopyStrValue(List(11, 14))
   * - Complex: "22t-25c-26.1t" -> ComplexStrValue(...)
   */
  private def parseStrValue(markerName: String, valueStr: String): Either[String, StrValue] = {
    val trimmed = valueStr.trim

    // Check for complex notation (contains letters like 't', 'c', 'q')
    if (trimmed.matches(".*[tcq].*") && trimmed.contains("-")) {
      parseComplexValue(trimmed)
    }
    // Check for multi-copy (contains dash or comma between numbers)
    else if (trimmed.contains("-") || trimmed.contains(",")) {
      parseMultiCopyValue(trimmed)
    }
    // Simple integer value
    else {
      Try(trimmed.toInt).toEither
        .left.map(_ => s"Not a valid integer: $trimmed")
        .map(SimpleStrValue.apply)
    }
  }

  /**
   * Parses multi-copy values like "11-14" or "11,14" into MultiCopyStrValue.
   */
  private def parseMultiCopyValue(valueStr: String): Either[String, MultiCopyStrValue] = {
    val separator = if (valueStr.contains(",")) "," else "-"
    val parts = valueStr.split(separator).map(_.trim)

    val values = parts.flatMap { p =>
      Try(p.toInt).toOption
    }.toList

    if (values.isEmpty) {
      Left(s"No valid integers found in: $valueStr")
    } else {
      Right(MultiCopyStrValue(values.sorted))
    }
  }

  /**
   * Parses complex multi-allelic values like "22t-25c-26.1t" into ComplexStrValue.
   */
  private def parseComplexValue(valueStr: String): Either[String, ComplexStrValue] = {
    // Pattern: number (with optional decimal) followed by optional letter(s)
    val allelePattern = """(\d+\.?\d*)([tcqTCQ]*)""".r
    val parts = valueStr.split("-").map(_.trim)

    val alleles = parts.flatMap { part =>
      allelePattern.findFirstMatchIn(part).map { m =>
        val repeats = m.group(1).toDouble
        val designation = Option(m.group(2)).filter(_.nonEmpty).map(_.toLowerCase)
        val count = designation match {
          case Some("c") => 2 // cis = both copies
          case Some("t") => 1 // trans = one copy
          case Some("q") => 4 // quad
          case _ => 1
        }
        StrAllele(repeats, count, designation)
      }
    }.toList

    if (alleles.isEmpty) {
      Left(s"Could not parse complex value: $valueStr")
    } else {
      Right(ComplexStrValue(alleles, Some(valueStr)))
    }
  }

  /** Checks if a string looks like a valid marker name */
  private def looksLikeMarkerName(name: String): Boolean = {
    val upper = name.toUpperCase.trim
    upper.startsWith("DYS") ||
    upper.startsWith("DYF") ||
    upper.startsWith("GATA") ||
    upper.startsWith("Y-GATA") ||
    upper.startsWith("YCAII") ||
    upper.startsWith("CDY") ||
    upper == "H4" ||
    upper.matches("Y[A-Z]+.*")
  }

  /** Normalizes marker names to a standard format */
  private def normalizeMarkerName(name: String): String = {
    // Remove quotes, trim, and standardize case
    val cleaned = name.replace("\"", "").trim

    // Handle common variations
    cleaned.toUpperCase match {
      case n if n.startsWith("DYS") => n
      case n if n.startsWith("DYF") => n
      case n if n.startsWith("GATA") => n
      case n if n.startsWith("Y-GATA") => n.replace("Y-GATA", "YGATA")
      case n if n.startsWith("YCAII") => n
      case n if n.startsWith("CDY") => n
      case n => n
    }
  }

  /** Infers which panel a marker belongs to based on its name */
  private def inferMarkerPanel(markerName: String): Option[String] = {
    val upper = markerName.toUpperCase
    // Core Y-12 markers (most common)
    val y12Markers = Set(
      "DYS393", "DYS390", "DYS19", "DYS391", "DYS385A", "DYS385B",
      "DYS426", "DYS388", "DYS439", "DYS389I", "DYS392", "DYS389II"
    )

    if (y12Markers.exists(m => upper.startsWith(m.replace("I", "").replace("II", "")))) {
      Some("Y12")
    } else {
      None // Most markers beyond Y12 vary by vendor
    }
  }

  /** Infers the panel name from marker count and provider */
  private def inferPanel(markerCount: Int, provider: Option[String]): String = {
    provider match {
      case Some("FTDNA") =>
        markerCount match {
          case n if n <= 12 => "Y12"
          case n if n <= 25 => "Y25"
          case n if n <= 37 => "Y37"
          case n if n <= 67 => "Y67"
          case n if n <= 111 => "Y111"
          case n if n <= 500 => "Y500"
          case _ => "Y700"
        }
      case Some("YSEQ") =>
        if (markerCount > 400) "YSEQ_ALPHA" else "CUSTOM"
      case _ =>
        s"CUSTOM_$markerCount"
    }
  }

  /** Splits a CSV line handling quoted values */
  private def splitCsvLine(line: String): Array[String] = {
    // Simple CSV parsing - handles quoted fields with commas
    val result = scala.collection.mutable.ArrayBuffer.empty[String]
    var current = new StringBuilder()
    var inQuotes = false

    line.foreach { c =>
      c match {
        case '"' =>
          inQuotes = !inQuotes
        case ',' if !inQuotes =>
          result += current.toString.trim
          current = new StringBuilder()
        case '\t' if !inQuotes =>
          // Also support tab-separated
          result += current.toString.trim
          current = new StringBuilder()
        case _ =>
          current += c
      }
    }
    result += current.toString.trim
    result.toArray
  }
}
