package com.decodingus.analysis

import com.decodingus.util.Logger
import org.broadinstitute.hellbender.Main

import java.io.{ByteArrayOutputStream, File, OutputStream, PrintStream}
import java.nio.file.{Files, Paths}

/**
 * Safely executes GATK tools and captures stdout/stderr.
 * Uses GATK's instanceMain() method which returns exit codes instead of calling System.exit().
 */
object GatkRunner {

  private val log = Logger("GatkRunner")

  case class GatkResult(exitCode: Int, stdout: String, stderr: String)

  /**
   * A PrintStream that intercepts lines and calls a callback for progress parsing.
   * Also accumulates all output for the final result.
   */
  private class ProgressCapturingStream(
    accumulator: ByteArrayOutputStream,
    onLine: String => Unit
  ) extends PrintStream(accumulator) {
    private val lineBuffer = new StringBuilder

    override def write(b: Int): Unit = {
      super.write(b)
      if (b == '\n') {
        onLine(lineBuffer.toString)
        lineBuffer.clear()
      } else {
        lineBuffer.append(b.toChar)
      }
    }

    override def write(buf: Array[Byte], off: Int, len: Int): Unit = {
      super.write(buf, off, len)
      val str = new String(buf, off, len)
      str.foreach { c =>
        if (c == '\n') {
          onLine(lineBuffer.toString)
          lineBuffer.clear()
        } else {
          lineBuffer.append(c)
        }
      }
    }
  }

  /**
   * Parses GATK/Picard ProgressLogger output to extract progress info.
   * Picard logs lines like: "INFO ... Processed 1,000,000 records"
   * GATK Walker logs lines like: "INFO ... chr1:12345678" showing current position
   *
   * @param line The log line to parse
   * @param totalRecords Optional total record count for percentage calculation
   * @param contigLengths Optional map of contig names to lengths for position-based progress
   * @return Optional (message, fraction) tuple
   */
  private def parseProgressLine(
    line: String,
    totalRecords: Option[Long] = None,
    contigLengths: Option[Map[String, Long]] = None
  ): Option[(String, Double)] = {
    // Picard-style: "Processed 1,000,000 records" or "Read 1,000,000 records"
    val recordsPattern = """(?:Processed|Read)\s+([\d,]+)\s+(?:records?|loci)""".r.unanchored
    // GATK Walker-style position: "chr1:12345678" or just position info
    val positionPattern = """(chr[\dXYM]+):(\d+)""".r.unanchored

    line match {
      case recordsPattern(countStr) =>
        val count = countStr.replace(",", "").toLong
        val fraction = totalRecords.map(total => count.toDouble / total).getOrElse(0.0)
        Some((s"Processed ${countStr} records", fraction.min(0.99)))

      case positionPattern(contig, posStr) if contigLengths.isDefined =>
        val pos = posStr.toLong
        contigLengths.flatMap { lengths =>
          lengths.get(contig).map { length =>
            val fraction = pos.toDouble / length
            (s"Processing $contig:$posStr", fraction.min(0.99))
          }
        }

      case _ => None
    }
  }

  // Cache directory for cloud file indexes
  private val IndexCacheDir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "indexes")

  /**
   * Checks if a path is a cloud/remote URL (gs://, s3://, http://, https://)
   */
  private def isCloudPath(path: String): Boolean = {
    path.startsWith("gs://") || path.startsWith("s3://") ||
    path.startsWith("http://") || path.startsWith("https://")
  }

  /**
   * Gets the index file extension for a BAM/CRAM file
   */
  private def getIndexExtension(path: String): String = {
    if (path.toLowerCase.endsWith(".cram")) ".crai" else ".bai"
  }

  /**
   * Ensures a BAM/CRAM index exists, creating one if necessary.
   * For local files: looks for/creates index next to the BAM/CRAM
   * For cloud URLs: looks for/creates index in ~/.decodingus/cache/indexes/
   *
   * @param bamPath Path to the BAM/CRAM file (local path or cloud URL)
   * @return Right(indexPath) if index exists or was created, Left(error) if creation failed
   */
  def ensureIndex(bamPath: String): Either[String, String] = {
    val isCram = bamPath.toLowerCase.endsWith(".cram")
    val indexExt = getIndexExtension(bamPath)

    if (isCloudPath(bamPath)) {
      // Cloud path - check/create index in cache directory
      ensureCloudIndex(bamPath, indexExt)
    } else {
      // Local path - check/create index next to the file
      ensureLocalIndex(bamPath, indexExt)
    }
  }

  private def ensureLocalIndex(bamPath: String, indexExt: String): Either[String, String] = {
    // Possible index locations for local files
    val possibleIndexPaths = Seq(
      bamPath + indexExt,                                    // file.bam.bai
      bamPath.stripSuffix(if (indexExt == ".crai") ".cram" else ".bam") + indexExt  // file.bai
    )

    possibleIndexPaths.find(p => new File(p).exists()) match {
      case Some(existingIndex) =>
        Right(existingIndex)
      case None =>
        // Need to create the index - GATK will put it next to the input file
        log.info("BAM index not found, creating with BuildBamIndex...")
        val args = Array(
          "BuildBamIndex",
          "-I", bamPath
        )
        run(args) match {
          case Right(_) =>
            // Find the created index
            possibleIndexPaths.find(p => new File(p).exists()) match {
              case Some(idx) =>
                log.info(s"Created index: $idx")
                Right(idx)
              case None =>
                Left("BuildBamIndex completed but index file not found")
            }
          case Left(error) =>
            Left(s"Failed to create BAM index: $error")
        }
    }
  }

  private def ensureCloudIndex(bamPath: String, indexExt: String): Either[String, String] = {
    // For cloud files, first check if index exists at the cloud location
    val cloudIndexPath = bamPath + indexExt

    // Create cache directory if needed
    if (!Files.exists(IndexCacheDir)) {
      Files.createDirectories(IndexCacheDir)
    }

    // Generate a cache filename based on the cloud URL
    val cacheFileName = bamPath.hashCode.toHexString + indexExt
    val cachedIndexPath = IndexCacheDir.resolve(cacheFileName)

    if (Files.exists(cachedIndexPath)) {
      log.debug(s"Using cached index: $cachedIndexPath")
      Right(cachedIndexPath.toString)
    } else {
      // Try to create index - GATK should handle cloud URLs
      // We specify output location in the cache
      log.info("Cloud BAM index not found, creating with BuildBamIndex...")
      val args = Array(
        "BuildBamIndex",
        "-I", bamPath,
        "-O", cachedIndexPath.toString
      )
      run(args) match {
        case Right(_) =>
          if (Files.exists(cachedIndexPath)) {
            log.info(s"Created index in cache: $cachedIndexPath")
            Right(cachedIndexPath.toString)
          } else {
            Left("BuildBamIndex completed but index file not found in cache")
          }
        case Left(error) =>
          Left(s"Failed to create BAM index: $error")
      }
    }
  }

  /**
   * Runs a GATK tool with the given arguments.
   * Uses instanceMain() to avoid System.exit() calls.
   *
   * @param args Command line arguments for GATK (tool name first, then options)
   * @return Either an error message (Left) or success (Right with exit code 0)
   */
  def run(args: Array[String]): Either[String, GatkResult] = {
    runWithProgress(args, None, None, None)
  }

  /**
   * Runs a GATK tool with progress callback support.
   * Parses GATK/Picard log output to extract progress information.
   *
   * @param args Command line arguments for GATK (tool name first, then options)
   * @param onProgress Optional callback receiving (message, fractionComplete) - fractions are 0.0 to 1.0
   * @param totalRecords Optional total record count for Picard-style progress (e.g., total reads)
   * @param contigLengths Optional contig lengths for position-based progress (e.g., Map("chr1" -> 248956422L))
   * @return Either an error message (Left) or success (Right with exit code 0)
   */
  def runWithProgress(
    args: Array[String],
    onProgress: Option[(String, Double) => Unit],
    totalRecords: Option[Long] = None,
    contigLengths: Option[Map[String, Long]] = None
  ): Either[String, GatkResult] = {
    val originalOut = System.out
    val originalErr = System.err

    val stdoutCapture = new ByteArrayOutputStream()
    val stderrCapture = new ByteArrayOutputStream()

    // Line handler that parses progress and calls the callback
    val lineHandler: String => Unit = onProgress match {
      case Some(callback) => line =>
        parseProgressLine(line, totalRecords, contigLengths).foreach {
          case (msg, fraction) => callback(msg, fraction)
        }
      case None => _ => ()
    }

    try {
      // Use progress-capturing streams if callback provided, otherwise simple capture
      val outStream = onProgress match {
        case Some(_) => new ProgressCapturingStream(stdoutCapture, lineHandler)
        case None => new PrintStream(stdoutCapture)
      }
      val errStream = onProgress match {
        case Some(_) => new ProgressCapturingStream(stderrCapture, lineHandler)
        case None => new PrintStream(stderrCapture)
      }

      System.setOut(outStream)
      System.setErr(errStream)

      // Use instanceMain which returns exit code instead of calling System.exit()
      val gatkMain = new Main()
      val exitCodeObj = gatkMain.instanceMain(args)
      val exitCode = exitCodeObj match {
        case i: java.lang.Integer => i.intValue()
        case _ => 0
      }

      val stdout = stdoutCapture.toString
      val stderr = stderrCapture.toString

      if (exitCode == 0) {
        onProgress.foreach(_(s"${args.headOption.getOrElse("GATK")} complete", 1.0))
        Right(GatkResult(exitCode, stdout, stderr))
      } else {
        val toolName = args.headOption.getOrElse("Unknown")
        Left(s"$toolName failed with exit code $exitCode.\n$stderr")
      }

    } catch {
      case e: Exception =>
        val stdout = stdoutCapture.toString
        val stderr = stderrCapture.toString
        val toolName = args.headOption.getOrElse("Unknown")
        Left(s"$toolName threw an exception: ${e.getMessage}\n$stderr")
    } finally {
      // Restore original streams
      System.setOut(originalOut)
      System.setErr(originalErr)
    }
  }
}
