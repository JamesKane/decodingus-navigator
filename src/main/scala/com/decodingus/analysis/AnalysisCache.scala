package com.decodingus.analysis

import com.decodingus.model.{ContigSummary, CoverageSummary, LibraryStats, WgsMetrics}
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*
import io.circe.{Decoder, Encoder}

import java.io.{BufferedInputStream, File, FileInputStream, PrintWriter}
import java.nio.file.{Files, Paths}
import java.security.MessageDigest
import scala.util.Using

object AnalysisCache {

  private val CacheDir = Paths.get(System.getProperty("user.home"), ".decodingus", "cache")

  // Encoders/Decoders for Circe
  implicit val libraryStatsEncoder: Encoder[LibraryStats] = deriveEncoder
  implicit val libraryStatsDecoder: Decoder[LibraryStats] = deriveDecoder
  implicit val wgsMetricsEncoder: Encoder[WgsMetrics] = deriveEncoder
  implicit val wgsMetricsDecoder: Decoder[WgsMetrics] = deriveDecoder
  implicit val contigSummaryEncoder: Encoder[ContigSummary] = deriveEncoder
  implicit val contigSummaryDecoder: Decoder[ContigSummary] = deriveDecoder
  implicit val coverageSummaryEncoder: Encoder[CoverageSummary] = deriveEncoder
  implicit val coverageSummaryDecoder: Decoder[CoverageSummary] = deriveDecoder

  /**
   * Calculates the SHA-256 hash of a file.
   * This reads the entire file, so it should be run on a background thread.
   */
  def calculateSha256(file: File): String = {
    val buffer = new Array[Byte](8192)
    val digest = MessageDigest.getInstance("SHA-256")
    val bis = new BufferedInputStream(new FileInputStream(file))
    try {
      var bytesRead = bis.read(buffer)
      while (bytesRead != -1) {
        digest.update(buffer, 0, bytesRead)
        bytesRead = bis.read(buffer)
      }
    } finally {
      bis.close()
    }
    digest.digest().map("%02x".format(_)).mkString
  }

  /**
   * Saves the CoverageSummary to the cache directory with the filename <sha256>.json.
   */
  def save(sha256: String, summary: CoverageSummary): Unit = {
    if (!Files.exists(CacheDir)) {
      Files.createDirectories(CacheDir)
    }
    val cacheFile = CacheDir.resolve(s"$sha256.json").toFile
    Using(new PrintWriter(cacheFile)) { writer =>
      writer.write(summary.asJson.noSpaces)
    }
  }

  /**
   * Retrieves a CoverageSummary from the cache if it exists.
   */
  def load(sha256: String): Option[CoverageSummary] = {
    val cacheFile = CacheDir.resolve(s"$sha256.json").toFile
    if (cacheFile.exists()) {
      val source = scala.io.Source.fromFile(cacheFile)
      try {
        decode[CoverageSummary](source.mkString) match {
          case Right(summary) => Some(summary)
          case Left(error) =>
            println(s"Failed to decode cached file: $error")
            None
        }
      } finally {
        source.close()
      }
    } else {
      None
    }
  }

  def exists(sha256: String): Boolean = {
    CacheDir.resolve(s"$sha256.json").toFile.exists()
  }
}
