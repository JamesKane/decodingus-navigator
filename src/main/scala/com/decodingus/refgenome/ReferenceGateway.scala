package com.decodingus.refgenome

import com.decodingus.analysis.GatkRunner
import com.decodingus.refgenome.config.{ReferenceConfig, ReferenceConfigService}
import sttp.client3.*

import java.io.IOException
import java.nio.file.{Files, Path, Paths}
import scala.sys.process.*

// NOTE: Do NOT import org.broadinstitute.hellbender.Main here!
// Use GatkRunner which handles Log4j initialization properly.

/**
 * Result type for reference resolution when user confirmation is needed.
 */
sealed trait ReferenceResolveResult

object ReferenceResolveResult {
  /** Reference was found locally and is ready to use */
  case class Available(path: Path) extends ReferenceResolveResult

  /** Reference not found, download is required - includes estimated size info */
  case class DownloadRequired(build: String, url: String, estimatedSizeMB: Int) extends ReferenceResolveResult

  /** Error occurred */
  case class Error(message: String) extends ReferenceResolveResult
}

class ReferenceGateway(onProgress: (Long, Long) => Unit) {
  private val cache = new ReferenceCache

  private val referenceUrls: Map[String, String] = ReferenceConfig.knownBuilds

  // Estimated download sizes in MB for user information
  private val estimatedSizes: Map[String, Int] = Map(
    "GRCh38" -> 3100, // ~3.1 GB uncompressed, downloads as uncompressed
    "GRCh37" -> 900, // ~900 MB compressed
    "CHM13v2" -> 900 // ~900 MB compressed
  )

  /**
   * Resolves a reference, automatically downloading if config allows.
   * This is the original behavior for backward compatibility.
   */
  def resolve(referenceBuild: String): Either[String, Path] = {
    cache.getPath(referenceBuild) match {
      case Some(path) =>
        println(s"Found reference $referenceBuild in cache: $path")
        validateAndCreateReferenceFiles(path)
      case None =>
        val config = ReferenceConfigService.load()
        val buildConfig = config.getOrDefault(referenceBuild)

        // Check if auto-download is enabled for this build
        if (buildConfig.autoDownload || !config.promptBeforeDownload) {
          referenceUrls.get(referenceBuild) match {
            case Some(url) =>
              downloadReference(referenceBuild, url).flatMap(validateAndCreateReferenceFiles)
            case None => Left(s"Unknown reference build: $referenceBuild")
          }
        } else {
          Left(s"Reference $referenceBuild not found locally. Configure a local path or enable download in Settings.")
        }
    }
  }

  /**
   * Checks if a reference is available without downloading.
   * Returns a result indicating availability or download requirement.
   */
  def checkAvailability(referenceBuild: String): ReferenceResolveResult = {
    cache.getPath(referenceBuild) match {
      case Some(path) =>
        validateAndCreateReferenceFiles(path) match {
          case Right(validPath) => ReferenceResolveResult.Available(validPath)
          case Left(error) => ReferenceResolveResult.Error(error)
        }
      case None =>
        referenceUrls.get(referenceBuild) match {
          case Some(url) =>
            val sizeMB = estimatedSizes.getOrElse(referenceBuild, 1000)
            ReferenceResolveResult.DownloadRequired(referenceBuild, url, sizeMB)
          case None =>
            ReferenceResolveResult.Error(s"Unknown reference build: $referenceBuild")
        }
    }
  }

  /**
   * Downloads a reference after user confirmation.
   * Call this after checkAvailability returns DownloadRequired and user approves.
   */
  def downloadAndResolve(referenceBuild: String): Either[String, Path] = {
    referenceUrls.get(referenceBuild) match {
      case Some(url) =>
        downloadReference(referenceBuild, url).flatMap(validateAndCreateReferenceFiles)
      case None =>
        Left(s"Unknown reference build: $referenceBuild")
    }
  }

  /**
   * Gets the URL for a reference build (for display to user).
   */
  def getDownloadUrl(referenceBuild: String): Option[String] = referenceUrls.get(referenceBuild)

  /**
   * Gets the estimated download size in MB.
   */
  def getEstimatedSizeMB(referenceBuild: String): Int = estimatedSizes.getOrElse(referenceBuild, 1000)

  private def downloadReference(referenceBuild: String, url: String): Either[String, Path] = {
    println(s"Downloading reference $referenceBuild from $url")
    val tempFileRaw = Files.createTempFile(s"ref-$referenceBuild", ".tmp")
    val tempFileGzipped = Files.createTempFile(s"ref-$referenceBuild", ".fa.gz")
    Files.deleteIfExists(tempFileGzipped) // Delete the empty .fa.gz file created by createTempFile

    val request = basicRequest.get(uri"$url").response(asFile(tempFileRaw.toFile))

    val backend = HttpURLConnectionBackend()
    val response = request.send(backend)

    response.body match {
      case Right(file) =>
        println("Download complete.")
        val sourcePathForCache = if (url.endsWith(".gz")) {
          // Already gzipped, just move the raw downloaded file to the .fa.gz temp path
          Files.move(file.toPath, tempFileGzipped)
          tempFileGzipped
        } else {
          // Not gzipped, apply bgzip
          println(s"Compressing $file with bgzip...")
          val command = s"bgzip -c ${file.toPath} > ${tempFileGzipped}"
          try {
            val exitCode = command.!
            if (exitCode != 0) {
              Files.deleteIfExists(file.toPath)
              Files.deleteIfExists(tempFileGzipped)
              return Left(s"Failed to bgzip $file. Exit code: $exitCode")
            }
            Files.deleteIfExists(file.toPath) // Delete the original uncompressed temp file
            tempFileGzipped
          } catch {
            case e: IOException =>
              Files.deleteIfExists(file.toPath)
              Files.deleteIfExists(tempFileGzipped)
              return Left(s"Failed to execute bgzip for $file: ${e.getMessage}")
          }
        }
        println("Caching reference.")
        val finalPath = cache.put(referenceBuild, sourcePathForCache)
        Right(finalPath)
      case Left(error) =>
        Files.deleteIfExists(tempFileRaw)
        Files.deleteIfExists(tempFileGzipped)
        Left(s"Failed to download reference: $error")
    }
  }

  private def validateAndCreateReferenceFiles(referencePath: Path): Either[String, Path] = {
    val faiPath = Paths.get(referencePath.toString + ".fai")
    val dictPath = Paths.get(referencePath.getParent.toString, referencePath.getFileName.toString.replace(".fa.gz", ".dict"))

    // Check and create .fai index
    if (!Files.exists(faiPath)) {
      println(s"Creating FASTA index for $referencePath...")
      val command = s"samtools faidx $referencePath"
      try {
        val exitCode = command.!
        if (exitCode != 0) {
          return Left(s"Failed to create FASTA index for $referencePath. Exit code: $exitCode")
        }
      } catch {
        case e: IOException =>
          return Left(s"Failed to execute samtools faidx for $referencePath: ${e.getMessage}")
      }
    }

    // Check and create .dict dictionary
    if (!Files.exists(dictPath)) {
      println(s"Creating sequence dictionary for $referencePath using GATK...")
      val args = Array(
        "CreateSequenceDictionary",
        "-R", referencePath.toAbsolutePath.toString,
        "-O", dictPath.toAbsolutePath.toString
      )
      GatkRunner.run(args) match {
        case Right(_) => // Success
        case Left(error) =>
          return Left(s"Failed to create sequence dictionary for $referencePath: $error")
      }
    }
    Right(referencePath)
  }
}
