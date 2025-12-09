package com.decodingus.ancestry.reference

import com.decodingus.ancestry.model.{AlleleFrequencyMatrix, AncestryPanelType, PCALoadings}
import com.decodingus.config.FeatureToggles

import java.io.{BufferedInputStream, FileOutputStream}
import java.net.{HttpURLConnection, URI}
import java.nio.file.{Files, Path}
import java.util.zip.GZIPInputStream
import scala.util.{Failure, Success, Try, Using}

/**
 * Result of checking ancestry reference data availability.
 */
sealed trait AncestryReferenceResult
object AncestryReferenceResult {
  /** Reference data is available locally */
  case class Available(
    sitesVcf: Path,
    alleleFreqs: AlleleFrequencyMatrix,
    pcaLoadings: PCALoadings
  ) extends AncestryReferenceResult

  /** Reference data needs to be downloaded */
  case class DownloadRequired(
    panelType: AncestryPanelType,
    referenceBuild: String,
    estimatedSizeMB: Int,
    downloadUrl: String
  ) extends AncestryReferenceResult

  /** Error accessing reference data */
  case class Error(message: String) extends AncestryReferenceResult
}

/**
 * Manages ancestry reference data download and caching.
 * Follows the same pattern as ReferenceGateway for genome references.
 */
class AncestryReferenceGateway(onProgress: (Long, Long) => Unit) {

  private val cache = new AncestryReferenceCache()

  // Base URL for reference data downloads
  private val baseUrl = "https://reference-data.decodingus.org/ancestry"

  // Estimated download sizes in MB
  private val panelSizes: Map[AncestryPanelType, Int] = Map(
    AncestryPanelType.Aims -> 15,
    AncestryPanelType.GenomeWide -> 250
  )

  /**
   * Check if ancestry reference data is available for a panel and reference build.
   * Does not download - use downloadAndResolve for that.
   */
  def checkAvailability(
    panelType: AncestryPanelType,
    referenceBuild: String
  ): AncestryReferenceResult = {
    if (cache.isPanelAvailable(panelType, referenceBuild)) {
      loadCachedPanel(panelType, referenceBuild)
    } else {
      val version = FeatureToggles.ancestryAnalysis.referenceVersion
      val panelName = panelType match {
        case AncestryPanelType.Aims => "aims"
        case AncestryPanelType.GenomeWide => "genome-wide"
      }
      AncestryReferenceResult.DownloadRequired(
        panelType = panelType,
        referenceBuild = referenceBuild,
        estimatedSizeMB = panelSizes.getOrElse(panelType, 100),
        downloadUrl = s"$baseUrl/$version/$panelName/$referenceBuild"
      )
    }
  }

  /**
   * Load cached panel data.
   */
  private def loadCachedPanel(
    panelType: AncestryPanelType,
    referenceBuild: String
  ): AncestryReferenceResult = {
    val sitesVcf = cache.getSitesVcfPath(panelType, referenceBuild)
    val alleleFreqPath = cache.getAlleleFreqPath(panelType)
    val pcaLoadingsPath = cache.getPcaLoadingsPath(panelType)

    (for {
      alleleFreqs <- AlleleFrequencyMatrix.load(alleleFreqPath)
      pcaLoadings <- PCALoadings.load(pcaLoadingsPath)
    } yield {
      AncestryReferenceResult.Available(sitesVcf, alleleFreqs, pcaLoadings)
    }).getOrElse {
      AncestryReferenceResult.Error("Failed to load cached panel data")
    }
  }

  /**
   * Resolve ancestry reference data, downloading if necessary.
   * Returns path to sites VCF and loaded reference matrices.
   */
  def resolve(
    panelType: AncestryPanelType,
    referenceBuild: String
  ): Either[String, (Path, AlleleFrequencyMatrix, PCALoadings)] = {
    checkAvailability(panelType, referenceBuild) match {
      case AncestryReferenceResult.Available(sitesVcf, alleleFreqs, pcaLoadings) =>
        Right((sitesVcf, alleleFreqs, pcaLoadings))

      case AncestryReferenceResult.DownloadRequired(_, _, _, _) =>
        Left("Ancestry reference data not available. Download required.")

      case AncestryReferenceResult.Error(msg) =>
        Left(msg)
    }
  }

  /**
   * Download and resolve ancestry reference data.
   * Call this after user confirms download.
   */
  def downloadAndResolve(
    panelType: AncestryPanelType,
    referenceBuild: String
  ): Either[String, (Path, AlleleFrequencyMatrix, PCALoadings)] = {
    val version = FeatureToggles.ancestryAnalysis.referenceVersion
    val panelName = panelType match {
      case AncestryPanelType.Aims => "aims"
      case AncestryPanelType.GenomeWide => "genome-wide"
    }

    // Download each required file
    val filesToDownload = List(
      (s"${referenceBuild}_sites.vcf.gz", cache.getSitesVcfPath(panelType, referenceBuild)),
      (s"${referenceBuild}_sites.vcf.gz.tbi", cache.getSitesVcfIndexPath(panelType, referenceBuild)),
      ("allele_freqs.bin", cache.getAlleleFreqPath(panelType)),
      ("pca_loadings.bin", cache.getPcaLoadingsPath(panelType))
    )

    val downloadResults = filesToDownload.map { case (filename, localPath) =>
      if (Files.exists(localPath)) {
        Right(localPath)
      } else {
        val url = s"$baseUrl/$version/$panelName/$filename"
        downloadFile(url, localPath)
      }
    }

    // Check for any download failures
    val failures = downloadResults.collect { case Left(err) => err }
    if (failures.nonEmpty) {
      Left(s"Download failed: ${failures.mkString(", ")}")
    } else {
      resolve(panelType, referenceBuild)
    }
  }

  /**
   * Download a file from URL to local path.
   */
  private def downloadFile(url: String, localPath: Path): Either[String, Path] = {
    Try {
      Files.createDirectories(localPath.getParent)
      val connection = URI.create(url).toURL.openConnection().asInstanceOf[HttpURLConnection]
      connection.setRequestMethod("GET")
      connection.setConnectTimeout(30000)
      connection.setReadTimeout(60000)

      val responseCode = connection.getResponseCode
      if (responseCode != 200) {
        throw new RuntimeException(s"HTTP $responseCode: ${connection.getResponseMessage}")
      }

      val contentLength = connection.getContentLengthLong
      var downloaded = 0L

      Using.resources(
        new BufferedInputStream(connection.getInputStream),
        new FileOutputStream(localPath.toFile)
      ) { (in, out) =>
        val buffer = new Array[Byte](8192)
        var bytesRead = 0
        while ({bytesRead = in.read(buffer); bytesRead != -1}) {
          out.write(buffer, 0, bytesRead)
          downloaded += bytesRead
          if (contentLength > 0) {
            onProgress(downloaded, contentLength)
          }
        }
      }

      localPath
    } match {
      case Success(path) => Right(path)
      case Failure(e) => Left(s"Download failed for $url: ${e.getMessage}")
    }
  }

  /**
   * Get estimated download size for a panel.
   */
  def getEstimatedSize(panelType: AncestryPanelType): Int = {
    panelSizes.getOrElse(panelType, 100)
  }

  /**
   * Check which builds have cached data for a panel.
   */
  def listCachedBuilds(panelType: AncestryPanelType): List[String] = {
    val knownBuilds = List("GRCh38", "GRCh37", "CHM13v2")
    knownBuilds.filter(build => cache.isPanelAvailable(panelType, build))
  }
}
