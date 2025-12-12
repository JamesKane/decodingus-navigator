package com.decodingus.refgenome

import htsjdk.samtools.liftover.LiftOver
import htsjdk.samtools.util.{Interval, Log}
import sttp.client3.*

import java.io.{BufferedInputStream, FileInputStream, FileOutputStream, PrintWriter}
import java.nio.file.{Files, Path}
import java.util.zip.GZIPInputStream
import scala.io.Source
import scala.util.Using

/**
 * Paths to all Y chromosome region files for a reference build.
 */
case class YRegionPaths(
  cytobands: Path,
  palindromes: Path,
  strs: Path,
  par: Path,
  xtr: Path,
  ampliconic: Path,
  centromeres: Option[Path] = None
)

/**
 * Gateway for downloading and caching Y chromosome region annotation files.
 *
 * Downloads from:
 * - ybrowse.org: cytobands, palindromes, STRs (GFF3 format, GRCh38)
 * - GIAB genome-stratifications: PAR, XTR, ampliconic (BED format, GRCh38)
 *
 * Supports liftover to GRCh37 and CHM13v2 using htsjdk LiftOver.
 *
 * @param onProgress Callback for progress updates (message, 0.0-1.0)
 */
class YRegionGateway(onProgress: (String, Double) => Unit = (_, _) => ()) {
  private val cache = new YRegionCache
  private val liftoverGateway = new LiftoverGateway((_, _) => ())

  // ybrowse.org GFF3 URLs (GRCh38 coordinates)
  private val ybrowseUrls: Map[String, String] = Map(
    "cytobands" -> "https://ybrowse.org/gbrowse2/gff/cytobands_hg38.gff3",
    "palindromes" -> "https://ybrowse.org/gbrowse2/gff/palindromes_hg38.gff3",
    "strs" -> "https://ybrowse.org/gbrowse2/gff/str_hg38.gff3"
  )

  // GIAB genome-stratifications BED URLs (GRCh38)
  // From: https://github.com/genome-in-a-bottle/genome-stratifications
  private val giabUrls: Map[String, String] = Map(
    "par" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/XY/GRCh38_chrY_PAR.bed",
    "xtr" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/XY/GRCh38_chrY_XTR.bed",
    "ampliconic" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/XY/GRCh38_chrY_ampliconic.bed"
  )

  // Builds that need liftover from GRCh38
  private val liftoverBuilds: Set[String] = Set("GRCh37", "CHM13v2")

  /**
   * Resolve all region files for a reference build.
   *
   * @param referenceBuild Target reference build (GRCh38, GRCh37, CHM13v2)
   * @return Either error message or paths to all region files
   */
  def resolveAll(referenceBuild: String): Either[String, YRegionPaths] = {
    for {
      cytobands <- resolve("cytobands", referenceBuild)
      palindromes <- resolve("palindromes", referenceBuild)
      strs <- resolve("strs", referenceBuild)
      par <- resolve("par", referenceBuild)
      xtr <- resolve("xtr", referenceBuild)
      ampliconic <- resolve("ampliconic", referenceBuild)
    } yield YRegionPaths(cytobands, palindromes, strs, par, xtr, ampliconic)
  }

  /**
   * Resolve a single region file for a reference build.
   *
   * @param regionType Type of region (cytobands, palindromes, strs, par, xtr, ampliconic)
   * @param referenceBuild Target reference build
   * @return Either error message or path to the region file
   */
  def resolve(regionType: String, referenceBuild: String): Either[String, Path] = {
    // Check cache first
    cache.getPath(regionType, referenceBuild) match {
      case Some(path) =>
        println(s"[YRegionGateway] Found $regionType for $referenceBuild in cache: $path")
        Right(path)
      case None =>
        if (referenceBuild == "GRCh38") {
          // Download directly for GRCh38
          downloadRegionFile(regionType, referenceBuild)
        } else if (liftoverBuilds.contains(referenceBuild)) {
          // First ensure we have GRCh38, then liftover
          resolve(regionType, "GRCh38").flatMap { grch38Path =>
            liftoverRegionFile(grch38Path, regionType, "GRCh38", referenceBuild)
          }
        } else {
          Left(s"Unsupported reference build: $referenceBuild")
        }
    }
  }

  /**
   * Download a region file from the appropriate source.
   */
  private def downloadRegionFile(regionType: String, referenceBuild: String): Either[String, Path] = {
    val url = ybrowseUrls.get(regionType).orElse(giabUrls.get(regionType))
    url match {
      case Some(u) => downloadFile(regionType, referenceBuild, u)
      case None => Left(s"Unknown region type: $regionType")
    }
  }

  /**
   * Download a file from URL and cache it.
   */
  private def downloadFile(regionType: String, referenceBuild: String, url: String): Either[String, Path] = {
    onProgress(s"Downloading $regionType for $referenceBuild...", 0.0)
    println(s"[YRegionGateway] Downloading $regionType from $url")

    val ext = if (YRegionCache.gff3Types.contains(regionType)) ".gff3" else ".bed"
    val tempFile = Files.createTempFile(s"yregion-$regionType-", ext)

    try {
      val backend = HttpURLConnectionBackend()

      // Handle both plain and gzipped files
      val isGzipped = url.endsWith(".gz")
      val request = basicRequest.get(uri"$url")
        .response(asFile(tempFile.toFile))
        .readTimeout(scala.concurrent.duration.Duration(60, "seconds"))

      val response = request.send(backend)

      response.body match {
        case Right(_) =>
          val finalPath = if (isGzipped) {
            onProgress("Decompressing...", 0.5)
            val decompressed = decompressGzip(tempFile, ext)
            Files.deleteIfExists(tempFile)
            decompressed
          } else {
            tempFile
          }

          val cachedPath = cache.put(regionType, referenceBuild, finalPath)
          onProgress(s"$regionType ready.", 1.0)
          println(s"[YRegionGateway] Cached $regionType at $cachedPath")
          Right(cachedPath)

        case Left(error) =>
          Files.deleteIfExists(tempFile)
          Left(s"Failed to download $regionType: $error")
      }
    } catch {
      case e: Exception =>
        Files.deleteIfExists(tempFile)
        Left(s"Error downloading $regionType: ${e.getMessage}")
    }
  }

  /**
   * Decompress a gzipped file.
   */
  private def decompressGzip(gzPath: Path, ext: String): Path = {
    val outputPath = Files.createTempFile("yregion-decompressed-", ext)
    Using.resources(
      new GZIPInputStream(new BufferedInputStream(new FileInputStream(gzPath.toFile))),
      new FileOutputStream(outputPath.toFile)
    ) { (gzIn, out) =>
      val buffer = new Array[Byte](8192)
      var len = gzIn.read(buffer)
      while (len > 0) {
        out.write(buffer, 0, len)
        len = gzIn.read(buffer)
      }
    }
    outputPath
  }

  /**
   * Liftover a region file from one build to another.
   */
  private def liftoverRegionFile(
    sourcePath: Path,
    regionType: String,
    fromBuild: String,
    toBuild: String
  ): Either[String, Path] = {
    onProgress(s"Lifting over $regionType from $fromBuild to $toBuild...", 0.0)
    println(s"[YRegionGateway] Lifting over $regionType from $fromBuild to $toBuild")

    for {
      chainFile <- liftoverGateway.resolve(fromBuild, toBuild)
      liftedPath <- if (YRegionCache.gff3Types.contains(regionType)) {
        liftoverGff3(sourcePath, chainFile, regionType, toBuild)
      } else {
        liftoverBed(sourcePath, chainFile, regionType, toBuild)
      }
    } yield {
      val cachedPath = cache.put(regionType, toBuild, liftedPath)
      onProgress(s"$regionType liftover complete.", 1.0)
      println(s"[YRegionGateway] Cached lifted $regionType at $cachedPath")
      cachedPath
    }
  }

  /**
   * Liftover a GFF3 file.
   */
  private def liftoverGff3(
    gff3Path: Path,
    chainPath: Path,
    regionType: String,
    targetBuild: String
  ): Either[String, Path] = {
    // Suppress verbose htsjdk logging
    Log.setGlobalLogLevel(Log.LogLevel.WARNING)

    val outputPath = Files.createTempFile(s"yregion-$regionType-$targetBuild-", ".gff3")

    try {
      val liftOver = new LiftOver(chainPath.toFile)
      var liftedCount = 0
      var failedCount = 0

      Using.resources(
        Source.fromFile(gff3Path.toFile),
        new PrintWriter(outputPath.toFile)
      ) { (source, writer) =>
        // Write header
        writer.println("##gff-version 3")
        writer.println(s"# Lifted over from GRCh38 to $targetBuild")

        for (line <- source.getLines() if !line.startsWith("#") && line.nonEmpty) {
          val fields = line.split("\t")
          if (fields.length >= 9) {
            val chrom = normalizeChromForChain(fields(0))
            val start = fields(3).toInt
            val end = fields(4).toInt

            // Create interval (htsjdk uses 1-based coordinates like GFF3)
            val interval = new Interval(chrom, start, end)
            val lifted = liftOver.liftOver(interval)

            if (lifted != null) {
              // Write lifted record
              val liftedFields = fields.clone()
              liftedFields(0) = lifted.getContig
              liftedFields(3) = lifted.getStart.toString
              liftedFields(4) = lifted.getEnd.toString
              writer.println(liftedFields.mkString("\t"))
              liftedCount += 1
            } else {
              failedCount += 1
            }
          }
        }
      }

      println(s"[YRegionGateway] GFF3 liftover complete: $liftedCount regions lifted, $failedCount failed")

      if (liftedCount > 0) Right(outputPath)
      else {
        Files.deleteIfExists(outputPath)
        Left("GFF3 liftover produced no valid output")
      }
    } catch {
      case e: Exception =>
        Files.deleteIfExists(outputPath)
        Left(s"GFF3 liftover failed: ${e.getMessage}")
    }
  }

  /**
   * Liftover a BED file.
   */
  private def liftoverBed(
    bedPath: Path,
    chainPath: Path,
    regionType: String,
    targetBuild: String
  ): Either[String, Path] = {
    // Suppress verbose htsjdk logging
    Log.setGlobalLogLevel(Log.LogLevel.WARNING)

    val outputPath = Files.createTempFile(s"yregion-$regionType-$targetBuild-", ".bed")

    try {
      val liftOver = new LiftOver(chainPath.toFile)
      var liftedCount = 0
      var failedCount = 0

      Using.resources(
        Source.fromFile(bedPath.toFile),
        new PrintWriter(outputPath.toFile)
      ) { (source, writer) =>
        for (line <- source.getLines() if !line.startsWith("#") && !line.startsWith("track") && !line.startsWith("browser") && line.nonEmpty) {
          val fields = line.split("\t")
          if (fields.length >= 3) {
            val chrom = normalizeChromForChain(fields(0))
            val start = fields(1).toInt  // BED is 0-based
            val end = fields(2).toInt

            // Create interval (htsjdk uses 1-based, so convert BED 0-based)
            val interval = new Interval(chrom, start + 1, end)
            val lifted = liftOver.liftOver(interval)

            if (lifted != null) {
              // Convert back to BED (0-based)
              val liftedStart = lifted.getStart - 1
              val liftedEnd = lifted.getEnd
              val liftedFields = fields.clone()
              liftedFields(0) = lifted.getContig
              liftedFields(1) = liftedStart.toString
              liftedFields(2) = liftedEnd.toString
              writer.println(liftedFields.mkString("\t"))
              liftedCount += 1
            } else {
              failedCount += 1
            }
          }
        }
      }

      println(s"[YRegionGateway] BED liftover complete: $liftedCount regions lifted, $failedCount failed")

      if (liftedCount > 0) Right(outputPath)
      else {
        Files.deleteIfExists(outputPath)
        Left("BED liftover produced no valid output")
      }
    } catch {
      case e: Exception =>
        Files.deleteIfExists(outputPath)
        Left(s"BED liftover failed: ${e.getMessage}")
    }
  }

  /**
   * Normalize chromosome name for chain file lookup.
   */
  private def normalizeChromForChain(chrom: String): String = {
    if (chrom.startsWith("chr")) chrom
    else if (chrom == "MT") "chrM"
    else s"chr$chrom"
  }

  /**
   * Load a YRegionAnnotator for a reference build.
   * Downloads all required region files and creates an annotator.
   *
   * @param referenceBuild Target reference build
   * @param callableLociPath Optional path to callable_loci.bed file
   * @return Either error message or configured annotator
   */
  def loadAnnotator(
    referenceBuild: String,
    callableLociPath: Option[Path] = None
  ): Either[String, YRegionAnnotator] = {
    for {
      paths <- resolveAll(referenceBuild)
      annotator <- createAnnotator(paths, referenceBuild, callableLociPath)
    } yield annotator
  }

  /**
   * Create an annotator from resolved paths.
   */
  private def createAnnotator(
    paths: YRegionPaths,
    referenceBuild: String,
    callableLociPath: Option[Path]
  ): Either[String, YRegionAnnotator] = {
    try {
      // Parse GFF3 files
      val cytobands = RegionFileParser.parseGff3(paths.cytobands)
        .map(records => YRegionAnnotator.gff3ToRegions(records, RegionType.Cytoband))
        .getOrElse(Nil)

      val palindromes = RegionFileParser.parseGff3(paths.palindromes)
        .map(records => YRegionAnnotator.gff3ToRegions(records, RegionType.Palindrome))
        .getOrElse(Nil)

      val strs = RegionFileParser.parseGff3(paths.strs)
        .map(records => YRegionAnnotator.gff3ToRegions(records, RegionType.STR))
        .getOrElse(Nil)

      // Parse BED files
      val pars = RegionFileParser.parseBed(paths.par)
        .map(records => YRegionAnnotator.bedToRegions(records, RegionType.PAR))
        .getOrElse(Nil)

      val xtrs = RegionFileParser.parseBed(paths.xtr)
        .map(records => YRegionAnnotator.bedToRegions(records, RegionType.XTR))
        .getOrElse(Nil)

      val ampliconic = RegionFileParser.parseBed(paths.ampliconic)
        .map(records => YRegionAnnotator.bedToRegions(records, RegionType.Ampliconic))
        .getOrElse(Nil)

      // Get hardcoded heterochromatin boundaries
      val heterochromatin = referenceBuild match {
        case "GRCh38" => YRegionAnnotator.grch38Heterochromatin
        case "GRCh37" => YRegionAnnotator.grch37Heterochromatin
        case _ => Nil
      }

      // Parse callable loci if provided
      val callablePositions: Option[Set[Long]] = callableLociPath.flatMap { path =>
        RegionFileParser.parseBed(path).toOption.map { records =>
          // Convert BED intervals to set of callable positions
          // For performance, we only track Y chromosome
          val yRecords = RegionFileParser.filterYChromosome(records, _.chrom)
          yRecords.flatMap { r =>
            val (start, end) = RegionFileParser.bedToOneBased(r.start, r.end)
            (start to end)
          }.toSet
        }
      }

      val annotator = YRegionAnnotator.fromRegions(
        cytobands = cytobands,
        palindromes = palindromes,
        strs = strs,
        pars = pars,
        xtrs = xtrs,
        ampliconic = ampliconic,
        heterochromatin = heterochromatin,
        callablePositions = callablePositions
      )

      println(s"[YRegionGateway] Created annotator with ${annotator.totalRegionCount} regions")
      Right(annotator)
    } catch {
      case e: Exception =>
        Left(s"Failed to create annotator: ${e.getMessage}")
    }
  }
}

object YRegionGateway {
  /**
   * Create a gateway with no progress callback.
   */
  def apply(): YRegionGateway = new YRegionGateway()

  /**
   * Create a gateway with a progress callback.
   */
  def apply(onProgress: (String, Double) => Unit): YRegionGateway = new YRegionGateway(onProgress)
}
