package com.decodingus.refgenome

import com.decodingus.util.Logger
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
  centromeres: Option[Path] = None,
  sequenceClass: Option[Path] = None  // T2T sequence class file (contains X-DEG, etc.)
)

/**
 * Gateway for downloading and caching Y chromosome region annotation files.
 *
 * Downloads from build-specific native sources when available:
 *
 * GRCh38:
 * - ybrowse.org: cytobands, palindromes, STRs (GFF3 format)
 * - GIAB genome-stratifications: PAR, XTR, ampliconic (BED format)
 *
 * CHM13v2.0 (hs1) - Native T2T annotations for best quality:
 * - T2T/CHM13 S3: amplicons, palindromes, cytobands, sequence class (BED format)
 * - GIAB genome-stratifications: PAR, XTR (BED format)
 *
 * GRCh37:
 * - Liftover from GRCh38 (no native files available)
 *
 * @param onProgress Callback for progress updates (message, 0.0-1.0)
 */
class YRegionGateway(onProgress: (String, Double) => Unit = (_, _) => ()) {
  private val log = Logger[YRegionGateway]
  private val cache = new YRegionCache
  private val liftoverGateway = new LiftoverGateway((_, _) => ())

  // ========== GRCh38 Sources ==========

  // ybrowse.org GFF3 URLs (GRCh38/hg38 coordinates)
  private val grch38YbrowseUrls: Map[String, String] = Map(
    "cytobands" -> "https://ybrowse.org/gbrowse2/gff/cytobands_hg38.gff3",
    "palindromes" -> "https://ybrowse.org/gbrowse2/gff/palindromes_hg38.gff3",
    "strs" -> "https://ybrowse.org/gbrowse2/gff/str_hg38.gff3"
  )

  // GIAB genome-stratifications BED URLs (GRCh38)
  // From: https://github.com/genome-in-a-bottle/genome-stratifications
  private val grch38GiabUrls: Map[String, String] = Map(
    "par" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/XY/GRCh38_chrY_PAR.bed",
    "xtr" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/XY/GRCh38_chrY_XTR.bed",
    "ampliconic" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/XY/GRCh38_chrY_ampliconic.bed"
  )

  // ========== GRCh37 Sources ==========

  // ybrowse.org GFF3 URLs (GRCh37/hg19 coordinates) - native files, no liftover needed!
  private val grch37YbrowseUrls: Map[String, String] = Map(
    "cytobands" -> "https://ybrowse.org/gbrowse2/gff/cytobands_hg19.gff3",
    "palindromes" -> "https://ybrowse.org/gbrowse2/gff/palindromes_hg19.gff3",
    "strs" -> "https://ybrowse.org/gbrowse2/gff/str_hg19.gff3"
  )

  // GIAB genome-stratifications BED URLs (GRCh37)
  // From: https://github.com/genome-in-a-bottle/genome-stratifications
  private val grch37GiabUrls: Map[String, String] = Map(
    "par" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh37/XY/GRCh37_chrY_PAR.bed",
    "xtr" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh37/XY/GRCh37_chrY_XTR.bed",
    "ampliconic" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh37/XY/GRCh37_chrY_ampliconic.bed"
  )

  // ========== CHM13v2.0 (hs1) Native Sources ==========

  // T2T CHM13 S3 URLs - Native annotations from the T2T consortium
  // From: https://github.com/marbl/CHM13
  // These are the gold standard for CHM13v2.0 with 30+ Mbp more Y sequence than GRCh38
  private val chm13T2tUrls: Map[String, String] = Map(
    "cytobands" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed",
    "palindromes" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_inverted_repeats_v1.bed",
    "ampliconic" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_amplicons_v1.bed",
    "sequence_class" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_chrXY_sequence_class_v1.bed",
    "azf_dyz" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_AZF_DYZ_v1.bed",
    "censat" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed"
  )

  // GIAB genome-stratifications BED URLs (CHM13v2.0)
  // From: https://github.com/genome-in-a-bottle/genome-stratifications/tree/master/CHM13v2.0/XY
  private val chm13GiabUrls: Map[String, String] = Map(
    "par" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/CHM13v2.0/XY/CHM13v2.0_chrY_PAR.bed",
    "xtr" -> "https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/CHM13v2.0/XY/CHM13v2.0_chrY_XTR.bed"
  )

  // Note: All builds now have native sources - liftover is only used as fallback for missing region types

  /**
   * Resolve all region files for a reference build.
   *
   * @param referenceBuild Target reference build (GRCh38, GRCh37, CHM13v2)
   * @return Either error message or paths to all region files
   */
  def resolveAll(referenceBuild: String): Either[String, YRegionPaths] = {
    val normalizedBuild = normalizeBuildName(referenceBuild)

    for {
      cytobands <- resolve("cytobands", referenceBuild)
      palindromes <- resolve("palindromes", referenceBuild)
      strs <- resolve("strs", referenceBuild)
      par <- resolve("par", referenceBuild)
      xtr <- resolve("xtr", referenceBuild)
      ampliconic <- resolve("ampliconic", referenceBuild)
    } yield {
      // sequence_class is only available natively for CHM13v2 (contains X-DEG regions)
      val sequenceClass = if (normalizedBuild == "CHM13v2") {
        resolve("sequence_class", referenceBuild).toOption
      } else {
        None
      }
      YRegionPaths(cytobands, palindromes, strs, par, xtr, ampliconic, sequenceClass = sequenceClass)
    }
  }

  /**
   * Resolve a single region file for a reference build.
   *
   * @param regionType Type of region (cytobands, palindromes, strs, par, xtr, ampliconic)
   * @param referenceBuild Target reference build
   * @return Either error message or path to the region file
   */
  def resolve(regionType: String, referenceBuild: String): Either[String, Path] = {
    // Normalize build name
    val normalizedBuild = normalizeBuildName(referenceBuild)

    // Check cache first
    cache.getPath(regionType, normalizedBuild) match {
      case Some(path) =>
        log.debug(s"Found $regionType for $normalizedBuild in cache: $path")
        Right(path)
      case None =>
        normalizedBuild match {
          case "GRCh38" =>
            // Download directly for GRCh38
            downloadGrch38RegionFile(regionType)

          case "GRCh37" =>
            // Download native GRCh37 files from ybrowse + GIAB
            downloadGrch37RegionFile(regionType)

          case "CHM13v2" =>
            // Download native CHM13v2.0 files from T2T + GIAB - best quality!
            downloadChm13RegionFile(regionType)

          case _ =>
            Left(s"Unsupported reference build: $referenceBuild")
        }
    }
  }

  /**
   * Normalize build name to canonical form.
   */
  private def normalizeBuildName(build: String): String = {
    build.toLowerCase match {
      case b if b.contains("chm13") || b == "hs1" || b.contains("t2t") => "CHM13v2"
      case b if b.contains("grch38") || b == "hg38" => "GRCh38"
      case b if b.contains("grch37") || b == "hg19" => "GRCh37"
      case _ => build
    }
  }

  /**
   * Download a GRCh38 region file from ybrowse or GIAB.
   */
  private def downloadGrch38RegionFile(regionType: String): Either[String, Path] = {
    val url = grch38YbrowseUrls.get(regionType).orElse(grch38GiabUrls.get(regionType))
    url match {
      case Some(u) =>
        val isGff3 = grch38YbrowseUrls.contains(regionType)
        downloadFile(regionType, "GRCh38", u, isGff3)
      case None => Left(s"Unknown region type for GRCh38: $regionType")
    }
  }

  /**
   * Download a GRCh37 region file from ybrowse or GIAB.
   * Native files available - no liftover needed!
   */
  private def downloadGrch37RegionFile(regionType: String): Either[String, Path] = {
    val url = grch37YbrowseUrls.get(regionType).orElse(grch37GiabUrls.get(regionType))
    url match {
      case Some(u) =>
        val isGff3 = grch37YbrowseUrls.contains(regionType)
        downloadFile(regionType, "GRCh37", u, isGff3)
      case None => Left(s"Unknown region type for GRCh37: $regionType")
    }
  }

  /**
   * Download a CHM13v2.0 region file from native T2T or GIAB sources.
   * CHM13v2.0 uses native annotations - no liftover needed!
   */
  private def downloadChm13RegionFile(regionType: String): Either[String, Path] = {
    // CHM13v2.0 doesn't have STRs in T2T annotations - liftover from GRCh38 as fallback
    if (regionType == "strs") {
      log.info("CHM13v2.0 STRs not available natively, lifting over from GRCh38")
      return resolve("strs", "GRCh38").flatMap { grch38Path =>
        liftoverRegionFile(grch38Path, "strs", "GRCh38", "CHM13v2")
      }
    }

    val url = chm13T2tUrls.get(regionType).orElse(chm13GiabUrls.get(regionType))
    url match {
      case Some(u) =>
        // All CHM13v2.0 files are BED format
        downloadFile(regionType, "CHM13v2", u, isGff3 = false)
      case None => Left(s"Unknown region type for CHM13v2.0: $regionType")
    }
  }

  /**
   * Download a file from URL and cache it.
   *
   * @param regionType Type of region
   * @param referenceBuild Target build
   * @param url Download URL
   * @param isGff3 True if file is GFF3 format, false for BED format
   */
  private def downloadFile(regionType: String, referenceBuild: String, url: String, isGff3: Boolean): Either[String, Path] = {
    onProgress(s"Downloading $regionType for $referenceBuild...", 0.0)
    log.info(s"Downloading $regionType from $url")

    val ext = if (isGff3) ".gff3" else ".bed"
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
          log.debug(s"Cached $regionType at $cachedPath")
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
    log.info(s"Lifting over $regionType from $fromBuild to $toBuild")

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
      log.debug(s"Cached lifted $regionType at $cachedPath")
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

      log.info(s"GFF3 liftover complete: $liftedCount regions lifted, $failedCount failed")

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

      log.info(s"BED liftover complete: $liftedCount regions lifted, $failedCount failed")

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
    val normalizedBuild = normalizeBuildName(referenceBuild)

    try {
      // CHM13v2.0 uses BED format for all files (from T2T consortium)
      // GRCh38/GRCh37 use GFF3 from ybrowse for cytobands/palindromes/strs
      val useBedFormat = normalizedBuild == "CHM13v2"

      val cytobands = if (useBedFormat) {
        RegionFileParser.parseBed(paths.cytobands)
          .map(records => YRegionAnnotator.bedToRegions(
            RegionFileParser.filterYChromosome(records, _.chrom),
            RegionType.Cytoband
          ))
          .getOrElse(Nil)
      } else {
        RegionFileParser.parseGff3(paths.cytobands)
          .map(records => YRegionAnnotator.gff3ToRegions(records, RegionType.Cytoband))
          .getOrElse(Nil)
      }

      val palindromes = if (useBedFormat) {
        RegionFileParser.parseBed(paths.palindromes)
          .map(records => YRegionAnnotator.bedToRegions(records, RegionType.Palindrome))
          .getOrElse(Nil)
      } else {
        RegionFileParser.parseGff3(paths.palindromes)
          .map(records => YRegionAnnotator.gff3ToRegions(records, RegionType.Palindrome))
          .getOrElse(Nil)
      }

      // STRs are always GFF3 (from ybrowse) - CHM13v2 falls back to liftover
      val strs = RegionFileParser.parseGff3(paths.strs)
        .map(records => YRegionAnnotator.gff3ToRegions(records, RegionType.STR))
        .getOrElse(Nil)

      // BED files (same format for all builds)
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
      // CHM13v2.0 has better heterochromatin definition from censat file (future enhancement)
      val heterochromatin = normalizedBuild match {
        case "GRCh38" => YRegionAnnotator.grch38Heterochromatin
        case "GRCh37" => YRegionAnnotator.grch37Heterochromatin
        case "CHM13v2" => YRegionAnnotator.chm13v2Heterochromatin
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

      // Parse X-degenerate regions from sequence_class file (CHM13v2 only)
      // The sequence_class file has entries like: chrY  start  end  X-DEG
      val xdegenerate = paths.sequenceClass.flatMap { seqClassPath =>
        RegionFileParser.parseBed(seqClassPath).toOption.map { records =>
          // Filter for X-DEG (X-degenerate) entries on Y chromosome
          val xdegRecords = records.filter { r =>
            r.name.exists(_.toUpperCase.contains("X-DEG")) &&
              RegionFileParser.filterYChromosome(List(r), _.chrom).nonEmpty
          }
          YRegionAnnotator.bedToRegions(xdegRecords, RegionType.XDegenerate)
        }
      }.getOrElse(Nil)

      val annotator = YRegionAnnotator.fromRegions(
        cytobands = cytobands,
        palindromes = palindromes,
        strs = strs,
        pars = pars,
        xtrs = xtrs,
        ampliconic = ampliconic,
        heterochromatin = heterochromatin,
        xdegenerate = xdegenerate,
        callablePositions = callablePositions
      )

      log.info(s"Created annotator with ${annotator.totalRegionCount} regions (${xdegenerate.size} X-degenerate)")
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
