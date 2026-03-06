package com.decodingus.refgenome

import sttp.client3.*

import java.io.{BufferedInputStream, FileInputStream, FileOutputStream, PrintWriter}
import java.nio.file.{Files, Path}
import java.util.zip.GZIPInputStream
import scala.io.Source
import scala.util.Using

/**
 * Gateway for downloading and caching HipSTR STR reference BED files.
 *
 * HipSTR reference files contain known Short Tandem Repeat (STR) regions
 * identified using Tandem Repeats Finder. These are used to annotate
 * indels found during haplogroup analysis.
 *
 * Source: https://github.com/HipSTR-Tool/HipSTR-references/
 */
class StrReferenceGateway(onProgress: (String, Double) => Unit = (_, _) => ()) {
  private val cache = new StrReferenceCache
  private val liftoverGateway = new LiftoverGateway((_, _) => ())

  // HipSTR provides GRCh38 directly; we'll liftover for other builds
  private val hipstrUrls: Map[String, String] = Map(
    "GRCh38" -> "https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/GRCh38.hipstr_reference.bed.gz"
  )

  // Builds that need liftover from GRCh38
  private val liftoverBuilds: Set[String] = Set("GRCh37", "CHM13v2")

  /**
   * Resolve STR reference BED file for the given reference build.
   * Downloads from HipSTR if GRCh38, or lifts over from GRCh38 for other builds.
   */
  def resolve(referenceBuild: String): Either[String, Path] = {
    cache.getPath(referenceBuild) match {
      case Some(path) =>
        println(s"[StrReferenceGateway] Found STR reference in cache: $path")
        Right(path)
      case None =>
        if (hipstrUrls.contains(referenceBuild)) {
          downloadAndDecompress(referenceBuild, hipstrUrls(referenceBuild))
        } else if (liftoverBuilds.contains(referenceBuild)) {
          // First ensure we have GRCh38, then liftover
          resolve("GRCh38").flatMap { grch38Path =>
            liftoverStrReference(grch38Path, "GRCh38", referenceBuild)
          }
        } else {
          Left(s"No STR reference available for build: $referenceBuild")
        }
    }
  }

  private def downloadAndDecompress(referenceBuild: String, url: String): Either[String, Path] = {
    onProgress(s"Downloading HipSTR reference for $referenceBuild...", 0.0)
    println(s"[StrReferenceGateway] Downloading STR reference from $url")

    val tempGz = Files.createTempFile(s"hipstr-$referenceBuild", ".bed.gz")
    val tempBed = Files.createTempFile(s"hipstr-$referenceBuild", ".bed")

    try {
      val request = basicRequest.get(uri"$url").response(asFile(tempGz.toFile))
      val backend = HttpURLConnectionBackend()
      val response = request.send(backend)

      response.body match {
        case Right(_) =>
          onProgress("Decompressing...", 0.5)
          // Decompress the gzip file
          Using.resources(
            new GZIPInputStream(new BufferedInputStream(new FileInputStream(tempGz.toFile))),
            new FileOutputStream(tempBed.toFile)
          ) { (gzIn, out) =>
            val buffer = new Array[Byte](8192)
            var len = gzIn.read(buffer)
            while (len > 0) {
              out.write(buffer, 0, len)
              len = gzIn.read(buffer)
            }
          }

          Files.deleteIfExists(tempGz)
          val finalPath = cache.put(referenceBuild, tempBed)
          onProgress("STR reference ready.", 1.0)
          println(s"[StrReferenceGateway] Cached STR reference at $finalPath")
          Right(finalPath)

        case Left(error) =>
          Files.deleteIfExists(tempGz)
          Files.deleteIfExists(tempBed)
          Left(s"Failed to download STR reference: $error")
      }
    } catch {
      case e: Exception =>
        Files.deleteIfExists(tempGz)
        Files.deleteIfExists(tempBed)
        Left(s"Error downloading STR reference: ${e.getMessage}")
    }
  }

  private def liftoverStrReference(sourcePath: Path, fromBuild: String, toBuild: String): Either[String, Path] = {
    onProgress(s"Lifting over STR reference from $fromBuild to $toBuild...", 0.0)
    println(s"[StrReferenceGateway] Lifting over STR reference from $fromBuild to $toBuild")

    for {
      chainFile <- liftoverGateway.resolve(fromBuild, toBuild)
      liftedPath <- performLiftover(sourcePath, chainFile, toBuild)
    } yield {
      val finalPath = cache.put(toBuild, liftedPath)
      onProgress("STR reference liftover complete.", 1.0)
      println(s"[StrReferenceGateway] Cached lifted STR reference at $finalPath")
      finalPath
    }
  }

  /**
   * Normalize chromosome name to match chain file conventions.
   * HipSTR BED uses "1", "2", etc. but UCSC chain files use "chr1", "chr2", etc.
   */
  private def normalizeChromForChain(chrom: String): String = {
    if (chrom.startsWith("chr")) chrom
    else if (chrom == "MT") "chrM"
    else s"chr$chrom"
  }

  /**
   * Perform liftover of BED file using htsjdk's LiftOver directly.
   * This is much faster than GATK's interval list approach for large BED files.
   */
  private def performLiftover(
                               bedPath: Path,
                               chainPath: Path,
                               targetBuild: String
                             ): Either[String, Path] = {
    import htsjdk.samtools.liftover.LiftOver
    import htsjdk.samtools.util.{Interval, Log}

    // Suppress verbose htsjdk logging for failed liftover regions
    Log.setGlobalLogLevel(Log.LogLevel.WARNING)

    val outputBed = Files.createTempFile(s"hipstr-$targetBuild", ".bed")

    try {
      println(s"[StrReferenceGateway] Loading chain file for liftover...")
      val liftOver = new LiftOver(chainPath.toFile)

      var liftedCount = 0
      var failedCount = 0

      println(s"[StrReferenceGateway] Lifting over STR regions (this may take a moment)...")
      Using.resources(
        Source.fromFile(bedPath.toFile),
        new PrintWriter(outputBed.toFile)
      ) { (source, writer) =>
        for (line <- source.getLines() if !line.startsWith("#") && line.nonEmpty) {
          val fields = line.split("\t")
          if (fields.length >= 5) {
            val chrom = fields(0)
            val start = fields(1).toInt // BED is 0-based
            val end = fields(2).toInt
            val period = fields(3)
            val numRepeats = fields(4)
            val name = if (fields.length > 5) fields(5) else ""

            // Normalize chromosome name for chain file lookup (add "chr" prefix if needed)
            val chainChrom = normalizeChromForChain(chrom)

            // Create interval (htsjdk uses 1-based coordinates)
            val interval = new Interval(chainChrom, start + 1, end)
            val lifted = liftOver.liftOver(interval)

            if (lifted != null) {
              // Convert back to BED (0-based)
              // Keep the output contig name as-is from the chain (will have chr prefix)
              val liftedChrom = lifted.getContig
              val liftedStart = lifted.getStart - 1
              val liftedEnd = lifted.getEnd
              writer.println(s"$liftedChrom\t$liftedStart\t$liftedEnd\t$period\t$numRepeats\t$name")
              liftedCount += 1
            } else {
              failedCount += 1
            }
          }
        }
      }

      println(s"[StrReferenceGateway] Liftover complete: $liftedCount regions lifted, $failedCount failed")

      if (liftedCount > 0) {
        Right(outputBed)
      } else {
        Files.deleteIfExists(outputBed)
        Left("Liftover produced no valid output")
      }
    } catch {
      case e: Exception =>
        Files.deleteIfExists(outputBed)
        Left(s"Liftover failed: ${e.getMessage}")
    }
  }

}
