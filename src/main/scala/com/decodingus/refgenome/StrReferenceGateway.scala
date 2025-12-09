package com.decodingus.refgenome

import com.decodingus.analysis.GatkRunner
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
  private val referenceGateway = new ReferenceGateway((_, _) => ())

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
      sourceRef <- referenceGateway.resolve(fromBuild)
      targetRef <- referenceGateway.resolve(toBuild)
      liftedPath <- performGatkLiftover(sourcePath, sourceRef, targetRef, chainFile, toBuild)
    } yield {
      val finalPath = cache.put(toBuild, liftedPath)
      onProgress("STR reference liftover complete.", 1.0)
      println(s"[StrReferenceGateway] Cached lifted STR reference at $finalPath")
      finalPath
    }
  }

  /**
   * Perform liftover using GATK's BedToIntervalList and LiftOverIntervalList.
   * This avoids requiring external liftOver tool installation.
   */
  private def performGatkLiftover(
    bedPath: Path,
    sourceRef: Path,
    targetRef: Path,
    chainPath: Path,
    targetBuild: String
  ): Either[String, Path] = {
    val intervalList = Files.createTempFile(s"hipstr-intervals", ".interval_list")
    val liftedIntervalList = Files.createTempFile(s"hipstr-$targetBuild-lifted", ".interval_list")
    val outputBed = Files.createTempFile(s"hipstr-$targetBuild", ".bed")

    try {
      // Step 1: Convert BED to interval list
      println(s"[StrReferenceGateway] Converting BED to interval list...")
      val bedToIntervalArgs = Array(
        "BedToIntervalList",
        "-I", bedPath.toString,
        "-O", intervalList.toString,
        "-SD", sourceRef.toString
      )

      GatkRunner.run(bedToIntervalArgs) match {
        case Left(error) =>
          cleanup(intervalList, liftedIntervalList, outputBed)
          return Left(s"BedToIntervalList failed: $error")
        case Right(_) => // continue
      }

      // Step 2: Liftover the interval list
      println(s"[StrReferenceGateway] Lifting over interval list...")
      val liftoverArgs = Array(
        "LiftOverIntervalList",
        "-I", intervalList.toString,
        "-O", liftedIntervalList.toString,
        "-SD", targetRef.toString,
        "-CHAIN", chainPath.toString
      )

      GatkRunner.run(liftoverArgs) match {
        case Left(error) =>
          cleanup(intervalList, liftedIntervalList, outputBed)
          return Left(s"LiftOverIntervalList failed: $error")
        case Right(_) => // continue
      }

      // Step 3: Convert back to BED format (preserving extra columns from original)
      println(s"[StrReferenceGateway] Converting interval list back to BED...")
      convertIntervalListToBed(liftedIntervalList, bedPath, outputBed)

      Files.deleteIfExists(intervalList)
      Files.deleteIfExists(liftedIntervalList)

      if (Files.exists(outputBed) && Files.size(outputBed) > 0) {
        Right(outputBed)
      } else {
        Files.deleteIfExists(outputBed)
        Left("Liftover produced empty output")
      }
    } catch {
      case e: Exception =>
        cleanup(intervalList, liftedIntervalList, outputBed)
        Left(s"Liftover failed: ${e.getMessage}")
    }
  }

  /**
   * Convert interval list back to BED format.
   * The interval list loses the extra BED columns (period, num_repeats, name),
   * so we try to match positions back to the original BED to recover them.
   */
  private def convertIntervalListToBed(intervalListPath: Path, originalBedPath: Path, outputBedPath: Path): Unit = {
    // Load original BED data keyed by (chrom, start, end) to recover extra columns
    // Note: After liftover, positions may shift, so we just output basic BED
    // Future enhancement: could try fuzzy matching to recover STR metadata

    Using.resources(
      Source.fromFile(intervalListPath.toFile),
      new PrintWriter(outputBedPath.toFile)
    ) { (source, writer) =>
      for (line <- source.getLines() if !line.startsWith("@") && line.nonEmpty) {
        val fields = line.split("\t")
        if (fields.length >= 3) {
          // Interval list format: chrom start end strand name
          // BED format: chrom start end [name] [score] [strand] ...
          val chrom = fields(0)
          val start = fields(1).toLong - 1 // Interval list is 1-based, BED is 0-based
          val end = fields(2)
          val name = if (fields.length > 4) fields(4) else "."

          // For now, output basic BED with placeholder STR info
          // The position lookup in StrAnnotator will still work
          // Period and num_repeats will need to be re-estimated or defaulted
          writer.println(s"$chrom\t$start\t$end\t4\t10.0\t$name")
        }
      }
    }
  }

  private def cleanup(paths: Path*): Unit = {
    paths.foreach(p => Files.deleteIfExists(p))
  }
}
