package com.decodingus.refgenome

import sttp.client3.*

import java.io.IOException
import java.nio.file.{Files, Path, Paths}
import sys.process._
import org.broadinstitute.hellbender.Main

class ReferenceGateway(onProgress: (Long, Long) => Unit) {
  private val cache = new ReferenceCache

  private val referenceUrls: Map[String, String] = Map(
    "GRCh38" -> "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
    "GRCh37" -> "https://storage.googleapis.com/genomics-public-data/references/hg19/v0/Homo_sapiens_assembly19.fasta.gz",
    "CHM13v2" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    // Add other references here
  )

  def resolve(referenceBuild: String): Either[String, Path] = {
    cache.getPath(referenceBuild) match {
      case Some(path) =>
        println(s"Found reference $referenceBuild in cache: $path")
        validateAndCreateReferenceFiles(path)
      case None =>
        referenceUrls.get(referenceBuild) match {
          case Some(url) =>
            downloadReference(referenceBuild, url).flatMap(validateAndCreateReferenceFiles)
          case None => Left(s"Unknown reference build: $referenceBuild")
        }
    }
  }

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
      println(s"Creating sequence dictionary for $referencePath using GATK library...")
      val args = Array(
        "CreateSequenceDictionary",
        "-R", referencePath.toAbsolutePath.toString,
        "-O", dictPath.toAbsolutePath.toString
      )
      try {
        Main.main(args)
      } catch {
        case e: Exception =>
          return Left(s"Failed to create sequence dictionary for $referencePath using GATK library: ${e.getMessage}")
      }
    }
    Right(referencePath)
  }
}
