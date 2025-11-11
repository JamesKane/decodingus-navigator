package com.decodingus.refgenome

import sttp.client3.*

import java.nio.file.{Files, Path}

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
        Right(path)
      case None =>
        referenceUrls.get(referenceBuild) match {
          case Some(url) => downloadReference(referenceBuild, url)
          case None => Left(s"Unknown reference build: $referenceBuild")
        }
    }
  }

  private def downloadReference(referenceBuild: String, url: String): Either[String, Path] = {
    println(s"Downloading reference $referenceBuild from $url")
    val tempFile = Files.createTempFile(s"ref-$referenceBuild", ".fa.gz")

    val request = basicRequest.get(uri"$url").response(asFile(tempFile.toFile))

    val backend = HttpURLConnectionBackend()
    val response = request.send(backend)

    response.body match {
      case Right(file) =>
        println("Download complete. Caching reference.")
        val finalPath = cache.put(referenceBuild, file.toPath)
        Right(finalPath)
      case Left(error) =>
        Left(s"Failed to download reference: $error")
    }
  }
}
