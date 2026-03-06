package com.decodingus.refgenome

import sttp.client3.*

import java.nio.file.{Files, Path}

class LiftoverGateway(onProgress: (Long, Long) => Unit) {
  private val cache = new LiftoverCache

  private val chainFileUrls: Map[(String, String), String] = Map(
    ("GRCh38", "GRCh37") -> "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    ("GRCh37", "GRCh38") -> "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    ("GRCh38", "CHM13v2") -> "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz",
    ("CHM13v2", "GRCh38") -> "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz"
    // Add other chain files here
  )

  def resolve(from: String, to: String): Either[String, Path] = {
    cache.getPath(from, to) match {
      case Some(path) =>
        println(s"Found chain file in cache: $path")
        Right(path)
      case None =>
        chainFileUrls.get((from, to)) match {
          case Some(url) => downloadChainFile(from, to, url)
          case None => Left(s"No chain file found for liftover from $from to $to")
        }
    }
  }

  private def downloadChainFile(from: String, to: String, url: String): Either[String, Path] = {
    println(s"Downloading chain file from $url")
    val tempFile = Files.createTempFile(s"chain-${from}To${to}", ".chain.gz")

    val request = basicRequest.get(uri"$url").response(asFile(tempFile.toFile))

    val backend = HttpURLConnectionBackend()
    val response = request.send(backend)

    response.body match {
      case Right(file) =>
        println("Download complete. Caching chain file.")
        val finalPath = cache.put(from, to, file.toPath)
        Right(finalPath)
      case Left(error) =>
        Left(s"Failed to download chain file: $error")
    }
  }
}
