package com.decodingus.refgenome

import java.nio.file.Path

class LiftoverCache extends FileCache {
  protected def cacheSubdir = "liftover"

  def getPath(from: String, to: String): Option[Path] =
    cachedPath(s"${from}To${to.capitalize}.over.chain.gz")

  def put(from: String, to: String, file: Path): Path =
    moveToCache(file, s"${from}To${to.capitalize}.over.chain.gz")
}
