package com.decodingus.refgenome

import java.nio.file.Path

/**
 * Cache for STR (Short Tandem Repeat) reference BED files.
 * Uses HipSTR reference files which contain known STR regions.
 */
class StrReferenceCache extends FileCache {
  protected def cacheSubdir = "str"

  def getPath(referenceBuild: String): Option[Path] =
    cachedPath(s"$referenceBuild.hipstr_reference.bed")

  def put(referenceBuild: String, file: Path): Path =
    moveToCache(file, s"$referenceBuild.hipstr_reference.bed")
}
