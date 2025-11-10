package com.decodingus.model

import scala.collection.mutable

case class LibraryStats(
  readCount: Int = 0,
  totalReadLength: Long = 0,
  pairedReads: Int = 0,
  totalInsertSize: Long = 0,
  pairedCount: Int = 0,
  lengthDistribution: mutable.Map[Int, Int] = mutable.Map(),
  insertSizeDistribution: mutable.Map[Long, Int] = mutable.Map(),
  aligner: String = "Unknown",
  referenceBuild: String = "Unknown",
  flowCells: mutable.Map[String, Int] = mutable.Map(),
  instruments: mutable.Map[String, Int] = mutable.Map(),
  platformCounts: mutable.Map[String, Int] = mutable.Map(),
  genomeSize: Long = 0,
  averageDepth: Double = 0.0
)