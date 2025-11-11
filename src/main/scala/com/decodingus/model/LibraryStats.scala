package com.decodingus.model

import scala.collection.mutable

case class LibraryStats(
                         readCount: Int = 0,
                         pairedReads: Int = 0,
                         lengthDistribution: mutable.Map[Int, Int] = mutable.Map(),
                         insertSizeDistribution: mutable.Map[Long, Int] = mutable.Map(),
                         aligner: String = "Unknown",
                         referenceBuild: String = "Unknown",
                         sampleName: String = "Unknown",
                         flowCells: mutable.Map[String, Int] = mutable.Map(),
                         instruments: mutable.Map[String, Int] = mutable.Map(),
                         mostFrequentInstrument: String = "Unknown",
                         inferredPlatform: String = "Unknown",
                         platformCounts: mutable.Map[String, Int] = mutable.Map()
                       )
