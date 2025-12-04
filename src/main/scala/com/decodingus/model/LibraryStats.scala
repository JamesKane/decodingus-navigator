package com.decodingus.model

case class LibraryStats(
                         readCount: Int = 0,
                         pairedReads: Int = 0,
                         lengthDistribution: Map[Int, Int] = Map(),
                         insertSizeDistribution: Map[Long, Int] = Map(),
                         aligner: String = "Unknown",
                         referenceBuild: String = "Unknown",
                         sampleName: String = "Unknown",
                         flowCells: Map[String, Int] = Map(),
                         instruments: Map[String, Int] = Map(),
                         mostFrequentInstrument: String = "Unknown",
                         inferredPlatform: String = "Unknown",
                         platformCounts: Map[String, Int] = Map()
                       )
