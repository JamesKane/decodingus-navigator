package com.decodingus.model

case class LibraryStats(
                         readCount: Int = 0,
                         pairedReads: Int = 0,
                         lengthDistribution: Map[Int, Int] = Map(),
                         insertSizeDistribution: Map[Long, Int] = Map(),
                         aligner: String = "Unknown",
                         referenceBuild: String = "Unknown",
                         sampleName: String = "Unknown",
                         libraryId: String = "Unknown", // @RG LB - GATK required, stable across re-alignments
                         platformUnit: Option[String] = None, // @RG PU - optional but best for fingerprinting
                         flowCells: Map[String, Int] = Map(),
                         instruments: Map[String, Int] = Map(),
                         mostFrequentInstrumentId: String = "Unknown",
                         mostFrequentInstrument: String = "Unknown",
                         inferredPlatform: String = "Unknown",
                         platformCounts: Map[String, Int] = Map()
                       ) {
  /**
   * Compute a fingerprint hash for identifying the same sequencing run
   * across different reference alignments.
   *
   * Uses SM (Sample Name) + Platform as the primary identifier since PU (Platform Unit)
   * is unreliable - most BAMs don't tag it correctly.
   *
   * Priority:
   * 1. SM + Platform - sample name on a specific platform
   * 2. Read stats - fallback when headers are incomplete
   */
  def computeRunFingerprint: String = {
    import java.security.MessageDigest

    val input = if (sampleName != "Unknown" && inferredPlatform != "Unknown") {
      // Primary: SM + Platform uniquely identifies a sequencing run
      s"SM:$sampleName:PL:$inferredPlatform"
    } else if (sampleName != "Unknown") {
      // Fallback: just sample name
      s"SM:$sampleName"
    } else {
      // Last resort: use read statistics
      val lengthHash = lengthDistribution.toSeq.sorted.hashCode()
      s"STATS:$readCount:$lengthHash:$pairedReads"
    }

    val md = MessageDigest.getInstance("SHA-256")
    md.update(input.getBytes("UTF-8"))
    md.digest().take(16).map("%02x".format(_)).mkString
  }

  /**
   * Confidence level of the fingerprint.
   */
  def fingerprintConfidence: String = {
    if (sampleName != "Unknown" && inferredPlatform != "Unknown") "HIGH"
    else if (sampleName != "Unknown") "MEDIUM"
    else "LOW"
  }
}
