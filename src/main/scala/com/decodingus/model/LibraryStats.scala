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
   * Priority:
   * 1. PU (Platform Unit) - most unique, includes flowcell+lane+barcode
   * 2. LB + SM - library + sample combination
   * 3. Read stats - fallback when headers are incomplete
   */
  def computeRunFingerprint: String = {
    import java.security.MessageDigest

    val input = platformUnit match {
      case Some(pu) if pu.nonEmpty =>
        // Best case: PU uniquely identifies the run
        s"PU:$pu"
      case _ if libraryId != "Unknown" && sampleName != "Unknown" =>
        // Good case: LB + SM combination (GATK required fields)
        s"LB:$libraryId:SM:$sampleName:PL:$inferredPlatform"
      case _ =>
        // Fallback: use read statistics
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
    platformUnit match {
      case Some(pu) if pu.nonEmpty => "HIGH"
      case _ if libraryId != "Unknown" && sampleName != "Unknown" => "MEDIUM"
      case _ => "LOW"
    }
  }
}
