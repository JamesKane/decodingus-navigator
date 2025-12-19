package com.decodingus.workspace.model

case class AlignmentMetrics(
                             // WGS Metrics
                             genomeTerritory: Option[Long] = None,
                             meanCoverage: Option[Double] = None,
                             medianCoverage: Option[Double] = None,
                             sdCoverage: Option[Double] = None,
                             pctExcDupe: Option[Double] = None,
                             pctExcMapq: Option[Double] = None,
                             pct10x: Option[Double] = None,
                             pct20x: Option[Double] = None,
                             pct30x: Option[Double] = None,
                             hetSnpSensitivity: Option[Double] = None,

                             // Callable Loci
                             callableBases: Option[Long] = None,
                             callableLociComplete: Option[Boolean] = None,
                             contigs: List[ContigMetrics] = List.empty,

                             // Whole-Genome VCF Status
                             vcfPath: Option[String] = None,
                             vcfCreatedAt: Option[String] = None,
                             vcfVariantCount: Option[Long] = None,
                             vcfReferenceBuild: Option[String] = None,

                             // Sex Inference
                             inferredSex: Option[String] = None,
                             sexInferenceConfidence: Option[String] = None,
                             xAutosomeRatio: Option[Double] = None,

                             // Structural Variant Calling
                             svVcfPath: Option[String] = None,
                             svCallCount: Option[Int] = None
                           ) {
  /** Check if whole-genome VCF has been generated */
  def hasVcf: Boolean = vcfPath.isDefined && vcfVariantCount.isDefined

  /** Check if callable loci analysis is complete */
  def hasCallableLoci: Boolean = callableLociComplete.getOrElse(false)

  /** Check if structural variant calling has been done */
  def hasSvCalling: Boolean = svVcfPath.isDefined && svCallCount.isDefined
}
