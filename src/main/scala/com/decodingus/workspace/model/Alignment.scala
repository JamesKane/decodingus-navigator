package com.decodingus.workspace.model

/**
 * An alignment of sequence data to a reference genome.
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.alignment).
 *
 * Independently managed for granular updates - metrics can be updated without touching parent records.
 *
 * @param atUri           The AT URI of this alignment record
 * @param meta            Record metadata for tracking changes and sync
 * @param sequenceRunRef  AT URI of the parent sequence run record
 * @param biosampleRef    AT URI of the grandparent biosample (denormalized for query efficiency)
 * @param referenceBuild  Reference genome build (GRCh38, GRCh37, T2T-CHM13, hg19, hg38)
 * @param aligner         Tool and version used for alignment (e.g., BWA-MEM 0.7.17)
 * @param variantCaller   Tool used for variant calling (e.g., GATK HaplotypeCaller 4.2)
 * @param files           Metadata about aligned data files (e.g., BAM, CRAM, VCF). Files remain local.
 * @param metrics         Quality control metrics for the alignment
 */
case class Alignment(
  atUri: Option[String],
  meta: RecordMeta,
  sequenceRunRef: String,
  biosampleRef: Option[String] = None,
  referenceBuild: String,
  aligner: String,
  variantCaller: Option[String] = None,
  files: List[FileInfo] = List.empty,
  metrics: Option[AlignmentMetrics] = None
)

object Alignment {
  /** Known reference build values */
  val KnownReferenceBuilds: Set[String] = Set(
    "GRCh38", "GRCh37", "T2T-CHM13", "hg19", "hg38"
  )
}
