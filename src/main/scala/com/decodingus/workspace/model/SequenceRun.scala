package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * A sequencing run representing one library preparation and sequencing session.
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.sequencerun).
 *
 * Can be independently managed - created, updated, or deleted without cascading changes.
 *
 * @param atUri              The AT URI of this sequence run record
 * @param meta               Record metadata for tracking changes and sync
 * @param biosampleRef       AT URI of the parent biosample record
 * @param platformName       Sequencing platform (ILLUMINA, PACBIO, NANOPORE, ION_TORRENT, BGI, ELEMENT, ULTIMA)
 * @param instrumentModel    Specific instrument model (e.g., NovaSeq 6000, Sequel II)
 * @param instrumentId       Unique instrument identifier extracted from @RG headers (for lab inference)
 * @param testType           Type of test (WGS, EXOME, TARGETED, RNA_SEQ, AMPLICON)
 * @param libraryLayout      Paired-end or Single-end (PAIRED, SINGLE)
 * @param totalReads         Total number of reads (from CollectAlignmentSummaryMetrics TOTAL_READS)
 * @param pfReads            Pass-filter reads (PF_READS)
 * @param pfReadsAligned     PF reads that aligned (PF_READS_ALIGNED)
 * @param pctPfReadsAligned  Percentage of PF reads aligned
 * @param readsPaired        Reads aligned in pairs (READS_ALIGNED_IN_PAIRS)
 * @param pctReadsPaired     Percentage of aligned reads that are paired
 * @param pctProperPairs     Percentage of reads aligned as proper pairs
 * @param readLength         Mean read length (MEAN_READ_LENGTH)
 * @param maxReadLength      Maximum read length (used for GATK --READ_LENGTH parameter)
 * @param meanInsertSize     Mean insert size of the library (from CollectInsertSizeMetrics)
 * @param medianInsertSize   Median insert size (more robust than mean)
 * @param stdInsertSize      Standard deviation of insert size
 * @param pairOrientation    Read pair orientation (FR, RF, TANDEM)
 * @param flowcellId         Flowcell identifier if available from headers
 * @param runDate            Date of the sequencing run if extractable
 * @param files              Metadata about raw data files (e.g., FASTQs). Files remain local.
 * @param alignmentRefs      AT URIs of alignment records derived from this sequence run
 */
case class SequenceRun(
  atUri: Option[String],
  meta: RecordMeta,
  biosampleRef: String,
  platformName: String,
  instrumentModel: Option[String] = None,
  instrumentId: Option[String] = None,
  testType: String,
  libraryLayout: Option[String] = None,
  totalReads: Option[Long] = None,
  pfReads: Option[Long] = None,
  pfReadsAligned: Option[Long] = None,
  pctPfReadsAligned: Option[Double] = None,
  readsPaired: Option[Long] = None,
  pctReadsPaired: Option[Double] = None,
  pctProperPairs: Option[Double] = None,
  readLength: Option[Int] = None,
  maxReadLength: Option[Int] = None,
  meanInsertSize: Option[Double] = None,
  medianInsertSize: Option[Double] = None,
  stdInsertSize: Option[Double] = None,
  pairOrientation: Option[String] = None,
  flowcellId: Option[String] = None,
  runDate: Option[LocalDateTime] = None,
  files: List[FileInfo] = List.empty,
  alignmentRefs: List[String] = List.empty
)

object SequenceRun {
  /** Known platform values */
  val KnownPlatforms: Set[String] = Set(
    "ILLUMINA", "PACBIO", "NANOPORE", "ION_TORRENT", "BGI", "ELEMENT", "ULTIMA"
  )

  /**
   * Known test type codes (aligned with TestTypes).
   * Use TestTypes.byCode() to get full TestTypeDefinition.
   */
  val KnownTestTypes: Set[String] = Set(
    // Whole genome - short read
    "WGS", "WGS_LOW_PASS",
    // Whole genome - long read
    "WGS_HIFI", "WGS_NANOPORE", "WGS_CLR",
    // Exome
    "WES",
    // Targeted Y-DNA
    "BIG_Y_500", "BIG_Y_700", "Y_ELITE", "Y_PRIME",
    // Targeted mtDNA
    "MT_FULL_SEQUENCE", "MT_PLUS", "MT_CR_ONLY",
    // Legacy compatibility
    "EXOME", "TARGETED", "RNA_SEQ", "AMPLICON", "Unknown"
  )

  /** Known library layout values */
  val KnownLibraryLayouts: Set[String] = Set("PAIRED", "SINGLE")

  /**
   * Get display name for a test type code.
   */
  def testTypeDisplayName(code: String): String = {
    import com.decodingus.genotype.model.TestTypes
    TestTypes.byCode(code).map(_.displayName).getOrElse(code)
  }

  /**
   * Check if test type supports Y-DNA haplogroup analysis.
   */
  def supportsYDna(testTypeCode: String): Boolean = {
    import com.decodingus.genotype.model.TestTypes
    TestTypes.byCode(testTypeCode).exists(_.supportsHaplogroupY)
  }

  /**
   * Check if test type supports mtDNA haplogroup analysis.
   */
  def supportsMtDna(testTypeCode: String): Boolean = {
    import com.decodingus.genotype.model.TestTypes
    TestTypes.byCode(testTypeCode).exists(_.supportsHaplogroupMt)
  }

  /**
   * Check if test type supports ancestry analysis.
   */
  def supportsAncestry(testTypeCode: String): Boolean = {
    import com.decodingus.genotype.model.TestTypes
    TestTypes.byCode(testTypeCode).exists(_.supportsAncestry)
  }
}
