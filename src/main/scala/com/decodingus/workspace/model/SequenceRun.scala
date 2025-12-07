package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * A sequencing run representing one library preparation and sequencing session.
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.sequencerun).
 *
 * Can be independently managed - created, updated, or deleted without cascading changes.
 *
 * @param atUri           The AT URI of this sequence run record
 * @param meta            Record metadata for tracking changes and sync
 * @param biosampleRef    AT URI of the parent biosample record
 * @param platformName    Sequencing platform (ILLUMINA, PACBIO, NANOPORE, ION_TORRENT, BGI, ELEMENT, ULTIMA)
 * @param instrumentModel Specific instrument model (e.g., NovaSeq 6000, Sequel II)
 * @param instrumentId    Unique instrument identifier extracted from @RG headers (for lab inference)
 * @param testType        Type of test (WGS, EXOME, TARGETED, RNA_SEQ, AMPLICON)
 * @param libraryLayout   Paired-end or Single-end (PAIRED, SINGLE)
 * @param totalReads      Total number of reads
 * @param readLength      Average read length
 * @param meanInsertSize  Mean insert size of the library
 * @param flowcellId      Flowcell identifier if available from headers
 * @param runDate         Date of the sequencing run if extractable
 * @param files           Metadata about raw data files (e.g., FASTQs). Files remain local.
 * @param alignmentRefs   AT URIs of alignment records derived from this sequence run
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
  readLength: Option[Int] = None,
  meanInsertSize: Option[Double] = None,
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

  /** Known test type values */
  val KnownTestTypes: Set[String] = Set(
    "WGS", "EXOME", "TARGETED", "RNA_SEQ", "AMPLICON"
  )

  /** Known library layout values */
  val KnownLibraryLayouts: Set[String] = Set("PAIRED", "SINGLE")
}
