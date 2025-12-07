package com.decodingus.workspace.model

/**
 * DEPRECATED: Use SequenceRun instead.
 *
 * This class is maintained for backward compatibility with UI components.
 * It represents the legacy embedded sequence data model. New code should use
 * SequenceRun which is a first-class record with its own AT URI.
 *
 * @see SequenceRun for the new first-class record type
 */
@deprecated("Use SequenceRun instead - this embedded model is being phased out", "2.0")
case class SequenceData(
  platformName: String,
  instrumentModel: Option[String],
  testType: String,
  libraryLayout: Option[String],
  totalReads: Option[Long],
  readLength: Option[Int],
  meanInsertSize: Option[Double],
  files: List[FileInfo],
  alignments: List[AlignmentData]
)

object SequenceData {
  /** Convert a SequenceRun to legacy SequenceData for UI compatibility */
  def fromSequenceRun(run: SequenceRun, alignmentRecords: List[Alignment]): SequenceData = {
    SequenceData(
      platformName = run.platformName,
      instrumentModel = run.instrumentModel,
      testType = run.testType,
      libraryLayout = run.libraryLayout,
      totalReads = run.totalReads,
      readLength = run.readLength,
      meanInsertSize = run.meanInsertSize,
      files = run.files,
      alignments = alignmentRecords.map(AlignmentData.fromAlignment)
    )
  }

  /** Convert legacy SequenceData to a SequenceRun record */
  def toSequenceRun(data: SequenceData, biosampleRef: String, meta: RecordMeta): (SequenceRun, List[Alignment]) = {
    val seqRunUri = s"local:sequencerun:${java.util.UUID.randomUUID()}"

    val alignmentPairs = data.alignments.zipWithIndex.map { case (alignData, idx) =>
      AlignmentData.toAlignment(alignData, seqRunUri, Some(biosampleRef), meta)
    }

    val sequenceRun = SequenceRun(
      atUri = Some(seqRunUri),
      meta = meta,
      biosampleRef = biosampleRef,
      platformName = data.platformName,
      instrumentModel = data.instrumentModel,
      testType = data.testType,
      libraryLayout = data.libraryLayout,
      totalReads = data.totalReads,
      readLength = data.readLength,
      meanInsertSize = data.meanInsertSize,
      files = data.files,
      alignmentRefs = alignmentPairs.flatMap(_.atUri)
    )

    (sequenceRun, alignmentPairs)
  }
}
