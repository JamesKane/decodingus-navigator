package com.decodingus.workspace.model

/**
 * DEPRECATED: Use Alignment instead.
 *
 * This class is maintained for backward compatibility with UI components.
 * It represents the legacy embedded alignment data model. New code should use
 * Alignment which is a first-class record with its own AT URI.
 *
 * @see Alignment for the new first-class record type
 */
@deprecated("Use Alignment instead - this embedded model is being phased out", "2.0")
case class AlignmentData(
  referenceBuild: String,
  aligner: String,
  files: List[FileInfo],
  metrics: Option[AlignmentMetrics]
)

object AlignmentData {
  /** Convert an Alignment record to legacy AlignmentData for UI compatibility */
  def fromAlignment(alignment: Alignment): AlignmentData = {
    AlignmentData(
      referenceBuild = alignment.referenceBuild,
      aligner = alignment.aligner,
      files = alignment.files,
      metrics = alignment.metrics
    )
  }

  /** Convert legacy AlignmentData to an Alignment record */
  def toAlignment(data: AlignmentData, sequenceRunRef: String, biosampleRef: Option[String], meta: RecordMeta): Alignment = {
    Alignment(
      atUri = Some(s"local:alignment:${java.util.UUID.randomUUID()}"),
      meta = meta,
      sequenceRunRef = sequenceRunRef,
      biosampleRef = biosampleRef,
      referenceBuild = data.referenceBuild,
      aligner = data.aligner,
      files = data.files,
      metrics = data.metrics
    )
  }
}
