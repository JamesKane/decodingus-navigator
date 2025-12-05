package com.decodingus.workspace.model

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
