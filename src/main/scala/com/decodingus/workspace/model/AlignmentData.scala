package com.decodingus.workspace.model

case class AlignmentData(
  referenceBuild: String,
  aligner: String,
  files: List[FileInfo],
  metrics: Option[AlignmentMetrics]
)
