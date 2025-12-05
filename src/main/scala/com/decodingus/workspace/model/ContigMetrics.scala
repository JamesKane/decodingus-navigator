package com.decodingus.workspace.model

case class ContigMetrics(
  contigName: String,
  callableBases: Long,
  meanCoverage: Option[Double],
  poorMappingQuality: Option[Long],
  lowCoverage: Option[Long],
  noCoverage: Option[Long]
)
