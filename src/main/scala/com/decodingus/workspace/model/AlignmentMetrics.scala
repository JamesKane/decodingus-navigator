package com.decodingus.workspace.model

case class AlignmentMetrics(
  genomeTerritory: Option[Long],
  meanCoverage: Option[Double],
  medianCoverage: Option[Double],
  sdCoverage: Option[Double],
  pctExcDupe: Option[Double],
  pctExcMapq: Option[Double],
  pct10x: Option[Double],
  pct20x: Option[Double],
  pct30x: Option[Double],
  hetSnpSensitivity: Option[Double],
  contigs: List[ContigMetrics]
)
