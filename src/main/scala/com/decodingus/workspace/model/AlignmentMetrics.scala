package com.decodingus.workspace.model

case class AlignmentMetrics(
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
  callableBases: Option[Long] = None,
  contigs: List[ContigMetrics] = List.empty
)
