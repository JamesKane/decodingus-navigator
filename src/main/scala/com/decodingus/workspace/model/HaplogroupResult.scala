package com.decodingus.workspace.model

case class HaplogroupResult(
  haplogroupName: String,
  score: Double,
  matchingSnps: Option[Int],
  mismatchingSnps: Option[Int],
  ancestralMatches: Option[Int],
  treeDepth: Option[Int],
  lineagePath: Option[List[String]]
)
