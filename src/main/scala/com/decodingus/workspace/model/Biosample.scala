package com.decodingus.workspace.model

import java.time.LocalDateTime

case class Biosample(
  sampleAccession: String,
  donorIdentifier: String,
  description: Option[String],
  centerName: Option[String],
  sex: Option[String],
  sequenceData: List[SequenceData],
  haplogroups: Option[HaplogroupAssignments],
  createdAt: Option[LocalDateTime]
)
