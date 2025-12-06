package com.decodingus.workspace.model

case class Project(
  projectName: String,
  atUri: Option[String],
  description: Option[String],
  administrator: String,
  members: List[String]
)
