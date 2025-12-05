package com.decodingus.workspace.model

case class Project(
  projectName: String,
  description: Option[String],
  administrator: String,
  members: List[String]
)
